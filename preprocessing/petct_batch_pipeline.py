#!/usr/bin/env python3
"""
PET/CT batch pipeline (one row = one subject)

  1. CT  : DICOM --(dcm2niix)--> ct   (CT_*.nii.gz)
  2. SUV : PET DICOM --(pysuv, SUVbw)--> suv  (SUV_*.nii.gz, native PET grid, float32)
  3. PET : suv --(resample to CT grid, identity transform)--> pet  (PET_*.nii.gz, float32)

Input CSV columns (required): subjid, ctdcm, petdcm, output_dir, ct, suv, pet
  ../prep/petctdcm.csv  (override with --csv)

Outputs
  - ct / suv / pet paths given in the CSV
  - ../logs/petct_error.csv   : failed subjects (step + message)
  - ../logs/petct_qc.csv      : SUV parameters + warnings for every processed subject

Design notes
  - dcm2niix runs in a per-subject temp folder inside output_dir; only the chosen CT
    file is moved to its final name, so no ROI / Eq / json leftovers ever reach output_dir.
  - dcm2niix ROI overlay outputs (_ROIn) are ignored with a warning.
    Any other extra output (_Eq uneven spacing/missing slice, _Tilt gantry tilt,
    additional series) stops that subject; no CT/SUV/PET is written for it.
  - No json sidecar (-b n).
  - Default: subjects whose ct, suv and pet outputs all exist are skipped.
    --overwrite reprocesses everything. A subject with only some outputs is always redone.
  - SUV computation delegated to pysuv (validated vs QIBA DRO; matches LIFEx v25.06).
"""
from __future__ import annotations

import argparse
import gc
import os
import re
import shutil
import subprocess
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import psutil
import pydicom
import pysuv
import SimpleITK as sitk
from tqdm import tqdm


def get_script_dir() -> Path:
    try:
        return Path(__file__).resolve().parent
    except NameError:
        return Path(sys.argv[0]).resolve().parent


# === Paths ===
SCRIPT_DIR      = get_script_dir()
CSV_INPUT_PATH  = SCRIPT_DIR.parent / "prep" / "petctdcm.csv"
LOGS_DIR        = SCRIPT_DIR.parent / "logs"
ERROR_LOG_PATH  = LOGS_DIR / "petct_error.csv"
QC_LOG_PATH     = LOGS_DIR / "petct_qc.csv"

REQUIRED_COLS   = ["subjid", "ctdcm", "petdcm", "output_dir", "ct", "suv", "pet"]
DCM2NIIX_TIMEOUT_S = 900
GB = 1024 ** 3


class StepError(Exception):
    """Error raised inside a pipeline step; carries the step name for logging."""
    def __init__(self, step: str, msg: str):
        super().__init__(msg)
        self.step = step


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def format_exception_msg(err: Exception) -> str:
    msg = str(err)
    if "Exception thrown in SimpleITK" in msg or "itk::error" in msg.lower():
        m = re.search(r"itk::ERROR:\s*(?:[A-Za-z0-9_]+\([^)]*\):\s*)?([^\n\r]+)", msg, re.IGNORECASE)
        if m:
            return f"SimpleITK Error: {m.group(1).strip()}"
    return " ".join(msg.split())


def nii_stem(path: Path) -> str:
    name = path.name
    for ext in (".nii.gz", ".nii"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return path.stem


def tmp_sibling(final_path: Path) -> Path:
    """Temp file in the same folder (same volume -> os.replace is atomic); keeps .nii.gz suffix."""
    return final_path.parent / f".tmp_{os.getpid()}_{final_path.name}"


def read_header(path: Path):
    try:
        return pydicom.dcmread(path, stop_before_pixels=True, force=True)
    except Exception:
        return None


def first_dicom_header(folder: Path, modality: str | None = None):
    for p in sorted(folder.rglob("*")):
        if not p.is_file():
            continue
        ds = read_header(p)
        if ds is None or "SOPClassUID" not in ds:
            continue
        if modality and getattr(ds, "Modality", None) != modality:
            continue
        return ds
    return None


# --------------------------------------------------------------------------- #
# Step 1: CT DICOM -> NIfTI (dcm2niix)
# --------------------------------------------------------------------------- #
def convert_ct(ctdcm: Path, ct_out: Path, warn: list[str]) -> None:
    step = "ct"
    if not ctdcm.is_dir():
        raise StepError(step, f"CT DICOM directory not found: {ctdcm}")

    ct_out.parent.mkdir(parents=True, exist_ok=True)
    base = nii_stem(ct_out)
    tmp_dir = ct_out.parent / f".tmp_dcm2niix_{os.getpid()}"
    shutil.rmtree(tmp_dir, ignore_errors=True)
    tmp_dir.mkdir(parents=True)

    try:
        cmd = ["dcm2niix", "-b", "n", "-z", "y", "-w", "1",
               "-f", base, "-o", str(tmp_dir), str(ctdcm)]
        try:
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, shell=False, timeout=DCM2NIIX_TIMEOUT_S)
        except subprocess.TimeoutExpired:
            raise StepError(step, f"dcm2niix timed out after {DCM2NIIX_TIMEOUT_S}s")
        if res.returncode != 0:
            msg = " ".join((res.stderr.strip() or res.stdout.strip()).split())
            raise StepError(step, f"dcm2niix failed (code {res.returncode}): {msg[:300]}")

        # Keep dcm2niix's own warnings (e.g. slice spacing details, tilt) even on success
        seen = set()
        for line in (res.stdout + "\n" + res.stderr).splitlines():
            line = " ".join(line.split())
            low = line.lower()
            if ("warning" in low or "duplicate" in low) and line not in seen:
                seen.add(line)
                warn.append(f"dcm2niix: {line[:250]}")

        produced = sorted(p for p in tmp_dir.glob("*.nii.gz") if p.stat().st_size > 0)
        if not produced:
            raise StepError(step, "dcm2niix produced no NIfTI output")

        # Policy: ROI overlays (_ROIn) -> warning only, ignored.
        #         Anything else besides the single expected file (_Eq, _Tilt, extra series,
        #         unexpected suffix) -> stop this subject. Nothing is written to output_dir.
        def suffix(p: Path) -> str:
            return nii_stem(p)[len(base):]

        roi = [p for p in produced if re.match(r"^_ROI\d*", suffix(p))]
        rest = [p for p in produced if p not in roi]
        if roi:
            warn.append(f"CT ROI overlay outputs ignored: {[p.name for p in roi]}")

        exact = [p for p in rest if suffix(p) == ""]
        others = [p for p in rest if suffix(p) != ""]
        if others:
            kinds = []
            if any("_Eq" in suffix(p) for p in others):
                kinds.append("uneven slice spacing/missing slices (_Eq)")
            if any("_Tilt" in suffix(p) for p in others):
                kinds.append("gantry tilt (_Tilt)")
            if not kinds:
                kinds.append("unexpected extra outputs")
            raise StepError(step, f"CT conversion stopped: {', '.join(kinds)}; "
                                  f"dcm2niix outputs: {[p.name for p in rest]}")
        if len(exact) != 1:
            raise StepError(step, f"CT conversion stopped: expected 1 output, got {[p.name for p in rest]}")
        chosen = exact[0]

        os.replace(chosen, ct_out)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# --------------------------------------------------------------------------- #
# Step 2: PET DICOM -> SUVbw NIfTI (pysuv)
# --------------------------------------------------------------------------- #
def collect_pet_slices(petdcm: Path, warn: list[str]) -> list[tuple[Path, pydicom.Dataset]]:
    step = "suv"
    if not petdcm.is_dir():
        raise StepError(step, f"PET DICOM directory not found: {petdcm}")

    slices = []
    for p in petdcm.rglob("*"):
        if not p.is_file():
            continue
        ds = read_header(p)
        if ds is None or getattr(ds, "Modality", None) != "PT":
            continue
        if "ImagePositionPatient" not in ds or "ImageOrientationPatient" not in ds:
            continue
        slices.append((p, ds))
    if not slices:
        raise StepError(step, f"No PET (Modality=PT) images in {petdcm}")

    # Drop exact duplicates (same SOPInstanceUID exported twice under different file names).
    # Files on disk are never renamed/moved; de-duplication is done in memory only.
    by_uid: dict[str, tuple[Path, pydicom.Dataset]] = {}
    dup_files = []
    for p, ds in sorted(slices, key=lambda t: t[0].name):
        uid = str(getattr(ds, "SOPInstanceUID", "")) or str(p)
        if uid in by_uid:
            dup_files.append(p.name)
        else:
            by_uid[uid] = (p, ds)
    if dup_files:
        warn.append(f"PET duplicate files ignored (same SOPInstanceUID): {len(dup_files)} "
                    f"e.g. {dup_files[:3]}")
    slices = list(by_uid.values())

    series: dict[str, list] = {}
    for p, ds in slices:
        series.setdefault(str(ds.SeriesInstanceUID), []).append((p, ds))
    if len(series) > 1:
        desc = [f"#{getattr(v[0][1], 'SeriesNumber', '?')} '{getattr(v[0][1], 'SeriesDescription', '')}' "
                f"({len(v)} files)" for v in series.values()]
        raise StepError(step, f"{len(series)} PET series in one folder: {'; '.join(desc)}")

    iop = np.array(slices[0][1].ImageOrientationPatient, dtype=float)
    normal = np.cross(iop[:3], iop[3:])
    slices.sort(key=lambda t: float(np.dot(np.array(t[1].ImagePositionPatient, dtype=float), normal)))
    return slices


def check_pet_consistency(slices, warn: list[str]) -> None:
    step = "suv"
    ref = slices[0][1]

    def same(attr_fn, name):
        vals = {attr_fn(ds) for _, ds in slices}
        if len(vals) > 1:
            raise StepError(step, f"Inconsistent {name} across PET slices: {sorted(map(str, vals))[:5]}")

    same(lambda d: (int(d.Rows), int(d.Columns)), "matrix size")
    same(lambda d: tuple(round(float(x), 5) for x in d.PixelSpacing), "PixelSpacing")
    same(lambda d: tuple(round(float(x), 5) for x in d.ImageOrientationPatient), "ImageOrientationPatient")
    same(lambda d: str(getattr(d, "Units", "")), "Units")
    same(lambda d: str(getattr(d, "DecayCorrection", "")), "DecayCorrection")
    same(lambda d: str(getattr(d, "PatientWeight", "")), "PatientWeight")

    def rp(d, key):
        try:
            return str(getattr(d.RadiopharmaceuticalInformationSequence[0], key, ""))
        except Exception:
            return ""
    same(lambda d: rp(d, "RadionuclideTotalDose"), "RadionuclideTotalDose")
    same(lambda d: rp(d, "RadiopharmaceuticalStartTime"), "RadiopharmaceuticalStartTime")

    # Uniform slice spacing (affine assumes it)
    if len(slices) > 1:
        iop = np.array(ref.ImageOrientationPatient, dtype=float)
        normal = np.cross(iop[:3], iop[3:])
        pos = np.array([np.dot(np.array(ds.ImagePositionPatient, dtype=float), normal) for _, ds in slices])
        d = np.diff(pos)
        if np.any(d <= 1e-4):
            raise StepError(step, "Duplicate PET slice positions detected")
        if not np.allclose(d, d.mean(), atol=0.01):
            raise StepError(step, f"Non-uniform PET slice spacing: min {d.min():.3f}, max {d.max():.3f} mm")

    n_expected = getattr(ref, "NumberOfSlices", None)
    if n_expected is not None and int(n_expected) != len(slices):
        warn.append(f"PET NumberOfSlices tag={int(n_expected)} but {len(slices)} files found")


def build_affine(slices) -> np.ndarray:
    ds0 = slices[0][1]
    iop = np.array(ds0.ImageOrientationPatient, dtype=float)
    row_cos, col_cos = iop[:3], iop[3:]
    dr, dc = (float(x) for x in ds0.PixelSpacing)
    ipp0 = np.array(ds0.ImagePositionPatient, dtype=float)
    if len(slices) > 1:
        ippN = np.array(slices[-1][1].ImagePositionPatient, dtype=float)
        slice_vec = (ippN - ipp0) / (len(slices) - 1)
    else:
        slice_vec = np.cross(row_cos, col_cos) * float(getattr(ds0, "SliceThickness", 1.0) or 1.0)

    aff = np.eye(4)
    aff[:3, 0] = row_cos * dc
    aff[:3, 1] = col_cos * dr
    aff[:3, 2] = slice_vec
    aff[:3, 3] = ipp0
    return np.diag([-1.0, -1.0, 1.0, 1.0]) @ aff   # LPS -> RAS


def convert_suv(petdcm: Path, suv_out: Path, warn: list[str], qc: dict) -> pydicom.Dataset:
    step = "suv"
    slices = collect_pet_slices(petdcm, warn)
    check_pet_consistency(slices, warn)
    ref_path, ref = slices[0]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            info = pysuv.get_suv_parameters(str(ref_path))
        except Exception as e:
            raise StepError(step, f"pysuv parameter extraction failed: {format_exception_msg(e)}")

        qc.update({
            "manufacturer": getattr(info, "manufacturer", None),
            "model": str(getattr(ref, "ManufacturerModelName", "")),
            "units": getattr(info, "units", None),
            "decay_correction": getattr(info, "decay_correction", None),
            "injected_dose_mbq": getattr(info, "injected_dose_mbq", None),
            "patient_weight_kg": getattr(info, "patient_weight_kg", None),
            "half_life_min": getattr(info, "half_life_min", None),
            "injection_time": getattr(info, "injection_time", None),
            "scan_time": getattr(info, "scan_time", None),
            "uptake_min": getattr(info, "elapsed_time_min", None),
            "pet_matrix": f"{int(ref.Columns)}x{int(ref.Rows)}x{len(slices)}",
        })

        # Allocate float32 directly; read pixel data one slice at a time (low memory)
        nx, ny, nz = int(ref.Columns), int(ref.Rows), len(slices)
        suv = np.empty((nx, ny, nz), dtype=np.float32)
        for k, (p, _) in enumerate(slices):
            try:
                s = pysuv.compute_suv(str(p), normalization="bw")
            except Exception as e:
                raise StepError(step, f"pysuv failed on {p.name}: {format_exception_msg(e)}")
            suv[:, :, k] = np.asarray(s, dtype=np.float32).T   # (rows, cols) -> (x, y)

    for w in caught:
        msg = " ".join(str(w.message).split())
        if msg not in warn:
            warn.append(f"pysuv: {msg}")

    if not np.isfinite(suv).all():
        raise StepError(step, "Non-finite SUV values produced")
    smax = float(suv.max())
    qc["suv_max_native"] = round(smax, 4)
    if smax <= 0:
        raise StepError(step, "SUV volume is all zero")

    # QC heuristics (warn only)
    up = qc.get("uptake_min")
    if up is not None and (up < 0 or up > 300):
        warn.append(f"Unusual uptake time: {up:.1f} min")
    wt = qc.get("patient_weight_kg")
    if wt is not None and (wt < 20 or wt > 250):
        warn.append(f"Unusual patient weight: {wt} kg")

    affine = build_affine(slices)
    img = nib.Nifti1Image(suv, affine)
    img.header.set_data_dtype(np.float32)
    img.header.set_slope_inter(1.0, 0.0)
    img.set_qform(affine, code=1)
    img.set_sform(affine, code=1)
    img.header.set_xyzt_units("mm")
    img.header["intent_name"] = b"SUVbw g/mL"
    img.header["descrip"] = b"SUVbw [g/mL], float32, computed with pysuv"

    suv_out.parent.mkdir(parents=True, exist_ok=True)
    tmp = tmp_sibling(suv_out)
    try:
        nib.save(img, tmp)
        os.replace(tmp, suv_out)
    finally:
        if tmp.exists():
            tmp.unlink()
    del suv, img
    return ref


# --------------------------------------------------------------------------- #
# Step 3: Resample SUV to CT grid
# --------------------------------------------------------------------------- #
def resample_to_ct(suv_path: Path, ct_path: Path, pet_out: Path, qc: dict) -> None:
    step = "resample"
    pet_img = sitk.ReadImage(str(suv_path))
    ct_img = sitk.ReadImage(str(ct_path))

    if pet_img.GetDimension() == 4:
        if pet_img.GetSize()[3] != 1:
            raise StepError(step, f"Unsupported multi-frame 4D PET (frames={pet_img.GetSize()[3]})")
        pet_img = pet_img[:, :, :, 0]
    if pet_img.GetDimension() != 3:
        raise StepError(step, f"Invalid PET dimension: {pet_img.GetDimension()}D")
    if ct_img.GetDimension() == 4 and ct_img.GetSize()[3] == 1:
        ct_img = ct_img[:, :, :, 0]
    if ct_img.GetDimension() != 3:
        raise StepError(step, f"Invalid CT dimension: {ct_img.GetDimension()}D")

    stat = sitk.MinimumMaximumImageFilter()
    stat.Execute(pet_img)
    orig_max = stat.GetMaximum()

    res = sitk.Resample(pet_img, ct_img, sitk.Transform(), sitk.sitkLinear, 0.0, sitk.sitkFloat32)

    checks = [
        ("Size", res.GetSize() == ct_img.GetSize()),
        ("Spacing", np.allclose(res.GetSpacing(), ct_img.GetSpacing(), atol=1e-4)),
        ("Origin", np.allclose(res.GetOrigin(), ct_img.GetOrigin(), atol=1e-4)),
        ("Direction", np.allclose(res.GetDirection(), ct_img.GetDirection(), atol=1e-4)),
    ]
    for name, ok in checks:
        if not ok:
            raise StepError(step, f"Geometric mismatch in {name} between resampled PET and CT")

    stat.Execute(res)
    res_max = stat.GetMaximum()
    qc["ct_matrix"] = "x".join(map(str, ct_img.GetSize()))
    qc["suv_max_resampled"] = round(float(res_max), 4)
    if orig_max > 0 and res_max <= 0:
        raise StepError(step, "Zero spatial overlap: resampled PET is entirely blank")

    del pet_img, ct_img
    pet_out.parent.mkdir(parents=True, exist_ok=True)
    tmp = tmp_sibling(pet_out)
    try:
        sitk.WriteImage(res, str(tmp), useCompression=True)
        if not tmp.is_file() or tmp.stat().st_size == 0:
            raise StepError(step, "Resampled output is missing or 0 bytes")
        os.replace(tmp, pet_out)
    finally:
        if tmp.exists():
            tmp.unlink()
    del res


# --------------------------------------------------------------------------- #
# Worker: one subject, all three steps
# --------------------------------------------------------------------------- #
def process_subject(row: dict) -> tuple[dict, dict | None]:
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)
    warn: list[str] = []
    qc: dict = {"subjid": row.get("subjid")}
    step = "input"

    try:
        for k in REQUIRED_COLS:
            v = row.get(k)
            if v is None or pd.isna(v) or str(v).strip() == "":
                raise StepError("input", f"Missing value in column '{k}'")

        ctdcm, petdcm = Path(str(row["ctdcm"]).strip()), Path(str(row["petdcm"]).strip())
        ct_out = Path(str(row["ct"]).strip())
        suv_out = Path(str(row["suv"]).strip())
        pet_out = Path(str(row["pet"]).strip())

        step = "ct"
        convert_ct(ctdcm, ct_out, warn)

        step = "suv"
        pet_ref = convert_suv(petdcm, suv_out, warn, qc)
        gc.collect()

        # Same scanner coordinate system is required for identity-transform resampling
        step = "check"
        ct_ref = first_dicom_header(ctdcm, modality="CT")
        if ct_ref is not None:
            for_ct = str(getattr(ct_ref, "FrameOfReferenceUID", ""))
            for_pet = str(getattr(pet_ref, "FrameOfReferenceUID", ""))
            if not for_ct or not for_pet:
                warn.append("FrameOfReferenceUID missing in CT or PET")
            elif for_ct != for_pet:
                warn.append("CT and PET FrameOfReferenceUID differ (identity resampling may be misaligned)")
            if str(getattr(ct_ref, "StudyInstanceUID", "")) != str(getattr(pet_ref, "StudyInstanceUID", "")):
                warn.append("CT and PET StudyInstanceUID differ")

        step = "resample"
        resample_to_ct(suv_out, ct_out, pet_out, qc)

        qc["status"] = "ok"
        qc["warnings"] = " | ".join(warn)
        return qc, None

    except Exception as e:
        err_step = e.step if isinstance(e, StepError) else step
        msg = str(e) if isinstance(e, StepError) else format_exception_msg(e)
        qc["status"] = f"failed:{err_step}"
        qc["warnings"] = " | ".join(warn)
        return qc, {**row, "step": err_step, "error_msg": msg}
    finally:
        gc.collect()


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def default_workers() -> int:
    """
    Per subject peak RAM ~3-4 GB (high-res PET 440x440x~500 + CT-grid float32 resample).
    i5-14400 (6P+4E, 16 threads), 32 GB RAM, SSD->SSD: 4 workers is the sweet spot.
    Also capped by currently available memory (4 GB per worker).
    """
    phys = psutil.cpu_count(logical=False) or os.cpu_count() or 2
    avail_gb = psutil.virtual_memory().available / GB
    return max(1, min(4, phys - 2, int(avail_gb // 4)))


def make_executor(n: int) -> ProcessPoolExecutor:
    # Fresh process per subject releases ITK/numpy memory fully (Python >= 3.11)
    try:
        return ProcessPoolExecutor(max_workers=n, max_tasks_per_child=1)
    except TypeError:
        return ProcessPoolExecutor(max_workers=n)


def main() -> None:
    ap = argparse.ArgumentParser(description="Batch PET/CT: CT dcm2niix + PET SUVbw + resample to CT")
    ap.add_argument("--csv", type=Path, default=CSV_INPUT_PATH, help=f"input CSV (default {CSV_INPUT_PATH})")
    ap.add_argument("--workers", type=int, default=None, help="number of worker processes (default: auto)")
    ap.add_argument("--overwrite", action="store_true",
                    help="reprocess every subject. Default: skip subjects whose ct, suv and pet "
                         "outputs all already exist (non-empty)")
    args = ap.parse_args()

    if not args.csv.exists():
        sys.exit(f"[ERROR] Input CSV not found: {args.csv}")
    if not shutil.which("dcm2niix"):
        sys.exit("[ERROR] 'dcm2niix' not found in PATH")

    df = pd.read_csv(args.csv, dtype={"subjid": str})
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        sys.exit(f"[ERROR] CSV missing columns: {missing}")

    n_workers = args.workers or default_workers()
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    all_rows = df.to_dict("records")

    def done(r: dict) -> bool:
        paths = [r.get(k) for k in ("ct", "suv", "pet")]
        if any(p is None or pd.isna(p) or str(p).strip() == "" for p in paths):
            return False
        return all(os.path.isfile(str(p).strip()) and os.path.getsize(str(p).strip()) > 0 for p in paths)

    if args.overwrite:
        rows, skipped = all_rows, []
    else:
        rows = [r for r in all_rows if not done(r)]
        skipped = [r for r in all_rows if done(r)]

    # Carry previous QC rows for skipped subjects so the QC table stays complete
    prev_qc = {}
    if skipped and QC_LOG_PATH.exists():
        try:
            old = pd.read_csv(QC_LOG_PATH, dtype={"subjid": str})
            prev_qc = {str(r["subjid"]): r for r in old.to_dict("records")}
        except Exception:
            prev_qc = {}

    print(f"Subjects: {len(all_rows)} | to process: {len(rows)} | skipped (outputs exist): {len(skipped)} | "
          f"overwrite: {args.overwrite} | workers: {n_workers} | "
          f"RAM available: {psutil.virtual_memory().available / GB:.1f} GB")

    qc_records, error_records = [], []
    for r in skipped:
        rec = dict(prev_qc.get(str(r.get("subjid")), {"subjid": r.get("subjid")}))
        prev_status = rec.get("status")
        rec["status"] = f"skipped (previous: {prev_status})" if isinstance(prev_status, str) else "skipped"
        qc_records.append(rec)

    with make_executor(n_workers) as ex:
        futures = {ex.submit(process_subject, r): r for r in rows}
        with tqdm(total=len(futures), desc="PET/CT", unit="subj") as pbar:
            for fut in as_completed(futures):
                row = futures[fut]
                try:
                    qc, err = fut.result()
                    qc_records.append(qc)
                    if err:
                        error_records.append(err)
                except Exception as e:   # worker crash (e.g. out of memory)
                    error_records.append({**row, "step": "worker", "error_msg": format_exception_msg(e)})
                    qc_records.append({"subjid": row.get("subjid"), "status": "failed:worker"})
                pbar.update(1)
                pbar.set_postfix(failed=len(error_records))

    qc_df = pd.DataFrame(qc_records)
    if not qc_df.empty:
        qc_df = qc_df.sort_values("subjid")
    qc_df.to_csv(QC_LOG_PATH, index=False)

    n_warn = int((qc_df.get("warnings", pd.Series(dtype=str)).fillna("") != "").sum())
    print(f"[i] QC table: {QC_LOG_PATH} ({n_warn} subjects with warnings)")

    if error_records:
        pd.DataFrame(error_records).to_csv(ERROR_LOG_PATH, index=False)
        print(f"[!] {len(error_records)} failed subjects -> {ERROR_LOG_PATH}")
    else:
        if ERROR_LOG_PATH.exists():
            ERROR_LOG_PATH.unlink()
        print(f"[+] All processed subjects completed successfully ({len(skipped)} skipped).")


if __name__ == "__main__":
    main()
