#!/usr/bin/env python3
"""
Split a folder of mixed DICOM files into one subfolder per series.

  - Groups by SeriesInstanceUID (the only reliable series key; SeriesNumber or
    SeriesDescription can repeat, e.g. two PET recons both called "PET WB").
  - Subfolder name follows the existing PETCTDCM convention:
      <SeriesNumber><SERIESDESCRIPTION alnum, upper><Modality>
      e.g. 'PET WB (OncoFreeze AI)', #10, PT -> 10PETWBONCOFREEZEAIPT
           'CT WB 3.0 Br38', #3, CT          -> 3CTWB30BR38CT
    Only if two different series map to the same name is '_<last 6 of UID>' added (reported).
  - Original file names are kept as-is.
  - Exact duplicates (same SOPInstanceUID, identical bytes) are written once (first file name kept).
  - Same SOPInstanceUID with DIFFERENT content is reported as a conflict and both files kept.
  - Only if two different files would land on the same name, '_1', '_2'... is appended.
  - Copies by default (source untouched); --move to move instead.

Usage
  python split_dicom_series.py <input_dir> <output_dir>
  python split_dicom_series.py F:/PETCTDCM/SUBJ01 F:/PETCTDCM_split/SUBJ01 --move
"""
from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import pydicom


def file_md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def folder_name(num, mod, desc) -> str:
    return f"{num}{re.sub(r'[^A-Za-z0-9]+', '', str(desc)).upper()}{mod}"


def main() -> None:
    ap = argparse.ArgumentParser(description="Split DICOM folder by SeriesInstanceUID")
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--move", action="store_true", help="move files instead of copying")
    a = ap.parse_args()
    if not a.input.is_dir():
        sys.exit(f"Input folder not found: {a.input}")

    series = defaultdict(dict)          # series_uid -> {sop_uid: [paths]}
    meta = {}
    skipped = 0
    for p in sorted(a.input.rglob("*")):
        if not p.is_file():
            continue
        try:
            ds = pydicom.dcmread(p, stop_before_pixels=True, force=True)
            suid, sop = str(ds.SeriesInstanceUID), str(ds.SOPInstanceUID)
        except Exception:
            skipped += 1
            continue
        series[suid].setdefault(sop, []).append(p)
        if suid not in meta:
            meta[suid] = (getattr(ds, "SeriesNumber", "NA"), getattr(ds, "Modality", "NA"),
                          getattr(ds, "SeriesDescription", ""))

    if not series:
        sys.exit("No DICOM files found")

    names = defaultdict(list)
    for suid, (num, mod, desc) in meta.items():
        names[folder_name(num, mod, desc)].append(suid)

    print(f"{'folder':<60} {'images':>7} {'dupes':>6} {'conflicts':>9}")
    for suid, sops in sorted(series.items(), key=lambda kv: str(meta[kv[0]][0])):
        num, mod, desc = meta[suid]
        base = folder_name(num, mod, desc)
        if len(names[base]) > 1:                      # two different series, same name
            base = f"{base}_{suid[-6:]}"
            print(f"[!] folder name clash, UID tail added: {base}")
        dest = a.output / base
        dest.mkdir(parents=True, exist_ok=True)
        n_dup = n_conf = 0
        for sop, paths in sops.items():
            groups = defaultdict(list)                   # content hash -> paths
            for p in paths:
                groups[file_md5(p) if len(paths) > 1 else ""].append(p)
            n_dup += sum(len(g) - 1 for g in groups.values())
            if len(groups) > 1:
                n_conf += 1
            for h, g in groups.items():
                target = dest / g[0].name                # keep original file name
                k = 1
                while target.exists():                   # name clash only (e.g. same name in subfolders)
                    target = dest / f"{g[0].stem}_{k}{g[0].suffix}"
                    k += 1
                (shutil.move if a.move else shutil.copy2)(g[0], target)
                if a.move:                               # drop the redundant identical copies
                    for extra in g[1:]:
                        extra.unlink()
        print(f"{dest.name:<60} {len(sops):>7} {n_dup:>6} {n_conf:>9}")

    if skipped:
        print(f"\n{skipped} non-DICOM files ignored")


if __name__ == "__main__":
    main()
