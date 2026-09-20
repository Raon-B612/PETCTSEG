import os
import sys
import re
import gc
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import psutil
import argparse
import pandas as pd
import numpy as np
import SimpleITK as sitk
from tqdm import tqdm

def get_script_dir():
    try:
        # Standard script execution
        return Path(__file__).resolve().parent
    except NameError:
        # Fallback for interactive/notebook execution
        return Path(sys.argv[0]).resolve().parent

# === Paths and Environment Setup ===
SCRIPT_DIR        = Path(get_script_dir())
CSV_INPUT_PATH    = SCRIPT_DIR.parent / "prep" / "resample_suv_prep.csv"
ERROR_OUTPUT_PATH = SCRIPT_DIR.parent / "logs" / "resampling_error.csv"

# === Resampling Core Function ===
def resample_pet_to_ct(pet_image, ct_image):
    """
    Resample PET image to match CT physical grid (Origin, Spacing, Direction, Size)
    using an Identity Transform to preserve native scanner coordinates.
    """
    identity_transform = sitk.Transform()

    return sitk.Resample(
        pet_image,
        ct_image,
        identity_transform,
        sitk.sitkLinear,
        0.0,
        sitk.sitkFloat32
    )

# === Exception Formatting Helper ===
def format_exception_msg(err: Exception) -> str:
    """
    Extract core actionable description from SimpleITK/ITK C++ exceptions
    and sanitize multi-line traceback into a single-line string for CSV integrity.
    """
    msg = str(err)
    if "Exception thrown in SimpleITK" in msg or "itk::ERROR" in msg.lower():
        match = re.search(r"itk::ERROR:\s*(?:[A-Za-z0-9_]+\([^)]*\):\s*)?([^\n\r]+)", msg, re.IGNORECASE)
        if match:
            return f"SimpleITK Error: {match.group(1).strip()}"
    return " ".join(msg.split())

# === Single Row Worker Function ===
def process_row(task_args):
    """
    Worker task executing resampling pipeline with memory reclamation and atomic file writes.
    """
    row_idx, raw_suv, raw_ct, raw_pet = task_args

    # Initialize handles for guaranteed garbage collection in finally block
    pet_image = None
    ct_image = None
    resampled_pet = None
    stat_filter = None
    temp_output_path = None

    try:
        # Limit SimpleITK internal threads to 1 per worker to prevent CPU contention
        sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)

        # 1. Validate raw inputs before string casting to catch NaN and empty values
        if (
            pd.isna(raw_suv) or pd.isna(raw_ct) or pd.isna(raw_pet) or
            str(raw_suv).strip() == '' or str(raw_ct).strip() == '' or str(raw_pet).strip() == ''
        ):
            return row_idx, "Missing or invalid file path (suv, ct, or pet) in CSV row"

        pet_image_path = str(raw_suv).strip()
        ct_image_path = str(raw_ct).strip()
        output_path = str(raw_pet).strip()

        # 2. Check source file existence on disk
        if not os.path.isfile(pet_image_path):
            return row_idx, f"PET source file not found: {pet_image_path}"
        if not os.path.isfile(ct_image_path):
            return row_idx, f"CT source file not found: {ct_image_path}"

        # 3. Read images
        pet_image = sitk.ReadImage(pet_image_path)
        ct_image = sitk.ReadImage(ct_image_path)

        # 4. Handle 4D PET dimensions (squeeze singleton 4th axis to match 3D CT)
        if pet_image.GetDimension() == 4:
            if pet_image.GetSize()[3] == 1:
                pet_image = pet_image[:, :, :, 0]
            else:
                return row_idx, f"Unsupported multi-frame 4D PET (frames={pet_image.GetSize()[3]})"
        elif pet_image.GetDimension() != 3:
            return row_idx, f"Invalid PET dimension: {pet_image.GetDimension()}D (expected 3D)"

        if ct_image.GetDimension() != 3:
            return row_idx, f"Invalid CT dimension: {ct_image.GetDimension()}D (expected 3D)"

        # 5. Measure input PET maximum intensity
        stat_filter = sitk.MinimumMaximumImageFilter()
        stat_filter.Execute(pet_image)
        orig_max_val = stat_filter.GetMaximum()

        # 6. Resample directly to float32 in memory
        resampled_pet = resample_pet_to_ct(pet_image, ct_image)

        # 7. In-memory geometric verification against CT reference
        checks = [
            ("Size", resampled_pet.GetSize() == ct_image.GetSize()),
            ("Spacing", np.allclose(resampled_pet.GetSpacing(), ct_image.GetSpacing(), atol=1e-4)),
            ("Origin", np.allclose(resampled_pet.GetOrigin(), ct_image.GetOrigin(), atol=1e-4)),
            ("Direction", np.allclose(resampled_pet.GetDirection(), ct_image.GetDirection(), atol=1e-4))
        ]
        for name, passed in checks:
            if not passed:
                return row_idx, f"Geometric mismatch in {name} between resampled PET and CT"

        # 8. Detect silent spatial misalignment (all-zero resampled volume)
        stat_filter.Execute(resampled_pet)
        resampled_max_val = stat_filter.GetMaximum()

        if orig_max_val > 0.0 and resampled_max_val <= 0.0:
            return row_idx, "Zero physical spatial overlap detected: resampled PET is entirely blank"

        # 9. Atomic disk write via temporary file to prevent corrupted partial outputs
        target_dir = Path(output_path).parent
        target_dir.mkdir(parents=True, exist_ok=True)
        temp_output_path = target_dir / f".tmp_{os.getpid()}_{row_idx}_{Path(output_path).name}"

        sitk.WriteImage(resampled_pet, str(temp_output_path), useCompression=True)

        # Validate temporary file integrity prior to replacement
        if not temp_output_path.is_file() or temp_output_path.stat().st_size == 0:
            if temp_output_path.exists():
                temp_output_path.unlink()
            return row_idx, f"Output file creation failed or file is 0 bytes: {output_path}"

        # Atomic rename to final destination
        os.replace(temp_output_path, output_path)
        temp_output_path = None

        return row_idx, None

    except Exception as e:
        # Clean up dangling temporary file on failure
        if temp_output_path is not None and os.path.exists(temp_output_path):
            try:
                os.remove(temp_output_path)
            except OSError:
                pass
        return row_idx, format_exception_msg(e)

    finally:
        # Explicit memory deallocation of large SimpleITK C++ buffers per task
        del pet_image, ct_image, resampled_pet, stat_filter
        gc.collect()

# === Main Pipeline Execution ===
def main():
    parser = argparse.ArgumentParser(description="Resample PET to the CT grid.")
    parser.add_argument(
        "--hdd",
        action="store_true",
        help="Enable HDD mode: strictly caps worker processes at 4 to prevent disk thrashing."
    )
    args = parser.parse_args()

    # Determine CPU core availability with safe fallback
    physical_cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2

    # Storage-aware worker allocation
    if args.hdd:
        num_workers = min(4, max(1, physical_cores // 2))
        storage_mode = "HDD (Throttled, Max 4 Workers)"
    else:
        num_workers = max(1, physical_cores // 2)
        storage_mode = "SSD (High Throughput, Full Concurrency)"

    print(f"Active Mode: {storage_mode} -> Running with {num_workers} workers")

    if not CSV_INPUT_PATH.exists():
        print(f"[ERROR] Input CSV not found: {CSV_INPUT_PATH}")
        sys.exit(1)

    # Prepare log output directory
    ERROR_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    # Reset DataFrame index to guarantee 1:1 positional and label integrity
    df = pd.read_csv(CSV_INPUT_PATH).reset_index(drop=True)
    progress_bar = tqdm(total=len(df), desc='Resampling', unit='image')
    error_records = []

    # Extract primitive values via fast C-level zip to bypass pandas overhead
    tasks = list(zip(df.index, df['suv'], df['ct'], df['pet']))

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_row, task): task[0] for task in tasks}

        for future in as_completed(futures):
            row_idx = futures[future]
            try:
                idx, error = future.result()
                if error:
                    err_row = df.loc[idx].to_dict()
                    err_row['error_msg'] = error
                    error_records.append(err_row)
            except Exception as e:
                err_row = df.loc[row_idx].to_dict()
                err_row['error_msg'] = format_exception_msg(e)
                error_records.append(err_row)
            finally:
                progress_bar.update(1)

    progress_bar.close()

    if error_records:
        error_df = pd.DataFrame(error_records)
        error_df.to_csv(ERROR_OUTPUT_PATH, index=False)
        print(f"[!] Saved {len(error_records)} failed cases to {ERROR_OUTPUT_PATH}")
    else:
        print("[+] All cases resampled and verified successfully.")

if __name__ == "__main__":
    main()