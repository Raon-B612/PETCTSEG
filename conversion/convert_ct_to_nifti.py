import os
import sys
import shutil
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import psutil
import pandas as pd
from tqdm import tqdm

def get_script_dir():
    try:
        # Standard script execution
        return Path(__file__).resolve().parent
    except NameError:
        # Fallback for interactive or notebook environments
        return Path(sys.argv[0]).resolve().parent

# === Centralized Path Configuration ===
SCRIPT_DIR        = Path(get_script_dir())
CSV_INPUT_PATH    = SCRIPT_DIR.parent / "prep" / "dcm2nii_prep.csv"
ERROR_OUTPUT_PATH = SCRIPT_DIR.parent / "logs" / "dcm2nii_error.csv"

# Allocate half of physical CPU cores with fallback to avoid TypeError on None
NUM_WORKERS = max(1, (psutil.cpu_count(logical=False) or os.cpu_count() or 2) - 1)

def execute_command(task):
    """
    Execute dcm2niix per DICOM directory with mandatory overwrite (-w 1).
    """
    row_idx, raw_subjid, raw_outdir, raw_ctdcm = task

    # 1. Validate raw inputs before string casting
    if (
        pd.isna(raw_subjid) or pd.isna(raw_outdir) or pd.isna(raw_ctdcm) or
        str(raw_subjid).strip() == '' or str(raw_outdir).strip() == '' or str(raw_ctdcm).strip() == ''
    ):
        return row_idx, "Missing or invalid values (subjid, output_dir, or ctdcm) in CSV row"

    subjid = str(raw_subjid).strip()
    outdir = Path(str(raw_outdir).strip())
    ctpath = Path(str(raw_ctdcm).strip())

    # 2. Check source DICOM directory existence
    if not ctpath.exists() or not ctpath.is_dir():
        return row_idx, f"Source DICOM directory not found: {ctpath}"

    # Ensure target output directory exists
    outdir.mkdir(parents=True, exist_ok=True)

    # 3. Construct command as argument list (safe against shell injection and spacing bugs)
    # -o: output directory, -f: filename prefix, -z y: gzip compression, -w 1: overwrite existing
    cmd = [
        "dcm2niix",
        "-o", str(outdir),
        "-f", f"CT_{subjid}",
        "-z", "y",
        "-w", "1",
        str(ctpath)
    ]

    try:
        # Execute external CLI process and capture stdout/stderr streams
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            shell=False
        )

        # Check process exit code
        if result.returncode != 0:
            err_msg = result.stderr.strip() or result.stdout.strip()
            clean_msg = " ".join(err_msg.split())
            return row_idx, f"dcm2niix failed (code {result.returncode}): {clean_msg[:300]}"

        # 4. Post-execution output integrity check
        expected_output = outdir / f"CT_{subjid}.nii.gz"
        if not expected_output.is_file() or expected_output.stat().st_size == 0:
            # Fallback check if dcm2niix added series/echo suffix
            generated_files = list(outdir.glob(f"CT_{subjid}*.nii.gz"))
            valid_generated = [f for f in generated_files if f.stat().st_size > 0]
            if not valid_generated:
                return row_idx, f"Output file missing or 0 bytes after conversion in {outdir}"

        return row_idx, None

    except Exception as e:
        return row_idx, f"Execution exception: {str(e)}"

def main():
    if not CSV_INPUT_PATH.exists():
        print(f"[ERROR] Input CSV not found: {CSV_INPUT_PATH}")
        sys.exit(1)

    # Verify dcm2niix executable exists in system PATH
    if not shutil.which("dcm2niix"):
        print("[ERROR] 'dcm2niix' executable not found in system PATH.")
        sys.exit(1)

    # Prepare log output directory
    ERROR_OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(CSV_INPUT_PATH).reset_index(drop=True)

    # Extract primitive values directly via fast C-level zip to bypass pandas overhead
    tasks = list(zip(df.index, df['subjid'], df['output_dir'], df['ctdcm']))

    print(f"Total entries: {len(df)}")
    print(f"Using {NUM_WORKERS} workers for dcm2niix conversion (overwrite mode enabled).")

    progress_bar = tqdm(total=len(tasks), desc='Converting DICOM to NIfTI', unit='scan')
    error_records = []

    # Using ThreadPoolExecutor as GIL is released during subprocess execution
    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
        futures = {executor.submit(execute_command, task): task[0] for task in tasks}

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
                err_row['error_msg'] = str(e)
                error_records.append(err_row)
            finally:
                progress_bar.update(1)

    progress_bar.close()

    # Save failed cases to error CSV
    if error_records:
        error_df = pd.DataFrame(error_records)
        error_df.to_csv(ERROR_OUTPUT_PATH, index=False)
        print(f"[!] Saved {len(error_records)} failed cases to {ERROR_OUTPUT_PATH}")
    else:
        print("[+] All DICOM datasets converted and verified successfully.")

if __name__ == '__main__':
    main()