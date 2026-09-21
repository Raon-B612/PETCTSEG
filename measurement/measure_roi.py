import os
import sys
import re
import gc
import argparse
import multiprocessing as mp
from pathlib import Path
from datetime import date
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import psutil
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

# === Centralized Path Configuration ===
SCRIPT_DIR               = Path(get_script_dir())
EXCEL_INPUT_PATH         = SCRIPT_DIR.parent / "prep" / "measure_prep.csv"
LOGS_DIR                 = SCRIPT_DIR.parent / "logs"
ERROR_LOG_STATIC_PATH    = LOGS_DIR / "measurement_error.csv"
ERROR_LOG_DYNAMIC_PATH   = LOGS_DIR / f"measurement_error_{date.today().isoformat()}.csv"

# === Exception Formatting Helper ===
def format_exception_msg(err: Exception) -> str:
    """
    Extract core actionable description from SimpleITK/ITK C++ exceptions
    and sanitize multi-line traceback into a single-line string for CSV logging.
    """
    msg = str(err)
    if "Exception thrown in SimpleITK" in msg or "itk::ERROR" in msg.lower():
        match = re.search(r"itk::ERROR:\s*(?:[A-Za-z0-9_]+\([^)]*\):\s*)?([^\n\r]+)", msg, re.IGNORECASE)
        if match:
            return f"SimpleITK Error: {match.group(1).strip()}"
    return " ".join(msg.split())

# === Worker IPC Counter State ===
_worker_active_counter = None

def _init_worker(counter):
    """Initializer to share the atomic counter with each child worker process."""
    global _worker_active_counter
    _worker_active_counter = counter

# === Group Worker Function ===
def process_pet_ct_group(group_data):
    """
    Processes all segmentations associated with a single PET/CT pair.
    Returns (success_flag, list_of_error_records) to eliminate multiprocessing file contention.
    """
    # 워커 시작 시 활성 프로세스 수 카운트 증가
    if _worker_active_counter is not None:
        with _worker_active_counter.get_lock():
            _worker_active_counter.value += 1

    ct_path, pet_path, tasks_in_group = group_data
    errors_encountered = []
    ct_img = None
    pet_img = None

    # Limit SimpleITK internal threads to 1 per worker to prevent CPU contention
    sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(1)

    try:
        # 1. Skip entire group if all expected output CSV files already exist and are non-empty
        all_completed = all(
            os.path.isfile(item['csv']) and os.path.getsize(item['csv']) > 0
            for item in tasks_in_group
        )
        if all_completed:
            return True, []

        # 2. Check source file existence for CT and PET
        if not os.path.isfile(ct_path):
            err_msg = f"CT file not found: {ct_path}"
            return False, [{**item, 'error_msg': err_msg} for item in tasks_in_group]

        if not os.path.isfile(pet_path):
            err_msg = f"PET file not found: {pet_path}"
            return False, [{**item, 'error_msg': err_msg} for item in tasks_in_group]

        # Load large 3D volumes once per (CT, PET) group
        ct_img = sitk.ReadImage(ct_path)
        pet_img = sitk.ReadImage(pet_path)

        # Handle 4D singleton axis if present in PET volume
        if pet_img.GetDimension() == 4 and pet_img.GetSize()[3] == 1:
            pet_img = pet_img[:, :, :, 0]
        elif pet_img.GetDimension() != 3:
            err_msg = f"Invalid PET dimension: {pet_img.GetDimension()}D (expected 3D)"
            return False, [{**item, 'error_msg': err_msg} for item in tasks_in_group]

        if ct_img.GetDimension() != 3:
            err_msg = f"Invalid CT dimension: {ct_img.GetDimension()}D (expected 3D)"
            return False, [{**item, 'error_msg': err_msg} for item in tasks_in_group]

        # Verify geometric dimension alignment between CT and PET
        if ct_img.GetSize() != pet_img.GetSize():
            err_msg = f"Grid size mismatch between CT {ct_img.GetSize()} and PET {pet_img.GetSize()}"
            return False, [{**item, 'error_msg': err_msg} for item in tasks_in_group]

        # Group tasks by unique segmentation file
        seg_dict = {}
        for item in tasks_in_group:
            seg_dict.setdefault(item['seg'], []).append(item)

        # Output schema definitions (maintaining both original and normalized metric columns)
        metric_columns = [
            'label', 'suv_mean', 'suv_sd', 'suv_median', 'suv_max', 'suv_min',
            'hu_mean', 'hu_sd', 'hu_median', 'hu_max', 'hu_min',
            'volume', 'volume_ml', 'area', 'mean_slice_area_mm2', 'tlg', 'tlg_raw'
        ]

        for seg_path, sub_items in seg_dict.items():
            # Skip this SEG if all associated target CSVs already exist
            if all(os.path.isfile(r['csv']) and os.path.getsize(r['csv']) > 0 for r in sub_items):
                continue

            if not os.path.isfile(seg_path):
                for r in sub_items:
                    errors_encountered.append({**r, 'error_msg': f"SEG file not found: {seg_path}"})
                continue

            seg_img = None
            ct_filter = None
            pet_filter = None
            shape_filter = None

            try:
                # Load and cast mask to UInt32 to prevent 'Pixel type not supported' exceptions
                raw_seg = sitk.ReadImage(seg_path)
                if raw_seg.GetDimension() == 4 and raw_seg.GetSize()[3] == 1:
                    raw_seg = raw_seg[:, :, :, 0]
                seg_img = sitk.Cast(raw_seg, sitk.sitkUInt32)
                del raw_seg

                # Verify spatial alignment between SEG and CT
                if seg_img.GetSize() != ct_img.GetSize():
                    err_msg = f"Grid size mismatch between SEG {seg_img.GetSize()} and CT {ct_img.GetSize()}"
                    for r in sub_items:
                        errors_encountered.append({**r, 'error_msg': err_msg})
                    continue

                # Initialize and execute SimpleITK feature extractors
                ct_filter = sitk.LabelIntensityStatisticsImageFilter()
                pet_filter = sitk.LabelIntensityStatisticsImageFilter()
                shape_filter = sitk.LabelShapeStatisticsImageFilter()

                ct_filter.Execute(seg_img, ct_img)
                pet_filter.Execute(seg_img, pet_img)
                shape_filter.Execute(seg_img)

                labels = shape_filter.GetLabels()
                z_spacing = seg_img.GetSpacing()[2]

                results_cache = []
                for label in labels:
                    bbox = shape_filter.GetBoundingBox(label)
                    z_len = bbox[5] * z_spacing
                    vol_mm3 = shape_filter.GetPhysicalSize(label)
                    vol_ml = vol_mm3 / 1000.0  # Normalized: 1 cm^3 (mL) = 1000 mm^3
                    
                    suv_mean = pet_filter.GetMean(label)

                    results_cache.append({
                        'label': label,
                        'suv_mean': suv_mean,
                        'suv_sd': pet_filter.GetStandardDeviation(label),
                        'suv_median': pet_filter.GetMedian(label),
                        'suv_max': pet_filter.GetMaximum(label),
                        'suv_min': pet_filter.GetMinimum(label),
                        'hu_mean': ct_filter.GetMean(label),
                        'hu_sd': ct_filter.GetStandardDeviation(label),
                        'hu_median': ct_filter.GetMedian(label),
                        'hu_max': ct_filter.GetMaximum(label),
                        'hu_min': ct_filter.GetMinimum(label),
                        'volume': vol_mm3,                                          # Preserved mm^3
                        'volume_ml': vol_ml,                                        # Standard mL
                        'area': (vol_mm3 / z_len * 10.0) if z_len > 0 else 0.0,      # Legacy formula
                        'mean_slice_area_mm2': (vol_mm3 / z_len) if z_len > 0 else 0.0,
                        'tlg': suv_mean * vol_ml,                                    # Standard TLG: SUV_mean * MTV (mL)
                        'tlg_raw': suv_mean * vol_mm3                                # Legacy raw TLG
                    })

                # Handle empty segmentation masks safely without KeyError
                if results_cache:
                    df_out = pd.DataFrame(results_cache).set_index('label')
                else:
                    df_out = pd.DataFrame(columns=metric_columns).set_index('label')

                # Atomic write to disk for each row using process-unique temporary files
                for r in sub_items:
                    out_csv = Path(r['csv'])
                    out_csv.parent.mkdir(parents=True, exist_ok=True)
                    tmp_csv = out_csv.parent / f".tmp_{os.getpid()}_{out_csv.name}"
                    
                    try:
                        df_out.to_csv(tmp_csv)
                        os.replace(tmp_csv, out_csv)
                    finally:
                        if tmp_csv.exists():
                            try:
                                tmp_csv.unlink()
                            except OSError:
                                pass

            except Exception as seg_err:
                for r in sub_items:
                    errors_encountered.append({**r, 'error_msg': format_exception_msg(seg_err)})
            finally:
                del seg_img, ct_filter, pet_filter, shape_filter

        return (len(errors_encountered) == 0), errors_encountered

    except Exception as e:
        err_str = format_exception_msg(e)
        return False, [{**item, 'error_msg': err_str} for item in tasks_in_group]

    finally:
        del ct_img, pet_img
        gc.collect()
        # 조기 종료(return), 성공, 예외 상황 상관없이 항상 카운트 감소
        if _worker_active_counter is not None:
            with _worker_active_counter.get_lock():
                _worker_active_counter.value -= 1

# === Main Pipeline Execution ===
def main():
    parser = argparse.ArgumentParser(description="Extract PET/CT quantitative metrics from segmentations.")
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

    if not EXCEL_INPUT_PATH.exists():
        print(f"[ERROR] Input file not found: {EXCEL_INPUT_PATH}")
        sys.exit(1)

    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    cols = ['ct', 'pet', 'seg', 'module', 'subjid', 'output_dir', 'csv']
    df = pd.read_csv(EXCEL_INPUT_PATH, usecols=cols, dtype={'subjid': 'str'})

    # Pre-validate input rows for NaN or blank values prior to multi-process dispatch
    valid_tasks = []
    immediate_errors = []
    required_keys = ['ct', 'pet', 'seg', 'output_dir', 'csv']

    for row in df.to_dict('records'):
        is_invalid = any(pd.isna(row.get(k)) or str(row.get(k)).strip() == '' for k in required_keys)
        if is_invalid:
            immediate_errors.append({**row, 'error_msg': 'Missing or empty path in input CSV row'})
        else:
            valid_tasks.append(row)

    # Group valid tasks by unique (CT, PET) file pairs
    grouped_tasks = {}
    for row in valid_tasks:
        key = (str(row['ct']).strip(), str(row['pet']).strip())
        grouped_tasks.setdefault(key, []).append(row)

    task_payloads = [
        (ct, pet, items)
        for (ct, pet), items in grouped_tasks.items()
    ]

    print(f"Total entries: {len(df)} (Valid: {len(valid_tasks)}, Invalid: {len(immediate_errors)})")
    print(f"Unique PET/CT pairs: {len(task_payloads)}")
    print(f"Active Mode: {storage_mode} -> Running with {num_workers} workers")

    all_error_records = list(immediate_errors)

    # 프로세스 간 공유할 스레드/프로세스 안전한 카운터
    active_counter = mp.Value('i', 0)

    with ProcessPoolExecutor(
        max_workers=num_workers,
        initializer=_init_worker,
        initargs=(active_counter,)
    ) as executor:
        futures = [executor.submit(process_pet_ct_group, payload) for payload in task_payloads]
        not_done = set(futures)

        with tqdm(total=len(task_payloads), ascii=' #', desc="Processing") as pbar:
            while not_done:
                # 2초 동안 대기하면서 완료된 작업 수거 (작업이 없어도 2초마다 반환되어 루프 돎)
                done, not_done = wait(not_done, timeout=2, return_when=FIRST_COMPLETED)

                for future in done:
                    try:
                        _, errors = future.result()
                        if errors:
                            all_error_records.extend(errors)
                    except Exception as fut_err:
                        all_error_records.append({'error_msg': format_exception_msg(fut_err)})
                    pbar.update(1)

                # 현재 작업 중인 실제 워커 수 조회 및 tqdm postfix 갱신
                with active_counter.get_lock():
                    current_active = active_counter.value

                pbar.set_postfix(active=f"{current_active}/{num_workers}")

    # Centralized single-thread logging with atomic write and append fallback
    if all_error_records:
        error_df = pd.DataFrame(all_error_records)
        for log_path in [ERROR_LOG_STATIC_PATH, ERROR_LOG_DYNAMIC_PATH]:
            log_file = Path(log_path)
            if log_file.exists():
                error_df.to_csv(log_file, mode='a', index=False, header=False)
            else:
                error_df.to_csv(log_file, index=False)
        print(f"[!] Processing finished with {len(all_error_records)} errors. Log saved to {ERROR_LOG_STATIC_PATH}")
    else:
        print("[+] All measurements completed and verified successfully.")

if __name__ == "__main__":
    main()
    