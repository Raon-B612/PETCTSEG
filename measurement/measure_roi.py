import os
import SimpleITK as sitk
import pandas as pd
import numpy as np
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import date
import psutil
from pathlib import Path

# === Centralized Path Configuration ===
BASE_DIR                     = Path("~/Documents/PETCTSEG/scripts").expanduser()
INPUT_FILENAME               = "prep/measure_prep.csv"
ERROR_LOG_STATIC_FILENAME    = "logs/measurement_error.csv"
ERROR_LOG_DYNAMIC_PREFIX     = "logs/measurement_error_"

# HDD parallelism optimization: Keep low to avoid disk thrashing
# Capping workers to ensure sequential-like I/O for large volumes
NUM_WORKERS                  = max(4, max(1, psutil.cpu_count(logical=False) // 2))

EXCEL_INPUT_PATH             = os.path.join(BASE_DIR, INPUT_FILENAME)
ERROR_LOG_STATIC_PATH        = os.path.join(BASE_DIR, ERROR_LOG_STATIC_FILENAME)
ERROR_LOG_DYNAMIC_PATH       = os.path.join(BASE_DIR, f"{ERROR_LOG_DYNAMIC_PREFIX}{date.today().isoformat()}.csv")

def save_error_row(row_dict):
    error_df = pd.DataFrame([row_dict])
    for path in [ERROR_LOG_STATIC_PATH, ERROR_LOG_DYNAMIC_PATH]:
        try:
            if os.path.exists(path):
                error_df.to_csv(path, mode='a', index=False, header=False)
            else:
                error_df.to_csv(path, index=False)
        except Exception as log_err:
            print(f"Failed to write error log to {path}: {log_err}")

def process_pet_ct_group(group_data):
    """
    Processes all segmentations associated with a single PET/CT pair.
    group_data: ( (ct_path, pet_path), dataframe_subset )
    """
    (ct_path, pet_path), subset_df = group_data
    sitk.ProcessObject_SetGlobalDefaultNumberOfThreads(1)

    # Pre-check: if all output files in this entire subset exist, skip loading PET/CT
    if all(os.path.exists(r['csv']) and os.path.getsize(r['csv']) > 0 for _, r in subset_df.iterrows()):
        return None

    try:
        # Load the large volumes only once per group
        ct_img = sitk.ReadImage(ct_path)
        pet_img = sitk.ReadImage(pet_path)

        # Process each unique SEG file within this PET/CT group
        for seg_path, seg_group in subset_df.groupby('seg'):
            
            # Skip if this specific SEG result already exists
            if all(os.path.exists(r['csv']) and os.path.getsize(r['csv']) > 0 for _, r in seg_group.iterrows()):
                continue

            seg_img = sitk.ReadImage(seg_path)

            # Initialize and execute filters
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
                vol = shape_filter.GetPhysicalSize(label)
                area = vol / z_len * 10 if z_len > 0 else 0
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
                    'volume': vol,
                    'area': area,
                    'tlg': suv_mean * vol
                })

            df_out = pd.DataFrame(results_cache).set_index('label')

            # Write out results for each row associated with this SEG
            for _, row in seg_group.iterrows():
                os.makedirs(row['output_dir'], exist_ok=True)
                df_out.to_csv(row['csv'])

    except Exception as e:
        first_row = subset_df.iloc[0]
        print(f"Error processing Subject {first_row['subjid']}: {e}")
        for _, row in subset_df.iterrows():
            save_error_row(row.to_dict())

def main():
    cols = ['ct', 'pet', 'seg', 'module', 'subjid', 'output_dir', 'csv']
    df = pd.read_csv(EXCEL_INPUT_PATH, usecols=cols, dtype={'subjid': 'str'})

    # Group by CT and PET to ensure large volumes are loaded once
    pet_ct_groups = list(df.groupby(['ct', 'pet']))
    
    print(f"Total entries: {len(df)}")
    print(f"Unique PET/CT pairs: {len(pet_ct_groups)}")
    print(f"Using {NUM_WORKERS} worker processes for optimized I/O")

    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        futures = [executor.submit(process_pet_ct_group, group) for group in pet_ct_groups]
        for _ in tqdm(as_completed(futures), total=len(pet_ct_groups), ascii=' #'):
            pass

if __name__ == "__main__":
    main()