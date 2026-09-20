import os
import SimpleITK as sitk
import pandas as pd
import numpy as np
import nibabel as nib
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import psutil

from pathlib import Path
import sys

def get_script_dir():
    try:
        # Standard script execution
        return Path(__file__).resolve().parent
    except NameError:
        # Fallback for environments without __file__ (e.g., interactive, Jupyter)
        return Path(sys.argv[0]).resolve().parent

# === 경로 및 설정 ===
# BASE_DIR = "C:/Users/Heartwork101/Documents/PETCTSEG/scripts"
BASE_DIR = get_script_dir()
INPUT_FILENAME = "prep/pet_resample_prep.csv"
ERROR_FILENAME = "logs/resampling_error.csv"
NUM_WORKERS = max(1, psutil.cpu_count(logical=False) // 2)
# NUM_WORKERS = 1

CSV_INPUT_PATH = os.path.join(BASE_DIR, INPUT_FILENAME)
ERROR_OUTPUT_PATH = os.path.join(BASE_DIR, ERROR_FILENAME)

# === Float32 → Float16 변환 ===
def convert_float32_to_float16(input_filepath, output_filepath):
    try:
        nifti_img = nib.load(input_filepath, mmap=False)
        data = nifti_img.get_fdata()
        data_float16 = data.astype(np.float16)
        nifti_img_float16 = nib.Nifti1Image(data_float16, affine=nifti_img.affine, header=nifti_img.header)
        nib.save(nifti_img_float16, output_filepath)
    except Exception as e:
        return f"Error converting {input_filepath} to float16: {e}"

# === Rigid 정합 기반 Resample 함수 ===
def resample_pet_image_with_rigid(pet_image, ct_image):
    initial_transform = sitk.CenteredTransformInitializer(
        ct_image,
        pet_image,
        sitk.Euler3DTransform(),  # rigid transform
        sitk.CenteredTransformInitializerFilter.GEOMETRY
    )

    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(ct_image)
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetTransform(initial_transform)
    return resampler.Execute(pet_image)

# === 한 줄 처리 함수 ===
def process_row(row):
    try:
        pet_image_path = row['suv']
        ct_image_path = row['ct']
        output_path = row['pet']

        if pd.isna(pet_image_path) or pd.isna(ct_image_path):
            return row.name, "Missing PET or CT image path"

        # 이미지 불러오기
        pet_image = sitk.ReadImage(pet_image_path)
        ct_image = sitk.ReadImage(ct_image_path)

        # Rigid 정합 후 리샘플링
        resampled_pet_image = resample_pet_image_with_rigid(pet_image, ct_image)
        resampled_pet_image = sitk.Cast(resampled_pet_image, sitk.sitkFloat32)

        # 저장 경로 확보
        output_dir = os.path.dirname(output_path)
        os.makedirs(output_dir, exist_ok=True)

        # 저장
        sitk.WriteImage(resampled_pet_image, output_path)

        # Float16 변환 (단, SPECT는 제외)
        if "TMJ" not in output_path:
            error = convert_float32_to_float16(output_path, output_path)
            if error:
                return row.name, error

        # === 파일 존재 여부 및 크기(0 byte) 확인 추가 ===
        # Check if file exists and size is greater than 0
        if not os.path.exists(output_path):
            return row.name, f"File was not created: {output_path}"
        
        if os.path.getsize(output_path) == 0:
            return row.name, f"File is 0 bytes: {output_path}"

        # 정합 결과 검증
        final_pet_image = sitk.ReadImage(output_path)
        checks = [
            ("Size", final_pet_image.GetSize() == ct_image.GetSize()),
            ("Spacing", np.allclose(final_pet_image.GetSpacing(), ct_image.GetSpacing(), atol=1e-5)),
            ("Origin", np.allclose(final_pet_image.GetOrigin(), ct_image.GetOrigin(), atol=1e-5)),
            ("Direction", np.allclose(final_pet_image.GetDirection(), ct_image.GetDirection(), atol=1e-5))
        ]
        for name, passed in checks:
            if not passed:
                return row.name, f"Mismatch in {name} between resampled PET and CT."

        return row.name, None

    except Exception as e:
        return row.name, str(e)

# === 메인 함수 ===
def main():
    print(f"Running in multi-threaded mode (n={NUM_WORKERS}) ...")
    
    df = pd.read_csv(CSV_INPUT_PATH)
    progress_bar = tqdm(total=len(df), desc='Processing', unit='image')
    error_rows = []
    
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        futures = {executor.submit(process_row, row): row for _, row in df.iterrows()}

        for future in as_completed(futures):
            row = futures[future]
            try:
                row_index, error = future.result()
                if error:
                    # Capture error message and store in error_rows
                    error_data = df.iloc[row_index].copy()
                    error_data['error_msg'] = error
                    error_rows.append(error_data)
            except Exception as e:
                error_data = row.copy()
                error_data['error_msg'] = str(e)
                error_rows.append(error_data)
            progress_bar.update(1)

    progress_bar.close()

    if error_rows:
        error_df = pd.DataFrame(error_rows)
        error_df.to_csv(ERROR_OUTPUT_PATH, index=False)
        print(f"Saved rows with errors to {ERROR_OUTPUT_PATH}")

if __name__ == "__main__":
    main()