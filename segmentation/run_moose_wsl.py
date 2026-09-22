import os
from pathlib import Path
import re
import socket
import time
from moosez import moose
import pandas as pd
import requests
import torch
from tqdm import tqdm

LHM_DATA_URL = "http://192.168.0.3:8085/data.json"
TEMP_THRESHOLD = 55.0
TEMP_CHECK_INTERVAL = 2


def get_cpu_temp_from_lhm(url=LHM_DATA_URL):
    """LHM에서 CPU 온도를 조회하며, 일시적 네트워크 지연 시 최대 3회 재시도합니다."""
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=3)
            r.raise_for_status()
            data = r.json()
            core_temps = []

            def extract_temp(value):
                s = str(value).strip()
                m = re.match(r"(\d{2})(?:\.\d+)?", s)
                return float(m.group(1)) if m else None

            def walk(node, parents=None):
                if parents is None:
                    parents = []
                if isinstance(node, dict):
                    text = str(node.get("Text", ""))
                    value = str(node.get("Value", ""))
                    children = node.get("Children", [])
                    current_path = parents + ([text] if text else [])
                    path_str = " > ".join(current_path)

                    if (
                        "Intel Core i9-14900K" in path_str
                        and "Temperatures" in path_str
                    ):
                        if (
                            text.startswith("P-Core #")
                            or text.startswith("E-Core #")
                        ) and "Distance to TjMax" not in text:
                            temp = extract_temp(value)
                            if temp is not None:
                                core_temps.append(temp)
                    for child in children:
                        walk(child, current_path)
                elif isinstance(node, list):
                    for item in node:
                        walk(item, parents)

            walk(data)
            if core_temps:
                return max(core_temps)
        except Exception:
            time.sleep(1)

    return None  # 모니터링 서버 응답 없을 시 None 반환하여 파이프라인 중단 방지


def wait_until_cpu_cool(
    threshold=TEMP_THRESHOLD, interval=TEMP_CHECK_INTERVAL
):
    """CPU 온도가 기준치 이하로 떨어질 때까지 대기합니다."""
    waiting_announced = False
    while True:
        temp = get_cpu_temp_from_lhm()
        if temp is None:
            # LHM 서버 미응답 시 과도한 대기를 방지하고 작업 진행
            break
        if temp < threshold:
            break
        if not waiting_announced:
            tqdm.write(
                f"[WAIT]  CPU Temp high ({temp:.1f}°C >= {threshold}°C). Cooling down..."
            )
            waiting_announced = True
        time.sleep(interval)


def check_file_exists_case_insensitive(filepath):
    """WSL/Linux 환경에서 대소문자 불일치(clin_CT vs clin_ct)로 인한 파일 탐색 실패를 방어합니다."""
    if os.path.exists(filepath):
        return filepath

    parent_dir = os.path.dirname(filepath)
    if not os.path.exists(parent_dir):
        return None

    target_name_lower = os.path.basename(filepath).lower()
    for existing_file in os.listdir(parent_dir):
        if existing_file.lower() == target_name_lower:
            return os.path.join(parent_dir, existing_file)

    return None


def validate_row_inputs(row):
    """입력 데이터 무결성 검증"""
    required_cols = ["ct", "module", "output_dir", "outfile", "seg"]
    missing_cols = [col for col in required_cols if col not in row.index]
    if missing_cols:
        raise KeyError(f"Missing required columns: {missing_cols}")

    infile = str(row["ct"]).strip()
    module_raw = str(row["module"]).strip()
    output_path = str(row["output_dir"]).strip()
    outfile = str(row["outfile"]).strip()
    segfile = str(row["seg"]).strip()

    if not infile or pd.isna(row["ct"]):
        raise ValueError("Empty input CT path")
    if not module_raw or pd.isna(row["module"]):
        raise ValueError("Empty module value")
    if not output_path or pd.isna(row["output_dir"]):
        raise ValueError("Empty output_dir value")
    # if not outfile or pd.isna(row["outfile"]):
    #     raise ValueError("Empty outfile value")
    if not segfile or pd.isna(row["seg"]):
            raise ValueError("Empty seg value")

    if not os.path.exists(infile):
        raise FileNotFoundError(f"Input CT not found: {infile}")
    if os.path.isdir(infile):
        raise IsADirectoryError(f"Input CT path is a directory: {infile}")

    return infile, module_raw, output_path, outfile, segfile


def process_row(row):
    try:
        infile, module_raw, output_path, outfile, segfile = validate_row_inputs(row)

        module = f"clin_ct_{re.sub(r'^moose_', '', module_raw)}"
        os.makedirs(output_path, exist_ok=True)

        # 대소문자 무관 파일 존재 확인
        existing_segfile = check_file_exists_case_insensitive(segfile)
        # existing_outfile = check_file_exists_case_insensitive(outfile)

        if existing_segfile:
            if os.path.getsize(existing_segfile) == 0:
                tqdm.write(
                    f"[CLEAN] Zero-size file removed -> {os.path.basename(existing_segfile)}"
                )
                os.remove(existing_segfile)
            else:
                tqdm.write(
                    f"[SKIP]  Output exists -> {os.path.basename(existing_segfile)}"
                )
                return

        # CPU 쿨다운 대기
        wait_until_cpu_cool()

        # 현재 행의 infile 경로에서 환자_날짜 폴더명 동적 추출
        subjid_long = Path(infile).parent.name
        tqdm.write(
            f"[RUN]   Moose task: {module} | Target: {subjid_long}"
        )

        # Segmentation 실행 (CUDA)
        moose(infile, module, output_path, "cuda")

        # 결과 검증
        final_outfile = check_file_exists_case_insensitive(outfile)
        if not final_outfile or not os.path.exists(final_outfile):
            raise FileNotFoundError(f"Generated output not found: {outfile}")

        size = os.path.getsize(final_outfile)
        if size == 0:
            raise RuntimeError(f"Generated output is zero-size: {outfile}")

        # outfile -> segfile 
        if os.path.abspath(final_outfile) != os.path.abspath(segfile):
            os.replace(final_outfile, segfile)

        tqdm.write(f"[OK]    Verified ({size:,} bytes) -> {os.path.basename(segfile)}\n")

        # PyTorch VRAM 캐시 해제 (연속 실행 시 OOM 방지)
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    except Exception as e:
        tqdm.write(f"[FAIL]  Target: {row.get('ct', 'NA')} | Error: {e}\n")


def main():
    csv_path = "/mnt/c/Users/Heartwork101/Documents/PETCTSEG/scripts/prep/moose_prep_wsl.csv"
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV file not found at: {csv_path}")

    df = pd.read_csv(csv_path)

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Total Progress"):
        try:
            process_row(row)
        except Exception as e:
            tqdm.write(f"[FATAL] Row Error: {row.get('ct', 'NA')} -> {e}")
            continue

    tqdm.write("[DONE]  All processing complete.")


if __name__ == "__main__":
    main()