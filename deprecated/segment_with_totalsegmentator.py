import csv
import time
import gc
import re
import os
import sys
import socket
import requests
from tqdm import tqdm
from pathlib import Path
from totalsegmentator.python_api import totalsegmentator

def get_local_ip():
    try:
        # Dummy connection to get the real local IP interface
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
        s.close()
        return ip
    except Exception:
        return "127.0.0.1"

# Automatically sets the LHM URL based on current machine IP
LOCAL_IP = get_local_ip()
LHM_DATA_URL = f"http://{LOCAL_IP}:8085/data.json"
TEMP_THRESHOLD = 55.0
TEMP_CHECK_INTERVAL = 2

def get_cpu_temp_from_lhm(url=LHM_DATA_URL):
    r = requests.get(url, timeout=5)
    r.raise_for_status()
    data = r.json()
    core_temps = []

    def extract_temp(value):
        s = str(value).strip()
        m = re.match(r"(\d{2})(?:\.\d+)?", s)
        if m:
            return float(m.group(1))
        return None

    def walk(node, parents=None):
        if parents is None:
            parents = []
        if isinstance(node, dict):
            text = str(node.get("Text", ""))
            value = str(node.get("Value", ""))
            children = node.get("Children", [])
            current_path = parents + ([text] if text else [])
            path_str = " > ".join(current_path)

            if "Intel Core i9-14900K" in path_str and "Temperatures" in path_str:
                if text.startswith("P-Core #") or text.startswith("E-Core #"):
                    if "Distance to TjMax" not in text:
                        temp = extract_temp(value)
                        if temp is not None:
                            core_temps.append(temp)
            for child in children:
                walk(child, current_path)
        elif isinstance(node, list):
            for item in node:
                walk(item, parents)

    walk(data)
    if not core_temps:
        raise RuntimeError("CPU core temps not found in LHM data")
    return max(core_temps)

def wait_until_cpu_cool(threshold=TEMP_THRESHOLD, interval=TEMP_CHECK_INTERVAL):
    print()
    while True:
        temp = get_cpu_temp_from_lhm()
        if temp < threshold:
            # print(f"...CPU temp OK: {temp:.1f} < {threshold}.")
            # print(f" ", end="")
            return
        # print(f"...CPU temp high: {temp:.1f} >= {threshold}, waiting {interval}s.")
        print(f".", end="")
        time.sleep(interval)

def run_total_segmentator(organ, infile, outfile):
    if os.path.exists(outfile):
        if os.path.getsize(outfile) == 0:
            print(f"\n...REMOVE: zero-size file -> {outfile}")
            os.remove(outfile)
        else:
            print(f"\n...SKIP: Output exists -> {outfile}")
            return
    
    wait_until_cpu_cool()

    # # Determine if higher order resampling is needed (corresponds to --higher_order_resampling)
    # use_higher_order = True if "HCV" in infile else False
    use_higher_order = False
    
    # infile_filename = os.path.basename(infile)
    print(f"Running task: {organ} on {infile}", end="")

    # Replaced subprocess with Python API call
    # All available parameters for totalsegmentator v2.x are listed below
    totalsegmentator(
        input=infile,
        output=outfile,
        ml=True,
        nr_thr_resamp=1,
        nr_thr_saving=1,
        fast=False,
        nora_tag="None",
        preview=False,
        task=organ,
        roi_subset=None,
        statistics=False,
        quiet=True,
        verbose=False,
        test=0,
        skip_saving=False,
        device="gpu",
        license_number=None,
        higher_order_resampling=use_higher_order, 
        # body_seg=False,
        # force_split=True,
        output_type="nifti",
        fastest=False,
        radiomics=False,
        crop_path=None,                 
        stats_aggregation="mean"
    )

    if not (os.path.exists(outfile) and os.path.getsize(outfile) > 0):
        print(f"...MISSING OUTPUT: {outfile}\n")

def get_script_dir():
    try:
        # Standard script execution
        return Path(__file__).resolve().parent
    except NameError:
        # Fallback for environments without __file__ (e.g., interactive, Jupyter)
        return Path(sys.argv[0]).resolve().parent

def main():
    # Load the CSV file into a pandas DataFrame
    SCRIPT_DIR = Path(get_script_dir())
    csv_file_path = SCRIPT_DIR.parent / "prep" / "tseg_prep_emma.csv"
    error_rows = []

    with open(csv_file_path, newline="", encoding="utf-8-sig") as csv_file:
        csv_reader = list(csv.DictReader(csv_file))
        for row in tqdm(csv_reader, desc="Processing", unit="file"):
            ct_path = row["ct"]
            seg_path = row["seg"]
            organ = re.sub(r"^tseg_", "", row["module"])
            
            try:
                run_total_segmentator(organ, ct_path, seg_path)
                # tqdm.write(f'|> Processed: {row["ct"]}\n')
                # tqdm.write(f'|> Saved to {row["seg"]}\n')
            except Exception as e:
                tqdm.write(f'!= Error processing {row["ct"]}: {e}\n')
                error_rows.append(row)
            finally:
                gc.collect()

    print("\nProcessing complete.")
    if error_rows:
        print("\nFailed rows:")
        for row in error_rows:
            print(row)

if __name__ == "__main__":
    main()