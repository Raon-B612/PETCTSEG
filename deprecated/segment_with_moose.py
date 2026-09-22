import re
import os
import time
import pandas as pd
import requests
from moosez import moose
from tqdm import tqdm
import socket

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
    """
    Fetch CPU core temperatures from LibreHardwareMonitor JSON API
    and return the maximum core temperature.
    """
    r = requests.get(url, timeout=5)
    r.raise_for_status()
    data = r.json()

    core_temps = []

    def extract_temp(value):
        """
        Extract numeric temperature value from a string.
        Examples:
            '45.0 °C' -> 45.0
            '9 °C'    -> 9.0
            '100 °C'  -> 100.0
        """
        if value is None:
            return None

        s = str(value).strip()
        m = re.search(r"(-?\d+(?:\.\d+)?)", s)
        if not m:
            return None

        try:
            temp = float(m.group(1))
        except ValueError:
            return None

        # Basic sanity check for temperature values
        if temp < 0 or temp > 150:
            return None

        return temp

    def walk(node, parents=None):
        """
        Recursively traverse JSON tree to find CPU core temperatures.
        """
        if parents is None:
            parents = []

        if isinstance(node, dict):
            text = str(node.get("Text", "")).strip()
            value = node.get("Value", "")
            children = node.get("Children", [])

            current_path = parents + ([text] if text else [])
            path_str = " > ".join(current_path)

            # Limit to the target CPU and temperature section
            if "Intel Core i9-14900K" in path_str and "Temperatures" in path_str:
                # Accept only per-core temperature sensors
                if re.match(r"^[PE]-Core #\d+$", text):
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
        raise RuntimeError("CPU core temperatures were not found in LibreHardwareMonitor data.json")

    return max(core_temps)


def wait_until_cpu_cool(threshold=TEMP_THRESHOLD, interval=TEMP_CHECK_INTERVAL):
    """
    Block execution until CPU temperature drops below threshold.
    """
    while True:
        temp = get_cpu_temp_from_lhm()
        if temp < threshold:
            # print(f"...CPU temp OK: {temp:.1f} < {threshold}")
            print(f" ", end="")
            return
        # print(f"...CPU temp high: {temp:.1f} >= {threshold}, waiting {interval}s")
        print(f".", end="")
        time.sleep(interval)


def validate_row_inputs(row):
    """
    Validate required row fields before processing.
    """
    required_cols = ["ct", "module", "output_dir", "outfile"]
    missing_cols = [col for col in required_cols if col not in row.index]
    if missing_cols:
        raise KeyError(f"Missing required columns: {missing_cols}")

    infile = row["ct"]
    module_raw = row["module"]
    output_path = row["output_dir"]
    outfile = row["outfile"]

    if pd.isna(infile) or str(infile).strip() == "":
        raise ValueError("Empty input CT path")
    if pd.isna(module_raw) or str(module_raw).strip() == "":
        raise ValueError("Empty module value")
    if pd.isna(output_path) or str(output_path).strip() == "":
        raise ValueError("Empty output_dir value")
    if pd.isna(outfile) or str(outfile).strip() == "":
        raise ValueError("Empty outfile value")

    infile = str(infile).strip()
    module_raw = str(module_raw).strip()
    output_path = str(output_path).strip()
    outfile = str(outfile).strip()

    if not os.path.exists(infile):
        raise FileNotFoundError(f"Input CT not found: {infile}")
    if os.path.isdir(infile):
        raise IsADirectoryError(f"Input CT path is a directory, not a file: {infile}")

    return infile, module_raw, output_path, outfile


def validate_output_file(outfile):
    """
    Validate generated output file.
    """
    if not os.path.exists(outfile):
        raise FileNotFoundError(f"Output file not found: {outfile}")

    if os.path.isdir(outfile):
        raise IsADirectoryError(f"Output path is a directory, not a file: {outfile}")

    size = os.path.getsize(outfile)
    if size == 0:
        raise RuntimeError(f"Output file is zero-size: {outfile}")

    return size


def process_row(row):
    """
    Process a single row from the CSV:
    - Validate input values
    - Check existing output
    - Wait for CPU cooling
    - Run moose segmentation
    - Validate output file
    """
    try:
        infile, module_raw, output_path, outfile = validate_row_inputs(row)
        
        # Normalize module name by removing 'moose_' prefix if present
        module = f"clin_ct_{re.sub(r'^moose_', '', module_raw)}"

        # Ensure output directory exists
        os.makedirs(output_path, exist_ok=True)

        # If output already exists, remove only if it is zero-size
        if os.path.exists(outfile):
            if os.path.getsize(outfile) == 0:
                print(f"REMOVE: zero-size file -> {outfile}")
                os.remove(outfile)
            else:
                print(f"SKIP: Output exists -> {outfile}")
                return

        print()
        # wait_until_cpu_cool()

        
        # infile_filename = os.path.basename(infile)
        print(f"|> Running Moose API for task: {module} on {infile}")
    
        # Run segmentation in CUDA mode
        moose(infile, module, output_path, "cuda")

        # Validate output after execution
        size = validate_output_file(outfile)
        print(f"|> Output verified ({size} bytes) -> {outfile}\n")

    except Exception as e:
        print(f"!= ERROR: Failed processing ct={row.get('ct', 'NA')}\n")
        # print(f"!= OUTFILE: {row.get('outfile', 'NA')}")
        # print(f"!= Reason: {e}")
        return


def main():
    """
    Main execution:
    - Load CSV
    - Iterate rows sequentially
    - Process each case with isolated error handling
    """
    df = pd.read_csv(
        "C:/Users/Heartwork101/Documents/PETCTSEG/scripts/prep/moose_prep_emma.csv"
    )

    for _, row in tqdm(df.iterrows(), total=len(df), desc="Processing rows"):
        try:
            process_row(row)
        except Exception as e:
            print(f"FATAL ROW ERROR: {row.get('ct', 'NA')} -> {e}")
            continue

    print("Processing complete.")


if __name__ == "__main__":
    main()
    