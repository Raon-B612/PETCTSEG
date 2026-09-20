import os
import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import psutil
from pathlib import Path

# Function to execute command
def execute_command(row):
    subjid, outdir, ctpath = row['subjid'], row['output_dir'], row['ctdcm']
    
    command = f'dcm2niix -o "{outdir}" -f CT_{subjid} -z y -w 1 "{ctpath}"'
    
    # Check if the output directory exists, if not, create it
    os.makedirs(outdir, exist_ok=True)
    
    # Execute the command
    subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


from pathlib import Path
import sys

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
    CSV_INPUT_PATH = SCRIPT_DIR.parent / "prep" / "dcm2nii_prep.csv"
        
    df = pd.read_csv(CSV_INPUT_PATH)
    
    # Calculate max workers as half of the available CPU cores
    num_cores = psutil.cpu_count(logical=False)
    num_workers = max(1, num_cores // 2)
    # num_workers = 1
    
    print(f"Using {num_workers} workers.")
    
    # Using ThreadPoolExecutor to run commands in parallel
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Create a progress bar instance
        progress_bar = tqdm(total=len(df), desc='Processing')
        
        # Submit tasks to the executor
        futures = [executor.submit(execute_command, row) for _, row in df.iterrows()]
        
        # Wait for the futures to complete and update progress bar accordingly
        for future in as_completed(futures):
            progress_bar.update(1)
        
        progress_bar.close()

if __name__ == '__main__':
    main()
