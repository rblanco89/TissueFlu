import sys
import json
import pandas as pd
import numpy as np
import os

def main(json_file):
    # Check input
    if not json_file.endswith('.json'):
        print("Error: Input file must be a .json file.")
        return
    
    # Derive base name
    base_name = os.path.splitext(os.path.basename(json_file))[0]
    
    # Load JSON
    with open(json_file, 'r', encoding='utf-8-sig') as f:
        data = json.load(f)
        cells = data.get('cells', [])
        if not cells:
            print("Error: No 'cells' key or empty list in JSON.")
            return
    
    # Normalize and convert to DataFrame
    df = pd.DataFrame(cells)
    df = df[['X', 'Y', 'Z']]  # ensure order
    
    # Save CSV
    csv_file = f"{base_name}.csv"
    df.to_csv(csv_file, index=False)
    print(f"Saved coordinates to {csv_file}")
    

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 json2csv.py <input_file.json>")
    else:
        main(sys.argv[1])
