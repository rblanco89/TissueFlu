import sys
import json
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
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
    
#    # Compute pairwise distances
#    positions = df[['X', 'Y', 'Z']].values
#    distance_matrix = cdist(positions, positions)
#    np.fill_diagonal(distance_matrix, np.inf)
#    min_distances = np.min(distance_matrix, axis=1)
#    
#    # Save distances
#    dist_file = f"{base_name}_MinDistances.csv"
#    pd.DataFrame({'min_distance': min_distances}).to_csv(dist_file, index=False)
#    print(f"Saved minimum distances to {dist_file}")
#    
#    # Plot histogram
#    plt.hist(min_distances, bins=30, edgecolor='black')
#    plt.xlabel('Minimum Distance to Another Cell')
#    plt.ylabel('Frequency')
#    plt.title('Distribution of Minimum Inter-Cell Distances')
#    plt.grid(True)
#    plt.tight_layout()
#    plt.savefig("minDistance_distribution.png")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 json2csv.py <input_file.json>")
    else:
        main(sys.argv[1])
