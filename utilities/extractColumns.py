import pandas as pd

# Load the CSV file
df = pd.read_csv('slice_3.csv')

# Keep only the 2nd and 3rd columns (X and Y)
xy = df.iloc[:, [1, 2]].copy()
xy.columns = ['X', 'Y']  # optional: rename for clarity

# Add Z column with zero values
xy.loc[:, 'Z'] = 0.0 

# Save to new file
xy.to_csv('slice_3_XYZ.csv', index=False)
