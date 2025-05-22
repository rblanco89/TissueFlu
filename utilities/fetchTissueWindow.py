import pandas as pd
import numpy as np

# Load the scaled coordinates
df = pd.read_csv('cell_coords_scaled.csv')

# Compute bounding box
x_min, x_max = df['X'].min(), df['X'].max()
y_min, y_max = df['Y'].min(), df['Y'].max()
z_min, z_max = df['Z'].min(), df['Z'].max()

# Define 10% box dimensions
x_len = x_max - x_min
y_len = y_max - y_min
z_len = z_max - z_min

# Define the origin of the window (e.g., bottom corner of the tissue)
# You can randomize this to get different windows
x0 = 0.25*x_len
y0 = 0.25*y_len
z0 = 0.25*z_len

# Window range
x1 = x0 + 0.5 * x_len
y1 = y0 + 0.5 * y_len
z1 = z0 + 0.5 * z_len

print(f"x0 = {x0}; x1 = {x1}")
print(f"y0 = {y0}; y1 = {y1}")
print(f"z0 = {z0}; z1 = {z1}")

# Filter points inside the window
window_df = df[(df['X'] >= x0) & (df['X'] <= x1) &
               (df['Y'] >= y0) & (df['Y'] <= y1) &
               (df['Z'] >= z0) & (df['Z'] <= z1)]
window_df.head()

# Save the result
window_df.to_csv('cell_coords_window.csv', index=False)
print(f"Selected {len(window_df)} cells from 10% spatial window.")
