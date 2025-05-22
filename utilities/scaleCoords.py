import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

# Load previously saved CSV
df = pd.read_csv('Cell_Coords_Normalized.csv')
coords = df[['X', 'Y', 'Z']].values

# Compute nearest neighbor distances
nn = NearestNeighbors(n_neighbors=2, algorithm='auto', metric='euclidean')
nn.fit(coords)
distances, _ = nn.kneighbors(coords)

# The second column contains distance to nearest *other* cell
min_distances = distances[:, 1]
global_min = np.min(min_distances)

# Scale coordinates so global minimum becomes 1
scaled_coords = coords / global_min

# Save scaled coordinates
df_scaled = pd.DataFrame(scaled_coords, columns=['X', 'Y', 'Z'])
df_scaled.to_csv('cell_coords_scaled.csv', index=False)
