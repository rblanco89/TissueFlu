from sklearn.neighbors import NearestNeighbors
import numpy as np
import pandas as pd

# Load previously saved CSV
df = pd.read_csv('Cell_Coords_Normalized.csv')
coords = df[['X', 'Y', 'Z']].values

# Fit NearestNeighbors (k=2 because first neighbor is itself)
nn = NearestNeighbors(n_neighbors=2, algorithm='auto', metric='euclidean')
nn.fit(coords)

# Compute distances to the two nearest neighbors
distances, indices = nn.kneighbors(coords)

# Take the 2nd column (distance to the closest *other* point)
min_distances = distances[:, 1]

# Save to CSV
pd.DataFrame({'min_distance': min_distances}).to_csv('Min_Cell_Distances.csv', index=False)
