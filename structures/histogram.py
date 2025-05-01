import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV file
df = pd.read_csv('Min_Cell_Distances.csv')

# Plot histogram
plt.hist(df['min_distance'], bins=30, edgecolor='black')
plt.xlabel('Minimum Distance to Another Cell')
plt.ylabel('Frequency')
plt.title('Distribution of Minimum Inter-Cell Distances')
plt.grid(True)
plt.tight_layout()
plt.show()
