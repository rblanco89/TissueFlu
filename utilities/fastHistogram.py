import numpy as np
import matplotlib.pyplot as plt

# Load your data
data = np.loadtxt('tmp.txt')

# Plot histogram
plt.figure(figsize=(8, 5))
plt.hist(data, bins=range(int(data.min()), int(data.max()) + 2), edgecolor='black', align='left')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of Poisson Random Numbers')
plt.grid(True)
plt.tight_layout()
plt.show()
