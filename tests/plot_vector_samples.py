import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("data/vector_samples.csv")

# 3D scatter
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection="3d")
ax.scatter(df.x, df.y, df.z, s=1, alpha=0.3)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Random unit vectors within cone")

plt.tight_layout()
plt.show()

# Check histograms
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
theta = np.arccos(df.z)
phi = np.arctan2(df.y, df.x)
axs[0].hist(np.cos(theta), bins=50, density=True)
axs[0].set_xlabel("cos(theta)")
axs[0].set_ylabel("Density")
axs[1].hist(phi, bins=50, density=True)
axs[1].set_xlabel("phi")
axs[1].set_ylabel("Density")
plt.tight_layout()
plt.show()
