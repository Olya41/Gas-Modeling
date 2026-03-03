import numpy as np

data = np.loadtxt("positions.txt")
NUM = data.shape[0]
N = (data.shape[1] - 1) // 3

with open("trajectory.xyz", "w") as f:
    for i in range(NUM):
        t = data[i, 0]
        coords = data[i, 1:].reshape(N, 3)
        f.write(f"{N}\n")
        f.write(f"Time={t:.4f}\n")
        for j in range(N):
            f.write(f"Ar {coords[j, 0]:.6f} {coords[j, 1]:.6f} {coords[j, 2]:.6f}\n")

print(f"Saved trajectory.xyz ({NUM} frames, {N} particles)")
