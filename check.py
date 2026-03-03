import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("velocities.txt")
N = (data.shape[1] - 1) // 3
last_row = data[-1]
t = last_row[0]
vel = last_row[1:].reshape(N, 3)

Ek = 0.5 * np.sum(vel**2)
T_kin = 2 * Ek / (3 * N)

v_all = vel.flatten()
v_arr = np.linspace(v_all.min() * 1.3, v_all.max() * 1.3, 200)
boltzmann_1d = (1 / (2 * np.pi * T_kin))**0.5 * np.exp(-v_arr**2 / (2 * T_kin))

plt.figure(figsize=(8, 5))
plt.hist(v_all, bins=30, density=True, alpha=0.7, label="v_x + v_y + v_z")
plt.plot(v_arr, boltzmann_1d, "r-", lw=2, label=f"Boltzmann (T = {T_kin:.3f})")
plt.xlabel("v")
plt.ylabel("p(v)")
plt.title(f"Velocity distribution, t = {t:.2f}")
plt.legend()
plt.tight_layout()
plt.savefig("vel_gist.png", dpi=150)
print(f"T_kin = {T_kin:.4f}, saved to vel_gist.png")

# Линеаризация: ln(P) vs v^2
counts, bin_edges = np.histogram(v_all, bins=30, density=True)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
mask = counts > 0
v2 = bin_centers[mask]**2
ln_p = np.log(counts[mask])

slope, intercept = np.polyfit(v2, ln_p, 1)
T_fit = -1 / (2 * slope)

plt.figure(figsize=(8, 5))
plt.scatter(v2, ln_p, label="Simulation")
plt.plot(v2, slope * v2 + intercept, "r-", lw=2, label=f"Fit: T = {T_fit:.3f}")
plt.xlabel(r"$v^2$")
plt.ylabel("ln p(v)")
plt.title(f"Linearization, t = {t:.2f}")
plt.legend()
plt.tight_layout()
plt.savefig("vel_linear.png", dpi=150)
print(f"T_fit = {T_fit:.4f}, saved to vel_linear.png")
