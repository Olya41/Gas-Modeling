from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Корень репозитория (запуск из analysis/ или из корня)
ROOT = Path(__file__).resolve().parent.parent


# Читаем NUM_FOR_HIST из params
def read_param(name, default=None):
    with open(ROOT / "params") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                if k.strip() == name:
                    return v.strip()
    return default

raw = read_param("NUM_FOR_HIST", "1000")
if raw == "NUM":
    raw = read_param("NUM", "1000")
NUM_FOR_HIST = int(raw)

data = np.loadtxt(ROOT / "output/data/velocities.txt")
N = (data.shape[1] - 1) // 3
total_rows = data.shape[0]

n_snapshots = min(NUM_FOR_HIST, total_rows)

rows = data[-n_snapshots:]
t_first = rows[0, 0]
t_last = rows[-1, 0]
vel_all_snapshots = rows[:, 1:].reshape(-1, 3)

Ek = 0.5 * np.sum(vel_all_snapshots**2) / n_snapshots
T_kin = 2 * Ek / (3 * N)

v_all = vel_all_snapshots.flatten()
v_arr = np.linspace(v_all.min() * 1.3, v_all.max() * 1.3, 200)
boltzmann_1d = (1 / (2 * np.pi * T_kin))**0.5 * np.exp(-v_arr**2 / (2 * T_kin))

plt.figure(figsize=(8, 5))
plt.hist(v_all, bins=60, density=True, alpha=0.7, label="v_x + v_y + v_z")
plt.plot(v_arr, boltzmann_1d, "r-", lw=2, label=f"Boltzmann (T_kin = {T_kin:.3f})")
plt.xlabel("v")
plt.ylabel("p(v)")
plt.title("Velocity distribution")
plt.legend()
plt.tight_layout()
plt.savefig(ROOT / "output/plots/vel_gist.png", dpi=150)
print(f"T_kin = {T_kin:.4f}, {n_snapshots} snapshots, t = [{t_first:.2f}, {t_last:.2f}]")

# Линеаризация: ln(P) vs v^2
counts, bin_edges = np.histogram(v_all, bins=60, density=True)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
mask = counts > 0
v2 = bin_centers[mask]**2
ln_p = np.log(counts[mask])

slope_kin = -1 / (2 * T_kin)
intercept_kin = 0.5 * np.log(1 / (2 * np.pi * T_kin))

plt.figure(figsize=(8, 5))
plt.scatter(v2, ln_p, label="Simulation")
plt.plot(v2, slope_kin * v2 + intercept_kin, "r-", lw=2, label=f"Theory: T_kin = {T_kin:.3f}")
plt.xlabel(r"$v^2$")
plt.ylabel("ln p(v)")
plt.title("Linearization")
plt.legend()
plt.tight_layout()
plt.savefig(ROOT / "output/plots/vel_linear.png", dpi=150)

# Дополнительный график с пунктиром на 2*max(Ek)/m и продленной теоретической линией
# Рассчитываем максимальную кинетическую энергию
Ek_max = np.max(np.sum(0.5 * rows[:, 1:].reshape(n_snapshots, N, 3)**2, axis=(1,2)))
v2_max = 2*Ek_max  # v² = 2E_max/m

# Создаем расширенный массив v² для теоретической линии до пунктира
v2_extended = np.linspace(0, v2_max, 200)
theory_line = slope_kin * v2_extended + intercept_kin

plt.figure(figsize=(8, 5))
plt.scatter(v2, ln_p, label="Simulation")
plt.plot(v2_extended, theory_line, "r-", lw=2, label=f"Theory: T_kin = {T_kin:.3f}")
plt.axvline(x=v2_max, color='red', linestyle='--', label=f'v² = 2E_max/m = {v2_max:.2f}')
plt.xlabel(r"$v^2$")
plt.ylabel("ln p(v)")
plt.title("Linearization with max energy line")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(ROOT / "output/plots/vel_linear_with_max.png", dpi=150)

# Energy vs time
edata = np.loadtxt(ROOT / "output/data/energies.txt")
times = edata[:, 0]
Ek_t  = edata[:, 1]
Ep_t  = edata[:, 2]
Etot  = edata[:, 3]

plt.figure(figsize=(10, 6))
plt.plot(times, Ek_t, label="Ek")
plt.plot(times, Ep_t, label="Ep")
plt.plot(times, Etot, label="Etot", lw=2, color="black")
plt.xlabel("t")
plt.ylabel("E")
plt.title("E(t)")
plt.legend()
plt.tight_layout()
plt.savefig(ROOT / "output/plots/energy.png", dpi=150)
print("Saved output/plots/energy.png")

# Total energy separately
plt.figure(figsize=(10, 6))
plt.plot(times, Etot, color='black')
plt.xlabel('t')
plt.ylabel('Etot')
plt.title('Total Energy E(t)')
plt.tight_layout()
plt.savefig(ROOT / "output/plots/total_energy.png", dpi=150)
print("Saved output/plots/total_energy.png")
