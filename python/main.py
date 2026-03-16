import numpy as np
from moving import init_lattice, moving

# --- Загрузка параметров ---
def load_params(filename="params"):
    p = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            key, _, val = line.partition('=')
            p[key.strip()] = val.strip()
    return p

par = load_params()
N          = int(par["N"])
dt         = float(par["dt"])
NUM        = int(par["NUM"])
rho        = float(par["rho"])
r_cut      = float(par["r_cut"])
energy_list = par["energy_log"].lower() == "true"

# --- Вычисляемые параметры ---
cell_size = (N / rho) ** (1/3)

# --- Начальные условия ---
coord_0 = init_lattice(N, cell_size)
velocities_0 = np.zeros((N, 3))

# --- Вывод параметров ---
print("=== MD Lennard-Jones ===")
print(f"N = {N}")
print(f"dt = {dt}")
print(f"NUM = {NUM} steps (total time {NUM * dt:.2f})")
print(f"rho = {rho}")
print(f"cell_size = {cell_size:.4f}")
print(f"r_cut = {r_cut}")
print(f"energy_list = {energy_list}")
print()

# --- Запуск ---
moving(coord_0, velocities_0, cell_size, NUM, dt, r_cut=r_cut, energy_list=energy_list)
