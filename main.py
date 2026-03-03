import numpy as np
from moving import init_lattice, moving

# --- Параметры симуляции ---
N = 200
dt = 0.01
NUM = 4000
rho = 0.3
r_cut = 2.5
energy_list = True

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
