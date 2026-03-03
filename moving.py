import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
from verlet import velocity_verlet_step, kinetic_energy, potential_energy
from forces import compute_accelerations


def init_lattice(N: int, cell_size: float) -> np.ndarray:
    n = int(np.ceil(N ** (1/3)))
    spacing = cell_size / n
    coords = []
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                if len(coords) < N:
                    coords.append([(ix + 0.5) * spacing,
                                   (iy + 0.5) * spacing,
                                   (iz + 0.5) * spacing])
    return np.array(coords)


def init_velocities(N: int, T_target: float = 1.0) -> np.ndarray:
    vel = np.random.normal(0, np.sqrt(T_target), (N, 3))
    vel -= vel.mean(axis=0)
    return vel


def moving(positions0: np.ndarray,
           velocities0: np.ndarray,
           cell_size: float,
           NUM: int = 1000,
           dt: float = 0.01,
           r_cut: float = 2.5,
           energy_list: bool = False,
           **kwargs) -> Tuple[np.ndarray, np.ndarray]:

    pos = positions0.copy()
    vel = velocities0.copy()
    acc = compute_accelerations(pos, cell_size, r_cut)
    n = len(pos)

    positions_traj = np.zeros((NUM + 1, n, 3))
    velocities_traj = np.zeros((NUM + 1, n, 3))
    energies = np.zeros((NUM + 1, 3)) if energy_list else None

    positions_traj[0] = pos
    velocities_traj[0] = vel

    if energy_list:
        ek = kinetic_energy(vel)
        pe = potential_energy(pos, cell_size, r_cut)
        T_kin = 2 * ek / (3 * n)
        energies[0] = [ek, pe, ek + pe]
        print(f"Step 0/{NUM}: Ek={ek:.4f}, Ep={pe:.4f}, Etot={ek+pe:.4f}, T={T_kin:.4f}")

    for t in range(1, NUM + 1):
        pos, vel, acc = velocity_verlet_step(pos, vel, acc, dt, cell_size, r_cut)
        positions_traj[t] = pos
        velocities_traj[t] = vel
        if energy_list:
            ek = kinetic_energy(vel)
            pe = potential_energy(pos, cell_size, r_cut)
            T_kin = 2 * ek / (3 * n)
            energies[t] = [ek, pe, ek + pe]
        if t % 100 == 0:
            if energy_list:
                print(f"Step {t}/{NUM}: Ek={ek:.4f}, Ep={pe:.4f}, Etot={ek+pe:.4f}, T={T_kin:.4f}")
            else:
                print(f"Step {t}/{NUM}")

    times = np.arange(NUM + 1) * dt
    pos_with_time = np.column_stack([times, positions_traj.reshape(NUM + 1, n * 3)])
    vel_with_time = np.column_stack([times, velocities_traj.reshape(NUM + 1, n * 3)])

    np.savetxt("positions.txt", pos_with_time)
    np.savetxt("velocities.txt", vel_with_time)

    if energies is not None:
        np.savetxt("energies.txt", np.column_stack([times, energies]),
                   header="t    Ek    Ep    Etot")
        plt.figure(figsize=(10, 6))
        plt.plot(times, energies[:, 0], label="Ek")
        plt.plot(times, energies[:, 1], label="Ep")
        plt.plot(times, energies[:, 2], label="Etot", lw=2, color="black")
        plt.xlabel("t")
        plt.ylabel("E")
        plt.title("E(t)")
        plt.legend()
        plt.tight_layout()
        plt.savefig("energy.png", dpi=150)
    else:
        open("energies.txt", "w").close()
        open("energy.png", "w").close()

    print("Saved to positions.txt, velocities.txt, energies.txt, energy.png")

    return positions_traj, velocities_traj, energies
