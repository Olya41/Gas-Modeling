
import numpy as np
from numba import njit
from forces import compute_accelerations, lennard_jones_potential
from typing import Tuple


def velocity_verlet_step(positions: np.ndarray,
                        velocities: np.ndarray,
                        accelerations: np.ndarray,
                        dt: float,
                        cell_size: float,
                        r_cut: float = 2.5,
                        **kwargs) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

    a_old = accelerations.copy()
    new_positions = positions + velocities * dt + 0.5 * a_old * dt**2
    new_accelerations = compute_accelerations(new_positions, cell_size, r_cut)
    new_velocities = velocities + 0.5 * (a_old + new_accelerations) * dt

    mask1 = new_positions > cell_size
    mask2 = new_positions < 0
    new_positions[mask1] -= cell_size
    new_positions[mask2] += cell_size

    return new_positions, new_velocities, new_accelerations


def position_verlet_step(positions: np.ndarray,
                        positions_prev: np.ndarray,
                        accelerations: np.ndarray,
                        dt: float,
                        cell_size: float,
                        **kwargs) -> Tuple[np.ndarray, np.ndarray]:

    disp = positions - positions_prev
    disp -= cell_size * np.round(disp / cell_size)

    new_positions = positions + disp + accelerations * dt**2

    new_disp = new_positions - positions
    velocities = (disp + new_disp) / (2 * dt)

    new_positions = new_positions % cell_size

    return new_positions, velocities


def kinetic_energy(velocities: np.ndarray) -> float:
    return 0.5 * np.sum(velocities ** 2)


@njit
def potential_energy(positions: np.ndarray, cell_size: float, r_cut: float = 2.5) -> float:
    N = len(positions)
    rc = min(r_cut, cell_size / 2)
    u_shift = lennard_jones_potential(rc)
    pe = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            r_vec = positions[j] - positions[i]
            for k in range(3):
                if r_vec[k] > cell_size / 2:
                    r_vec[k] -= cell_size
                if r_vec[k] < -cell_size / 2:
                    r_vec[k] += cell_size
            r = np.sqrt(r_vec[0]**2 + r_vec[1]**2 + r_vec[2]**2)
            if r < rc:
                pe += lennard_jones_potential(r) - u_shift
    return pe
