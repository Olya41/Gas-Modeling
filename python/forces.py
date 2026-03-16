import numpy as np
from numba import njit


@njit
def lennard_jones_potential(r: float) -> float:
    return 4 * ((1 / r) ** 12 - (1 / r) ** 6)


@njit
def lennard_jones_force(r: float) -> float:
    return 24 * (2 * (1 / r) ** 12 - (1 / r) ** 6) / r


@njit
def compute_accelerations(coords: np.ndarray,
                          cell_size: float,
                          r_cut: float = 2.5) -> np.ndarray:
    N = len(coords)
    accelerations = np.zeros_like(coords)
    rc = min(r_cut, cell_size / 2)

    for i in range(N):
        for j in range(i + 1, N):
            r_vec = coords[j] - coords[i]

            for k in range(3):
                if r_vec[k] > cell_size / 2:
                    r_vec[k] -= cell_size
                if r_vec[k] < -cell_size / 2:
                    r_vec[k] += cell_size

            r = np.sqrt(r_vec[0]**2 + r_vec[1]**2 + r_vec[2]**2)
            if r < rc:
                f = lennard_jones_force(r)
                for k in range(3):
                    fk = f * r_vec[k] / r
                    accelerations[i, k] -= fk
                    accelerations[j, k] += fk

    return accelerations
