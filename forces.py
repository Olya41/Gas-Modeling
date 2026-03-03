import numpy as np


def lennard_jones_potential(r: float) -> float:
    return 4 * ((1 / r) ** 12 - (1 / r) ** 6)


def lennard_jones_force(r: float) -> float:
    return 24 *(2 * (1 / r) ** 12 - (1 / r) ** 6) / r


def compute_accelerations(coords: np.ndarray,
                          cell_size,
                          r_cut: float = 2.5) -> np.ndarray:

    N = len(coords)
    accelerations = np.zeros_like(coords)
    rc = min(r_cut, cell_size / 2)

    for i in range(N):
        for j in range(i + 1, N):
            r_vec = (coords[j] - coords[i]).copy()

            for k in range(len(r_vec)):
                if r_vec[k] > cell_size / 2:
                    r_vec[k] -= cell_size
                if r_vec[k] < -cell_size / 2:
                    r_vec[k] += cell_size

            r = np.linalg.norm(r_vec)
            if r < rc:
                force_vec = lennard_jones_force(r) * r_vec / r
                accelerations[i] -= force_vec
                accelerations[j] += force_vec

    return accelerations











