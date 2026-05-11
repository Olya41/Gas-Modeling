"""
Ядро детектирования перицентров и суммирования длин свободного пробега
(логика как в free_run.compute_free_paths, без привязки к глобальным траекториям).
"""
from __future__ import annotations

from itertools import combinations

import numpy as np


def compute_free_path_lengths(
    positions: np.ndarray,
    velocities: np.ndarray,
    cell_size: float,
    dt_sim: float,
    threshold: float,
) -> list[float]:
    """
    Ищет интервалы r_ij < threshold, перицентры внутри полных интервалов,
    для каждой частицы суммирует ∫|v|dt между последовательными перицентрами.

    positions: (M, N, 3), velocities: (M, N, 3), threshold — порог центр–центр в ед. σ.
    """
    m_steps, n, _ = positions.shape
    if m_steps != velocities.shape[0] or velocities.shape[1] != n:
        raise ValueError("positions и velocities должны совпасть по (M, N)")
    particle_peri: list[list[int]] = [[] for _ in range(n)]

    for i, j in combinations(range(n), 2):
        dr = positions[:, j, :] - positions[:, i, :]
        dr -= cell_size * np.round(dr / cell_size)
        dists = np.linalg.norm(dr, axis=1)

        below = dists < threshold
        if not np.any(below):
            continue

        padded = np.empty(m_steps + 2, dtype=bool)
        padded[0] = False
        padded[1:-1] = below
        padded[-1] = False
        diff = np.diff(padded.astype(np.int8))
        starts = np.where(diff == 1)[0]
        ends = np.where(diff == -1)[0]

        for s, e in zip(starts, ends):
            if s == 0 or e == m_steps:
                continue
            peri = s + int(np.argmin(dists[s:e]))
            particle_peri[i].append(peri)
            particle_peri[j].append(peri)

    free_paths: list[float] = []
    for ii in range(n):
        peri_sorted = sorted(set(particle_peri[ii]))
        if len(peri_sorted) < 2:
            continue
        for k in range(len(peri_sorted) - 1):
            a, b = peri_sorted[k], peri_sorted[k + 1]
            if b <= a + 1:
                continue
            speeds = np.linalg.norm(velocities[a:b, ii, :], axis=1)
            free_paths.append(float(np.sum(speeds) * dt_sim))

    return free_paths


def compute_lambda_sim(
    positions: np.ndarray,
    velocities: np.ndarray,
    cell_size: float,
    dt_sim: float,
    threshold: float,
    rho: float,
) -> tuple[float, float, int, float]:
    """
    Средняя длина свободного пробега λ(sim), σ по сегментам, их число,
    и λ(theory) = 1/(√2·ρ·π·d²), d = threshold.
    При отсутствии сегментов — (nan, nan, 0, theory).
    """
    paths = compute_free_path_lengths(
        positions, velocities, cell_size, dt_sim, threshold
    )
    coll_section = np.pi * threshold**2
    lambda_theory = 1.0 / max(np.sqrt(2) * rho * coll_section, 1e-30)
    if not paths:
        return float("nan"), float("nan"), 0, float(lambda_theory)
    arr = np.asarray(paths, dtype=np.float64)
    return float(np.mean(arr)), float(np.std(arr)), int(arr.size), float(lambda_theory)
