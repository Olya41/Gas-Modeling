"""
Зависимость λ(sim) от радиуса столкновения d (порог r_ij < d), d = (0.90…1.10) σ, шаг 0.02 σ.

Кадров берётся мало для быстрого прогона (см. N_ROWS_TAIL).

  python3 analysis/lambda_collision_sweep.py

Результат: output/data/lambda_vs_collision_radius.csv, output/plots/lambda_vs_collision_radius.png
"""
from __future__ import annotations

import io
import subprocess

import matplotlib.pyplot as plt
import numpy as np

from free_path_core import compute_lambda_sim
from sim_io import ROOT, read_param


# Уменьшенное окно траектории (последние строки) — главный рычаг времени счёта
N_ROWS_TAIL = 2000

# d / σ
THRESHOLDS_SIGMA = np.linspace(0.90, 1.10, 11)


def read_tail(path, n: int) -> np.ndarray:
    out = subprocess.run(
        ["tail", "-n", str(n), str(path)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return np.loadtxt(io.StringIO(out))


def main() -> None:
    rho = float(read_param("rho", "0.8"))
    dt_param = float(read_param("dt", "0.01"))

    pos_path = ROOT / "output/data/positions.txt"
    vel_path = ROOT / "output/data/velocities.txt"
    with open(pos_path) as f:
        first_row = f.readline().split()
        total_lines = 1 + sum(1 for _ in f)
    n = (len(first_row) - 1) // 3
    n_rows = min(N_ROWS_TAIL, total_lines)

    pos_raw = read_tail(pos_path, n_rows)
    vel_raw = read_tail(vel_path, n_rows)
    times = pos_raw[:, 0]
    m = len(times)
    positions = pos_raw[:, 1:].reshape(m, n, 3)
    velocities = vel_raw[:, 1:].reshape(m, n, 3)
    dt_sim = float(times[1] - times[0]) if m > 1 else dt_param
    cell_size = (n / rho) ** (1.0 / 3.0)

    print(
        f"λ(d): {m} кадров (хвост из {total_lines}), N = {n}, dt = {dt_sim:.5g}, "
        f"L = {cell_size:.5g}"
    )

    rows: list[tuple[float, float, float, int, float]] = []
    for d in THRESHOLDS_SIGMA:
        lam_mean, lam_std, n_seg, lam_th = compute_lambda_sim(
            positions, velocities, cell_size, dt_sim, d, rho
        )
        rows.append((d, lam_mean, lam_std, n_seg, lam_th))
        print(
            f"  d = {d:.2f} σ  →  λ(sim) = {lam_mean:.4g} ± {lam_std:.3g},  "
            f"сегментов = {n_seg},  λ(th) = {lam_th:.4g}"
        )

    out_csv = ROOT / "output/data/lambda_vs_collision_radius.csv"
    out_plot = ROOT / "output/plots/lambda_vs_collision_radius.png"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_plot.parent.mkdir(parents=True, exist_ok=True)

    header = "d_sigma,lambda_sim,lambda_std,n_segments,lambda_theory"
    with open(out_csv, "w", encoding="utf-8") as f:
        f.write(header + "\n")
        for d, lam_mean, lam_std, n_seg, lam_th in rows:
            f.write(
                f"{d:.6f},{lam_mean:.8f},{lam_std:.8f},{n_seg},{lam_th:.8f}\n"
            )
    print(f"\nCSV: {out_csv}")

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    d_arr = np.array([r[0] for r in rows])
    lam_sim_arr = np.array([r[1] for r in rows])
    lam_th_arr = np.array([r[4] for r in rows])
    ax.plot(d_arr, lam_sim_arr, "C0-o", lw=1.8, ms=6, label=r"$\lambda_{\mathrm{sim}}$")
    ax.plot(
        d_arr,
        lam_th_arr,
        "C1--",
        lw=1.5,
        label=r"$\lambda_{\mathrm{theory}} = 1/(\sqrt{2}\rho\pi d^2)$",
    )
    ax.set_xlabel(r"радиус столкновения $d$ ($\sigma$)")
    ax.set_ylabel(r"$\lambda$")
    ax.set_title(
        f"Средний свободный пробег vs порог сближения\n"
        f"{m} кадров, N = {n}, ρ = {rho}"
    )
    ax.grid(True, alpha=0.35)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_plot, dpi=150)
    plt.close(fig)
    print(f"график: {out_plot}")


if __name__ == "__main__":
    main()
