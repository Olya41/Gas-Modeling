"""
Радиальная плотность частиц вокруг одной «центральной» частицы
(в смысле: ближайшей к геометрическому центру ячейки [0, L)³, τ=0)
с УСРЕДНЕНИЕМ ПО ВРЕМЕНИ. Расстояния — минимальные образы (ПГУ), r < L/2.

  python3 analysis/density_plots.py

Параметры: см. `params` — rho, msd_t0 / dens_t0, dens_n_bins, dens_stride (опц.).
→ output/plots/density_from_center.png
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from sim_io import ROOT, read_param


def _read_positions(path) -> tuple[np.ndarray, np.ndarray, int]:
    d = np.loadtxt(path)
    times = d[:, 0]
    n = (d.shape[1] - 1) // 3
    pos = d[:, 1:].reshape(d.shape[0], n, 3)
    return pos, times, n


def _min_image_separation(r0: np.ndarray, r1: np.ndarray, l_box: float) -> np.ndarray:
    dr = r1 - r0
    return dr - l_box * np.round(dr / l_box)


def _min_image_from_point(
    p0: np.ndarray, pos_all: np.ndarray, l_box: float
) -> np.ndarray:
    """Мин. образы: (N,3) от p0 (3,) к каждой из pos_all."""
    dr = pos_all - p0[None, :]
    return dr - l_box * np.round(dr / l_box)


def _pick_central_index(pos_frame: np.ndarray, l_box: float) -> int:
    """Индекс частицы, ближайшей (ПГУ) к центру куба (L/2, L/2, L/2)."""
    c3 = np.full((1, 3), 0.5 * l_box)
    d = _min_image_separation(pos_frame, c3, l_box)
    return int(np.argmin(np.sum(d * d, axis=1)))


def compute_radial_density(
    root=ROOT,
    n_bins: int | None = None,
    time_stride: int | None = None,
    t0: float | None = None,
) -> dict:
    """
    n(r) = (число пар в оболочке, усреднённое по t) / V_оболочки, остальные N−1
    соседей вокруг фиксированного центра.
    """
    rho = float(read_param("rho", "0.8"))
    if n_bins is None:
        _b = read_param("dens_n_bins", None)
        n_bins = int(_b) if _b is not None else 80
    if time_stride is None:
        _s = read_param("dens_stride", None)
        time_stride = int(_s) if _s is not None else 1
    if time_stride < 1:
        time_stride = 1

    pos, times, n = _read_positions(root / "output/data/positions.txt")
    m = pos.shape[0]
    t_span = float(times[-1] - times[0]) if m > 1 else 0.0
    msd_t0 = read_param("msd_t0", None)
    if t0 is not None:
        t_start = float(t0)
    else:
        _d0 = read_param("dens_t0", None)
        if _d0 is not None:
            t_start = float(_d0)
        elif msd_t0 is not None:
            t_start = float(msd_t0)
        else:
            t_start = 0.1 * t_span
    i0 = int(np.searchsorted(times, t_start, side="left"))
    i0 = min(max(i0, 0), m - 1)

    l_box = (n / rho) ** (1.0 / 3.0)
    r_max = 0.5 * l_box * 0.999
    r_edges = np.linspace(0.0, r_max, n_bins + 1)
    r_cent = 0.5 * (r_edges[1:] + r_edges[:-1])
    v_shells = (4.0 * np.pi / 3.0) * (r_edges[1:] ** 3 - r_edges[:-1] ** 3)
    v_shells = np.where(v_shells < 1e-30, np.nan, v_shells)

    central_idx = _pick_central_index(pos[i0], l_box)

    hist = np.zeros(n_bins, dtype=np.float64)
    n_used = 0
    mask_others = np.ones(n, dtype=bool)
    mask_others[central_idx] = False
    for t in range(i0, m, time_stride):
        p0 = pos[t, central_idx]
        dr = _min_image_from_point(p0, pos[t], l_box)
        r = np.linalg.norm(dr, axis=1)
        rs = r[mask_others]
        rs = rs[rs < r_max]
        c, _ = np.histogram(rs, bins=r_edges)
        hist += c.astype(np.float64)
        n_used += 1
    if n_used < 1:
        raise SystemExit("Нет кадров для усреднения (dens_t0 слишком велико?).")
    mean_counts = hist / n_used
    n_density = mean_counts / v_shells
    rho_mean = n / (l_box**3)

    return {
        "N": n,
        "rho": rho,
        "L": l_box,
        "i_start": i0,
        "t_start": float(times[i0]),
        "i_end": m - 1,
        "t_end": float(times[-1]),
        "n_frames_used": n_used,
        "time_stride": time_stride,
        "central_index": central_idx,
        "n_bins": n_bins,
        "r_edges": r_edges,
        "r_center": r_cent,
        "n_density": n_density,
        "rho_mean": rho_mean,
    }


def plot_density(
    d: dict | None = None,
    out_path=None,
    root=ROOT,
) -> None:
    d = d or compute_radial_density(root=root)
    out_path = out_path or (root / "output/plots/density_from_center.png")

    r = d["r_center"]
    n_r = d["n_density"]
    rho0 = d["rho_mean"]
    g_r = n_r / rho0
    n_frames = d["n_frames_used"]

    fs_label, fs_title, fs_tick = 14, 15, 12

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    m = np.isfinite(n_r) & (r > 0)
    ax1.plot(r[m], n_r[m], "C0-", lw=1.8, label=r"$\langle n(r)\rangle$")
    ax1.axhline(rho0, color="0.4", ls="--", lw=1.2, label=r"$\rho_0 = N/L^3$")
    ax1.set_ylabel(r"локальная плотность $n(r)$", fontsize=fs_label)
    ax1.grid(True, alpha=0.35)
    i0 = d["central_index"]
    ax1.set_title(
        rf"плотность на расстоянии от центр. частицы $i = {i0}$ "
        rf"({d['n_frames_used']} снимков, stride = {d['time_stride']})",
        fontsize=fs_title,
    )
    ax1.legend(loc="best", fontsize=10)
    ax1.tick_params(axis="both", labelsize=fs_tick)

    ax2.plot(r[m], g_r[m], "C1-", lw=1.8, label=r"$\langle n(r)\rangle / \rho_0$")
    ax2.axhline(1.0, color="0.4", ls="--", lw=1.0)
    ax2.set_xlabel(r"расстояние $r$ (ПГУ), ед. $\sigma$", fontsize=fs_label)
    ax2.set_ylabel(r"$n(r)/\rho_0$", fontsize=fs_label)
    ax2.grid(True, alpha=0.35)
    ax2.legend(loc="best", fontsize=10)
    ax2.tick_params(axis="both", labelsize=fs_tick)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"сохранено: {out_path}")


def main() -> None:
    d = compute_radial_density()
    n = d["N"]
    l_box = d["L"]
    print("--- радиальная плотность вокруг одной «центральной» частицы ---")
    print(f"  N = {n},  rho = {d['rho']},  L = (N/rho)^{{1/3}} = {l_box:.6f}")
    print(
        f"  кадр начала: i = {d['i_start']},  t = {d['t_start']:.4f}  "
        f"(dens_t0 / msd_t0 / 0.1·T, см. код)"
    )
    print(
        f"  конец: t = {d['t_end']:.4f},  снимков в среднем: {d['n_frames_used']},  "
        f"dens_stride = {d['time_stride']}"
    )
    print(
        f"  центр. частица: индекс i = {d['central_index']}  "
        f"(ближайшая (ПГУ) к (L/2, L/2, L/2) на t_start)"
    )
    print(
        f"  радиус сетки: 0…{0.5 * l_box * 0.999:.4f}  (< L/2),  bins = {d['n_bins']}"
    )
    print(f"  rho_0 = N/L^3 = {d['rho_mean']:.6f}")
    plot_density(d=d)
    m = np.isfinite(d["n_density"])
    k = m & (d["r_center"] > 0)
    r_mid = float(np.median(d["r_center"][k]))
    n_mid = float(np.interp(r_mid, d["r_center"][k], d["n_density"][k]))
    print(f"  пример: n(r≈{r_mid:.2f}) / rho_0 ≈ {n_mid / d['rho_mean']:.3f}")


if __name__ == "__main__":
    main()
