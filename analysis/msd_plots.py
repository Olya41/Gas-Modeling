"""
MSD по распрямлённой траектории после прогрева (msd_t0).

Быстрый перерисовывание графиков без тяжёлого free_run (столкновения):
  python3 analysis/msd_plots.py
из корня репозитория (или с PYTHONPATH=…).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import matplotlib.pyplot as plt

from sim_io import ROOT, read_param


def _unwrap_positions(
    pos_wrapped: np.ndarray, cell_size: float
) -> np.ndarray:
    M_full = pos_wrapped.shape[0]
    out = pos_wrapped.copy()
    for t in range(1, M_full):
        delta = pos_wrapped[t] - pos_wrapped[t - 1]
        delta -= cell_size * np.round(delta / cell_size)
        out[t] = out[t - 1] + delta
    return out


def compute_msd(
    root: Path = ROOT,
) -> dict[str, Any]:
    """Загрузка траектории, распрямление, τ, MSD, D, печатные величины."""
    rho = float(read_param("rho", "0.8"))
    dt = float(read_param("dt", "0.01"))

    pos_full = np.loadtxt(root / "output/data/positions.txt")
    times_full = pos_full[:, 0]
    M_full = len(times_full)
    N = (pos_full.shape[1] - 1) // 3
    pos_wrapped_full = pos_full[:, 1:].reshape(M_full, N, 3)

    cell_size = (N / rho) ** (1.0 / 3.0)
    pos_unwrap_full = _unwrap_positions(pos_wrapped_full, cell_size)

    _span = float(times_full[-1]) - float(times_full[0])
    _msd_t0_raw = read_param("msd_t0", None)
    msd_t0 = float(_msd_t0_raw) if _msd_t0_raw is not None else 0.1 * _span
    i_msd0 = int(np.searchsorted(times_full, msd_t0, side="left"))
    i_msd0 = min(i_msd0, M_full - 1)
    t_msd_ref = float(times_full[i_msd0])

    displacement = pos_unwrap_full - pos_unwrap_full[i_msd0]
    msd_all = np.mean(np.sum(displacement**2, axis=2), axis=1)

    vel_full = np.loadtxt(root / "output/data/velocities.txt")
    M_vel = vel_full.shape[0]
    n_msd = M_full - i_msd0
    if M_vel < M_full:
        n_msd = min(n_msd, M_vel - i_msd0)
    if n_msd < 3:
        raise SystemExit(
            f"Слишком мало кадров для MSD после t_ref: n_msd = {n_msd}. "
            f"Уменьшите msd_t0 или удлините прогон."
        )

    tau = times_full[i_msd0 : i_msd0 + n_msd] - t_msd_ref
    msd = msd_all[i_msd0 : i_msd0 + n_msd]

    slope_msd, intercept_msd = np.polyfit(tau, msd, 1)
    D_from_msd = slope_msd / 6.0

    V = vel_full[i_msd0 : i_msd0 + n_msd, 1:].reshape(n_msd, N, 3)
    mean_speed = float(np.mean(np.linalg.norm(V, axis=2)))
    lambda_from_msd = 3.0 * D_from_msd / mean_speed

    dt_msd = float(times_full[1] - times_full[0]) if M_full > 1 else dt

    return {
        "N": N,
        "M_full": M_full,
        "n_msd": n_msd,
        "t_msd_ref": t_msd_ref,
        "t_first": float(times_full[0]),
        "t_last": float(times_full[-1]),
        "dt_msd": dt_msd,
        "msd_t0": msd_t0,
        "i_msd0": i_msd0,
        "tau": tau,
        "msd": msd,
        "slope_msd": float(slope_msd),
        "intercept_msd": float(intercept_msd),
        "D_from_msd": D_from_msd,
        "mean_speed": mean_speed,
        "lambda_from_msd": lambda_from_msd,
    }


def plot_msd(
    out_path: Path | None = None,
    root: Path = ROOT,
    zoom_tau: float | None = None,
) -> None:
    """
    out_path: куда сохранить (по умолчанию output/plots/msd.png).
    zoom_tau: ширина нижнего графика; None → params msd_zoom_tau, иначе 10.
    """
    if zoom_tau is None:
        _z = read_param("msd_zoom_tau", None)
        zoom_tau = float(_z) if _z is not None else 10.0

    d = compute_msd(root=root)
    out_path = out_path or (root / "output/plots/msd.png")

    print(
        f"MSD (после прогрева): {d['M_full']} снимков,  "
        f"t ∈ [{d['t_first']:.2f}, {d['t_last']:.2f}],  dt = {d['dt_msd']:.4f},  "
        f"msd_t0 = {d['msd_t0']:.4f} → t_ref = {d['t_msd_ref']:.4f},  "
        f"точек MSD = {d['n_msd']}"
    )

    tau = d["tau"]
    msd = d["msd"]
    t_msd_ref = d["t_msd_ref"]
    slope_msd = d["slope_msd"]
    intercept_msd = d["intercept_msd"]
    msd_t0 = d["msd_t0"]
    D_from_msd = d["D_from_msd"]
    mean_speed = d["mean_speed"]
    lambda_from_msd = d["lambda_from_msd"]
    tau_max = float(tau[-1])

    fs_label, fs_title, fs_tick = 15, 16, 13

    fig, (ax_full, ax_zoom) = plt.subplots(2, 1, figsize=(12, 10), sharex=False)

    ax_full.plot(tau, msd, lw=1.8, color="C0", label="MSD")
    msd_fit = slope_msd * tau + intercept_msd
    ax_full.plot(
        tau,
        msd_fit,
        "k--",
        lw=1.8,
        label=rf"лин. весь интервал: ${slope_msd:.4f}\,\tau + {intercept_msd:.2f}$",
    )
    ax_full.set_ylabel(
        r"$\langle |r(t_{\mathrm{ref}}+\tau) - r(t_{\mathrm{ref}})|^2 \rangle$",
        fontsize=fs_label,
    )
    ax_full.tick_params(axis="both", labelsize=fs_tick)
    ax_full.grid(True, alpha=0.35)
    ax_full.set_title(
        rf"MSD после прогрева: $t_{{\mathrm{{ref}}}} = {t_msd_ref:.4f}$",
        fontsize=fs_title,
    )
    ax_full.set_xlabel(r"$\tau = t - t_{\mathrm{ref}}$", fontsize=fs_label)
    ax_full.legend(loc="upper left", fontsize=11)

    # Нижняя панель: MSD + парабола (фит по [0,2], рисуем до τ=5)
    t_zoom_hi = min(zoom_tau, tau_max)
    i1 = int(np.searchsorted(tau, t_zoom_hi, side="right"))
    i1 = max(2, min(i1, len(tau)))

    ax_zoom.plot(tau[:i1], msd[:i1], lw=2.0, color="C0", label="MSD")

    t_par_hi = min(2.0, tau_max)  # интервал для polyfit
    t_par_draw_hi = min(5.0, t_zoom_hi)  # рисуем ту же квадрику до 5 (не шире окна)
    m_par = (tau >= 0.0) & (tau <= t_par_hi)
    if np.count_nonzero(m_par) >= 3:
        tp, mp = tau[m_par], msd[m_par]
        c2, c1, c0 = np.polyfit(tp, mp, 2)
        tau_par = np.linspace(0.0, t_par_draw_hi, 200)
        y_par = ((c2 * tau_par) + c1) * tau_par + c0
        ax_zoom.plot(
            tau_par,
            y_par,
            color="C2",
            lw=1.8,
            alpha=0.95,
            label=rf"парабола фит $[0,\,2]$, рис. $[0,\,{t_par_draw_hi:.0f}]$: "
            rf"${c2:.3f}\tau^2+{c1:.3f}\tau+{c0:.2f}$",
        )
    else:
        c2 = c1 = c0 = float("nan")

    ax_zoom.set_ylabel(
        r"$\langle |r(t_{\mathrm{ref}}+\tau) - r(t_{\mathrm{ref}})|^2 \rangle$",
        fontsize=fs_label,
    )
    ax_zoom.tick_params(axis="both", labelsize=fs_tick)
    ax_zoom.grid(True, alpha=0.35)
    ax_zoom.set_xlim(0.0, t_zoom_hi)
    _y_list = [msd[:i1]]
    if np.count_nonzero(m_par) >= 3:
        _y_list.append(y_par)
    _y = np.concatenate(_y_list)
    _ymin, _ymax = float(_y.min()), float(_y.max())
    if _ymax > _ymin:
        _pad = 0.04 * (_ymax - _ymin)
        ax_zoom.set_ylim(_ymin - _pad, _ymax + _pad)
    ax_zoom.set_title(
        fr"Увеличение: $\tau \in [0, {t_zoom_hi:.0f}]$",
        fontsize=fs_title,
    )
    ax_zoom.set_xlabel(r"$\tau = t - t_{\mathrm{ref}}$", fontsize=fs_label)
    ax_zoom.legend(loc="upper left", fontsize=10)

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=200)
    plt.close()

    print("--- MSD: линейная аппроксимация по τ после t_ref ---")
    print(f"t_ref = {t_msd_ref:.6f}  (порог msd_t0 = {msd_t0:.6f})")
    print(f"d(MSD)/dτ (наклон) = {slope_msd:.6f}")
    print(f"D = (d MSD/dτ) / 6 = {D_from_msd:.6f}")
    print(f"<|v|> на интервале MSD = {mean_speed:.6f}")
    print(f"lambda = 3 D / <|v|> = {lambda_from_msd:.6f}")
    if np.count_nonzero(m_par) >= 3:
        print(
            f"парабола (фит [0, {t_par_hi:.1f}]): "
            f"MSD ≈ {c2:.6f} τ² + {c1:.6f} τ + {c0:.6f}"
        )
        print(
            f"  на графике: та же кривая на τ ∈ [0, {t_par_draw_hi:.1f}]"
        )
    print(f"сохранено: {out_path}")
    _print_msd_comparisons(
        slope_msd,
        float(d["D_from_msd"]),
        c2,
        np.count_nonzero(m_par) >= 3,
    )


def _print_msd_comparisons(
    k: float,
    D: float,
    c2: float,
    parabola_ok: bool,
) -> None:
    """k = dMSD/dτ, D = k/6; lambda1 = 3D/√(3T), lambda2 = 3D/√(c2)."""
    T_want = float(read_param("T_want", "1.0"))
    tau_tr = float(read_param("msd_tau_trans", "3.0"))
    v1 = float(np.sqrt(3.0 * T_want))
    v2 = float(np.sqrt(c2)) if parabola_ok and c2 > 0 else float("nan")

    print()
    print("===ВЫЧИСЛЕНИЯ===")
    print(f"  v = sqrt(3T) = {v1:.6f}")
    if parabola_ok and np.isfinite(v2):
        print(f"  v = sqrt(c2) = {v2:.6f}")
    else:
        print("  v = sqrt(c2) =  —")
    print(f"  D = k/6 = {D:.6f}")
    lam1 = 3.0 * D / v1
    print(f"  lambda1 = 3D / sqrt(3T) = {lam1:.6f}")
    if parabola_ok and np.isfinite(v2) and v2 > 0:
        lam2 = 3.0 * D / v2
        print(f"  lambda2 = 3D / sqrt(c2) = {lam2:.6f}")
    v_tau = v1 * tau_tr
    print(f"  v·τ* = sqrt(3T)·τ* = {v_tau:.4f}")


def main() -> None:
    plot_msd()


if __name__ == "__main__":
    main()
