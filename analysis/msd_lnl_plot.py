"""
ln(MSD) vs ln(τ): весь ряд τ после t_ref; на графике кусок с -1 ≤ ln τ ≤ 5.

  python3 analysis/msd_lnl_plot.py

→ output/plots/msd_lnl.png
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from msd_plots import compute_msd
from sim_io import ROOT


def _ln_ln_line(tau: np.ndarray, msd: np.ndarray, mask: np.ndarray) -> tuple[float, float] | None:
    """y = a·x + b в коорд. (ln τ, ln MSD); None если < 2 точек."""
    if int(np.count_nonzero(mask)) < 2:
        return None
    t = np.asarray(tau, dtype=float)[mask]
    m = np.asarray(msd, dtype=float)[mask]
    if np.any(t <= 0) or np.any(m <= 0):
        return None
    x, y = np.log(t), np.log(m)
    a, b = np.polyfit(x, y, 1)
    return (float(a), float(b))


def main() -> None:
    d = compute_msd(ROOT)
    tau = d["tau"]
    msd = d["msd"]

    # весь ряд по τ после t_ref
    tz, mz = tau, msd

    fs_label, fs_title, fs_tick = 14, 15, 12
    ln_lo, ln_hi = -1.0, 5.0

    m_first = (tz > 0) & (tz <= 2.0) & (mz > 0)
    # вторая прямая: 3…5 по оси ln τ (как на рисунке), не «последние 2 с»
    ln_seg_lo, ln_seg_hi = 3.0, 5.0
    lntz = np.full_like(tz, -np.inf, dtype=float)
    pos = tz > 0
    lntz[pos] = np.log(tz[pos])
    m_last = (pos) & (mz > 0) & (lntz >= ln_seg_lo) & (lntz <= ln_seg_hi)
    fit1 = _ln_ln_line(tz, mz, m_first)
    fit2 = _ln_ln_line(tz, mz, m_last)
    x_int, y_int = None, None
    if fit1 is not None and fit2 is not None:
        a1, b1 = fit1
        a2, b2 = fit2
        den = a1 - a2
        if abs(den) > 1e-12:
            x_int = (b2 - b1) / den
            y_int = a1 * x_int + b1

    fig, ax = plt.subplots(figsize=(9, 6))
    okd = (tz > 0) & (mz > 0)
    if np.count_nonzero(okd) > 0:
        lx = np.log(tz[okd])
        ly = np.log(mz[okd])
        w = (lx >= ln_lo) & (lx <= ln_hi)
        if np.count_nonzero(w) > 0:
            ax.plot(
                lx[w],
                ly[w],
                color="C0",
                lw=2.0,
                label="MSD (данные)",
            )
    x_grid = np.linspace(ln_lo, ln_hi, 200)
    if fit1 is not None:
        a1, b1 = fit1
        ax.plot(
            x_grid,
            a1 * x_grid + b1,
            "C1--",
            lw=1.5,
            label=r"прям. по $0<\tau\leq 2$ с",
        )
    if fit2 is not None:
        a2, b2 = fit2
        ax.plot(
            x_grid,
            a2 * x_grid + b2,
            "C3--",
            lw=1.5,
            label=rf"прям. по $3 \leq \ln\tau \leq 5$ (шкала графика)",
        )
    if x_int is not None and y_int is not None:
        tau_star = float(np.exp(x_int))
        ax.plot(
            [x_int],
            [y_int],
            "k*",
            ms=12,
            zorder=5,
            label=rf"пересеч., $\tau^* \approx {tau_star:.3f}\,\mathrm{{с}}$",
        )
    ax.set_xlabel(r"$\ln\tau$", fontsize=fs_label)
    ax.set_ylabel(r"$\ln\,\mathrm{MSD}$", fontsize=fs_label)
    ax.set_title(
        r"$\ln\,\mathrm{MSD}$ — $\ln\tau$",
        fontsize=fs_title,
    )
    ax.set_xlim(ln_lo, ln_hi)
    ax.tick_params(axis="both", labelsize=fs_tick)
    ax.grid(True, alpha=0.35)
    ax.legend(loc="lower right", fontsize=10)
    fig.tight_layout()

    out = ROOT / "output/plots/msd_lnl.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"сохранено: {out}")
    if x_int is not None and y_int is not None:
        print("--- прямые: ln MSD = a·ln τ + b ---")
        if fit1 is not None:
            print(f"  0<τ≤2 (с):         a1={fit1[0]:.6f},  b1={fit1[1]:.6f}")
        if fit2 is not None:
            print(
                f"  3≤ln τ≤5 (по рис.):  a2={fit2[0]:.6f},  b2={fit2[1]:.6f}"
            )
        print(
            f"  пересечение:  ln τ = {x_int:.6f},  ln MSD = {y_int:.6f}  |  "
            f"τ = {np.exp(x_int):.6f} с,  MSD = {np.exp(y_int):.6f}"
        )
    elif fit1 is not None and fit2 is not None:
        print("--- прямые параллельны, пересечения нет ---")
    else:
        print(
            "--- нет пары прямых (нужны ≥2 точек при 0<τ≤2 с и в 3≤lnτ≤5) ---"
        )


if __name__ == "__main__":
    main()
