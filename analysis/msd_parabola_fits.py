"""
Отдельный график: MSD и три квадратичных фита на τ ∈ [0,1], [0,2], [0,3] (единицы как в t_ref).

  python3 analysis/msd_parabola_fits.py

Сохраняет output/plots/msd_parabola_fits.png
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from msd_plots import compute_msd
from sim_io import ROOT

# верхние границы окон фита, с
FITS = (1.0, 2.0, 3.0)
COLORS = ("C0", "C1", "C2")
# до куда рисовать ось и экстраполяцию парабол
TAU_PLOT_MAX = 4.5


def main() -> None:
    d = compute_msd(ROOT)
    tau = d["tau"]
    msd = d["msd"]
    t_ref = d["t_msd_ref"]
    tau_max = float(tau[-1])

    # ось: до 4.5 с (или по данным, если прогон короче)
    x_max = min(TAU_PLOT_MAX, max(tau_max, 0.1))
    # данные на графике — только в этом окне по τ
    m_show = tau <= x_max
    tw, msw = tau[m_show], msd[m_show]

    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(tw, msw, color="0.2", lw=2.0, label="MSD (данные)", zorder=1)

    for t_hi, color in zip(FITS, COLORS):
        if t_hi > tau_max + 1e-12:
            print(f"[0, {t_hi:.0f}] с: нет τ до этой границы, фит пропущен")
            continue
        m = (tau >= 0.0) & (tau <= t_hi)
        npts = int(np.count_nonzero(m))
        if npts < 3:
            print(f"[0, {t_hi:.0f}] с: слишком мало точек ({npts}), пропуск")
            continue
        c2, c1, c0 = np.polyfit(tau[m], msd[m], 2)
        tau_line = np.linspace(0.0, min(TAU_PLOT_MAX, tau_max), 300)
        y_line = ((c2 * tau_line) + c1) * tau_line + c0
        ax.plot(
            tau_line,
            y_line,
            color=color,
            lw=2.0,
            label=rf"фит $[0,\,{t_hi:.0f}]$ с: ${c2:.3f}\tau^2+{c1:.3f}\tau+{c0:.2f}$",
        )
        print(
            f"фит [0, {t_hi:.0f}] с: MSD ≈ {c2:.6f} τ² + {c1:.6f} τ + {c0:.6f}  "
            f"({npts} точек)"
        )

    ax.set_xlim(0.0, x_max)
    ax.set_xlabel(r"$\tau = t - t_{\mathrm{ref}}$ (с)", fontsize=13)
    ax.set_ylabel(
        r"$\langle |r(t_{\mathrm{ref}}+\tau) - r(t_{\mathrm{ref}})|^2 \rangle$",
        fontsize=13,
    )
    ax.set_title(
        rf"Параболы по данным на $[0,1]$, $[0,2]$, $[0,3]$ с;  "
        r"$t_{\mathrm{ref}}$" + f" = {t_ref:.4f} с",
        fontsize=14,
    )
    ax.grid(True, alpha=0.35)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()

    out = ROOT / "output/plots/msd_parabola_fits.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"сохранено: {out}")


if __name__ == "__main__":
    main()
