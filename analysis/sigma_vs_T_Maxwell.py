"""
r_min vs T через Максвелловское распределение скоростей.

Для каждой T сэмплируем N_SAMPLES пар частиц с независимыми
максвелловскими скоростями (каждая компонента ~ N(0, sqrt(T))),
вычисляем энергию в ЦОМ E = |v_rel|² / 4  (лобовой удар, m=1, μ=½),
находим r_min из V(r_min) = E аналитически и усредняем.

Среднее <|v_rel|²> = 6T (складываются дисперсии двух частиц по 3 осям),
поэтому <E> = 3T/2 — это больше, чем модельное 3T/4 при v_rel = sqrt(3T).

  python3 analysis/sigma_vs_T_Maxwell.py

→ output/plots/sigma_vs_T_Maxwell.png
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from count_sigma import r_lj_turn_from_E, r_V_equals_three_fourths_T
from sim_io import ROOT, read_param

N_SAMPLES = 30_000
N_T = 80
T_LO, T_HI = 0.05, 5.0


def main() -> None:
    T_want = float(read_param("T_want", "1.0"))
    T_arr = np.linspace(T_LO, T_HI, N_T)
    rng = np.random.default_rng(42)

    r_mean = np.empty(N_T)
    r_std = np.empty(N_T)

    for k, T in enumerate(T_arr):
        # Максвелл: каждая компонента ~ N(0, sqrt(T)), m=1
        v1 = rng.normal(0.0, np.sqrt(T), (N_SAMPLES, 3))
        v2 = rng.normal(0.0, np.sqrt(T), (N_SAMPLES, 3))
        v_rel2 = np.sum((v2 - v1) ** 2, axis=1)   # |v_rel|²; <v_rel²> = 6T
        E = 0.25 * v_rel2                           # E = μ v_rel²/2, μ=1/2
        r_vals = r_lj_turn_from_E(E)                # аналитика, векторно
        valid = np.isfinite(r_vals)
        r_mean[k] = np.mean(r_vals[valid]) if np.any(valid) else np.nan
        r_std[k] = np.std(r_vals[valid]) if np.any(valid) else np.nan

    r_model = r_V_equals_three_fourths_T(T_arr)   # старая кривая, E = 3T/4

    # r_min(T) при среднем <E> = 3T/2 как дополнительный ориентир
    r_mean_E = r_lj_turn_from_E(1.5 * T_arr)

    # --- график ---
    fig, ax = plt.subplots(figsize=(9, 6))

    ax.fill_between(
        T_arr,
        r_mean - r_std,
        r_mean + r_std,
        alpha=0.18,
        color="C0",
        label=r"$\langle r_{\min}\rangle \pm \sigma$ (Максвелл)",
    )
    ax.plot(
        T_arr,
        r_mean,
        "C0-",
        lw=2.2,
        label=rf"$\langle r_{{\min}}\rangle$ (Максвелл, $N={N_SAMPLES}$ пар на точку)",
    )
    ax.plot(
        T_arr,
        r_model,
        "C1--",
        lw=1.8,
        label=(
            r"$v_{\mathrm{rel}}=\sqrt{3T}$, $E=\frac{3T}{4}$:  "
            r"$r_{\min}=\left(\frac{2}{1+\sqrt{1+3T/4}}\right)^{1/6}$"
        ),
    )
    ax.plot(
        T_arr,
        r_mean_E,
        "C2:",
        lw=2.0,
        label=r"$E=\langle E\rangle=\frac{3T}{2}$:  $r_{\min}=\left(\frac{2}{1+\sqrt{1+3T/2}}\right)^{1/6}$",
    )

    r_at_Twant = float(r_V_equals_three_fourths_T(np.asarray([T_want]))[0])
    ax.scatter(
        [T_want],
        [r_at_Twant],
        s=50,
        c="magenta",
        zorder=5,
        edgecolors="k",
        linewidths=0.5,
        label=rf"params: $T={T_want}$,  $r_{{\min}}={r_at_Twant:.4f}$",
    )

    ax.set_xlabel(r"$T$", fontsize=14)
    ax.set_ylabel(r"$r_{\min}$ (в ед. $\sigma$)", fontsize=14)
    ax.set_title(r"$r_{\min}$ vs $T$: Максвелловское распределение скоростей", fontsize=14)
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=9, loc="upper right")
    fig.tight_layout()

    out = ROOT / "output/plots/sigma_vs_T_Maxwell.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"сохранено: {out}")
    print(f"<E>/T при Максвелл: 3/2 = {1.5:.4f}  (против 3/4 = {0.75:.4f} в модели)")


if __name__ == "__main__":
    main()
