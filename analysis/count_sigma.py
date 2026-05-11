"""
Модельный лобовой «удар» пары в ЛД: для каждой T — характерная относительная
скорость v_rel = sqrt(3T) (как 3D RMS при m=1), в СОМ E = ½μ v_rel², μ=½.
Классически: в ближайшей остановке E = V(r_min) на ветви 0<r<1 (LJ «стенка»);
численно: корень 4(r^{-12}−r^{-6}) = E, даёт r_min(T).

На графике только связь V=3T/4 и v_rel=sqrt(3T) в ЦОМ; маркер — T_want из params.

  python3 analysis/count_sigma.py

→ output/plots/sigma_vs_T.png
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from sim_io import ROOT, read_param


def _V_lj(r: float) -> float:
    r = float(r)
    if r <= 0:
        return float("inf")
    inv2 = 1.0 / (r * r)
    inv6 = inv2 * inv2 * inv2
    return 4.0 * (inv6 * inv6 - inv6)


def r_min_from_energy(e_kin: float) -> float:
    """
    Для 1D лобового LJ в COM: E = ½ μ v_rel² = V(r_min) в момент остановки (радиальная
    скорость 0) на **внутренней** ветви V(r) > 0, r < σ.
    """
    e_kin = float(e_kin)
    if e_kin <= 0:
        return float("nan")

    def f(r: float) -> float:
        return _V_lj(r) - e_kin

    r_lo, r_hi = 0.11, 0.9999
    fl, fh = f(r_lo), f(r_hi)
    if fl * fh > 0:
        if fh > 0:
            while r_lo > 1.01e-3 and f(r_lo) * f(r_hi) > 0:
                r_lo *= 0.7
        else:
            return float("nan")
        fl, fh = f(r_lo), f(r_hi)
    if fl * fh > 0:
        return float("nan")
    for _ in range(80):
        m = 0.5 * (r_lo + r_hi)
        if f(m) * f(r_lo) > 0:
            r_lo = m
        else:
            r_hi = m
    return float(0.5 * (r_lo + r_hi))


def v_from_T(T: float) -> float:
    """m=1, характерная 3D-скорость: v_rms = sqrt(3T)."""
    return float(np.sqrt(3.0 * max(T, 1e-12)))


def e_kin_com(v_rel: float) -> float:
    """E = ½ μ v²,  μ=1/2,  m₁=m₂=1."""
    return 0.25 * (v_rel * v_rel)


def r_lj_turn_from_E(E_arr: np.ndarray) -> np.ndarray:
    """
    Решение 4(u²−u)=E при u=r^{−6} (поворот, r<1): r^6 = 2/(1+√(1+E)).

    При E→0+ внутреннее стремление к r→1 слева через r^6 = 2(√(1+E)−1)/E.
    """
    E = np.asarray(E_arr, dtype=float)
    inner = np.where(E > 0.0, 2.0 * (np.sqrt(1.0 + E) - 1.0) / E, np.nan)
    return np.power(inner, 1.0 / 6.0)


def r_V_equals_three_fourths_T(T_arr: np.ndarray) -> np.ndarray:
    """4/r¹² − 4/r⁶ = 3T/4 ⇒ E = 3T/4 (как E_кин в модели v_rel = √(3T), μ = ½)."""
    return r_lj_turn_from_E(0.75 * np.asarray(T_arr, dtype=float))


def main() -> None:
    T_want = float(read_param("T_want", "1.0"))

    n = 100
    T_lo, T_hi = 0.05, 5.0
    T_list = np.linspace(T_lo, T_hi, n)

    r_curve = r_V_equals_three_fourths_T(T_list)
    r_at_Twant = float(r_V_equals_three_fourths_T(np.asarray([T_want]))[0])

    fs = 14
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(
        T_list,
        r_curve,
        "C0-",
        lw=2.2,
        label=(
            r"$\frac{4}{r^{12}}-\frac{4}{r^{6}}=\frac{3T}{4}$"
            "\n"
            r"$r_{\min}=\left(\frac{8\left(\sqrt{1+\frac{3T}{4}}-1\right)}{3T}\right)^{1/6}$"
        ),
    )
    ax.scatter(
        [T_want],
        [r_at_Twant],
        s=45,
        c="magenta",
        zorder=5,
        edgecolors="k",
        linewidths=0.45,
        label=rf"params: $T = {T_want}$,  $r_{{\min}} = {r_at_Twant:.4f}$",
    )
    ax.set_xlabel(r"$T$", fontsize=fs)
    ax.set_ylabel(r"$r_{\min}$ (в ед. $\sigma$)", fontsize=fs)
    ax.set_title(r"$r_{\min}$ vs $T$", fontsize=14)
    ax.grid(True, alpha=0.4)
    ax.legend(loc="best", fontsize=11, alignment="left")
    fig.tight_layout()
    out = ROOT / "output/plots/sigma_vs_T.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"сохранено: {out}")
    print(
        f"T_want = {T_want}  →  v_rel = {v_from_T(T_want):.6f}  →  "
        f"r_min = {r_at_Twant:.6f}"
    )


if __name__ == "__main__":
    main()
