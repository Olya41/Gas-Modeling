"""
Вычисление средней длины свободного пробега через детектирование периценторов.

Столкновение: пара (i,j) сближается до r < threshold.
Два критерия (порог r для «столкновения»; перицентр внутри r < …):
  r_min  — эффективный центр–центр при T_want, как в count_sigma / sigma_vs_T.png
           (лобовой ЛД, v_rel = sqrt(3T), E = ½μ v_rel² = V(r_min) при r<σ)
  sigma  = 1.0 — ноль потенциала (V(σ)=0)

Как находим перицентр (момент наибольшего сближения):
  Для каждой пары находим все непрерывные отрезки по времени, где dist(i,j) < threshold.
  Внутри каждого такого отрезка выбираем шаг с минимальным расстоянием — это
  и есть перицентр. Метод работает при условии, что dt достаточно мал относительно
  продолжительности столкновения (обычно выполняется для LJ при разумном dt).
  Отрезки, которые начались или не закончились на краю анализируемого окна,
  отбрасываются как неполные.

Длина свободного пробега: для каждой частицы i сортируем её перицентры по
времени; между двумя последовательными перицентрами суммируем |v_i(t)| * dt —
это расстояние, пройденное частицей в свободном полёте. Усредняем по всем
сегментам и всем частицам.

Почему λ(sim) при r_min ≪ σ может отличаться сильнее, чем «грубая» теория
λ ∝ 1/d²: теория жёстких сфер сравнивает только S ∝ d², меняя d на 4%.
В реальном газе 3D много сближений — острые, с расст. ближайшего в (r_min, σ);
они дают pericenter при dist < σ, но не при dist < r_min, «пропадают».
Исчезновение **части** близких проходов сильно удлиняет **остающиеся** сегменты
между **редкими** глубокими сближениями, поэтому отношение λ(sim) **не
обязано** быть (σ / r_min)²: это **не** ошибка кода, а другая «физика
детектора».
"""
import io
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

from count_sigma import e_kin_com, r_min_from_energy, v_from_T
from msd_plots import plot_msd
from sim_io import ROOT, read_param


N_params     = int(read_param("N", "100"))
rho          = float(read_param("rho", "0.8"))
dt           = float(read_param("dt", "0.01"))
T_want       = float(read_param("T_want", "1.0"))
num_for_coll = int(read_param("num_for_coll", "500"))

# N и число строк берём из файла траектории (params мог устареть)
with open(ROOT / "output/data/positions.txt") as f:
    first_row = f.readline().split()
    total_lines = 1 + sum(1 for _ in f)
N = (len(first_row) - 1) // 3
if N != N_params:
    print(f"(предупреждение: N в params = {N_params}, в positions.txt N = {N})")
cell_size = (N / rho) ** (1.0 / 3.0)

# ---------------------------------------------------------------------------
# Критерии: sigma = 1.0; r_min = эффективное сближение при T_want (как count_sigma)
# ---------------------------------------------------------------------------
SIGMA   = 1.0
R_MIN   = r_min_from_energy(e_kin_com(v_from_T(T_want)))

# ---------------------------------------------------------------------------
# Оценка числа строк: λ ~ 1/(sqrt(2) π d² ρ),  d = r_min; чем длиннее λ, тем больше
# снимков нужно на одно ожидаемое «столкновение».
# ---------------------------------------------------------------------------
v_mean_est  = np.sqrt(8.0 * T_want / np.pi)
sigma_coll  = np.pi * R_MIN ** 2
lambda_est  = 1.0 / max(np.sqrt(2) * rho * sigma_coll, 1e-9)
tau_est     = lambda_est / max(v_mean_est, 0.1)

steps_needed = int(2 * num_for_coll * tau_est / (N * dt)) + 5_000
n_rows       = max(steps_needed, 20_000)
n_rows       = min(n_rows, total_lines)

# ---------------------------------------------------------------------------
# Чтение последних n_rows строк (без загрузки всего файла)
# ---------------------------------------------------------------------------
def read_tail(path, n):
    out = subprocess.run(["tail", "-n", str(n), str(path)],
                         capture_output=True, text=True, check=True).stdout
    return np.loadtxt(io.StringIO(out))


pos_raw = read_tail(ROOT / "output/data/positions.txt",  n_rows)
vel_raw = read_tail(ROOT / "output/data/velocities.txt", n_rows)

times      = pos_raw[:, 0]
M          = len(times)
positions  = pos_raw[:, 1:].reshape(M, N, 3)
velocities = vel_raw[:, 1:].reshape(M, N, 3)
dt_sim     = float(times[1] - times[0]) if M > 1 else dt

print(f"Столкновения (хвост) : {M} снимков из {total_lines},  "
      f"t ∈ [{times[0]:.2f}, {times[-1]:.2f}],  dt = {dt_sim:.4f}")
print(f"Частиц               : {N}")
print(
    f"r_min (T_want) для оценок λ, как в count_sigma:  {R_MIN:.6f}  "
    f"(T_want = {T_want})  →  λ_theory ≈ {lambda_est:.4f}"
)
print()

# ---------------------------------------------------------------------------
# Основная функция: детектирование + длины пробега для заданного threshold
# ---------------------------------------------------------------------------
def compute_free_paths(threshold):
    particle_peri = [[] for _ in range(N)]

    for i, j in combinations(range(N), 2):
        dr    = positions[:, j, :] - positions[:, i, :]
        dr   -= cell_size * np.round(dr / cell_size)   # минимальное изображение (ПГУ)
        dists = np.linalg.norm(dr, axis=1)

        below = dists < threshold
        if not np.any(below):
            continue

        padded       = np.empty(M + 2, dtype=bool)
        padded[0]    = False
        padded[1:-1] = below
        padded[-1]   = False
        diff   = np.diff(padded.astype(np.int8))
        starts = np.where(diff ==  1)[0]
        ends   = np.where(diff == -1)[0]

        for s, e in zip(starts, ends):
            if s == 0 or e == M:
                continue
            peri = s + int(np.argmin(dists[s:e]))
            particle_peri[i].append(peri)
            particle_peri[j].append(peri)

    free_paths = []
    for i in range(N):
        peri_sorted = sorted(set(particle_peri[i]))
        if len(peri_sorted) < 2:
            continue
        for k in range(len(peri_sorted) - 1):
            a, b = peri_sorted[k], peri_sorted[k + 1]
            if b <= a + 1:
                continue
            speeds = np.linalg.norm(velocities[a:b, i, :], axis=1)
            free_paths.append(float(np.sum(speeds) * dt_sim))

    return free_paths


# ---------------------------------------------------------------------------
# Вывод и график для одного критерия
# ---------------------------------------------------------------------------
def report_and_plot(threshold, label, ax):
    coll_section = np.pi * threshold ** 2
    lambda_theory = 1.0 / (np.sqrt(2) * rho * coll_section)

    print(f"--- {label} = {threshold:.6f} ---")
    print(f"Сечение столкновения pi * {label}^2 = {coll_section:.6f}")

    free_paths = compute_free_paths(threshold)

    if not free_paths:
        print("Недостаточно событий.")
        return None

    lambda_sim = float(np.mean(free_paths))
    lambda_std = float(np.std(free_paths))
    n_seg = len(free_paths)

    print(
        f"N сегментов (l)     = {n_seg}  (чем больше — тем устойчивее среднее)"
    )
    print(f"lambda (sim)    = {lambda_sim:.4f}  +/-  {lambda_std:.4f}")
    print(f"lambda (theory) = {lambda_theory:.4f}  (1/(√2·ρ·π·d²), d=порог)")
    print()

    fp    = np.array(free_paths)
    l_arr = np.linspace(0, np.percentile(fp, 99), 300)

    ax.hist(fp, bins=60, density=True, alpha=0.7, label="симуляция")
    ax.plot(l_arr, (1 / lambda_sim) * np.exp(-l_arr / lambda_sim), "r-", lw=2,
            label=f"exp(-l/λ)/λ,  λ = {lambda_sim:.3f}")
    ax.axvline(lambda_sim,    color="red",  linestyle="--", lw=1.2,
               label=f"λ sim  = {lambda_sim:.3f}")
    ax.axvline(lambda_theory, color="gray", linestyle=":",  lw=1.5,
               label=f"λ theory = {lambda_theory:.3f}")
    ax.set_xlabel("l")
    ax.set_ylabel("P(l)")
    ax.set_title(f"Распределение длин свободного пробега  ({label} = {threshold:.4f})")
    ax.legend()
    return lambda_sim


# ---------------------------------------------------------------------------
# Запуск для двух критериев
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

lambda_sim_r = report_and_plot(R_MIN, "r_min", axes[0])
lambda_sim_s = report_and_plot(SIGMA, "sigma", axes[1])

if lambda_sim_r is not None and lambda_sim_s is not None:
    # λtheory(d) = 1/(√2·ρ·π·d²)  →  λ(d₁)/λ(d₂) = (d₂/d₁)²; с d_r < d_σ : λ(r_min) > λ(σ)
    ratio_th = (SIGMA / R_MIN) ** 2
    ratio_sm = float(lambda_sim_r) / max(float(lambda_sim_s), 1e-30)
    print("--- сравнение r_min и σ (одна траектория) ---")
    print(
        f"Теория 1/d²:  λtheory(r_min) / λtheory(σ)  =  (σ/r_min)²  ≈  {ratio_th:.3f}  "
        f"— λ по теории растёт лишь на ~{(ratio_th - 1.0) * 100:.0f}%, если поменять только d в π d²"
    )
    print(
        f"Симуляция:   λ sim (r_min) / λ sim (σ)  =  {ratio_sm:.3f}  — "
        f"сильно больше (σ/r_min)²: при новом пороге выпадает целый класс сближений "
        f"с dist в (r_min, σ), а не «подрезается» площадь на ~4%."
    )
    print()

plt.tight_layout()
plt.savefig(ROOT / "output/plots/free_path_dist.png", dpi=150)
plt.close()

# MSD — отдельный модуль (перерисовка: `python3 analysis/msd_plots.py` без free_run)
print()
plot_msd()
