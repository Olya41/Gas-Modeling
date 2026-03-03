import numpy as np
import matplotlib.pyplot as plt

# --- Чтение positions.txt ---
data = np.loadtxt("positions.txt")
times = data[:, 0]
N = (data.shape[1] - 1) // 3
positions_wrapped = data[:, 1:].reshape(len(times), N, 3)
cell_size = np.max(positions_wrapped[0])  # позиции лежат в [0, cell_size)

# --- "Распрямление" траекторий (unwrap periodic boundaries) ---
# Если между соседними шагами координата прыгнула на > cell_size/2,
# значит частица перешла через границу ячейки.
# dx > cell_size/2  => частица ушла через левую границу (прыгок ~0 -> ~cell_size) => вычитаем cell_size
# dx < -cell_size/2 => частица ушла через правую границу (прыжок ~cell_size -> ~0) => прибавляем cell_size

positions_unwrapped = positions_wrapped.copy()

for t in range(1, len(times)):
    delta = positions_wrapped[t] - positions_wrapped[t - 1]
    # Корректируем скачки
    delta -= cell_size * np.round(delta / cell_size)
    positions_unwrapped[t] = positions_unwrapped[t - 1] + delta

# --- Сохранение unlimited_positions.txt ---
unlimited_data = np.column_stack([times, positions_unwrapped.reshape(len(times), N * 3)])
np.savetxt("unlimited_positions.txt", unlimited_data)
print(f"Saved unlimited_positions.txt ({unlimited_data.shape[0]} rows, {unlimited_data.shape[1]} cols)")

# --- MSD (среднеквадратичное смещение) для оценки коэффициента диффузии ---
displacement = positions_unwrapped - positions_unwrapped[0]  # (steps, N, 3)
msd = np.mean(np.sum(displacement**2, axis=2), axis=1)        # <|r(t)-r(0)|^2>

# --- Оценка коэффициента диффузии D из MSD = 6Dt ---
# Берём последние 3/4 измерений
fit_start = len(times) // 4
slope, intercept = np.polyfit(times[fit_start:], msd[fit_start:], 1)

D = slope / 6

# --- Средняя длина свободного пробега ---
# l = 3D / <v>, средняя скорость берётся за тот же интервал
vel_data = np.loadtxt("velocities.txt")
velocities = vel_data[:, 1:].reshape(len(times), N, 3)
speeds = np.sqrt(np.sum(velocities**2, axis=2))  # (steps, N)
mean_speed = np.mean(speeds[fit_start:])
mean_free_path = 3 * D / mean_speed

# --- График ---
fit_times = times[fit_start:]
plt.figure(figsize=(8, 5))
plt.plot(times, msd, label="MSD")
plt.plot(fit_times, slope * fit_times + intercept, '--', label=f"fit: 6Dt, D={D:.4f}")
plt.xlabel("t")
plt.ylabel(r"$\langle |r(t) - r(0)|^2 \rangle$")
plt.title("MSD (среднеквадратичное смещение)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("mean_free_path.png", dpi=150)

# --- Вывод ---
print(f"Коэффициент наклона прямой: {slope:.4f}")
print(f"Коэффициент диффузии D = coeff/6 = {D:.4f}")
print(f"Средняя скорость v = {mean_speed:.4f}")
print(f"Средняя длина свободного пробега λ = 3D/v = {mean_free_path:.4f}")
