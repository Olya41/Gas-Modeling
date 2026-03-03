#!/bin/bash
set -e

echo "=== Сборка ==="
cmake -B build && cmake --build build

echo ""
echo "=== Симуляция (C++) ==="
./build/md_simulation

echo ""
echo "=== График энергий (plot_energy.py) ==="
python3 plot_energy.py

echo ""
echo "=== Распределение скоростей (check.py) ==="
python3 check.py

echo ""
echo "=== MSD и длина свободного пробега (free_run.py) ==="
python3 free_run.py

echo ""
echo "=== Траектория XYZ (plot_animation.py) ==="
python3 plot_animation.py

echo ""
echo "=== Готово ==="
echo "Выходные файлы:"
echo "  positions.txt, velocities.txt, energies.txt  — сырые данные"
echo "  unlimited_positions.txt                       — unwrapped траектории"
echo "  energy.png                                    — график энергий"
echo "  vel_gist.png, vel_linear.png                  — распределение скоростей"
echo "  mean_free_path.png                            — MSD и диффузия"
echo "  trajectory.xyz                                — траектория для визуализации"
