#include "forces.h"
#include <cmath>
#include <algorithm>

double lj_potential(double r) {
    double r6 = 1.0 / (r * r * r * r * r * r);
    return 4.0 * (r6 * r6 - r6);
}

double lj_force(double r) {
    double r6 = 1.0 / (r * r * r * r * r * r);
    return 24.0 * (2.0 * r6 * r6 - r6) / r;
}

void compute_accelerations(const std::vector<Vec3>& coords,
                           double cell_size,
                           double r_cut,
                           std::vector<Vec3>& acc) {
    int N = static_cast<int>(coords.size());
    double rc = std::min(r_cut, cell_size / 2.0);
    double rc2 = rc * rc;
    double half_cell = cell_size / 2.0;

    for (auto& a : acc) a = {0.0, 0.0, 0.0};

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = coords[j][0] - coords[i][0];
            double dy = coords[j][1] - coords[i][1];
            double dz = coords[j][2] - coords[i][2];

            if (dx >  half_cell) dx -= cell_size;
            if (dx < -half_cell) dx += cell_size;
            if (dy >  half_cell) dy -= cell_size;
            if (dy < -half_cell) dy += cell_size;
            if (dz >  half_cell) dz -= cell_size;
            if (dz < -half_cell) dz += cell_size;

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rc2) {
                double r = std::sqrt(r2);
                double f = lj_force(r) / r;
                double fx = f * dx;
                double fy = f * dy;
                double fz = f * dz;

                acc[i][0] -= fx;
                acc[i][1] -= fy;
                acc[i][2] -= fz;
                acc[j][0] += fx;
                acc[j][1] += fy;
                acc[j][2] += fz;
            }
        }
    }
}
