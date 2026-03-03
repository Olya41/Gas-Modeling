#include "verlet.h"
#include <cmath>
#include <algorithm>

void velocity_verlet_step(std::vector<Vec3>& pos,
                          std::vector<Vec3>& vel,
                          std::vector<Vec3>& acc,
                          double dt,
                          double cell_size,
                          double r_cut) {
    int N = static_cast<int>(pos.size());
    double half_dt2 = 0.5 * dt * dt;
    double half_dt  = 0.5 * dt;

    // positions update: r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt^2
    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) {
            pos[i][d] += vel[i][d] * dt + acc[i][d] * half_dt2;
        }
    }

    // save old accelerations, compute new
    std::vector<Vec3> acc_old = acc;
    compute_accelerations(pos, cell_size, r_cut, acc);

    // velocities update: v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt
    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) {
            vel[i][d] += half_dt * (acc_old[i][d] + acc[i][d]);
        }
    }

    // wrap positions into [0, cell_size)
    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) {
            if (pos[i][d] >= cell_size) pos[i][d] -= cell_size;
            if (pos[i][d] < 0.0)       pos[i][d] += cell_size;
        }
    }
}

double kinetic_energy(const std::vector<Vec3>& vel) {
    double ek = 0.0;
    for (const auto& v : vel) {
        ek += v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    }
    return 0.5 * ek;
}

double potential_energy(const std::vector<Vec3>& pos,
                        double cell_size,
                        double r_cut) {
    int N = static_cast<int>(pos.size());
    double rc = std::min(r_cut, cell_size / 2.0);
    double rc2 = rc * rc;
    double half_cell = cell_size / 2.0;
    double u_shift = lj_potential(rc);
    double pe = 0.0;

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = pos[j][0] - pos[i][0];
            double dy = pos[j][1] - pos[i][1];
            double dz = pos[j][2] - pos[i][2];

            if (dx >  half_cell) dx -= cell_size;
            if (dx < -half_cell) dx += cell_size;
            if (dy >  half_cell) dy -= cell_size;
            if (dy < -half_cell) dy += cell_size;
            if (dz >  half_cell) dz -= cell_size;
            if (dz < -half_cell) dz += cell_size;

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < rc2) {
                double r = std::sqrt(r2);
                pe += lj_potential(r) - u_shift;
            }
        }
    }
    return pe;
}
