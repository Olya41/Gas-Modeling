#include "simulation.h"
#include "verlet.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

std::vector<Vec3> init_lattice(int N, double cell_size) {
    int n = static_cast<int>(std::ceil(std::cbrt(static_cast<double>(N))));
    double spacing = cell_size / n;
    std::vector<Vec3> coords;
    coords.reserve(N);
    for (int ix = 0; ix < n; ++ix)
        for (int iy = 0; iy < n; ++iy)
            for (int iz = 0; iz < n; ++iz) {
                if (static_cast<int>(coords.size()) >= N) break;
                coords.push_back({(ix + 0.5) * spacing,
                                  (iy + 0.5) * spacing,
                                  (iz + 0.5) * spacing});
            }
    return coords;
}

std::vector<Vec3> init_velocities(int N, double T_target) {
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, std::sqrt(T_target));

    std::vector<Vec3> vel(N);
    Vec3 mean = {0.0, 0.0, 0.0};
    for (int i = 0; i < N; ++i) {
        vel[i] = {dist(rng), dist(rng), dist(rng)};
        mean[0] += vel[i][0];
        mean[1] += vel[i][1];
        mean[2] += vel[i][2];
    }
    mean[0] /= N;
    mean[1] /= N;
    mean[2] /= N;
    for (int i = 0; i < N; ++i) {
        vel[i][0] -= mean[0];
        vel[i][1] -= mean[1];
        vel[i][2] -= mean[2];
    }
    return vel;
}

static void write_row(std::ofstream& f, double t,
                       const std::vector<Vec3>& data, int N) {
    f << std::scientific << std::setprecision(14) << t;
    for (int i = 0; i < N; ++i)
        for (int d = 0; d < 3; ++d)
            f << ' ' << data[i][d];
    f << '\n';
}

void run_simulation(const SimParams& p) {
    double cell_size = std::cbrt(static_cast<double>(p.N) / p.rho);

    std::cout << "=== MD Lennard-Jones (C++) ===" << std::endl;
    std::cout << "N = " << p.N << std::endl;
    std::cout << "dt = " << p.dt << std::endl;
    std::cout << "NUM = " << p.NUM << " steps (total time "
              << p.NUM * p.dt << ")" << std::endl;
    std::cout << "rho = " << p.rho << std::endl;
    std::cout << "cell_size = " << cell_size << std::endl;
    std::cout << "r_cut = " << p.r_cut << std::endl;
    std::cout << "energy_log = " << (p.energy_log ? "true" : "false")
              << std::endl << std::endl;

    auto pos = init_lattice(p.N, cell_size);
    auto vel = std::vector<Vec3>(p.N, {0.0, 0.0, 0.0});
    std::vector<Vec3> acc(p.N, {0.0, 0.0, 0.0});
    compute_accelerations(pos, cell_size, p.r_cut, acc);

    std::ofstream f_pos("positions.txt");
    std::ofstream f_vel("velocities.txt");
    std::ofstream f_ene("energies.txt");

    f_pos << std::scientific << std::setprecision(14);
    f_vel << std::scientific << std::setprecision(14);
    f_ene << std::scientific << std::setprecision(14);

    if (p.energy_log)
        f_ene << "# t    Ek    Ep    Etot" << std::endl;

    // step 0
    double t = 0.0;
    write_row(f_pos, t, pos, p.N);
    write_row(f_vel, t, vel, p.N);

    if (p.energy_log) {
        double ek = kinetic_energy(vel);
        double pe = potential_energy(pos, cell_size, p.r_cut);
        double T_kin = 2.0 * ek / (3.0 * p.N);
        f_ene << t << ' ' << ek << ' ' << pe << ' ' << ek + pe << '\n';
        std::printf("Step 0/%d: Ek=%.4f, Ep=%.4f, Etot=%.4f, T=%.4f\n",
                    p.NUM, ek, pe, ek + pe, T_kin);
    }

    for (int step = 1; step <= p.NUM; ++step) {
        velocity_verlet_step(pos, vel, acc, p.dt, cell_size, p.r_cut);
        t = step * p.dt;

        write_row(f_pos, t, pos, p.N);
        write_row(f_vel, t, vel, p.N);

        if (p.energy_log) {
            double ek = kinetic_energy(vel);
            double pe = potential_energy(pos, cell_size, p.r_cut);
            double T_kin = 2.0 * ek / (3.0 * p.N);
            f_ene << t << ' ' << ek << ' ' << pe << ' ' << ek + pe << '\n';
            if (step % 100 == 0)
                std::printf("Step %d/%d: Ek=%.4f, Ep=%.4f, Etot=%.4f, T=%.4f\n",
                            step, p.NUM, ek, pe, ek + pe, T_kin);
        } else {
            if (step % 100 == 0)
                std::printf("Step %d/%d\n", step, p.NUM);
        }
    }

    f_pos.close();
    f_vel.close();
    f_ene.close();

    std::cout << "Saved to positions.txt, velocities.txt, energies.txt"
              << std::endl;
}
