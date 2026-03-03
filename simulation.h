#pragma once
#include "forces.h"
#include <string>

struct SimParams {
    int    N          = 200;
    int    NUM        = 4000;
    double dt         = 0.01;
    double rho        = 0.3;
    double r_cut      = 2.5;
    bool   energy_log = true;
};

std::vector<Vec3> init_lattice(int N, double cell_size);
std::vector<Vec3> init_velocities(int N, double T_target = 1.0);

void run_simulation(const SimParams& p);
