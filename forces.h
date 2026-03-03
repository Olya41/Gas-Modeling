#pragma once
#include <vector>
#include <array>

using Vec3 = std::array<double, 3>;

double lj_potential(double r);
double lj_force(double r);

void compute_accelerations(const std::vector<Vec3>& coords,
                           double cell_size,
                           double r_cut,
                           std::vector<Vec3>& acc);
