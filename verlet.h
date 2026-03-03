#pragma once
#include "forces.h"

void velocity_verlet_step(std::vector<Vec3>& pos,
                          std::vector<Vec3>& vel,
                          std::vector<Vec3>& acc,
                          double dt,
                          double cell_size,
                          double r_cut);

double kinetic_energy(const std::vector<Vec3>& vel);

double potential_energy(const std::vector<Vec3>& pos,
                        double cell_size,
                        double r_cut);
