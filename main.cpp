#include "simulation.h"

int main() {
    SimParams p;
    p.N          = 100;
    p.dt         = 0.01;
    p.NUM        = 2000;
    p.rho        = 0.3;
    p.r_cut      = 2.5;
    p.energy_log = true;

    run_simulation(p);
    return 0;
}
