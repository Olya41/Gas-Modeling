#include "simulation.h"
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <cmath>

static SimParams load_params(const std::string& filename = "params") {
    SimParams p;
    std::ifstream f(filename);
    std::string line;
    std::map<std::string, std::string> kv;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto pos = line.find('=');
        if (pos == std::string::npos) continue;
        auto trim = [](std::string s) {
            s.erase(0, s.find_first_not_of(" \t"));
            auto end = s.find_last_not_of(" \t");
            return end == std::string::npos ? "" : s.substr(0, end + 1);
        };
        kv[trim(line.substr(0, pos))] = trim(line.substr(pos + 1));
    }
    if (kv.count("N"))          p.N          = std::stoi(kv["N"]);
    if (kv.count("NUM"))        p.NUM        = std::stoi(kv["NUM"]);
    if (kv.count("dt"))         p.dt         = std::stod(kv["dt"]);
    if (kv.count("rho"))        p.rho        = std::stod(kv["rho"]);
    if (kv.count("r_cut"))      p.r_cut      = std::stod(kv["r_cut"]);
    if (kv.count("energy_log")) p.energy_log = (kv["energy_log"] == "true");
    if (kv.count("T_want"))     p.T_want     = std::stod(kv["T_want"]);
    if (kv.count("renew_freq")) p.renew_freq = std::stoi(kv["renew_freq"]);
    if (kv.count("crit_num"))   p.crit_num   = std::stoi(kv["crit_num"]);
    return p;
}

int main() {
    SimParams p = load_params();
    double cell_size = std::cbrt((double)p.N / p.rho);

    std::cout << "=== MD Lennard-Jones (C++) ===" << std::endl;
    std::cout << "N          = " << p.N << std::endl;
    std::cout << "dt         = " << p.dt << std::endl;
    std::cout << "NUM        = " << p.NUM
              << "  (total time " << p.NUM * p.dt << ")" << std::endl;
    std::cout << "rho        = " << p.rho << std::endl;
    std::cout << "cell_size  = " << cell_size << std::endl;
    std::cout << "r_cut      = " << p.r_cut << std::endl;
    std::cout << "energy_log = " << std::boolalpha << p.energy_log << std::endl;
    std::cout << "T_want     = " << p.T_want << std::endl;
    std::cout << "renew_freq = " << p.renew_freq << std::endl;
    std::cout << "crit_num   = " << p.crit_num << std::endl;
    std::cout << std::endl;

    run_simulation(p);
    return 0;
}
