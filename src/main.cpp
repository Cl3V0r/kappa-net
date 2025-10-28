#include "polymer.hpp"
#include <iostream>

int main() {
    int num_beads = 100;
    double kappa = 15.0;
    double beta = 1.0;
    double rotation_angle = 0.5;
    int sweeps = 1000000;

    PolymerSimulation sim(num_beads, kappa, beta, rotation_angle, sweeps);
    sim.run();

    std::cout << "Simulation complete for kappa = " << kappa << std::endl;
    return 0;
}
