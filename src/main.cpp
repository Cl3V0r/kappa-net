#include "polymer.hpp"
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

int main() {
    // ------------------------------
    // Simulation parameters
    // ------------------------------
    int num_beads = 100;
    double beta = 1.0;
    double rotation_angle = 0.5;
    int sweeps = 1000000;
    unsigned long seed = 42;

    // κ values to sweep over
    std::vector<double> kappa_values;
    for (double k = 5.0; k <= 20.0; k += 1.0)
        kappa_values.push_back(k);

    // ------------------------------
    // Ensure output folder exists
    // ------------------------------
    fs::path out_dir = "data/raw/";
    fs::create_directories(out_dir);

    // ------------------------------
    // Loop over kappa
    // ------------------------------
    for (double kappa : kappa_values) {
        std::cout << "Running simulation for kappa = " << kappa << " ..." << std::endl;

        PolymerSimulation sim(num_beads, kappa, beta, rotation_angle, sweeps, seed);
        sim.run();

        std::string filename = out_dir.string() + "config_kappa_" + std::to_string(int(kappa)) + ".csv";
        sim.save_configuration(filename, sweeps);

        std::cout << "Saved configuration to " << filename << std::endl;
    }

    std::cout << "All simulations complete!" << std::endl;
    return 0;
}
