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

    // For ML: increase sweeps to get many configurations
    int sweeps = 5000000;       // total MC sweeps
    int save_every = 1000;      // save every 1000 sweeps

    unsigned long seed = 42;

    // Sweep over multiple kappa values
    std::vector<double> kappa_values;
    for (double k = 5.0; k <= 20.0; k += 0.1)
        kappa_values.push_back(k);

    // Ensure output folder exists
    fs::path out_dir = "data/raw/";
    fs::create_directories(out_dir);

    // ------------------------------
    // Loop over kappa values
    // ------------------------------
    for (double kappa : kappa_values) {
        std::cout << "Running simulation for kappa = " << kappa << " ..." << std::endl;

        PolymerSimulation sim(num_beads, kappa, beta, rotation_angle, sweeps, seed);

        // Provide save_every to the simulation
        sim.set_save_every(save_every);

        sim.run();
        std::string kappa_str = std::to_string(round(kappa * 10) / 10); // round to 1 decimal place
        std::string filename = out_dir.string() + "config_kappa_" + kappa_str + ".csv";
        sim.save_configuration(filename, sweeps);

        std::cout << "Saved configurations to " << filename << std::endl;
    }

    std::cout << "All simulations complete!" << std::endl;
    return 0;
}
