#include "polymer.hpp"
#include <catch2/catch_all.hpp>
#include <cmath>
#include <vector>
#include <numeric>
#include <iostream>

TEST_CASE("Measured tangent correlation decays exponentially", "[correlation]") {
    int num_beads = 40;
    double kappa = 8.0;
    double beta = 1.0;
    double rotation_angle = 0.3;
    int sweeps = 5000000;   // more sweeps for smoother average

    PolymerSimulation sim(num_beads, kappa, beta, rotation_angle, sweeps);
    sim.run();

    const auto &corr = sim.get_tangent_correlation();
    const auto &counts = sim.get_tangent_counts();

    double Lp = persistence_length(beta, kappa, 1.0);
    REQUIRE(Lp > 0.0);

    for (int j = 1; j < num_beads; ++j) { // check only first few
        if (j >= (int)corr.size() || counts[j] == 0) continue;

        double measured = corr[j];
        double expected = std::exp(-static_cast<double>(j) / Lp);

        // j-dependent tolerance: more relaxed for higher j
        double tolerance = 0.1 + 0.02 * j;

        INFO("j = " << j
             << ", measured = " << measured
             << ", expected = " << expected
             << ", relative error = " << std::fabs(measured - expected) / expected
             << ", tolerance = " << tolerance);
        REQUIRE(std::fabs(measured - expected) / expected < tolerance);
    }
}
