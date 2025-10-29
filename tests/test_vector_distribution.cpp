#include <catch2/catch_all.hpp>
#include "vector.hpp"
#include <gsl/gsl_rng.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <filesystem>
#include <iostream>

TEST_CASE("Generated vectors are uniformly distributed in cone", "[vector_sampling]") {
    const double angle = 0.5;  // radians
    const int N = 100000;
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, 1234);

    Vec v0(0, 0, 1);
    std::vector<double> cos_thetas, phis;
    cos_thetas.reserve(N);
    phis.reserve(N);

    const char* write_csv_env = std::getenv("WRITE_CSV");
    bool write_csv = (write_csv_env && std::string(write_csv_env) == "1");
    std::ofstream csv_file;

    if (write_csv) {
        std::string base_dir = std::string(PROJECT_SOURCE_DIR) + "/data";
        std::filesystem::create_directories(base_dir);
        csv_file.open(base_dir + "/vector_samples.csv");
        csv_file << "x,y,z\n";
    }


    for (int i = 0; i < N; ++i) {
        Vec r = get_random_vector_fast(v0, angle, rng);
        double cos_theta = r.z;
        double phi = std::atan2(r.y, r.x);
        if (phi < 0) phi += 2 * M_PI;
        cos_thetas.push_back(cos_theta);
        phis.push_back(phi);

        if (write_csv) {
            csv_file << r.x << "," << r.y << "," << r.z << "\n";
        }

    }


    if (write_csv) {
        csv_file.close();
        WARN("Wrote 100k samples to tests/data/vector_samples.csv for visualization.");
    }


    gsl_rng_free(rng);

    // KS test for cos(theta)
    std::sort(cos_thetas.begin(), cos_thetas.end());
    double cos_min = std::cos(angle);
    double ks_max = 0.0;
    for (int i = 0; i < N; ++i) {
        double F_emp = double(i + 1) / N;
        double F_th = (cos_thetas[i] - cos_min) / (1.0 - cos_min);
        ks_max = std::max(ks_max, std::fabs(F_emp - F_th));
    }
    INFO("KS statistic for cos(theta): " << ks_max);
    REQUIRE(ks_max < 0.02);

    // KS test for phi uniformity
    std::sort(phis.begin(), phis.end());
    double ks_max_phi = 0.0;
    for (int i = 0; i < N; ++i) {
        double F_emp = double(i + 1) / N;
        double F_th = phis[i] / (2.0 * M_PI);
        ks_max_phi = std::max(ks_max_phi, std::fabs(F_emp - F_th));
    }
    INFO("KS statistic for phi: " << ks_max_phi);
    REQUIRE(ks_max_phi < 0.02);
}

