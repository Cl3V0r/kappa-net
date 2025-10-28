#ifndef POLYMER_HPP
#define POLYMER_HPP

#include <vector>
#include <gsl/gsl_rng.h>
#include <algorithm>
#include "vector.hpp"

class PolymerSimulation {
public:
    PolymerSimulation(int num_beads, double kappa, double beta,
                      double rotation_angle, int sweeps,
                      unsigned long seed = 0);

    void run();
    void set_save_every(int n) { save_every_ = n; }
    void save_configuration(const std::string &filename, int step);
    double compute_energy() const;
    double compute_variance() const;

    // --- new: access to measured correlations ---
    const std::vector<double> &get_tangent_correlation() const { return tangent_corr_; }
    const std::vector<int>    &get_tangent_counts() const      { return tangent_count_; }

private:
    int num_beads_;
    double kappa_;
    double beta_;
    double rotation_angle_;
    int sweeps_;
    int save_every_ = 1000;   // default saving interval

    gsl_rng *rng_;
    std::vector<Vec> tangents_;
    std::vector<double> tangent_corr_;
    std::vector<int> tangent_count_;

    void initialize_chain();
    void metropolis_step();
    void measure_tangent_correlation();
    double compute_local_energy(int index, const Vec &new_tangent) const;
};

double persistence_length(double beta, double kappa, double b0);

#endif // POLYMER_HPP
