#include "polymer.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>

PolymerSimulation::PolymerSimulation(int num_beads, double kappa, double beta,
                                     double rotation_angle, int sweeps,
                                     unsigned long seed)
    : num_beads_(num_beads), kappa_(kappa), beta_(beta),
      rotation_angle_(rotation_angle), sweeps_(sweeps)
{
    rng_ = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng_, seed);

    tangents_.resize(num_beads_);
    tangent_corr_.assign(num_beads_, 0.0);
    tangent_count_.assign(num_beads_, 0);
}

void PolymerSimulation::initialize_chain() {
    Vec current(1, 0, 0);
    for (int i = 0; i < num_beads_; ++i) {
        current = get_random_vector(current, rotation_angle_, rng_);
        tangents_[i] = current;
    }
}

double PolymerSimulation::compute_local_energy(int index, const Vec &new_tangent) const {
    double energy = 0.0;
    if (index != 0)
        energy += tangents_[index - 1] * new_tangent * (-1.0) + 1;
    if (index != num_beads_ - 1)
        energy += new_tangent * tangents_[index + 1] * (-1.0) + 1;
    return energy;
}

void PolymerSimulation::metropolis_step() {
    int which = gsl_rng_uniform(rng_) * num_beads_;
    Vec new_tangent = get_random_vector(tangents_[which], rotation_angle_, rng_);

    double dE_old = compute_local_energy(which, tangents_[which]);
    double dE_new = compute_local_energy(which, new_tangent);
    double dE = kappa_ * (dE_old - dE_new);

    if (dE > 0 || gsl_rng_uniform(rng_) < exp(dE)) {
        tangents_[which] = new_tangent;
    }
}

void PolymerSimulation::measure_tangent_correlation() {
    for (int i = 0; i < num_beads_; ++i)
        for (int j = i + 1; j < num_beads_; ++j) {
            tangent_corr_[j - i] += tangents_[j] * tangents_[i];
            tangent_count_[j - i]++;
        }
}

void PolymerSimulation::run() {
    initialize_chain();

    // Equilibration
    for (int step = 0; step < 100000; ++step)
        metropolis_step();

    // Production
    for (int step = 0; step < sweeps_; ++step) {
        metropolis_step();

        if (step % num_beads_ == 0)
            measure_tangent_correlation();

        // Save configuration every save_every_ sweeps
        if (step >= 10000 && step % save_every_ == 0) {
            std::string filename = "data/raw/config_kappa_" + std::to_string(int(kappa_)) + ".csv";
            save_configuration(filename, step);
        }
    }

    // Normalize correlations
    for (size_t i = 0; i < tangent_corr_.size(); ++i)
        if (tangent_count_[i] > 0)
            tangent_corr_[i] /= static_cast<double>(tangent_count_[i]);
}


void PolymerSimulation::save_configuration(const std::string &filename, int step) {
    std::ofstream data(filename, std::ios_base::app);  // append mode
    for (int i = 0; i < num_beads_ - 1; ++i) {
        double dot = scalar_product(tangents_[i + 1], tangents_[i]);
        data << std::acos(std::clamp(dot, -1.0, 1.0)) << " ";
    }
    data << kappa_ << "\n";  // store kappa as the last column
}

double PolymerSimulation::compute_energy() const {
    double E = 0;
    for (int i = 0; i < num_beads_ - 1; ++i)
        E += 1.0 - scalar_product(tangents_[i + 1], tangents_[i]);
    return E * kappa_ / double(num_beads_ - 1);
}

double PolymerSimulation::compute_variance() const {
    double meanE = compute_energy();
    double var = 0;
    for (int i = 0; i < num_beads_ - 1; ++i) {
        double term = (1.0 - scalar_product(tangents_[i + 1], tangents_[i])) * kappa_;
        var += pow(term, 2);
    }
    var /= double(num_beads_ - 1);
    var -= pow(meanE, 2);
    return var;
}

double persistence_length(double beta, double kappa, double b0) {
    return -b0 / log((1 / tanh(beta * kappa) - 1 / (beta * kappa)));
}
