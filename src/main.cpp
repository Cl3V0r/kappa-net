#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
gsl_rng *rng; // random number generator

// 3D vector class
class Vec {
public:
    double x, y, z;
    double length;

    Vec(double inx, double iny, double inz) { x = inx; y = iny; z = inz; }
    Vec() { x = y = z = 0; }

    inline Vec operator+(const Vec &v) {
        Vec temp;
        temp.length = this->length;
        temp.x = x + v.x;
        temp.y = y + v.y;
        temp.z = z + v.z;
        return temp;
    }

    inline bool operator==(const Vec &v) {
        return v.x == x && v.y == y && v.z == z;
    }

    inline double operator*(const Vec &v) { // dot product
        return (v.x * x + v.y * y + v.z * z);
    }

    inline Vec operator*(const double &c) { // scalar multiplication
        Vec temp;
        temp.length = this->length * fabs(c);
        temp.x = c * x;
        temp.y = c * y;
        temp.z = c * z;
        return temp;
    }

    inline void operator=(const Vec &v) {
        x = v.x;
        y = v.y;
        z = v.z;
        length = v.length;
    }

    double magnitude() { return sqrt(*this * *this); }
};

// Generate a random unit vector within a certain angle from v
Vec get_random_vector(const Vec v, double angle) {
    Vec result;
    while (true) {
        gsl_ran_dir_3d(rng, &(result.x), &(result.y), &(result.z)); // random unit vector
        if (result * v > cos(angle) && result * v > 0)
            return result;
    }
}

// Cross product
Vec cross_product(const Vec a, const Vec b) {
    Vec temp;
    temp.x = a.y * b.z - b.y * a.z;
    temp.y = a.z * b.x - b.z * a.x;
    temp.z = a.x * b.y - b.x * a.y;
    return temp;
}

// Scalar (dot) product
double scalar_product(const Vec a, const Vec b) {
    return (b.x * a.x + b.y * a.y + b.z * a.z);
}

// Persistence length function (unused here)
double persistence_length(double beta, double kappa, double b0) {
    return -b0 / log((1 / tanh(beta * kappa) - 1 / (beta * kappa)));
}

int main() {
    // --- Random number initialization ---
    int seed = 0;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, seed);

    // --- Simulation parameters ---
    int num_beads = 100;
    double beta = 1.0;
    double rotation_angle = 0.25;

    // --- Helper variables ---
    vector<double> tangent_corr(num_beads);
    vector<int> tangent_count(num_beads);
    Vec x_axis(1, 0, 0);
    Vec new_tangent = x_axis;
    double kappa;
    int num_sweeps = 1000000;
    int num_measurements = 1;
    int accepted_moves;
    double dE, dE_old, dE_new;
    double variance, angle, energy;

    // --- Loop over different kappa values ---
    for (int j = 0; j < num_measurements; j++) {
        kappa = 15.0;
        accepted_moves = 0;
        rotation_angle = 0.5;

        // --- Initialize chain with random tangent vectors ---
        vector<Vec> tangents;
        for (int i = 0; i < num_beads; i++) {
            new_tangent = get_random_vector(new_tangent, rotation_angle);
            tangents.push_back(new_tangent);
        }

        // --- Monte Carlo simulation ---
        for (int i = 0; i < num_sweeps; i++) {
            // --- Proposed move ---
            int which = gsl_rng_uniform(rng) * (num_beads);
            Vec new_tangent = get_random_vector(tangents[which], rotation_angle);

            // --- Compute old and new local energies ---
            dE_old = 0;
            if (which != 0)
                dE_old += tangents[which - 1] * tangents[which] * (-1.0) + 1;
            if (which != num_beads - 1)
                dE_old += tangents[which] * tangents[which + 1] * (-1.0) + 1;

            dE_new = 0;
            if (which != 0)
                dE_new += tangents[which - 1] * new_tangent * (-1.0) + 1;
            if (which != num_beads - 1)
                dE_new += new_tangent * tangents[which + 1] * (-1.0) + 1;

            dE = kappa * (dE_old - dE_new);

            // --- Metropolis acceptance criterion ---
            if (dE > 0 || gsl_rng_uniform(rng) < exp(dE)) {
                tangents[which] = new_tangent;
                accepted_moves++;
                dE = dE_new;
            } else {
                dE = dE_old;
            }

            // --- Measure tangent correlations ---
            for (int j = 0; j < num_beads; j++)
                for (int k = j + 1; k < num_beads; k++) {
                    tangent_corr[k - j] += tangents[k] * tangents[j];
                    tangent_count[k - j]++;
                }

            // --- Save configurations periodically ---
            if (i >= 10000 && i % 1000 == 0) {
                ofstream data("data/config-kappa-15.csv", ios_base::app);
                ofstream data_energy("data/energy-15.dat", ios_base::app);
                ofstream data_variance("data/variance-15.dat", ios_base::app);
                ofstream accept("data/accept15.dat", ios_base::app);

                accept << accepted_moves / double(i) << " " << kappa << endl;
                accept.close();

                energy = 0;
                variance = 0;

                for (int i = 0; i < num_beads - 1; i++) {
                    angle = scalar_product(tangents[i + 1], tangents[i]);
                    energy += 1.0 - angle;
                    variance += pow((1.0 - angle) * kappa, 2);
                    data << acos(angle) << " ";
                }

                variance /= double(num_beads - 1);
                energy *= kappa / double(num_beads - 1);
                variance -= pow(energy, 2);

                data_variance << variance << endl;
                data_variance.close();
                data_energy << energy << endl;
                data_energy.close();

                data << kappa << endl;
                data.close();
            }

            // --- Adjust rotation angle to maintain acceptance rate ---
            if (accepted_moves / double(i) > 0.65 && i < 200)
                rotation_angle += 0.01;
            else if (accepted_moves / double(i) < 0.35 && i < 200)
                rotation_angle -= 0.01;

            // --- Progress output ---
            cout << j << " of " << num_measurements
                 << " with acceptance rate: " << accepted_moves / double(i)
                 << " kappa " << kappa << "\r";
        }
    }
}
