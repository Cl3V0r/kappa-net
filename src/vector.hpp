#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <cmath>
#include <string>
#include <gsl/gsl_rng.h>  // Needed for gsl_rng*

class Vec {
public:
    double x, y, z;
    double length;

    Vec(double inx, double iny, double inz);
    Vec();
    Vec(const Vec &other) = default;

    Vec operator+(const Vec &v) const;
    bool operator==(const Vec &v) const;
    double operator*(const Vec &v) const;
    Vec operator*(const double &c) const;
    Vec &operator=(const Vec &v);
    double magnitude() const;
};

// Utility functions
Vec get_random_vector(const Vec &v, double angle, gsl_rng *rng);
Vec get_random_vector_fast(const Vec &v, double angle, gsl_rng *rng);
Vec cross_product(const Vec &a, const Vec &b);
double scalar_product(const Vec &a, const Vec &b);

#endif // VECTOR_HPP

