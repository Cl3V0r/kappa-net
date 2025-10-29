#include "vector.hpp"
#include <gsl/gsl_randist.h>
#include <cmath>

Vec::Vec(double inx, double iny, double inz) : x(inx), y(iny), z(inz), length(0) {}
Vec::Vec() : x(0), y(0), z(0), length(0) {}

Vec Vec::operator+(const Vec &v) const {
    Vec result;
    result.length = this->length;
    result.x = x + v.x;
    result.y = y + v.y;
    result.z = z + v.z;
    return result;
}

bool Vec::operator==(const Vec &v) const {
    return (x == v.x && y == v.y && z == v.z);
}

double Vec::operator*(const Vec &v) const {  // dot product
    return (x * v.x + y * v.y + z * v.z);
}

Vec Vec::operator*(const double &c) const {  // scalar multiplication
    Vec result;
    result.length = this->length * std::fabs(c);
    result.x = c * x;
    result.y = c * y;
    result.z = c * z;
    return result;
}

Vec &Vec::operator=(const Vec &v) {
    if (this != &v) {
        x = v.x;
        y = v.y;
        z = v.z;
        length = v.length;
    }
    return *this;
}

double Vec::magnitude() const {
    return std::sqrt((*this) * (*this));
}

// Generate a random unit vector within a cone of "angle" from v
Vec get_random_vector(const Vec &v, double angle, gsl_rng *rng) {
    Vec result;
    while (true) {
        gsl_ran_dir_3d(rng, &result.x, &result.y, &result.z);  // random unit vector
        double dot = result * v;
        if (dot > std::cos(angle) && dot > 0)
            return result;
    }
}

Vec get_random_vector_fast(const Vec &v, double angle, gsl_rng *rng) {
    // Normalize input vector v
    Vec axis = v;
    double norm = axis.magnitude();
    if (norm == 0.0)
        return Vec(1, 0, 0);
    axis = axis * (1.0 / norm);

    // Uniform sampling in cone
    double cos_angle = std::cos(angle);
    double cos_theta = gsl_rng_uniform(rng) * (1.0 - cos_angle) + cos_angle;
    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    double phi = gsl_rng_uniform(rng) * 2.0 * M_PI;

    // Construct orthonormal basis (u, w)
    Vec u;
    if (std::fabs(axis.x) > 0.1 || std::fabs(axis.y) > 0.1) {
        u = Vec(-axis.y, axis.x, 0.0);
    } else {
        u = Vec(0.0, -axis.z, axis.y);
    }
    double u_norm = u.magnitude();
    if (u_norm == 0.0) u = Vec(1, 0, 0);
    else u = u * (1.0 / u_norm);
    Vec w = cross_product(axis, u);

    // Combine in cone coordinates
    Vec result = (u * (sin_theta * std::cos(phi)) +
                  w * (sin_theta * std::sin(phi)) +
                  axis * cos_theta);
    return result * (1.0 / result.magnitude());
}

Vec cross_product(const Vec &a, const Vec &b) {
    Vec result;
    result.x = a.y * b.z - b.y * a.z;
    result.y = a.z * b.x - b.z * a.x;
    result.z = a.x * b.y - b.x * a.y;
    return result;
}

double scalar_product(const Vec &a, const Vec &b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}
