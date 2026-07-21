#pragma once

#include <array>

namespace axifem::numeric {

using Matrix4 = std::array<double, 16>;

struct Q1MagneticElementMatrices {
    Matrix4 stiffness;
    Matrix4 sigma_mass;
};

// Return row-major 4x4 matrices for the node order
// (ra, za), (rb, za), (rb, zb), (ra, zb). The DOFs are nodal A_phi values.
Matrix4 ComputeQ1InverseVandermonde(
    double ra, double rb, double za, double zb);
Matrix4 ComputeQ1MagneticStiffness(
    double ra, double rb, double za, double zb, double mu);
Matrix4 ComputeQ1SigmaMass(
    double ra, double rb, double za, double zb, double sigma);
Q1MagneticElementMatrices ComputeQ1MagneticElementMatrices(
    double ra, double rb, double za, double zb, double mu, double sigma);

}  // namespace axifem::numeric
