// Axisymmetric volume PEEC of a revolved conductor section.
//
// Under an axisymmetric (loop-voltage) drive each meridian triangle of the
// section carries one azimuthal ring current (P0).  The circuit is
//   (R + j omega M) I = V,   R_k = 2 pi r_k / (sigma A_k),
// with M the Galerkin ring partial inductance: the coaxial-filament mutual
// averaged over both triangles.  Self and near pairs subtract
// mu0 r (ln(8 r) - 2 - ln|x - y|) and integrate its logarithm in closed form.
//
// Plain C++ on the radia::host scheduling boundary (rad_host_parallel.h), so
// the same kernel serves the pybind11 module and a MATLAB MEX gateway.
#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace radia {
namespace peec_axisym {

// Meridian triangles: points (r, z) in metres, row-major n_points x 2;
// triangles, row-major n_tri x 3, zero-based point indices.
struct RingSection {
    const double* points = nullptr;
    std::size_t n_points = 0;
    const int* triangles = nullptr;
    std::size_t n_tri = 0;
};

struct RingInductance {
    std::size_t n = 0;
    std::vector<double> M;          // n x n, row-major, symmetric [H]
    std::vector<double> area;       // meridian area of each triangle [m^2]
    std::vector<double> centroid;   // n x 2 (r, z) [m]
};

struct RingSolution {
    std::size_t n = 0;
    std::complex<double> impedance;              // V / I_total at V = 1 [ohm]
    std::vector<std::complex<double>> current;   // ring currents at V = 1 [A]
    std::vector<double> resistance;              // R_k [ohm]
    std::vector<double> cell_loss;               // time-average loss per 1 A peak terminal current [W]
    double true_relative_residual = 0.0;
};

// Pair classes by centroid distance d and sizes s = sqrt(2 A):
// d < near_factor (s_a + s_b): closed-form near rule (7 x 6 points + log term);
// d < moderate_factor (s_a + s_b): 3 x 3 product rule; otherwise centroids.
RingInductance AssembleRingInductance(const RingSection& section,
                                      double near_factor = 3.0,
                                      double moderate_factor = 10.0);

// Solves (R + j omega M) I = 1 with LAPACK and checks the true residual
// against 1e-6 (raises otherwise).
RingSolution SolveLoopVoltage(const RingInductance& L, double frequency_hz, double sigma);

// Coaxial circular filament mutual inductance (Maxwell, Carlson integrals) [H].
double FilamentMutual(double r1, double z1, double r2, double z2);

}  // namespace peec_axisym
}  // namespace radia
