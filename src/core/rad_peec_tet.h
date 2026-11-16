// Tetrahedral volume PEEC: lowest-order face (RT0) currents in a conductor.
//
// On a flat tetrahedron every lowest-order H(div) basis function is affine,
//   w(x) = c (x - v_opp),
// with v_opp the vertex opposite its face.  The caller supplies c and the
// DOF number of each (tetrahedron, local face) from the finite-element space
// that owns orientation and scaling (NGSolve), so this kernel never chooses
// signs.  The partial inductance is the Galerkin form
//   L_ij = mu0/(4 pi) int int w_i(x) . w_j(y) / |x - y| dx dy.
// Near pairs integrate the inner tetrahedron in closed form (face formulas
// for int 1/R and int R) at Gauss points of the outer one; moderate pairs use
// a product Gauss rule; far pairs the centroids.  The circuit
//   [R + j omega L, D^T; D, 0] [I; phi] = [b; 0]
// is solved densely with LAPACK and its true residual is checked.
//
// Plain C++ on the radia::host scheduling boundary (rad_host_parallel.h).
#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace radia {
namespace peec_tet {

struct TetMesh {
    const double* points = nullptr;     // n_points x 3, row-major [m]
    std::size_t n_points = 0;
    const int* tets = nullptr;          // n_tet x 4 point indices
    std::size_t n_tet = 0;
    const double* face_coef = nullptr;  // n_tet x 4: c of the face opposite local vertex k
    const int* face_dof = nullptr;      // n_tet x 4: DOF of that face, -1 when not an unknown
    std::size_t n_dof = 0;
};

// Pair classes by centroid distance d and bounding radii s:
// d < near_factor (s_a + s_b): closed-form inner, Gauss outer (near_order^3 points);
// d < moderate_factor (s_a + s_b): 2^3 x 2^3 product rule; otherwise centroids.
std::vector<double> AssembleInductance(const TetMesh& mesh, double near_factor = 3.0,
                                       double moderate_factor = 10.0, int near_order = 3);

// Closed-form potentials of a uniform unit density on one tetrahedron at x:
//   p0 = int 1/|x - y| dy,  p1 = int y/|x - y| dy.
void TetPotential(const double x[3], const double v[4][3], double& p0, double p1[3]);

struct Triplets {
    std::vector<int> row, col;
    std::vector<double> value;
};

struct KKTSolution {
    std::vector<std::complex<double>> current;    // DOF currents
    std::vector<std::complex<double>> potential;  // constraint multipliers
    double true_relative_residual = 0.0;
};

// Solves [R + j omega L, D^T; D, 0] [I; phi] = [b; 0] (dense LU), raising when
// the true relative residual exceeds 1e-6 or the dense system and its factor
// would need more than max_bytes.  L is n_dof x n_dof row-major,
// R is n_dof x n_dof and D is n_con x n_dof, both as triplets.
KKTSolution SolveKKT(std::size_t n_dof, const std::vector<double>& L, const Triplets& R,
                     std::size_t n_con, const Triplets& D,
                     const std::vector<std::complex<double>>& b, double omega, double max_bytes);

}  // namespace peec_tet
}  // namespace radia
