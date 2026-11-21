// PRIMA model-order reduction of a dense branch PEEC port model (CPU).
//
// Full model (modified nodal analysis, x = [I; V]):
//   (G + s E) x = B u,   y = B^T x
//   E = [[L, 0], [0, 0]],  G = [[R, A^T], [-A, 0]],  B = [0; P]
// with branch inductance L (dense, symmetric), branch resistance R
// (diagonal, > 0), reduced incidence A (ground node removed) and port
// columns P (+1 at the positive node, -1 at the negative node).  The
// input u is the port current, the output y the port voltage
// (impedance orientation).
//
// Reduction: block Arnoldi on (G + s0 E)^{-1} E from (G + s0 E)^{-1} B,
// real s0 >= 0, then the congruence E_q = V^T E V, G_q = V^T G V,
// B_q = V^T B.  G + G^T = diag(2R, 0) >= 0 and E >= 0 are kept by the
// congruence, so the reduced model is passive (Odabasioglu, Celik and
// Pileggi, IEEE TCAD 17(8), 1998).
//
// The shifted saddle system is solved densely through Cholesky
// factors of K = R + s0 L and of the Schur complement A K^{-1} A^T
// (LAPACK dpotrf/dpotrs).  Every solve checks its true residual
// against RELATIVE_LIMIT and raises.  No sparse direct fallback.
#pragma once

#include <cstddef>
#include <vector>

namespace radia {
namespace peec_prima {

constexpr double RELATIVE_LIMIT = 1e-6;

struct PortModelInput {
    std::size_t n_branch = 0;
    std::vector<double> L;            // n_branch x n_branch, row-major, symmetric
    std::vector<double> R;            // n_branch, > 0
    std::vector<int> branch_nodes;    // n_branch x 2, current from [0] to [1]
    std::size_t n_nodes = 0;
    std::vector<int> port_plus;       // per port
    std::vector<int> port_minus;      // per port; all share ground = port_minus[0]
};

struct ReducedModel {
    std::size_t q = 0;                // reduced order actually built
    std::size_t n_ports = 0;
    std::size_t n_full = 0;           // N + number of non-ground nodes
    std::vector<double> V;            // n_full x q, row-major, orthonormal columns
    std::vector<double> E;            // q x q, symmetric
    std::vector<double> G;            // q x q
    std::vector<double> B;            // q x n_ports
    std::vector<int> deflated;        // Krylov candidates dropped (index of source column)
    double max_solve_residual = 0.0;  // worst relative true residual of the shifted solves
    double min_E_eig_over_max = 0.0;  // reduced energy spectrum, checked >= -roundoff
};

// Raises std::invalid_argument on malformed input and std::runtime_error on
// a failed factorisation, a solve above RELATIVE_LIMIT or a reduced model
// with negative magnetic energy.
ReducedModel Reduce(const PortModelInput& in, std::size_t order, double s0);

}  // namespace peec_prima
}  // namespace radia
