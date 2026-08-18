// Batched application of Dragt-Finn factorized Lie maps to 6D ensembles.
//
// This is the tracking-side hot path of the accelerator Lie-map route: the
// map is built once per orbit (Python, radia.accelerator_lie_topopt) and the
// particle ensemble is then advanced through the factor chain
//     exp(:f5:) -> exp(:f4:) -> exp(:f3:) -> R
// with per-generator implicit-midpoint flows (symplectic at finite
// amplitude).  The kernel mirrors apply_dragt_finn_map /
// apply_homogeneous_lie_generator exactly: same predictor, same Newton
// residual and 6x6 solve, same infinity-norm convergence criterion, so the
// batch result matches the single-state Python reference to roundoff.

#ifndef RAD_LIE_MAP_BATCH_H
#define RAD_LIE_MAP_BATCH_H

#include <cstddef>

namespace rad_lie {

// Advance n_states canonical 6D states through the Dragt-Finn factors.
//
// R:        (6,6) row-major linear part.
// f3:       (6,6,6) cubic generator tensor.
// f4:       (6,6,6,6) quartic generator tensor.
// f5:       (6,6,6,6,6) quintic generator tensor, or nullptr for a
//           third-order factorization (no f5 factor).
// poisson:  (6,6) row-major canonical Poisson matrix.
// states_in / states_out: (n_states, 6) row-major; out may alias in.
//
// Gradient/Hessian conventions (identical to the Python reference):
//   deg 3: G_i = T_ijk z_j z_k / 2        H_ij = T_ijk z_k
//   deg 4: G_i = T_ijkl z_j z_k z_l / 6   H_ij = T_ijkl z_k z_l / 2
//   deg 5: G_i = T_ijklm z_j..z_m / 24    H_ij = T_ijklm z_k z_l z_m / 6
//
// Throws std::invalid_argument on malformed controls and
// std::runtime_error when an implicit-midpoint Newton flow fails to
// converge for any particle.
void ApplyDragtFinnMapBatch(
    const double* R,
    const double* f3,
    const double* f4,
    const double* f5,
    const double* poisson,
    const double* states_in,
    double* states_out,
    std::size_t n_states,
    int substeps,
    double newton_tolerance,
    int maximum_newton_iterations);

}  // namespace rad_lie

#endif  // RAD_LIE_MAP_BATCH_H
