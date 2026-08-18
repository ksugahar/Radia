// Native construction kernel for the fourth-order Dragt-Finn Lie map.
//
// This is the value path (parameter_jacobians=False) of
// radia.accelerator_lie_topopt._fourth_order_lie_map_from_vector_potential_
// polynomials: per-segment degree-5 Hamiltonian jets built from (Ay, As)
// polynomial coefficient arrays on the shared 462-monomial basis, the
// nonautonomous stage-jet RK4 flow of the factorial R/T/U/V map ODE, and
// the sequential segment composition.  Python keeps the Dragt-Finn
// factorization, gates the result, and packages it; MATLAB reaches the
// same kernel through a standalone MEX ABI (plain arrays in/out, no
// Python objects) per the capability-parity policy.
//
// The forward-AD (topopt jacobian) mode stays in Python: it is a co-valid
// path for an input class this kernel does not cover, not a PoC.

#ifndef RAD_LIE_MAP_KERNEL_H
#define RAD_LIE_MAP_KERNEL_H

#include <cstddef>

namespace rad_lie {

// Integrate the fourth-order map tensors of a segmented s-polynomial
// vector-potential Hamiltonian.
//
// Ay, As:   (n_segments, s_order_count, d+1, d+1) row-major coefficient
//           arrays; coefficient [k][i][j] multiplies zeta^k x^i y^j with
//           zeta in [-1, 1] across the segment.  s_order_count == 1 is the
//           constant-jet (autonomous) case.
// lengths:  (n_segments) segment lengths in metres.
// curvature:(n_segments, curvature_columns) zeta-polynomial reference
//           curvature per segment; curvature_columns == 1 is constant.
// poisson:  (6,6) canonical Poisson matrix (row-major).
// longitudinal_covariant: 1 = As already covariant, 0 = physical tangent
//           projection (multiplied by the metric 1 + h x here).
//
// Outputs: R (6,6), T (6,6,6), U (6,6,6,6), V (6,6,6,6,6) row-major map
// tensors, linear_out (n_segments, 6) worst-stage Hamiltonian linear term
// per segment, worst_linear_out the global maximum |H1|.
//
// Throws std::invalid_argument on malformed inputs, on a reference-orbit
// gate violation ("reference orbit is not a Hamiltonian fixed trajectory
// ..."), and on the step-cap overflow, mirroring the Python messages.
void LieMapTensorsFromSpolyArrays(
    const double* Ay,
    const double* As,
    std::size_t n_segments,
    std::size_t s_order_count,
    std::size_t transverse_order_count,   // d+1, 2..6
    const double* lengths,
    const double* curvature,
    std::size_t curvature_columns,
    double magnetic_rigidity,
    double curvature_sign,
    double reference_beta,
    int longitudinal_covariant,
    const double* poisson,
    double maximum_step_m,
    long long maximum_steps,
    double reference_orbit_tolerance,
    double* R_out,
    double* T_out,
    double* U_out,
    double* V_out,
    double* linear_out,
    double* worst_linear_out);

}  // namespace rad_lie

#endif  // RAD_LIE_MAP_KERNEL_H
