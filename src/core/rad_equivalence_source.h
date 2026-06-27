/*-------------------------------------------------------------------------
*
* File name:      rad_equivalence_source.h
*
* Project:        RADIA
*
* Description:    Equivalence-theorem near-field source kernels
*                 (Schelkunoff/Love).  CST Near-Field-Source equivalent
*                 for the Radia/NGSolve stack.
*
*                 Given (E, H) on a closed surface around all primary
*                 sources, reconstruct the exterior field at any obs
*                 point via the Stratton-Chu surface integral.  See
*                 docs/equivalence_source/CPP_DESIGN.md for theory and
*                 radia_mcp.fem.equivalence_source_knowledge for the
*                 user-level recipes.
*
*                 Phase A (this file): magnetostatic reduction (omega = 0)
*                 only -- 1/R kernel.  Phase B adds the time-harmonic
*                 dyadic Green's function (exp(-jkR)/R + (1/k^2) grad-grad
*                 psi correction).
*
* First release:  2026-05-26 (Phase A)
*
-------------------------------------------------------------------------*/

#ifndef __RAD_EQUIVALENCE_SOURCE_H
#define __RAD_EQUIVALENCE_SOURCE_H

namespace radia { namespace eqsrc {

/**
 * Magnetostatic H reconstruction at observation points from per-face
 * surface sources on a closed triangulated surface.
 *
 * Schelkunoff/Love scalar-potential + vector-potential representation,
 * differentiated inside the integral, gives the static reduction
 *
 *     H(r) = (1/(4 pi)) integral_dOmega
 *              { grad(1/R) x J_s  -  (n . H_s) grad(1/R) } dS
 *
 *     where R = |r - r'|,  grad(1/R) = -(r - r') / R^3
 *           J_s = n x H_s        (equivalent surface current)
 *           (n . H_s) = rho_m / mu_0  (equivalent surface charge / mu_0)
 *
 * Per-face quadrature: single centroid sample (one quadrature point
 * per face).  Higher-order Gauss quadrature is a roadmap item.
 *
 * Sign convention is fixed by the empirical magnetic-dipole unit test:
 * m = m_z z-hat on-axis gives Hz > 0 outside.  See the Phase 1 golden
 * test in validation_test/equivalence_source/phase1_static_coil.py for the
 * end-to-end verification (0.83% PASS at 14160 surface panels).
 *
 * Memory layout: all arrays are C-contiguous row-major
 * (CLAUDE.md row-major policy).
 *
 * Parallelisation: NGSolve TaskManager over observation points
 * (CLAUDE.md TaskManager policy).
 *
 * @param centroids  (n_faces, 3) row-major: face centroid positions [m]
 * @param normals    (n_faces, 3) row-major: OUTWARD unit normals
 * @param areas      (n_faces,)   row-major: face areas [m^2]
 * @param H_surf     (n_faces, 3) row-major: H phasor REAL part at
 *                                            centroid [A/m].  For pure
 *                                            magnetostatic the imag part
 *                                            is zero by definition; this
 *                                            kernel ignores it.
 * @param n_faces    number of surface panels
 * @param obs        (n_obs, 3)   row-major: observation points [m].
 *                                            Must be OUTSIDE the closed
 *                                            surface for the result to
 *                                            be physically meaningful;
 *                                            evaluation INSIDE returns
 *                                            ~0 (the equivalence-theorem
 *                                            null-field property).
 * @param n_obs      number of obs points
 * @param H_out      (n_obs, 3)   row-major: reconstructed H output
 *                                            [A/m] (overwritten).
 * @param n_threads  0 = TaskManager default; > 0 = request that many.
 */
void EvaluateStaticH(
    const double* centroids,
    const double* normals,
    const double* areas,
    const double* H_surf,
    int n_faces,
    const double* obs,
    int n_obs,
    double* H_out,
    int n_threads);

/**
 * Time-harmonic (E, H) reconstruction at observation points from per-face
 * complex (E_s, H_s) on a closed triangulated surface.  Full Stratton-Chu
 * with the DYADIC Green's function (I + grad-grad/k^2) psi -- correctly
 * reproduces deep-near-field behaviour where the scalar form undershoots
 * by up to a factor of 3.
 *
 * Formula (per docs/equivalence_source/CPP_DESIGN.md sec 3.2):
 *
 *   psi(r,r')  = exp(-jkR) / (4 pi R),  k = omega sqrt(mu_0 eps_0)
 *   h          = jkR (complex)
 *   G_bar . v  = psi {(h^2+h+1)/h^2 v  -  (h^2+3h+3)/h^2 (v.R_hat) R_hat}
 *   grad_psi   = -(jk + 1/R) psi (r - r') / R
 *
 *   E(r) = integral { -jw mu_0 G_bar . J_s  +  grad_psi x M_s } dS'
 *   H(r) = integral {  jw eps_0 G_bar . M_s  +  grad_psi x J_s } dS'
 *
 *   J_s = n x H_s,  M_s = n x E_s.
 *
 * The dyadic Green-function form already contains the longitudinal
 * surface-charge contribution.  Do not add separate scalar-charge
 * terms to this kernel; that was the source of the historical
 * Hertzian-dipole Phase 2 mismatch.
 *
 * Sign convention: same as EvaluateStaticH (verified against analytic
 * magnetic dipole + Hertzian dipole far field).
 *
 * Memory layout: all arrays C-contiguous row-major.
 * Parallel: NGSolve TaskManager over obs points.
 *
 * @param centroids  (n_faces, 3) row-major  -- face centroids [m]
 * @param normals    (n_faces, 3) row-major  -- OUTWARD unit normals
 * @param areas      (n_faces,)   row-major  -- face areas [m^2]
 * @param E_re,E_im  (n_faces, 3) row-major  -- E phasor at centroid [V/m]
 * @param H_re,H_im  (n_faces, 3) row-major  -- H phasor at centroid [A/m]
 * @param n_faces    number of surface panels
 * @param obs        (n_obs, 3) row-major    -- observation points [m]
 * @param n_obs      number of obs points
 * @param omega      angular frequency [rad/s]  (must be > 0; for static
 *                   problems use EvaluateStaticH instead)
 * @param E_re_out, E_im_out  (n_obs, 3) row-major -- reconstructed E
 * @param H_re_out, H_im_out  (n_obs, 3) row-major -- reconstructed H
 * @param n_threads  0 = TaskManager default; > 0 = request that many
 */
void EvaluateHarmonic(
    const double* centroids,
    const double* normals,
    const double* areas,
    const double* E_re, const double* E_im,
    const double* H_re, const double* H_im,
    int n_faces,
    const double* obs,
    int n_obs,
    double omega,
    double* E_re_out, double* E_im_out,
    double* H_re_out, double* H_im_out,
    int n_threads);

}}  // namespace radia::eqsrc

#endif  // __RAD_EQUIVALENCE_SOURCE_H
