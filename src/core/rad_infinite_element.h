/* rad_infinite_element.h -- static / low-frequency INFINITE-ELEMENT (IE) DtN surface operator.
 *
 * The (static) infinite element closes an exterior Laplace problem on a truncation surface Gamma by
 * expanding the decaying exterior field as (surface FE) x (radial decay basis) and statically
 * condensing the radial DOFs onto the surface trace.  On a SPHERE this is IDENTICAL to the Kelvin
 * transformation (same exterior polynomial space; see examples/.../act7_28_ie_vs_kelvin_fair_dtn) --
 * this module is the C++ port of the Python prototype acts7_32/7_33 (the spherical IE), the
 * foundation a later non-spherical (surface-conforming) IE builds on.
 *
 * THE ASSEMBLY (faithful to act7_32 ie_surface_operator):
 *   exterior energy = tensor product   A^{kl} = R1_kl * Mtil + R0_kl * Ktil   (block per radial pair)
 *     R1_kl = int_a^inf rho_k' rho_l' r^2 dr      (radial STIFFNESS, P x P)
 *     R0_kl = int_a^inf rho_k  rho_l      dr      (radial MASS,      P x P)
 *     Mtil  = unit-sphere surface mass    (caller-supplied; = physical M^S for a=1)
 *     Ktil  = unit-sphere Laplace-Beltrami(caller-supplied; gradients' a-factors cancel)
 *   then statically condense the P-1 radial "bubble" levels onto the trace level (level 0) ->
 *   the DtN surface stiffness S_Gamma (N x N), which the caller adds to the interior FE system.
 *   For a single mode n: Mtil->1, Ktil->n(n+1), and eig(S_Gamma, Mtil) -> the analytic ladder (n+1)/a.
 *
 * RADIAL BASIS -- orthogonalized, NEVER naive monomials (the act7_28 conditioning lesson).  With
 * t = a/r in (0,1]:  N_1(t)=t (the vertex / trace-carrying function, value 1 at the surface), and
 * N_k(t) = (L_k(2t-1) - L_{k-2}(2t-1))/(2k-1) for k>=2 (integrated-Legendre bubbles, value 0 at both
 * t=0 and t=1).  span{N_1..N_P} = {(a/r)^1..(a/r)^P} -- the SAME exterior space as the monomials but
 * well-conditioned AND with a CLEAN trace (only N_1 nonzero at the surface -> trace vector g = e_1).
 *
 * Conventions (CLAUDE.md): row-major [target][source]; +N physical sign; pure C++ + MKL (no NGSolve --
 * the surface mass M^S and Laplace-Beltrami K^S are supplied by the caller).
 */
#ifndef RAD_INFINITE_ELEMENT_H
#define RAD_INFINITE_ELEMENT_H

#include <vector>

namespace rad_ie {

/* Gauss-Legendre nodes/weights on [0,1] (n points).  Exposed for the radial-quadrature golden test. */
void GaussLegendre01(int n, std::vector<double>& x, std::vector<double>& w);

/* Orthogonal nodal radial decay operators (a = sphere radius).  Outputs row-major P x P:
 *   R1 = radial stiffness, R0 = radial mass; and the trace vector g (length P, = e_1 for this basis).
 * nq = number of Gauss-Legendre points on [0,1] used for the radial integrals. */
void RadialOperators(int P, double a, std::vector<double>& R1, std::vector<double>& R0,
                     std::vector<double>& g, int nq = 160);

/* Condensed DtN surface operator S (row-major N x N) from the unit-sphere surface mass Mtil and
 * Laplace-Beltrami Ktil (both row-major N x N, symmetric).  Builds the P-level tensor blocks
 * A^{kl} = R1_kl*Mtil + R0_kl*Ktil and statically condenses the bubble levels (1..P-1) onto the
 * trace level (0; nodal basis g=e_1).  a = sphere radius (sets the radial integrals).
 * Returns 0 on success, nonzero LAPACK info on a failed condensation solve. */
int DtNSurfaceOperator(const double* Mtil, const double* Ktil, int N, int P, double a,
                       std::vector<double>& S, int nq = 160);

}  // namespace rad_ie

#endif  // RAD_INFINITE_ELEMENT_H
