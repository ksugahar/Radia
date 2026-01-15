/*-------------------------------------------------------------------------
*
* File name:      rad_elliptic_integral.h
*
* Project:        RADIA
*
* Description:    Elliptic integral functions for analytical coil field
*                 computation (K(k), E(k) complete elliptic integrals)
*                 and arc coil near-field B-field calculation.
*
* References:
*   [1] J.C. Maxwell, "A Treatise on Electricity and Magnetism," Vol. 2,
*       Art. 701-706, Oxford: Clarendon Press, 1873.
*   [2] J.C. Simpson et al., "Simple Analytic Expressions for the Magnetic
*       Field of a Circular Current Loop," NASA/TM-2013-217919, 2001.
*   [3] C. Hastings et al., "Approximations for Digital Computers,"
*       Princeton University Press, 1955.
*   [4] A. Kameari, EMPY_Field library (original implementation)
*   [5] A. Kameari, "Calculation of Transient 3D Eddy Current using
*       Edge-Elements", IEEE Trans. Magn., Vol. 26, No. 2, pp. 466-469, 1990.
*       - Arc coil B-field with incomplete elliptic integrals
*       - Analytical cross-section integration (distances to 4 corners)
*   [6] M. Abramowitz, I.A. Stegun, "Handbook of Mathematical Functions,"
*       National Bureau of Standards, Chapter 17: Elliptic Integrals, 1964.
*   [7] R. Piessens, E. de Doncker-Kapenga, C.W. Ueberhuber, D.K. Kahaner,
*       "QUADPACK: A Subroutine Package for Automatic Integration,"
*       Springer-Verlag, 1983. - Gauss-Kronrod 7-15 adaptive quadrature
*
* Author(s):      A. Kameari (original EMPY), adapted for Radia
*
* First release:  2025
*
-------------------------------------------------------------------------*/

#ifndef __RAD_ELLIPTIC_INTEGRAL_H
#define __RAD_ELLIPTIC_INTEGRAL_H

namespace RadElliptic {

/**
 * Complete elliptic integrals of the first and second kind using
 * Hastings polynomial approximation (CELIDD method).
 *
 * Input:
 *   k2: k^2 where k is the modulus (0 <= k2 < 1)
 *
 * Output:
 *   K: Complete elliptic integral of the first kind K(k)
 *   E: Complete elliptic integral of the second kind E(k)
 *
 * Returns:
 *   0: Success
 *   1: Error (k2 out of range)
 */
int CompleteEllipticIntegrals(double k2, double& K, double& E);

/**
 * Solid angle subtended by a circular current loop at an observation point.
 * Used for magnetic scalar potential computation: Phi = I * Omega / (4*pi)
 *
 * Input:
 *   R:  radial distance from axis (cylindrical coordinates)
 *   Z:  axial distance from coil center
 *   CR: coil radius
 *   CZ: coil center z-coordinate (0 if centered at origin)
 *
 * Output:
 *   Omega: solid angle in steradians [sr]
 *
 * Note: The solid angle is discontinuous across the loop plane (z=CZ)
 *       inside the loop (R < CR). This is physically correct.
 *       On-axis (R=0): Omega = 2*pi * (1 - |z|/sqrt(CR^2 + z^2)) * sign(z)
 *       Off-axis: Uses elliptic integral formula.
 *
 * References:
 *   [1] J.D. Jackson, "Classical Electrodynamics," 3rd ed., Ch. 5.6
 *   [2] MIT 6.013 Electromagnetics, Section 8.3
 */
double CircularLoopSolidAngle(double R, double Z, double CR, double CZ);

/**
 * B-field and A-field computation for a circular current loop using
 * analytical elliptic integral formulas.
 *
 * Input:
 *   R:  radial distance from axis (cylindrical coordinates)
 *   Z:  axial distance from coil center
 *   CR: coil radius
 *   CZ: coil center z-coordinate (0 if centered at origin)
 *   CI: current (A)
 *
 * Output:
 *   BR: radial component of B-field (T)
 *   BZ: axial component of B-field (T)
 *
 * Note: On-axis (R=0) case is handled analytically.
 */
void CircularLoopBField(double R, double Z, double CR, double CZ, double CI,
                        double& BR, double& BZ);

/**
 * Vector potential A_phi for a circular current loop using
 * analytical elliptic integral formulas.
 *
 * Input:
 *   R:  radial distance from axis (cylindrical coordinates)
 *   Z:  axial distance from coil center
 *   CR: coil radius
 *   CZ: coil center z-coordinate (0 if centered at origin)
 *   CI: current (A)
 *
 * Output:
 *   Aphi: azimuthal component of vector potential (T*m)
 *
 * Note: On-axis (R=0), Aphi = 0 by symmetry.
 */
double CircularLoopAPhi(double R, double Z, double CR, double CZ, double CI);

/**
 * Perfect elliptic integral computation for complete integrals K and E.
 * Uses polynomial approximation valid for k2 near 0 and near 1.
 *
 * Input:
 *   k2:  k^2 (modulus squared), 0 <= k2 < 1
 *   kd2: 1 - k2 (complementary modulus squared)
 *
 * Output:
 *   K: Complete elliptic integral of first kind
 *   E: Complete elliptic integral of second kind
 *
 * Returns:
 *   1 if k2 >= EPS_k2 (normal computation)
 *   0 if k2 < EPS_k2 (near-zero case, K and E not meaningful)
 */
int PerfectEllipticIntegral(double k2, double kd2, double* K, double* E);

/**
 * Incomplete elliptic integrals F(phi, m) and E(phi, m).
 * Uses arithmetic-geometric mean algorithm.
 *
 * Reference:
 *   [1] A. Kameari, EMPY implementation (ell_fe function)
 *   [2] M. Abramowitz, I.A. Stegun, "Handbook of Mathematical Functions",
 *       Chapter 17: Elliptic Integrals
 *
 * Input:
 *   phi: amplitude (radians)
 *   m:   parameter m = k^2 (modulus squared)
 *
 * Output:
 *   F: Incomplete elliptic integral of first kind F(phi, m)
 *   E: Incomplete elliptic integral of second kind E(phi, m)
 */
void IncompleteEllipticIntegrals(double phi, double m, double& F, double& E);

/**
 * Incomplete elliptic integrals between two angles.
 * Computes F(phi2,m) - F(phi1,m) and E(phi2,m) - E(phi1,m).
 *
 * Input:
 *   phi1, phi2: amplitude range (radians)
 *   m:          parameter m = k^2 (modulus squared)
 *
 * Output:
 *   F: F(phi2,m) - F(phi1,m)
 *   E: E(phi2,m) - E(phi1,m)
 */
void IncompleteEllipticIntegralsDiff(double phi1, double phi2, double m,
                                      double& F, double& E);

/**
 * Arc coil B-field using incomplete elliptic integrals (far-field/line current).
 *
 * Reference:
 *   [1] A. Kameari, "Calculation of Transient 3D Eddy Current using
 *       Edge-Elements", IEEE Trans. Magn., Vol. 26, No. 2, 1990
 *
 * This implements the analytical formula for the magnetic field of an arc coil
 * when the observation point is far from the coil cross-section.
 * The coil is approximated as a line current along the arc.
 *
 * Input:
 *   rho:        radial distance from axis in local frame (m)
 *   z:          axial distance from coil plane in local frame (m)
 *   rc:         coil radius (m)
 *   phi1, phi2: arc angles relative to observation point (rad)
 *   current:    total current (A)
 *   cross_area: cross-section area (m^2)
 *
 * Output:
 *   B_rho, B_theta, B_z: B-field components in cylindrical coords (T)
 */
void ArcCoilBFieldAnalytical(double rho, double z, double rc,
                              double phi1, double phi2,
                              double current, double cross_area,
                              double& B_rho, double& B_theta, double& B_z);

/**
 * Arc coil field computation with distance-based method selection.
 *
 * This is the main entry point for arc coil field calculation.
 * It automatically selects between:
 *   1. Analytical (line current) formula for far-field
 *   2. Returns false for near-field (caller should use numerical integration)
 *
 * Input:
 *   obs_rho, obs_z: observation point in local cylindrical coords (m)
 *   rc:             coil center radius (m)
 *   phi1, phi2:     arc angles (rad)
 *   a:              half radial width (m)
 *   b:              half axial width (m)
 *   current:        total current (A)
 *
 * Output:
 *   B_rho, B_theta, B_z: B-field components (T)
 *
 * Returns:
 *   true:  analytical formula was used (far-field)
 *   false: numerical integration needed (near-field), B components set to 0
 */
bool ArcCoilBFieldAdaptive(double obs_rho, double obs_z, double rc,
                            double phi1, double phi2,
                            double a, double b, double current,
                            double& B_rho, double& B_theta, double& B_z);

/**
 * Parameter structure for arc coil field integrand.
 *
 * This structure holds all geometric parameters needed for the
 * analytical cross-section integration at each azimuthal angle.
 *
 * Reference:
 *   A. Kameari, EMPY_Field implementation (field_param structure)
 */
struct ArcFieldParams {
    double RHO;       // Observation point radial distance
    double RHOSQ;     // RHO^2
    double ZSQ;       // z^2
    double ASQ;       // rc^2 (coil radius squared)
    double AA;        // rc (coil radius)
    double Z1;        // z (observation point axial distance)
    double ARZSQ;     // (rc - rho)^2 + z^2
    double TAR;       // 2 * rho * rc
    double AREA;      // Cross-section area
    double DFT;       // Distance threshold for line approximation

    double AL;        // rc - a (inner coil radius)
    double AU;        // rc + a (outer coil radius)
    double ZZL;       // z + b
    double ZZU;       // z - b
    double ALSQ;      // (rc - a - rho)^2
    double AUSQ;      // (rc + a - rho)^2
    double ZLSQ;      // (z + b)^2
    double ZUSQ;      // (z - b)^2
    double ARZ[4];    // Squared distances to 4 corners at phi=0
    double ARL;       // 2 * rho * (rc - a)
    double ARU;       // 2 * rho * (rc + a)

    int calcA;        // Flag: compute vector potential A
    int calcB;        // Flag: compute magnetic field B
    int symmetricArc; // Flag: using symmetry (only positive half)
};

/**
 * Arc coil field integrand with analytical cross-section integration.
 *
 * For each azimuthal angle ALPHA, computes the contribution to B and A fields
 * by analytically integrating over the rectangular cross-section using
 * distances to 4 corners (R[0-3]).
 *
 * The analytical formula for cross-section integration is based on:
 *   - 1/r integral over rectangle using logarithm and arctangent
 *   - Distances to 4 corners: R[i] = sqrt((r_corner - r_obs)^2 + dz^2)
 *
 * Reference:
 *   [5] A. Kameari, "Calculation of Transient 3D Eddy Current using
 *       Edge-Elements", IEEE Trans. Magn., Vol. 26, No. 2, 1990.
 *       See also EMPY_Field implementation (arc_field_integrand function)
 *
 * Input:
 *   ALPHA: azimuthal angle (radians)
 *   param: pointer to ArcFieldParams structure
 *
 * Output:
 *   F: array of integrand values
 *      - If calcB && !symmetricArc: F[0]=B_rho, F[1]=B_theta, F[2]=B_z
 *      - If calcB && symmetricArc:  F[0]=B_rho, F[1]=B_z
 *      - Additional A-field terms follow if calcA is set
 */
void ArcFieldIntegrand(double ALPHA, double* F, const ArcFieldParams* param);

/**
 * Arc coil near-field integration using analytical cross-section formula
 * with adaptive 1D integration over azimuthal direction.
 *
 * This computes the B-field of an arc coil when the observation point is
 * near the coil cross-section. Instead of 2D Gauss quadrature over the
 * cross-section, it uses:
 *   1. Analytical formula for cross-section integral at each phi [5]
 *   2. Adaptive Gauss-Kronrod 7-15 integration over azimuth [7]
 *
 * References:
 *   [5] A. Kameari, "Calculation of Transient 3D Eddy Current using
 *       Edge-Elements", IEEE Trans. Magn., Vol. 26, No. 2, 1990.
 *       - Analytical cross-section integration formula
 *   [7] R. Piessens et al., "QUADPACK: A Subroutine Package for
 *       Automatic Integration," Springer-Verlag, 1983.
 *       - Gauss-Kronrod 7-15 adaptive quadrature (QAG algorithm)
 *
 * Input:
 *   obs_rho: observation point radial distance from axis (m)
 *   obs_z:   observation point axial distance from coil center (m)
 *   rc:      coil center radius (m)
 *   a:       half radial width (m)
 *   b:       half axial width (m)
 *   phi1, phi2: arc angles (rad)
 *   current: total current (A)
 *
 * Output:
 *   B_rho, B_theta, B_z: B-field components in cylindrical coords (T)
 */
void ArcCoilBFieldNearField(double obs_rho, double obs_z, double rc,
                             double a, double b,
                             double phi1, double phi2,
                             double current,
                             double& B_rho, double& B_theta, double& B_z);

} // namespace RadElliptic

#endif // __RAD_ELLIPTIC_INTEGRAL_H
