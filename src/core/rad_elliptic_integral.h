/*-------------------------------------------------------------------------
*
* File name:      rad_elliptic_integral.h
*
* Project:        RADIA
*
* Description:    Elliptic integral functions for analytical coil field
*                 computation (K(k), E(k) complete elliptic integrals).
*
* References:
*   [1] J.C. Maxwell, "A Treatise on Electricity and Magnetism," Vol. 2,
*       Art. 701-706, Oxford: Clarendon Press, 1873.
*   [2] J.C. Simpson et al., "Simple Analytic Expressions for the Magnetic
*       Field of a Circular Current Loop," NASA/TM-2013-217919, 2001.
*   [3] C. Hastings et al., "Approximations for Digital Computers,"
*       Princeton University Press, 1955.
*   [4] M. Abramowitz, I.A. Stegun, "Handbook of Mathematical Functions,"
*       National Bureau of Standards, Chapter 17: Elliptic Integrals, 1964.
*
* Author(s):      Radia Development Team
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
 *   [1] M. Abramowitz, I.A. Stegun, "Handbook of Mathematical Functions",
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

} // namespace RadElliptic

#endif // __RAD_ELLIPTIC_INTEGRAL_H
