/*-------------------------------------------------------------------------
*
* File name:      rad_elliptic_integral.cpp
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
*   [2] J.C. Simpson, J.E. Lane, C.D. Immer, R.C. Youngquist,
*       "Simple Analytic Expressions for the Magnetic Field of a
*       Circular Current Loop," NASA/TM-2013-217919, 2001.
*   [3] C. Hastings, J.T. Hayward, J.P. Wong, "Approximations for
*       Digital Computers," Princeton University Press, 1955.
*   [4] W.J. Cody, "Chebyshev Approximations for the Complete Elliptic
*       Integrals K and E," Math. Comp. 19(92), pp. 105-112, 1965.
*   [5] A. Kameari, EMPY_Field library (original implementation)
*   [6] A. Kameari, "Calculation of Transient 3D Eddy Current using
*       Edge-Elements", IEEE Trans. Magn., Vol. 26, No. 2, pp. 466-469, 1990.
*       - Arc coil B-field with incomplete elliptic integrals
*       - Analytical cross-section integration (distances to 4 corners)
*   [7] M. Abramowitz, I.A. Stegun, "Handbook of Mathematical Functions,"
*       National Bureau of Standards, Chapter 17: Elliptic Integrals, 1964.
*   [8] R. Piessens, E. de Doncker-Kapenga, C.W. Ueberhuber, D.K. Kahaner,
*       "QUADPACK: A Subroutine Package for Automatic Integration,"
*       Springer-Verlag, 1983. - Gauss-Kronrod 7-15 adaptive quadrature
*
* Author(s):      A. Kameari (original EMPY), adapted for Radia
*
* First release:  2025
*
-------------------------------------------------------------------------*/

#include "rad_elliptic_integral.h"
#include <cmath>
#include <algorithm>

namespace RadElliptic {

// Physical constants
static const double PI = 3.141592653589793238;
static const double MU0 = 1.2566370614359173e-6;  // mu_0 = 4*pi*1e-7 H/m
static const double MU0_DIV_4PI = 1.0e-7;          // mu_0 / (4*pi) for Biot-Savart

// Tolerance for zero detection
static const double ZERO_TOL = 1.0e-10;

//-------------------------------------------------------------------------
// Complete elliptic integrals K(k) and E(k) using Hastings polynomial
// approximation (CELIDD method from EMPY)
//-------------------------------------------------------------------------
int CompleteEllipticIntegrals(double k2, double& K, double& E)
{
    // Hastings polynomial coefficients for K(k)
    static const double A0 = 1.38629436112;
    static const double A1 = 0.09666344259;
    static const double A2 = 0.03590092383;
    static const double A3 = 0.03742563713;
    static const double A4 = 0.01451196212;

    static const double B0 = 0.5;
    static const double B1 = 0.12498593597;
    static const double B2 = 0.06880248576;
    static const double B3 = 0.03328355346;
    static const double B4 = 0.00441787012;

    // Hastings polynomial coefficients for E(k)
    static const double C1 = 0.44325141463;
    static const double C2 = 0.06260601220;
    static const double C3 = 0.04757383546;
    static const double C4 = 0.01736506451;

    static const double D1 = 0.24998368310;
    static const double D2 = 0.09200180037;
    static const double D3 = 0.04069697526;
    static const double D4 = 0.00526449639;

    // Check validity: 0 <= k2 < 1
    if (k2 < 0.0 || k2 >= 1.0) {
        K = 0.0;
        E = 0.0;
        return 1;  // Error
    }

    // Compute using Hastings polynomial approximation
    double X1 = 1.0 - k2;  // complementary modulus squared
    double X2 = X1 * X1;
    double X3 = X1 * X2;
    double X4 = X2 * X2;
    double XI = 1.0 / X1;
    double XL = log(XI);

    // K(k) = A(X1) + B(X1) * ln(1/X1)
    double AA = A0 + A1 * X1 + A2 * X2 + A3 * X3 + A4 * X4;
    double BB = B0 + B1 * X1 + B2 * X2 + B3 * X3 + B4 * X4;
    K = AA + BB * XL;

    // E(k) = C(X1) + D(X1) * ln(1/X1)
    static const double C0 = 1.0;
    double CC = C0 + C1 * X1 + C2 * X2 + C3 * X3 + C4 * X4;
    double DD = D1 * X1 + D2 * X2 + D3 * X3 + D4 * X4;
    E = CC + DD * XL;

    return 0;  // Success
}

//-------------------------------------------------------------------------
// Perfect elliptic integral (from EMPY EllipticIntegral.cpp)
// Uses higher-precision polynomial approximation for k2 near 1
//-------------------------------------------------------------------------
int PerfectEllipticIntegral(double k2, double kd2, double* K, double* E)
{
    static const double EPS_k2 = 1.0e-6;
    static const double EPS_kd2 = 0.44408920985006e-15;

    static const double CK1[10] = {
        0.30072519903687e-03, 0.39684709020990e-02,
        0.10795990490592e-01, 0.10589953620989e-01, 0.75193867218084e-02,
        0.89266462945565e-02, 0.14942029142282e-01, 0.30885173001900e-01,
        0.96573590301742e-01, 0.13862943611199e+01
    };
    static const double CK2[10] = {
        0.66631752464607e-04, 0.17216147097987e-02,
        0.92811603829686e-02, 0.20690240005101e-01, 0.29503729348689e-01,
        0.37335546682286e-01, 0.48827155048118e-01, 0.70312495459547e-01,
        0.12499999999764e+00, 0.50000000000000e+00
    };
    static const double CE1[10] = {
        0.32519201550639e-03, 0.43025377747931e-02,
        0.11785841008734e-01, 0.11841925995501e-01, 0.90355277375409e-02,
        0.11716766944658e-01, 0.21836131405487e-01, 0.56805223329308e-01,
        0.44314718058337e+00, 0.0
    };
    static const double CE2[10] = {
        0.72031696345716e-04, 0.18645379184063e-02,
        0.10087958494375e-01, 0.22660309891604e-01, 0.32811069172721e-01,
        0.42672510126592e-01, 0.58592707184265e-01, 0.93749995116367e-01,
        0.24999999999746e+00, 0.0
    };

    // if k2 < EPS_k2, return 0 (K and E not meaningful)
    if (k2 < EPS_k2) {
        return 0;
    }

    double sumK1, sumK2, sumE1, sumE2;
    double log_kd2;

    if (kd2 > 0.0) {
        log_kd2 = log(kd2);
        if (kd2 >= EPS_kd2) {
            // Normal case: polynomial approximation
            sumK1 = CK1[0];
            sumK2 = CK2[0];
            for (int i = 1; i < 10; ++i) {
                sumK1 = sumK1 * kd2 + CK1[i];
                sumK2 = sumK2 * kd2 + CK2[i];
            }

            sumE2 = 0.0;
            sumE1 = sumE2;
            for (int i = 0; i < 9; ++i) {
                sumE1 = kd2 * (sumE1 + CE1[i]);
                sumE2 = kd2 * (sumE2 + CE2[i]);
            }

            *K = sumK1 - log_kd2 * sumK2;
            *E = 1.0 + (sumE1 - log_kd2 * sumE2);
        } else {
            // kd2 very small (k2 very close to 1)
            *K = CK1[9] - log_kd2 * CK2[9];
            *E = 1.0;
        }
    } else {
        // kd2 = 0 (k2 = 1): singular case
        *K = 13.0;  // Large value (actually infinite)
        *E = 1.0;
    }

    return 1;
}

//-------------------------------------------------------------------------
// Circular loop B-field using elliptic integrals (from EMPY bring0_)
//
// Radia now uses SI units (meters) internally, matching ELF.
// All coordinates (R, Z, CR, CZ) are in meters, CI is in Amperes.
// Output is in Tesla.
//
// Formula based on:
//   BR = (mu_0 * I / (2*pi)) * (z-zc) / (R * sqrt((R+CR)^2 + (z-zc)^2)) *
//        (-K + (R^2 + CR^2 + (z-zc)^2) / ((R-CR)^2 + (z-zc)^2) * E)
//   BZ = (mu_0 * I / (2*pi)) * 1 / sqrt((R+CR)^2 + (z-zc)^2) *
//        (K - (R^2 - CR^2 + (z-zc)^2) / ((R-CR)^2 + (z-zc)^2) * E)
//-------------------------------------------------------------------------
void CircularLoopBField(double R, double Z, double CR, double CZ, double CI,
                        double& BR, double& BZ)
{
    // SI formula: AA = mu_0/pi * I = 4e-7 * I [T*m/A * A = T*m]
    // Then AA / (2*sqrt(S+P)) gives [T] when S+P is in m^2
    double AA = 4.0e-7 * CI;

    // On-axis case (R = 0)
    double R_check = R + 1.0;
    if (R_check == 1.0) {  // R is essentially zero
        double CR2 = CR * CR;
        double dz = CZ - Z;
        double S = CR2 + dz * dz;
        S = sqrt(S) * S;  // S^(3/2)
        BR = 0.0;
        BZ = AA * PI / 2.0 / S * CR2;
        return;
    }

    // General case
    double dz = CZ - Z;
    double S = CR * CR + R * R + dz * dz;
    double P = 2.0 * CR * R;
    double RK2 = 2.0 * P / (S + P);

    double ELPK, ELPE;
    int IER;
    IER = CompleteEllipticIntegrals(RK2, ELPK, ELPE);

    if (IER != 0) {
        BR = 0.0;
        BZ = 0.0;
        return;
    }

    // S_mi_P = (R - CR)^2 + dz^2
    double S_mi_P = S - P;

    AA = AA / (2.0 * sqrt(S + P));
    BR = AA * (Z - CZ) * (-ELPK + S / S_mi_P * ELPE) / R;
    BZ = AA * (ELPK - (S - 2.0 * CR * CR) / S_mi_P * ELPE);
}

//-------------------------------------------------------------------------
// Circular loop vector potential A_phi using elliptic integrals (from EMPY psi_)
//
// IMPORTANT: This function works in Radia's internal units (mm).
// All coordinates (R, Z, CR, CZ) are in mm, CI is in Amperes.
// Output is in Tesla*mm (vector potential in Radia units).
//
// Formula:
//   A_phi = (mu_0 * I / (4*pi)) * sqrt(CR/R) * 2/k * ((1 - k^2/2)*K - E)
// where k^2 = 4*CR*R / ((CR+R)^2 + (z-zc)^2)
//
// Note: The function returns A_phi (vector potential component), NOT psi.
//       psi = R * A_phi (magnetic flux function)
//-------------------------------------------------------------------------
double CircularLoopAPhi(double R, double Z, double CR, double CZ, double CI)
{
    // On-axis case (R = 0): A_phi = 0 by symmetry
    double R_check = R + 1.0;
    if (R_check == 1.0) {
        return 0.0;
    }

    // SI formula for vector potential of circular loop:
    //   A_phi = (mu_0 * I / pi) * sqrt(CR/R) / k * ((1 - k^2/2)*K - E)
    // where k^2 = 4*CR*R / ((CR+R)^2 + dz^2)
    //
    // For coordinates in mm:
    // - CR, R, dz are all in mm
    // - A_phi in SI has units of T*m
    // - Radia expects A in T*mm (since coordinates are in mm)
    // - Since sqrt(CR/R) is dimensionless and k is dimensionless,
    //   the units come from mu_0 * I / pi which is in T*m
    // - To get T*mm output, multiply by 1e3
    //   BUT: the formula also involves sqrt(P) = sqrt(CR*R) which is in mm
    //   so we need to be careful about the actual formula used
    //
    // Using the flux function approach (from EMPY):
    //   psi = (mu_0 * I / pi) * sqrt(P) / k * ((1 - k^2/2)*K - E)
    //   where P = CR * R (in mm^2), so sqrt(P) is in mm
    // Radia now uses SI units (meters) internally, matching ELF.
    //
    // Derivation for SI units:
    //   A_phi(SI) [T*m] = (mu_0/pi) * I * sqrt(CR/R) * (1/k) * f(K,E)
    //   where f(K,E) = (1 - k^2/2)*K - E is dimensionless
    //   sqrt(CR/R) is dimensionless (ratio of radii)
    //   mu_0/pi = 4e-7 [T*m/A], I [A]
    //   So A_phi(SI) = 4e-7 * I * sqrt(CR/R) * (1/k) * f(K,E) [T*m]

    double UM = 4.0e-7 * CI;  // mu_0 * I / pi for SI units (T*m)
    double dz = CZ - Z;
    double S = CR * CR + R * R + dz * dz;
    double P = CR * R;
    double RK2 = 4.0 * P / (S + 2.0 * P);

    double ELPK, ELPE;
    int IER;
    IER = CompleteEllipticIntegrals(RK2, ELPK, ELPE);

    if (IER != 0) {
        return 0.0;
    }

    // A_phi = (mu_0 * I / pi) * sqrt(CR/R) * (1/k) * ((1 - k^2/2)*K - E)
    double RK = sqrt(RK2);
    double A_phi = UM * sqrt(CR / R) / RK * ((1.0 - RK2 / 2.0) * ELPK - ELPE);

    return A_phi;
}

//-------------------------------------------------------------------------
// Circular loop solid angle for magnetic scalar potential calculation
//
// Radia now uses SI units (meters) internally, matching ELF.
// All coordinates (R, Z, CR, CZ) are in meters.
// Output is in steradians [sr].
//
// The magnetic scalar potential is related to solid angle by:
//   Phi = I * Omega / (4*pi)   [A]
//
// References:
//   [1] J.D. Jackson, "Classical Electrodynamics," 3rd ed., Wiley, 1998,
//       Chapter 5.6.
//   [2] MIT OpenCourseWare, 6.013 Electromagnetics, Section 8.3,
//       "The Scalar Magnetic Potential"
//   [3] J.T. Conway, "Analytical solutions for the self-inductance of a
//       circular disk coil using complete elliptic integrals," J. Phys. A,
//       Vol. 34, pp. 3687-3695, 2001.
//
// Formula:
//   On-axis (R=0): Omega = 2*pi * (1 - z/sqrt(CR^2 + z^2)) * sign(z)
//   Off-axis:      Omega = sign(z) * (2*pi/k) * ((2-k^2)*K - 2*E)
//                  where k^2 = 4*CR*R / ((CR+R)^2 + z^2)
//
// Note: The solid angle has a discontinuity of 2*pi when crossing
//       the loop plane (z=0) inside the loop (R < CR).
//       This is physically correct - the scalar potential is multivalued.
//-------------------------------------------------------------------------
double CircularLoopSolidAngle(double R, double Z, double CR, double CZ)
{
    double dz = Z - CZ;  // z relative to coil center

    // Handle z = 0 case (on the loop plane)
    // Inside loop (R < CR): Omega = +-2*pi depending on approach direction
    // Outside loop (R > CR): Omega = 0
    // On loop (R = CR): Omega = +-pi (half the discontinuity)
    double dz_check = dz + 1.0;
    if (dz_check == 1.0) {  // dz is essentially zero
        // On the loop plane, return 0 (limit from either side averages to 0 outside)
        // Inside the loop, this is a singular/multivalued point
        // For numerical stability, return 0
        return 0.0;
    }

    // Sign of z determines the sign of the solid angle
    double sign_z = (dz > 0.0) ? 1.0 : -1.0;
    double z_abs = fabs(dz);

    // On-axis case (R = 0 or R << CR)
    // Use tolerance based on coil radius: R < CR * 1e-6 is "on-axis"
    double on_axis_tol = CR * 1.0e-6;
    if (on_axis_tol < 1.0e-10) on_axis_tol = 1.0e-10;  // Minimum tolerance

    if (R < on_axis_tol) {
        double dist = sqrt(CR * CR + dz * dz);
        // Omega = 2*pi * (1 - |z|/dist) * sign(z)
        // For z > 0 (above loop): positive solid angle (loop seen from front)
        // For z < 0 (below loop): negative solid angle (loop seen from back)
        double Omega = 2.0 * PI * (1.0 - z_abs / dist) * sign_z;
        return Omega;
    }

    // General off-axis case using elliptic integrals
    // k^2 = 4*CR*R / ((CR+R)^2 + z^2)
    double sum_r = CR + R;
    double denom = sum_r * sum_r + dz * dz;
    double k2 = 4.0 * CR * R / denom;

    // Avoid k2 >= 1 which would make elliptic integrals singular
    if (k2 >= 1.0) {
        k2 = 1.0 - 1.0e-10;
    }

    double ELPK, ELPE;
    int IER = CompleteEllipticIntegrals(k2, ELPK, ELPE);

    if (IER != 0) {
        return 0.0;
    }

    // Omega = sign(z) * (2*pi/k) * ((2-k^2)*K - 2*E)
    // This can also be written as:
    // Omega = sign(z) * (4*pi/k) * ((1-k^2/2)*K - E)
    //
    // Alternative form using Heuman's Lambda function:
    // Omega = 2*pi * (1 - Lambda0(beta, k)) * sign(z)
    // where sin(beta) = z / sqrt((R-CR)^2 + z^2)
    //
    // We use the direct elliptic integral form:
    double k = sqrt(k2);
    double Omega = sign_z * (2.0 * PI / k) * ((2.0 - k2) * ELPK - 2.0 * ELPE);

    return Omega;
}

//-------------------------------------------------------------------------
// Incomplete elliptic integrals F(phi, m) and E(phi, m)
// using arithmetic-geometric mean algorithm
//
// Reference:
//   [1] A. Kameari, EMPY implementation (ell_fe function)
//   [2] M. Abramowitz, I.A. Stegun, "Handbook of Mathematical Functions",
//       Chapter 17: Elliptic Integrals
//
// Parameters:
//   phi: amplitude (radians)
//   m: parameter m = k^2 (modulus squared)
//   F: output - incomplete elliptic integral of first kind
//   E: output - incomplete elliptic integral of second kind
//
// The integrals are defined as:
//   F(phi, m) = integral_0^phi 1/sqrt(1 - m*sin^2(t)) dt
//   E(phi, m) = integral_0^phi sqrt(1 - m*sin^2(t)) dt
//-------------------------------------------------------------------------
void IncompleteEllipticIntegrals(double phi, double m, double& F, double& E)
{
    const double HALF_PI = PI / 2.0;

    // Handle m >= 1 case (degenerate)
    if (m >= 1.0) {
        double s = std::sin(phi);
        F = 0.5 * std::log((1.0 + s) / (1.0 - s));
        int per = static_cast<int>(std::floor((phi + HALF_PI) / PI));
        double sign_val = (std::abs(per % 2) < 0.5) ? 1.0 : -1.0;
        E = std::sin(phi) * sign_val + 2.0 * per;
        return;
    }

    // Arithmetic-geometric mean algorithm
    double sgn = (phi >= 0.0) ? 1.0 : -1.0;
    phi = std::fabs(phi);

    double a = 1.0;
    double b = std::sqrt(1.0 - m);
    double twon = 1.0;
    double eok = 1.0 - 0.5 * m;
    double cs = 0.0;

    // AGM iteration (maximum ~8 passes for double precision)
    for (int iter = 0; iter < 20; ++iter) {
        double c = 0.5 * (a - b);
        if (c < 1.0e-15) break;

        double phase = (phi + HALF_PI) / PI;
        double cycle = std::floor(phase);
        if (std::fabs(cycle - phase) < 1.0e-10) {
            phi *= 1.0 + 1.0e-15;
        }
        phi += std::atan((b / a) * std::tan(phi)) + PI * cycle;

        cs += c * std::sin(phi);
        eok -= twon * c * c;

        twon *= 2.0;
        double am = a - c;
        if (am == a) break;
        b = std::sqrt(a * b);
        a = am;
    }

    F = phi / (twon * a) * sgn;
    // Bug fix: E formula was E = eok * F * sgn + sgn * cs (WRONG)
    // Correct formula: E = eok * F + sgn * cs
    // The extra "* sgn" after F was causing wrong sign for negative angles
    E = eok * F + sgn * cs;
}

//-------------------------------------------------------------------------
// Incomplete elliptic integrals between two angles
// Computes F(phi2,m) - F(phi1,m) and E(phi2,m) - E(phi1,m)
//
// Parameters:
//   phi1, phi2: amplitude range (radians)
//   m: parameter m = k^2 (modulus squared)
//   F: output - F(phi2,m) - F(phi1,m)
//   E: output - E(phi2,m) - E(phi1,m)
//-------------------------------------------------------------------------
void IncompleteEllipticIntegralsDiff(double phi1, double phi2, double m,
                                      double& F, double& E)
{
    double F2 = 0.0, E2 = 0.0;
    double F1 = 0.0, E1 = 0.0;

    if (phi2 != 0.0) {
        IncompleteEllipticIntegrals(phi2, m, F2, E2);
    }

    // Optimization for symmetric case
    if (phi2 == -phi1) {
        F = 2.0 * F2;
        E = 2.0 * E2;
        return;
    }

    if (phi1 != 0.0) {
        IncompleteEllipticIntegrals(phi1, m, F1, E1);
    }

    F = F2 - F1;
    E = E2 - E1;
}

//-------------------------------------------------------------------------
// Arc coil B-field using incomplete elliptic integrals (far-field approximation)
//
// Reference:
//   [1] A. Kameari, "Calculation of Transient 3D Eddy Current using
//       Edge-Elements", IEEE Trans. Magn., Vol. 26, No. 2, 1990
//   [2] A. Kameari, EMPY_Field implementation (field_ARC_line function)
//
// This implements the analytical formula for the magnetic field of an arc coil
// in the far-field regime (when observation point is far from coil cross-section).
// The coil is treated as a line current along the arc.
//
// Parameters:
//   rho: radial distance from axis in local frame (m)
//   z: axial distance from coil plane in local frame (m)
//   rc: coil radius (m)
//   phi1, phi2: arc start and end angles relative to observation point (rad)
//   current: total current (A)
//   cross_area: cross-section area (m^2)
//   B_rho, B_theta, B_z: output B-field components in cylindrical coords (T)
//
// The observation point is at (rho, 0, z) in the local frame where phi=0
// points toward the observation point's azimuthal position.
//-------------------------------------------------------------------------
void ArcCoilBFieldAnalytical(double rho, double z, double rc,
                              double phi1, double phi2,
                              double current, double cross_area,
                              double& B_rho, double& B_theta, double& B_z)
{
    const double ZOT = 1.0e-10;

    // Current density factor: J * mu_0 / (4*pi)
    double J_density = current / cross_area;
    double XI = J_density * MU0_DIV_4PI;

    //=========================================================================
    // Special case: On-axis observation point (rho = 0 or very small)
    //
    // When rho is small (e.g., due to SmallPositive offset in the caller),
    // the standard elliptic integral formula degenerates. Use the direct
    // on-axis analytical formula instead.
    //
    // For on-axis points, B_z is proportional to the arc angle span.
    // This is the key fix for partial arc calculations.
    //=========================================================================
    // Tolerance accounts for SmallPositive offset in caller (sqrt(1e-10) ≈ 1e-5)
    // and small numerical values. Use absolute tolerance to catch on-axis cases.
    const double RHO_TOLERANCE = 1.0e-4;  // 0.1mm absolute tolerance
    if (rho < RHO_TOLERANCE) {
        // For on-axis points, use direct formula
        // B_z = (mu_0 * I * R^2 / (2 * (R^2 + z^2)^1.5)) * (dphi / 2*pi)
        // But we have cross_area and current density, so:
        // B_z = XI * cross_area * rc^2 * dphi / (R^2 + z^2)^1.5
        //     where dphi = phi2 - phi1

        double ZSQ = z * z;
        double RAZSQ = rc * rc + ZSQ;  // rc^2 + z^2
        double RAZ = std::sqrt(RAZSQ);

        double dphi = phi2 - phi1;  // Arc angle span

        // On-axis B_z formula (from Biot-Savart for arc)
        B_z = XI * cross_area * rc * rc * dphi / (RAZSQ * RAZ);

        // B_rho and B_theta: Use angular dependence
        double SINA = std::sin(phi2) - std::sin(phi1);
        double COSA = std::cos(phi1) - std::cos(phi2);

        // For on-axis points, B_rho and B_theta are small and depend on angular terms
        // These represent the transverse field components that cancel for a full circle
        double BRT = XI * cross_area * rc / (RAZSQ * RAZ);
        B_rho = BRT * SINA;
        B_theta = BRT * COSA;

        return;
    }

    // Standard off-axis calculation using elliptic integrals
    double RHOSQ = rho * rho;
    double ZSQ = z * z;
    double RZSQ = RHOSQ + ZSQ;
    double ASQ = rc * rc;
    double ASRZS = ASQ + RZSQ;
    double APR = rc + rho;
    double APRSQ = APR * APR;
    double RAZSQ = APRSQ + ZSQ;
    double RAZ = std::sqrt(RAZSQ);

    // Quantities for far-field approximation (line current)
    double ARZSQ = (rc - rho) * (rc - rho) + ZSQ;  // min distance squared
    double TAR = 2.0 * rho * rc;

    // Modulus for elliptic integrals
    // CASQ = k^2 = 4*rho*rc / ((rho+rc)^2 + z^2)
    double CASQ = TAR + TAR;
    if (RAZSQ > 0.0) {
        CASQ /= RAZSQ;
    }

    // Angular quantities (shifted by pi for the formula)
    double OMEGA1 = 0.5 * (phi1 - PI);
    double OMEGA2 = 0.5 * (phi2 - PI);

    double SNOMG1 = std::sin(OMEGA1);
    double CSOMG1 = std::cos(OMEGA1);
    double SNOMG2 = std::sin(OMEGA2);
    double CSOMG2 = std::cos(OMEGA2);
    double SN1SQ = SNOMG1 * SNOMG1;
    double SN2SQ = SNOMG2 * SNOMG2;

    // Compute incomplete elliptic integrals
    double FEI = 0.0, EEI = 0.0;
    IncompleteEllipticIntegralsDiff(OMEGA1, OMEGA2, CASQ, FEI, EEI);

    B_rho = 0.0;
    B_theta = 0.0;
    B_z = 0.0;

    static const double EPS_k2 = 1.0e-6;

    if (CASQ >= EPS_k2) {
        // Normal case: use elliptic integral formula
        double D1 = std::sqrt(1.0 - CASQ * SN1SQ);
        double D2 = std::sqrt(1.0 - CASQ * SN2SQ);

        // B-field calculation
        double D1_inv = 1.0 / D1;
        double D2_inv = 1.0 / D2;
        double SNCS1 = SNOMG1 * CSOMG1;
        double SNCS2 = SNOMG2 * CSOMG2;
        double SL = EEI - CASQ * (SNCS2 * D2_inv - SNCS1 * D1_inv);

        double C = cross_area / RAZ;

        if (std::fabs(z) > ZOT) {
            double D = (z / rho) * C;
            B_theta = XI * D * (D1_inv - D2_inv);
            B_rho = XI * D * (SL * (ASRZS / ARZSQ) - FEI);
        }

        if (std::fabs(1.0 - RHOSQ / ASQ) > ZOT) {
            SL *= (RZSQ - ASQ) / ARZSQ;
        }
        B_z = XI * C * (FEI - SL);

    } else {
        // Small modulus case: simplified formula
        double C = (2.0 * rc / RAZ) * cross_area;
        double ARO = C * (SN2SQ - SN1SQ);

        B_z = XI * (C / RAZSQ) * APR * EEI;
        double G = z / RAZSQ;
        B_theta = -XI * G * ARO;
    }
}

//-------------------------------------------------------------------------
// Arc coil field computation with distance-based method selection
//
// This is the main entry point for arc coil field calculation.
// It automatically selects between:
//   1. Analytical (line current) formula for far-field
//   2. Numerical integration for near-field
//
// Parameters:
//   obs_rho, obs_z: observation point in local cylindrical coords (m)
//   rc: coil center radius (m)
//   phi1, phi2: arc angles (rad)
//   a: half radial width (m)
//   b: half axial width (m)
//   current: total current (A)
//   B_rho, B_theta, B_z: output B-field components (T)
//
// Returns: true if analytical formula was used, false if numerical integration
//-------------------------------------------------------------------------
bool ArcCoilBFieldAdaptive(double obs_rho, double obs_z, double rc,
                            double phi1, double phi2,
                            double a, double b, double current,
                            double& B_rho, double& B_theta, double& B_z)
{
    // Cross-section area
    double cross_area = 4.0 * a * b;

    // Minimum distance to coil cross-section corners
    double ARZSQ = (rc - obs_rho) * (rc - obs_rho) + obs_z * obs_z;
    double TAR = 2.0 * obs_rho * rc;

    // Check if phi range includes phi=0 (closest point)
    double DMIN = ARZSQ;
    if (phi1 * phi2 > 0.0) {
        // Arc doesn't cross phi=0, minimum distance is at arc endpoints
        double cos_max = std::max(std::cos(phi1), std::cos(phi2));
        DMIN += TAR * (1.0 - cos_max);
    }

    // Distance threshold for line current approximation
    // Use analytical formula when distance > LINE_APP_LIMIT * max(a, b)
    // Reduced from 3.0 to 2.0 to use more accurate far-field formula more often
    // The near-field formula (Kameari) has numerical issues for large aspect ratio coils
    const double LINE_APP_LIMIT = 2.0;  // Was 3.0, reduced to workaround near-field bugs
    double DFT = LINE_APP_LIMIT * std::max(a, b);
    DFT = DFT * DFT;

    if (DMIN >= DFT) {
        // Far-field: use analytical (line current) formula
        ArcCoilBFieldAnalytical(obs_rho, obs_z, rc, phi1, phi2,
                                 current, cross_area,
                                 B_rho, B_theta, B_z);
        return true;
    }

    // Near-field: return false to indicate numerical integration is needed
    B_rho = 0.0;
    B_theta = 0.0;
    B_z = 0.0;
    return false;
}

//-------------------------------------------------------------------------
// Helper function: ATAN4 - computes arctan combination for field integrals
// Reference: EMPY_Field (used in arc_field_integrand)
//-------------------------------------------------------------------------
static double ATAN4(double PP0, double PP1, double PP2, double PP3)
{
    return std::atan(PP2) - std::atan(PP0) + std::atan(PP3) - std::atan(PP1);
}

//-------------------------------------------------------------------------
// Safe log division: LOGDIV(A, B) = log(A/B) with zero handling
//-------------------------------------------------------------------------
static double LOGDIV(double A, double B)
{
    const double FEPS = 1.0e-50;
    if (std::fabs(A) > FEPS && std::fabs(B) > FEPS) {
        return std::log(A / B);
    }
    return 0.0;
}

//-------------------------------------------------------------------------
// Arc coil field integrand with analytical cross-section integration
//
// References:
//   [1] A. Kameari, "Calculation of Transient 3D Eddy Current using
//       Edge-Elements", IEEE Trans. Magn., Vol. 26, No. 2, pp. 466-469, 1990.
//       - Analytical cross-section integration using distances to 4 corners
//   [2] A. Kameari, EMPY_Field library (arc_field_integrand function)
//       - Original implementation of this algorithm
//
// For each azimuthal angle ALPHA (phi), this function computes the
// contribution to B and A fields by analytically integrating over the
// rectangular cross-section using distances to 4 corners (R[0-3]).
//
// Parameters (in param structure):
//   RHO: observation point radial distance from axis
//   AL, AU: inner/outer coil radius (rc-a, rc+a)
//   ZZL, ZZU: axial distances (z-b, z+b)
//   ARZ[4]: squared distances to 4 corners at phi=0
//   ARL, ARU: 2*rho*(rc-a), 2*rho*(rc+a)
//   ALSQ, AUSQ, ZLSQ, ZUSQ: squared terms
//   calcA, calcB: flags for A or B calculation
//   symmetricArc: 1 if using symmetry (only positive half)
//
// Output F array:
//   When symmetricArc=0 and calcB=1: F[0]=B_rho term, F[1]=B_theta term, F[2]=B_z term
//   When symmetricArc=1 and calcB=1: F[0]=B_rho term (symmetric), F[1]=B_z term
//   Similar pattern for A field terms
//-------------------------------------------------------------------------
void ArcFieldIntegrand(double ALPHA, double* F, const ArcFieldParams* param)
{
    double rp = param->RHO;
    double AL = param->AL;      // rc - a (inner radius)
    double AU = param->AU;      // rc + a (outer radius)
    double ARL = param->ARL;    // 2*rp*(rc-a)
    double ARU = param->ARU;    // 2*rp*(rc+a)
    double ALSQ = param->ALSQ;  // (rc-a-rp)^2
    double AUSQ = param->AUSQ;  // (rc+a-rp)^2
    double ZLSQ = param->ZLSQ;  // (z+b)^2
    double ZUSQ = param->ZUSQ;  // (z-b)^2
    double ZZL = param->ZZL;    // z + b
    double ZZU = param->ZZU;    // z - b

    double sinPhi = std::sin(ALPHA);
    double cosPhi = std::cos(ALPHA);
    double RSA = rp * sinPhi;
    double rpCosPhi = rp * cosPhi;

    // Geometric terms
    double GXL = AL - rpCosPhi;      // rc-a - rp*cos(phi)
    double GXU = AU - rpCosPhi;      // rc+a - rp*cos(phi)

    // (1 - cos(phi)) term - use Taylor expansion for small angles
    double COSA1 = 1.0 - cosPhi;
    if (COSA1 == 0.0) {
        COSA1 = ALPHA * ALPHA / 2.0;
    }

    double ARCL = ARL * COSA1;       // 2*rp*(rc-a)*(1-cos)
    double ARCU = ARU * COSA1;       // 2*rp*(rc+a)*(1-cos)
    double CLSQ = ARCL + ALSQ;       // (rc-a)^2 + rp^2 - 2*(rc-a)*rp*cos
    double CUSQ = ARCU + AUSQ;       // (rc+a)^2 + rp^2 - 2*(rc+a)*rp*cos

    // Distances from observation point to 4 corners of cross-section
    // R[i] = sqrt((r-rp)^2 + 2*rp*r*(1-cos) + (zp-z)^2) for corners
    double R[4];
    R[0] = std::sqrt(param->ARZ[0] + ARCL);  // Corner: r=rc-a, z=zc-b
    R[1] = std::sqrt(param->ARZ[1] + ARCL);  // Corner: r=rc-a, z=zc+b
    R[2] = std::sqrt(param->ARZ[2] + ARCU);  // Corner: r=rc+a, z=zc-b
    R[3] = std::sqrt(param->ARZ[3] + ARCU);  // Corner: r=rc+a, z=zc+b

    double RL = R[2] - R[0];  // R(outer,lower) - R(inner,lower)
    double RU = R[3] - R[1];  // R(outer,upper) - R(inner,upper)

    // Log terms for radial integration
    double AXU = ((GXU * R[3] + CUSQ) + ZUSQ) / R[3];
    double AXL = ((GXL * R[1] + CLSQ) + ZUSQ) / R[1];
    double ALXU = std::log(AXU / AXL);

    AXU = ((GXU * R[2] + CUSQ) + ZLSQ) / R[2];
    AXL = ((GXL * R[0] + CLSQ) + ZLSQ) / R[0];
    double ALXL = std::log(AXU / AXL);

    double ZAZA = ZZU * ALXU - ZZL * ALXL;

    // Additional log terms for B_z
    double RZU = ((ZUSQ + ZZU * R[1]) + CLSQ) / R[1];
    double RZL = ((ZLSQ + ZZL * R[0]) + CLSQ) / R[0];
    double RLZ = std::fabs(RZU / RZL);

    RZU = ((ZUSQ + ZZU * R[3]) + CUSQ) / R[3];
    RZL = ((ZLSQ + ZZL * R[2]) + CUSQ) / R[2];
    double RUZ = std::fabs(RZU / RZL);

    int index = 0;

    // B-field calculation
    if (param->calcB) {
        double PP[4];
        PP[0] = (ZZL * GXL) / (RSA * R[0]);
        PP[1] = (ZZU * GXL) / (RSA * R[1]);
        PP[2] = (ZZL * GXU) / (RSA * R[2]);
        PP[3] = (ZZU * GXU) / (RSA * R[3]);

        double BRT = (RU - RL) + rpCosPhi * (ALXU - ALXL);
        double F3 = RSA * ATAN4(PP[0], PP[1], PP[2], PP[3]) - rpCosPhi * LOGDIV(RLZ, RUZ) - ZAZA;

        F[index++] = cosPhi * BRT;  // B_rho term

        if (param->symmetricArc == 0) {
            F[index++] = sinPhi * BRT;  // B_theta term
            F[index++] = F3;             // B_z term
        } else {
            F[index++] = F3;  // B_z term (B_theta is antisymmetric, cancels)
        }
    }

    // A-field calculation (if needed)
    if (param->calcA) {
        double ZLRS = ZZL * RSA;
        double ZURS = ZZU * RSA;
        double PP[4];

        if (ZLRS != 0.0) {
            PP[0] = (CLSQ + GXL * R[0]) / ZLRS;
            PP[2] = (CUSQ + GXU * R[2]) / ZLRS;
        } else {
            PP[0] = 0.0;
            PP[2] = 0.0;
        }

        if (ZURS != 0.0) {
            PP[1] = (CLSQ + GXL * R[1]) / ZURS;
            PP[3] = (CUSQ + GXU * R[3]) / ZURS;
        } else {
            PP[1] = 0.0;
            PP[3] = 0.0;
        }

        double ART = 0.5 * (ZZU * RU - ZZL * RL)
                   + rpCosPhi * (ZAZA - RSA * ATAN4(PP[0], PP[1], PP[2], PP[3]))
                   + ((rpCosPhi * GXU + 0.5 * CUSQ) * std::log(RUZ)
                    - (rpCosPhi * GXL + 0.5 * CLSQ) * std::log(RLZ));

        if (param->symmetricArc == 0) {
            F[index++] = sinPhi * ART;   // A_rho term
            F[index++] = cosPhi * ART;   // A_theta term
        } else {
            F[index++] = cosPhi * ART;   // A_theta term (A_rho is antisymmetric)
        }
    }
}

//-------------------------------------------------------------------------
// Adaptive 1D integration using Gauss-Kronrod quadrature with depth limit
//
// Reference:
//   [1] P. Wynn, "On a Device for Computing the e_m(S_n) Transformation",
//       MTAC 10, pp. 91-96, 1956.
//   [2] R. Piessens et al., QUADPACK: A Subroutine Package for Automatic
//       Integration, Springer-Verlag, 1983.
//   [3] A. Kameari, EMPY_Field library - bounded integration algorithm
//
// Parameters:
//   ncal: number of function values to compute
//   func: integrand function
//   a, b: integration bounds
//   eps: relative tolerance
//   result: output array of size ncal
//   param: parameter structure passed to func
//   depth: current recursion depth (default 0)
//
// Returns: number of function evaluations
//
// IMPORTANT: Recursion depth is bounded to KMAX (7) to prevent exponential
// growth of function evaluations when convergence is slow. This matches
// EMPY_Field behavior where maximum evaluations is O(2^KMAX) = O(128).
//-------------------------------------------------------------------------
static int Adaptive1DIntegrationImpl(int ncal,
                                      void (*func)(double, double*, const ArcFieldParams*),
                                      double a, double b, double eps,
                                      double* result, const ArcFieldParams* param,
                                      int depth)
{
    const int KMAX = 7;  // Maximum recursion depth (matches EMPY kmax1D)
    const double EPSINT = 1.0e-10;

    // Gauss-Kronrod 7-15 rule nodes and weights
    static const double xgk[8] = {
        0.0,
        0.40584515137739716690660641,
        0.74153118559939443986386477,
        0.94910791234275852452618968,
        0.20778495500789846760068940,
        0.58608723546769113507090340,
        0.86486442335976907278971278,
        0.99145537112081263920685390
    };

    static const double wgk[8] = {
        0.20948214108472782801299917,
        0.20443294007529889241416199,
        0.16900472663926790282658342,
        0.10479001032225018383987633,
        0.20443294007529889241416199,
        0.16900472663926790282658342,
        0.10479001032225018383987633,
        0.02293532201052922496373200
    };

    static const double wg[4] = {
        0.41795918367346938775510204,
        0.38183005050511894495036978,
        0.27970539148927666790146777,
        0.12948496616886969327061143
    };

    // Initialize results
    for (int i = 0; i < ncal; ++i) {
        result[i] = 0.0;
    }

    // Integration bounds
    double c = 0.5 * (a + b);
    double h = 0.5 * (b - a);

    if (std::fabs(h) < EPSINT) {
        return 0;
    }

    // Function values at Gauss-Kronrod points
    double fv[8 * 5];  // Max 5 function components, 8 points

    // Evaluate at center
    func(c, &fv[0], param);

    int neval = 1;

    // Evaluate at other points (symmetric)
    for (int j = 1; j < 8; ++j) {
        double x1 = c - h * xgk[j];
        double x2 = c + h * xgk[j];

        double f1[5], f2[5];
        func(x1, f1, param);
        func(x2, f2, param);
        neval += 2;

        for (int i = 0; i < ncal; ++i) {
            fv[j * ncal + i] = f1[i] + f2[i];
        }
    }

    // Compute Gauss-Kronrod and Gauss estimates
    double resK[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double resG[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

    for (int i = 0; i < ncal; ++i) {
        resK[i] = fv[0 * ncal + i] * wgk[0];
        resG[i] = fv[0 * ncal + i] * wg[0];

        for (int j = 1; j < 4; ++j) {
            resK[i] += fv[j * ncal + i] * wgk[j];
            resG[i] += fv[j * ncal + i] * wg[j];
        }
        for (int j = 4; j < 8; ++j) {
            resK[i] += fv[j * ncal + i] * wgk[j];
        }

        resK[i] *= h;
        resG[i] *= h;
    }

    // Estimate error
    double err_est = 0.0;
    for (int i = 0; i < ncal; ++i) {
        double diff = std::fabs(resK[i] - resG[i]);
        if (diff > err_est) err_est = diff;
    }

    double abs_est = 0.0;
    for (int i = 0; i < ncal; ++i) {
        if (std::fabs(resK[i]) > abs_est) abs_est = std::fabs(resK[i]);
    }

    // Check convergence
    if (err_est < eps * abs_est || err_est < 1.0e-14) {
        for (int i = 0; i < ncal; ++i) {
            result[i] = resK[i];
        }
        return neval;
    }

    // Check recursion depth limit (CRITICAL: prevents exponential blowup)
    if (depth >= KMAX) {
        // Maximum depth reached - return best estimate
        for (int i = 0; i < ncal; ++i) {
            result[i] = resK[i];
        }
        return neval;
    }

    // Recursive subdivision with depth tracking
    double mid = c;
    double res_left[5], res_right[5];

    int n_left = Adaptive1DIntegrationImpl(ncal, func, a, mid, eps, res_left, param, depth + 1);
    int n_right = Adaptive1DIntegrationImpl(ncal, func, mid, b, eps, res_right, param, depth + 1);

    for (int i = 0; i < ncal; ++i) {
        result[i] = res_left[i] + res_right[i];
    }

    return neval + n_left + n_right;
}

// Public interface wrapper
int Adaptive1DIntegration(int ncal,
                          void (*func)(double, double*, const ArcFieldParams*),
                          double a, double b, double eps,
                          double* result, const ArcFieldParams* param)
{
    return Adaptive1DIntegrationImpl(ncal, func, a, b, eps, result, param, 0);
}

//-------------------------------------------------------------------------
// Arc coil near-field integration using analytical cross-section formula
// with adaptive 1D integration over azimuthal direction
//
// References:
//   [1] A. Kameari, "Calculation of Transient 3D Eddy Current using
//       Edge-Elements", IEEE Trans. Magn., Vol. 26, No. 2, pp. 466-469, 1990.
//       - Analytical cross-section integration formula
//   [2] R. Piessens, E. de Doncker-Kapenga, C.W. Ueberhuber, D.K. Kahaner,
//       "QUADPACK: A Subroutine Package for Automatic Integration,"
//       Springer-Verlag, 1983. - Gauss-Kronrod 7-15 adaptive quadrature
//   [3] A. Kameari, EMPY_Field library (arc_integration function)
//       - Original implementation of this algorithm
//
// This computes the B-field of an arc coil when the observation point is
// near the coil cross-section. Instead of 2D Gauss quadrature over the
// cross-section, it uses:
//   1. Analytical formula for cross-section integral at each phi [1]
//   2. Adaptive Gauss-Kronrod 7-15 integration over azimuth [2]
//
// Parameters:
//   obs_rho: observation point radial distance from axis (m)
//   obs_z: observation point axial distance from coil center (m)
//   rc: coil center radius (m)
//   a: half radial width (m)
//   b: half axial width (m)
//   phi1, phi2: arc angles (rad)
//   current: total current (A)
//   B_rho, B_theta, B_z: output B-field components (T)
//-------------------------------------------------------------------------
void ArcCoilBFieldNearField(double obs_rho, double obs_z, double rc,
                             double a, double b,
                             double phi1, double phi2,
                             double current,
                             double& B_rho, double& B_theta, double& B_z)
{
    // Set up parameter structure for integrand
    ArcFieldParams param;

    double cross_area = 4.0 * a * b;
    double J_density = current / cross_area;
    double XI = J_density * MU0_DIV_4PI;

    // Coil geometry
    double AL = rc - a;    // Inner radius
    double AU = rc + a;    // Outer radius
    double ZZL = obs_z + b;  // z + b (lower axial boundary relative to obs)
    double ZZU = obs_z - b;  // z - b (upper axial boundary relative to obs)

    //=========================================================================
    // Special case: On-axis observation point (rho = 0)
    //
    // When rho=0, the standard integrand has division by RSA = rho*sin(phi) = 0.
    // Use direct analytical formula instead (from EMPY_Field EM_COIL_Arc.cpp).
    //
    // For on-axis points, the azimuthal integration gives a simple result
    // because the geometry is axially symmetric. The field only depends on
    // the arc angle span (phi2 - phi1).
    //=========================================================================
    // Tolerance accounts for SmallPositive offset in caller (sqrt(1e-10) ≈ 1e-5)
    const double RHO_TOLERANCE = 1.0e-4;  // 0.1mm absolute tolerance
    if (obs_rho < RHO_TOLERANCE) {
        // For rho=0, ARZ values simplify (no rho terms)
        double ALSQ = AL * AL;  // Inner radius squared
        double AUSQ = AU * AU;  // Outer radius squared
        double ZLSQ = ZZL * ZZL;
        double ZUSQ = ZZU * ZZU;

        // Distances from on-axis point to 4 corners of cross-section
        double R0 = std::sqrt(ALSQ + ZLSQ);  // r=AL, z=+b
        double R1 = std::sqrt(ALSQ + ZUSQ);  // r=AL, z=-b
        double R2 = std::sqrt(AUSQ + ZLSQ);  // r=AU, z=+b
        double R3 = std::sqrt(AUSQ + ZUSQ);  // r=AU, z=-b

        double RL = R2 - R0;  // Outer - Inner at z=+b
        double RU = R3 - R1;  // Outer - Inner at z=-b

        // Log terms for radial integration (on-axis simplification)
        double ALXU = std::log((AU + R3) / (AL + R1));
        double ALXL = std::log((AU + R2) / (AL + R0));

        // B-field components for on-axis point
        // B_rho and B_theta depend on angular variation
        double SINA = std::sin(phi2) - std::sin(phi1);
        double COSA = std::cos(phi1) - std::cos(phi2);
        double BRT = RU - RL;

        B_rho = XI * SINA * BRT;
        B_theta = XI * COSA * BRT;

        // B_z uses the analytical formula (key fix for partial arcs)
        double BZWAVE = (phi2 - phi1) * (ZZL * ALXL - ZZU * ALXU);
        B_z = XI * BZWAVE;

        return;
    }

    // Standard off-axis calculation using adaptive integration
    param.RHO = obs_rho;
    param.RHOSQ = obs_rho * obs_rho;
    param.ZSQ = obs_z * obs_z;
    param.ASQ = rc * rc;
    param.AA = rc;
    param.Z1 = obs_z;

    param.AL = AL;
    param.AU = AU;
    param.ZZL = ZZL;
    param.ZZU = ZZU;

    // Squared distances
    param.ALSQ = (param.AL - obs_rho) * (param.AL - obs_rho);
    param.AUSQ = (param.AU - obs_rho) * (param.AU - obs_rho);
    param.ZLSQ = param.ZZL * param.ZZL;
    param.ZUSQ = param.ZZU * param.ZZU;

    // Corner distances at phi=0
    param.ARZ[0] = param.ALSQ + param.ZLSQ;  // (rc-a-rp)^2 + (z+b)^2
    param.ARZ[1] = param.ALSQ + param.ZUSQ;  // (rc-a-rp)^2 + (z-b)^2
    param.ARZ[2] = param.AUSQ + param.ZLSQ;  // (rc+a-rp)^2 + (z+b)^2
    param.ARZ[3] = param.AUSQ + param.ZUSQ;  // (rc+a-rp)^2 + (z-b)^2

    // Factors for azimuthal variation
    param.ARL = 2.0 * obs_rho * param.AL;  // 2*rp*(rc-a)
    param.ARU = 2.0 * obs_rho * param.AU;  // 2*rp*(rc+a)

    // Other geometric quantities
    param.ARZSQ = (rc - obs_rho) * (rc - obs_rho) + obs_z * obs_z;
    param.TAR = 2.0 * obs_rho * rc;
    param.AREA = cross_area;

    param.calcA = 0;
    param.calcB = 1;
    param.symmetricArc = 0;

    // Integration tolerance
    const double ArcIntegralAccuracy = 1.0e-6;

    // Handle symmetric case for efficiency
    if (phi1 * phi2 >= 0.0) {
        // Arc doesn't cross phi=0
        double result[5];
        Adaptive1DIntegration(3, ArcFieldIntegrand, phi1, phi2, ArcIntegralAccuracy, result, &param);

        B_rho = XI * result[0];
        B_theta = XI * result[1];
        B_z = XI * result[2];
    } else {
        // Arc crosses phi=0 - split integration for better accuracy
        double result1[5], result2[5];

        // First part: phi1 to 0
        if (phi1 < -0.035) {  // About 2 degrees
            Adaptive1DIntegration(3, ArcFieldIntegrand, phi1, -0.035, ArcIntegralAccuracy, result1, &param);
            double result1b[5];
            Adaptive1DIntegration(3, ArcFieldIntegrand, -0.035, 0.0, ArcIntegralAccuracy, result1b, &param);
            for (int i = 0; i < 3; ++i) result1[i] += result1b[i];
        } else {
            Adaptive1DIntegration(3, ArcFieldIntegrand, phi1, 0.0, ArcIntegralAccuracy, result1, &param);
        }

        // Second part: 0 to phi2
        if (phi2 > 0.035) {
            Adaptive1DIntegration(3, ArcFieldIntegrand, 0.0, 0.035, ArcIntegralAccuracy, result2, &param);
            double result2b[5];
            Adaptive1DIntegration(3, ArcFieldIntegrand, 0.035, phi2, ArcIntegralAccuracy, result2b, &param);
            for (int i = 0; i < 3; ++i) result2[i] += result2b[i];
        } else {
            Adaptive1DIntegration(3, ArcFieldIntegrand, 0.0, phi2, ArcIntegralAccuracy, result2, &param);
        }

        B_rho = XI * (result1[0] + result2[0]);
        B_theta = XI * (result1[1] + result2[1]);
        B_z = XI * (result1[2] + result2[2]);
    }
}

} // namespace RadElliptic
