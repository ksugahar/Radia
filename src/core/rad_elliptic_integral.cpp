/*-------------------------------------------------------------------------
*
* File name:      rad_elliptic_integral.cpp
*
* Project:        RADIA
*
* Description:    Elliptic integral functions for analytical coil field
*                 computation (K(k), E(k) complete elliptic integrals)
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
*
* Author(s):      A. Kameari (original EMPY), adapted for Radia
*
* First release:  2025
*
-------------------------------------------------------------------------*/

#include "rad_elliptic_integral.h"
#include <cmath>

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

} // namespace RadElliptic
