/*-------------------------------------------------------------------------
*
* File name:      rad_elliptic_integral.cpp
*
* Project:        RADIA
*
* Description:    Elliptic integral functions for analytical coil field
*                 computation (K(k), E(k) complete elliptic integrals).
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
*   [5] M. Abramowitz, I.A. Stegun, "Handbook of Mathematical Functions,"
*       National Bureau of Standards, Chapter 17: Elliptic Integrals, 1964.
*
* Author(s):      Radia Development Team
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
// Radia now uses SI units (meters) internally, matching ELF.
// All coordinates (R, Z, CR, CZ) are in meters, CI is in Amperes.
// Output is in Tesla*m (vector potential in SI units).
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
// All coordinates (R, Z, CR, CZ) are in meters (SI).
// Output is in steradians [sr].
//
// The magnetic scalar potential is related to solid angle by:
//   Phi = I * Omega / (4*pi)   [A]
//
// Method: Omega = integral over disk bounded by loop of z*dA/|r|^3.
//   Inner integral (radial) is computed analytically.
//   Outer integral (azimuthal) uses 32-point Gauss-Legendre quadrature.
//
// Note: The solid angle has a discontinuity of 4*pi when crossing
//       the loop plane (z=0) inside the loop (R < CR).
//
// References:
//   [1] J.D. Jackson, "Classical Electrodynamics," 3rd ed., Ch. 5.6
//   [2] J.T. Conway, J. Phys. A, Vol. 34, pp. 3687-3695, 2001.
//-------------------------------------------------------------------------
double CircularLoopSolidAngle(double R, double Z, double CR, double CZ)
{
    double dz = Z - CZ;

    if (fabs(dz) < 1.0e-15) {
        return 0.0;
    }

    // On-axis analytical formula
    double on_axis_tol = CR * 1.0e-8;
    if (on_axis_tol < 1.0e-15) on_axis_tol = 1.0e-15;

    if (R < on_axis_tol) {
        double dist = sqrt(CR * CR + dz * dz);
        double sign_z = (dz > 0.0) ? 1.0 : -1.0;
        return sign_z * 2.0 * PI * (1.0 - fabs(dz) / dist);
    }

    // General off-axis: numerical azimuthal quadrature with analytical radial integral
    //
    // Omega = integral_0^{2*pi} I(phi) dphi
    // where I(phi) = z * integral_0^CR r_d/(r_d^2 - 2*r_d*R*cos(phi) + R^2 + z^2)^{3/2} dr_d
    //
    // Antiderivative: integral r/(r^2+Br+C)^{3/2} dr = -(C + Br/2) / (D * sqrt(r^2+Br+C))
    // where B = -2*R*cos(phi), C = R^2+z^2, D = C - B^2/4 = R^2*sin^2(phi) + z^2
    //
    // I(phi) = dz * [-(C - CR*R*cos(phi))/(D*R_a) + C/(D*R_0)]
    //        = dz/D * [C/R_0 - (C - CR*R*cos(phi))/R_a]
    // where R_a = sqrt(CR^2 - 2*CR*R*cos(phi) + R^2 + z^2)
    //       R_0 = sqrt(R^2 + z^2)
    //       C   = R^2 + z^2 = R_0^2

    // 32-point Gauss-Legendre quadrature on [0, 2*pi]
    static const int N_GAUSS = 32;
    static const double gp32[] = {
        -0.9972638618494816, -0.9856115115452684, -0.9647622555875064, -0.9349060759377397,
        -0.8963211557660521, -0.8493676137325700, -0.7944837959679424, -0.7321821187402897,
        -0.6630442669302152, -0.5877157572407623, -0.5068999089322294, -0.4213512761306353,
        -0.3318686022821276, -0.2392873622521371, -0.1444719615827965, -0.0483076656877383,
         0.0483076656877383,  0.1444719615827965,  0.2392873622521371,  0.3318686022821276,
         0.4213512761306353,  0.5068999089322294,  0.5877157572407623,  0.6630442669302152,
         0.7321821187402897,  0.7944837959679424,  0.8493676137325700,  0.8963211557660521,
         0.9349060759377397,  0.9647622555875064,  0.9856115115452684,  0.9972638618494816
    };
    static const double gw32[] = {
        0.0070186100094701,  0.0162743947309057,  0.0253920653092621,  0.0342738629130214,
        0.0428358980222267,  0.0509980592623762,  0.0586840934785355,  0.0658222227763618,
        0.0723457941088485,  0.0781938957870703,  0.0833119242269467,  0.0876520930044038,
        0.0911738786957639,  0.0938443990808046,  0.0956387200792749,  0.0965400885147278,
        0.0965400885147278,  0.0956387200792749,  0.0938443990808046,  0.0911738786957639,
        0.0876520930044038,  0.0833119242269467,  0.0781938957870703,  0.0723457941088485,
        0.0658222227763618,  0.0586840934785355,  0.0509980592623762,  0.0428358980222267,
        0.0342738629130214,  0.0253920653092621,  0.0162743947309057,  0.0070186100094701
    };

    double C = R * R + dz * dz;
    double R_0 = sqrt(C);

    double Omega = 0.0;
    for (int j = 0; j < N_GAUSS; ++j) {
        // Map [-1,1] to [0, 2*pi]: phi = pi + pi*x
        double phi = PI + PI * gp32[j];
        double cos_phi = cos(phi);

        double D = R * R * (1.0 - cos_phi * cos_phi) + dz * dz;  // R^2*sin^2(phi) + z^2
        if (D < 1.0e-30) D = 1.0e-30;

        double R_a = sqrt(CR * CR - 2.0 * CR * R * cos_phi + C);

        // I(phi) = dz/D * [C/R_0 - (C - CR*R*cos(phi))/R_a]
        double I_phi = (dz / D) * (C / R_0 - (C - CR * R * cos_phi) / R_a);

        Omega += gw32[j] * I_phi;
    }

    // Scale by pi (half-interval width for mapping [-1,1] -> [0,2*pi])
    Omega *= PI;

    return Omega;
}

//-------------------------------------------------------------------------
// Incomplete elliptic integrals F(phi, m) and E(phi, m)
// using arithmetic-geometric mean algorithm
//
// Reference:
//   [1] M. Abramowitz, I.A. Stegun, "Handbook of Mathematical Functions",
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
    // Correct formula: E = eok * F + sgn * cs
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

} // namespace RadElliptic
