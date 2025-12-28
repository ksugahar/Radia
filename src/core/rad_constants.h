/*-------------------------------------------------------------------------
*
* File name:      rad_constants.h
*
* Project:        RADIA
*
* Description:    Unified physical and mathematical constants
*
* Author(s):      Claude Code (AI-assisted refactoring)
*
* First release:  2025-12-28
*
* Copyright (C):  European Synchrotron Radiation Facility, France
*
*-------------------------------------------------------------------------
*
* This header provides unified constants for magnetic field computation.
* All magnetic interaction calculations should use these constants to ensure
* consistency between different element types (tetrahedral, hexahedral, etc.)
*
* Mathematical Framework:
* -----------------------
* The magnetostatic problem is solved using the boundary integral equation:
*
*   H_int = H_ext + N * M
*
* where:
*   H_int = internal magnetic field (A/m)
*   H_ext = external/applied magnetic field (A/m)
*   N     = demagnetization tensor (interaction matrix)
*   M     = magnetization (A/m)
*
* For linear materials: M = chi * H_int = chi * (H_ext + N * M)
* Rearranging: (I - chi*N) * M = chi * H_ext
*
* Or equivalently (ELF formulation):
*   (N_demag + 1/chi * I) * M = H_ext
* where N_demag = -N (demagnetization convention)
*
* The interaction matrix elements are computed from solid angle integrals:
*
*   N_ij = (1/4pi) * integral[ (r - r') / |r - r'|^3 ] dS'
*
* The factor 1/(4*pi) appears naturally from the Green's function for the
* Laplacian in 3D: G(r,r') = 1 / (4*pi*|r-r'|)
*
*-------------------------------------------------------------------------*/

#ifndef __RAD_CONSTANTS_H
#define __RAD_CONSTANTS_H

#include <cmath>

namespace RadConst {

//-------------------------------------------------------------------------
// Mathematical Constants
//-------------------------------------------------------------------------

/// Pi (high precision)
constexpr double PI = 3.14159265358979323846264338327950288;

/// 2*Pi
constexpr double TWO_PI = 2.0 * PI;

/// 4*Pi
constexpr double FOUR_PI = 4.0 * PI;

/// 1/(4*Pi) - appears in solid angle integration and Green's function
/// This is THE canonical factor for magnetic interaction matrix elements
constexpr double INV_FOUR_PI = 1.0 / FOUR_PI;

/// 1/(2*Pi)
constexpr double INV_TWO_PI = 1.0 / TWO_PI;

/// 1/Pi
constexpr double INV_PI = 1.0 / PI;

//-------------------------------------------------------------------------
// Physical Constants (SI units)
//-------------------------------------------------------------------------

/// Permeability of free space: mu_0 = 4*pi*10^-7 H/m
constexpr double MU_0 = FOUR_PI * 1.0e-7;

/// 1/mu_0 = 1/(4*pi*10^-7) = 7.9577471545947668e5 m/H
constexpr double INV_MU_0 = 1.0 / MU_0;

/// mu_0 / (4*pi) = 10^-7 H/m
/// Useful for Biot-Savart law: B = (mu_0/4pi) * I * dl x r / r^3
constexpr double MU_0_OVER_FOUR_PI = 1.0e-7;

//-------------------------------------------------------------------------
// Numerical Constants
//-------------------------------------------------------------------------

/// Small epsilon for numerical comparisons and avoiding division by zero
constexpr double EPSILON = 1.0e-14;

/// Tolerance for geometric comparisons (coordinates)
constexpr double GEOM_TOL = 1.0e-12;

/// Tolerance for determining if a point is on a surface
constexpr double SURFACE_TOL = 1.0e-10;

/// Default tolerance for solid angle computation
constexpr double SOLID_ANGLE_TOL = 1.0e-12;

//-------------------------------------------------------------------------
// Field Computation Constants
//-------------------------------------------------------------------------

/// Sign convention for surface charge contribution to H field
/// H = (1/4pi) * integral[ sigma * (r-r')/|r-r'|^3 ] dS
/// For outward normal, sigma > 0 gives field pointing outward
constexpr double SURFACE_CHARGE_H_FACTOR = INV_FOUR_PI;

/// Sign convention for dipole contribution to H field
/// H = (1/4pi) * [ 3*(m.r)*r/r^5 - m/r^3 ]
constexpr double DIPOLE_H_FACTOR = INV_FOUR_PI;

/// Factor for B field from magnetization: B = mu_0 * (H + M)
/// Inside a uniformly magnetized region, B = mu_0*M + mu_0*H_demag
constexpr double B_FROM_M_FACTOR = MU_0;

//-------------------------------------------------------------------------
// Demagnetization Tensor Constants
//-------------------------------------------------------------------------

/// Trace of demagnetization tensor for any closed shape: Tr(N) = -1
/// (Convention: N gives H_demag = N * M, where H_demag opposes M)
constexpr double DEMAG_TENSOR_TRACE = -1.0;

/// Demagnetization factor for infinite slab perpendicular to field: N = -1
constexpr double DEMAG_SLAB_PERP = -1.0;

/// Demagnetization factor for infinite slab parallel to field: N = 0
constexpr double DEMAG_SLAB_PARA = 0.0;

/// Demagnetization factor for sphere: N = -1/3 (each axis)
constexpr double DEMAG_SPHERE = -1.0 / 3.0;

/// Demagnetization factor for infinite cylinder perpendicular to axis: N = -1/2
constexpr double DEMAG_CYLINDER_PERP = -1.0 / 2.0;

/// Demagnetization factor for infinite cylinder along axis: N = 0
constexpr double DEMAG_CYLINDER_PARA = 0.0;

//-------------------------------------------------------------------------
// Inline Helper Functions
//-------------------------------------------------------------------------

/// Check if a value is effectively zero
inline bool IsZero(double val, double tol = EPSILON) {
    return std::abs(val) < tol;
}

/// Check if two values are effectively equal
inline bool AreEqual(double a, double b, double tol = EPSILON) {
    return std::abs(a - b) < tol;
}

/// Safe division avoiding division by zero
inline double SafeDivide(double num, double denom, double default_val = 0.0) {
    return (std::abs(denom) > EPSILON) ? (num / denom) : default_val;
}

/// Clamp a value to a range
inline double Clamp(double val, double min_val, double max_val) {
    return (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
}

} // namespace RadConst

#endif // __RAD_CONSTANTS_H
