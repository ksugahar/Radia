#ifndef RAD_FMM3D_H
#define RAD_FMM3D_H

/**
 * @file rad_fmm3d.h
 * @brief FMM3D wrapper for fast multipole method acceleration
 *
 * Provides O(N) field computation using FMM3D library instead of O(N^2) direct.
 * Uses dipole approximation: each element represented as a point dipole at its center.
 *
 * Algorithm:
 *   H_total = H_far(FMM) + H_correction(near) + H_external
 *
 * where:
 *   H_far: FMM3D dipole field from all elements (O(N log N))
 *   H_correction: Exact MSC - dipole approximation for nearby elements
 *   H_external: Background field
 */

#include "gmvect.h"
#include <vector>

namespace RadFMM3D {

/**
 * @brief Check if FMM3D library is available
 *
 * @return true if FMM3D is linked and available
 */
bool IsAvailable();

/**
 * @brief Compute H field from dipole sources using FMM3D
 *
 * Uses lfmm3d_t_d_g_: dipole sources to target gradient.
 * The gradient of scalar potential gives H field:
 *   H = -grad(phi_m) = -grad(sum_i (m_i . r_i) / (4*pi*|r_i|^3))
 *
 * FMM3D computes gradient of potential as (grad_x, grad_y, grad_z),
 * which equals -H field directly.
 *
 * @param n_src Number of dipole sources (elements)
 * @param src_pos Source positions [n_src * 3] (x1,y1,z1,x2,y2,z2,...)
 * @param dipole_mom Dipole moments [n_src * 3] (mx1,my1,mz1,mx2,my2,mz2,...)
 * @param n_targ Number of target points
 * @param targ_pos Target positions [n_targ * 3]
 * @param H_out Output H field [n_targ * 3] (Hx1,Hy1,Hz1,...)
 * @param eps FMM tolerance (default 1e-6)
 * @param ier Output error code (0 = success)
 */
void ComputeDipoleField(
    int n_src, const double* src_pos, const double* dipole_mom,
    int n_targ, const double* targ_pos, double* H_out,
    double eps, int* ier);

/**
 * @brief Compute scalar potential phi_m from dipole sources using FMM3D
 *
 * Uses lfmm3d_t_d_p_: dipole sources to target potential.
 * phi_m = sum_i (m_i . r_i) / (4*pi*|r_i|^3)
 *
 * Note: The sign convention means H = -grad(phi_m)
 *
 * @param n_src Number of dipole sources
 * @param src_pos Source positions [n_src * 3]
 * @param dipole_mom Dipole moments [n_src * 3]
 * @param n_targ Number of target points
 * @param targ_pos Target positions [n_targ * 3]
 * @param phi_out Output scalar potential [n_targ]
 * @param eps FMM tolerance (default 1e-6)
 * @param ier Output error code (0 = success)
 */
void ComputeScalarPotential(
    int n_src, const double* src_pos, const double* dipole_mom,
    int n_targ, const double* targ_pos, double* phi_out,
    double eps, int* ier);

/**
 * @brief Direct O(N^2) dipole field computation (for verification)
 *
 * Computes dipole field directly without FMM for small problems or verification.
 * H = (1/4*pi) * sum_i [ 3(m_i . r_i)r_i - m_i*|r_i|^2 ] / |r_i|^5
 *
 * @param n_src Number of dipole sources
 * @param src_pos Source positions [n_src * 3]
 * @param dipole_mom Dipole moments [n_src * 3]
 * @param n_targ Number of target points
 * @param targ_pos Target positions [n_targ * 3]
 * @param H_out Output H field [n_targ * 3]
 */
void ComputeDipoleFieldDirect(
    int n_src, const double* src_pos, const double* dipole_mom,
    int n_targ, const double* targ_pos, double* H_out);

/**
 * @brief Compute dipole field at single target point (inline)
 *
 * For near-field correction calculations.
 *
 * @param m Dipole moment vector (mx, my, mz)
 * @param r Vector from source to target (rx, ry, rz)
 * @return H field vector at target
 */
TVector3d DipoleFieldSingle(const TVector3d& m, const TVector3d& r);

/**
 * @brief FMM parameters structure
 */
struct FMMParams {
    double eps;          ///< FMM tolerance (default 1e-6)
    double near_thresh;  ///< Near field threshold (multiplier of element size, default 3.0)

    FMMParams() : eps(1e-6), near_thresh(3.0) {}
};

/**
 * @brief Global FMM parameters
 */
extern FMMParams gFMMParams;

/**
 * @brief Set FMM tolerance
 */
void SetEpsilon(double eps);

/**
 * @brief Set near field threshold
 */
void SetNearThreshold(double thresh);

/**
 * @brief Get current FMM tolerance
 */
double GetEpsilon();

/**
 * @brief Get current near field threshold
 */
double GetNearThreshold();

} // namespace RadFMM3D

#endif // RAD_FMM3D_H
