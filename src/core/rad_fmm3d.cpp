/**
 * @file rad_fmm3d.cpp
 * @brief FMM3D wrapper implementation for fast multipole method acceleration
 */

#include "rad_fmm3d.h"
#include <cmath>
#include <cstring>
#include <vector>

// FMM3D library interface (Fortran, via C linkage)
// We declare the functions we need directly to avoid complex.h issues
#ifdef RADIA_FMM3D_ENABLED
extern "C" {
    // Laplace FMM: dipole sources to target gradient (what we need for H field)
    // t = target only, d = dipole sources, g = compute gradient
    void lfmm3d_t_d_g_(double *eps, int64_t *nsource,
                       double *source, double *dipvec, int64_t *nt, double *targ,
                       double *pottarg, double *gradtarg, int64_t *ier);

    // Laplace FMM: dipole sources to target potential (for scalar potential phi)
    void lfmm3d_t_d_p_(double *eps, int64_t *nsource,
                       double *source, double *dipvec, int64_t *nt, double *targ,
                       double *pottarg, int64_t *ier);
}
#endif

namespace RadFMM3D {

// Global FMM parameters
FMMParams gFMMParams;

//-----------------------------------------------------------------------------
// Check if FMM3D is available
//-----------------------------------------------------------------------------
bool IsAvailable()
{
#ifdef RADIA_FMM3D_ENABLED
    return true;
#else
    return false;
#endif
}

//-----------------------------------------------------------------------------
// Set FMM tolerance
//-----------------------------------------------------------------------------
void SetEpsilon(double eps)
{
    if (eps > 0.0 && eps < 1.0) {
        gFMMParams.eps = eps;
    }
}

//-----------------------------------------------------------------------------
// Set near field threshold
//-----------------------------------------------------------------------------
void SetNearThreshold(double thresh)
{
    if (thresh > 0.0) {
        gFMMParams.near_thresh = thresh;
    }
}

//-----------------------------------------------------------------------------
// Get current FMM tolerance
//-----------------------------------------------------------------------------
double GetEpsilon()
{
    return gFMMParams.eps;
}

//-----------------------------------------------------------------------------
// Get current near field threshold
//-----------------------------------------------------------------------------
double GetNearThreshold()
{
    return gFMMParams.near_thresh;
}

//-----------------------------------------------------------------------------
// Compute dipole field at single target point
//-----------------------------------------------------------------------------
TVector3d DipoleFieldSingle(const TVector3d& m, const TVector3d& r)
{
    // H = (1/4*pi) * [ 3(m . r)r / r^5 - m / r^3 ]
    //   = (1/4*pi) * [ 3(m . r)r - m*r^2 ] / r^5

    double r2 = r.x * r.x + r.y * r.y + r.z * r.z;
    if (r2 < 1e-30) {
        // Avoid division by zero for coincident points
        return TVector3d(0.0, 0.0, 0.0);
    }

    double r_mag = std::sqrt(r2);
    double r3 = r2 * r_mag;
    double r5 = r3 * r2;

    double m_dot_r = m.x * r.x + m.y * r.y + m.z * r.z;
    double coef3 = 3.0 * m_dot_r / r5;

    // Factor of 1/(4*pi) included
    const double inv_4pi = 0.07957747154594767;  // 1/(4*pi)

    TVector3d H;
    H.x = inv_4pi * (coef3 * r.x - m.x / r3);
    H.y = inv_4pi * (coef3 * r.y - m.y / r3);
    H.z = inv_4pi * (coef3 * r.z - m.z / r3);

    return H;
}

//-----------------------------------------------------------------------------
// Direct O(N^2) dipole field computation
//-----------------------------------------------------------------------------
void ComputeDipoleFieldDirect(
    int n_src, const double* src_pos, const double* dipole_mom,
    int n_targ, const double* targ_pos, double* H_out)
{
    // Initialize output to zero
    std::memset(H_out, 0, n_targ * 3 * sizeof(double));

    // Loop over all targets
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int t = 0; t < n_targ; ++t) {
        TVector3d pt;
        pt.x = targ_pos[t * 3 + 0];
        pt.y = targ_pos[t * 3 + 1];
        pt.z = targ_pos[t * 3 + 2];

        TVector3d H_sum(0.0, 0.0, 0.0);

        // Sum contributions from all sources
        for (int s = 0; s < n_src; ++s) {
            TVector3d m;
            m.x = dipole_mom[s * 3 + 0];
            m.y = dipole_mom[s * 3 + 1];
            m.z = dipole_mom[s * 3 + 2];

            TVector3d src;
            src.x = src_pos[s * 3 + 0];
            src.y = src_pos[s * 3 + 1];
            src.z = src_pos[s * 3 + 2];

            // Vector from source to target
            TVector3d r;
            r.x = pt.x - src.x;
            r.y = pt.y - src.y;
            r.z = pt.z - src.z;

            TVector3d H_contrib = DipoleFieldSingle(m, r);
            H_sum.x += H_contrib.x;
            H_sum.y += H_contrib.y;
            H_sum.z += H_contrib.z;
        }

        H_out[t * 3 + 0] = H_sum.x;
        H_out[t * 3 + 1] = H_sum.y;
        H_out[t * 3 + 2] = H_sum.z;
    }
}

//-----------------------------------------------------------------------------
// Compute H field from dipole sources using FMM3D
//-----------------------------------------------------------------------------
void ComputeDipoleField(
    int n_src, const double* src_pos, const double* dipole_mom,
    int n_targ, const double* targ_pos, double* H_out,
    double eps, int* ier)
{
#ifdef RADIA_FMM3D_ENABLED
    // FMM3D uses Fortran 64-bit integers
    int64_t ns = static_cast<int64_t>(n_src);
    int64_t nt = static_cast<int64_t>(n_targ);
    int64_t err = 0;

    // Allocate arrays for FMM3D (column-major order for Fortran)
    // FMM3D expects sources as [3, nsource] and dipoles as [3, nsource]
    // We need to transpose from row-major to column-major

    std::vector<double> src_fortran(n_src * 3);
    std::vector<double> dip_fortran(n_src * 3);
    std::vector<double> targ_fortran(n_targ * 3);
    std::vector<double> pot_out(n_targ);      // Scalar potential (not used but required)
    std::vector<double> grad_out(n_targ * 3); // Gradient = -H field

    // Transpose source positions: [n,3] -> [3,n] (column-major)
    for (int i = 0; i < n_src; ++i) {
        src_fortran[i]             = src_pos[i * 3 + 0];  // x
        src_fortran[n_src + i]     = src_pos[i * 3 + 1];  // y
        src_fortran[2 * n_src + i] = src_pos[i * 3 + 2];  // z
    }

    // Transpose dipole moments: [n,3] -> [3,n]
    for (int i = 0; i < n_src; ++i) {
        dip_fortran[i]             = dipole_mom[i * 3 + 0];  // mx
        dip_fortran[n_src + i]     = dipole_mom[i * 3 + 1];  // my
        dip_fortran[2 * n_src + i] = dipole_mom[i * 3 + 2];  // mz
    }

    // Transpose target positions: [n,3] -> [3,n]
    for (int i = 0; i < n_targ; ++i) {
        targ_fortran[i]              = targ_pos[i * 3 + 0];  // x
        targ_fortran[n_targ + i]     = targ_pos[i * 3 + 1];  // y
        targ_fortran[2 * n_targ + i] = targ_pos[i * 3 + 2];  // z
    }

    // Call FMM3D: dipole sources to target gradient
    // lfmm3d_t_d_g_: targets only, dipole sources, compute gradient
    lfmm3d_t_d_g_(
        &eps,
        &ns,
        src_fortran.data(),
        dip_fortran.data(),
        &nt,
        targ_fortran.data(),
        pot_out.data(),
        grad_out.data(),
        &err
    );

    *ier = static_cast<int>(err);

    // FMM3D returns gradient of scalar potential
    // For dipoles: phi(r) = (1/4pi) * sum_i (m_i . (r - r_i)) / |r - r_i|^3
    // grad(phi) = -H (by convention)
    // So H = -grad_out

    // Transpose gradient back: [3,n] -> [n,3] and negate for H
    for (int i = 0; i < n_targ; ++i) {
        H_out[i * 3 + 0] = -grad_out[i];               // Hx = -grad_x
        H_out[i * 3 + 1] = -grad_out[n_targ + i];      // Hy = -grad_y
        H_out[i * 3 + 2] = -grad_out[2 * n_targ + i];  // Hz = -grad_z
    }

#else
    // FMM3D not available, fall back to direct computation
    ComputeDipoleFieldDirect(n_src, src_pos, dipole_mom, n_targ, targ_pos, H_out);
    *ier = 0;
#endif
}

//-----------------------------------------------------------------------------
// Compute scalar potential from dipole sources using FMM3D
//-----------------------------------------------------------------------------
void ComputeScalarPotential(
    int n_src, const double* src_pos, const double* dipole_mom,
    int n_targ, const double* targ_pos, double* phi_out,
    double eps, int* ier)
{
#ifdef RADIA_FMM3D_ENABLED
    int64_t ns = static_cast<int64_t>(n_src);
    int64_t nt = static_cast<int64_t>(n_targ);
    int64_t err = 0;

    // Transpose arrays for Fortran column-major order
    std::vector<double> src_fortran(n_src * 3);
    std::vector<double> dip_fortran(n_src * 3);
    std::vector<double> targ_fortran(n_targ * 3);

    for (int i = 0; i < n_src; ++i) {
        src_fortran[i]             = src_pos[i * 3 + 0];
        src_fortran[n_src + i]     = src_pos[i * 3 + 1];
        src_fortran[2 * n_src + i] = src_pos[i * 3 + 2];
    }

    for (int i = 0; i < n_src; ++i) {
        dip_fortran[i]             = dipole_mom[i * 3 + 0];
        dip_fortran[n_src + i]     = dipole_mom[i * 3 + 1];
        dip_fortran[2 * n_src + i] = dipole_mom[i * 3 + 2];
    }

    for (int i = 0; i < n_targ; ++i) {
        targ_fortran[i]              = targ_pos[i * 3 + 0];
        targ_fortran[n_targ + i]     = targ_pos[i * 3 + 1];
        targ_fortran[2 * n_targ + i] = targ_pos[i * 3 + 2];
    }

    // Call FMM3D: dipole sources to target potential
    lfmm3d_t_d_p_(
        &eps,
        &ns,
        src_fortran.data(),
        dip_fortran.data(),
        &nt,
        targ_fortran.data(),
        phi_out,
        &err
    );

    *ier = static_cast<int>(err);

#else
    // FMM3D not available, compute directly
    // phi = (1/4pi) * sum_i (m_i . r_i) / |r_i|^3
    const double inv_4pi = 0.07957747154594767;

    std::memset(phi_out, 0, n_targ * sizeof(double));

    for (int t = 0; t < n_targ; ++t) {
        double px = targ_pos[t * 3 + 0];
        double py = targ_pos[t * 3 + 1];
        double pz = targ_pos[t * 3 + 2];

        double phi = 0.0;

        for (int s = 0; s < n_src; ++s) {
            double mx = dipole_mom[s * 3 + 0];
            double my = dipole_mom[s * 3 + 1];
            double mz = dipole_mom[s * 3 + 2];

            double sx = src_pos[s * 3 + 0];
            double sy = src_pos[s * 3 + 1];
            double sz = src_pos[s * 3 + 2];

            double rx = px - sx;
            double ry = py - sy;
            double rz = pz - sz;

            double r2 = rx * rx + ry * ry + rz * rz;
            if (r2 < 1e-30) continue;

            double r = std::sqrt(r2);
            double r3 = r2 * r;

            double m_dot_r = mx * rx + my * ry + mz * rz;
            phi += m_dot_r / r3;
        }

        phi_out[t] = inv_4pi * phi;
    }

    *ier = 0;
#endif
}

} // namespace RadFMM3D
