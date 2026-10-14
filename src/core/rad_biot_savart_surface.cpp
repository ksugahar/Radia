/*-------------------------------------------------------------------------
*
* File name:      rad_biot_savart_surface.cpp
*
* Project:        RADIA
*
* Description:    Implementation of Biot-Savart B/A from a triangulated
*                 surface with piecewise-constant complex J_s.
*                 See rad_biot_savart_surface.h for API.
*
-------------------------------------------------------------------------*/

#include "rad_biot_savart_surface.h"

#include <cmath>
#include <cstring>
#include <vector>

#include <core/taskmanager.hpp>

namespace radia { namespace bs {

namespace {
constexpr double MU_0_OVER_4PI = 1.0e-7;
constexpr double EPS_DEGEN     = 1e-30;

// 3-point Gauss quadrature on triangle (mid-edge), barycentric coords +
// weights summing to 1.  Exact for linear, ~1% on 1/r at >= 1 edge away.
constexpr double Q_BARY[3][3] = {
    {0.5, 0.5, 0.0},
    {0.0, 0.5, 0.5},
    {0.5, 0.0, 0.5}
};
constexpr double Q_W = 1.0 / 3.0;

struct TriCache {
    double v[3][3];      // 3 vertices x 3 coords
    double qpt[3][3];    // 3 quadrature points x 3 coords
    double area;
    double Jre[3];
    double Jim[3];
    bool   valid;
};

void build_tri_cache(const double* verts, int n_t,
                     const double* J_re, const double* J_im,
                     std::vector<TriCache>& cache)
{
    cache.resize(n_t);
    for (int t = 0; t < n_t; ++t) {
        TriCache& tc = cache[t];
        const double* p = &verts[t * 9];
        for (int i = 0; i < 3; ++i)
            for (int d = 0; d < 3; ++d)
                tc.v[i][d] = p[i * 3 + d];

        // Triangle area via |((v2-v1) x (v3-v1))| / 2
        double e12[3] = { tc.v[1][0]-tc.v[0][0],
                          tc.v[1][1]-tc.v[0][1],
                          tc.v[1][2]-tc.v[0][2] };
        double e13[3] = { tc.v[2][0]-tc.v[0][0],
                          tc.v[2][1]-tc.v[0][1],
                          tc.v[2][2]-tc.v[0][2] };
        double cx = e12[1]*e13[2] - e12[2]*e13[1];
        double cy = e12[2]*e13[0] - e12[0]*e13[2];
        double cz = e12[0]*e13[1] - e12[1]*e13[0];
        double cross_mag = std::sqrt(cx*cx + cy*cy + cz*cz);
        tc.area = 0.5 * cross_mag;
        if (tc.area < EPS_DEGEN) { tc.valid = false; continue; }
        tc.valid = true;

        // Pre-compute quadrature points in world coords.
        for (int q = 0; q < 3; ++q) {
            for (int d = 0; d < 3; ++d) {
                tc.qpt[q][d] = Q_BARY[q][0] * tc.v[0][d]
                             + Q_BARY[q][1] * tc.v[1][d]
                             + Q_BARY[q][2] * tc.v[2][d];
            }
        }

        for (int d = 0; d < 3; ++d) {
            tc.Jre[d] = J_re[t * 3 + d];
            tc.Jim[d] = J_im[t * 3 + d];
        }
    }
}

}  // anonymous namespace

void AFromTrianglesComplex(
    const double* verts, int n_t,
    const double* J_re,
    const double* J_im,
    const double* obs, int n_obs,
    double* A_re_out,
    double* A_im_out,
    int n_threads)
{
    std::memset(A_re_out, 0, sizeof(double) * 3 * (size_t)n_obs);
    std::memset(A_im_out, 0, sizeof(double) * 3 * (size_t)n_obs);

    std::vector<TriCache> cache;
    build_tri_cache(verts, n_t, J_re, J_im, cache);

    int requested_threads = (n_threads > 0)
        ? n_threads
        : ngcore::TaskManager::GetMaxThreads();
    if (requested_threads <= 0) requested_threads = 1;
    ngcore::RegionTaskManager rtm(requested_threads);

    ngcore::ParallelFor(ngcore::IntRange(n_obs), [&](size_t j_sz) {
        int j = static_cast<int>(j_sz);
        double ox = obs[j * 3 + 0];
        double oy = obs[j * 3 + 1];
        double oz = obs[j * 3 + 2];
        double A_re[3] = {0.0, 0.0, 0.0};
        double A_im[3] = {0.0, 0.0, 0.0};

        for (int t = 0; t < n_t; ++t) {
            const TriCache& tc = cache[t];
            if (!tc.valid) continue;

            // Sum_q w_q * 1/|P_obs - P_q|
            double inv_r_sum = 0.0;
            for (int q = 0; q < 3; ++q) {
                double dx = ox - tc.qpt[q][0];
                double dy = oy - tc.qpt[q][1];
                double dz = oz - tc.qpt[q][2];
                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < EPS_DEGEN) continue;   // singular: skip
                inv_r_sum += Q_W / std::sqrt(r2);
            }
            // Triangle contribution: pref * J_s_T * area * sum_q
            double pref = MU_0_OVER_4PI * tc.area * inv_r_sum;
            A_re[0] += pref * tc.Jre[0];
            A_re[1] += pref * tc.Jre[1];
            A_re[2] += pref * tc.Jre[2];
            A_im[0] += pref * tc.Jim[0];
            A_im[1] += pref * tc.Jim[1];
            A_im[2] += pref * tc.Jim[2];
        }

        A_re_out[j*3 + 0] = A_re[0];
        A_re_out[j*3 + 1] = A_re[1];
        A_re_out[j*3 + 2] = A_re[2];
        A_im_out[j*3 + 0] = A_im[0];
        A_im_out[j*3 + 1] = A_im[1];
        A_im_out[j*3 + 2] = A_im[2];
    });
}


void BFromTrianglesComplex(
    const double* verts, int n_t,
    const double* J_re,
    const double* J_im,
    const double* obs, int n_obs,
    double* B_re_out,
    double* B_im_out,
    int n_threads)
{
    std::memset(B_re_out, 0, sizeof(double) * 3 * (size_t)n_obs);
    std::memset(B_im_out, 0, sizeof(double) * 3 * (size_t)n_obs);

    std::vector<TriCache> cache;
    build_tri_cache(verts, n_t, J_re, J_im, cache);

    int requested_threads = (n_threads > 0)
        ? n_threads
        : ngcore::TaskManager::GetMaxThreads();
    if (requested_threads <= 0) requested_threads = 1;
    ngcore::RegionTaskManager rtm(requested_threads);

    ngcore::ParallelFor(ngcore::IntRange(n_obs), [&](size_t j_sz) {
        int j = static_cast<int>(j_sz);
        double ox = obs[j * 3 + 0];
        double oy = obs[j * 3 + 1];
        double oz = obs[j * 3 + 2];
        double B_re[3] = {0.0, 0.0, 0.0};
        double B_im[3] = {0.0, 0.0, 0.0};

        for (int t = 0; t < n_t; ++t) {
            const TriCache& tc = cache[t];
            if (!tc.valid) continue;

            // sum_q w_q * (P_obs - P_q) / |P_obs - P_q|^3
            double R[3] = {0.0, 0.0, 0.0};
            for (int q = 0; q < 3; ++q) {
                double dx = ox - tc.qpt[q][0];
                double dy = oy - tc.qpt[q][1];
                double dz = oz - tc.qpt[q][2];
                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < EPS_DEGEN) continue;
                double inv_r3 = Q_W / (r2 * std::sqrt(r2));
                R[0] += dx * inv_r3;
                R[1] += dy * inv_r3;
                R[2] += dz * inv_r3;
            }
            double pref = MU_0_OVER_4PI * tc.area;
            // J_s_T x R component-wise; complex via Re/Im of J only (R is real)
            // (J x R)_x = Jy*Rz - Jz*Ry
            // (J x R)_y = Jz*Rx - Jx*Rz
            // (J x R)_z = Jx*Ry - Jy*Rx
            B_re[0] += pref * (tc.Jre[1]*R[2] - tc.Jre[2]*R[1]);
            B_re[1] += pref * (tc.Jre[2]*R[0] - tc.Jre[0]*R[2]);
            B_re[2] += pref * (tc.Jre[0]*R[1] - tc.Jre[1]*R[0]);
            B_im[0] += pref * (tc.Jim[1]*R[2] - tc.Jim[2]*R[1]);
            B_im[1] += pref * (tc.Jim[2]*R[0] - tc.Jim[0]*R[2]);
            B_im[2] += pref * (tc.Jim[0]*R[1] - tc.Jim[1]*R[0]);
        }

        B_re_out[j*3 + 0] = B_re[0];
        B_re_out[j*3 + 1] = B_re[1];
        B_re_out[j*3 + 2] = B_re[2];
        B_im_out[j*3 + 0] = B_im[0];
        B_im_out[j*3 + 1] = B_im[1];
        B_im_out[j*3 + 2] = B_im[2];
    });
}

}}  // namespace radia::bs
