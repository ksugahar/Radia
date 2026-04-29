/*-------------------------------------------------------------------------
*
* File name:      rad_biot_savart_filaments.cpp
*
* Project:        RADIA
*
* Description:    Implementation of fast finite-segment Biot-Savart H.
*                 See rad_biot_savart_filaments.h for API.
*
-------------------------------------------------------------------------*/

#include "rad_biot_savart_filaments.h"

#include <cmath>
#include <cstring>
#include <algorithm>

#include <core/taskmanager.hpp>

namespace radia { namespace bs {

namespace {
constexpr double INV_4PI = 1.0 / (4.0 * 3.14159265358979323846);
constexpr double MU_0_OVER_4PI = 1.0e-7;     // mu_0 / (4*pi) [T*m/A]
constexpr double EPS_DEGEN = 1e-30;
constexpr double EPS_LOG = 1e-300;            // safe lower bound for log argument
}

void HFromSegmentsComplex(
    const double* segments, int n_seg,
    const double* obs,      int n_obs,
    const double* I_re,
    const double* I_im,
    double* H_re_out,
    double* H_im_out,
    int n_threads)
{
    // Pre-compute per-segment unit direction, length, and skip flag.
    // Output buffers initialized to zero.
    std::memset(H_re_out, 0, sizeof(double) * 3 * (size_t)n_obs);
    std::memset(H_im_out, 0, sizeof(double) * 3 * (size_t)n_obs);

    struct SegCache {
        double e_l[3];   // unit vector along segment
        double p1[3];
        double p2[3];
        double Ire, Iim;
        bool valid;
    };
    std::vector<SegCache> cache(n_seg);
    for (int k = 0; k < n_seg; ++k) {
        const double* p1 = &segments[k * 6 + 0];
        const double* p2 = &segments[k * 6 + 3];
        double dx = p2[0] - p1[0];
        double dy = p2[1] - p1[1];
        double dz = p2[2] - p1[2];
        double L = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (L < EPS_DEGEN) {
            cache[k].valid = false;
            continue;
        }
        double inv_L = 1.0 / L;
        cache[k].valid = true;
        cache[k].p1[0] = p1[0]; cache[k].p1[1] = p1[1]; cache[k].p1[2] = p1[2];
        cache[k].p2[0] = p2[0]; cache[k].p2[1] = p2[1]; cache[k].p2[2] = p2[2];
        cache[k].e_l[0] = dx * inv_L;
        cache[k].e_l[1] = dy * inv_L;
        cache[k].e_l[2] = dz * inv_L;
        cache[k].Ire = I_re[k];
        cache[k].Iim = I_im[k];
    }

    int requested_threads = (n_threads > 0)
        ? n_threads
        : ngcore::TaskManager::GetMaxThreads();
    if (requested_threads <= 0) requested_threads = 1;
    ngcore::RegionTaskManager rtm(requested_threads);

    // Parallel over obs points.  Each obs is independent.
    ngcore::ParallelFor(ngcore::IntRange(n_obs), [&](size_t j_sz) {
        int j = static_cast<int>(j_sz);
        double ox = obs[j * 3 + 0];
        double oy = obs[j * 3 + 1];
        double oz = obs[j * 3 + 2];
        double H_re[3] = {0.0, 0.0, 0.0};
        double H_im[3] = {0.0, 0.0, 0.0};

        for (int k = 0; k < n_seg; ++k) {
            const SegCache& sc = cache[k];
            if (!sc.valid) continue;
            // r1 = obs - p1
            double r1x = ox - sc.p1[0];
            double r1y = oy - sc.p1[1];
            double r1z = oz - sc.p1[2];
            // r2 = obs - p2
            double r2x = ox - sc.p2[0];
            double r2y = oy - sc.p2[1];
            double r2z = oz - sc.p2[2];
            // cross = e_l x r1
            double cx = sc.e_l[1] * r1z - sc.e_l[2] * r1y;
            double cy = sc.e_l[2] * r1x - sc.e_l[0] * r1z;
            double cz = sc.e_l[0] * r1y - sc.e_l[1] * r1x;
            double d2 = cx*cx + cy*cy + cz*cz;
            if (d2 < EPS_DEGEN) continue;
            double d = std::sqrt(d2);
            double r1_mag = std::sqrt(r1x*r1x + r1y*r1y + r1z*r1z);
            double r2_mag = std::sqrt(r2x*r2x + r2y*r2y + r2z*r2z);
            if (r1_mag < EPS_DEGEN || r2_mag < EPS_DEGEN) continue;
            double cos_a1 = (r1x*sc.e_l[0] + r1y*sc.e_l[1] + r1z*sc.e_l[2]) / r1_mag;
            double cos_a2 = (r2x*sc.e_l[0] + r2y*sc.e_l[1] + r2z*sc.e_l[2]) / r2_mag;
            double geom  = INV_4PI / d * (cos_a1 - cos_a2);
            // e_perp = cross / d
            double inv_d = 1.0 / d;
            double ex = cx * inv_d;
            double ey = cy * inv_d;
            double ez = cz * inv_d;
            // H_k = I * geom * e_perp.  Split complex.
            double w_re = sc.Ire * geom;
            double w_im = sc.Iim * geom;
            H_re[0] += w_re * ex;
            H_re[1] += w_re * ey;
            H_re[2] += w_re * ez;
            H_im[0] += w_im * ex;
            H_im[1] += w_im * ey;
            H_im[2] += w_im * ez;
        }

        H_re_out[j*3 + 0] = H_re[0];
        H_re_out[j*3 + 1] = H_re[1];
        H_re_out[j*3 + 2] = H_re[2];
        H_im_out[j*3 + 0] = H_im[0];
        H_im_out[j*3 + 1] = H_im[1];
        H_im_out[j*3 + 2] = H_im[2];
    });
}

void AFromSegmentsComplex(
    const double* segments, int n_seg,
    const double* obs,      int n_obs,
    const double* I_re,
    const double* I_im,
    double* A_re_out,
    double* A_im_out,
    int n_threads)
{
    std::memset(A_re_out, 0, sizeof(double) * 3 * (size_t)n_obs);
    std::memset(A_im_out, 0, sizeof(double) * 3 * (size_t)n_obs);

    struct SegCacheA {
        double e_l[3];   // unit vector along segment
        double p1[3];
        double p2[3];
        double L;
        double Ire, Iim;
        bool valid;
    };
    std::vector<SegCacheA> cache(n_seg);
    for (int k = 0; k < n_seg; ++k) {
        const double* p1 = &segments[k * 6 + 0];
        const double* p2 = &segments[k * 6 + 3];
        double dx = p2[0] - p1[0];
        double dy = p2[1] - p1[1];
        double dz = p2[2] - p1[2];
        double L = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (L < EPS_DEGEN) {
            cache[k].valid = false;
            continue;
        }
        double inv_L = 1.0 / L;
        cache[k].valid = true;
        cache[k].p1[0] = p1[0]; cache[k].p1[1] = p1[1]; cache[k].p1[2] = p1[2];
        cache[k].p2[0] = p2[0]; cache[k].p2[1] = p2[1]; cache[k].p2[2] = p2[2];
        cache[k].e_l[0] = dx * inv_L;
        cache[k].e_l[1] = dy * inv_L;
        cache[k].e_l[2] = dz * inv_L;
        cache[k].L = L;
        cache[k].Ire = I_re[k];
        cache[k].Iim = I_im[k];
    }

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

        for (int k = 0; k < n_seg; ++k) {
            const SegCacheA& sc = cache[k];
            if (!sc.valid) continue;
            // r1 = obs - p1
            double r1x = ox - sc.p1[0];
            double r1y = oy - sc.p1[1];
            double r1z = oz - sc.p1[2];
            // r2 = obs - p2
            double r2x = ox - sc.p2[0];
            double r2y = oy - sc.p2[1];
            double r2z = oz - sc.p2[2];
            double r1_mag = std::sqrt(r1x*r1x + r1y*r1y + r1z*r1z);
            double r2_mag = std::sqrt(r2x*r2x + r2y*r2y + r2z*r2z);
            // t = (obs - p1) . e_l  (signed projection on wire axis)
            double t  = r1x*sc.e_l[0] + r1y*sc.e_l[1] + r1z*sc.e_l[2];
            // numer = (L - t) + r2_mag,  denom = -t + r1_mag
            // Both are non-negative (numer = r2_mag + along-wire dist from p2 to obs's projection;
            // denom = r1_mag - signed-dist past p1) and become 0 only when obs is on the
            // wire's extension past p2 (numer) or past p1 (denom).
            double numer = (sc.L - t) + r2_mag;
            double denom = -t + r1_mag;
            if (numer < EPS_LOG) numer = EPS_LOG;
            if (denom < EPS_LOG) denom = EPS_LOG;
            double log_factor = std::log(numer / denom);
            // A_k = (mu_0/(4*pi)) * I_k * log_factor * e_l
            double w_re = MU_0_OVER_4PI * sc.Ire * log_factor;
            double w_im = MU_0_OVER_4PI * sc.Iim * log_factor;
            A_re[0] += w_re * sc.e_l[0];
            A_re[1] += w_re * sc.e_l[1];
            A_re[2] += w_re * sc.e_l[2];
            A_im[0] += w_im * sc.e_l[0];
            A_im[1] += w_im * sc.e_l[1];
            A_im[2] += w_im * sc.e_l[2];
        }

        A_re_out[j*3 + 0] = A_re[0];
        A_re_out[j*3 + 1] = A_re[1];
        A_re_out[j*3 + 2] = A_re[2];
        A_im_out[j*3 + 0] = A_im[0];
        A_im_out[j*3 + 1] = A_im[1];
        A_im_out[j*3 + 2] = A_im[2];
    });
}

}}  // namespace radia::bs
