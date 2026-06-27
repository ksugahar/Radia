/*-------------------------------------------------------------------------
*
* File name:      rad_equivalence_source.cpp
*
* Project:        RADIA
*
* Description:    Implementation of Schelkunoff/Love static H
*                 reconstruction.  See rad_equivalence_source.h for API.
*
*                 Replaces the Python Stratton-Chu loop in
*                 src/radia/equivalence_source.py:evaluate_static_H().
*                 Expected speedup ~50-100x vs numpy at N=1e3-1e4 faces
*                 (CLAUDE.md TaskManager + tight per-face inner loop).
*
*                 Sign convention verified against analytic magnetic
*                 dipole on a sphere mesh; see phase1_static_coil.py
*                 (Phase 1 golden, 0.83% PASS).
*
-------------------------------------------------------------------------*/

#include "rad_equivalence_source.h"

#include <cmath>
#include <complex>
#include <cstring>
#include <vector>

#include <core/taskmanager.hpp>

namespace radia { namespace eqsrc {

namespace {

constexpr double PI_         = 3.14159265358979323846;
constexpr double INV_FOUR_PI = 1.0 / (4.0 * PI_);
constexpr double R_MIN       = 1.0e-15;     // singular-distance guard
constexpr double MU_0        = 4.0e-7 * PI_;
constexpr double EPS_0       = 8.8541878128e-12;

// Per-face pre-computed values: centroid + outward normal + area
// + the two surface-source terms (J_s and n.H).  Hot-loop only needs
// these, no recomputation per obs point.
struct FaceCache {
    double cx, cy, cz;          // centroid
    double nx, ny, nz;          // outward unit normal
    double area;
    double Jsx, Jsy, Jsz;       // J_s = n x H_surf  [A/m]
    double n_dot_H;             // n . H_surf  (= rho_m / mu_0)
};

void build_face_cache(
    const double* centroids, const double* normals,
    const double* areas, const double* H_surf,
    int n_faces, std::vector<FaceCache>& cache)
{
    cache.resize(n_faces);
    for (int t = 0; t < n_faces; ++t) {
        FaceCache& fc = cache[t];
        fc.cx = centroids[t * 3 + 0];
        fc.cy = centroids[t * 3 + 1];
        fc.cz = centroids[t * 3 + 2];
        fc.nx = normals[t * 3 + 0];
        fc.ny = normals[t * 3 + 1];
        fc.nz = normals[t * 3 + 2];
        fc.area = areas[t];
        const double Hx = H_surf[t * 3 + 0];
        const double Hy = H_surf[t * 3 + 1];
        const double Hz = H_surf[t * 3 + 2];
        // J_s = n x H_surf
        fc.Jsx = fc.ny * Hz - fc.nz * Hy;
        fc.Jsy = fc.nz * Hx - fc.nx * Hz;
        fc.Jsz = fc.nx * Hy - fc.ny * Hx;
        // n . H_surf  (used as rho_m / mu_0 coefficient)
        fc.n_dot_H = fc.nx * Hx + fc.ny * Hy + fc.nz * Hz;
    }
}

}  // anonymous namespace


void EvaluateStaticH(
    const double* centroids,
    const double* normals,
    const double* areas,
    const double* H_surf,
    int n_faces,
    const double* obs,
    int n_obs,
    double* H_out,
    int n_threads)
{
    std::memset(H_out, 0, sizeof(double) * 3 * (size_t)n_obs);

    if (n_faces <= 0 || n_obs <= 0) return;

    std::vector<FaceCache> cache;
    build_face_cache(centroids, normals, areas, H_surf, n_faces, cache);

    int requested_threads = (n_threads > 0)
        ? n_threads
        : ngcore::TaskManager::GetMaxThreads();
    if (requested_threads <= 0) requested_threads = 1;
    ngcore::RegionTaskManager rtm(requested_threads);

    ngcore::ParallelFor(ngcore::IntRange(n_obs), [&](size_t j_sz) {
        const int j = static_cast<int>(j_sz);
        const double ox = obs[j * 3 + 0];
        const double oy = obs[j * 3 + 1];
        const double oz = obs[j * 3 + 2];
        double Hx_acc = 0.0;
        double Hy_acc = 0.0;
        double Hz_acc = 0.0;

        for (int t = 0; t < n_faces; ++t) {
            const FaceCache& fc = cache[t];
            // R_vec = obs - centroid
            const double Rx = ox - fc.cx;
            const double Ry = oy - fc.cy;
            const double Rz = oz - fc.cz;
            const double R2 = Rx * Rx + Ry * Ry + Rz * Rz;
            if (R2 < R_MIN * R_MIN) continue;   // skip singular face (obs on it)
            const double R = std::sqrt(R2);
            const double R3 = R2 * R;
            // grad(1/R) = -R_vec / R^3  (w.r.t. obs)
            const double gx = -Rx / R3;
            const double gy = -Ry / R3;
            const double gz = -Rz / R3;
            // Term 1: grad(1/R) x J_s  (vector)
            const double t1x = gy * fc.Jsz - gz * fc.Jsy;
            const double t1y = gz * fc.Jsx - gx * fc.Jsz;
            const double t1z = gx * fc.Jsy - gy * fc.Jsx;
            // Term 2: -(n . H_surf) grad(1/R)  (vector)
            const double t2x = -fc.n_dot_H * gx;
            const double t2y = -fc.n_dot_H * gy;
            const double t2z = -fc.n_dot_H * gz;
            // Per-face contribution: (t1 + t2) * dS
            Hx_acc += (t1x + t2x) * fc.area;
            Hy_acc += (t1y + t2y) * fc.area;
            Hz_acc += (t1z + t2z) * fc.area;
        }
        // Multiply by 1/(4 pi) once at the end
        H_out[j * 3 + 0] = INV_FOUR_PI * Hx_acc;
        H_out[j * 3 + 1] = INV_FOUR_PI * Hy_acc;
        H_out[j * 3 + 2] = INV_FOUR_PI * Hz_acc;
    });
}


// =====================================================================
// EvaluateHarmonic -- full Stratton-Chu with DYADIC Green's function
// =====================================================================
//
//   G_bar . v = psi {(h^2+h+1)/h^2 v - (h^2+3h+3)/h^2 (v.R_hat) R_hat}
//
//   where h = jkR (complex), psi = exp(-jkR)/(4 pi R).
//
//   E(r) = integral {-jw mu_0 G_bar.J_s + grad_psi x M_s}dS
//   H(r) = integral {+jw eps_0 G_bar.M_s + grad_psi x J_s}dS
//
//   J_s = n x H_s,  M_s = n x E_s.
//
// The dyadic Green's function captures the (1/k^2)*grad-grad-psi term that
// the scalar form omits.  This is the fix for the former 66% Phase 2
// undershoot at R_obs/lambda << 1.
//

namespace {

// Per-face cache for the harmonic kernel.
struct HarmFaceCache {
    double cx, cy, cz;                      // centroid
    double nx, ny, nz;                      // outward unit normal
    double area;
    std::complex<double> Js[3];             // n x H_s
    std::complex<double> Ms[3];             // n x E_s
};

void build_harm_cache(
    const double* centroids, const double* normals,
    const double* areas,
    const double* E_re, const double* E_im,
    const double* H_re, const double* H_im,
    int n_faces, std::vector<HarmFaceCache>& cache)
{
    using C = std::complex<double>;
    cache.resize(n_faces);
    for (int t = 0; t < n_faces; ++t) {
        HarmFaceCache& fc = cache[t];
        fc.cx = centroids[t*3+0]; fc.cy = centroids[t*3+1]; fc.cz = centroids[t*3+2];
        fc.nx = normals[t*3+0]; fc.ny = normals[t*3+1]; fc.nz = normals[t*3+2];
        fc.area = areas[t];
        const C Ex(E_re[t*3+0], E_im[t*3+0]);
        const C Ey(E_re[t*3+1], E_im[t*3+1]);
        const C Ez(E_re[t*3+2], E_im[t*3+2]);
        const C Hx(H_re[t*3+0], H_im[t*3+0]);
        const C Hy(H_re[t*3+1], H_im[t*3+1]);
        const C Hz(H_re[t*3+2], H_im[t*3+2]);
        // J_s = n x H_s
        fc.Js[0] = fc.ny * Hz - fc.nz * Hy;
        fc.Js[1] = fc.nz * Hx - fc.nx * Hz;
        fc.Js[2] = fc.nx * Hy - fc.ny * Hx;
        // M_s = n x E_s.  The dyadic Green-function form already
        // contains the longitudinal surface-charge contribution, so do
        // not add separate -(n.E)grad_psi / -(n.H)grad_psi terms here.
        fc.Ms[0] = fc.ny * Ez - fc.nz * Ey;
        fc.Ms[1] = fc.nz * Ex - fc.nx * Ez;
        fc.Ms[2] = fc.nx * Ey - fc.ny * Ex;
    }
}

}  // anonymous namespace


void EvaluateHarmonic(
    const double* centroids,
    const double* normals,
    const double* areas,
    const double* E_re, const double* E_im,
    const double* H_re, const double* H_im,
    int n_faces,
    const double* obs,
    int n_obs,
    double omega,
    double* E_re_out, double* E_im_out,
    double* H_re_out, double* H_im_out,
    int n_threads)
{
    using C = std::complex<double>;

    std::memset(E_re_out, 0, sizeof(double) * 3 * (size_t)n_obs);
    std::memset(E_im_out, 0, sizeof(double) * 3 * (size_t)n_obs);
    std::memset(H_re_out, 0, sizeof(double) * 3 * (size_t)n_obs);
    std::memset(H_im_out, 0, sizeof(double) * 3 * (size_t)n_obs);

    if (n_faces <= 0 || n_obs <= 0) return;
    if (omega <= 0.0) return;   // caller should use EvaluateStaticH

    std::vector<HarmFaceCache> cache;
    build_harm_cache(centroids, normals, areas,
                      E_re, E_im, H_re, H_im, n_faces, cache);

    const double k_w = omega * std::sqrt(MU_0 * EPS_0);   // wave number, real
    const C jwmu(0.0, omega * MU_0);                       // jw mu_0
    const C jweps(0.0, omega * EPS_0);                     // jw eps_0
    const C jk(0.0, k_w);                                  // jk

    int requested_threads = (n_threads > 0)
        ? n_threads
        : ngcore::TaskManager::GetMaxThreads();
    if (requested_threads <= 0) requested_threads = 1;
    ngcore::RegionTaskManager rtm(requested_threads);

    ngcore::ParallelFor(ngcore::IntRange(n_obs), [&](size_t j_sz) {
        const int j = static_cast<int>(j_sz);
        const double ox = obs[j*3+0];
        const double oy = obs[j*3+1];
        const double oz = obs[j*3+2];
        C Ex_acc(0,0), Ey_acc(0,0), Ez_acc(0,0);
        C Hx_acc(0,0), Hy_acc(0,0), Hz_acc(0,0);

        for (int t = 0; t < n_faces; ++t) {
            const HarmFaceCache& fc = cache[t];
            const double Rx = ox - fc.cx;
            const double Ry = oy - fc.cy;
            const double Rz = oz - fc.cz;
            const double R2 = Rx*Rx + Ry*Ry + Rz*Rz;
            if (R2 < R_MIN * R_MIN) continue;
            const double R = std::sqrt(R2);
            const double inv_R = 1.0 / R;
            // R-hat
            const double rhx = Rx * inv_R, rhy = Ry * inv_R, rhz = Rz * inv_R;
            // psi = exp(-jkR) / (4 pi R)
            const C exp_mjkR = std::exp(C(0.0, -k_w * R));
            const C psi = exp_mjkR * (INV_FOUR_PI * inv_R);
            // grad_psi = -(jk + 1/R) psi R_hat
            const C scal_g = -(jk + inv_R) * psi;
            const C gx = scal_g * rhx, gy = scal_g * rhy, gz = scal_g * rhz;

            // Dyadic kernel pre-factors (h = jkR)
            const C h = C(0.0, k_w * R);   // jkR
            const C h2 = h * h;
            const C alpha = (h2 + h + 1.0) / h2;   // multiplies v
            const C beta  = (h2 + 3.0*h + 3.0) / h2;  // multiplies (v.R_hat) R_hat
            //   G_bar . v = psi (alpha v - beta (v.R_hat) R_hat)

            // ============= E formula =============
            // E_per_face = [-jw mu_0 G_bar . J_s + grad_psi x M_s] dS
            // (1) -jw mu_0 G_bar . J_s
            const C Jdotr = fc.Js[0]*rhx + fc.Js[1]*rhy + fc.Js[2]*rhz;
            const C Gx_J = psi * (alpha * fc.Js[0] - beta * Jdotr * rhx);
            const C Gy_J = psi * (alpha * fc.Js[1] - beta * Jdotr * rhy);
            const C Gz_J = psi * (alpha * fc.Js[2] - beta * Jdotr * rhz);
            const C E1x = -jwmu * Gx_J;
            const C E1y = -jwmu * Gy_J;
            const C E1z = -jwmu * Gz_J;
            // (2) grad_psi x M_s
            const C E2x = gy * fc.Ms[2] - gz * fc.Ms[1];
            const C E2y = gz * fc.Ms[0] - gx * fc.Ms[2];
            const C E2z = gx * fc.Ms[1] - gy * fc.Ms[0];
            const double dS = fc.area;
            Ex_acc += (E1x + E2x) * dS;
            Ey_acc += (E1y + E2y) * dS;
            Ez_acc += (E1z + E2z) * dS;

            // ============= H formula =============
            // H_per_face = [+jw eps_0 G_bar . M_s + grad_psi x J_s] dS
            const C Mdotr = fc.Ms[0]*rhx + fc.Ms[1]*rhy + fc.Ms[2]*rhz;
            const C Gx_M = psi * (alpha * fc.Ms[0] - beta * Mdotr * rhx);
            const C Gy_M = psi * (alpha * fc.Ms[1] - beta * Mdotr * rhy);
            const C Gz_M = psi * (alpha * fc.Ms[2] - beta * Mdotr * rhz);
            const C H1x = jweps * Gx_M;
            const C H1y = jweps * Gy_M;
            const C H1z = jweps * Gz_M;
            // grad_psi x J_s
            const C H2x = gy * fc.Js[2] - gz * fc.Js[1];
            const C H2y = gz * fc.Js[0] - gx * fc.Js[2];
            const C H2z = gx * fc.Js[1] - gy * fc.Js[0];
            Hx_acc += (H1x + H2x) * dS;
            Hy_acc += (H1y + H2y) * dS;
            Hz_acc += (H1z + H2z) * dS;
        }

        E_re_out[j*3+0] = Ex_acc.real();  E_im_out[j*3+0] = Ex_acc.imag();
        E_re_out[j*3+1] = Ey_acc.real();  E_im_out[j*3+1] = Ey_acc.imag();
        E_re_out[j*3+2] = Ez_acc.real();  E_im_out[j*3+2] = Ez_acc.imag();
        H_re_out[j*3+0] = Hx_acc.real();  H_im_out[j*3+0] = Hx_acc.imag();
        H_re_out[j*3+1] = Hy_acc.real();  H_im_out[j*3+1] = Hy_acc.imag();
        H_re_out[j*3+2] = Hz_acc.real();  H_im_out[j*3+2] = Hz_acc.imag();
    });
}

}}  // namespace radia::eqsrc
