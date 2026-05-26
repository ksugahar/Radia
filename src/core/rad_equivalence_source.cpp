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
#include <cstring>
#include <vector>

#include <core/taskmanager.hpp>

namespace radia { namespace eqsrc {

namespace {

constexpr double INV_FOUR_PI = 1.0 / (4.0 * 3.14159265358979323846);
constexpr double R_MIN       = 1.0e-15;     // singular-distance guard

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

}}  // namespace radia::eqsrc
