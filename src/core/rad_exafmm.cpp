/*-------------------------------------------------------------------------
 *
 * File name:      rad_exafmm.cpp
 *
 * Project:        RADIA
 *
 * Description:    ExaFMM Integration for fast dipole field computation
 *
 * Author(s):      Claude Code (AI-assisted development)
 *
 * First release:  2025-12-30
 *
 * Copyright (C):  European Synchrotron Radiation Facility, France
 *
 *-------------------------------------------------------------------------
 *
 * This module provides:
 * 1. ExaFMM-t library integration (when RADIA_USE_EXAFMM is defined)
 * 2. Direct C++ dipole computation (always available as fallback)
 *
 * The magnetic scalar potential from a dipole m at position r' is:
 *   phi_m(r) = (1/4*pi) * (m . (r-r')) / |r-r'|^3
 *
 * The H-field is:
 *   H = -grad(phi_m) = (1/4*pi) * [3*(m.r)*r/r^5 - m/r^3]
 *
 *-------------------------------------------------------------------------*/

#include "rad_exafmm.h"
#include "rad_constants.h"
#include <cmath>
#include <vector>
#include <cstdint>

#include "rad_parallel.h"

// ExaFMM-t headers must be included OUTSIDE of any namespace
// to avoid std:: namespace conflicts with MSVC
#ifdef RADIA_USE_EXAFMM
#include "exafmm/dipole.h"
#include "exafmm/build_tree.h"
#include "exafmm/build_list.h"
#endif

namespace RadExaFMM {

//-------------------------------------------------------------------------
// Direct C++ implementation of dipole field computation
// This is always available and serves as fallback when ExaFMM is not linked
//-------------------------------------------------------------------------

/**
 * Compute H-field from a single dipole at one target point
 *
 * phi_m = (1/4*pi) * (m . r) / r^3
 * H = -grad(phi_m) = (1/4*pi) * [3*(m.r)*r/r^5 - m/r^3]
 */
inline void ComputeSingleDipoleH(
    double sx, double sy, double sz,
    double mx, double my, double mz,
    double tx, double ty, double tz,
    double& pot, double& gx, double& gy, double& gz,
    double thresh = 1e-15
) {
    // Vector from source to target
    double rx = tx - sx;
    double ry = ty - sy;
    double rz = tz - sz;

    double r2 = rx*rx + ry*ry + rz*rz;

    if (r2 < thresh * thresh) {
        return;  // Too close - skip
    }

    double r = std::sqrt(r2);
    double r3 = r2 * r;
    double r5 = r3 * r2;

    // m . r
    double m_dot_r = mx*rx + my*ry + mz*rz;

    // Scalar potential: phi = (m.r) / r^3 (1/4*pi applied at end)
    pot += m_dot_r / r3;

    // H = -grad(phi) = [3*(m.r)*r/r^5 - m/r^3]
    double coeff = 3.0 * m_dot_r / r5;
    gx += coeff * rx - mx / r3;
    gy += coeff * ry - my / r3;
    gz += coeff * rz - mz / r3;
}

/**
 * Compute vector potential A from a single dipole at one target point
 *
 * A = (mu_0/4*pi) * (m x r) / r^3
 *
 * where m is the magnetic dipole moment [A*m^2] and r is the vector
 * from the dipole to the observation point.
 *
 * Units: A has units of [T*m] = [kg/(A*s^2)] when:
 *   - mu_0 = 4*pi*10^-7 [H/m]
 *   - m is in [A*m^2]
 *   - r is in [m]
 */
inline void ComputeSingleDipoleA(
    double sx, double sy, double sz,
    double mx, double my, double mz,
    double tx, double ty, double tz,
    double& Ax, double& Ay, double& Az,
    double thresh = 1e-15
) {
    // Vector from source to target
    double rx = tx - sx;
    double ry = ty - sy;
    double rz = tz - sz;

    double r2 = rx*rx + ry*ry + rz*rz;

    if (r2 < thresh * thresh) {
        return;  // Too close - skip
    }

    double r = std::sqrt(r2);
    double r3 = r2 * r;

    // A = (m x r) / r^3 (mu_0/4*pi applied at end)
    // m x r = [my*rz - mz*ry, mz*rx - mx*rz, mx*ry - my*rx]
    Ax += (my*rz - mz*ry) / r3;
    Ay += (mz*rx - mx*rz) / r3;
    Az += (mx*ry - my*rx) / r3;
}

FMMResult ComputeDipoleFieldDirect(
    const double* sources,
    const double* dipoles,
    int64_t nsource,
    const double* targets,
    int64_t ntarget,
    double thresh
) {
    FMMResult result;
    result.error_code = 0;

    if (nsource <= 0 || ntarget <= 0) {
        result.error_code = -1;
        return result;
    }

    // Allocate output arrays
    result.pot.resize(ntarget, 0.0);
    result.gradx.resize(ntarget, 0.0);
    result.grady.resize(ntarget, 0.0);
    result.gradz.resize(ntarget, 0.0);

    // Direct O(N*M) computation with TaskManager parallelization
    ngcore::ParallelFor(ngcore::IntRange(ntarget), [&](size_t i) {
        double tx = targets[i * 3 + 0];
        double ty = targets[i * 3 + 1];
        double tz = targets[i * 3 + 2];

        double pot = 0.0;
        double gx = 0.0, gy = 0.0, gz = 0.0;

        for (int64_t j = 0; j < nsource; ++j) {
            double sx = sources[j * 3 + 0];
            double sy = sources[j * 3 + 1];
            double sz = sources[j * 3 + 2];

            double mx = dipoles[j * 3 + 0];
            double my = dipoles[j * 3 + 1];
            double mz = dipoles[j * 3 + 2];

            ComputeSingleDipoleH(sx, sy, sz, mx, my, mz, tx, ty, tz,
                                 pot, gx, gy, gz, thresh);
        }

        // Apply 1/(4*pi) factor
        result.pot[i] = pot * RadConst::INV_FOUR_PI;
        result.gradx[i] = gx * RadConst::INV_FOUR_PI;
        result.grady[i] = gy * RadConst::INV_FOUR_PI;
        result.gradz[i] = gz * RadConst::INV_FOUR_PI;
    });

    return result;
}

/**
 * Compute vector potential A from dipoles (direct O(N*M) method)
 *
 * A = (mu_0/4*pi) * sum_j [ m_j x (r - r_j) / |r - r_j|^3 ]
 *
 * @param sources   Dipole positions (3*nsource doubles)
 * @param dipoles   Dipole moments (3*nsource doubles) in A*m^2
 * @param nsource   Number of source dipoles
 * @param targets   Target positions (3*ntarget doubles)
 * @param ntarget   Number of target points
 * @param thresh    Distance threshold for skipping (default 1e-15)
 * @return          FMMResult with Ax, Ay, Az in gradx, grady, gradz
 */
FMMResult ComputeDipoleVectorPotentialDirect(
    const double* sources,
    const double* dipoles,
    int64_t nsource,
    const double* targets,
    int64_t ntarget,
    double thresh
) {
    FMMResult result;
    result.error_code = 0;

    if (nsource <= 0 || ntarget <= 0) {
        result.error_code = -1;
        return result;
    }

    // Allocate output arrays (Ax, Ay, Az stored in gradx, grady, gradz)
    result.pot.resize(ntarget, 0.0);    // Not used for A
    result.gradx.resize(ntarget, 0.0);  // Ax
    result.grady.resize(ntarget, 0.0);  // Ay
    result.gradz.resize(ntarget, 0.0);  // Az

    // Direct O(N*M) computation with TaskManager parallelization
    ngcore::ParallelFor(ngcore::IntRange(ntarget), [&](size_t i) {
        double tx = targets[i * 3 + 0];
        double ty = targets[i * 3 + 1];
        double tz = targets[i * 3 + 2];

        double Ax = 0.0, Ay = 0.0, Az = 0.0;

        for (int64_t j = 0; j < nsource; ++j) {
            double sx = sources[j * 3 + 0];
            double sy = sources[j * 3 + 1];
            double sz = sources[j * 3 + 2];

            double mx = dipoles[j * 3 + 0];
            double my = dipoles[j * 3 + 1];
            double mz = dipoles[j * 3 + 2];

            ComputeSingleDipoleA(sx, sy, sz, mx, my, mz, tx, ty, tz,
                                 Ax, Ay, Az, thresh);
        }

        // Apply mu_0/(4*pi) = 1e-7 factor
        result.gradx[i] = Ax * RadConst::MU_0_OVER_FOUR_PI;
        result.grady[i] = Ay * RadConst::MU_0_OVER_FOUR_PI;
        result.gradz[i] = Az * RadConst::MU_0_OVER_FOUR_PI;
    });

    return result;
}

#ifdef RADIA_USE_EXAFMM

//-------------------------------------------------------------------------
// ExaFMM-t library integration (when available)
//-------------------------------------------------------------------------

FMMResult ComputeDipoleField(
    double eps,
    const double* sources,
    const double* dipoles,
    int64_t nsource,
    const double* targets,
    int64_t ntarget
) {
    using namespace exafmm_t;

    FMMResult result;
    result.error_code = 0;

    if (nsource <= 0 || ntarget <= 0) {
        result.error_code = -1;
        return result;
    }

    result.pot.resize(ntarget, 0.0);
    result.gradx.resize(ntarget, 0.0);
    result.grady.resize(ntarget, 0.0);
    result.gradz.resize(ntarget, 0.0);

    // Threshold for using FMM vs direct (FMM has overhead for small N)
    const int64_t FMM_THRESHOLD = 1000;

    if (nsource < FMM_THRESHOLD || ntarget < FMM_THRESHOLD) {
        // Use direct computation for small problems
        return ComputeDipoleFieldDirect(sources, dipoles, nsource, targets, ntarget);
    }

    try {
        // Create source bodies with dipole moments
        Bodies<real_t> src_bodies(nsource);
        for (int64_t i = 0; i < nsource; i++) {
            src_bodies[i].ibody = static_cast<int>(i);
            src_bodies[i].X[0] = sources[i * 3 + 0];
            src_bodies[i].X[1] = sources[i * 3 + 1];
            src_bodies[i].X[2] = sources[i * 3 + 2];
            // For dipole kernel, we need vector source - store as 3 separate values
            // ExaFMM's standard kernels use scalar q, so we use gradient_P2P directly
            src_bodies[i].q = 1.0;  // Dummy - actual dipole handled in P2P
        }

        // Create target bodies
        Bodies<real_t> trg_bodies(ntarget);
        for (int64_t i = 0; i < ntarget; i++) {
            trg_bodies[i].ibody = static_cast<int>(i);
            trg_bodies[i].X[0] = targets[i * 3 + 0];
            trg_bodies[i].X[1] = targets[i * 3 + 1];
            trg_bodies[i].X[2] = targets[i * 3 + 2];
            trg_bodies[i].p = 0;
            trg_bodies[i].F[0] = 0;
            trg_bodies[i].F[1] = 0;
            trg_bodies[i].F[2] = 0;
        }

        // FMM parameters
        int P = 6;           // Expansion order (6-10 typical)
        int ncrit = 64;      // Max bodies per leaf

        // Create DipoleFmm instance
        DipoleFmm fmm(P, ncrit);

        // Build tree
        NodePtrs<real_t> leafs, nonleafs;
        Nodes<real_t> nodes;

        get_bounds(src_bodies, trg_bodies, fmm.x0, fmm.r0);
        nodes = build_tree(src_bodies, trg_bodies, leafs, nonleafs, fmm);

        // Build interaction lists
        init_rel_coord();
        build_list(nodes, fmm);
        fmm.M2L_setup(nonleafs);

        // Precompute matrices
        fmm.precompute();

        // Set source values (dipole moments) in leaf nodes
        for (auto* leaf : leafs) {
            int nsrc = leaf->isrcs.size();
            leaf->src_value.resize(nsrc * 3);
            for (int j = 0; j < nsrc; j++) {
                int orig_idx = leaf->isrcs[j];
                leaf->src_value[3*j+0] = dipoles[orig_idx * 3 + 0];
                leaf->src_value[3*j+1] = dipoles[orig_idx * 3 + 1];
                leaf->src_value[3*j+2] = dipoles[orig_idx * 3 + 2];
            }
            // Allocate target values (phi, Hx, Hy, Hz per target)
            int ntrg = leaf->itrgs.size();
            leaf->trg_value.resize(ntrg * 4, 0.0);
        }

        // FMM evaluation
        fmm.upward_pass(nodes, leafs);
        fmm.downward_pass(nodes, leafs);

        // Extract results back to original order
        for (auto* leaf : leafs) {
            for (size_t j = 0; j < leaf->itrgs.size(); j++) {
                int orig_idx = leaf->itrgs[j];
                result.pot[orig_idx]   = leaf->trg_value[4*j+0];
                result.gradx[orig_idx] = leaf->trg_value[4*j+1];
                result.grady[orig_idx] = leaf->trg_value[4*j+2];
                result.gradz[orig_idx] = leaf->trg_value[4*j+3];
            }
        }

    } catch (...) {
        // Fallback to direct computation on any error
        return ComputeDipoleFieldDirect(sources, dipoles, nsource, targets, ntarget);
    }

    return result;
}

bool IsAvailable() {
    return true;
}

const char* GetVersion() {
    return "ExaFMM-t (BSD-3-Clause) + Direct fallback";
}

#else // !RADIA_USE_EXAFMM

FMMResult ComputeDipoleField(
    double eps,
    const double* sources,
    const double* dipoles,
    int64_t nsource,
    const double* targets,
    int64_t ntarget
) {
    (void)eps;  // Suppress unused parameter warning
    return ComputeDipoleFieldDirect(sources, dipoles, nsource, targets, ntarget);
}

bool IsAvailable() {
    return false;  // ExaFMM library not linked
}

const char* GetVersion() {
    return "Direct dipole computation (ExaFMM not enabled)";
}

#endif // RADIA_USE_EXAFMM

//-------------------------------------------------------------------------
// Unified wrapper for vector potential A computation
// (Works with or without ExaFMM)
//-------------------------------------------------------------------------

FMMResult ComputeDipoleVectorPotential(
    double eps,
    const double* sources,
    const double* dipoles,
    int64_t nsource,
    const double* targets,
    int64_t ntarget
) {
    (void)eps;  // Currently unused (FMM for A not implemented yet)
    return ComputeDipoleVectorPotentialDirect(sources, dipoles, nsource, targets, ntarget);
}

} // namespace RadExaFMM
