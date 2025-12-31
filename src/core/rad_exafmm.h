/*-------------------------------------------------------------------------
 *
 * File name:      rad_exafmm.h
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
 * License:        BSD-3-Clause (ExaFMM-t is BSD-3-Clause licensed)
 *
 *-------------------------------------------------------------------------
 *
 * This module provides:
 * 1. ExaFMM-t library integration for O(N log N) field computation
 * 2. Direct C++ dipole computation as fallback (O(N*M))
 *
 * The magnetic scalar potential from a dipole m at position r' is:
 *   phi_m(r) = (1/4*pi) * (m . (r-r')) / |r-r'|^3
 *
 * The H-field is:
 *   H = -grad(phi_m) = (1/4*pi) * [3*(m.r)*r/r^5 - m/r^3]
 *
 *-------------------------------------------------------------------------*/

#ifndef RAD_EXAFMM_H
#define RAD_EXAFMM_H

#include <vector>
#include <cstdint>

namespace RadExaFMM {

/**
 * Dipole field computation result
 */
struct FMMResult {
    std::vector<double> pot;      // Scalar potential at targets
    std::vector<double> gradx;    // H-field x component
    std::vector<double> grady;    // H-field y component
    std::vector<double> gradz;    // H-field z component
    int error_code;               // 0 = success, -1 = invalid input
};

/**
 * Compute H-field at target points from magnetic dipoles
 *
 * Uses ExaFMM library if available (RADIA_USE_EXAFMM defined),
 * otherwise falls back to direct O(N*M) computation.
 *
 * @param eps         FMM tolerance (e.g., 1e-6). Ignored in direct computation.
 * @param sources     Source positions [x1,y1,z1, x2,y2,z2, ...] (nsource*3 values)
 * @param dipoles     Dipole moments [mx1,my1,mz1, mx2,my2,mz2, ...] (nsource*3 values)
 * @param nsource     Number of sources
 * @param targets     Target positions [x1,y1,z1, x2,y2,z2, ...] (ntarget*3 values)
 * @param ntarget     Number of targets
 * @return            FMMResult with potential and gradient (H-field)
 */
FMMResult ComputeDipoleField(
    double eps,
    const double* sources,
    const double* dipoles,
    int64_t nsource,
    const double* targets,
    int64_t ntarget
);

/**
 * Direct O(N*M) computation of dipole field (H-field and scalar potential)
 *
 * Always available. OpenMP parallelized.
 *
 * @param sources     Source positions [x1,y1,z1, ...]
 * @param dipoles     Dipole moments [mx1,my1,mz1, ...]
 * @param nsource     Number of sources
 * @param targets     Target positions [x1,y1,z1, ...]
 * @param ntarget     Number of targets
 * @param thresh      Distance threshold (skip contribution if |r| < thresh)
 * @return            FMMResult with potential (phi) and gradient (Hx, Hy, Hz)
 */
FMMResult ComputeDipoleFieldDirect(
    const double* sources,
    const double* dipoles,
    int64_t nsource,
    const double* targets,
    int64_t ntarget,
    double thresh = 1e-15
);

/**
 * Direct O(N*M) computation of vector potential A from dipoles
 *
 * A = (mu_0/4*pi) * sum_j [ m_j x (r - r_j) / |r - r_j|^3 ]
 *
 * @param sources     Source positions [x1,y1,z1, ...] in meters
 * @param dipoles     Dipole moments [mx1,my1,mz1, ...] in A*m^2
 * @param nsource     Number of sources
 * @param targets     Target positions [x1,y1,z1, ...] in meters
 * @param ntarget     Number of targets
 * @param thresh      Distance threshold (skip contribution if |r| < thresh)
 * @return            FMMResult with Ax, Ay, Az in gradx, grady, gradz [T*m]
 */
FMMResult ComputeDipoleVectorPotentialDirect(
    const double* sources,
    const double* dipoles,
    int64_t nsource,
    const double* targets,
    int64_t ntarget,
    double thresh = 1e-15
);

/**
 * Compute vector potential A at target points from magnetic dipoles
 *
 * Uses direct O(N*M) computation with OpenMP parallelization.
 *
 * @param eps         FMM tolerance (currently unused, for future FMM support)
 * @param sources     Source positions [x1,y1,z1, x2,y2,z2, ...] in meters
 * @param dipoles     Dipole moments [mx1,my1,mz1, ...] in A*m^2
 * @param nsource     Number of sources
 * @param targets     Target positions [x1,y1,z1, x2,y2,z2, ...] in meters
 * @param ntarget     Number of targets
 * @return            FMMResult with Ax, Ay, Az in gradx, grady, gradz [T*m]
 */
FMMResult ComputeDipoleVectorPotential(
    double eps,
    const double* sources,
    const double* dipoles,
    int64_t nsource,
    const double* targets,
    int64_t ntarget
);

/**
 * Check if ExaFMM library is available
 * @return true if ExaFMM is available, false if using direct computation only
 */
bool IsAvailable();

/**
 * Get version/status string
 */
const char* GetVersion();

} // namespace RadExaFMM

#endif // RAD_EXAFMM_H
