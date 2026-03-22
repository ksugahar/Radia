/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file constants.hpp
 * @brief Named numeric constants for SparseSolv library
 *
 * Centralizes threshold values used internally by solvers and preconditioners.
 * SolverConfig default values (tolerance, max_iterations, shift_parameter, etc.)
 * are NOT included here as they are part of the public API.
 */

#ifndef SPARSESOLV_CORE_CONSTANTS_HPP
#define SPARSESOLV_CORE_CONSTANTS_HPP

namespace sparsesolv {
namespace constants {

/// Threshold for detecting numerical breakdown in CG (pAp ~ 0)
/// and SGS-MRTR first-iteration (v_w ~ 0)
constexpr double BREAKDOWN_THRESHOLD = 1e-30;

/// Threshold for detecting denominator collapse in SGS-MRTR
/// general iterations (denom ~ 0)
constexpr double DENOMINATOR_BREAKDOWN = 1e-60;

/// Minimum absolute diagonal value for safe inversion in preconditioners
/// (Jacobi, SGS, IC diagonal checks)
constexpr double MIN_DIAGONAL_TOLERANCE = 1e-15;

} // namespace constants
} // namespace sparsesolv

#endif // SPARSESOLV_CORE_CONSTANTS_HPP
