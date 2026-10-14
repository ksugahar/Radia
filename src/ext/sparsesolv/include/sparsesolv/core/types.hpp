/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file types.hpp
 * @brief Basic type definitions for SparseSolv library
 */

#ifndef SPARSESOLV_CORE_TYPES_HPP
#define SPARSESOLV_CORE_TYPES_HPP

#include <cmath>
#include <cstdint>
#include <complex>
#include <stdexcept>
#include <vector>
#include <string>

namespace sparsesolv {

// Index type for sparse matrix indices
using index_t = std::int32_t;

// Complex type alias
using complex_t = std::complex<double>;

/**
 * @brief Result of an iterative solver
 */
struct SolverResult {
    bool converged = false;           ///< Whether the solver converged
    int iterations = 0;               ///< Number of iterations performed
    int best_iteration = 0;           ///< Iteration of the returned iterate (0 = initial guess)
    double final_residual = 0.0;      ///< Relative recursive residual of the returned iterate
                                      ///< (scaled system when diagonal scaling is on)
    double true_residual = 0.0;       ///< ||b - A x|| / ||b|| of the returned x, original system
    double actual_shift = 0.0;        ///< IC shift used (0 when no IC factor was applied)
    std::vector<double> residual_history;  ///< [initial, iteration 1, ...] (optional)

    /// Check if the solve was successful
    explicit operator bool() const { return converged; }
};

/**
 * @brief Divergence detection strategy
 */
enum class DivergenceCheck {
    None,             ///< No divergence check, run until max iterations
    StagnationCount   ///< Stop if residual stagnates for too many iterations
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_TYPES_HPP
