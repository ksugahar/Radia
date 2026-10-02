/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file solver_config.hpp
 * @brief Unified configuration for all SparseSolv solvers
 */

#ifndef SPARSESOLV_CORE_SOLVER_CONFIG_HPP
#define SPARSESOLV_CORE_SOLVER_CONFIG_HPP

#include "types.hpp"

namespace sparsesolv {

/**
 * @brief Configuration structure for iterative solvers
 *
 * This struct centralizes all solver parameters that were previously
 * scattered across multiple setter methods. Use this to configure
 * tolerance, iteration limits, preconditioning parameters, and
 * convergence behavior.
 *
 * Example:
 * @code
 * SolverConfig config;
 * config.tolerance = 1e-10;
 * config.max_iterations = 1000;
 * config.shift_parameter = 1.05;  // IC shift; start value when auto_shift
 * @endcode
 */
struct SolverConfig {
    //--------------------------------------------------
    // Convergence criteria
    //--------------------------------------------------

    /// Relative tolerance: stop when ||r|| / ||b|| < tolerance (strict), with
    /// the recursive residual of the (diagonally scaled, when enabled) system.
    /// Must be positive.
    double tolerance = 1e-8;

    /// Maximum number of iterations; 0 means 2 * n. Must not be negative.
    int max_iterations = 0;

    //--------------------------------------------------
    // Preconditioning parameters
    //--------------------------------------------------

    /// Multiplicative shift of the IC diagonal, alpha * a_ii, applied to rows
    /// with Re(a_ii) > 0.  Fixed value, or the start value when auto_shift is
    /// on.  Must be >= 1.
    double shift_parameter = 1.0;

    /// Solve the symmetrically scaled system S A S y = S b, x = S y, with
    /// S = diag(1/sqrt|a_ii|); the stopping test then uses the scaled system.
    bool diagonal_scaling = true;

    //--------------------------------------------------
    // Auto-shift parameters for IC decomposition
    //--------------------------------------------------

    /// Restart the factorization with shift + shift_increment while a pivot
    /// of a row with Re(a_ii) > 0 has Re(d_i) < min_diagonal_threshold * |a_ii|
    /// and shift < max_shift_value.  A pivot still below it at the limit is
    /// an error.
    bool auto_shift = true;

    /// Additive shift step of the auto-shift search
    double shift_increment = 0.01;

    /// The auto-shift search stops increasing the shift once it reaches this value
    double max_shift_value = 5.0;

    /// Pivot threshold relative to |a_ii| of the factored (scaled) matrix
    double min_diagonal_threshold = 1e-6;

    //--------------------------------------------------
    // Divergence detection
    //--------------------------------------------------

    /// Strategy for detecting divergence
    DivergenceCheck divergence_check = DivergenceCheck::StagnationCount;

    /// A residual below best * divergence_threshold resets the stagnation counter
    double divergence_threshold = 10.0;

    /// Stop when more than this many consecutive iterations neither improve
    /// the best residual nor stay below best * divergence_threshold
    int divergence_count = 10;

    //--------------------------------------------------
    // Result saving options
    //--------------------------------------------------

    /// Return the iterate with the smallest residual (the initial guess
    /// included) instead of the last iterate
    bool save_best_result = true;

    /// Save residual history for analysis/debugging
    bool save_residual_history = false;

    //--------------------------------------------------
    // Complex inner product
    //--------------------------------------------------

    /// Use conjugated inner product (a^H * b) for Hermitian systems.
    /// Default false uses unconjugated (a^T * b) for complex-symmetric systems.
    /// Has no effect for real-valued problems.
    bool conjugate = false;

    //--------------------------------------------------
    // ABMC ordering parameters (for parallel triangular solves)
    //--------------------------------------------------

    /// Enable ABMC (Algebraic Block Multi-Color) ordering.
    /// When enabled, the preconditioner reorders the matrix to enable
    /// parallel triangular solves using a two-level hierarchy:
    /// colors (sequential) -> blocks (parallel) -> rows (sequential).
    bool use_abmc = false;

    /// Number of rows per block for ABMC ordering (block size).
    /// Larger blocks reduce parallelism overhead but may decrease
    /// the degree of parallelism. Typical values: 2-16.
    int abmc_block_size = 4;

    /// Number of colors for ABMC graph coloring.
    /// More colors allow finer-grained parallelism but increase
    /// the number of sequential synchronization points.
    /// The actual number may be increased if the graph requires it.
    int abmc_num_colors = 4;

    /// When true, CG runs entirely in ABMC-reordered space (SpMV uses
    /// reordered matrix). When false (default), CG uses the original
    /// matrix for SpMV and only the preconditioner operates in
    /// reordered space. False is usually faster because it preserves
    /// the FEM mesh ordering cache locality for SpMV.
    bool abmc_reorder_spmv = false;

    /// Enable RCM (Reverse Cuthill-McKee) preprocessing before ABMC.
    /// RCM reduces bandwidth, improving cache locality for both SpMV
    /// and triangular solves. ABMC then operates on the bandwidth-reduced
    /// matrix.
    bool abmc_use_rcm = false;

    /// Iteration limit for a system of size n (0 means 2 * n)
    int iteration_limit(index_t n) const {
        return max_iterations > 0 ? max_iterations : 2 * static_cast<int>(n);
    }

    /// Reject option values that have no defined meaning
    void validate() const {
        if (!(tolerance > 0.0) || !std::isfinite(tolerance))
            throw std::invalid_argument("SparseSolv: tolerance must be positive and finite");
        if (max_iterations < 0)
            throw std::invalid_argument("SparseSolv: max_iterations must be >= 0 (0 means 2*n)");
        if (!(shift_parameter >= 1.0) || !std::isfinite(shift_parameter))
            throw std::invalid_argument("SparseSolv: IC shift must be >= 1 and finite");
        if (!(shift_increment > 0.0) || !(max_shift_value >= 1.0)
            || !(min_diagonal_threshold > 0.0))
            throw std::invalid_argument("SparseSolv: invalid auto-shift parameters");
        if (!(divergence_threshold > 0.0) || divergence_count < 0)
            throw std::invalid_argument(
                "SparseSolv: divergence_threshold must be > 0 and divergence_count >= 0");
    }
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_SOLVER_CONFIG_HPP
