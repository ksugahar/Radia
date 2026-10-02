/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file sparsesolv.hpp
/// @brief Main header — includes all SparseSolv components

#ifndef SPARSESOLV_HPP
#define SPARSESOLV_HPP

// Core components
#include "core/types.hpp"
#include "core/constants.hpp"
#include "core/solver_config.hpp"
#include "core/sparse_matrix_view.hpp"
#include "core/preconditioner.hpp"
#include "core/abmc_ordering.hpp"
#include "core/norms.hpp"

// Preconditioners
#include "preconditioners/ic_preconditioner.hpp"
// Solvers
#include "solvers/iterative_solver.hpp"
#include "solvers/cg_solver.hpp"
#include "solvers/cocr_solver.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace sparsesolv {

struct Version {
    static constexpr int major = 2;
    static constexpr int minor = 7;
    static constexpr int patch = 0;

    static std::string string() {
        return std::to_string(major) + "." +
               std::to_string(minor) + "." +
               std::to_string(patch);
    }
};

namespace detail {

/// conjugate=True has no meaning for ICCG (complex-symmetric IC factor) or COCR
inline void reject_conjugate(const SolverConfig& config, const char* why) {
    if (config.conjugate) {
        throw std::invalid_argument(std::string(why) + "; conjugate=True is not supported");
    }
}

/// ||b - A x|| / ||b|| on the given system; ||A x|| when b is exactly zero
template<typename Scalar>
double relative_true_residual(const SparseMatrixView<Scalar>& A, const Scalar* b,
                              const Scalar* x, index_t n) {
    std::vector<Scalar> r(n);
    A.multiply(x, r.data());
    parallel_for(n, [&](index_t i) { r[i] = b[i] - r[i]; });
    const double rn = robust_norm(r.data(), n);
    const double bn = robust_norm(b, n);
    return bn == 0.0 ? rn : rn / bn;
}

/**
 * @brief Solve a (possibly scaled) system: normalize b, test the initial
 *        guess, and only then build a preconditioner through `inner`
 *
 * b is divided by an exact power of two near max|b_i| (x0 likewise and x is
 * multiplied back), which leaves every relative residual unchanged and keeps
 * the iteration's norms clear of underflow and overflow.  An exactly zero b
 * returns x = 0, and an initial guess with ||b - A x0|| / ||b|| < tol returns
 * without iterating; neither builds the IC factor.
 */
template<typename Scalar, typename Inner>
SolverResult solve_system(const SparseMatrixView<Scalar>& A, const Scalar* b, Scalar* x,
                          index_t n, const SolverConfig& config, Inner&& inner) {
    const double bmax = max_component(b, n);
    if (!std::isfinite(bmax)) {
        throw std::invalid_argument("SparseSolv: right-hand side is not finite");
    }
    SolverResult result;
    if (bmax == 0.0) {
        std::fill(x, x + n, Scalar(0));
        result.converged = true;
        if (config.save_residual_history) result.residual_history.push_back(0.0);
        return result;
    }
    // Divide b and x0 by 2^k on their exponents (exact; no reciprocal that
    // could overflow for a subnormal max|b_i|)
    const int k = scale_exponent(bmax);
    std::vector<Scalar> bw(n), xw(n), r(n);
    parallel_for(n, [&](index_t i) {
        bw[i] = scale_by_pow2(b[i], -k);
        xw[i] = scale_by_pow2(x[i], -k);
    });
    if (!std::isfinite(max_component(xw.data(), n))) {
        throw std::invalid_argument(
            "SparseSolv: the initial guess divided by the right-hand-side scale is not finite");
    }
    A.multiply(xw.data(), r.data());
    parallel_for(n, [&](index_t i) { r[i] = bw[i] - r[i]; });
    const double rel0 = robust_norm(r.data(), n) / robust_norm(bw.data(), n);
    if (!std::isfinite(rel0)) {
        throw std::invalid_argument("SparseSolv: initial residual is not finite");
    }
    if (rel0 < config.tolerance) {
        result.converged = true;
        result.final_residual = rel0;
        if (config.save_residual_history) result.residual_history.push_back(rel0);
        return result;  // x is the initial guess
    }
    result = inner(A, bw.data(), xw.data(), config);
    parallel_for(n, [&](index_t i) { x[i] = scale_by_pow2(xw[i], k); });
    return result;
}

/**
 * @brief Validate, scale when configured, and solve through solve_system
 *
 * With diagonal_scaling, S = diag(1/sqrt|a_ii|) and the inner solver works on
 * (S A S) y = S b from y0 = S^{-1} x0; the result is x = S y.  Its stopping
 * test therefore uses the scaled system.  The inner configuration has
 * diagonal_scaling off.  Option and matrix validation (including the scale
 * factors) precede the zero right-hand side and initial-guess returns.
 * true_residual is always measured on A x = b.
 */
template<typename Scalar, typename Inner>
SolverResult solve_scaled(const SparseMatrixView<Scalar>& A, const Scalar* b, Scalar* x,
                          index_t n, const SolverConfig& config, Inner&& inner) {
    config.validate();
    if (A.rows() != n || A.cols() != n) {
        throw std::invalid_argument("SparseSolv: matrix size does not match the vectors");
    }
    if (!std::isfinite(max_component(b, n))) {
        throw std::invalid_argument("SparseSolv: right-hand side is not finite");
    }
    if (!std::isfinite(max_component(x, n))) {
        throw std::invalid_argument("SparseSolv: initial guess is not finite");
    }
    SolverConfig inner_config = config;
    inner_config.diagonal_scaling = false;
    SolverResult result;
    if (!config.diagonal_scaling) {
        result = solve_system(A, b, x, n, inner_config, inner);
    } else {
        const index_t* rp = A.row_ptr();
        const index_t* ci = A.col_idx();
        const Scalar* av = A.values();
        // s_i = 1/sqrt|a_ii|; a row with a zero diagonal takes
        // s_i = 1 / max_{j<i} |a_ij s_j| over its lower off-diagonal entries.
        std::vector<double> s(n);
        for (index_t i = 0; i < n; ++i) {
            double d = 0.0;
            bool found = false;
            for (index_t k = rp[i]; k < rp[i + 1]; ++k) {
                if (ci[k] == i) { d = std::abs(av[k]); found = true; break; }
            }
            if (!found || !std::isfinite(d)) {
                throw std::invalid_argument(
                    "SparseSolv diagonal scaling: missing or non-finite diagonal at row "
                    + std::to_string(i));
            }
            if (d > 0.0) {
                s[i] = 1.0 / std::sqrt(d);
                continue;
            }
            double m = 0.0;
            for (index_t k = rp[i]; k < rp[i + 1]; ++k) {
                if (ci[k] < i) m = std::max(m, std::abs(av[k]) * s[ci[k]]);
            }
            if (!(m > 0.0) || !std::isfinite(m)) {
                throw std::invalid_argument(
                    "SparseSolv diagonal scaling: zero diagonal without a lower "
                    "off-diagonal entry at row " + std::to_string(i));
            }
            s[i] = 1.0 / m;
        }
        std::vector<Scalar> values(rp[n]);
        parallel_for(n, [&](index_t i) {
            for (index_t k = rp[i]; k < rp[i + 1]; ++k) {
                // (a_ij s_i) s_j: the product s_i s_j alone overflows for tiny diagonals
                values[k] = av[k] * static_cast<Scalar>(s[i]) * static_cast<Scalar>(s[ci[k]]);
            }
        });
        std::vector<Scalar> bs(n), ys(n);
        parallel_for(n, [&](index_t i) {
            bs[i] = b[i] * static_cast<Scalar>(s[i]);
            ys[i] = x[i] / static_cast<Scalar>(s[i]);
        });
        SparseMatrixView<Scalar> scaled(n, n, rp, ci, values.data());
        result = solve_system(scaled, bs.data(), ys.data(), n, inner_config, inner);
        parallel_for(n, [&](index_t i) {
            x[i] = ys[i] * static_cast<Scalar>(s[i]);
        });
    }
    result.true_residual = relative_true_residual(A, b, x, n);
    return result;
}

/// ICCG on the (already scaled) system with the configured IC ordering
template<typename Scalar>
SolverResult iccg_unscaled(const SparseMatrixView<Scalar>& A, const Scalar* b, Scalar* x,
                           index_t size, const SolverConfig& config) {
    ICPreconditioner<Scalar> precond(config.shift_parameter);
    precond.set_config(config);
    precond.setup(A);

    CGSolver<Scalar> solver;
    solver.set_config(config);

    SolverResult result;
    if (precond.has_reordered_matrix()) {
        // Path 1: ABMC with abmc_reorder_spmv=true
        // CG runs entirely in ABMC-reordered space (legacy behavior)
        const auto& ord = precond.abmc_ordering();

        std::vector<Scalar> b_reord(size);
        for (index_t i = 0; i < size; ++i) b_reord[ord[i]] = b[i];

        std::vector<Scalar> x_reord(size);
        for (index_t i = 0; i < size; ++i) x_reord[ord[i]] = x[i];

        auto A_reord = precond.reordered_matrix_view();
        ICPrecondReorderedAdapter<Scalar> adapter(precond);
        result = solver.solve(A_reord, b_reord.data(), x_reord.data(),
                              size, &adapter);

        for (index_t i = 0; i < size; ++i) x[i] = x_reord[ord[i]];
    } else if (precond.has_rcm_matrix()) {
        // Path 2: RCM+ABMC split mode
        // SpMV uses RCM-reordered matrix (better cache locality than original)
        // Preconditioner handles RCM->ABMC permutation internally
        const auto& rcm_ord = precond.rcm_perm();

        std::vector<Scalar> b_rcm(size);
        for (index_t i = 0; i < size; ++i) b_rcm[rcm_ord[i]] = b[i];

        std::vector<Scalar> x_rcm(size);
        for (index_t i = 0; i < size; ++i) x_rcm[rcm_ord[i]] = x[i];

        auto A_rcm = precond.rcm_matrix_view();
        ICPrecondRCMABMCAdapter<Scalar> rcm_adapter(precond);
        result = solver.solve(A_rcm, b_rcm.data(), x_rcm.data(),
                              size, &rcm_adapter);

        for (index_t i = 0; i < size; ++i) x[i] = x_rcm[rcm_ord[i]];
    } else {
        // Path 3: Standard (default)
        // SpMV uses original matrix A (preserves FEM mesh cache locality)
        // ABMC permutation handled internally by preconditioner::apply_abmc()
        result = solver.solve(A, b, x, size, &precond);
    }
    result.actual_shift = precond.actual_shift();
    return result;
}

/// Unpreconditioned CG or COCR on the (already scaled) system
template<typename Solver, typename Scalar>
SolverResult krylov_unscaled(const SparseMatrixView<Scalar>& A, const Scalar* b, Scalar* x,
                             index_t size, const SolverConfig& config) {
    Solver solver;
    solver.set_config(config);
    return solver.solve(A, b, x, size, nullptr);
}

/// IC preconditioned COCR on the (already scaled) system
template<typename Scalar>
SolverResult iccocr_unscaled(const SparseMatrixView<Scalar>& A, const Scalar* b, Scalar* x,
                             index_t size, const SolverConfig& config) {
    ICPreconditioner<Scalar> precond(config.shift_parameter);
    precond.set_config(config);
    precond.setup(A);
    COCRSolver<Scalar> solver;
    solver.set_config(config);
    SolverResult result = solver.solve(A, b, x, size, &precond);
    result.actual_shift = precond.actual_shift();
    return result;
}

} // namespace detail

/**
 * @brief Solve A x = b with IC preconditioned CG
 *
 * Defaults: diagonal scaling S A S, auto IC shift from 1.0 in steps of 0.01
 * (pivot test Re(d) < 1e-6 |a_ii|, limit 5), stop when the recursive relative
 * residual of the scaled system is below tolerance, stagnation stop after
 * more than divergence_count iterations without progress, and the best
 * iterate (initial guess included) is returned.  x holds the initial guess.
 */
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_iccg(
    const SparseMatrixView<Scalar>& A,
    const Scalar* b,
    Scalar* x,
    index_t size,
    const SolverConfig& config = SolverConfig()
) {
    detail::reject_conjugate(config, "ICCG: the IC factor is complex symmetric");
    return detail::solve_scaled(A, b, x, size, config,
        [](const SparseMatrixView<Scalar>& M, const Scalar* rhs, Scalar* sol,
           const SolverConfig& cfg) {
            return detail::iccg_unscaled(M, rhs, sol, M.rows(), cfg);
        });
}

/// Solve A x = b with unpreconditioned CG (same stopping and scaling rules as ICCG)
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_cg(
    const SparseMatrixView<Scalar>& A, const Scalar* b, Scalar* x, index_t size,
    const SolverConfig& config = SolverConfig()
) {
    return detail::solve_scaled(A, b, x, size, config,
        [](const SparseMatrixView<Scalar>& M, const Scalar* rhs, Scalar* sol,
           const SolverConfig& cfg) {
            return detail::krylov_unscaled<CGSolver<Scalar>>(M, rhs, sol, M.rows(), cfg);
        });
}

/// Convenience function: Solve Ax=b using ICCG with std::vector
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_iccg(
    const SparseMatrixView<Scalar>& A,
    const std::vector<Scalar>& b,
    std::vector<Scalar>& x,
    const SolverConfig& config = SolverConfig()
) {
    if (x.size() != b.size()) {
        x.resize(b.size());
    }
    return solve_iccg(A, b.data(), x.data(), static_cast<index_t>(b.size()),
                      config);
}

/// Convenience function: Solve Ax=b using COCR (IC preconditioned)
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_cocr(
    const SparseMatrixView<Scalar>& A,
    const Scalar* b,
    Scalar* x,
    index_t size,
    const SolverConfig& config = SolverConfig()
) {
    detail::reject_conjugate(config, "COCR uses unconjugated products");
    return detail::solve_scaled(A, b, x, size, config,
        [](const SparseMatrixView<Scalar>& M, const Scalar* rhs, Scalar* sol,
           const SolverConfig& cfg) {
            return detail::iccocr_unscaled(M, rhs, sol, M.rows(), cfg);
        });
}

/// Solve A x = b with unpreconditioned COCR (same stopping and scaling rules as ICCG)
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_cocr_unpreconditioned(
    const SparseMatrixView<Scalar>& A, const Scalar* b, Scalar* x, index_t size,
    const SolverConfig& config = SolverConfig()
) {
    detail::reject_conjugate(config, "COCR uses unconjugated products");
    return detail::solve_scaled(A, b, x, size, config,
        [](const SparseMatrixView<Scalar>& M, const Scalar* rhs, Scalar* sol,
           const SolverConfig& cfg) {
            return detail::krylov_unscaled<COCRSolver<Scalar>>(M, rhs, sol, M.rows(), cfg);
        });
}

/// Convenience function: Solve Ax=b using COCR with std::vector
template<typename Scalar = double>
[[nodiscard]] inline SolverResult solve_cocr(
    const SparseMatrixView<Scalar>& A,
    const std::vector<Scalar>& b,
    std::vector<Scalar>& x,
    const SolverConfig& config = SolverConfig()
) {
    if (x.size() != b.size()) {
        x.resize(b.size());
    }
    return solve_cocr(A, b.data(), x.data(), static_cast<index_t>(b.size()), config);
}

} // namespace sparsesolv

#endif // SPARSESOLV_HPP
