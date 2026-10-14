/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file iterative_solver.hpp
 * @brief Base class for iterative linear solvers
 */

#ifndef SPARSESOLV_SOLVERS_ITERATIVE_SOLVER_HPP
#define SPARSESOLV_SOLVERS_ITERATIVE_SOLVER_HPP

#include "../core/types.hpp"
#include "../core/constants.hpp"
#include "../core/norms.hpp"
#include "../core/solver_config.hpp"
#include "../core/sparse_matrix_view.hpp"
#include "../core/preconditioner.hpp"
#include "../core/parallel.hpp"
#include <cmath>
#include <complex>
#include <algorithm>
#include <stdexcept>
#include <type_traits>

namespace sparsesolv {

/**
 * @brief Abstract base class for iterative linear solvers
 *
 * This class provides a common interface for iterative methods such as CG,
 * COCR, etc. It handles:
 * - Configuration management
 * - Convergence checking
 * - Divergence detection
 * - Residual history tracking
 * - Best result saving
 *
 * To implement a new solver:
 * 1. Override do_iterate() to implement the specific iteration scheme
 * 2. Optionally override name() to return a descriptive name
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class IterativeSolver {
public:
    virtual ~IterativeSolver() = default;

    /**
     * @brief Solve the linear system Ax = b
     *
     * @param A The system matrix
     * @param b Right-hand side vector
     * @param x Solution vector (initial guess on input, solution on output)
     * @param size System size
     * @param precond Preconditioner (optional, nullptr for no preconditioning)
     * @return SolverResult containing convergence info
     */
    SolverResult solve(
        const SparseMatrixView<Scalar>& A,
        const Scalar* b,
        Scalar* x,
        index_t size,
        const Preconditioner<Scalar>* precond = nullptr
    ) {
        config_.validate();
        size_ = size;
        A_ = &A;
        b_ = b;
        x_ = x;
        precond_ = precond;
        iteration_limit_ = config_.iteration_limit(size);

        allocate_work_vectors();

        // A zero right-hand side or an initial guess that already meets the
        // tolerance returns without iterating.
        if (prepare_iteration()) {
            return build_result(true, 0);
        }
        return do_iterate();
    }

    /**
     * @brief Solve with std::vector interface
     */
    SolverResult solve(
        const SparseMatrixView<Scalar>& A,
        const std::vector<Scalar>& b,
        std::vector<Scalar>& x,
        const Preconditioner<Scalar>* precond = nullptr
    ) {
        if (x.size() != b.size()) {
            x.resize(b.size());
        }
        return solve(A, b.data(), x.data(), static_cast<index_t>(b.size()), precond);
    }

    /// Set solver configuration
    void set_config(const SolverConfig& config) { config_ = config; }

    /// Get solver configuration
    const SolverConfig& config() const { return config_; }

    /// Get mutable configuration reference
    SolverConfig& config() { return config_; }

    /// Get the name of the solver
    virtual std::string name() const = 0;

protected:
    /**
     * @brief Implement the specific iteration scheme
     *
     * This method should perform the iteration until convergence or
     * max iterations. Use the helper methods for convergence checking.
     *
     * @return SolverResult with final state
     */
    virtual SolverResult do_iterate() = 0;

    /**
     * @brief Allocate work vectors needed by the solver
     *
     * Override this to allocate additional work vectors beyond r_, z_, p_, Ap_
     */
    virtual void allocate_work_vectors() {
        r_.resize(size_);
        z_.resize(size_);
        p_.resize(size_);
        Ap_.resize(size_);
    }

    /**
     * @brief Compute the initial residual and register the initial guess
     *
     * The initial guess is iteration 0 and the first best-result candidate.
     * @return true when no iteration is needed (zero right-hand side, or the
     *         initial relative residual is already below the tolerance)
     */
    bool prepare_iteration() {
        residual_history_.clear();
        bad_count_ = 0;
        best_iteration_ = 0;
        last_residual_ = 0.0;

        const double norm_b = robust_norm(b_, size_);
        if (!std::isfinite(norm_b)) {
            throw std::invalid_argument("SparseSolv: right-hand side is not finite");
        }
        if (norm_b == 0.0) {
            // A x = 0 has the solution x = 0
            std::fill(x_, x_ + size_, Scalar(0));
            normalizer_ = 1.0;
            best_residual_ = 0.0;
            if (config_.save_best_result) best_x_.assign(x_, x_ + size_);
            if (config_.save_residual_history) residual_history_.push_back(0.0);
            return true;
        }
        normalizer_ = norm_b;

        // r = b - A*x
        A_->multiply(x_, r_.data());
        parallel_for(size_, [&](index_t i) {
            r_[i] = b_[i] - r_[i];
        });
        const double rel0 = robust_norm(r_.data(), size_) / normalizer_;
        if (!std::isfinite(rel0)) {
            throw std::invalid_argument("SparseSolv: initial residual is not finite");
        }
        if (config_.save_residual_history) residual_history_.push_back(rel0);

        best_residual_ = rel0;
        last_residual_ = rel0;
        if (config_.save_best_result) best_x_.assign(x_, x_ + size_);
        return rel0 < config_.tolerance;
    }

    /**
     * @brief Apply preconditioner: z = M^{-1} * r
     */
    void apply_preconditioner() {
        if (precond_) {
            precond_->apply(r_.data(), z_.data(), size_);
        } else {
            // No preconditioning: z = r
            std::copy(r_.begin(), r_.end(), z_.begin());
        }
    }

    /**
     * @brief Apply preconditioner and compute dot(r, z) in one pass
     *
     * Computes z = M^{-1}*r and returns dot(r, z), fused to avoid
     * a separate kernel launch for the dot product.
     */
    Scalar apply_preconditioner_fused_dot() {
        if (precond_) {
            return precond_->apply_fused_dot(
                r_.data(), r_.data(), z_.data(), size_, config_.conjugate);
        } else {
            // No preconditioning: z = r, return dot(r, r)
            std::copy(r_.begin(), r_.end(), z_.begin());
            return dot_product(r_.data(), z_.data(), size_);
        }
    }

    /// Outcome of the per-iteration test
    enum class Step { Continue, Converged, Stop };

    /**
     * @brief Record iteration `iter` (1-based) and decide whether to stop
     *
     * rel = ||r|| / ||b|| with the recursive residual.  A new minimum
     * becomes the best iterate and resets the stagnation counter; a residual
     * below best * divergence_threshold also resets it; any other residual
     * increments it.  The solve converges when rel < tolerance and stops
     * when the counter exceeds divergence_count or rel is not finite.
     */
    Step check_convergence(double norm_r, int iter) {
        const double rel = norm_r / normalizer_;
        last_residual_ = rel;
        if (config_.save_residual_history) residual_history_.push_back(rel);
        if (!std::isfinite(rel)) return Step::Stop;

        if (rel < best_residual_) {
            best_residual_ = rel;
            best_iteration_ = iter;
            if (config_.save_best_result) std::copy(x_, x_ + size_, best_x_.begin());
            bad_count_ = 0;
        } else if (rel < best_residual_ * config_.divergence_threshold) {
            bad_count_ = 0;
        } else {
            ++bad_count_;
        }

        if (rel < config_.tolerance) return Step::Converged;
        if (config_.divergence_check == DivergenceCheck::StagnationCount
            && bad_count_ > config_.divergence_count) {
            return Step::Stop;
        }
        return Step::Continue;
    }

    /**
     * @brief Build the result and, with save_best_result, restore the best iterate
     */
    SolverResult build_result(bool converged, int iterations) {
        SolverResult result;
        result.converged = converged;
        result.iterations = iterations;
        if (config_.save_best_result) {
            std::copy(best_x_.begin(), best_x_.end(), x_);
            result.best_iteration = best_iteration_;
            result.final_residual = best_residual_;
        } else {
            result.best_iteration = iterations;
            result.final_residual = last_residual_;
        }
        result.residual_history = std::move(residual_history_);
        return result;
    }


    /**
     * @brief Compute dot product of two vectors
     *
     * When config_.conjugate is false (default):
     *   Unconjugated dot product (a^T * b) for complex-symmetric systems (A^T = A).
     *
     * When config_.conjugate is true:
     *   Conjugated dot product (a^H * b) for Hermitian systems (A^H = A).
     *   Uses std::conj(a[i]) * b[i].
     *
     * For real types, conjugation has no effect.
     */
    Scalar dot_product(const Scalar* a, const Scalar* b, index_t size) const {
        if constexpr (std::is_same_v<Scalar, double>) {
            // Real: conjugation has no effect
            return parallel_reduce_sum<Scalar>(size, [a, b](index_t i) {
                return a[i] * b[i];
            });
        } else {
            if (config_.conjugate) {
                // Hermitian: a^H * b = sum( conj(a[i]) * b[i] )
                return parallel_reduce_sum<Scalar>(size, [a, b](index_t i) {
                    return std::conj(a[i]) * b[i];
                });
            } else {
                // Complex-symmetric: a^T * b = sum( a[i] * b[i] )
                return parallel_reduce_sum<Scalar>(size, [a, b](index_t i) {
                    return a[i] * b[i];
                });
            }
        }
    }

    // Configuration
    SolverConfig config_;

    // Problem data (set during solve)
    index_t size_ = 0;
    const SparseMatrixView<Scalar>* A_ = nullptr;
    const Scalar* b_ = nullptr;
    Scalar* x_ = nullptr;
    const Preconditioner<Scalar>* precond_ = nullptr;

    // Work vectors
    std::vector<Scalar> r_;   // Residual
    std::vector<Scalar> z_;   // Preconditioned residual
    std::vector<Scalar> p_;   // Search direction
    std::vector<Scalar> Ap_;  // Matrix-vector product A*p

    // Convergence tracking
    int iteration_limit_ = 0;
    double normalizer_ = 1.0;       // ||b||
    double last_residual_ = 0.0;    // relative residual of the last iterate
    std::vector<double> residual_history_;

    // Best result tracking (the residual is tracked even when the iterate is not saved)
    std::vector<Scalar> best_x_;
    double best_residual_ = 0.0;
    int best_iteration_ = 0;

    // Stagnation counter
    int bad_count_ = 0;
};

} // namespace sparsesolv

#endif // SPARSESOLV_SOLVERS_ITERATIVE_SOLVER_HPP
