/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file cocr_solver.hpp
 * @brief COCR (Conjugate Orthogonal Conjugate Residual) solver
 *        for complex-symmetric systems (A^T = A)
 *
 * Reference: Sogabe & Zhang (2007), "A COCR method for solving complex
 * symmetric linear systems", J. Comput. Appl. Math., 199(2), 297-303.
 */

#ifndef SPARSESOLV_SOLVERS_COCR_SOLVER_HPP
#define SPARSESOLV_SOLVERS_COCR_SOLVER_HPP

#include "iterative_solver.hpp"
#include "../core/constants.hpp"

namespace sparsesolv {

/**
 * @brief Preconditioned COCR solver for complex-symmetric systems
 *
 * For A^T = A (complex-symmetric, NOT Hermitian), COCR is the optimal
 * short-recurrence Krylov solver. It uses unconjugated inner products
 * (x^T y, not x^H y) and minimizes ||A r~||_2 for smoother convergence
 * than COCG.
 *
 * Cost per iteration: 1 MatVec + 1 Preconditioner apply (same as CG).
 * Memory: 6 vectors (vs 4 for CG).
 *
 * Algorithm:
 * 1. r_0 = b - A*x_0
 * 2. rt_0 = M^{-1} * r_0
 * 3. p_0 = rt_0, q_0 = A*p_0
 * 4. rho_0 = rt_0^T * q_0  (unconjugated)
 * 5. For k = 0, 1, 2, ...:
 *    a. qt_k = M^{-1} * q_k
 *    b. mu_k = qt_k^T * q_k
 *    c. alpha_k = rho_k / mu_k
 *    d. x_{k+1} = x_k + alpha_k * p_k
 *    e. r_{k+1} = r_k - alpha_k * q_k
 *    f. rt_{k+1} = rt_k - alpha_k * qt_k  (vector update, NOT precond!)
 *    g. Check convergence
 *    h. t_{k+1} = A * rt_{k+1}
 *    i. rho_{k+1} = rt_{k+1}^T * t_{k+1}
 *    j. beta_k = rho_{k+1} / rho_k
 *    k. p_{k+1} = rt_{k+1} + beta_k * p_k
 *    l. q_{k+1} = t_{k+1} + beta_k * q_k
 *
 * @note For COCG (Conjugate Orthogonal Conjugate Gradient), use
 *       CGSolver with conjugate=false. COCG is mathematically identical.
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class COCRSolver : public IterativeSolver<Scalar> {
public:
    std::string name() const override { return "COCR"; }

protected:
    void allocate_work_vectors() override {
        IterativeSolver<Scalar>::allocate_work_vectors();
        qt_.resize(this->size_);
        t_.resize(this->size_);
    }

    SolverResult do_iterate() override {
        const index_t n = this->size_;
        auto& r = this->r_;       // residual
        auto& rt = this->z_;      // preconditioned residual (reuse z_)
        auto& p = this->p_;       // search direction
        auto& q = this->Ap_;      // A*p (reuse Ap_)
        Scalar* x = this->x_;
        const auto& config = this->config_;

        // rt = M^{-1} * r  (initial preconditioned residual)
        this->apply_preconditioner();

        // p = rt
        std::copy(rt.begin(), rt.end(), p.begin());

        // q = A * p
        this->A_->multiply(p.data(), q.data());

        // rho = rt^T * q  (unconjugated inner product)
        Scalar rho = unconjugated_dot(rt.data(), q.data(), n);

        // Cache CSR pointers for fused SpMV+dot
        const auto* A_rowptr = this->A_->row_ptr();
        const auto* A_colidx = this->A_->col_idx();
        const auto* A_vals = this->A_->values();

        // Main iteration loop
        for (int iter = 0; iter < config.max_iterations; ++iter) {
            // qt = M^{-1} * q  (preconditioner apply)
            if (this->precond_) {
                this->precond_->apply(q.data(), qt_.data(), n);
            } else {
                std::copy(q.begin(), q.end(), qt_.begin());
            }

            // mu = qt^T * q  (unconjugated)
            Scalar mu = unconjugated_dot(qt_.data(), q.data(), n);

            // Breakdown check
            if (std::abs(mu) < constants::BREAKDOWN_THRESHOLD) {
                double norm_r = this->compute_norm(r.data(), n);
                if (this->check_convergence(norm_r, iter)) {
                    return this->build_result(true, iter + 1, norm_r);
                }
                return this->build_result(false, iter, norm_r);
            }

            Scalar alpha = rho / mu;

            // x += alpha * p
            // r -= alpha * q
            // rt -= alpha * qt
            // norm_r = ||r||
            double norm_r_sq = parallel_reduce_sum<double>(n, [&](index_t i) -> double {
                x[i] += alpha * p[i];
                r[i] -= alpha * q[i];
                rt[i] -= alpha * qt_[i];
                return std::norm(r[i]);  // |r[i]|^2
            });
            double norm_r = std::sqrt(norm_r_sq);

            // Check convergence
            if (this->check_convergence(norm_r, iter)) {
                return this->build_result(true, iter + 1, norm_r);
            }

            // Check divergence
            if (this->is_diverging()) {
                return this->build_result(false, iter + 1, norm_r);
            }

            // Fused: t = A * rt AND rho_new = rt^T * t
            Scalar rho_new = parallel_reduce_sum<Scalar>(n, [&](index_t i) -> Scalar {
                Scalar s = Scalar(0);
                for (index_t k = A_rowptr[i]; k < A_rowptr[i + 1]; ++k)
                    s += A_vals[k] * rt[A_colidx[k]];
                t_[i] = s;
                return rt[i] * s;  // unconjugated
            });

            // Breakdown check
            if (std::abs(rho) < constants::BREAKDOWN_THRESHOLD) {
                return this->build_result(false, iter + 1, norm_r);
            }

            Scalar beta = rho_new / rho;
            rho = rho_new;

            // p = rt + beta * p
            // q = t + beta * q
            parallel_for(n, [&](index_t i) {
                p[i] = rt[i] + beta * p[i];
                q[i] = t_[i] + beta * q[i];
            });
        }

        // Max iterations reached
        double final_norm = this->compute_norm(r.data(), n);
        return this->build_result(false, config.max_iterations, final_norm);
    }

private:
    /**
     * @brief Unconjugated dot product: sum(a[i] * b[i])
     *
     * COCR always uses unconjugated inner products regardless of
     * config_.conjugate setting. This is fundamental to the algorithm
     * for complex-symmetric systems.
     */
    Scalar unconjugated_dot(const Scalar* a, const Scalar* b, index_t size) const {
        return parallel_reduce_sum<Scalar>(size, [a, b](index_t i) -> Scalar {
            return a[i] * b[i];  // no conjugation
        });
    }

    std::vector<Scalar> qt_;  // M^{-1} * q
    std::vector<Scalar> t_;   // A * rt
};

} // namespace sparsesolv

#endif // SPARSESOLV_SOLVERS_COCR_SOLVER_HPP
