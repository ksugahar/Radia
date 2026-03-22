/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * @file cg_solver.hpp
 * @brief Conjugate Gradient (CG) solver
 */

#ifndef SPARSESOLV_SOLVERS_CG_SOLVER_HPP
#define SPARSESOLV_SOLVERS_CG_SOLVER_HPP

#include "iterative_solver.hpp"
#include "../core/constants.hpp"

namespace sparsesolv {

/**
 * @brief Preconditioned Conjugate Gradient (PCG) solver
 *
 * Solves symmetric positive definite linear systems Ax = b using the
 * preconditioned conjugate gradient method.
 *
 * The algorithm:
 * 1. r_0 = b - A*x_0
 * 2. z_0 = M^{-1} * r_0
 * 3. p_0 = z_0
 * 4. For k = 0, 1, 2, ...:
 *    a. alpha_k = (r_k, z_k) / (p_k, A*p_k)
 *    b. x_{k+1} = x_k + alpha_k * p_k
 *    c. r_{k+1} = r_k - alpha_k * A*p_k
 *    d. Check convergence
 *    e. z_{k+1} = M^{-1} * r_{k+1}
 *    f. beta_k = (r_{k+1}, z_{k+1}) / (r_k, z_k)
 *    g. p_{k+1} = z_{k+1} + beta_k * p_k
 *
 * Requirements:
 * - Matrix A must be symmetric positive definite
 * - Preconditioner M should also be SPD
 *
 * Example:
 * @code
 * CGSolver<double> solver;
 * solver.config().tolerance = 1e-10;
 * solver.config().max_iterations = 1000;
 *
 * ICPreconditioner<double> precond(1.05);
 * precond.setup(A);
 *
 * auto result = solver.solve(A, b, x, n, &precond);
 * if (result.converged) {
 *     std::cout << "Converged in " << result.iterations << " iterations\n";
 * }
 * @endcode
 *
 * @tparam Scalar The scalar type (double or complex<double>)
 */
template<typename Scalar = double>
class CGSolver : public IterativeSolver<Scalar> {
public:
    std::string name() const override { return "CG"; }

protected:
    SolverResult do_iterate() override {
        const index_t n = this->size_;
        auto& r = this->r_;
        auto& z = this->z_;
        auto& p = this->p_;
        auto& Ap = this->Ap_;
        Scalar* x = this->x_;
        const auto& config = this->config_;

        // Fused: z = M^{-1} * r AND rz_old = dot(r, z)
        Scalar rz_old = this->apply_preconditioner_fused_dot();

        // p = z
        std::copy(z.begin(), z.end(), p.begin());

        // Cache CSR pointers for fused SpMV+dot
        const auto* A_rowptr = this->A_->row_ptr();
        const auto* A_colidx = this->A_->col_idx();
        const auto* A_vals = this->A_->values();

        // Main iteration loop
        for (int iter = 0; iter < config.max_iterations; ++iter) {
            // Fused: Ap = A*p AND pAp = (p, Ap) in single pass
            // Avoids re-reading p[] and Ap[] for separate dot product
            Scalar pAp = parallel_reduce_sum<Scalar>(n, [&](index_t i) -> Scalar {
                Scalar s = Scalar(0);
                for (index_t k = A_rowptr[i]; k < A_rowptr[i + 1]; ++k)
                    s += A_vals[k] * p[A_colidx[k]];
                Ap[i] = s;
                if constexpr (!std::is_same_v<Scalar, double>) {
                    return config.conjugate ? std::conj(p[i]) * s : p[i] * s;
                } else {
                    return p[i] * s;
                }
            });

            // Avoid division by zero
            if (std::abs(pAp) < constants::BREAKDOWN_THRESHOLD) {
                // Numerical breakdown - check if already converged
                double norm_r = this->compute_norm(r.data(), n);
                if (this->check_convergence(norm_r, iter)) {
                    return this->build_result(true, iter + 1, norm_r);
                }
                return this->build_result(false, iter, norm_r);
            }

            Scalar alpha = rz_old / pAp;

            // Fused: x += alpha*p, r -= alpha*Ap, norm_r = ||r|| in single pass
            // Avoids re-reading r[] for separate norm computation
            double norm_r_sq = parallel_reduce_sum<double>(n, [&](index_t i) -> double {
                x[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
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

            // Fused: z = M^{-1} * r AND rz_new = dot(r, z) in one pass
            // Avoids a separate kernel launch for the dot product
            Scalar rz_new = this->apply_preconditioner_fused_dot();
            Scalar beta = rz_new / rz_old;
            rz_old = rz_new;

            // p = z + beta * p
            parallel_for(n, [&](index_t i) {
                p[i] = z[i] + beta * p[i];
            });
        }

        // Max iterations reached
        double final_norm = this->compute_norm(r.data(), n);
        return this->build_result(false, config.max_iterations, final_norm);
    }
};

} // namespace sparsesolv

#endif // SPARSESOLV_SOLVERS_CG_SOLVER_HPP
