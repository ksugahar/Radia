/*
 * rad_gmres.h
 *
 * Templated restarted GMRES(m) with RIGHT preconditioning, for real (double) systems.
 *
 * Shares the radia::blas BLAS helpers from rad_bicgstab.h and uses the SAME matvec/precond
 * lambda interface, so it is a drop-in alternative to radia::bicgstab::Solve for the
 * collocation-MMMM moment solver:
 *
 *   auto matvec  = [&](const double* x, double* y){ ... y = A x ... };
 *   auto precond = [&](const double* x, double* y){ ... y = M^-1 x ... };
 *   auto res = radia::gmres::Solve(n, matvec, precond, rhs, sol, tol, max_iter, restart);
 *
 * GMRES is needed (over BiCGSTAB) for NON-NORMAL / INDEFINITE systems -- in particular the KKT
 * saddle system  [[A, Q],[Q^T, 0]]  used for the loop-free collocation-MMMM solve (A kept INTACT;
 * Q = field-null loop basis).  BiCGSTAB breaks down on that saddle (verified); GMRES converges.
 *
 * Algorithm: Saad-Schultz 1986, restarted GMRES(m), modified Gram-Schmidt Arnoldi + Givens
 * rotations, right preconditioning (so the Hessenberg residual is the TRUE residual ||b - A x||).
 *
 * Part of Radia project.
 */

#ifndef RAD_GMRES_H
#define RAD_GMRES_H

#include "rad_bicgstab.h"   // radia::blas::{dot,nrm2,axpy,scal,copy}
#include <vector>
#include <cmath>
#include <algorithm>

namespace radia {
namespace gmres {

struct Result {
    int iterations;   // total inner iterations (matvecs) performed
    double residual;  // final relative residual ||b - A x|| / ||b||
    bool converged;   // whether tol was reached
};

/**
 * @brief Restarted GMRES(m) with right preconditioning.  Solves A x = b.
 * @param n        system dimension
 * @param matvec   callable void(const double* x, double* y): y = A x
 * @param precond  callable void(const double* x, double* y): y = M^-1 x
 * @param rhs      right-hand side (size n)
 * @param sol      solution (size n, initial guess on input, result on output)
 * @param tol      relative tolerance on ||b - A x|| / ||b||
 * @param max_iter maximum total inner iterations (across restarts)
 * @param restart  restart length m (Krylov subspace size before restart)
 */
template<typename MatVecFunc, typename PrecFunc>
Result Solve(int n,
             MatVecFunc matvec,
             PrecFunc precond,
             const double* rhs,
             double* sol,
             double tol,
             int max_iter,
             int restart = 40)
{
    using radia::blas::dot;
    using radia::blas::nrm2;
    using radia::blas::axpy;
    using radia::blas::scal;
    using radia::blas::copy;

    Result result = {0, 1.0, false};
    if(n <= 0) return result;

    const int m = std::max(1, std::min(restart, std::max(1, max_iter)));
    std::vector<std::vector<double>> V((size_t)m + 1, std::vector<double>((size_t)n));
    std::vector<double> H((size_t)(m + 1) * (size_t)m, 0.0);  // H[i + j*(m+1)] = h_{i,j}
    std::vector<double> cs((size_t)m, 0.0), sn((size_t)m, 0.0), g((size_t)m + 1, 0.0);
    std::vector<double> w((size_t)n), z((size_t)n), tmp((size_t)n);

    double bnorm = nrm2(n, rhs);
    if(bnorm < 1e-30) bnorm = 1.0;

    int total_iter = 0;
    while(total_iter < max_iter && !result.converged)
    {
        // r0 = b - A x  (V[0] = normalized residual)
        matvec(sol, w.data());
        copy(n, rhs, V[0].data());
        axpy(n, -1.0, w.data(), V[0].data());
        double beta = nrm2(n, V[0].data());
        result.residual = beta / bnorm;
        if(result.residual < tol) { result.converged = true; break; }
        scal(n, 1.0 / beta, V[0].data());

        std::fill(g.begin(), g.end(), 0.0);
        g[0] = beta;

        int j = 0;
        for(; j < m && total_iter < max_iter; ++j, ++total_iter)
        {
            // right preconditioning: z = M^-1 v_j ; w = A z
            precond(V[(size_t)j].data(), z.data());
            matvec(z.data(), w.data());

            // modified Gram-Schmidt Arnoldi
            for(int i = 0; i <= j; ++i)
            {
                double hij = dot(n, w.data(), V[(size_t)i].data());
                H[(size_t)i + (size_t)j * (m + 1)] = hij;
                axpy(n, -hij, V[(size_t)i].data(), w.data());
            }
            double hjp = nrm2(n, w.data());
            H[(size_t)(j + 1) + (size_t)j * (m + 1)] = hjp;
            if(hjp > 1e-30)
            {
                copy(n, w.data(), V[(size_t)j + 1].data());
                scal(n, 1.0 / hjp, V[(size_t)j + 1].data());
            }

            // apply previous Givens rotations to the new column j
            for(int i = 0; i < j; ++i)
            {
                double h_i  = H[(size_t)i + (size_t)j * (m + 1)];
                double h_i1 = H[(size_t)(i + 1) + (size_t)j * (m + 1)];
                H[(size_t)i + (size_t)j * (m + 1)]       =  cs[(size_t)i] * h_i + sn[(size_t)i] * h_i1;
                H[(size_t)(i + 1) + (size_t)j * (m + 1)] = -sn[(size_t)i] * h_i + cs[(size_t)i] * h_i1;
            }

            // new Givens rotation to eliminate H[j+1, j]
            double h_jj = H[(size_t)j + (size_t)j * (m + 1)];
            double h_j1 = H[(size_t)(j + 1) + (size_t)j * (m + 1)];
            double denom = std::sqrt(h_jj * h_jj + h_j1 * h_j1);
            if(denom < 1e-30) denom = 1e-30;
            cs[(size_t)j] = h_jj / denom;
            sn[(size_t)j] = h_j1 / denom;
            H[(size_t)j + (size_t)j * (m + 1)] = cs[(size_t)j] * h_jj + sn[(size_t)j] * h_j1;
            H[(size_t)(j + 1) + (size_t)j * (m + 1)] = 0.0;

            // update the RHS rotation vector g
            double gt = cs[(size_t)j] * g[(size_t)j] + sn[(size_t)j] * g[(size_t)j + 1];
            g[(size_t)j + 1] = -sn[(size_t)j] * g[(size_t)j] + cs[(size_t)j] * g[(size_t)j + 1];
            g[(size_t)j] = gt;

            result.residual = std::abs(g[(size_t)j + 1]) / bnorm;
            result.iterations = total_iter + 1;
            if(result.residual < tol) { ++j; result.converged = true; break; }
            if(std::isnan(result.residual) || std::isinf(result.residual)) { ++j; break; }
        }

        // back-substitution: solve H(0:k-1, 0:k-1) y = g(0:k-1)
        const int k = j;
        std::vector<double> y((size_t)std::max(1, k), 0.0);
        for(int i = k - 1; i >= 0; --i)
        {
            double s = g[(size_t)i];
            for(int l = i + 1; l < k; ++l) s -= H[(size_t)i + (size_t)l * (m + 1)] * y[(size_t)l];
            double d = H[(size_t)i + (size_t)i * (m + 1)];
            y[(size_t)i] = (std::abs(d) > 1e-30) ? s / d : 0.0;
        }

        // x <- x + M^-1 (V(:,0:k-1) y)   (right preconditioning)
        std::fill(tmp.begin(), tmp.end(), 0.0);
        for(int i = 0; i < k; ++i) axpy(n, y[(size_t)i], V[(size_t)i].data(), tmp.data());
        precond(tmp.data(), z.data());
        axpy(n, 1.0, z.data(), sol);

        if(std::isnan(result.residual) || std::isinf(result.residual)) break;
    }

    return result;
}

} // namespace gmres
} // namespace radia

#endif // RAD_GMRES_H
