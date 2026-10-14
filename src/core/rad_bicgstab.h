/*
 * rad_bicgstab.h
 *
 * Templated BiCGSTAB solver for real (double) and complex (std::complex<double>) types.
 *
 * Shares MKL BLAS infrastructure with the Radia linear solver paths.
 * BLAS dispatch via function overloading:
 *   double           -> cblas_ddot, cblas_dnrm2, cblas_daxpy, cblas_dscal, cblas_dcopy, cblas_dgemv
 *   complex<double>  -> cblas_zdotc_sub, cblas_dznrm2, cblas_zaxpy, cblas_zscal, cblas_zcopy, cblas_zgemv
 *
 * Usage (PEEC complex):
 *   auto matvec = [&](const Z* x, Z* y) { blas::gemv(n, n, Z(1,0), A, n, x, Z(0,0), y); };
 *   auto precond = [&](const Z* x, Z* y) { for(int i=0;i<n;i++) y[i] = diag_inv[i]*x[i]; };
 *   auto result = bicgstab::Solve<Z>(n, matvec, precond, rhs, sol, 1e-10, 1000);
 *
 * Usage (real-valued magnetic solve):
 *   auto matvec = [&](const double* x, double* y) { blas::gemv(n, n, -1.0, K, n, x, 0.0, y); ... };
 *   auto result = bicgstab::Solve<double>(n, matvec, precond, rhs, sol, tol, max_iter);
 *
 * Part of Radia project
 */

#ifndef RAD_BICGSTAB_H
#define RAD_BICGSTAB_H

#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

#ifdef HAVE_LAPACK
#include "mkl_cblas.h"
#endif

namespace radia {

// ============================================================================
// BLAS dispatch: overloaded functions for double and complex<double>
// ============================================================================

namespace blas {

// ---------- dot product ----------

#ifdef HAVE_LAPACK
inline double dot(int n, const double* x, const double* y) {
    return cblas_ddot(n, x, 1, y, 1);
}

inline std::complex<double> dot(int n, const std::complex<double>* x,
                                 const std::complex<double>* y) {
    std::complex<double> result;
    cblas_zdotc_sub(n, x, 1, y, 1, &result);  // Hermitian: conj(x) . y
    return result;
}
#else
inline double dot(int n, const double* x, const double* y) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += x[i] * y[i];
    return sum;
}

inline std::complex<double> dot(int n, const std::complex<double>* x,
                                 const std::complex<double>* y) {
    std::complex<double> sum(0, 0);
    for (int i = 0; i < n; ++i) sum += std::conj(x[i]) * y[i];
    return sum;
}
#endif

// ---------- Euclidean norm ----------

#ifdef HAVE_LAPACK
inline double nrm2(int n, const double* x) {
    return cblas_dnrm2(n, x, 1);
}

inline double nrm2(int n, const std::complex<double>* x) {
    return cblas_dznrm2(n, x, 1);
}
#else
inline double nrm2(int n, const double* x) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += x[i] * x[i];
    return std::sqrt(sum);
}

inline double nrm2(int n, const std::complex<double>* x) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += std::norm(x[i]);  // |x[i]|^2
    return std::sqrt(sum);
}
#endif

// ---------- axpy: y += alpha * x ----------

#ifdef HAVE_LAPACK
inline void axpy(int n, double alpha, const double* x, double* y) {
    cblas_daxpy(n, alpha, x, 1, y, 1);
}

inline void axpy(int n, std::complex<double> alpha,
                 const std::complex<double>* x, std::complex<double>* y) {
    cblas_zaxpy(n, &alpha, x, 1, y, 1);
}
#else
inline void axpy(int n, double alpha, const double* x, double* y) {
    for (int i = 0; i < n; ++i) y[i] += alpha * x[i];
}

inline void axpy(int n, std::complex<double> alpha,
                 const std::complex<double>* x, std::complex<double>* y) {
    for (int i = 0; i < n; ++i) y[i] += alpha * x[i];
}
#endif

// ---------- scal: x *= alpha ----------

#ifdef HAVE_LAPACK
inline void scal(int n, double alpha, double* x) {
    cblas_dscal(n, alpha, x, 1);
}

inline void scal(int n, std::complex<double> alpha, std::complex<double>* x) {
    cblas_zscal(n, &alpha, x, 1);
}
#else
inline void scal(int n, double alpha, double* x) {
    for (int i = 0; i < n; ++i) x[i] *= alpha;
}

inline void scal(int n, std::complex<double> alpha, std::complex<double>* x) {
    for (int i = 0; i < n; ++i) x[i] *= alpha;
}
#endif

// ---------- copy: dst = src ----------

#ifdef HAVE_LAPACK
inline void copy(int n, const double* src, double* dst) {
    cblas_dcopy(n, src, 1, dst, 1);
}

inline void copy(int n, const std::complex<double>* src, std::complex<double>* dst) {
    cblas_zcopy(n, src, 1, dst, 1);
}
#else
inline void copy(int n, const double* src, double* dst) {
    for (int i = 0; i < n; ++i) dst[i] = src[i];
}

inline void copy(int n, const std::complex<double>* src, std::complex<double>* dst) {
    for (int i = 0; i < n; ++i) dst[i] = src[i];
}
#endif

// ---------- gemv: y = alpha * A * x + beta * y (row-major) ----------

#ifdef HAVE_LAPACK
inline void gemv(int m, int n, double alpha, const double* A, int lda,
                 const double* x, double beta, double* y) {
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, A, lda, x, 1, beta, y, 1);
}

inline void gemv(int m, int n, std::complex<double> alpha,
                 const std::complex<double>* A, int lda,
                 const std::complex<double>* x,
                 std::complex<double> beta, std::complex<double>* y) {
    cblas_zgemv(CblasRowMajor, CblasNoTrans, m, n, &alpha, A, lda, x, 1, &beta, y, 1);
}
#else
inline void gemv(int m, int n, double alpha, const double* A, int lda,
                 const double* x, double beta, double* y) {
    for (int i = 0; i < m; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) sum += A[i * lda + j] * x[j];
        y[i] = alpha * sum + beta * y[i];
    }
}

inline void gemv(int m, int n, std::complex<double> alpha,
                 const std::complex<double>* A, int lda,
                 const std::complex<double>* x,
                 std::complex<double> beta, std::complex<double>* y) {
    for (int i = 0; i < m; ++i) {
        std::complex<double> sum(0, 0);
        for (int j = 0; j < n; ++j) sum += A[i * lda + j] * x[j];
        y[i] = alpha * sum + beta * y[i];
    }
}
#endif

// ---------- gemm: C = alpha * A * B + beta * C (row-major) ----------

#ifdef HAVE_LAPACK
inline void gemm(CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                 int m, int n, int k,
                 std::complex<double> alpha,
                 const std::complex<double>* A, int lda,
                 const std::complex<double>* B, int ldb,
                 std::complex<double> beta,
                 std::complex<double>* C, int ldc) {
    cblas_zgemm(CblasRowMajor, transA, transB, m, n, k,
                &alpha, A, lda, B, ldb, &beta, C, ldc);
}
#endif

} // namespace blas

// ============================================================================
// Type traits helpers
// ============================================================================

namespace bicgstab_detail {

template<typename T> inline T zero_val();
template<> inline double zero_val<double>() { return 0.0; }
template<> inline std::complex<double> zero_val<std::complex<double>>() { return {0.0, 0.0}; }

template<typename T> inline T one_val();
template<> inline double one_val<double>() { return 1.0; }
template<> inline std::complex<double> one_val<std::complex<double>>() { return {1.0, 0.0}; }

inline double abs_val(double x) { return std::abs(x); }
inline double abs_val(std::complex<double> x) { return std::abs(x); }

} // namespace bicgstab_detail

// ============================================================================
// BiCGSTAB solver
// ============================================================================

namespace bicgstab {

struct Result {
    int iterations;     // Number of iterations performed
    double residual;    // Final relative residual ||r|| / ||b||
    bool converged;     // Whether tol was reached
};

/**
 * @brief Templated BiCGSTAB with preconditioner
 *
 * Solves A*x = b where A is defined implicitly by matvec.
 *
 * Algorithm: van der Vorst 1992, BiCGSTAB with right preconditioning.
 *
 * @tparam T Scalar type (double or std::complex<double>)
 * @tparam MatVecFunc Callable: void(const T* x, T* y) computing y = A*x
 * @tparam PrecFunc  Callable: void(const T* x, T* y) computing y = M^{-1}*x
 *
 * @param n       System dimension
 * @param matvec  Matrix-vector product function
 * @param precond Preconditioner function
 * @param rhs     Right-hand side (size n)
 * @param sol     Solution vector (size n, initial guess on input)
 * @param tol     Relative tolerance for convergence
 * @param max_iter Maximum iterations
 * @return Result struct with iteration count, residual, convergence flag
 */
template<typename T, typename MatVecFunc, typename PrecFunc>
Result Solve(int n,
             MatVecFunc matvec,
             PrecFunc precond,
             const T* rhs,
             T* sol,
             double tol,
             int max_iter) {
    using namespace bicgstab_detail;

    Result result = {0, 1.0, false};
    if (n <= 0) return result;

    std::vector<T> r(n), r0(n), p(n), v(n), s(n), t(n);
    std::vector<T> p_hat(n), s_hat(n);

    // r = b - A*x0
    matvec(sol, v.data());
    blas::copy(n, rhs, r.data());
    blas::axpy(n, T(-1), v.data(), r.data());
    blas::copy(n, r.data(), r0.data());

    T rho = one_val<T>(), alpha = one_val<T>(), omega = one_val<T>();
    std::fill(p.begin(), p.end(), zero_val<T>());
    std::fill(v.begin(), v.end(), zero_val<T>());

    double rhs_norm = blas::nrm2(n, rhs);
    if (rhs_norm < 1e-30) rhs_norm = 1.0;

    result.residual = blas::nrm2(n, r.data()) / rhs_norm;
    if (result.residual < tol) {
        result.converged = true;
        return result;
    }

    for (int iter = 1; iter <= max_iter; ++iter) {
        T rho_old = rho;
        rho = blas::dot(n, r0.data(), r.data());

        if (abs_val(rho) < 1e-30) {
            result.residual = blas::nrm2(n, r.data()) / rhs_norm;
            result.converged = result.residual < tol;
            result.iterations = iter;
            break;
        }

        if (iter == 1) {
            blas::copy(n, r.data(), p.data());
        } else {
            if (abs_val(rho_old * omega) < 1e-30) {
                result.residual = blas::nrm2(n, r.data()) / rhs_norm;
                result.converged = result.residual < tol;
                result.iterations = iter;
                break;
            }
            T beta = (rho / rho_old) * (alpha / omega);
            blas::axpy(n, -omega, v.data(), p.data());  // p -= omega * v
            blas::scal(n, beta, p.data());               // p *= beta
            blas::axpy(n, one_val<T>(), r.data(), p.data()); // p += r
        }

        // p_hat = M^{-1} * p
        precond(p.data(), p_hat.data());

        // v = A * p_hat
        matvec(p_hat.data(), v.data());

        T r0_dot_v = blas::dot(n, r0.data(), v.data());
        if (abs_val(r0_dot_v) < 1e-30) {
            result.residual = blas::nrm2(n, r.data()) / rhs_norm;
            result.converged = result.residual < tol;
            result.iterations = iter;
            break;
        }
        alpha = rho / r0_dot_v;

        // s = r - alpha * v
        blas::copy(n, r.data(), s.data());
        blas::axpy(n, -alpha, v.data(), s.data());

        double s_norm = blas::nrm2(n, s.data());
        if (s_norm / rhs_norm < tol) {
            blas::axpy(n, alpha, p_hat.data(), sol);
            result.residual = s_norm / rhs_norm;
            result.iterations = iter;
            result.converged = true;
            break;
        }

        // s_hat = M^{-1} * s
        precond(s.data(), s_hat.data());

        // t = A * s_hat
        matvec(s_hat.data(), t.data());

        T t_dot_s = blas::dot(n, t.data(), s.data());
        T t_dot_t = blas::dot(n, t.data(), t.data());
        if (abs_val(t_dot_t) < 1e-30) {
            blas::axpy(n, alpha, p_hat.data(), sol);
            result.residual = s_norm / rhs_norm;
            result.converged = result.residual < tol;
            result.iterations = iter;
            break;
        }
        omega = t_dot_s / t_dot_t;

        // x += alpha * p_hat + omega * s_hat
        blas::axpy(n, alpha, p_hat.data(), sol);
        blas::axpy(n, omega, s_hat.data(), sol);

        // r = s - omega * t
        blas::copy(n, s.data(), r.data());
        blas::axpy(n, -omega, t.data(), r.data());

        double r_norm = blas::nrm2(n, r.data());
        result.residual = r_norm / rhs_norm;
        result.iterations = iter;

        if (result.residual < tol) {
            result.converged = true;
            break;
        }

        // Detect numerical blowup
        if (std::isnan(result.residual) || std::isinf(result.residual)
            || result.residual > 1e15) {
            break;
        }

        if (abs_val(omega) < 1e-30) break;
    }

    return result;
}

/**
 * @brief Solve with identity preconditioner (no preconditioning)
 */
template<typename T, typename MatVecFunc>
Result SolveUnpreconditioned(int n,
                              MatVecFunc matvec,
                              const T* rhs,
                              T* sol,
                              double tol,
                              int max_iter) {
    auto identity = [](const T* x, T* y_out) {
        // This is never called with correct n because we capture nothing
        // Use the version below instead
    };
    // Use n-aware identity
    auto identity_precond = [n](const T* x, T* y_out) {
        blas::copy(n, x, y_out);
    };
    return Solve<T>(n, matvec, identity_precond, rhs, sol, tol, max_iter);
}

/**
 * @brief Solve dense system A*x = b with diagonal Jacobi preconditioner
 *
 * Convenience wrapper for dense matrices.
 *
 * @param n         System dimension
 * @param A         Dense matrix (n x n, row-major)
 * @param rhs       Right-hand side (size n)
 * @param sol       Solution (size n, initial guess on input)
 * @param tol       Convergence tolerance
 * @param max_iter  Maximum iterations
 * @return Result struct
 */
template<typename T>
Result DenseSolve(int n, const T* A, const T* rhs, T* sol,
                  double tol, int max_iter) {
    using namespace bicgstab_detail;

    // Dense matvec: y = A * x
    auto matvec = [A, n](const T* x, T* y) {
        blas::gemv(n, n, one_val<T>(), A, n, x, zero_val<T>(), y);
    };

    // Diagonal Jacobi preconditioner: y = diag(A)^{-1} * x
    std::vector<T> diag_inv(n);
    for (int i = 0; i < n; ++i) {
        T d = A[i * n + i];
        diag_inv[i] = (abs_val(d) > 1e-15) ? (one_val<T>() / d) : one_val<T>();
    }

    auto precond = [&diag_inv](const T* x, T* y) {
        int nn = static_cast<int>(diag_inv.size());
        for (int i = 0; i < nn; ++i) y[i] = diag_inv[i] * x[i];
    };

    return Solve<T>(n, matvec, precond, rhs, sol, tol, max_iter);
}

/**
 * @brief Invert dense matrix via column-by-column BiCGSTAB
 *
 * Solves A * X[:,j] = I[:,j] for each column j.
 * Useful for computing Y_branch = Z_eff^{-1} in PEEC MNA.
 *
 * @param n         Matrix dimension
 * @param A         Dense matrix (n x n, row-major, NOT modified)
 * @param X_inv     Output inverse (n x n, row-major)
 * @param tol       Convergence tolerance
 * @param max_iter  Maximum iterations per column
 * @return          Maximum residual across all columns
 */
template<typename T>
double DenseInvert(int n, const T* A, T* X_inv, double tol, int max_iter) {
    using namespace bicgstab_detail;

    // Dense matvec
    auto matvec = [A, n](const T* x, T* y) {
        blas::gemv(n, n, one_val<T>(), A, n, x, zero_val<T>(), y);
    };

    // Diagonal Jacobi
    std::vector<T> diag_inv(n);
    for (int i = 0; i < n; ++i) {
        T d = A[i * n + i];
        diag_inv[i] = (abs_val(d) > 1e-15) ? (one_val<T>() / d) : one_val<T>();
    }

    auto precond = [&diag_inv](const T* x, T* y) {
        int nn = static_cast<int>(diag_inv.size());
        for (int i = 0; i < nn; ++i) y[i] = diag_inv[i] * x[i];
    };

    double max_residual = 0.0;

    // Solve for each column of the identity
    std::vector<T> rhs(n, zero_val<T>());
    std::vector<T> sol(n, zero_val<T>());

    for (int j = 0; j < n; ++j) {
        // RHS = e_j (j-th unit vector)
        std::fill(rhs.begin(), rhs.end(), zero_val<T>());
        rhs[j] = one_val<T>();

        // Initial guess: diagonal inverse as starting point
        std::fill(sol.begin(), sol.end(), zero_val<T>());
        sol[j] = diag_inv[j];

        Result res = Solve<T>(n, matvec, precond, rhs.data(), sol.data(), tol, max_iter);

        // Store column j of the inverse (row-major: X_inv[i*n + j])
        for (int i = 0; i < n; ++i) {
            X_inv[i * n + j] = sol[i];
        }

        if (res.residual > max_residual) max_residual = res.residual;
    }

    return max_residual;
}

} // namespace bicgstab
} // namespace radia

#endif // RAD_BICGSTAB_H
