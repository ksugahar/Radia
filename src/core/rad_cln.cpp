/**
 * @file rad_cln.cpp
 * @brief CLN (Cauer Ladder Network) model order reduction implementation
 *
 * Uses Intel MKL for linear algebra operations.
 */

#include "rad_cln.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <stdexcept>

#include "rad_parallel.h"

// Intel MKL
#include <mkl.h>

// HACApK for ACA+ (optional, for aca_compress)
#ifdef RADIA_HACAPK_ENABLED
#include "../ext/HACApK/cHACApK_cpp.h"
#endif

namespace radia {
namespace cln {



std::vector<double> build_tridiagonal(const double* diag, int n) {
    std::vector<double> T(n * n, 0.0);

    // Build U matrix: U[i,i] = 1, U[i,i+1] = -1
    // T = U^T * diag(diag) * U
    // Which gives:
    //   T[i,i] = diag[i] + diag[i+1] (for i < n-1)
    //   T[n-1,n-1] = diag[n-1]
    //   T[i,i+1] = T[i+1,i] = -diag[i+1]

    for (int i = 0; i < n; ++i) {
        if (i < n - 1) {
            T[i * n + i] = diag[i] + diag[i + 1];
            T[i * n + (i + 1)] = -diag[i + 1];
            T[(i + 1) * n + i] = -diag[i + 1];
        } else {
            T[i * n + i] = diag[i];
        }
    }

    return T;
}

std::complex<double> compute_cln_impedance(
    const double* R_diag,
    const double* L_tridiag,
    int n,
    double freq
) {
    const double PI = 3.14159265358979323846;
    std::complex<double> s(0.0, 2.0 * PI * freq);

    // Build Z = R_diag + s * L_tridiag (complex matrix)
    std::vector<std::complex<double>> Z(n * n);
    for (int i = 0; i < n * n; ++i) {
        Z[i] = R_diag[i] + s * L_tridiag[i];
    }

    // Voltage source: V = [1, 0, 0, ...]
    std::vector<std::complex<double>> V(n, 0.0);
    V[0] = 1.0;

    // Solve Z * I = V using LAPACK zgesv
    std::vector<int> ipiv(n);

    // Note: LAPACK expects column-major, but our matrices are row-major
    // For symmetric/Hermitian matrices, row-major = column-major
    // So we can use zgesv directly

    // Convert to MKL complex type
    std::vector<MKL_Complex16> Z_mkl(n * n);
    std::vector<MKL_Complex16> V_mkl(n);
    for (int i = 0; i < n * n; ++i) {
        Z_mkl[i].real = Z[i].real();
        Z_mkl[i].imag = Z[i].imag();
    }
    for (int i = 0; i < n; ++i) {
        V_mkl[i].real = V[i].real();
        V_mkl[i].imag = V[i].imag();
    }

    int info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, n, 1, Z_mkl.data(), n, ipiv.data(), V_mkl.data(), 1);
    if (info != 0) {
        return std::complex<double>(0.0, 0.0);  // Error case
    }

    // Impedance = 1 / I[0]
    std::complex<double> I0(V_mkl[0].real, V_mkl[0].imag);
    return 1.0 / I0;
}

void compute_cln_impedance_sweep(
    const double* R_diag,
    const double* L_tridiag,
    int n,
    const double* freqs,
    int n_freqs,
    std::complex<double>* Z_out
) {
    // TaskManager parallelization for frequency sweep
    ngcore::ParallelFor(ngcore::IntRange(n_freqs), [&](size_t f) {
        Z_out[f] = compute_cln_impedance(R_diag, L_tridiag, n, freqs[f]);
    });
}

std::vector<double> transform_coupling(
    const double* Q,
    const double* M_LS,
    int n_loop,
    int n_reduced,
    int n_star
) {
    // M_LS_reduced = Q^T * M_LS
    // Q: (n_loop x n_reduced), M_LS: (n_loop x n_star)
    // Result: (n_reduced x n_star)

    std::vector<double> M_LS_reduced(n_reduced * n_star);

    // Use BLAS dgemm: C = alpha * A^T * B + beta * C
    // A = Q (n_loop x n_reduced), B = M_LS (n_loop x n_star)
    // C = M_LS_reduced (n_reduced x n_star)
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                n_reduced, n_star, n_loop,
                1.0, Q, n_reduced, M_LS, n_star,
                0.0, M_LS_reduced.data(), n_star);

    return M_LS_reduced;
}

std::vector<double> transform_port_vector(
    const double* Q,
    const double* v,
    int n_loop,
    int n_reduced
) {
    // v_reduced = Q^T * v
    // Q: (n_loop x n_reduced), v: (n_loop)
    // Result: (n_reduced)

    std::vector<double> v_reduced(n_reduced);

    // Use BLAS dgemv: y = alpha * A^T * x + beta * y
    cblas_dgemv(CblasRowMajor, CblasTrans,
                n_loop, n_reduced,
                1.0, Q, n_reduced, v, 1,
                0.0, v_reduced.data(), 1);

    return v_reduced;
}

std::complex<double> compute_loop_star_impedance(
    const LoopStarCLN& cln,
    double freq,
    const double* port_loop
) {
    const double PI = 3.14159265358979323846;
    std::complex<double> s(0.0, 2.0 * PI * freq);
    std::complex<double> s_inv = 1.0 / s;

    int n_L = cln.n_reduced;
    int n_S = cln.n_star;
    int n_total = n_L + n_S;

    // Build full system matrix Z
    // Z = [Z_LL  Z_LS]
    //     [Z_SL  Z_SS]
    std::vector<MKL_Complex16> Z(n_total * n_total);

    // Z_LL = R_diag + s * L_tridiag  (n_L x n_L)
    for (int i = 0; i < n_L; ++i) {
        for (int j = 0; j < n_L; ++j) {
            int idx_full = i * n_total + j;
            int idx_L = i * n_L + j;
            std::complex<double> val = cln.R_diag[idx_L] + s * cln.L_tridiag[idx_L];
            Z[idx_full].real = val.real();
            Z[idx_full].imag = val.imag();
        }
    }

    // Z_LS = s * M_LS_reduced  (n_L x n_S)
    for (int i = 0; i < n_L; ++i) {
        for (int j = 0; j < n_S; ++j) {
            int idx_full = i * n_total + (n_L + j);
            int idx_LS = i * n_S + j;
            std::complex<double> val = s * cln.M_LS_reduced[idx_LS];
            Z[idx_full].real = val.real();
            Z[idx_full].imag = val.imag();
        }
    }

    // Z_SL = Z_LS^T = s * M_LS_reduced^T  (n_S x n_L)
    for (int i = 0; i < n_S; ++i) {
        for (int j = 0; j < n_L; ++j) {
            int idx_full = (n_L + i) * n_total + j;
            int idx_LS = j * n_S + i;  // Transpose
            std::complex<double> val = s * cln.M_LS_reduced[idx_LS];
            Z[idx_full].real = val.real();
            Z[idx_full].imag = val.imag();
        }
    }

    // Z_SS = P / s  (n_S x n_S)
    for (int i = 0; i < n_S; ++i) {
        for (int j = 0; j < n_S; ++j) {
            int idx_full = (n_L + i) * n_total + (n_L + j);
            int idx_P = i * n_S + j;
            std::complex<double> val = cln.P[idx_P] * s_inv;
            Z[idx_full].real = val.real();
            Z[idx_full].imag = val.imag();
        }
    }

    // Transform port vector
    auto v_reduced = transform_port_vector(cln.Q.data(), port_loop, cln.n_loop, n_L);

    // Build voltage vector V = [v_reduced; 0]
    std::vector<MKL_Complex16> V(n_total);
    for (int i = 0; i < n_L; ++i) {
        V[i].real = v_reduced[i];
        V[i].imag = 0.0;
    }
    for (int i = 0; i < n_S; ++i) {
        V[n_L + i].real = 0.0;
        V[n_L + i].imag = 0.0;
    }

    // Solve Z * I = V
    std::vector<int> ipiv(n_total);
    int info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, n_total, 1, Z.data(), n_total, ipiv.data(), V.data(), 1);

    if (info != 0) {
        return std::complex<double>(0.0, 0.0);  // Error
    }

    // Impedance = 1 / (v_reduced^T * I_L)
    // where I_L = V[0:n_L] (after solve, V contains I)
    std::complex<double> dot_product(0.0, 0.0);
    for (int i = 0; i < n_L; ++i) {
        std::complex<double> I_i(V[i].real, V[i].imag);
        dot_product += v_reduced[i] * I_i;
    }

    return 1.0 / dot_product;
}

// ============================================================
// ACA+ Low-Rank Approximation
// ============================================================

std::vector<double> ACAResult::reconstruct() const {
    // P_approx = U * V^T
    std::vector<double> P_approx(n * n, 0.0);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
                n, n, k,
                1.0, U.data(), k, V.data(), k,
                0.0, P_approx.data(), n);
    return P_approx;
}

// Global variable for ACA+ callback (thread-local for safety)
namespace {
    thread_local const double* g_aca_matrix = nullptr;
    thread_local int g_aca_matrix_n = 0;
}

#ifdef RADIA_HACAPK_ENABLED
// Callback function for HACApK ACA+
extern "C" double aca_entry_func(int i, int j, int /* i_bemv */) {
    // HACApK uses 1-based indices
    int i0 = i - 1;
    int j0 = j - 1;
    if (i0 >= 0 && i0 < g_aca_matrix_n && j0 >= 0 && j0 < g_aca_matrix_n) {
        return g_aca_matrix[i0 * g_aca_matrix_n + j0];
    }
    return 0.0;
}
#endif

ACAResult aca_compress(
    const double* P,
    int n,
    double eps,
    int kmax
) {
    ACAResult result;
    result.n = n;
    result.converged = false;

    if (n <= 0) {
        result.k = 0;
        result.compression_ratio = 0.0;
        return result;
    }

    if (kmax < 0) {
        kmax = n;
    }

#ifdef RADIA_HACAPK_ENABLED
    // Use HACApK ACA+ algorithm
    g_aca_matrix = P;
    g_aca_matrix_n = n;

    // Allocate output matrices (maximum rank = kmax)
    result.U.resize(n * kmax, 0.0);
    result.V.resize(n * kmax, 0.0);

    // Create identity permutation (no reordering)
    std::vector<int> lod(n);
    for (int i = 0; i < n; ++i) {
        lod[i] = i + 1;  // 1-based
    }

    // ACA+ parameters
    std::vector<double> param(100, 0.0);
    param[71] = eps;   // ACA tolerance
    param[72] = 1.0;   // ACA_EPS factor

    // Compute Frobenius norm for scaling
    double znrmmat = 0.0;
    for (int i = 0; i < n * n; ++i) {
        znrmmat += P[i] * P[i];
    }
    znrmmat = std::sqrt(znrmmat);
    if (znrmmat < 1e-30) {
        znrmmat = 1.0;
    }

    // Call ACA+ algorithm
    int actual_rank = HACApK_acaplus_cpp(
        result.U.data(),  // U: n x kmax
        result.V.data(),  // V: n x kmax
        param.data(),
        n,               // ndl (rows)
        n,               // ndt (cols)
        1,               // nstrtl (1-based start)
        1,               // nstrtt (1-based start)
        lod.data(),
        0,               // i_bemv (dummy)
        kmax,            // maximum rank
        eps,             // tolerance
        znrmmat,         // matrix norm
        eps,             // pACA_EPS
        aca_entry_func
    );

    // Clean up
    g_aca_matrix = nullptr;
    g_aca_matrix_n = 0;

    if (actual_rank > 0) {
        result.k = actual_rank;
        result.converged = true;

        // Resize to actual rank
        std::vector<double> U_compact(n * actual_rank);
        std::vector<double> V_compact(n * actual_rank);

        // Copy columns (HACApK stores in column-major for U, V)
        // U is stored as n x kmax (row-major), take first actual_rank columns
        for (int i = 0; i < n; ++i) {
            for (int k_idx = 0; k_idx < actual_rank; ++k_idx) {
                U_compact[i * actual_rank + k_idx] = result.U[i * kmax + k_idx];
                V_compact[i * actual_rank + k_idx] = result.V[i * kmax + k_idx];
            }
        }

        result.U = std::move(U_compact);
        result.V = std::move(V_compact);
    } else {
        // ACA+ failed or matrix is zero
        result.k = 0;
    }

#else
    // Fallback: Use SVD for low-rank approximation (if HACApK not available)
    // This is slower but always works

    // Copy P for SVD (will be overwritten)
    std::vector<double> A(P, P + n * n);

    // SVD: A = U * S * V^T
    std::vector<double> S(n);
    std::vector<double> U(n * n);
    std::vector<double> Vt(n * n);
    std::vector<double> superb(n);

    int info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A',
                               n, n, A.data(), n,
                               S.data(), U.data(), n, Vt.data(), n,
                               superb.data());

    if (info == 0) {
        // Determine rank based on tolerance
        double S_max = S[0];
        if (S_max < 1e-30) {
            S_max = 1.0;
        }

        int actual_rank = 0;
        for (int i = 0; i < n && i < kmax; ++i) {
            if (S[i] / S_max > eps) {
                actual_rank = i + 1;
            }
        }

        if (actual_rank > 0) {
            result.k = actual_rank;
            result.converged = true;

            // U_k = U[:, 0:k] * diag(sqrt(S[0:k]))
            // V_k = V[:, 0:k] * diag(sqrt(S[0:k]))
            result.U.resize(n * actual_rank);
            result.V.resize(n * actual_rank);

            for (int i = 0; i < n; ++i) {
                for (int k_idx = 0; k_idx < actual_rank; ++k_idx) {
                    double sqrt_s = std::sqrt(S[k_idx]);
                    result.U[i * actual_rank + k_idx] = U[i * n + k_idx] * sqrt_s;
                    result.V[i * actual_rank + k_idx] = Vt[k_idx * n + i] * sqrt_s;  // Transpose
                }
            }
        } else {
            result.k = 0;
        }
    } else {
        // SVD failed
        result.k = 0;
    }
#endif

    result.compression_ratio = (n > 0) ? static_cast<double>(result.k) / n : 0.0;
    return result;
}

std::complex<double> compute_compressed_impedance(
    const LoopStarCLNCompressed& cln,
    double freq,
    const double* port_loop
) {
    const double PI = 3.14159265358979323846;
    std::complex<double> s(0.0, 2.0 * PI * freq);
    std::complex<double> s_inv = 1.0 / s;

    int n_L = cln.n_reduced;
    int n_S = cln.k_star;  // Compressed star dimension
    int n_total = n_L + n_S;

    // Build full system matrix Z in mode space
    // Z = [Z_LL       Z_LS_mode  ]
    //     [Z_SL_mode  Z_SS_mode  ]
    std::vector<MKL_Complex16> Z(n_total * n_total);

    // Z_LL = R_diag + s * L_tridiag  (n_L x n_L)
    for (int i = 0; i < n_L; ++i) {
        for (int j = 0; j < n_L; ++j) {
            int idx_full = i * n_total + j;
            int idx_L = i * n_L + j;
            std::complex<double> val = cln.R_diag[idx_L] + s * cln.L_tridiag[idx_L];
            Z[idx_full].real = val.real();
            Z[idx_full].imag = val.imag();
        }
    }

    // Z_LS_mode = s * M_LS_mode  (n_L x n_S)
    for (int i = 0; i < n_L; ++i) {
        for (int j = 0; j < n_S; ++j) {
            int idx_full = i * n_total + (n_L + j);
            int idx_LS = i * n_S + j;
            std::complex<double> val = s * cln.M_LS_mode[idx_LS];
            Z[idx_full].real = val.real();
            Z[idx_full].imag = val.imag();
        }
    }

    // Z_SL_mode = Z_LS_mode^T  (n_S x n_L)
    for (int i = 0; i < n_S; ++i) {
        for (int j = 0; j < n_L; ++j) {
            int idx_full = (n_L + i) * n_total + j;
            int idx_LS = j * n_S + i;  // Transpose
            std::complex<double> val = s * cln.M_LS_mode[idx_LS];
            Z[idx_full].real = val.real();
            Z[idx_full].imag = val.imag();
        }
    }

    // Z_SS_mode = P_mode / s  (n_S x n_S)
    for (int i = 0; i < n_S; ++i) {
        for (int j = 0; j < n_S; ++j) {
            int idx_full = (n_L + i) * n_total + (n_L + j);
            int idx_P = i * n_S + j;
            std::complex<double> val = cln.P_mode[idx_P] * s_inv;
            Z[idx_full].real = val.real();
            Z[idx_full].imag = val.imag();
        }
    }

    // Transform port vector: v_reduced = Q^T * port_loop
    auto v_reduced = transform_port_vector(cln.Q.data(), port_loop, cln.n_loop, n_L);

    // Build voltage vector V = [v_reduced; 0]
    std::vector<MKL_Complex16> V(n_total);
    for (int i = 0; i < n_L; ++i) {
        V[i].real = v_reduced[i];
        V[i].imag = 0.0;
    }
    for (int i = 0; i < n_S; ++i) {
        V[n_L + i].real = 0.0;
        V[n_L + i].imag = 0.0;
    }

    // Solve Z * I = V
    std::vector<int> ipiv(n_total);
    int info = LAPACKE_zgesv(LAPACK_ROW_MAJOR, n_total, 1, Z.data(), n_total, ipiv.data(), V.data(), 1);

    if (info != 0) {
        return std::complex<double>(0.0, 0.0);  // Error
    }

    // Impedance = 1 / (v_reduced^T * I_L)
    std::complex<double> dot_product(0.0, 0.0);
    for (int i = 0; i < n_L; ++i) {
        std::complex<double> I_i(V[i].real, V[i].imag);
        dot_product += v_reduced[i] * I_i;
    }

    return 1.0 / dot_product;
}

}  // namespace cln
}  // namespace radia
