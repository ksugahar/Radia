/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/// @file dense_matrix.hpp
/// @brief Small dense matrix with LU inverse

#ifndef SPARSESOLV_CORE_DENSE_MATRIX_HPP
#define SPARSESOLV_CORE_DENSE_MATRIX_HPP

#include "types.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef SPARSESOLV_USE_MKL
#include <mkl_cblas.h>
#include <mkl_lapacke.h>
#endif

namespace sparsesolv {

namespace detail {
    /// Absolute value that works for both real and complex
    template<typename T>
    inline double dense_abs(T val) { return std::abs(val); }

#ifdef SPARSESOLV_USE_MKL
    // ---- GEMM wrappers (row-major) ----
    inline void mkl_gemm(index_t m, index_t n, index_t k,
                          double alpha, const double* A, index_t lda,
                          const double* B, index_t ldb,
                          double beta, double* C, index_t ldc) {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    inline void mkl_gemm(index_t m, index_t n, index_t k,
                          complex_t alpha, const complex_t* A, index_t lda,
                          const complex_t* B, index_t ldb,
                          complex_t beta, complex_t* C, index_t ldc) {
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
    }

    // ---- LU factorization + inverse wrappers ----
    inline lapack_int mkl_getrf(index_t n, double* data, index_t lda,
                                 lapack_int* ipiv) {
        return LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, data, lda, ipiv);
    }

    inline lapack_int mkl_getrf(index_t n, complex_t* data, index_t lda,
                                 lapack_int* ipiv) {
        return LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n,
                              reinterpret_cast<lapack_complex_double*>(data),
                              lda, ipiv);
    }

    inline lapack_int mkl_getri(index_t n, double* data, index_t lda,
                                 const lapack_int* ipiv) {
        return LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, data, lda, ipiv);
    }

    inline lapack_int mkl_getri(index_t n, complex_t* data, index_t lda,
                                 const lapack_int* ipiv) {
        return LAPACKE_zgetri(LAPACK_ROW_MAJOR, n,
                              reinterpret_cast<lapack_complex_double*>(data),
                              lda, ipiv);
    }
#endif
}

/// Small dense matrix (row-major) with LU factorization
template<typename Scalar>
class DenseMatrix {
public:
    DenseMatrix() : rows_(0), cols_(0) {}

    DenseMatrix(index_t rows, index_t cols)
        : rows_(rows), cols_(cols), data_(rows * cols, Scalar(0)) {}

    index_t rows() const { return rows_; }
    index_t cols() const { return cols_; }

    Scalar& operator()(index_t i, index_t j) {
        return data_[i * cols_ + j];
    }

    const Scalar& operator()(index_t i, index_t j) const {
        return data_[i * cols_ + j];
    }

    /// Scale row i by factor
    void scale_row(index_t i, Scalar factor) {
        for (index_t j = 0; j < cols_; ++j)
            data_[i * cols_ + j] *= factor;
    }

    /// Scale column j by factor
    void scale_col(index_t j, Scalar factor) {
        for (index_t i = 0; i < rows_; ++i)
            data_[i * cols_ + j] *= factor;
    }

    /// Negate all entries
    void negate() {
        for (auto& v : data_) v = -v;
    }

    /// Matrix multiplication: C = A * B
    static DenseMatrix multiply(const DenseMatrix& A, const DenseMatrix& B) {
        DenseMatrix C(A.rows_, B.cols_);
#ifdef SPARSESOLV_USE_MKL
        detail::mkl_gemm(A.rows_, B.cols_, A.cols_,
                          Scalar(1), A.data_.data(), A.cols_,
                          B.data_.data(), B.cols_,
                          Scalar(0), C.data_.data(), C.cols_);
#else
        for (index_t i = 0; i < A.rows_; ++i) {
            for (index_t p = 0; p < A.cols_; ++p) {
                Scalar a_ip = A(i, p);
                for (index_t j = 0; j < B.cols_; ++j) {
                    C(i, j) += a_ip * B(p, j);
                }
            }
        }
#endif
        return C;
    }

    /// Multiply-add: C += A * B
    static void multiply_add(const DenseMatrix& A, const DenseMatrix& B,
                              DenseMatrix& C) {
#ifdef SPARSESOLV_USE_MKL
        detail::mkl_gemm(A.rows_, B.cols_, A.cols_,
                          Scalar(1), A.data_.data(), A.cols_,
                          B.data_.data(), B.cols_,
                          Scalar(1), C.data_.data(), C.cols_);
#else
        for (index_t i = 0; i < A.rows_; ++i) {
            for (index_t p = 0; p < A.cols_; ++p) {
                Scalar a_ip = A(i, p);
                for (index_t j = 0; j < B.cols_; ++j) {
                    C(i, j) += a_ip * B(p, j);
                }
            }
        }
#endif
    }

    /// Matrix-vector multiply: y = A * x
    void matvec(const Scalar* x, Scalar* y) const {
        for (index_t i = 0; i < rows_; ++i) {
            Scalar sum = Scalar(0);
            for (index_t j = 0; j < cols_; ++j) {
                sum += data_[i * cols_ + j] * x[j];
            }
            y[i] = sum;
        }
    }

    /// Compute inverse in-place via LU decomposition
    void invert() {
        if (rows_ != cols_)
            throw std::runtime_error("DenseMatrix::invert: matrix must be square ("
                + std::to_string(rows_) + " x " + std::to_string(cols_) + ")");
        const index_t n = rows_;
        if (n == 0) return;

#ifdef SPARSESOLV_USE_MKL
        std::vector<lapack_int> ipiv(n);
        lapack_int info = detail::mkl_getrf(n, data_.data(), n, ipiv.data());
        if (info != 0)
            throw std::runtime_error("DenseMatrix::invert: dgetrf failed (info=" +
                                     std::to_string(info) + ")");
        info = detail::mkl_getri(n, data_.data(), n, ipiv.data());
        if (info != 0)
            throw std::runtime_error("DenseMatrix::invert: dgetri failed (info=" +
                                     std::to_string(info) + ")");
#else
        // LU factorize
        auto piv = lu_factorize();

        // Build permutation matrix P where PA = LU, i.e. P[k, piv[k]] = 1
        DenseMatrix inv(n, n);
        for (index_t k = 0; k < n; ++k) {
            inv(k, piv[k]) = Scalar(1);
        }

        // Forward substitution: solve L * Y = P * I
        for (index_t j = 0; j < n; ++j) {
            for (index_t i = 1; i < n; ++i) {
                for (index_t k = 0; k < i; ++k) {
                    inv(i, j) -= (*this)(i, k) * inv(k, j);
                }
            }
        }

        // Backward substitution: solve U * X = Y
        for (index_t j = 0; j < n; ++j) {
            for (index_t i = n; i-- > 0;) {
                for (index_t k = i + 1; k < n; ++k) {
                    inv(i, j) -= (*this)(i, k) * inv(k, j);
                }
                inv(i, j) /= (*this)(i, i);
            }
        }

        // Copy result back
        data_ = std::move(inv.data_);
#endif
    }

private:
    /// In-place LU with partial pivoting. Returns pivot permutation.
    std::vector<index_t> lu_factorize() {
        const index_t n = rows_;
        std::vector<index_t> piv(n);
        for (index_t i = 0; i < n; ++i) piv[i] = i;

        for (index_t k = 0; k < n; ++k) {
            // Find pivot: max |A[i,k]| for i >= k
            double max_val = 0;
            index_t max_row = k;
            for (index_t i = k; i < n; ++i) {
                double val = detail::dense_abs((*this)(i, k));
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }

            if (max_val < 1e-15) {
                throw std::runtime_error(
                    "DenseMatrix::lu_factorize: singular or near-singular matrix");
            }

            // Swap rows k and max_row
            if (max_row != k) {
                std::swap(piv[k], piv[max_row]);
                for (index_t j = 0; j < n; ++j) {
                    std::swap((*this)(k, j), (*this)(max_row, j));
                }
            }

            // Eliminate below
            Scalar diag = (*this)(k, k);
            for (index_t i = k + 1; i < n; ++i) {
                Scalar factor = (*this)(i, k) / diag;
                (*this)(i, k) = factor;  // Store L factor
                for (index_t j = k + 1; j < n; ++j) {
                    (*this)(i, j) -= factor * (*this)(k, j);
                }
            }
        }
        return piv;
    }

    index_t rows_, cols_;
    std::vector<Scalar> data_;
};

} // namespace sparsesolv

#endif // SPARSESOLV_CORE_DENSE_MATRIX_HPP
