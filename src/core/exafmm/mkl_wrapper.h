#ifndef mkl_wrapper_h
#define mkl_wrapper_h

/**
 * @file mkl_wrapper.h
 * @brief Direct MKL BLAS/LAPACK calls for ExaFMM-t
 *
 * This file provides inline wrappers that directly call Intel MKL functions.
 * No fallback implementations - MKL is required.
 */

// Prevent Windows min/max macro conflicts with std::min/std::max
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#include <complex>
#include <algorithm>  // for std::min, std::max
#include "exafmm_t.h"
#include "timer.h"

// Intel MKL headers (MKL is required)
#include <mkl_blas.h>
#include <mkl_lapack.h>

using std::complex;

namespace exafmm_t {

  //==========================================================================
  // BLAS Level 2: Matrix-Vector Operations
  //==========================================================================

  //! Real matrix-vector multiply (row major): y = A * x
  inline void gemv(int m, int n, real_t* A, real_t* x, real_t* y) {
    char trans = 'T';  // Row major requires transpose
    real_t alpha = 1.0, beta = 0.0;
    MKL_INT incx = 1, incy = 1;
    MKL_INT mm = m, nn = n;
    dgemv(&trans, &nn, &mm, &alpha, A, &nn, x, &incx, &beta, y, &incy);
    add_flop((long long)(2*m*n));
  }

  //! Complex matrix-vector multiply (row major): y = A * x
  inline void gemv(int m, int n, complex_t* A, complex_t* x, complex_t* y) {
    char trans = 'T';
    MKL_Complex16 alpha = {1., 0.}, beta = {0., 0.};
    MKL_INT incx = 1, incy = 1;
    MKL_INT mm = m, nn = n;
    zgemv(&trans, &nn, &mm, &alpha, reinterpret_cast<MKL_Complex16*>(A), &nn,
          reinterpret_cast<MKL_Complex16*>(x), &incx, &beta,
          reinterpret_cast<MKL_Complex16*>(y), &incy);
  }

  //==========================================================================
  // BLAS Level 3: Matrix-Matrix Operations
  //==========================================================================

  //! Real matrix-matrix multiply (row major): C = A * B
  inline void gemm(int m, int n, int k, real_t* A, real_t* B, real_t* C) {
    char transA = 'N', transB = 'N';
    real_t alpha = 1.0, beta = 0.0;
    MKL_INT mm = m, nn = n, kk = k;
    dgemm(&transA, &transB, &nn, &mm, &kk, &alpha, B, &nn, A, &kk, &beta, C, &nn);
  }

  //! Complex matrix-matrix multiply (row major): C = A * B
  inline void gemm(int m, int n, int k, complex_t* A, complex_t* B, complex_t* C) {
    char transA = 'N', transB = 'N';
    MKL_Complex16 alpha = {1., 0.}, beta = {0., 0.};
    MKL_INT mm = m, nn = n, kk = k;
    zgemm(&transA, &transB, &nn, &mm, &kk, &alpha,
          reinterpret_cast<MKL_Complex16*>(B), &nn,
          reinterpret_cast<MKL_Complex16*>(A), &kk, &beta,
          reinterpret_cast<MKL_Complex16*>(C), &nn);
  }

  //==========================================================================
  // LAPACK: SVD
  //==========================================================================

  //! Real SVD (row major): A = U * S * VT
  inline void svd(int m, int n, real_t* A, real_t* S, real_t* U, real_t* VT) {
    char JOBU = 'S', JOBVT = 'S';
    MKL_INT INFO;
    MKL_INT LWORK = std::max(3*std::min(m,n)+std::max(m,n), 5*std::min(m,n));
    LWORK = std::max(LWORK, (MKL_INT)1);
    int kk = std::min(m, n);
    RealVec tS(kk, 0.);
    RealVec WORK(LWORK);
    MKL_INT mm = m, nn = n;
    dgesvd(&JOBU, &JOBVT, &nn, &mm, A, &nn, &tS[0], VT, &nn, U, (MKL_INT*)&kk, &WORK[0], &LWORK, &INFO);
    // Copy singular values from 1D layout (tS) to 2D layout (S)
    for(int i=0; i<kk; i++) {
      S[i*n+i] = tS[i];
    }
  }

  //! Complex SVD (row major): A = U * S * VH
  inline void svd(int m, int n, complex_t* A, real_t* S, complex_t* U, complex_t* VH) {
    char JOBU = 'S', JOBVT = 'S';
    MKL_INT INFO;
    MKL_INT LWORK = std::max(3*std::min(m,n)+std::max(m,n), 5*std::min(m,n));
    LWORK = std::max(LWORK, (MKL_INT)1);
    int kk = std::min(m, n);
    RealVec tS(kk, 0.);
    ComplexVec WORK(LWORK);
    RealVec RWORK(5*kk);
    MKL_INT mm = m, nn = n;
    zgesvd(&JOBU, &JOBVT, &nn, &mm,
           reinterpret_cast<MKL_Complex16*>(A), &nn, &tS[0],
           reinterpret_cast<MKL_Complex16*>(VH), &nn,
           reinterpret_cast<MKL_Complex16*>(U), (MKL_INT*)&kk,
           reinterpret_cast<MKL_Complex16*>(&WORK[0]), &LWORK, &RWORK[0], &INFO);
    // Copy singular values from 1D layout (tS) to 2D layout (S)
    for(int i=0; i<kk; i++) {
      S[i*n+i] = tS[i];
    }
  }

  //==========================================================================
  // Matrix Utility Functions
  //==========================================================================

  //! Transpose real matrix
  inline RealVec transpose(RealVec& vec, int m, int n) {
    RealVec temp(vec.size());
    for(int i=0; i<m; i++) {
      for(int j=0; j<n; j++) {
        temp[j*m+i] = vec[i*n+j];
      }
    }
    return temp;
  }

  //! Transpose complex matrix
  inline ComplexVec transpose(ComplexVec& vec, int m, int n) {
    ComplexVec temp(vec.size());
    for(int i=0; i<m; i++) {
      for(int j=0; j<n; j++) {
        temp[j*m+i] = vec[i*n+j];
      }
    }
    return temp;
  }

  //! Conjugate transpose complex matrix
  inline ComplexVec conjugate_transpose(ComplexVec& vec, int m, int n) {
    ComplexVec temp(vec.size());
    for(int i=0; i<m; i++) {
      for(int j=0; j<n; j++) {
        temp[j*m+i] = std::conj(vec[i*n+j]);
      }
    }
    return temp;
  }

}  // end namespace exafmm_t
#endif
