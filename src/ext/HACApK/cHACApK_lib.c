/*
!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 1.3.0                                             *
!                                                                     *
!   License                                                           *
!     This file is part of HACApK.                                    *
!     HACApK is a free software, you can use it under the terms       *
!     of The MIT License (MIT). See LICENSE file and User's guide     *
!     for more details.                                               *
!                                                                     *
!   ppOpen-HPC project:                                               *
!     Open Source Infrastructure for Development and Execution of     *
!     Large-Scale Scientific Applications on Post-Peta-Scale          *
!     Supercomputers with Automatic Tuning (AT).                      *
!                                                                     *
!   Sponsorship:                                                      *
!     Japan Science and Technology Agency (JST), Basic Research       *
!     Programs: CREST, Development of System Software Technologies    *
!     for post-Peta Scale High Performance Computing.                 *
!                                                                     *
!   Copyright (c) 2015 <Akihiro Ida and Takeshi Iwashita>             *
!                                                                     *
!=====================================================================*
*/
#include "cHACApK_lib.h"
#include <math.h>

/* BLAS dgemv for optimized matrix-vector operations */
#ifdef HAVE_LAPACK
/* Intel MKL CBLAS interface */
#ifdef __cplusplus
extern "C" {
#endif

/* CBLAS constants */
typedef enum { CblasRowMajor=101, CblasColMajor=102 } CBLAS_ORDER;
typedef enum { CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113 } CBLAS_TRANSPOSE;

/* CBLAS dgemv: y = alpha*A*x + beta*y */
extern void cblas_dgemv(CBLAS_ORDER order, CBLAS_TRANSPOSE TransA,
                        int M, int N, double alpha, const double *A, int lda,
                        const double *X, int incX, double beta, double *Y, int incY);

/* CBLAS dcopy: y = x (with strides) */
extern void cblas_dcopy(int N, const double *X, int incX, double *Y, int incY);

/* CBLAS dnrm2: Euclidean norm */
extern double cblas_dnrm2(int N, const double *X, int incX);

#ifdef __cplusplus
}
#endif
#endif /* HAVE_LAPACK */

//***cHACApK_unrm_d
// Compute the Euclidean norm of a vector
// Optimized with BLAS dnrm2 when available
double cHACApK_unrm_d(int n, double *vec) {
#ifdef HAVE_LAPACK
  return cblas_dnrm2(n, vec, 1);
#else
  double norm = 0.0;
  int i;
  for (i = 0; i < n; i++) {
    norm += vec[i] * vec[i];
  }
  return sqrt(norm);
#endif
}

//***cHACApK_extract_col
// Extract column from column-major matrix with stride access
// Optimized with BLAS dcopy when available
// src: column-major matrix (m x n), dst: output vector (k elements)
// Extracts: dst[i] = src[row + i*ldsrc] for i=0..k-1
void cHACApK_extract_col(double *dst, const double *src, int row, int ldsrc, int k) {
#ifdef HAVE_LAPACK
  /* Use BLAS dcopy with stride for optimized strided access */
  if (k > 0) {
    cblas_dcopy(k, src + row, ldsrc, dst, 1);
  }
#else
  int i;
  for (i = 0; i < k; i++) {
    dst[i] = src[row + i * ldsrc];
  }
#endif
}

//***cHACApK_maxabsvalloc_d
// Find maximum absolute value and its location
void cHACApK_maxabsvalloc_d(double *vec, double *maxval, int *loc, int n) {
  int i;
  *maxval = 0.0;
  *loc = 0;
  for (i = 0; i < n; i++) {
    double absval = fabs(vec[i]);
    if (absval > *maxval) {
      *maxval = absval;
      *loc = i;
    }
  }
}

//***cHACApK_minabsvalloc_d
// Find minimum absolute value and its location
void cHACApK_minabsvalloc_d(double *vec, double *minval, int *loc, int n) {
  int i;
  *minval = fabs(vec[0]);
  *loc = 0;
  for (i = 1; i < n; i++) {
    double absval = fabs(vec[i]);
    if (absval < *minval) {
      *minval = absval;
      *loc = i;
    }
  }
}

//***cHACApK_adotsub_dsm
// Subtract the dot product of (zaa column k) and zz from vec
// vec(1:ndl) = vec(1:ndl) - zaa(1:ndl, 1:k) * zz(1:k)
// Optimized with BLAS dgemv when available
void cHACApK_adotsub_dsm(double *vec, double *zaa, double *zz, int ndl, int k, int ldaa) {
#ifdef HAVE_LAPACK
  /* Use BLAS dgemv for optimized matrix-vector multiplication:
   * vec = 1.0 * vec + (-1.0) * zaa * zz
   * = vec - zaa * zz
   *
   * cblas_dgemv(order, trans, M, N, alpha, A, lda, x, incx, beta, y, incy)
   * Computes: y = alpha * A * x + beta * y
   *
   * Here: y = vec (M x 1), A = zaa (M x N = ndl x k), x = zz (N x 1)
   *       alpha = -1.0, beta = 1.0
   *       => vec = -1.0 * zaa * zz + 1.0 * vec = vec - zaa * zz
   */
  if (k > 0) {
    cblas_dgemv(CblasColMajor, CblasNoTrans,
                ndl, k, -1.0, zaa, ldaa, zz, 1, 1.0, vec, 1);
  }
#else
  /* Fallback: simple loop implementation */
  int i, j;
  for (i = 0; i < ndl; i++) {
    double sum = 0.0;
    for (j = 0; j < k; j++) {
      sum += zaa[i + j * ldaa] * zz[j];
    }
    vec[i] -= sum;
  }
#endif
}

//***cHACApK_med3
int cHACApK_med3(
  int nl,
  int nr,
  int nlr2)
{
  if(nl < nr) {
    if (nr < nlr2) { return nr; } else if (nlr2 < nl) { return nl; } else { return nlr2; }
  } else {
    if (nlr2 < nr) { return nr;  } else if (nl < nlr2) { return nl; } else { return nlr2; }
  }
}
