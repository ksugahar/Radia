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
*/#ifndef CHACAPK_LIB_H_INCLUDED
#define CHACAPK_LIB_H_INCLUDED

/* Vector norm computation (BLAS dnrm2 optimized) */
extern double cHACApK_unrm_d(int n, double *vec);

/* Extract column with stride access (BLAS dcopy optimized) */
extern void cHACApK_extract_col(double *dst, const double *src, int row, int ldsrc, int k);

/* Find maximum absolute value and location */
extern void cHACApK_maxabsvalloc_d(double *vec, double *maxval, int *loc, int n);

/* Find minimum absolute value and location */
extern void cHACApK_minabsvalloc_d(double *vec, double *minval, int *loc, int n);

/* Dot product subtraction: vec -= zaa * zz */
extern void cHACApK_adotsub_dsm(double *vec, double *zaa, double *zz, int ndl, int k, int ldaa);

/* Median of three */
extern int cHACApK_med3(
  int nl,
  int nr,
  int nlr2);

#endif // CHACAPK_LIB_H_INCLUDED
