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

//***cHACApK_unrm_d
// Compute the Euclidean norm of a vector
double cHACApK_unrm_d(int n, double *vec) {
  double norm = 0.0;
  int i;
  for (i = 0; i < n; i++) {
    norm += vec[i] * vec[i];
  }
  return sqrt(norm);
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
// vec(1:ndl) = vec(1:ndl) - zaa(1:ndl, 1:k-1) * zz(1:k-1)
void cHACApK_adotsub_dsm(double *vec, double *zaa, double *zz, int ndl, int k, int ldaa) {
  int i, j;
  for (i = 0; i < ndl; i++) {
    double sum = 0.0;
    for (j = 0; j < k; j++) {
      sum += zaa[i + j * ldaa] * zz[j];
    }
    vec[i] -= sum;
  }
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
