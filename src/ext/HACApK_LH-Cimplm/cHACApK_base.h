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
!C**************************************************************************
!C  This file includes basic routines for H-matrices
!C  created by Akihiro Ida at Kyoto University on May 2012
!C  added functions related to ACA+ to HACApK1.0.0 on Nov. 2016
!C  corrected the allocation for st_ctl%lthr on Nov. 2016
!C  added a function related to HACApK_view to HACApK1.1.0 on May 2017
!C  added a function related to writing H-matrix to HACApK1.2.0 on May 2017
!C  added functions related to Block clustering to HACApK1.2.0 on May 2017
!C  translated to C language by Akihiro Ida and Kazuya Goto
!C**************************************************************************
*/
#ifndef CHACAPK_BASE_H_INCLUDED
#define CHACAPK_BASE_H_INCLUDED

#include <stdint.h>

//
// typedefs
//

typedef struct st_cHACApK_cluster *st_cHACApK_cluster;
typedef struct st_cHACApK_leafmtx *st_cHACApK_leafmtx;
typedef struct st_cHACApK_leafmtxp *st_cHACApK_leafmtxp;
typedef struct st_cHACApK_lcontrol *st_cHACApK_lcontrol;

//
// definition of structs
//

//*** struct st_cHACApK_cluster
struct st_cHACApK_cluster {
  int ndim;
  int nstrt, nsize, ndpth, nnson, nmbr;
  int ndscd;
  double *bmin;
  double *bmax;
  double zwdth;
  st_cHACApK_cluster *pc_sons;
};

//*** struct st_cHACApK_leafmtx
struct st_cHACApK_leafmtx {
  int ltmtx;  // kind of the matrix; 1:rk 2:full
  int kt;
  int nstrtl,ndl;
  int nstrtt,ndt;
  double *a1, *a2;

  int nlf; // number of leaves(sub-matrices) in the MPI process
  st_cHACApK_leafmtx *st_lf;
};

//*** struct st_cHACApK_leafmtxp
struct st_cHACApK_leafmtxp {
  int nd; // number of unknowns of whole matrix
  int nlf; // number of leaves(sub-matrices) in the MPI process
  int nlfkt; // number of low-rank sub matrices in the MPI process
  int ktmax;
  int nbl; //number of blocks for MPI assignment
  int nlfalt; //number of leaves in row(column) of whole matrix
  int nlfl,nlft;  // number of leaves in row and column in the MPI process
  int ndlfs,ndtfs;  // vector sizes in the MPI process
  st_cHACApK_leafmtx *st_lf;
  int64_t **lnlfl2g_t;  // 2D array
  int *lbstrtl, *lbstrtt; // Start points of each block in row and column
  int *lbndl, *lbndt; // vector sizes of each block in row and column
  int *lbndlfs, *lbndtfs; // vector sizes of each MPI process in row and column
  int *lbl2t; // bit vector for recieving data on each MPI process
};

//*** struct st_cHACApK_lcontrol
struct st_cHACApK_lcontrol {
  int *lod, *lsp, *lnp, *lthr, *lpmd;
  double *param, *time;
  int lf_umpi;
};

//
// prototype of functions
//

extern void cHACApK_generate_frame_blrleaf(
  st_cHACApK_leafmtxp st_leafmtxp,
  int ibemv,
  st_cHACApK_lcontrol st_ctl,
  double **gmid_t,  // 2D array [ndim+1][nofc+1]
  int lnmtx[4+1],
  int nofc,
  int nffc,
  int ndim);

extern void cHACApK_setcutthread(
  int *lthr,
  st_cHACApK_leafmtxp st_leafmtxp,
  st_cHACApK_lcontrol st_ctl,
  int64_t mem8,
  int nthr,
  int ktp);

extern void cHACApK_count_blrnmb(
  st_cHACApK_cluster st_cltl,
  st_cHACApK_cluster st_cltt,
  double *param,
  int *lpmd,
  int *lnmtx,
  int nofc,
  int nffc,
  int *p_ndpth);

extern void cHACApK_count_blrleaf(
  st_cHACApK_leafmtx *st_leafmtx,
  st_cHACApK_cluster st_cltl,
  st_cHACApK_cluster st_cltt,
  double *param,
  int *lpmd,
  int *lnmtx,
  int nofc,
  int nffc,
  int *p_ndpth);

extern void cHACApK_generate_blrleaf(
  st_cHACApK_leafmtx *st_leafmtx,
  st_cHACApK_cluster st_cltl,
  st_cHACApK_cluster st_cltt,
  double *param,
  int *lpmd,
  int *lnmtx,
  int nofc,
  int nffc,
  int *p_nlf,
  int *p_ndpth);

extern void cHACApK_count_lntmx(
  st_cHACApK_cluster st_cltl,
  st_cHACApK_cluster st_cltt,
  double *param,
  int *lpmd,
  int *lnmtx,
  int nofc,
  int nffc,
  int *p_ndpth);

extern void cHACApK_generate_leafmtx(
  st_cHACApK_leafmtx *st_leafmtx,
  st_cHACApK_cluster st_cltl,
  st_cHACApK_cluster st_cltt,
  double *param,
  int *lpmd,
  int *lnmtx,
  int nofc,
  int nffc,
  int *p_nlf,
  int *p_ndpth);

extern void cHACApK_sort_leafmtx(
  st_cHACApK_leafmtx *st_leafmtx,
  int nlf);

extern void cHACApK_qsort_col_leafmtx(
  st_cHACApK_leafmtx *st_leafmtx,
  int nlf_s,
  int nlf_e);

extern void cHACApK_qsort_row_leafmtx(
  st_cHACApK_leafmtx *st_leafmtx,
  int nlf_s,
  int nlf_e);

extern void cHACApK_free_st_clt(
  st_cHACApK_cluster st_clt);

extern st_cHACApK_cluster cHACApK_generate_cluster(
  int *p_nmbr,
  int ndpth,
  int nstrt,
  int nsize,
  int ndim,
  int nson);

extern void cHACApK_bndbox(
  st_cHACApK_cluster st_clt,
  double **zgmid_t, // 2D array [st_clt->ndim+1][nofc+1]
  int *lod,
  int nofc);

extern void cHACApK_generate_cbitree(
  st_cHACApK_cluster *st_clt,
  double **zgmid_t, // 2D array [ndim+1][md+1]
  double *param,
  int *lpmd,
  int *lod,
  int *p_ndpth,
  int ndscd,
  int nsrt,
  int nd,
  int md,
  int ndim,
  int *p_nclst);

#endif // CHACAPK_BASE_H_INCLUDED
