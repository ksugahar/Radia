/*
 * C++ Compatible Wrapper for HACApK
 *
 * This header provides C++ compatible type definitions for HACApK.
 * The original cHACApK_base.h uses C-style typedefs that conflict
 * with C++ name lookup rules.
 *
 * Usage: Include this header in C++ code instead of cHACApK_base.h
 *
 * Copyright (c) 2025 Radia Project
 * License: MIT
 */

#ifndef CHACAPK_CPP_H_INCLUDED
#define CHACAPK_CPP_H_INCLUDED

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Structure definitions (must match cHACApK_base.h exactly)
 * Using struct tags directly to avoid typedef conflicts
 */

/* Forward declarations */
struct st_cHACApK_cluster_s;
struct st_cHACApK_leafmtx_s;
struct st_cHACApK_leafmtxp_s;
struct st_cHACApK_lcontrol_s;

/* Cluster structure for hierarchical partitioning */
struct st_cHACApK_cluster_s {
    int ndim;
    int nstrt, nsize, ndpth, nnson, nmbr;
    int ndscd;
    double *bmin;
    double *bmax;
    double zwdth;
    struct st_cHACApK_cluster_s *pc_sons;
};

/* Leaf matrix structure (low-rank or dense block) */
struct st_cHACApK_leafmtx_s {
    int ltmtx;  /* kind of the matrix; 1:rk 2:full */
    int kt;     /* rank for low-rank blocks */
    int nstrtl, ndl;  /* row start and dimension */
    int nstrtt, ndt;  /* column start and dimension */
    double *a1, *a2;  /* U and V for low-rank, or A for dense */

    int nlf; /* number of leaves(sub-matrices) in the MPI process */
    struct st_cHACApK_leafmtx_s **st_lf;  /* Array of pointers (matches C typedef) */
};

/* Leaf matrix collection structure */
struct st_cHACApK_leafmtxp_s {
    int nd;     /* number of unknowns of whole matrix */
    int nlf;    /* number of leaves(sub-matrices) in the MPI process */
    int nlfkt;  /* number of low-rank sub matrices in the MPI process */
    int ktmax;  /* maximum rank */
    int nbl;    /* number of blocks for MPI assignment */
    int nlfalt; /* number of leaves in row(column) of whole matrix */
    int nlfl, nlft;   /* number of leaves in row and column in the MPI process */
    int ndlfs, ndtfs; /* vector sizes in the MPI process */
    struct st_cHACApK_leafmtx_s **st_lf;  /* Array of pointers (matches C typedef) */
    int64_t **lnlfl2g_t;  /* 2D array */
    int *lbstrtl, *lbstrtt; /* Start points of each block in row and column */
    int *lbndl, *lbndt;     /* vector sizes of each block in row and column */
    int *lbndlfs, *lbndtfs; /* vector sizes of each MPI process in row and column */
    int *lbl2t; /* bit vector for receiving data on each MPI process */
};

/* Control structure for HACApK */
struct st_cHACApK_lcontrol_s {
    int *lod, *lsp, *lnp, *lthr, *lpmd;
    double *param, *time;
    int lf_umpi;
};

/* C++ compatible type aliases */
typedef struct st_cHACApK_cluster_s  HACApK_cluster;
typedef struct st_cHACApK_leafmtx_s  HACApK_leafmtx;
typedef struct st_cHACApK_leafmtxp_s HACApK_leafmtxp;
typedef struct st_cHACApK_lcontrol_s HACApK_lcontrol;

/*
 * Entry function callback type
 * Called by HACApK to get matrix element A(i,j)
 * i, j are 1-based indices
 */
typedef double (*HACApK_entry_func)(int i, int j, int i_bemv);

/*
 * HACApK function prototypes for C++ usage
 * These wrap the original C functions
 */

/* Initialize control structure */
void HACApK_init(HACApK_lcontrol *st_ctl, int nd, int print_level);

/* Generate cluster tree */
void HACApK_generate_cbitree_cpp(
    HACApK_cluster **p_st_clt,
    double **zgmid_t,  /* 2D array [ndim+1][md+1] */
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

/* Compute bounding boxes */
void HACApK_bndbox_cpp(
    HACApK_cluster *st_clt,
    double **zgmid_t,
    int *lod,
    int nofc);

/* Count leaf blocks */
void HACApK_count_lntmx_cpp(
    HACApK_cluster *st_cltl,
    HACApK_cluster *st_cltt,
    double *param,
    int *lpmd,
    int *lnmtx,
    int nofc,
    int nffc,
    int *p_ndpth);

/* Generate leaf matrix structure */
void HACApK_generate_leafmtx_cpp(
    HACApK_leafmtx **p_st_leafmtx,
    HACApK_cluster *st_cltl,
    HACApK_cluster *st_cltt,
    double *param,
    int *lpmd,
    int *lnmtx,
    int nofc,
    int nffc,
    int *p_nlf,
    int *p_ndpth);

/* Fill leaf matrices using ACA+ compression */
void HACApK_fill_leafmtx_cpp(
    HACApK_leafmtx *st_lf,
    int i_bemv,
    double *param,
    double znrmmat,
    int *lpmd,
    int *lnmtx,
    int *lodl,
    int *lodt,
    int nd,
    int nlf,
    int *lnps,
    int *lnpe,
    int *lthr,
    HACApK_entry_func entry_func);

/* Sort leaf matrices for efficient traversal */
void HACApK_sort_leafmtx_cpp(HACApK_leafmtx *st_leafmtx, int nlf);

/* Free cluster tree */
void HACApK_free_cluster_cpp(HACApK_cluster *st_clt);

/* Free leaf matrix */
void HACApK_free_leafmtx_cpp(HACApK_leafmtx *st_lf, int nlf);

/* Free control structure */
void HACApK_free_control_cpp(HACApK_lcontrol *st_ctl);

/* ACA+ low-rank approximation */
int HACApK_acaplus_cpp(
    double *zaa,  /* Output: U matrix (ndl x kmax) */
    double *zab,  /* Output: V matrix (ndt x kmax) */
    double *param,
    int ndl,
    int ndt,
    int nstrtl,
    int nstrtt,
    int *lod,
    int i_bemv,
    int kmax,
    double eps,
    double znrmmat,
    double pACA_EPS,
    HACApK_entry_func entry_func);

/*
 * High-level C wrapper functions for Radia integration
 * These use void* to avoid C/C++ struct type compatibility issues
 */

/* Allocate/free opaque structures */
void* HACApK_alloc_leafmtxp(void);
void* HACApK_alloc_lcontrol(void);
void HACApK_free_leafmtxp(void *ptr);
void HACApK_free_lcontrol(void *ptr);

/* Build complete H-matrix from entry function (cHACApK_entry_ij callback) */
int HACApK_build_hmatrix_wrapper(
    void *leafmtxp,            /* Allocated by HACApK_alloc_leafmtxp() */
    void *ctl,                 /* Allocated by HACApK_alloc_lcontrol() */
    double *coordinates,       /* [n_elem * ndim], row-major */
    int n_elem,
    int nffc,                  /* DOF per element (3 for tet, 6 for hex) */
    int ndim,                  /* Spatial dimension (3) */
    double eps,                /* ACA+ tolerance */
    int leaf_size,             /* Minimum cluster size */
    double eta,                /* Admissibility parameter */
    int print_level);

/* Matrix-vector product: y = A * x using H-matrix */
void HACApK_matvec_wrapper(
    void *leafmtxp,
    void *ctl,
    const double *x,
    double *y,
    int nd);

/* Free all H-matrix resources (call before HACApK_free_leafmtxp/lcontrol) */
void HACApK_free_hmatrix_wrapper(
    void *leafmtxp,
    void *ctl);

/* Accessor functions for opaque structure fields */
int HACApK_leafmtxp_get_nd(void *ptr);
int HACApK_leafmtxp_get_nlf(void *ptr);
int HACApK_leafmtxp_get_nlfkt(void *ptr);
int HACApK_leafmtxp_get_ktmax(void *ptr);
int* HACApK_lcontrol_get_lod(void *ptr);

/*
 * Update dense diagonal blocks for nonlinear iteration
 * This is the key function for efficient nonlinear material support.
 * Only dense diagonal blocks (ltmtx==2 && nstrtl==nstrtt) are recomputed.
 * Low-rank off-diagonal blocks remain unchanged.
 *
 * This follows ELF_MAGIC's HACApK_update_diagonal_omp pattern.
 */
void HACApK_update_diagonal_wrapper(
    void *leafmtxp,            /* st_cHACApK_leafmtxp* */
    void *ctl,                 /* st_cHACApK_lcontrol* */
    HACApK_entry_func entry_func);

/**
 * Fast diagonal update: only update true diagonal elements (i==j)
 *
 * This is much faster than HACApK_update_diagonal_wrapper because:
 * - Only updates diagonal entries, not entire blocks
 * - Uses pre-computed N_ii values (passed as diag_N array)
 * - Just adds new inv_chi[i] to stored N_ii values
 *
 * @param leafmtxp  H-matrix leaf pointer
 * @param ctl       H-matrix control structure
 * @param diag_N    Pre-computed diagonal N values (size = ndof)
 * @param inv_chi   New inverse susceptibility values (size = ndof)
 * @param ndof      Total degrees of freedom
 */
void HACApK_update_diagonal_fast_wrapper(
    void *leafmtxp,
    void *ctl,
    const double *diag_N,
    const double *inv_chi,
    int ndof);

#ifdef __cplusplus
}
#endif

#endif /* CHACAPK_CPP_H_INCLUDED */
