/*
 * C++ Compatible Wrapper for HACApK
 *
 * This header provides C++ compatible interface for HACApK.
 * It includes cHACApK_base.h which now uses standard C/C++ compatible
 * struct definitions (st_cHACApK_*_t types with separate pointer typedefs).
 *
 * Usage: Include this header in C++ code for HACApK integration.
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

/* Include the base HACApK header which defines all structures */
#include "cHACApK_base.h"

/*
 * C++ compatible type aliases
 * These provide cleaner names for C++ usage
 */
typedef st_cHACApK_cluster_t  HACApK_cluster;
typedef st_cHACApK_leafmtx_t  HACApK_leafmtx;
typedef st_cHACApK_leafmtxp_t HACApK_leafmtxp;
typedef st_cHACApK_lcontrol_t HACApK_lcontrol;

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

/* Build H-matrix with variable DOF per element (for mixed hex+tetra meshes) */
int HACApK_build_hmatrix_varDOF_wrapper(
    void *leafmtxp,            /* Allocated by HACApK_alloc_leafmtxp() */
    void *ctl,                 /* Allocated by HACApK_alloc_lcontrol() */
    double *coordinates,       /* [n_elem * ndim], row-major */
    int n_elem,
    int *dof_offset,           /* [n_elem + 1], cumulative DOF offset (0-based) */
    int total_dof,             /* Total DOF count = dof_offset[n_elem] */
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

/* Transpose matvec: y = A^T x (mirror of HACApK_matvec_wrapper with l/t roles swapped). */
void HACApK_matvec_transpose_wrapper(
    void *leafmtxp,
    void *ctl,
    const double *x,
    double *y,
    int nd);

/* Symmetric matvec: y = G_sym x, G_sym EXACTLY symmetric (built from the upper-triangular leaves;
 * the lower triangle is the exact transpose of the upper).  Valid only for a symmetric cluster tree
 * (rows == cols, one geometry -- e.g. the charge Gram), where the leaf partition is symmetric. */
void HACApK_matvec_sym_wrapper(
    void *leafmtxp,
    void *ctl,
    const double *x,
    double *y,
    int nd);

/* Row-major batch [nrhs][nd].  Traverses every symmetric H-matrix leaf once
 * and uses BLAS-3 GEMM across right-hand sides. */
void HACApK_matvec_sym_many_wrapper(
    void *leafmtxp,
    void *ctl,
    const double *x,
    double *y,
    int nd,
    int nrhs);

/* Symmetric row-major batch restricted to an active principal submatrix.
 * active_prefix is the inclusive-prefix count in HACApK's permuted ordering
 * (length nd+1, active_prefix[0]=0).  A leaf is skipped exactly when either
 * its row or column range contains no active entry. */
void HACApK_matvec_sym_many_masked_wrapper(
    void *leafmtxp,
    void *ctl,
    const double *x,
    double *y,
    int nd,
    int nrhs,
    const int *active_prefix);

/* Prepared symmetric batch apply.  active_prefix and diagonal_scale are
 * optional.  Scaling is fused into permutation/reduction so callers can
 * evaluate S*A*S without allocating two row-major scratch batches. */
void HACApK_matvec_sym_many_prepared_wrapper(
    void *leafmtxp,
    void *ctl,
    const double *x,
    double *y,
    int nd,
    int nrhs,
    const int *active_prefix,
    const double *diagonal_scale);

/* Optional matvec profiler. Enabled by RADIA_HDIV_HMATVEC_STATS=1.
 * values[0..7] = total_s, zero_s, permute_s, leaf_s, reduce_s, meta_s,
 *                lowrank_flop_est, dense_flop_est.
 * counts[0..19] = calls, lowrank_leaves, dense_leaves,
 * mirrored_upper_leaves, diagonal_leaves, skipped_lower_leaves, last_nd,
 * last_nthr, lowrank_upper_leaves, dense_upper_leaves,
 * inactive_skipped_leaves, lowrank_directions, dense_directions, gemm_calls,
 * lowrank_rank_sum, lowrank_rank_max, and cumulative direction counts for
 * ranks <= 4, 8, 16, and 32.
 */
void HACApK_matvec_stats_reset(void);
void HACApK_matvec_stats_get(double *values, int n_values,
                             int64_t *counts, int n_counts);

/* Free all H-matrix resources (call before HACApK_free_leafmtxp/lcontrol) */
void HACApK_free_hmatrix_wrapper(
    void *leafmtxp,
    void *ctl);

/* Reset all HACApK global state (call between solves) */
void HACApK_reset_global_state(void);

/* Accessor functions for opaque structure fields */
int HACApK_leafmtxp_get_nd(void *ptr);
int HACApK_leafmtxp_get_nlf(void *ptr);
int HACApK_leafmtxp_get_nlfkt(void *ptr);
int HACApK_leafmtxp_get_ktmax(void *ptr);
int* HACApK_lcontrol_get_lod(void *ptr);

/**
 * Calculate actual H-matrix memory usage (ELF-compatible)
 *
 * Returns memory in bytes by iterating over all leaf blocks:
 * - Low-rank (ltmtx=1): (ndl + ndt) * kt * sizeof(double)
 * - Dense (ltmtx=2):    ndl * ndt * sizeof(double)
 *
 * @param leafmtxp        H-matrix leaf pointer
 * @param hmat_bytes_out  Output: actual H-matrix memory [bytes]
 * @param dense_bytes_out Output: equivalent dense memory [bytes]
 */
void HACApK_get_memory_stats(void *leafmtxp,
                              int64_t *hmat_bytes_out,
                              int64_t *dense_bytes_out);

/* Accessor functions for current lod array during H-matrix build
 * These are set by HACApK_build_hmatrix_varDOF_wrapper before fill and
 * cleared after fill, allowing the C++ callback to access permutation info.
 * lod[permuted_1based] = original_1based
 */
int* HACApK_get_current_lod(void);
int HACApK_get_current_lod_size(void);

/*
 * Update dense diagonal blocks for nonlinear iteration
 * This is the key function for efficient nonlinear material support.
 * Only dense diagonal blocks (ltmtx==2 && nstrtl==nstrtt) are recomputed.
 * Low-rank off-diagonal blocks remain unchanged.
 *
 * This follows ELF_MAGIC's HACApK_update_diagonal_omp pattern.
 */
void HACApK_update_diagonal_wrapper(
    void *leafmtxp,            /* st_cHACApK_leafmtxp_t* */
    void *ctl,                 /* st_cHACApK_lcontrol_t* */
    HACApK_entry_func entry_func);

/*
 * Fast diagonal update using pre-computed N_ii values
 * Only updates true diagonal entries (i==j) using: a1[il + ndt*il] = diag_N[orig_i-1] - inv_chi[orig_i-1]
 * This is O(ndof) instead of O(block_size^2 * n_diag_blocks)
 */
void HACApK_update_diagonal_fast_wrapper(
    void *leafmtxp,            /* st_cHACApK_leafmtxp_t* */
    void *ctl,                 /* st_cHACApK_lcontrol_t* */
    const double *diag_N,      /* Pre-computed N_ii diagonal elements [ndof] */
    const double *inv_chi,     /* Inverse susceptibility diagonal [ndof] */
    int ndof);                 /* Total DOF count */

/**
 * Set entry function for cHACApK_entry_ij callback
 * Must be called before HACApK_build_hmatrix_wrapper
 */
void HACApK_set_entry_func(HACApK_entry_func func);

/**
 * Clear entry function callback
 */
void HACApK_clear_entry_func(void);

#ifdef __cplusplus
}
#endif

#endif /* CHACAPK_CPP_H_INCLUDED */
