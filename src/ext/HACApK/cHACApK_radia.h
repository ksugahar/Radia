/*
 * HACApK for Radia - Simplified H-matrix with ACA+ for single-process execution
 *
 * Based on HACApK (ppOpen-HPC project, MIT License)
 * Copyright (c) 2015 Akihiro Ida and Takeshi Iwashita
 *
 * Simplified for Radia by removing MPI support and keeping OpenMP only.
 * Logging unified with Radia conventions.
 *
 * Copyright (c) 2025 Radia Project
 */

#ifndef CHACAPK_RADIA_H_INCLUDED
#define CHACAPK_RADIA_H_INCLUDED

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Forward declarations for opaque structures
 */
typedef struct st_cHACApK_cluster_impl *st_cHACApK_cluster;
typedef struct st_cHACApK_leafmtx_impl *st_cHACApK_leafmtx;
typedef struct st_cHACApK_leafmtxp_impl *st_cHACApK_leafmtxp;
typedef struct st_cHACApK_lcontrol_impl *st_cHACApK_lcontrol;

/*
 * Configuration parameters for H-matrix construction
 */
typedef struct {
    double aca_eps;      /* ACA+ compression tolerance (default: 1e-4) */
    int    leaf_size;    /* Minimum cluster size (default: 32) */
    double eta;          /* Admissibility parameter (default: 2.0) */
    int    max_rank;     /* Maximum ACA rank (default: 200) */
    int    print_level;  /* 0=silent, 1=standard, 2=debug */
    int    num_threads;  /* Number of OpenMP threads (0=auto) */
} HACApK_config;

/*
 * Statistics from H-matrix construction
 */
typedef struct {
    int n_lowrank;       /* Number of low-rank (ACA compressed) blocks */
    int n_dense;         /* Number of dense blocks */
    int max_rank;        /* Maximum rank among low-rank blocks */
    int n_leaves;        /* Total number of leaf blocks */
    int n_dof;           /* Total degrees of freedom */
    double compression;  /* Compression ratio (memory used / full matrix) */
    double build_time;   /* Time to build H-matrix (seconds) */
} HACApK_stats;

/*
 * Callback function type for computing matrix elements
 *
 * Parameters:
 *   i     - Row index (1-based)
 *   j     - Column index (1-based)
 *   bemv  - User-defined parameter (passed from build function)
 *
 * Returns:
 *   Matrix element A(i,j)
 */
typedef double (*HACApK_entry_func)(int i, int j, int bemv);

/*
 * Initialize default configuration
 */
void HACApK_init_config(HACApK_config *config);

/*
 * Allocate H-matrix structures
 * Returns NULL on failure
 */
st_cHACApK_leafmtxp HACApK_alloc_leafmtxp(void);
st_cHACApK_lcontrol HACApK_alloc_lcontrol(void);

/*
 * Free H-matrix structures
 */
void HACApK_free_leafmtxp(st_cHACApK_leafmtxp ptr);
void HACApK_free_lcontrol(st_cHACApK_lcontrol ptr);

/*
 * Build H-matrix from coordinates and entry callback
 *
 * Parameters:
 *   leafmtxp    - Allocated leafmtxp structure
 *   ctl         - Allocated control structure
 *   coordinates - Element center coordinates [n_elem * ndim], row-major
 *   n_elem      - Number of elements
 *   nffc        - DOF per element (3 for standard, 6 for hex MSC)
 *   ndim        - Spatial dimension (3)
 *   config      - Configuration parameters
 *   entry_func  - Callback to compute matrix elements
 *   bemv        - User parameter passed to entry_func
 *   stats       - Output statistics (can be NULL)
 *
 * Returns:
 *   0 on success, non-zero on error
 */
int HACApK_build(
    st_cHACApK_leafmtxp leafmtxp,
    st_cHACApK_lcontrol ctl,
    const double *coordinates,
    int n_elem,
    int nffc,
    int ndim,
    const HACApK_config *config,
    HACApK_entry_func entry_func,
    int bemv,
    HACApK_stats *stats);

/*
 * Matrix-vector product: y = A * x
 *
 * Parameters:
 *   leafmtxp - Built H-matrix
 *   ctl      - Control structure
 *   x        - Input vector [n_dof]
 *   y        - Output vector [n_dof]
 *   n_dof    - Vector size (leafmtxp->nd)
 */
void HACApK_matvec(
    st_cHACApK_leafmtxp leafmtxp,
    st_cHACApK_lcontrol ctl,
    const double *x,
    double *y,
    int n_dof);

/*
 * Get DOF permutation array (cluster reordering)
 * Returns pointer to internal array, do not free
 */
const int* HACApK_get_permutation(st_cHACApK_lcontrol ctl);

/*
 * Get total DOF count
 */
int HACApK_get_ndof(st_cHACApK_leafmtxp leafmtxp);

/*
 * Utility functions
 */

/* Compute Euclidean norm of vector */
double HACApK_norm2(const double *vec, int n);

/* Find max absolute value and its location */
void HACApK_maxabs_loc(const double *vec, int n, double *maxval, int *loc);

/* Find min absolute value and its location */
void HACApK_minabs_loc(const double *vec, int n, double *minval, int *loc);

/* vec -= A * x (for ACA+ orthogonalization) */
void HACApK_gemv_sub(double *vec, const double *A, const double *x, int nrow, int ncol, int lda);

/* Median of three (for quicksort) */
int HACApK_med3(int a, int b, int c);

#ifdef __cplusplus
}
#endif

#endif /* CHACAPK_RADIA_H_INCLUDED */
