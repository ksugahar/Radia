/*
 * C++ Compatible Wrapper Implementation for HACApK
 *
 * This file provides C wrapper functions that bridge C++ code to the
 * original HACApK C implementation. It compiles as pure C to avoid
 * the typedef struct * conflicts that occur in C++.
 *
 * Copyright (c) 2025 Radia Project
 * License: MIT
 */

#include "cHACApK_base.h"
#include "cHACApK_calc_entry_ij.h"
#include "cHACApK_lib.h"
#include "mpi_stub.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <mkl.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

#include "rad_hacapk_parallel.h"

/*=========================================================================
 * BLAS declarations (using Intel MKL LAPACK conventions)
 * These are provided by the system BLAS library linked via CMake
 *=========================================================================*/
#ifdef __cplusplus
extern "C" {
#endif

/* BLAS Level 2: Matrix-vector multiply
 * y = alpha * A * x + beta * y  (trans='N')
 * y = alpha * A^T * x + beta * y  (trans='T')
 */
void dgemv_(const char *trans, const int *m, const int *n,
            const double *alpha, const double *a, const int *lda,
            const double *x, const int *incx,
            const double *beta, double *y, const int *incy);
void dgemm_(const char *transa, const char *transb,
            const int *m, const int *n, const int *k,
            const double *alpha, const double *a, const int *lda,
            const double *b, const int *ldb, const double *beta,
            double *c, const int *ldc);

#ifdef __cplusplus
}
#endif

/*
 * Global lod pointer for C++ callback access during H-matrix build
 * Set before cHACApK_fill_leafmtx_hyp, cleared after
 */
static int* g_hacapk_lod = NULL;
static int g_hacapk_lod_size = 0;

/* Accessor functions for C++ code */
int* HACApK_get_current_lod(void) {
    return g_hacapk_lod;
}

int HACApK_get_current_lod_size(void) {
    return g_hacapk_lod_size;
}

/*
 * Note: This file is compiled as C (not C++)
 * The struct types here use the original HACApK typedefs.
 *
 * The wrapper functions take void* parameters that are cast to the
 * correct HACApK types internally.
 */

/* Helper to get number of threads (via TaskManager) */
static int get_num_threads(void) {
    return hacapk_get_num_threads();
}

/*=========================================================================
 * Persistent work arrays for matvec (avoid repeated allocation)
 * ELF-style: allocate once, reuse across multiple matvec calls
 *=========================================================================*/
static double *g_x_perm = NULL;      /* Permuted input vector */
static double *g_y_perm = NULL;      /* Permuted output vector */
static double **g_y_thread = NULL;   /* Thread-local y arrays */
static double **g_tmp_vec = NULL;    /* Thread-local tmp vectors */
static int g_matvec_nd = 0;          /* Size of allocated arrays */
static int g_matvec_nthr = 0;        /* Number of threads */
static int g_matvec_ktmax = 0;       /* Max rank for tmp_vec */
static double *g_batch_x_perm = NULL;    /* column-major [nd,nrhs] */
static double **g_batch_y_thread = NULL; /* thread-local [nd,nrhs] */
static double **g_batch_tmp = NULL;      /* thread-local [ktmax,nrhs] */
static int g_batch_nd = 0;
static int g_batch_nrhs = 0;
static int g_batch_nthr = 0;
static int g_batch_ktmax = 0;

typedef struct {
    int64_t calls;
    int64_t lowrank_leaves;
    int64_t dense_leaves;
    int64_t mirrored_upper_leaves;
    int64_t diagonal_leaves;
    int64_t skipped_lower_leaves;
    int64_t last_nd;
    int64_t last_nthr;
    double total_s;
    double zero_s;
    double permute_s;
    double leaf_s;
    double reduce_s;
    double meta_s;
    double lowrank_flop_est;
    double dense_flop_est;
} HACApKMatvecStats;

static HACApKMatvecStats g_matvec_stats;
static int g_matvec_stats_enabled_cache = -1;
static int g_matvec_mkl_threads_cache = -2;

static int matvec_stats_enabled(void) {
    if (g_matvec_stats_enabled_cache < 0) {
        const char *env = getenv("RADIA_HDIV_HMATVEC_STATS");
        g_matvec_stats_enabled_cache = (env && env[0] && env[0] != '0') ? 1 : 0;
    }
    return g_matvec_stats_enabled_cache;
}

static double matvec_wall_time(void) {
#ifdef _WIN32
    static LARGE_INTEGER freq;
    static int initialized = 0;
    LARGE_INTEGER now;
    if (!initialized) {
        QueryPerformanceFrequency(&freq);
        initialized = 1;
    }
    QueryPerformanceCounter(&now);
    return (double)now.QuadPart / (double)freq.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1.0e-9 * (double)ts.tv_nsec;
#endif
}

static int matvec_mkl_threads(void) {
    if (g_matvec_mkl_threads_cache == -2) {
        const char *env = getenv("RADIA_HACAPK_MATVEC_MKL_THREADS");
        if (env && env[0]) {
            g_matvec_mkl_threads_cache = atoi(env);
        } else {
            g_matvec_mkl_threads_cache = 1;
        }
    }
    return g_matvec_mkl_threads_cache;
}

static void matvec_mkl_begin(void) {
    int nth = matvec_mkl_threads();
    if (nth > 0) mkl_set_num_threads_local(nth);
}

static void matvec_mkl_end(void) {
    if (matvec_mkl_threads() > 0) mkl_set_num_threads_local(0);
}

void HACApK_matvec_stats_reset(void) {
    memset(&g_matvec_stats, 0, sizeof(g_matvec_stats));
    g_matvec_stats_enabled_cache = -1;
    g_matvec_mkl_threads_cache = -2;
}

void HACApK_matvec_stats_get(double *values, int n_values,
                             int64_t *counts, int n_counts) {
    int i;
    if (values) {
        for (i = 0; i < n_values; ++i) values[i] = 0.0;
        if (n_values > 0) values[0] = g_matvec_stats.total_s;
        if (n_values > 1) values[1] = g_matvec_stats.zero_s;
        if (n_values > 2) values[2] = g_matvec_stats.permute_s;
        if (n_values > 3) values[3] = g_matvec_stats.leaf_s;
        if (n_values > 4) values[4] = g_matvec_stats.reduce_s;
        if (n_values > 5) values[5] = g_matvec_stats.meta_s;
        if (n_values > 6) values[6] = g_matvec_stats.lowrank_flop_est;
        if (n_values > 7) values[7] = g_matvec_stats.dense_flop_est;
    }
    if (counts) {
        for (i = 0; i < n_counts; ++i) counts[i] = 0;
        if (n_counts > 0) counts[0] = g_matvec_stats.calls;
        if (n_counts > 1) counts[1] = g_matvec_stats.lowrank_leaves;
        if (n_counts > 2) counts[2] = g_matvec_stats.dense_leaves;
        if (n_counts > 3) counts[3] = g_matvec_stats.mirrored_upper_leaves;
        if (n_counts > 4) counts[4] = g_matvec_stats.diagonal_leaves;
        if (n_counts > 5) counts[5] = g_matvec_stats.skipped_lower_leaves;
        if (n_counts > 6) counts[6] = g_matvec_stats.last_nd;
        if (n_counts > 7) counts[7] = g_matvec_stats.last_nthr;
    }
}

/* Initialize persistent matvec buffers */
static void init_matvec_buffers(int nd, int nthr, int ktmax) {
    int i;
    if (g_matvec_nd != nd || g_matvec_nthr != nthr || g_matvec_ktmax < ktmax) {
        /* Free old buffers */
        if (g_x_perm) { free(g_x_perm); g_x_perm = NULL; }
        if (g_y_perm) { free(g_y_perm); g_y_perm = NULL; }
        if (g_y_thread) {
            for (i = 0; i < g_matvec_nthr; i++) {
                if (g_y_thread[i]) free(g_y_thread[i]);
            }
            free(g_y_thread);
            g_y_thread = NULL;
        }
        if (g_tmp_vec) {
            for (i = 0; i < g_matvec_nthr; i++) {
                if (g_tmp_vec[i]) free(g_tmp_vec[i]);
            }
            free(g_tmp_vec);
            g_tmp_vec = NULL;
        }

        /* Allocate new buffers */
        g_x_perm = (double*)malloc(sizeof(double) * nd);
        g_y_perm = (double*)malloc(sizeof(double) * nd);
        g_y_thread = (double**)malloc(sizeof(double*) * nthr);
        g_tmp_vec = (double**)malloc(sizeof(double*) * nthr);
        for (i = 0; i < nthr; i++) {
            g_y_thread[i] = (double*)malloc(sizeof(double) * nd);
            g_tmp_vec[i] = (double*)malloc(sizeof(double) * (ktmax > 0 ? ktmax : 1));
        }

        g_matvec_nd = nd;
        g_matvec_nthr = nthr;
        g_matvec_ktmax = ktmax;
    }
}

static void free_batch_matvec_buffers(void) {
    int i;
    if (g_batch_x_perm) { free(g_batch_x_perm); g_batch_x_perm = NULL; }
    if (g_batch_y_thread) {
        for (i = 0; i < g_batch_nthr; ++i) free(g_batch_y_thread[i]);
        free(g_batch_y_thread); g_batch_y_thread = NULL;
    }
    if (g_batch_tmp) {
        for (i = 0; i < g_batch_nthr; ++i) free(g_batch_tmp[i]);
        free(g_batch_tmp); g_batch_tmp = NULL;
    }
    g_batch_nd = g_batch_nrhs = g_batch_nthr = g_batch_ktmax = 0;
}

static void init_batch_matvec_buffers(int nd, int nrhs, int nthr, int ktmax) {
    int i;
    if (g_batch_nd == nd && g_batch_nrhs == nrhs &&
        g_batch_nthr == nthr && g_batch_ktmax >= ktmax) return;
    free_batch_matvec_buffers();
    g_batch_x_perm = (double*)malloc(sizeof(double)*(size_t)nd*nrhs);
    g_batch_y_thread = (double**)malloc(sizeof(double*)*nthr);
    g_batch_tmp = (double**)malloc(sizeof(double*)*nthr);
    for (i = 0; i < nthr; ++i) {
        g_batch_y_thread[i] = (double*)malloc(sizeof(double)*(size_t)nd*nrhs);
        g_batch_tmp[i] = (double*)malloc(
            sizeof(double)*(size_t)(ktmax > 0 ? ktmax : 1)*nrhs);
    }
    g_batch_nd = nd; g_batch_nrhs = nrhs; g_batch_nthr = nthr;
    g_batch_ktmax = ktmax;
}

/* Free persistent matvec buffers (call when H-matrix is destroyed) */
static void free_matvec_buffers(void) {
    int i;
    if (g_x_perm) { free(g_x_perm); g_x_perm = NULL; }
    if (g_y_perm) { free(g_y_perm); g_y_perm = NULL; }
    if (g_y_thread) {
        for (i = 0; i < g_matvec_nthr; i++) {
            if (g_y_thread[i]) free(g_y_thread[i]);
        }
        free(g_y_thread);
        g_y_thread = NULL;
    }
    if (g_tmp_vec) {
        for (i = 0; i < g_matvec_nthr; i++) {
            if (g_tmp_vec[i]) free(g_tmp_vec[i]);
        }
        free(g_tmp_vec);
        g_tmp_vec = NULL;
    }
    g_matvec_nd = 0;
    g_matvec_nthr = 0;
    g_matvec_ktmax = 0;
    free_batch_matvec_buffers();
}

/* Public function to reset all HACApK global state */
void HACApK_reset_global_state(void) {
    /* Free persistent matvec buffers */
    free_matvec_buffers();

    /* Clear lod state */
    g_hacapk_lod = NULL;
    g_hacapk_lod_size = 0;

    /* Clear optional profiling state. */
    HACApK_matvec_stats_reset();
}

/*=========================================================================
 * High-level wrapper: Build complete H-matrix
 *=========================================================================*/

int HACApK_build_hmatrix_wrapper(
    void *leafmtxp_void,           /* st_cHACApK_leafmtxp* */
    void *ctl_void,                /* st_cHACApK_lcontrol* */
    double *coordinates,           /* [n_elem * ndim], row-major */
    int n_elem,
    int nffc,                      /* DOF per element (3 for tet, 6 for hex) */
    int ndim,                      /* Spatial dimension (3) */
    double eps,                    /* ACA+ tolerance */
    int leaf_size,                 /* Minimum cluster size */
    double eta,                    /* Admissibility parameter */
    int print_level)
{
    st_cHACApK_leafmtxp leafmtxp = (st_cHACApK_leafmtxp)leafmtxp_void;
    st_cHACApK_lcontrol ctl = (st_cHACApK_lcontrol)ctl_void;
    st_cHACApK_cluster st_clt;
    st_cHACApK_leafmtx *st_leafmtx;
    double **gmid_t;
    int *lodfc;
    int nd, nofc;
    int lnmtx[5];
    int ndpth, nclst, nlf;
    int i_bemv = 0;
    int il, ig, is, ip;
    double znrmmat;
    int nthr;

    nofc = n_elem;
    nd = nofc * nffc;
    nthr = get_num_threads();

    if (print_level > 0) {
        printf("[HACApK] Building H-matrix: n_elem=%d, nffc=%d, nd=%d\n", n_elem, nffc, nd);
        printf("[HACApK] Parameters: eps=%.2e, leaf_size=%d, eta=%.2f\n", eps, leaf_size, eta);
    }

    /* Initialize leaf matrix pointer structure */
    leafmtxp->nd = nd;
    leafmtxp->nlf = 0;
    leafmtxp->nlfkt = 0;
    leafmtxp->ktmax = 0;
    leafmtxp->st_lf = NULL;

    /* Allocate and set control parameters */
    ctl->param = (double*)calloc(101, sizeof(double));
    ctl->lpmd = (int*)calloc(100, sizeof(int));
    ctl->time = (double*)calloc(100, sizeof(double));
    ctl->lod = (int*)calloc(nd + 1, sizeof(int));
    ctl->lthr = (int*)calloc(nthr + 10, sizeof(int));
    ctl->lsp = NULL;
    ctl->lnp = NULL;

    if (!ctl->param || !ctl->lpmd || !ctl->lod || !ctl->lthr) {
        fprintf(stderr, "[HACApK] Error: Memory allocation failed\n");
        return -1;
    }

    /* Set HACApK parameters */
    ctl->param[1] = (double)print_level;   /* Print level */
    ctl->param[21] = (double)leaf_size;    /* Minimum leaf size */
    ctl->param[22] = 1.0;                   /* Max leaf factor */
    ctl->param[41] = 1.0;                   /* npgl (sqrt of MPI size) */
    ctl->param[42] = 0;                     /* Block size (auto) */
    ctl->param[43] = 1.0;                   /* Block division factor */
    ctl->param[51] = eta;                   /* Admissibility (eta) */
    ctl->param[60] = 1;                     /* ACA mode (1=ACA, 2=ACA+) */
    ctl->param[61] = 1;                     /* ACA norm (MREM) */
    ctl->param[62] = (double)200;          /* Max rank initial */
    ctl->param[63] = (double)200;          /* Max rank */
    ctl->param[64] = 1;                     /* Min rank (ELF uses 1) */
    ctl->param[71] = eps;                   /* ACA tolerance */
    /* param[72]: ACA_EPS multiplier (HACApK standard: 1.0e-3, LatticeH: 1.0e-9) */
    /* ELF uses standard HACApK (not LatticeH), so use 1.0e-3 for compatibility */
    ctl->param[72] = 1.0e-3;

    /* MPI stub setup (single process) */
    ctl->lpmd[1] = MPI_COMM_WORLD;   /* Communicator */
    ctl->lpmd[2] = 1;                 /* Number of MPI processes */
    ctl->lpmd[3] = 0;                 /* MPI rank */
    ctl->lpmd[4] = 0;                 /* MPI log */
    ctl->lpmd[20] = nthr;             /* Number of threads */

    /* Initialize DOF ordering (identity permutation) */
    for (il = 1; il <= nd; il++) {
        ctl->lod[il] = il;
    }

    /* Allocate 2D coordinate array for HACApK [ndim+1][nofc+1] */
    gmid_t = (double**)malloc(sizeof(double*) * (ndim + 1));
    for (il = 0; il <= ndim; il++) {
        gmid_t[il] = (double*)malloc(sizeof(double) * (nofc + 1));
    }

    /* Fill coordinates (1-based indexing) */
    for (il = 0; il < nofc; il++) {
        for (ig = 0; ig < ndim; ig++) {
            gmid_t[ig + 1][il + 1] = coordinates[il * ndim + ig];
        }
    }

    /* Allocate temporary element ordering array */
    lodfc = (int*)malloc(sizeof(int) * (nofc + 1));
    for (il = 1; il <= nofc; il++) {
        lodfc[il] = il;
    }

    /*=========================================================================
     * Generate cluster tree
     *=========================================================================*/
    ndpth = 0;
    nclst = 0;

    cHACApK_generate_cbitree(&st_clt, gmid_t, ctl->param, ctl->lpmd, lodfc,
                              &ndpth, 0, 1, nofc, nofc, ndim, &nclst);

    if (print_level > 0) {
        printf("[HACApK] Cluster tree: nclst=%d, ndpth=%d\n", nclst, ndpth);
    }

    /* Compute bounding boxes */
    cHACApK_bndbox(st_clt, gmid_t, lodfc, nofc);

    /* Map element ordering to DOF ordering */
    for (il = 1; il <= nofc; il++) {
        for (ig = 1; ig <= nffc; ig++) {
            is = ig + (il - 1) * nffc;
            ctl->lod[is] = (lodfc[il] - 1) * nffc + ig;
        }
    }

    /*=========================================================================
     * Count leaf matrices
     *=========================================================================*/
    for (il = 0; il < 5; il++) lnmtx[il] = 0;
    ndpth = 0;
    cHACApK_count_lntmx(st_clt, st_clt, ctl->param, ctl->lpmd, lnmtx, nofc, nffc, &ndpth);

    nlf = lnmtx[1] + lnmtx[2];  /* Total leaves = low-rank + dense */
    if (print_level > 0) {
        printf("[HACApK] Leaf count: low-rank=%d, dense=%d, total=%d\n",
               lnmtx[1], lnmtx[2], nlf);
    }

    /*=========================================================================
     * Generate leaf matrix structure
     *=========================================================================*/
    st_leafmtx = (st_cHACApK_leafmtx*)calloc(nlf + 1, sizeof(st_cHACApK_leafmtx));
    if (!st_leafmtx) {
        fprintf(stderr, "[HACApK] Error: Memory allocation for st_leafmtx failed\n");
        cHACApK_free_st_clt(st_clt);
        free(lodfc);
        for (il = 0; il <= ndim; il++) free(gmid_t[il]);
        free(gmid_t);
        return -1;
    }

    /* Initialize leaf structures (indices 1 to nlf, 0 is unused sentinel) */
    st_leafmtx[0] = NULL;  /* Explicit sentinel */
    for (ip = 1; ip <= nlf; ip++) {
        st_leafmtx[ip] = (st_cHACApK_leafmtx)calloc(1, sizeof(st_cHACApK_leafmtx_t));
        if (!st_leafmtx[ip]) {
            fprintf(stderr, "[HACApK] Error: Memory allocation for leaf %d failed\n", ip);
            return -1;
        }
    }

    /* Generate leaf matrix structure */
    for (il = 0; il < 5; il++) lnmtx[il] = 0;
    ndpth = 0;
    nlf = 0;
    if (print_level > 0) {
        printf("[HACApK] Generating leaf matrix structure...\n");
        fflush(stdout);
    }
    cHACApK_generate_leafmtx(st_leafmtx, st_clt, st_clt, ctl->param, ctl->lpmd,
                              lnmtx, nofc, nffc, &nlf, &ndpth);

    if (print_level > 0) {
        printf("[HACApK] Generated %d leaf matrices\n", nlf);
        fflush(stdout);
    }

    /* Sort leaves for efficient traversal */
    if (print_level > 0) {
        printf("[HACApK] Sorting leaf matrices...\n");
        fflush(stdout);
    }
    cHACApK_sort_leafmtx(st_leafmtx, nlf);
    if (print_level > 0) {
        printf("[HACApK] Sort completed\n");
        fflush(stdout);
    }

    /*=========================================================================
     * Fill leaf matrices with values using ACA+
     *=========================================================================*/

    /* Compute matrix norm estimate for ACA convergence */
    znrmmat = 1.0;  /* Use 1.0 as default; could compute Frobenius norm estimate */

    /* Set thread work distribution
     * lthr[ith] = start index for thread ith (1-based leaf indices)
     * lthr[ith+1] - 1 = end index for thread ith
     * lthr[0] = 1 (first leaf)
     * lthr[nthr] = nlf + 1 (marks end)
     */
    {
        int *lnps = (int*)calloc(nthr + 1, sizeof(int));
        int *lnpe = (int*)calloc(nthr + 1, sizeof(int));
        int nlf_per_thread = nlf / nthr;
        int remainder = nlf % nthr;
        int start = 1;

        /* Initialize lthr array */
        ctl->lthr[0] = 1;  /* First thread starts at leaf 1 */

        for (il = 0; il < nthr; il++) {
            int this_thread_nlf = nlf_per_thread + (il < remainder ? 1 : 0);
            lnps[il] = start;
            lnpe[il] = start + this_thread_nlf - 1;
            start = lnpe[il] + 1;
            ctl->lthr[il + 1] = start;  /* Next thread's start = current end + 1 */
        }

        /* Fill leaf matrices using ACA+ */
        if (print_level > 0) {
            printf("[HACApK] Starting ACA+ fill: nlf=%d, kparam=%d, nd=%d\n",
                   nlf, (int)ctl->param[63], nd);
            fflush(stdout);
        }
        cHACApK_fill_leafmtx_hyp(st_leafmtx, i_bemv, ctl->param, znrmmat,
                                  ctl->lpmd, lnmtx, ctl->lod, ctl->lod, nd, nlf,
                                  lnps, lnpe, ctl->lthr);
        if (print_level > 0) {
            printf("[HACApK] ACA+ fill completed\n");
            fflush(stdout);
        }

        free(lnps);
        free(lnpe);
    }

    /*=========================================================================
     * Store results in output structures
     *=========================================================================*/
    leafmtxp->nd = nd;
    leafmtxp->nlf = nlf;
    leafmtxp->st_lf = st_leafmtx;

    /* Count low-rank blocks and compute max rank after ACA+ fill */
    leafmtxp->nlfkt = 0;
    leafmtxp->ktmax = 0;
    for (ip = 1; ip <= nlf; ip++) {
        if (st_leafmtx[ip]->ltmtx == 1) {
            leafmtxp->nlfkt++;
            if (st_leafmtx[ip]->kt > leafmtxp->ktmax) {
                leafmtxp->ktmax = st_leafmtx[ip]->kt;
            }
        }
    }

    if (print_level > 0) {
        printf("[HACApK] H-matrix built: nlf=%d, nlfkt=%d, ktmax=%d\n",
               leafmtxp->nlf, leafmtxp->nlfkt, leafmtxp->ktmax);
    }

    /* Phase 4: keep the cluster-tree root alive for downstream block-tree
     * builders (cHACApK_build_block_tree -> cHACApK_hlu_decomp). The cluster
     * tree is owned by leafmtxp from here on, freed in HACApK_free_leafmtxp. */
    leafmtxp->st_clt_root = st_clt;

    /* Cleanup temporary arrays (st_clt now owned by leafmtxp). */
    free(lodfc);
    for (il = 0; il <= ndim; il++) free(gmid_t[il]);
    free(gmid_t);

    return 0;
}

/*=========================================================================
 * Recursive helper for generating leaf matrices with variable DOF
 * Following ELF's HACApK_generate_leafmtx_varDOF pattern exactly
 *
 * Key insight: DOF indices are computed by cumulative sum in permuted order
 * - perm_dof_start_l = sum of DOFs for elements 1..(elem_start-1) in Morton order
 * - ndl = sum of DOFs for elements elem_start..elem_end in Morton order
 *=========================================================================*/
static void generate_leafmtx_varDOF_recursive(
    st_cHACApK_leafmtx *st_leafmtx,
    st_cHACApK_cluster st_cltl,
    st_cHACApK_cluster st_cltt,
    double *param,
    int *lnmtx,
    int nofc,
    int *dof_offset,
    int *lodfc,
    int *p_nlf,
    int nlf_max)
{
    int il, it, id;
    int elem_start_l, elem_end_l, elem_start_t, elem_end_t;
    int perm_dof_start_l, perm_dof_start_t;
    int ndl, ndt;
    int nleaf;
    double eta, zdistlt, zs;
    double avg_dof;
    int total_dof = dof_offset[nofc];

    /* Null check for cluster pointers */
    if (!st_cltl || !st_cltt) return;

    /* Get element ranges for clusters */
    elem_start_l = st_cltl->nstrt;
    elem_end_l = st_cltl->nstrt + st_cltl->nsize - 1;
    elem_start_t = st_cltt->nstrt;
    elem_end_t = st_cltt->nstrt + st_cltt->nsize - 1;

    /* Debug: print cluster ranges if print_level > 1 */
    if (param[1] > 1) {
        printf("[varDOF] Cluster l: nstrt=%d, nsize=%d, nnson=%d, elem_range=[%d,%d]\n",
               st_cltl->nstrt, st_cltl->nsize, st_cltl->nnson, elem_start_l, elem_end_l);
        printf("[varDOF] Cluster t: nstrt=%d, nsize=%d, nnson=%d, elem_range=[%d,%d]\n",
               st_cltt->nstrt, st_cltt->nsize, st_cltt->nnson, elem_start_t, elem_end_t);
    }

    /* Validate element ranges - handle 0-based vs 1-based */
    if (st_cltl->nsize <= 0 || st_cltt->nsize <= 0) return;
    if (elem_start_l < 1 || elem_start_l > nofc) return;
    if (elem_start_t < 1 || elem_start_t > nofc) return;
    if (elem_end_l < 1 || elem_end_l > nofc) return;
    if (elem_end_t < 1 || elem_end_t > nofc) return;

    /*=========================================================================
     * Calculate permuted DOF start index (1-based) for ROW cluster
     * ELF pattern: Sum DOFs of all elements before elem_start_l in permuted order
     *=========================================================================*/
    perm_dof_start_l = 1;
    for (il = 1; il < elem_start_l; il++) {
        int elem = lodfc[il] - 1;  /* 0-based element index */
        perm_dof_start_l += dof_offset[elem + 1] - dof_offset[elem];
    }

    /* Total DOFs in row cluster */
    ndl = 0;
    for (il = elem_start_l; il <= elem_end_l; il++) {
        int elem = lodfc[il] - 1;
        ndl += dof_offset[elem + 1] - dof_offset[elem];
    }

    /*=========================================================================
     * Calculate permuted DOF start index (1-based) for COLUMN cluster
     *=========================================================================*/
    perm_dof_start_t = 1;
    for (it = 1; it < elem_start_t; it++) {
        int elem = lodfc[it] - 1;
        perm_dof_start_t += dof_offset[elem + 1] - dof_offset[elem];
    }

    /* Total DOFs in column cluster */
    ndt = 0;
    for (it = elem_start_t; it <= elem_end_t; it++) {
        int elem = lodfc[it] - 1;
        ndt += dof_offset[elem + 1] - dof_offset[elem];
    }

    /* Skip empty blocks */
    if (ndl == 0 || ndt == 0) return;

    /* Compute average DOF and leaf size threshold */
    avg_dof = (double)total_dof / (double)nofc;
    nleaf = (int)((param[21] + 1) * avg_dof);
    if (nleaf < 3) nleaf = 3;

    /* Get admissibility parameter */
    eta = param[51];

    /* Compute distance between clusters (need valid bmin/bmax) */
    zs = 0.0;
    if (st_cltl->bmin && st_cltl->bmax && st_cltt->bmin && st_cltt->bmax) {
        for (id = 0; id < st_cltl->ndim; id++) {
            if (st_cltl->bmax[id] < st_cltt->bmin[id]) {
                double diff = st_cltt->bmin[id] - st_cltl->bmax[id];
                zs += diff * diff;
            } else if (st_cltt->bmax[id] < st_cltl->bmin[id]) {
                double diff = st_cltl->bmin[id] - st_cltt->bmax[id];
                zs += diff * diff;
            }
        }
    }
    zdistlt = sqrt(zs);

    /*=========================================================================
     * Check admissibility: clusters are well-separated if
     * min(diam_l, diam_t) <= eta * dist(l, t)
     *=========================================================================*/
    if ((st_cltl->zwdth <= eta * zdistlt || st_cltt->zwdth <= eta * zdistlt) &&
        (ndl >= nleaf && ndt >= nleaf)) {
        /* Admissible: create low-rank block */
        int nlf = *p_nlf + 1;
        if (nlf <= nlf_max) {
            st_cHACApK_leafmtx lf = st_leafmtx[nlf];
            if (lf) {
                lf->nstrtl = perm_dof_start_l;
                lf->ndl = ndl;
                lf->nstrtt = perm_dof_start_t;
                lf->ndt = ndt;
                lf->kt = 0;
                lf->ltmtx = 1;  /* Low-rank */
                *p_nlf = nlf;
                lnmtx[1]++;
            }
        }
    } else {
        /* Not admissible: check if we should create dense block or recurse */
        int nnsonl = st_cltl->nnson;
        int nnsont = st_cltt->nnson;

        if (nnsonl == 0 || nnsont == 0 || ndl <= nleaf || ndt <= nleaf) {
            /* Create dense block */
            int nlf = *p_nlf + 1;
            if (nlf <= nlf_max) {
                st_cHACApK_leafmtx lf = st_leafmtx[nlf];
                if (lf) {
                    lf->nstrtl = perm_dof_start_l;
                    lf->ndl = ndl;
                    lf->nstrtt = perm_dof_start_t;
                    lf->ndt = ndt;
                    lf->ltmtx = 2;  /* Dense */
                    *p_nlf = nlf;
                    lnmtx[2]++;
                }
            }
        } else {
            /* Recurse into children - check pc_sons is valid */
            /* NOTE: pc_sons is 1-indexed in HACApK! */
            if (st_cltl->pc_sons && st_cltt->pc_sons) {
                for (il = 1; il <= nnsonl; il++) {
                    for (it = 1; it <= nnsont; it++) {
                        generate_leafmtx_varDOF_recursive(
                            st_leafmtx,
                            st_cltl->pc_sons[il],  /* Already a pointer, not &pc_sons[il] */
                            st_cltt->pc_sons[it],
                            param, lnmtx, nofc, dof_offset, lodfc, p_nlf, nlf_max);
                    }
                }
            }
        }
    }
}

/*=========================================================================
 * Build H-matrix with variable DOF per element (for mixed hex+tetra meshes)
 * Following ELF_MAGIC's HACApK_generate_varDOF_omp pattern
 *=========================================================================*/
int HACApK_build_hmatrix_varDOF_wrapper(
    void *leafmtxp_void,
    void *ctl_void,
    double *coordinates,
    int n_elem,
    int *dof_offset,
    int total_dof,
    int ndim,
    double eps,
    int leaf_size,
    double eta,
    int print_level)
{
    st_cHACApK_leafmtxp leafmtxp = (st_cHACApK_leafmtxp)leafmtxp_void;
    st_cHACApK_lcontrol ctl = (st_cHACApK_lcontrol)ctl_void;
    st_cHACApK_cluster st_clt;
    st_cHACApK_leafmtx *st_leafmtx;
    double **gmid_t;
    int *lodfc;
    int nofc;
    int lnmtx[5];
    int ndpth, nclst, nlf, nlf_max;
    int i_bemv = 0;
    int il, ig, ip;
    double znrmmat;
    int nthr;

    nofc = n_elem;
    nthr = get_num_threads();

    if (print_level > 0) {
        printf("[HACApK] Building H-matrix (varDOF): n_elem=%d, total_dof=%d\n", n_elem, total_dof);
        printf("[HACApK] Parameters: eps=%.2e, leaf_size=%d, eta=%.2f\n", eps, leaf_size, eta);
    }

    /* Initialize leaf matrix pointer structure */
    leafmtxp->nd = total_dof;
    leafmtxp->nlf = 0;
    leafmtxp->nlfkt = 0;
    leafmtxp->ktmax = 0;
    leafmtxp->st_lf = NULL;

    /* Allocate and set control parameters */
    ctl->param = (double*)calloc(101, sizeof(double));
    ctl->lpmd = (int*)calloc(100, sizeof(int));
    ctl->time = (double*)calloc(100, sizeof(double));
    ctl->lod = (int*)calloc(total_dof + 1, sizeof(int));
    ctl->lthr = (int*)calloc(nthr + 10, sizeof(int));
    ctl->lsp = NULL;
    ctl->lnp = NULL;

    if (!ctl->param || !ctl->lpmd || !ctl->lod || !ctl->lthr) {
        fprintf(stderr, "[HACApK] Error: Memory allocation failed\n");
        return -1;
    }

    /* Set HACApK parameters */
    ctl->param[1] = (double)print_level;
    ctl->param[21] = (double)leaf_size;
    ctl->param[22] = 1.0;
    ctl->param[41] = 1.0;
    ctl->param[42] = 0;
    ctl->param[43] = 1.0;
    ctl->param[51] = eta;
    ctl->param[60] = 1;  /* ACA mode (1=ACA, 2=ACA+) */
    ctl->param[61] = 1;
    ctl->param[62] = (double)200;
    ctl->param[63] = (double)200;
    ctl->param[64] = 1;
    ctl->param[71] = eps;
    ctl->param[72] = 1.0e-3;

    /* MPI stub setup */
    ctl->lpmd[1] = MPI_COMM_WORLD;
    ctl->lpmd[2] = 1;
    ctl->lpmd[3] = 0;
    ctl->lpmd[4] = 0;
    ctl->lpmd[20] = nthr;

    /* Allocate 2D coordinate array for HACApK [ndim+1][nofc+1] */
    gmid_t = (double**)malloc(sizeof(double*) * (ndim + 1));
    for (il = 0; il <= ndim; il++) {
        gmid_t[il] = (double*)malloc(sizeof(double) * (nofc + 1));
    }

    /* Fill coordinates (1-based indexing) */
    for (il = 0; il < nofc; il++) {
        for (ig = 0; ig < ndim; ig++) {
            gmid_t[ig + 1][il + 1] = coordinates[il * ndim + ig];
        }
    }

    /* Allocate temporary element ordering array */
    lodfc = (int*)malloc(sizeof(int) * (nofc + 1));
    for (il = 1; il <= nofc; il++) {
        lodfc[il] = il;
    }

    /*=========================================================================
     * Generate cluster tree (element-based, same as uniform DOF)
     *=========================================================================*/
    ndpth = 0;
    nclst = 0;

    cHACApK_generate_cbitree(&st_clt, gmid_t, ctl->param, ctl->lpmd, lodfc,
                              &ndpth, 0, 1, nofc, nofc, ndim, &nclst);

    if (print_level > 0) {
        printf("[HACApK] Cluster tree: nclst=%d, ndpth=%d\n", nclst, ndpth);
    }

    /* Compute bounding boxes */
    cHACApK_bndbox(st_clt, gmid_t, lodfc, nofc);

    /*=========================================================================
     * Build DOF permutation array (lod)
     * For each element in Morton order, append its DOFs to lod
     *=========================================================================*/
    {
        int dof_idx = 1;  /* 1-based DOF index for lod */
        for (il = 1; il <= nofc; il++) {
            int elem = lodfc[il] - 1;  /* 0-based element index */
            int elem_dof = dof_offset[elem + 1] - dof_offset[elem];
            for (ig = 0; ig < elem_dof; ig++) {
                /* DOF global index (1-based) = dof_offset[elem] + ig + 1 */
                ctl->lod[dof_idx] = dof_offset[elem] + ig + 1;
                dof_idx++;
            }
        }
    }

    /*=========================================================================
     * Estimate maximum number of leaves
     * Conservative estimate: 4 * nofc (typical for H-matrix)
     *=========================================================================*/
    nlf_max = 4 * nofc;
    if (nlf_max < 100) nlf_max = 100;

    if (print_level > 0) {
        printf("[HACApK] Allocating up to %d leaf matrices\n", nlf_max);
    }

    /*=========================================================================
     * Allocate leaf matrix structure
     *=========================================================================*/
    st_leafmtx = (st_cHACApK_leafmtx*)calloc(nlf_max + 1, sizeof(st_cHACApK_leafmtx));
    if (!st_leafmtx) {
        fprintf(stderr, "[HACApK] Error: Memory allocation for st_leafmtx failed\n");
        cHACApK_free_st_clt(st_clt);
        free(lodfc);
        for (il = 0; il <= ndim; il++) free(gmid_t[il]);
        free(gmid_t);
        return -1;
    }

    st_leafmtx[0] = NULL;
    for (ip = 1; ip <= nlf_max; ip++) {
        st_leafmtx[ip] = (st_cHACApK_leafmtx)calloc(1, sizeof(st_cHACApK_leafmtx_t));
        if (!st_leafmtx[ip]) {
            fprintf(stderr, "[HACApK] Error: Memory allocation for leaf %d failed\n", ip);
            return -1;
        }
    }

    /*=========================================================================
     * Generate leaf matrices using recursive varDOF approach (ELF pattern)
     *=========================================================================*/
    for (il = 0; il < 5; il++) lnmtx[il] = 0;
    nlf = 0;

    generate_leafmtx_varDOF_recursive(
        st_leafmtx, st_clt, st_clt,
        ctl->param, lnmtx, nofc, dof_offset, lodfc, &nlf, nlf_max);

    if (print_level > 0) {
        printf("[HACApK] Generated %d leaf matrices (varDOF): lowrank=%d, dense=%d\n",
               nlf, lnmtx[1], lnmtx[2]);
    }

    /* Sort leaves for efficient traversal */
    cHACApK_sort_leafmtx(st_leafmtx, nlf);

    /*=========================================================================
     * Fill leaf matrices with values using ACA+
     *=========================================================================*/
    znrmmat = 1.0;

    /* Set global lod for C++ callback access during matrix fill */
    g_hacapk_lod = ctl->lod;
    g_hacapk_lod_size = total_dof;

    {
        int *lnps = (int*)calloc(nthr + 1, sizeof(int));
        int *lnpe = (int*)calloc(nthr + 1, sizeof(int));
        int nlf_per_thread = nlf / nthr;
        int remainder = nlf % nthr;
        int start = 1;

        ctl->lthr[0] = 1;

        for (il = 0; il < nthr; il++) {
            int this_thread_nlf = nlf_per_thread + (il < remainder ? 1 : 0);
            lnps[il] = start;
            lnpe[il] = start + this_thread_nlf - 1;
            start = lnpe[il] + 1;
            ctl->lthr[il + 1] = start;
        }

        cHACApK_fill_leafmtx_hyp(st_leafmtx, i_bemv, ctl->param, znrmmat,
                                  ctl->lpmd, lnmtx, ctl->lod, ctl->lod, total_dof, nlf,
                                  lnps, lnpe, ctl->lthr);

        free(lnps);
        free(lnpe);
    }

    /* Clear global lod after matrix fill */
    g_hacapk_lod = NULL;
    g_hacapk_lod_size = 0;

    /*=========================================================================
     * Store results
     *=========================================================================*/
    leafmtxp->nd = total_dof;
    leafmtxp->nlf = nlf;
    leafmtxp->st_lf = st_leafmtx;

    leafmtxp->nlfkt = 0;
    leafmtxp->ktmax = 0;
    for (ip = 1; ip <= nlf; ip++) {
        if (st_leafmtx[ip]->ltmtx == 1) {
            leafmtxp->nlfkt++;
            if (st_leafmtx[ip]->kt > leafmtxp->ktmax) {
                leafmtxp->ktmax = st_leafmtx[ip]->kt;
            }
        }
    }

    if (print_level > 0) {
        printf("[HACApK] H-matrix built (varDOF): nlf=%d, nlfkt=%d, ktmax=%d\n",
               leafmtxp->nlf, leafmtxp->nlfkt, leafmtxp->ktmax);
    }

    /* Phase 4: keep the cluster-tree root alive (see varDOF analogue above). */
    leafmtxp->st_clt_root = st_clt;

    /* Cleanup (st_clt now owned by leafmtxp). */
    free(lodfc);
    for (il = 0; il <= ndim; il++) free(gmid_t[il]);
    free(gmid_t);

    return 0;
}

/*=========================================================================
 * Matvec callback functions for TaskManager parallelization
 * These are used by hacapk_parallel_for() and hacapk_parallel_job()
 *=========================================================================*/

/* Zero one thread-local y array */
static void matvec_zero_y(int idx, void *data) {
    typedef struct { int nd; } zero_y_ctx;
    zero_y_ctx *ctx = (zero_y_ctx*)data;
    memset(g_y_thread[idx], 0, sizeof(double) * ctx->nd);
}

/* Permute one element of input vector */
static void matvec_permute_input(int idx, void *data) {
    typedef struct { const double *x; int *lod; } permute_ctx;
    permute_ctx *ctx = (permute_ctx*)data;
    /* idx is 0-based, lod is 1-based */
    g_x_perm[idx] = ctx->x[ctx->lod[idx + 1] - 1];
}

static void matvec_collect_leaf_profile(
    st_cHACApK_leafmtx *st_lf, int nlf, int symmetric,
    int64_t *lowrank_leaves, int64_t *dense_leaves,
    int64_t *mirrored_upper_leaves, int64_t *diagonal_leaves,
    int64_t *skipped_lower_leaves,
    double *lowrank_flop_est, double *dense_flop_est)
{
    int ip;
    *lowrank_leaves = 0;
    *dense_leaves = 0;
    *mirrored_upper_leaves = 0;
    *diagonal_leaves = 0;
    *skipped_lower_leaves = 0;
    *lowrank_flop_est = 0.0;
    *dense_flop_est = 0.0;
    for (ip = 1; ip <= nlf; ++ip) {
        st_cHACApK_leafmtx leaf = st_lf[ip];
        int upper = 0;
        double block_flop;
        if (!leaf) continue;
        if (symmetric) {
            if (leaf->nstrtl > leaf->nstrtt) {
                ++(*skipped_lower_leaves);
                continue;
            }
            if (leaf->nstrtl == leaf->nstrtt) ++(*diagonal_leaves);
            upper = (leaf->nstrtl < leaf->nstrtt);
        }
        if (leaf->ltmtx == 1) {
            int kt = leaf->kt;
            if (!leaf->a1 || !leaf->a2 || kt <= 0) continue;
            ++(*lowrank_leaves);
            block_flop = 2.0 * (double)kt * (double)(leaf->ndl + leaf->ndt);
            *lowrank_flop_est += block_flop;
            if (upper) {
                ++(*mirrored_upper_leaves);
                *lowrank_flop_est += block_flop;
            }
        } else {
            if (!leaf->a1) continue;
            ++(*dense_leaves);
            block_flop = 2.0 * (double)leaf->ndl * (double)leaf->ndt;
            *dense_flop_est += block_flop;
            if (upper) {
                ++(*mirrored_upper_leaves);
                *dense_flop_est += block_flop;
            }
        }
    }
}

/* Thread function for H-matrix matvec: each thread processes leaf blocks
 * using round-robin distribution and thread-local y accumulation */
static void matvec_thread_func(int tid, int nthr, void *data) {
    typedef struct {
        st_cHACApK_leafmtx *st_lf;
        int nlf;
    } matvec_job_ctx;
    matvec_job_ctx *ctx = (matvec_job_ctx*)data;
    st_cHACApK_leafmtx *st_lf = ctx->st_lf;
    int nlf = ctx->nlf;

    const double d_one = 1.0;
    const double d_zero = 0.0;
    const int i_one = 1;

    double *y_local = g_y_thread[tid];
    double *tmp_vec = g_tmp_vec[tid];

    /* Round-robin distribution of leaf blocks across threads */
    int ip;
    matvec_mkl_begin();
    for (ip = tid + 1; ip <= nlf; ip += nthr) {
        st_cHACApK_leafmtx leaf = st_lf[ip];
        if (!leaf) continue;

        int ndl = leaf->ndl;
        int ndt = leaf->ndt;
        int nstrtl = leaf->nstrtl;
        int nstrtt = leaf->nstrtt;
        double *a1 = leaf->a1;
        double *a2 = leaf->a2;

        if (leaf->ltmtx == 1) {
            /* Low-rank block: y += U * (V^T * x) */
            int kt = leaf->kt;
            if (!a1 || !a2 || kt <= 0) continue;

            dgemv_("T", &ndt, &kt, &d_one, a1, &ndt,
                   &g_x_perm[nstrtt - 1], &i_one,
                   &d_zero, tmp_vec, &i_one);

            dgemv_("N", &ndl, &kt, &d_one, a2, &ndl,
                   tmp_vec, &i_one,
                   &d_one, &y_local[nstrtl - 1], &i_one);
        } else {
            /* Dense block: y += A * x */
            if (!a1) continue;
            dgemv_("T", &ndt, &ndl, &d_one, a1, &ndt,
                   &g_x_perm[nstrtt - 1], &i_one,
                   &d_one, &y_local[nstrtl - 1], &i_one);
        }
    }
    matvec_mkl_end();
}

/* TRANSPOSE H-matrix matvec thread func: y += A^T x.  Each leaf block(l,t) (rows l = nstrtl..,
 * cols t = nstrtt..) contributes its TRANSPOSE -- i.e. x[l] -> y[t].  Mirror of matvec_thread_func
 * with the l/t roles and the dgemv "N"/"T" flags swapped.  Used to (a) probe matrix symmetry and
 * (b) form the symmetrized apply 0.5*(G x + G^T x).  Reads g_x_perm, writes g_y_thread[tid]. */
static void matvec_transpose_thread_func(int tid, int nthr, void *data) {
    typedef struct { st_cHACApK_leafmtx *st_lf; int nlf; } matvec_job_ctx;
    matvec_job_ctx *ctx = (matvec_job_ctx*)data;
    st_cHACApK_leafmtx *st_lf = ctx->st_lf;
    int nlf = ctx->nlf;
    const double d_one = 1.0;
    const double d_zero = 0.0;
    const int i_one = 1;
    double *y_local = g_y_thread[tid];
    double *tmp_vec = g_tmp_vec[tid];
    int ip;
    matvec_mkl_begin();
    for (ip = tid + 1; ip <= nlf; ip += nthr) {
        st_cHACApK_leafmtx leaf = st_lf[ip];
        if (!leaf) continue;
        int ndl = leaf->ndl, ndt = leaf->ndt;
        int nstrtl = leaf->nstrtl, nstrtt = leaf->nstrtt;
        double *a1 = leaf->a1, *a2 = leaf->a2;
        if (leaf->ltmtx == 1) {
            /* block(l,t) = a2 a1^T (a1: ndt x kt, a2: ndl x kt).  transpose = a1 a2^T: y[t] += a1 (a2^T x[l]) */
            int kt = leaf->kt;
            if (!a1 || !a2 || kt <= 0) continue;
            dgemv_("T", &ndl, &kt, &d_one, a2, &ndl,        /* tmp(kt) = a2^T x[l] */
                   &g_x_perm[nstrtl - 1], &i_one, &d_zero, tmp_vec, &i_one);
            dgemv_("N", &ndt, &kt, &d_one, a1, &ndt,        /* y[t] += a1 tmp */
                   tmp_vec, &i_one, &d_one, &y_local[nstrtt - 1], &i_one);
        } else {
            /* dense block(l,t) = a1^T (a1: ndt x ndl col-major).  transpose = a1: y[t] += a1 x[l] */
            if (!a1) continue;
            dgemv_("N", &ndt, &ndl, &d_one, a1, &ndt,
                   &g_x_perm[nstrtl - 1], &i_one, &d_one, &y_local[nstrtt - 1], &i_one);
        }
    }
    matvec_mkl_end();
}

/* SYMMETRIC H-matrix matvec thread func: y = G_sym x with G_sym EXACTLY symmetric, built from the
 * UPPER-triangular leaves only.  For a symmetric cluster tree (rows == cols, one geometry) the leaf
 * partition is symmetric, so every above-diagonal block is covered by an upper leaf (nstrtl < nstrtt);
 * that leaf supplies BOTH its own block (x[t]->y[l]) AND the mirror as its transpose (x[l]->y[t]),
 * so the lower triangle is the EXACT transpose of the upper regardless of the (independently-ACA'd)
 * lower leaves -- which makes the operator machine-symmetric and removes the ACA-asymmetry failure mode
 * for CG/MINRES.
 * Strictly-lower leaves (nstrtl > nstrtt) are SKIPPED (their pair is handled by the upper mirror);
 * diagonal leaves (nstrtl == nstrtt, the dense self/near block filled exactly from the symmetric
 * kernel) are applied once. */
static void matvec_sym_thread_func(int tid, int nthr, void *data) {
    typedef struct { st_cHACApK_leafmtx *st_lf; int nlf; } matvec_job_ctx;
    matvec_job_ctx *ctx = (matvec_job_ctx*)data;
    st_cHACApK_leafmtx *st_lf = ctx->st_lf;
    int nlf = ctx->nlf;
    const double d_one = 1.0;
    const double d_zero = 0.0;
    const int i_one = 1;
    double *y_local = g_y_thread[tid];
    double *tmp_vec = g_tmp_vec[tid];
    int ip;
    matvec_mkl_begin();
    for (ip = tid + 1; ip <= nlf; ip += nthr) {
        st_cHACApK_leafmtx leaf = st_lf[ip];
        if (!leaf) continue;
        int ndl = leaf->ndl, ndt = leaf->ndt;
        int nstrtl = leaf->nstrtl, nstrtt = leaf->nstrtt;
        double *a1 = leaf->a1, *a2 = leaf->a2;
        if (nstrtl > nstrtt) continue;          /* strictly-lower: covered by the upper mirror */
        int upper = (nstrtl < nstrtt);
        if (leaf->ltmtx == 1) {
            int kt = leaf->kt;
            if (!a1 || !a2 || kt <= 0) continue;
            /* forward x[t] -> y[l]:  tmp = a1^T x[t];  y[l] += a2 tmp */
            dgemv_("T", &ndt, &kt, &d_one, a1, &ndt,
                   &g_x_perm[nstrtt - 1], &i_one, &d_zero, tmp_vec, &i_one);
            dgemv_("N", &ndl, &kt, &d_one, a2, &ndl,
                   tmp_vec, &i_one, &d_one, &y_local[nstrtl - 1], &i_one);
            if (upper) {                         /* mirror transpose x[l] -> y[t]:  tmp = a2^T x[l]; y[t] += a1 tmp */
                dgemv_("T", &ndl, &kt, &d_one, a2, &ndl,
                       &g_x_perm[nstrtl - 1], &i_one, &d_zero, tmp_vec, &i_one);
                dgemv_("N", &ndt, &kt, &d_one, a1, &ndt,
                       tmp_vec, &i_one, &d_one, &y_local[nstrtt - 1], &i_one);
            }
        } else {
            if (!a1) continue;
            /* forward dense block(l,t)=a1^T: x[t] -> y[l] */
            dgemv_("T", &ndt, &ndl, &d_one, a1, &ndt,
                   &g_x_perm[nstrtt - 1], &i_one, &d_one, &y_local[nstrtl - 1], &i_one);
            if (upper) {                         /* mirror transpose = a1: x[l] -> y[t] */
                dgemv_("N", &ndt, &ndl, &d_one, a1, &ndt,
                       &g_x_perm[nstrtl - 1], &i_one, &d_one, &y_local[nstrtt - 1], &i_one);
            }
        }
    }
    matvec_mkl_end();
}

typedef struct {
    st_cHACApK_leafmtx *st_lf;
    int nlf, nd, nrhs;
    const int *active_prefix;
} matvec_many_job_ctx;

static void matvec_sym_many_thread_func(int tid, int nthr, void *data) {
    matvec_many_job_ctx *ctx = (matvec_many_job_ctx*)data;
    const double one = 1.0, zero = 0.0;
    const char trans = 'T', no_trans = 'N';
    double *y_local = g_batch_y_thread[tid];
    double *tmp = g_batch_tmp[tid];
    int ip;
    matvec_mkl_begin();
    for (ip = tid + 1; ip <= ctx->nlf; ip += nthr) {
        st_cHACApK_leafmtx leaf = ctx->st_lf[ip];
        int upper, ndl, ndt, nstrtl, nstrtt;
        double *a1, *a2;
        if (!leaf || leaf->nstrtl > leaf->nstrtt) continue;
        upper = leaf->nstrtl < leaf->nstrtt;
        ndl = leaf->ndl; ndt = leaf->ndt;
        nstrtl = leaf->nstrtl; nstrtt = leaf->nstrtt;
        if (ctx->active_prefix &&
            (ctx->active_prefix[nstrtl-1+ndl] ==
                 ctx->active_prefix[nstrtl-1] ||
             ctx->active_prefix[nstrtt-1+ndt] ==
                 ctx->active_prefix[nstrtt-1]))
            continue;
        a1 = leaf->a1; a2 = leaf->a2;
        if (leaf->ltmtx == 1) {
            int kt = leaf->kt;
            if (!a1 || !a2 || kt <= 0) continue;
            dgemm_(&trans, &no_trans, &kt, &ctx->nrhs, &ndt,
                   &one, a1, &ndt, &g_batch_x_perm[nstrtt-1], &ctx->nd,
                   &zero, tmp, &kt);
            dgemm_(&no_trans, &no_trans, &ndl, &ctx->nrhs, &kt,
                   &one, a2, &ndl, tmp, &kt, &one,
                   &y_local[nstrtl-1], &ctx->nd);
            if (upper) {
                dgemm_(&trans, &no_trans, &kt, &ctx->nrhs, &ndl,
                       &one, a2, &ndl, &g_batch_x_perm[nstrtl-1], &ctx->nd,
                       &zero, tmp, &kt);
                dgemm_(&no_trans, &no_trans, &ndt, &ctx->nrhs, &kt,
                       &one, a1, &ndt, tmp, &kt, &one,
                       &y_local[nstrtt-1], &ctx->nd);
            }
        } else {
            if (!a1) continue;
            dgemm_(&trans, &no_trans, &ndl, &ctx->nrhs, &ndt,
                   &one, a1, &ndt, &g_batch_x_perm[nstrtt-1], &ctx->nd,
                   &one, &y_local[nstrtl-1], &ctx->nd);
            if (upper)
                dgemm_(&no_trans, &no_trans, &ndt, &ctx->nrhs, &ndl,
                       &one, a1, &ndt, &g_batch_x_perm[nstrtl-1], &ctx->nd,
                       &one, &y_local[nstrtt-1], &ctx->nd);
        }
    }
    matvec_mkl_end();
}

static void matvec_many_zero_y(int tid, void *data) {
    matvec_many_job_ctx *ctx = (matvec_many_job_ctx*)data;
    memset(g_batch_y_thread[tid], 0,
           sizeof(double)*(size_t)ctx->nd*ctx->nrhs);
}

typedef struct {
    const double *x;
    double *y;
    int *lod;
    int nd, nrhs, nthr;
} matvec_many_io_ctx;

static void matvec_many_permute(int idx, void *data) {
    matvec_many_io_ctx *ctx = (matvec_many_io_ctx*)data;
    const int rhs = idx / ctx->nd;
    const int pos = idx - rhs*ctx->nd;
    g_batch_x_perm[pos + (size_t)ctx->nd*rhs] =
        ctx->x[(size_t)rhs*ctx->nd + ctx->lod[pos+1]-1];
}

static void matvec_many_reduce(int idx, void *data) {
    matvec_many_io_ctx *ctx = (matvec_many_io_ctx*)data;
    const int rhs = idx / ctx->nd;
    const int pos = idx - rhs*ctx->nd;
    double sum = 0.0;
    int tid;
    for (tid = 0; tid < ctx->nthr; ++tid)
        sum += g_batch_y_thread[tid][pos + (size_t)ctx->nd*rhs];
    ctx->y[(size_t)rhs*ctx->nd + ctx->lod[pos+1]-1] = sum;
}

/* Reduce thread-local y arrays and apply inverse permutation for one element */
static void matvec_reduce_output(int idx, void *data) {
    typedef struct { double *y; int *lod; int nthr; } reduce_ctx;
    reduce_ctx *ctx = (reduce_ctx*)data;
    int il = idx + 1;  /* Convert 0-based to 1-based */
    double sum = 0.0;
    int t;
    for (t = 0; t < ctx->nthr; t++) {
        sum += g_y_thread[t][idx];
    }
    ctx->y[ctx->lod[il] - 1] = sum;
}

/*=========================================================================
 * Matrix-vector product: y = A * x using H-matrix
 * Optimized version using BLAS dgemv, persistent buffers, no atomic ops
 *=========================================================================*/

void HACApK_matvec_wrapper(
    void *leafmtxp_void,
    void *ctl_void,
    const double *x,
    double *y,
    int nd)
{
    st_cHACApK_leafmtxp leafmtxp = (st_cHACApK_leafmtxp)leafmtxp_void;
    st_cHACApK_lcontrol ctl = (st_cHACApK_lcontrol)ctl_void;
    int *lod = ctl->lod;
    int nlf = leafmtxp->nlf;
    st_cHACApK_leafmtx *st_lf = leafmtxp->st_lf;
    int il, ip, t;
    int ktmax = leafmtxp->ktmax;
    int nthr = 1;

    /* BLAS constants */
    const double d_one = 1.0;
    const double d_zero = 0.0;
    const int i_one = 1;

    /* Validate pointers */
    if (!leafmtxp || !ctl || !lod || !st_lf || !x || !y) {
        return;
    }

    nthr = hacapk_get_num_threads();

    /* Initialize persistent buffers (only allocates if size changed) */
    init_matvec_buffers(nd, nthr, ktmax);

    /* Zero thread-local y arrays (parallelized via TaskManager) */
    {
        typedef struct { int nd; } zero_y_ctx;
        zero_y_ctx zctx = { nd };
        hacapk_parallel_for(nthr, matvec_zero_y, &zctx);
    }

    /* Pre-permute input vector (parallelized via TaskManager) */
    {
        typedef struct { const double *x; int *lod; } permute_ctx;
        permute_ctx pctx = { x, lod };
        hacapk_parallel_for(nd, matvec_permute_input, &pctx);
    }

    /* H-matrix matrix-vector product using TaskManager
     * Each thread processes leaf blocks using round-robin distribution
     * and accumulates into thread-local y arrays for lock-free operation
     */
    {
        typedef struct {
            st_cHACApK_leafmtx *st_lf;
            int nlf;
        } matvec_job_ctx;
        matvec_job_ctx mctx = { st_lf, nlf };
        hacapk_parallel_job(matvec_thread_func, &mctx);
    }

    /* Reduce thread-local y arrays and apply inverse permutation */
    {
        typedef struct { double *y; int *lod; int nthr; } reduce_ctx;
        reduce_ctx rctx = { y, lod, nthr };
        hacapk_parallel_for(nd, matvec_reduce_output, &rctx);
    }
}

/* Shared driver for the transpose / symmetric matvec variants: same buffer setup, input permute,
 * and output reduce as HACApK_matvec_wrapper, but with a caller-supplied leaf thread function. */
static void hacapk_matvec_run(
    void *leafmtxp_void, void *ctl_void, const double *x, double *y, int nd,
    void (*thread_func)(int, int, void *), int symmetric_profile)
{
    st_cHACApK_leafmtxp leafmtxp = (st_cHACApK_leafmtxp)leafmtxp_void;
    st_cHACApK_lcontrol ctl = (st_cHACApK_lcontrol)ctl_void;
    if (!leafmtxp || !ctl || !ctl->lod || !leafmtxp->st_lf || !x || !y) return;
    int *lod = ctl->lod;
    int nlf = leafmtxp->nlf;
    st_cHACApK_leafmtx *st_lf = leafmtxp->st_lf;
    int ktmax = leafmtxp->ktmax;
    int nthr = hacapk_get_num_threads();
    int stats_on = matvec_stats_enabled();
    double t_total0 = 0.0, t0 = 0.0;
    double zero_s = 0.0, permute_s = 0.0, leaf_s = 0.0, reduce_s = 0.0, meta_s = 0.0;
    int64_t lowrank_leaves = 0, dense_leaves = 0, mirrored_upper_leaves = 0;
    int64_t diagonal_leaves = 0, skipped_lower_leaves = 0;
    double lowrank_flop_est = 0.0, dense_flop_est = 0.0;

    if (stats_on) {
        t_total0 = matvec_wall_time();
        t0 = t_total0;
        matvec_collect_leaf_profile(st_lf, nlf, symmetric_profile,
                                    &lowrank_leaves, &dense_leaves,
                                    &mirrored_upper_leaves, &diagonal_leaves,
                                    &skipped_lower_leaves,
                                    &lowrank_flop_est, &dense_flop_est);
        meta_s = matvec_wall_time() - t0;
    }

    init_matvec_buffers(nd, nthr, ktmax);
    if (stats_on) t0 = matvec_wall_time();
    { typedef struct { int nd; } zero_y_ctx; zero_y_ctx zc = { nd };
      hacapk_parallel_for(nthr, matvec_zero_y, &zc); }
    if (stats_on) { zero_s = matvec_wall_time() - t0; t0 = matvec_wall_time(); }
    { typedef struct { const double *x; int *lod; } permute_ctx; permute_ctx pc = { x, lod };
      hacapk_parallel_for(nd, matvec_permute_input, &pc); }
    if (stats_on) { permute_s = matvec_wall_time() - t0; t0 = matvec_wall_time(); }
    { typedef struct { st_cHACApK_leafmtx *st_lf; int nlf; } matvec_job_ctx;
      matvec_job_ctx mc = { st_lf, nlf };
      hacapk_parallel_job(thread_func, &mc); }
    if (stats_on) { leaf_s = matvec_wall_time() - t0; t0 = matvec_wall_time(); }
    { typedef struct { double *y; int *lod; int nthr; } reduce_ctx; reduce_ctx rc = { y, lod, nthr };
      hacapk_parallel_for(nd, matvec_reduce_output, &rc); }
    if (stats_on) {
        reduce_s = matvec_wall_time() - t0;
        g_matvec_stats.calls += 1;
        g_matvec_stats.total_s += matvec_wall_time() - t_total0;
        g_matvec_stats.zero_s += zero_s;
        g_matvec_stats.permute_s += permute_s;
        g_matvec_stats.leaf_s += leaf_s;
        g_matvec_stats.reduce_s += reduce_s;
        g_matvec_stats.meta_s += meta_s;
        g_matvec_stats.lowrank_leaves += lowrank_leaves;
        g_matvec_stats.dense_leaves += dense_leaves;
        g_matvec_stats.mirrored_upper_leaves += mirrored_upper_leaves;
        g_matvec_stats.diagonal_leaves += diagonal_leaves;
        g_matvec_stats.skipped_lower_leaves += skipped_lower_leaves;
        g_matvec_stats.lowrank_flop_est += lowrank_flop_est;
        g_matvec_stats.dense_flop_est += dense_flop_est;
        g_matvec_stats.last_nd = nd;
        g_matvec_stats.last_nthr = nthr;
    }
}

/* y = A^T x (transpose H-matvec). */
void HACApK_matvec_transpose_wrapper(
    void *leafmtxp_void, void *ctl_void, const double *x, double *y, int nd)
{
    hacapk_matvec_run(leafmtxp_void, ctl_void, x, y, nd, matvec_transpose_thread_func, 0);
}

/* y = G_sym x with G_sym EXACTLY symmetric (upper-triangular leaves define both triangles). */
void HACApK_matvec_sym_wrapper(
    void *leafmtxp_void, void *ctl_void, const double *x, double *y, int nd)
{
    hacapk_matvec_run(leafmtxp_void, ctl_void, x, y, nd, matvec_sym_thread_func, 1);
}

static void hacapk_matvec_sym_many_run(
    void *leafmtxp_void, void *ctl_void, const double *x, double *y,
    int nd, int nrhs, const int *active_prefix)
{
    st_cHACApK_leafmtxp leafmtxp = (st_cHACApK_leafmtxp)leafmtxp_void;
    st_cHACApK_lcontrol ctl = (st_cHACApK_lcontrol)ctl_void;
    int nthr;
    matvec_many_job_ctx job;
    matvec_many_io_ctx io;
    if (!leafmtxp || !ctl || !ctl->lod || !leafmtxp->st_lf ||
        !x || !y || nd < 1 || nrhs < 1) return;
    nthr = hacapk_get_num_threads();
    init_batch_matvec_buffers(nd, nrhs, nthr, leafmtxp->ktmax);
    job.st_lf = leafmtxp->st_lf; job.nlf = leafmtxp->nlf;
    job.nd = nd; job.nrhs = nrhs;
    job.active_prefix = active_prefix;
    io.x = x; io.y = y; io.lod = ctl->lod; io.nd = nd;
    io.nrhs = nrhs; io.nthr = nthr;
    hacapk_parallel_for(nthr, matvec_many_zero_y, &job);
    hacapk_parallel_for(nd*nrhs, matvec_many_permute, &io);
    hacapk_parallel_job(matvec_sym_many_thread_func, &job);
    hacapk_parallel_for(nd*nrhs, matvec_many_reduce, &io);
}

void HACApK_matvec_sym_many_wrapper(
    void *leafmtxp_void, void *ctl_void, const double *x, double *y,
    int nd, int nrhs)
{
    hacapk_matvec_sym_many_run(
        leafmtxp_void, ctl_void, x, y, nd, nrhs, NULL);
}

void HACApK_matvec_sym_many_masked_wrapper(
    void *leafmtxp_void, void *ctl_void, const double *x, double *y,
    int nd, int nrhs, const int *active_prefix)
{
    if (!active_prefix || active_prefix[0] != 0) return;
    hacapk_matvec_sym_many_run(
        leafmtxp_void, ctl_void, x, y, nd, nrhs, active_prefix);
}

/*=========================================================================
 * Free H-matrix resources
 *=========================================================================*/

void HACApK_free_hmatrix_wrapper(
    void *leafmtxp_void,
    void *ctl_void)
{
    st_cHACApK_leafmtxp leafmtxp = (st_cHACApK_leafmtxp)leafmtxp_void;
    st_cHACApK_lcontrol ctl = (st_cHACApK_lcontrol)ctl_void;
    int ip;

    if (leafmtxp && leafmtxp->st_lf) {
        for (ip = 1; ip <= leafmtxp->nlf; ip++) {
            if (leafmtxp->st_lf[ip]) {
                if (leafmtxp->st_lf[ip]->a1) free(leafmtxp->st_lf[ip]->a1);
                if (leafmtxp->st_lf[ip]->a2) free(leafmtxp->st_lf[ip]->a2);
                free(leafmtxp->st_lf[ip]);
            }
        }
        free(leafmtxp->st_lf);
        leafmtxp->st_lf = NULL;
    }

    if (ctl) {
        if (ctl->lod) { free(ctl->lod); ctl->lod = NULL; }
        if (ctl->lsp) { free(ctl->lsp); ctl->lsp = NULL; }
        if (ctl->lnp) { free(ctl->lnp); ctl->lnp = NULL; }
        if (ctl->lthr) { free(ctl->lthr); ctl->lthr = NULL; }
        if (ctl->lpmd) { free(ctl->lpmd); ctl->lpmd = NULL; }
        if (ctl->param) { free(ctl->param); ctl->param = NULL; }
        if (ctl->time) { free(ctl->time); ctl->time = NULL; }
    }
}

/*=========================================================================
 * Allocate structures (called from C++)
 *=========================================================================*/

void* HACApK_alloc_leafmtxp(void) {
    return calloc(1, sizeof(st_cHACApK_leafmtxp_t));
}

void* HACApK_alloc_lcontrol(void) {
    return calloc(1, sizeof(st_cHACApK_lcontrol_t));
}

void HACApK_free_leafmtxp(void *ptr) {
    if (!ptr) return;
    st_cHACApK_leafmtxp lp = (st_cHACApK_leafmtxp)ptr;
    /* Phase 4: free the preserved cluster-tree root if any. */
    if (lp->st_clt_root) {
        cHACApK_free_st_clt(lp->st_clt_root);
        lp->st_clt_root = NULL;
    }
    free(ptr);
}

void HACApK_free_lcontrol(void *ptr) {
    if (ptr) free(ptr);
}

/*=========================================================================
 * Accessor functions for C++ to read struct fields
 *=========================================================================*/

int HACApK_leafmtxp_get_nd(void *ptr) {
    st_cHACApK_leafmtxp lp = (st_cHACApK_leafmtxp)ptr;
    return lp ? lp->nd : 0;
}

int HACApK_leafmtxp_get_nlf(void *ptr) {
    st_cHACApK_leafmtxp lp = (st_cHACApK_leafmtxp)ptr;
    return lp ? lp->nlf : 0;
}

int HACApK_leafmtxp_get_nlfkt(void *ptr) {
    st_cHACApK_leafmtxp lp = (st_cHACApK_leafmtxp)ptr;
    return lp ? lp->nlfkt : 0;
}

int HACApK_leafmtxp_get_ktmax(void *ptr) {
    st_cHACApK_leafmtxp lp = (st_cHACApK_leafmtxp)ptr;
    return lp ? lp->ktmax : 0;
}

int* HACApK_lcontrol_get_lod(void *ptr) {
    st_cHACApK_lcontrol ctl = (st_cHACApK_lcontrol)ptr;
    return ctl ? ctl->lod : NULL;
}

/*=========================================================================
 * Calculate actual H-matrix memory usage (ELF-compatible)
 *
 * Returns memory in bytes by iterating over all leaf blocks:
 * - Low-rank (ltmtx=1): (ndl + ndt) * kt * sizeof(double)
 * - Dense (ltmtx=2):    ndl * ndt * sizeof(double)
 *
 * This matches ELF's HACApK_get_stats implementation exactly.
 *=========================================================================*/
void HACApK_get_memory_stats(void *leafmtxp_void,
                              int64_t *hmat_bytes_out,
                              int64_t *dense_bytes_out) {
    st_cHACApK_leafmtxp lp = (st_cHACApK_leafmtxp)leafmtxp_void;
    int64_t hmat_bytes = 0;
    int64_t dense_bytes = 0;
    int ip;

    if (!lp || !lp->st_lf) {
        if (hmat_bytes_out) *hmat_bytes_out = 0;
        if (dense_bytes_out) *dense_bytes_out = 0;
        return;
    }

    /* Iterate over all leaf blocks (1-indexed) */
    for (ip = 1; ip <= lp->nlf; ip++) {
        st_cHACApK_leafmtx leaf = lp->st_lf[ip];
        if (!leaf) continue;

        int ndl = leaf->ndl;
        int ndt = leaf->ndt;

        if (leaf->ltmtx == 1) {
            /* Low-rank block: stores U (ndl x kt) and V (ndt x kt) */
            int kt = leaf->kt;
            hmat_bytes += (int64_t)(ndl + ndt) * kt * sizeof(double);
        } else {
            /* Dense block: stores full ndl x ndt matrix */
            hmat_bytes += (int64_t)ndl * ndt * sizeof(double);
        }

        /* Equivalent dense for this block */
        dense_bytes += (int64_t)ndl * ndt * sizeof(double);
    }

    if (hmat_bytes_out) *hmat_bytes_out = hmat_bytes;
    if (dense_bytes_out) *dense_bytes_out = dense_bytes;
}

/*=========================================================================
 * Update dense diagonal blocks for nonlinear iteration
 *
 * This function follows ELF_MAGIC's HACApK_update_diagonal_omp pattern:
 * - Only dense diagonal blocks (ltmtx==2 && nstrtl==nstrtt) are recomputed
 * - Low-rank off-diagonal blocks remain unchanged (geometry doesn't change)
 * - Entry function callback is called for each element of diagonal blocks
 *
 * This is essential for nonlinear material iteration where chi changes
 * but the geometry (and thus off-diagonal low-rank approximations) stays fixed.
 *=========================================================================*/

void HACApK_update_diagonal_wrapper(
    void *leafmtxp_void,
    void *ctl_void,
    double (*entry_func)(int i, int j, int i_bemv))
{
    st_cHACApK_leafmtxp leafmtxp = (st_cHACApK_leafmtxp)leafmtxp_void;
    st_cHACApK_lcontrol ctl = (st_cHACApK_lcontrol)ctl_void;
    int ip, il, it;
    int i_bemv = 0;
    int *lod;
    int n_diag_updated = 0;

    if (!leafmtxp || !ctl || !entry_func) return;

    lod = ctl->lod;
    if (!lod) return;

    /* Iterate over all leaf blocks (parallelized via TaskManager) */
    /* Note: n_diag_updated is informational only, exact count not critical */
    for (ip = 1; ip <= leafmtxp->nlf; ip++) {
        st_cHACApK_leafmtx leaf = leafmtxp->st_lf[ip];

        if (!leaf) continue;

        /* Check if this is a dense diagonal block */
        /* ltmtx==2: dense block, nstrtl==nstrtt: diagonal position */
        if (leaf->ltmtx == 2 && leaf->nstrtl == leaf->nstrtt) {
            int ndl = leaf->ndl;
            int ndt = leaf->ndt;
            int nstrtl = leaf->nstrtl;
            int nstrtt = leaf->nstrtt;
            double *a1 = leaf->a1;

            /* Recompute all entries in this dense diagonal block */
            /* Storage: a1[it + ndt * il] = A[il, it] (row-major) */
            for (il = 0; il < ndl; il++) {
                int ill = lod[nstrtl + il];  /* Global row index (1-based) */
                for (it = 0; it < ndt; it++) {
                    int itt = lod[nstrtt + it];  /* Global column index (1-based) */
                    a1[it + ndt * il] = entry_func(ill, itt, i_bemv);
                }
            }
            n_diag_updated++;
        }
    }
}

/*=========================================================================
 * Fast diagonal update: only update true diagonal elements (i==j)
 *
 * This is MUCH faster than HACApK_update_diagonal_wrapper because:
 * - Only iterates over diagonal entries in diagonal blocks
 * - Uses pre-computed N_ii values instead of calling entry_func
 * - O(ndof) instead of O(block_size^2 * n_diag_blocks)
 *
 * For a 1000-element problem (6000 DOF), this reduces from:
 * - Old: ~180 diagonal blocks * 32*32 entries * expensive entry_func call
 * - New: 6000 diagonal entries * simple array lookup
 *=========================================================================*/

void HACApK_update_diagonal_fast_wrapper(
    void *leafmtxp_void,
    void *ctl_void,
    const double *diag_N,
    const double *inv_chi,
    int ndof)
{
    st_cHACApK_leafmtxp leafmtxp = (st_cHACApK_leafmtxp)leafmtxp_void;
    st_cHACApK_lcontrol ctl = (st_cHACApK_lcontrol)ctl_void;
    int ip, il;
    int *lod;

    if (!leafmtxp || !ctl || !diag_N || !inv_chi) return;

    lod = ctl->lod;
    if (!lod) return;

    /* Iterate over all leaf blocks */
    for (ip = 1; ip <= leafmtxp->nlf; ip++) {
        st_cHACApK_leafmtx leaf = leafmtxp->st_lf[ip];

        if (!leaf) continue;

        /* Only process dense diagonal blocks (ltmtx==2 && nstrtl==nstrtt) */
        if (leaf->ltmtx == 2 && leaf->nstrtl == leaf->nstrtt) {
            int ndl = leaf->ndl;
            int ndt = leaf->ndt;
            int nstrtl = leaf->nstrtl;
            double *a1 = leaf->a1;

            /* Only update true diagonal entries (il == it) */
            /* This is the key optimization: O(ndl) instead of O(ndl * ndt) */
            for (il = 0; il < ndl && il < ndt; il++) {
                int global_idx = lod[nstrtl + il] - 1;  /* Convert 1-based to 0-based */

                if (global_idx >= 0 && global_idx < ndof) {
                    /* A_ii = -K/(4pi) + 1/chi = -diag_N + inv_chi */
                    /* diag_N contains K/(4pi) (positive), we need -K/(4pi) */
                    /* Storage: a1[it + ndt * il] where it==il for diagonal */
                    a1[il + ndt * il] = -diag_N[global_idx] + inv_chi[global_idx];
                }
            }
        }
    }
}
