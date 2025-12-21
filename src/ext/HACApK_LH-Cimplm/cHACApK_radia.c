/*
 * HACApK for Radia - Simplified H-matrix with ACA+ for single-process execution
 *
 * Based on HACApK (ppOpen-HPC project, MIT License)
 * Copyright (c) 2015 Akihiro Ida and Takeshi Iwashita
 *
 * Simplified for Radia by removing MPI support and keeping OpenMP only.
 * This file provides the simplified C interface defined in cHACApK_radia.h.
 *
 * Copyright (c) 2025 Radia Project
 */

#include "cHACApK_radia.h"
#include "cHACApK_base.h"
#include "cHACApK_calc_entry_ij.h"
#include "hacapk_log.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _WIN32
#include <time.h>
#else
#include <sys/time.h>
#endif

/*
 * Internal structures (match those in cHACApK_base.h)
 */
struct st_cHACApK_cluster_impl {
    int ndim;
    int nstrt, nsize, ndpth, nnson, nmbr;
    int ndscd;
    double *bmin;
    double *bmax;
    double zwdth;
    struct st_cHACApK_cluster_impl *pc_sons;
};

struct st_cHACApK_leafmtx_impl {
    int ltmtx;  /* 1:low-rank 2:dense */
    int kt;     /* rank */
    int nstrtl, ndl;  /* row start and size */
    int nstrtt, ndt;  /* column start and size */
    double *a1, *a2;  /* low-rank factors or dense matrix */

    int nlf;    /* number of sub-leaves */
    struct st_cHACApK_leafmtx_impl *st_lf;
};

struct st_cHACApK_leafmtxp_impl {
    int nd;     /* total DOF */
    int nlf;    /* number of leaves */
    int nlfkt;  /* number of low-rank leaves */
    int ktmax;  /* maximum rank */
    int nbl;    /* number of blocks */
    int nlfalt;
    int nlfl, nlft;
    int ndlfs, ndtfs;
    struct st_cHACApK_leafmtx_impl *st_lf;
    int64_t **lnlfl2g_t;
    int *lbstrtl, *lbstrtt;
    int *lbndl, *lbndt;
    int *lbndlfs, *lbndtfs;
    int *lbl2t;
};

struct st_cHACApK_lcontrol_impl {
    int *lod;       /* DOF permutation */
    int *lsp;
    int *lnp;
    int *lthr;      /* thread work distribution */
    int *lpmd;      /* control parameters */
    double *param;  /* numerical parameters */
    double *time;   /* timing information */
    int lf_umpi;
};

/*
 * Global state for matrix entry callback
 */
static HACApK_entry_func g_entry_func = NULL;
static int g_bemv = 0;

/*
 * Get current time in seconds
 */
static double get_wtime(void) {
#ifdef _WIN32
    return (double)clock() / (double)CLOCKS_PER_SEC;
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
#endif
}

/*
 * Initialize default configuration
 */
void HACApK_init_config(HACApK_config *config) {
    if (config == NULL) return;

    config->aca_eps = 1.0e-4;
    config->leaf_size = 32;
    config->eta = 2.0;
    config->max_rank = 200;
    config->print_level = 0;  /* Silent by default for Radia */
    config->num_threads = 0;  /* Auto-detect */
}

/*
 * Allocate H-matrix structures
 */
st_cHACApK_leafmtxp HACApK_alloc_leafmtxp(void) {
    struct st_cHACApK_leafmtxp_impl *ptr =
        (struct st_cHACApK_leafmtxp_impl *)calloc(1, sizeof(struct st_cHACApK_leafmtxp_impl));
    return ptr;
}

st_cHACApK_lcontrol HACApK_alloc_lcontrol(void) {
    struct st_cHACApK_lcontrol_impl *ptr =
        (struct st_cHACApK_lcontrol_impl *)calloc(1, sizeof(struct st_cHACApK_lcontrol_impl));
    return ptr;
}

/*
 * Free H-matrix structures
 */
void HACApK_free_leafmtxp(st_cHACApK_leafmtxp ptr) {
    struct st_cHACApK_leafmtxp_impl *p = (struct st_cHACApK_leafmtxp_impl *)ptr;
    if (p == NULL) return;

    /* Free leaf matrices */
    if (p->st_lf != NULL) {
        for (int i = 1; i <= p->nlf; i++) {
            if (p->st_lf[i] != NULL) {
                if (p->st_lf[i]->a1 != NULL) free(p->st_lf[i]->a1);
                if (p->st_lf[i]->a2 != NULL) free(p->st_lf[i]->a2);
                free(p->st_lf[i]);
            }
        }
        free(p->st_lf);
    }

    /* Free auxiliary arrays */
    if (p->lnlfl2g_t != NULL) {
        for (int i = 1; i <= p->nlfl; i++) {
            if (p->lnlfl2g_t[i] != NULL) free(p->lnlfl2g_t[i]);
        }
        free(p->lnlfl2g_t);
    }
    if (p->lbstrtl != NULL) free(p->lbstrtl);
    if (p->lbstrtt != NULL) free(p->lbstrtt);
    if (p->lbndl != NULL) free(p->lbndl);
    if (p->lbndt != NULL) free(p->lbndt);
    if (p->lbndlfs != NULL) free(p->lbndlfs);
    if (p->lbndtfs != NULL) free(p->lbndtfs);
    if (p->lbl2t != NULL) free(p->lbl2t);

    free(p);
}

void HACApK_free_lcontrol(st_cHACApK_lcontrol ptr) {
    struct st_cHACApK_lcontrol_impl *p = (struct st_cHACApK_lcontrol_impl *)ptr;
    if (p == NULL) return;

    if (p->lod != NULL) free(p->lod);
    if (p->lsp != NULL) free(p->lsp);
    if (p->lnp != NULL) free(p->lnp);
    if (p->lthr != NULL) free(p->lthr);
    if (p->lpmd != NULL) free(p->lpmd);
    if (p->param != NULL) free(p->param);
    if (p->time != NULL) free(p->time);

    free(p);
}

/*
 * Callback bridge for cHACApK_entry_ij
 * This function is called by the HACApK library to get matrix elements
 */
double cHACApK_entry_ij(int i, int j, int bemv) {
    if (g_entry_func != NULL) {
        return g_entry_func(i, j, bemv);
    }
    return 0.0;
}

/*
 * Build H-matrix from coordinates and entry callback
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
    HACApK_stats *stats)
{
    struct st_cHACApK_leafmtxp_impl *p_leafmtxp = (struct st_cHACApK_leafmtxp_impl *)leafmtxp;
    struct st_cHACApK_lcontrol_impl *p_ctl = (struct st_cHACApK_lcontrol_impl *)ctl;

    if (p_leafmtxp == NULL || p_ctl == NULL || coordinates == NULL || entry_func == NULL) {
        HACAPK_LOG_ERROR_MSG("Invalid arguments to HACApK_build");
        return -1;
    }

    double t_start = get_wtime();

    /* Set global callback */
    g_entry_func = entry_func;
    g_bemv = bemv;

    /* Determine number of threads */
    int nthr = config->num_threads;
    if (nthr <= 0) {
#ifdef _OPENMP
        nthr = omp_get_max_threads();
#else
        nthr = 1;
#endif
    }

    /* Total DOF */
    int nd = n_elem * nffc;
    p_leafmtxp->nd = nd;

    HACAPK_LOG_INFO_MSG("Building H-matrix: n_elem=%d, nffc=%d, nd=%d, nthr=%d",
                        n_elem, nffc, nd, nthr);

    /* Allocate control arrays */
    /* Parameter array (100 elements, 1-indexed in original code) */
    p_ctl->param = (double *)calloc(100, sizeof(double));
    p_ctl->lpmd = (int *)calloc(100, sizeof(int));
    p_ctl->lod = (int *)calloc(nd + 1, sizeof(int));
    p_ctl->lthr = (int *)calloc(nthr + 2, sizeof(int));
    p_ctl->time = (double *)calloc(20, sizeof(double));

    if (p_ctl->param == NULL || p_ctl->lpmd == NULL ||
        p_ctl->lod == NULL || p_ctl->lthr == NULL) {
        HACAPK_LOG_ERROR_MSG("Memory allocation failed for control structures");
        return -2;
    }

    /* Set parameters */
    /* param[1]: print level */
    p_ctl->param[1] = config->print_level;

    /* param[21]: minimum leaf size */
    p_ctl->param[21] = config->leaf_size;

    /* param[51]: admissibility parameter eta */
    p_ctl->param[51] = config->eta;

    /* param[61]: ACA mode (1=relative, 2=absolute, 3=scaled) */
    p_ctl->param[61] = 1;

    /* param[62]: expected rank (for memory allocation) */
    p_ctl->param[62] = config->max_rank;

    /* param[63]: maximum rank */
    p_ctl->param[63] = config->max_rank;

    /* param[64]: minimum rank */
    p_ctl->param[64] = 1;

    /* param[71]: ACA tolerance */
    p_ctl->param[71] = config->aca_eps;

    /* param[72]: ACA_EPS multiplier */
    p_ctl->param[72] = 1.0;

    /* param[41], [42], [43]: block clustering parameters */
    p_ctl->param[41] = 1;  /* npgl = 1 for single process */
    p_ctl->param[42] = 0;  /* auto block size */
    p_ctl->param[43] = 1;

    /* lpmd: MPI parameters (all single-process values) */
    p_ctl->lpmd[1] = 0;    /* icomm (communicator handle) */
    p_ctl->lpmd[2] = 1;    /* nrank (number of processes) */
    p_ctl->lpmd[3] = 0;    /* mpinr (this process rank) */
    p_ctl->lpmd[4] = 0;    /* mpilog (log level) */
    p_ctl->lpmd[20] = nthr; /* number of threads */

    /* Initialize DOF ordering */
    for (int i = 1; i <= nd; i++) {
        p_ctl->lod[i] = i;
    }

    /* Create coordinate array in HACApK format */
    /* gmid_t[dim][elem] - transposed coordinate array */
    double **gmid_t = (double **)malloc((ndim + 1) * sizeof(double *));
    if (gmid_t == NULL) {
        HACAPK_LOG_ERROR_MSG("Memory allocation failed for coordinate array");
        return -3;
    }

    for (int d = 0; d <= ndim; d++) {
        gmid_t[d] = (double *)malloc((n_elem + 1) * sizeof(double));
        if (gmid_t[d] == NULL) {
            HACAPK_LOG_ERROR_MSG("Memory allocation failed for coordinate array dimension %d", d);
            for (int k = 0; k < d; k++) free(gmid_t[k]);
            free(gmid_t);
            return -3;
        }
    }

    /* Copy coordinates (1-indexed) */
    for (int i = 0; i < n_elem; i++) {
        for (int d = 0; d < ndim; d++) {
            gmid_t[d + 1][i + 1] = coordinates[i * ndim + d];
        }
    }

    /* Build H-matrix using cHACApK functions */
    int lnmtx[5] = {0, 0, 0, 0, 0};

    /* Call the frame generation function */
    /* This builds the cluster tree and leaf matrix structure */
    cHACApK_generate_frame_blrleaf(
        (st_cHACApK_leafmtxp)p_leafmtxp,
        bemv,
        (st_cHACApK_lcontrol)p_ctl,
        gmid_t,
        lnmtx,
        n_elem,
        nffc,
        ndim);

    /* Fill leaf matrices with ACA+ compression */
    int nlf = p_leafmtxp->nlf;
    double znrmmat = 1.0;  /* Matrix norm estimate */

    if (nlf > 0) {
        int *lnps = (int *)malloc((nthr + 1) * sizeof(int));
        int *lnpe = (int *)malloc((nthr + 1) * sizeof(int));

        if (lnps != NULL && lnpe != NULL) {
            /* Simple work distribution */
            int nlf_per_thr = nlf / nthr;
            for (int i = 0; i < nthr; i++) {
                lnps[i] = i * nlf_per_thr + 1;
                lnpe[i] = (i == nthr - 1) ? nlf : (i + 1) * nlf_per_thr;
            }

            cHACApK_fill_leafmtx_hyp(
                p_leafmtxp->st_lf,
                bemv,
                p_ctl->param,
                znrmmat,
                p_ctl->lpmd,
                lnmtx,
                p_ctl->lod,
                p_ctl->lod,
                nd,
                nlf,
                lnps,
                lnpe,
                p_ctl->lthr);

            free(lnps);
            free(lnpe);
        }
    }

    /* Free coordinate array */
    for (int d = 0; d <= ndim; d++) {
        free(gmid_t[d]);
    }
    free(gmid_t);

    double t_end = get_wtime();

    /* Compute statistics */
    if (stats != NULL) {
        stats->n_lowrank = lnmtx[1];
        stats->n_dense = lnmtx[2];
        stats->n_leaves = nlf;
        stats->n_dof = nd;
        stats->max_rank = p_leafmtxp->ktmax;
        stats->build_time = t_end - t_start;

        /* Compute compression ratio */
        int64_t full_size = (int64_t)nd * nd;
        int64_t compressed_size = 0;
        for (int i = 1; i <= nlf; i++) {
            struct st_cHACApK_leafmtx_impl *lf = p_leafmtxp->st_lf[i];
            if (lf->ltmtx == 1) {
                /* Low-rank: (ndl + ndt) * kt */
                compressed_size += (int64_t)(lf->ndl + lf->ndt) * lf->kt;
            } else {
                /* Dense: ndl * ndt */
                compressed_size += (int64_t)lf->ndl * lf->ndt;
            }
        }
        stats->compression = (double)compressed_size / (double)full_size;
    }

    HACAPK_LOG_INFO_MSG("H-matrix built: %d leaves (%d low-rank, %d dense), time=%.3f sec",
                        nlf, lnmtx[1], lnmtx[2], t_end - t_start);

    return 0;
}

/*
 * Matrix-vector product: y = A * x
 */
void HACApK_matvec(
    st_cHACApK_leafmtxp leafmtxp,
    st_cHACApK_lcontrol ctl,
    const double *x,
    double *y,
    int n_dof)
{
    struct st_cHACApK_leafmtxp_impl *p = (struct st_cHACApK_leafmtxp_impl *)leafmtxp;
    struct st_cHACApK_lcontrol_impl *c = (struct st_cHACApK_lcontrol_impl *)ctl;

    if (p == NULL || c == NULL || x == NULL || y == NULL) return;

    /* Initialize output */
    memset(y, 0, n_dof * sizeof(double));

    /* Apply permutation to input */
    double *xp = (double *)malloc(n_dof * sizeof(double));
    if (xp == NULL) return;

    for (int i = 1; i <= n_dof; i++) {
        xp[i - 1] = x[c->lod[i] - 1];
    }

    /* Temporary output in permuted order */
    double *yp = (double *)calloc(n_dof, sizeof(double));
    if (yp == NULL) {
        free(xp);
        return;
    }

    /* Matrix-vector multiplication over all leaves */
    int nlf = p->nlf;

#pragma omp parallel for schedule(dynamic) if(nlf > 10)
    for (int ip = 1; ip <= nlf; ip++) {
        struct st_cHACApK_leafmtx_impl *lf = p->st_lf[ip];
        int ltmtx = lf->ltmtx;
        int ndl = lf->ndl;
        int ndt = lf->ndt;
        int nstrtl = lf->nstrtl;
        int nstrtt = lf->nstrtt;

        if (ltmtx == 1) {
            /* Low-rank matrix: y += U * (V^T * x)
             * Storage from cHACApK_base.c:
             *   a1 = zab (ndt x kt) = V matrix (column factor)
             *   a2 = zaa (ndl x kt) = U matrix (row factor)
             * Low-rank: A ≈ U * V^T where U is ndl x kt, V is ndt x kt
             */
            int kt = lf->kt;
            double *V = lf->a1;  /* V stored as ndt x kt, column-major */
            double *U = lf->a2;  /* U stored as ndl x kt, column-major */

            /* Compute z = V^T * x (size kt) */
            double *z = (double *)calloc(kt, sizeof(double));
            if (z == NULL) continue;

            for (int k = 0; k < kt; k++) {
                double sum = 0.0;
                for (int j = 0; j < ndt; j++) {
                    sum += V[j + ndt * k] * xp[nstrtt - 1 + j];
                }
                z[k] = sum;
            }

            /* Compute y += U * z */
#pragma omp critical
            {
                for (int i = 0; i < ndl; i++) {
                    double sum = 0.0;
                    for (int k = 0; k < kt; k++) {
                        sum += U[i + ndl * k] * z[k];
                    }
                    yp[nstrtl - 1 + i] += sum;
                }
            }

            free(z);
        } else {
            /* Dense matrix: y += A * x
             * Storage from cHACApK_base.c line 1053:
             *   a1[it + ndt * il] = A[il, it]
             * This is row-major storage: a1[row * ncol + col]
             */
            double *a = lf->a1;

#pragma omp critical
            {
                for (int i = 0; i < ndl; i++) {
                    double sum = 0.0;
                    for (int j = 0; j < ndt; j++) {
                        sum += a[j + ndt * i] * xp[nstrtt - 1 + j];
                    }
                    yp[nstrtl - 1 + i] += sum;
                }
            }
        }
    }

    /* Apply inverse permutation to output */
    for (int i = 1; i <= n_dof; i++) {
        y[c->lod[i] - 1] = yp[i - 1];
    }

    free(xp);
    free(yp);
}

/*
 * Get DOF permutation array
 */
const int* HACApK_get_permutation(st_cHACApK_lcontrol ctl) {
    struct st_cHACApK_lcontrol_impl *c = (struct st_cHACApK_lcontrol_impl *)ctl;
    if (c == NULL) return NULL;
    return c->lod;
}

/*
 * Get total DOF count
 */
int HACApK_get_ndof(st_cHACApK_leafmtxp leafmtxp) {
    struct st_cHACApK_leafmtxp_impl *p = (struct st_cHACApK_leafmtxp_impl *)leafmtxp;
    if (p == NULL) return 0;
    return p->nd;
}

/*
 * Utility functions (wrappers around cHACApK_lib functions)
 */

double HACApK_norm2(const double *vec, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}

void HACApK_maxabs_loc(const double *vec, int n, double *maxval, int *loc) {
    *maxval = 0.0;
    *loc = 0;
    for (int i = 0; i < n; i++) {
        double absval = fabs(vec[i]);
        if (absval > *maxval) {
            *maxval = absval;
            *loc = i;
        }
    }
}

void HACApK_minabs_loc(const double *vec, int n, double *minval, int *loc) {
    if (n <= 0) {
        *minval = 0.0;
        *loc = 0;
        return;
    }
    *minval = fabs(vec[0]);
    *loc = 0;
    for (int i = 1; i < n; i++) {
        double absval = fabs(vec[i]);
        if (absval < *minval) {
            *minval = absval;
            *loc = i;
        }
    }
}

void HACApK_gemv_sub(double *vec, const double *A, const double *x, int nrow, int ncol, int lda) {
    for (int i = 0; i < nrow; i++) {
        double sum = 0.0;
        for (int j = 0; j < ncol; j++) {
            sum += A[i + j * lda] * x[j];
        }
        vec[i] -= sum;
    }
}

int HACApK_med3(int a, int b, int c) {
    if (a < b) {
        if (b < c) return b;
        else if (c < a) return a;
        else return c;
    } else {
        if (c < b) return b;
        else if (a < c) return a;
        else return c;
    }
}
