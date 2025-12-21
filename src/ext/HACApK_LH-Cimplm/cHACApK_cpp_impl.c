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

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * Note: This file is compiled as C (not C++)
 * The struct types here use the original HACApK typedefs.
 *
 * The wrapper functions take void* parameters that are cast to the
 * correct HACApK types internally.
 */

/* Helper to get number of threads */
static int get_num_threads(void) {
/* TEMPORARY: Force single thread for debugging */
#if 0 && defined(_OPENMP)
    return omp_get_max_threads();
#else
    return 1;
#endif
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
    ctl->param[60] = 2;                     /* ACA mode (2=ACA+) */
    ctl->param[61] = 1;                     /* ACA norm (MREM) */
    ctl->param[62] = (double)200;          /* Max rank initial */
    ctl->param[63] = (double)200;          /* Max rank */
    ctl->param[71] = eps;                   /* ACA tolerance */
    ctl->param[72] = 1.0e-3 * eps;          /* Relative tolerance */

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
        st_leafmtx[ip] = (st_cHACApK_leafmtx)calloc(1, sizeof(struct st_cHACApK_leafmtx));
        if (!st_leafmtx[ip]) {
            fprintf(stderr, "[HACApK] Error: Memory allocation for leaf %d failed\n", ip);
            return -1;
        }
    }

    /* Generate leaf matrix structure */
    for (il = 0; il < 5; il++) lnmtx[il] = 0;
    ndpth = 0;
    nlf = 0;
    cHACApK_generate_leafmtx(st_leafmtx, st_clt, st_clt, ctl->param, ctl->lpmd,
                              lnmtx, nofc, nffc, &nlf, &ndpth);

    if (print_level > 0) {
        printf("[HACApK] Generated %d leaf matrices\n", nlf);
    }

    /* Sort leaves for efficient traversal */
    cHACApK_sort_leafmtx(st_leafmtx, nlf);

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
        cHACApK_fill_leafmtx_hyp(st_leafmtx, i_bemv, ctl->param, znrmmat,
                                  ctl->lpmd, lnmtx, ctl->lod, ctl->lod, nd, nlf,
                                  lnps, lnpe, ctl->lthr);

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

    /* Cleanup temporary arrays */
    cHACApK_free_st_clt(st_clt);
    free(lodfc);
    for (il = 0; il <= ndim; il++) free(gmid_t[il]);
    free(gmid_t);

    return 0;
}

/*=========================================================================
 * Matrix-vector product: y = A * x using H-matrix
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
    double *x_perm, *y_perm;
    int il, it, ip, k;

    /* Allocate permuted vectors */
    x_perm = (double*)malloc(sizeof(double) * nd);
    y_perm = (double*)calloc(nd, sizeof(double));

    /* Pre-permute input vector */
    for (il = 1; il <= nd; il++) {
        x_perm[il - 1] = x[lod[il] - 1];
    }

    /* H-matrix matrix-vector product (single-threaded) */
    for (ip = 1; ip <= nlf; ip++) {
        st_cHACApK_leafmtx leaf = st_lf[ip];
        int ndl = leaf->ndl;
        int ndt = leaf->ndt;
        int nstrtl = leaf->nstrtl;
        int nstrtt = leaf->nstrtt;
        double *a1 = leaf->a1;
        double *a2 = leaf->a2;

        if (leaf->ltmtx == 1) {
            /* Low-rank block: y += U * (V^T * x)
             * Storage: a1 = V (ndt x kt), a2 = U (ndl x kt)
             */
            int kt = leaf->kt;
            double *V = a1;
            double *U = a2;
            double *tmp = (double*)calloc(kt, sizeof(double));

            /* Compute tmp = V^T * x_perm */
            for (k = 0; k < kt; k++) {
                double sum = 0.0;
                for (it = 0; it < ndt; it++) {
                    sum += V[it + ndt * k] * x_perm[nstrtt - 1 + it];
                }
                tmp[k] = sum;
            }

            /* Compute y_perm += U * tmp */
            for (il = 0; il < ndl; il++) {
                double sum = 0.0;
                for (k = 0; k < kt; k++) {
                    sum += U[il + ndl * k] * tmp[k];
                }
                y_perm[nstrtl - 1 + il] += sum;
            }

            free(tmp);
        } else {
            /* Dense block: y += A * x
             * Storage: a1[it + ndt * il] = A[il, it] (row-major)
             */
            for (il = 0; il < ndl; il++) {
                double sum = 0.0;
                for (it = 0; it < ndt; it++) {
                    sum += a1[it + ndt * il] * x_perm[nstrtt - 1 + it];
                }
                y_perm[nstrtl - 1 + il] += sum;
            }
        }
    }

    /* Inverse permutation: y(lod(i)) = y_perm(i) */
    for (il = 1; il <= nd; il++) {
        y[lod[il] - 1] = y_perm[il - 1];
    }

    free(x_perm);
    free(y_perm);
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
    return calloc(1, sizeof(struct st_cHACApK_leafmtxp));
}

void* HACApK_alloc_lcontrol(void) {
    return calloc(1, sizeof(struct st_cHACApK_lcontrol));
}

void HACApK_free_leafmtxp(void *ptr) {
    if (ptr) free(ptr);
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

    /* Iterate over all leaf blocks (single-threaded for debugging) */
    /* #pragma omp parallel for schedule(dynamic, 16) reduction(+:n_diag_updated) */
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

#if 0  /* Debug output - disabled */
    printf("[HACApK] Updated %d dense diagonal blocks out of %d total leaves\n", n_diag_updated, leafmtxp->nlf);
#endif
}
