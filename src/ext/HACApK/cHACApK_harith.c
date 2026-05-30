/* cHACApK_harith.c  -- H-matrix arithmetic (Phase 1/2 minimal: dense-leaf only).
 *
 * Block-recursive H-LU on the block-tree view (cHACApK_block_tree.h).
 * Leaves are assumed DENSE (ltmtx == 2). Low-rank leaves return an
 * error code for the caller to handle (skip / densify / fall-back).
 *
 * Layout assumption: HACApK dense leaf stores a1 as COLUMN-MAJOR
 * (ndl x ndt). This matches HACApK's call-pattern with LAPACK / BLAS.
 *
 * No pivoting (matches H2Lib's lrdecomp_hmatrix). For high-mu_r non-
 * symmetric MSC at large scale we will need to add diagonal-leaf
 * pivoting in Phase 3 (the LAPACKE_dgetrf row swaps stay local to the
 * leaf if we cap the pivot at the leaf boundary).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mkl.h>

#include "cHACApK_harith.h"

/* ---------- HACApK <-> internal storage convention --------------- *
 *
 * HACApK leaf storage convention (from cHACApK_base.c fill_leafmtx):
 *   ltmtx=1 (rk):  a1 is V of shape (ndt x kt) column-major
 *                  a2 is U of shape (ndl x kt) column-major
 *                  matrix M = U V^T = a2 a1^T  (ndl x ndt)
 *   ltmtx=2 (dense): a1 is M stored ROW-MAJOR with leading dim ndt
 *                    i.e., a1[col + row*ndt] = M[row, col]
 *                    Equivalently, a1 is M^T column-major with ldA = ndt.
 *
 * Our internal H-LU convention (matches build_deep_tree self-test):
 *   ltmtx=1 (rk):  a1 is U of shape (ndl x kt) column-major
 *                  a2 is V of shape (ndt x kt) column-major
 *                  matrix M = U V^T = a1 a2^T  (ndl x ndt)
 *   ltmtx=2 (dense): a1 is M column-major with leading dim ndl
 *                    i.e., a1[row + col*ndl] = M[row, col]
 *
 * Conversion is in-place:
 *   - rk:  swap a1 <-> a2 pointers  (O(1))
 *   - dense: transpose a1  (O(ndl * ndt) per leaf, requires temp buffer
 *            for non-square; in-place transpose for square via swap loop) */

void cHACApK_convert_leafmtxp_to_internal(struct st_cHACApK_leafmtxp_t *lp)
{
    if (!lp || !lp->st_lf) return;
    for (int ip = 1; ip <= lp->nlf; ip++) {
        st_cHACApK_leafmtx_t *lf = lp->st_lf[ip];
        if (!lf) continue;
        int ndl = lf->ndl, ndt = lf->ndt;
        if (lf->ltmtx == 1) {
            /* rk: swap a1 (V) <-> a2 (U) so my code finds U at a1. */
            double *tmp = lf->a1; lf->a1 = lf->a2; lf->a2 = tmp;
        } else if (lf->ltmtx == 2) {
            /* dense: transpose a1 from HACApK row-major (ldA=ndt) to
             * internal column-major (ldA=ndl). */
            double *new_a1 = (double*)malloc(sizeof(double) * (size_t)ndl * (size_t)ndt);
            if (!new_a1) continue;  /* allocation failure: skip (caller will see wrong results) */
            for (int row = 0; row < ndl; row++)
                for (int col = 0; col < ndt; col++)
                    new_a1[row + (size_t)col * (size_t)ndl] =
                        lf->a1[col + (size_t)row * (size_t)ndt];
            free(lf->a1);
            lf->a1 = new_a1;
        }
    }
}

void cHACApK_convert_leafmtxp_to_hacapk(struct st_cHACApK_leafmtxp_t *lp)
{
    if (!lp || !lp->st_lf) return;
    for (int ip = 1; ip <= lp->nlf; ip++) {
        st_cHACApK_leafmtx_t *lf = lp->st_lf[ip];
        if (!lf) continue;
        int ndl = lf->ndl, ndt = lf->ndt;
        if (lf->ltmtx == 1) {
            double *tmp = lf->a1; lf->a1 = lf->a2; lf->a2 = tmp;
        } else if (lf->ltmtx == 2) {
            /* transpose internal col-major (ldA=ndl) back to HACApK row-major (ldA=ndt) */
            double *new_a1 = (double*)malloc(sizeof(double) * (size_t)ndl * (size_t)ndt);
            if (!new_a1) continue;
            for (int row = 0; row < ndl; row++)
                for (int col = 0; col < ndt; col++)
                    new_a1[col + (size_t)row * (size_t)ndt] =
                        lf->a1[row + (size_t)col * (size_t)ndl];
            free(lf->a1);
            lf->a1 = new_a1;
        }
    }
}


/* ---------- module-level stats (single-threaded, single-decomp) ---- */
static cHACApK_hlu_stats_t g_stats;

const cHACApK_hlu_stats_t *cHACApK_hlu_last_stats(void) { return &g_stats; }

static void stats_reset(void) { memset(&g_stats, 0, sizeof(g_stats)); }

/* ---------- leaf utilities --------------------------------------- */

static inline int leaf_is_dense(const st_cHACApK_block_node_t *n)
{ return n && n->leaf_mtx && n->leaf_kind == 2; }

static inline int leaf_is_rk(const st_cHACApK_block_node_t *n)
{ return n && n->leaf_mtx && n->leaf_kind == 1; }

static inline double *leaf_dense_data(const st_cHACApK_block_node_t *n)
{ return n->leaf_mtx->a1; }

static inline int leaf_rows(const st_cHACApK_block_node_t *n)
{ return n->leaf_mtx->ndl; }

static inline int leaf_cols(const st_cHACApK_block_node_t *n)
{ return n->leaf_mtx->ndt; }

/* ---------- rk leaf utilities ---------------------------------- *
 * HACApK rk leaf storage (ltmtx=1):
 *   a1 = U,  shape ndl x kt, column-major
 *   a2 = V,  shape ndt x kt, column-major
 *   matrix represented = U * V^T  (ndl x ndt)
 * The rank kt is stored in leaf_mtx->kt. */

static inline int leaf_rk_rank(const st_cHACApK_block_node_t *n)
{ return n->leaf_mtx->kt; }

static inline double *leaf_rk_U(const st_cHACApK_block_node_t *n)
{ return n->leaf_mtx->a1; }

static inline double *leaf_rk_V(const st_cHACApK_block_node_t *n)
{ return n->leaf_mtx->a2; }


/* ---------- dense-leaf LAPACK primitives ------------------------- */

/* In-place NO-PIVOT dense LU (Doolittle): A overwritten with L (unit-lower,
 * below the diagonal) and U (on and above the diagonal). For Phase 1
 * minimal we restrict to diagonally-dominant matrices (validated by the
 * self-test); Phase 3 will add diagonal-leaf pivoting for general MSC.
 *
 * This explicit loop avoids the LAPACKE_dgetrf pivot tracking that would
 * have to be propagated through every off-diagonal trsm and the solve
 * forward sweep -- significant complexity for the recursive layout.
 * Matches H2Lib's lrdecomp_hmatrix which is also pivot-free. */
static int dense_lu_inplace(double *A, int n)
{
    for (int i = 0; i < n; i++) {
        double pivot = A[i + i * n];
        if (pivot == 0.0) return CHACAPK_HARITH_ERR_LAPACK;
        double inv_p = 1.0 / pivot;
        /* L column: A[j,i] /= A[i,i] for j > i */
        for (int j = i + 1; j < n; j++) A[j + i * n] *= inv_p;
        /* Update trailing block: A[j,k] -= A[j,i] * A[i,k] for j,k > i */
        for (int k = i + 1; k < n; k++) {
            double aik = A[i + k * n];
            for (int j = i + 1; j < n; j++) {
                A[j + k * n] -= A[j + i * n] * aik;
            }
        }
    }
    g_stats.n_dense_lu++;
    return CHACAPK_HARITH_OK;
}
/* dummy stubs (no longer used; preserved-comment for diff readability) */
static int *find_ipiv(void *leaf_key) { (void)leaf_key; return NULL; }
static void clear_ipiv_registry(void) { }

/* L_diag is an n_diag x n_diag dense LU factor (lower in-place + upper).
 * Solve L_diag X = B  (X overwrites B), B is m x n. */
static void dense_tri_solve_left_lower(
    const double *L_diag, int n_diag,
    double *B, int m, int n)
{
    /* dtrsm: B := alpha * L^{-1} * B, with L unit-lower */
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                m, n, 1.0, L_diag, n_diag, B, m);
    g_stats.n_dense_trsm++;
}

/* Solve X U_diag = B (X overwrites B), U_diag has upper triangle of the LU. */
static void dense_tri_solve_right_upper(
    const double *U_diag, int n_diag,
    double *B, int m, int n)
{
    cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
                m, n, 1.0, U_diag, n_diag, B, m);
    g_stats.n_dense_trsm++;
}

/* ---------- rk truncation via SVD recompression ------------------ *
 *
 * Given U (m x k_in) and V (n x k_in) defining an m x n matrix M = U V^T,
 * find U_new (m x k_new) and V_new (n x k_new) with k_new <= k_in such that
 * U_new V_new^T ~= M.  The truncation is determined by:
 *   tol_rel: keep singular values sigma_i > tol_rel * sigma_max
 *   k_max:   hard cap on the new rank
 *
 * Algorithm (standard, ~H2Lib's trunc_rkmatrix):
 *   1. QR of U:  U = Q_U R_U     (Q_U: m x k_in, R_U: k_in x k_in upper)
 *   2. QR of V:  V = Q_V R_V     (Q_V: n x k_in, R_V: k_in x k_in upper)
 *   3. S = R_U R_V^T             (k_in x k_in)
 *   4. SVD: S = Us Sigma Vs^T    (LAPACK returns V^T already)
 *   5. k_new = #{sigma_i > tol_rel * sigma_max}, capped at k_max
 *   6. U_new = Q_U Us[:,:k_new] diag(sigma[:k_new])
 *      V_new = Q_V (Vs^T)[:k_new,:]^T     (= Q_V Vs[:,:k_new])
 *
 * Returns CHACAPK_HARITH_OK on success and writes the new rk factors into
 * *U_new_out, *V_new_out (caller frees), with *k_new_out giving the rank.
 * On any LAPACK error, returns CHACAPK_HARITH_ERR_LAPACK and the outputs
 * are left NULL / 0.  Memory: O(m k_in + n k_in + k_in^2) temporary. */
static int rkleaf_recompress(
    const double *U, const double *V, int m, int n, int k_in,
    double tol_rel, int k_max,
    double **U_new_out, double **V_new_out, int *k_new_out)
{
    *U_new_out = NULL;
    *V_new_out = NULL;
    *k_new_out = 0;
    if (k_in <= 0 || m <= 0 || n <= 0) return CHACAPK_HARITH_ERR_NULL;

    /* Workspace copies (LAPACK overwrites inputs). */
    double *U_work = (double*)malloc(sizeof(double) * (size_t)m * (size_t)k_in);
    double *V_work = (double*)malloc(sizeof(double) * (size_t)n * (size_t)k_in);
    double *R_U    = (double*)calloc((size_t)k_in * (size_t)k_in, sizeof(double));
    double *R_V    = (double*)calloc((size_t)k_in * (size_t)k_in, sizeof(double));
    double *tau    = (double*)malloc(sizeof(double) * (size_t)k_in);
    if (!U_work || !V_work || !R_U || !R_V || !tau) {
        free(U_work); free(V_work); free(R_U); free(R_V); free(tau);
        return CHACAPK_HARITH_ERR_NULL;
    }
    memcpy(U_work, U, sizeof(double) * (size_t)m * (size_t)k_in);
    memcpy(V_work, V, sizeof(double) * (size_t)n * (size_t)k_in);

    /* QR of U_work -> Q stored implicitly, R in upper triangle. */
    int info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m, k_in, U_work, m, tau);
    if (info != 0) goto lapack_err;
    for (int j = 0; j < k_in; j++)
        for (int i = 0; i <= j; i++) R_U[i + j*k_in] = U_work[i + j*m];
    info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, m, k_in, k_in, U_work, m, tau);
    if (info != 0) goto lapack_err;
    /* Now U_work = Q_U (m x k_in). */

    /* QR of V_work. */
    info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, n, k_in, V_work, n, tau);
    if (info != 0) goto lapack_err;
    for (int j = 0; j < k_in; j++)
        for (int i = 0; i <= j; i++) R_V[i + j*k_in] = V_work[i + j*n];
    info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, n, k_in, k_in, V_work, n, tau);
    if (info != 0) goto lapack_err;
    /* Now V_work = Q_V (n x k_in). */

    free(tau); tau = NULL;

    /* S = R_U R_V^T  (k_in x k_in). */
    double *S = (double*)malloc(sizeof(double) * (size_t)k_in * (size_t)k_in);
    if (!S) { free(U_work); free(V_work); free(R_U); free(R_V);
              return CHACAPK_HARITH_ERR_NULL; }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                k_in, k_in, k_in, 1.0,
                R_U, k_in, R_V, k_in, 0.0, S, k_in);
    free(R_U); free(R_V);

    /* SVD: S = Us Sigma Vs^T. */
    double *sigma  = (double*)malloc(sizeof(double) * (size_t)k_in);
    double *Us     = (double*)malloc(sizeof(double) * (size_t)k_in * (size_t)k_in);
    double *Vt     = (double*)malloc(sizeof(double) * (size_t)k_in * (size_t)k_in);
    double *superb = (double*)malloc(sizeof(double) * (size_t)(k_in > 1 ? k_in - 1 : 1));
    if (!sigma || !Us || !Vt || !superb) {
        free(sigma); free(Us); free(Vt); free(superb);
        free(S); free(U_work); free(V_work);
        return CHACAPK_HARITH_ERR_NULL;
    }
    info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'S', k_in, k_in, S, k_in,
                          sigma, Us, k_in, Vt, k_in, superb);
    free(superb); free(S);
    if (info != 0) { free(sigma); free(Us); free(Vt);
                     free(U_work); free(V_work);
                     return CHACAPK_HARITH_ERR_LAPACK; }

    /* Determine new rank. */
    int k_new = 0;
    double sv_max = sigma[0];
    if (sv_max > 0.0) {
        double thresh = tol_rel * sv_max;
        for (int i = 0; i < k_in; i++) {
            if (sigma[i] > thresh) k_new++;
            else break;
        }
        if (k_new > k_max) k_new = k_max;
        if (k_new < 1)     k_new = 1;
    } else {
        k_new = 1;  /* rank-1 to keep something */
    }

    /* U_new = Q_U Us[:,:k_new] diag(sigma[:k_new])
     *       = Q_U (Us scaled by sigma columns)  (m x k_new) */
    double *U_new = (double*)malloc(sizeof(double) * (size_t)m * (size_t)k_new);
    double *V_new = (double*)malloc(sizeof(double) * (size_t)n * (size_t)k_new);
    double *Us_scaled = (double*)malloc(sizeof(double) * (size_t)k_in * (size_t)k_new);
    if (!U_new || !V_new || !Us_scaled) {
        free(U_new); free(V_new); free(Us_scaled);
        free(sigma); free(Us); free(Vt);
        free(U_work); free(V_work);
        return CHACAPK_HARITH_ERR_NULL;
    }
    for (int j = 0; j < k_new; j++) {
        double s = sigma[j];
        for (int i = 0; i < k_in; i++)
            Us_scaled[i + j*k_in] = Us[i + j*k_in] * s;
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                m, k_new, k_in, 1.0,
                U_work, m, Us_scaled, k_in, 0.0, U_new, m);
    free(Us_scaled);

    /* V_new = Q_V Vs[:,:k_new] = Q_V (Vt[:k_new,:])^T.
     * Pass Vt with ldb = k_in; CblasTrans + N=k_new accesses only the first
     * k_new rows of Vt as the "B" matrix (B is logically N x K = k_new x k_in
     * post-Trans). */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                n, k_new, k_in, 1.0,
                V_work, n, Vt, k_in, 0.0, V_new, n);

    free(sigma); free(Us); free(Vt);
    free(U_work); free(V_work);

    *U_new_out = U_new;
    *V_new_out = V_new;
    *k_new_out = k_new;
    return CHACAPK_HARITH_OK;

lapack_err:
    free(U_work); free(V_work); free(R_U); free(R_V); free(tau);
    return CHACAPK_HARITH_ERR_LAPACK;
}


/* ---------- Phase 3.6 helpers: dense materialization + distribution ---- *
 *
 * For mixed leaf+internal cases in h_addmul / htrsm we fall back to a
 * "materialize + flat BLAS + distribute back" path. Cost is O(m*n) for
 * materialization plus the natural cost of the BLAS operation. For real
 * HACApK trees where mixed cases occur at intermediate levels (not the
 * root), the cost is bounded by the largest mixed-case sub-block.
 *
 * Also handles the dense(A) * dense(B) -> rk(C) Phase 3.5 extension
 * case (form C_full = U_c V_c^T + alpha A B, SVD-truncate to rk).
 */

/* SVD-truncate an m x n dense matrix to rk factors U (m x k_new) V (n x k_new)
 * with k_new determined by relative tolerance against sigma_max, capped at
 * k_max. Returns CHACAPK_HARITH_OK on success and writes new factors into
 * *U_out, *V_out (caller frees), *k_out the rank. */
static int dense_to_rk_truncate(
    const double *dense, int m, int n,
    double tol_rel, int k_max,
    double **U_out, double **V_out, int *k_out)
{
    *U_out = NULL; *V_out = NULL; *k_out = 0;
    if (m <= 0 || n <= 0) return CHACAPK_HARITH_ERR_NULL;
    int kmin = (m < n) ? m : n;
    if (k_max <= 0 || k_max > kmin) k_max = kmin;

    double *work   = (double*)malloc(sizeof(double) * (size_t)m * (size_t)n);
    double *sigma  = (double*)malloc(sizeof(double) * (size_t)kmin);
    double *Us     = (double*)malloc(sizeof(double) * (size_t)m * (size_t)kmin);
    double *Vt     = (double*)malloc(sizeof(double) * (size_t)kmin * (size_t)n);
    double *superb = (double*)malloc(sizeof(double) * (size_t)(kmin > 1 ? kmin - 1 : 1));
    if (!work || !sigma || !Us || !Vt || !superb) {
        free(work); free(sigma); free(Us); free(Vt); free(superb);
        return CHACAPK_HARITH_ERR_NULL;
    }
    memcpy(work, dense, sizeof(double) * (size_t)m * (size_t)n);

    int info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, work, m,
                              sigma, Us, m, Vt, kmin, superb);
    free(work); free(superb);
    if (info != 0) { free(sigma); free(Us); free(Vt); return CHACAPK_HARITH_ERR_LAPACK; }

    int k_new = 0;
    double sv_max = sigma[0];
    if (sv_max > 0.0) {
        double thresh = tol_rel * sv_max;
        for (int i = 0; i < kmin; i++) {
            if (sigma[i] > thresh) k_new++;
            else break;
        }
        if (k_new > k_max) k_new = k_max;
        if (k_new < 1)     k_new = 1;
    } else {
        k_new = 1;
    }

    double *U_new = (double*)malloc(sizeof(double) * (size_t)m * (size_t)k_new);
    double *V_new = (double*)malloc(sizeof(double) * (size_t)n * (size_t)k_new);
    if (!U_new || !V_new) {
        free(sigma); free(Us); free(Vt); free(U_new); free(V_new);
        return CHACAPK_HARITH_ERR_NULL;
    }
    /* U_new = U_s[:,:k_new] * diag(sigma[:k_new])  (m x k_new col-major). */
    for (int j = 0; j < k_new; j++) {
        double s = sigma[j];
        for (int i = 0; i < m; i++)
            U_new[i + (size_t)j*(size_t)m] = Us[i + (size_t)j*(size_t)m] * s;
    }
    /* V_new = (Vt[:k_new, :])^T  (n x k_new col-major). */
    for (int j = 0; j < n; j++)
        for (int i = 0; i < k_new; i++)
            V_new[j + (size_t)i*(size_t)n] = Vt[i + (size_t)j*(size_t)kmin];

    free(sigma); free(Us); free(Vt);
    *U_out = U_new; *V_out = V_new; *k_out = k_new;
    return CHACAPK_HARITH_OK;
}

/* Materialize a block-tree node as a dense buffer (col-major, leading dim
 * = row_cluster->nsize). Works on dense leaf, rk leaf, or internal node.
 * Returns heap-allocated buffer (caller frees) or NULL on error. */
static double *materialize_node_as_dense(const st_cHACApK_block_node_t *node)
{
    if (!node) return NULL;
    int m = node->row_cluster->nsize;
    int n = node->col_cluster->nsize;
    double *out = (double*)calloc((size_t)m * (size_t)n, sizeof(double));
    if (!out) return NULL;

    if (leaf_is_dense(node)) {
        memcpy(out, leaf_dense_data(node), sizeof(double) * (size_t)m * (size_t)n);
        return out;
    }
    if (leaf_is_rk(node)) {
        int kt = leaf_rk_rank(node);
        if (kt > 0) {
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                        m, n, kt, 1.0,
                        leaf_rk_U(node), m,
                        leaf_rk_V(node), n,
                        0.0, out, m);
        }
        return out;
    }
    /* Internal: recurse on children, copy each into its slot in out. */
    int nr = node->nrsons, nc = node->ncsons;
    int row_base = node->row_cluster->nstrt;
    int col_base = node->col_cluster->nstrt;
    for (int j_s = 0; j_s < nc; j_s++) {
        for (int i_s = 0; i_s < nr; i_s++) {
            const st_cHACApK_block_node_t *child = node->sons[i_s + j_s * nr];
            int cm = child->row_cluster->nsize;
            int cn = child->col_cluster->nsize;
            int row_off = child->row_cluster->nstrt - row_base;
            int col_off = child->col_cluster->nstrt - col_base;
            double *cd = materialize_node_as_dense(child);
            if (!cd) { free(out); return NULL; }
            for (int jj = 0; jj < cn; jj++) {
                memcpy(&out[row_off + (size_t)(col_off + jj) * (size_t)m],
                       &cd[(size_t)jj * (size_t)cm],
                       sizeof(double) * (size_t)cm);
            }
            free(cd);
        }
    }
    return out;
}

/* Add a dense buffer D (column-major, leading dim D_ld) of size m x n
 * to a block-tree node C (which must span the same m x n region).
 * - Dense leaf C: in-place add into C's data.
 * - Rk leaf C: form C_full = U_c V_c^T, add D, SVD-truncate back to rk.
 * - Internal C: recurse with sub-views of D into each child. */
static int add_dense_to_node(
    const double *D, int D_ld,
    int m, int n,
    st_cHACApK_block_node_t *C)
{
    if (!D || !C) return CHACAPK_HARITH_ERR_NULL;

    if (leaf_is_dense(C)) {
        if (leaf_rows(C) != m || leaf_cols(C) != n) return CHACAPK_HARITH_ERR_TOPOLOGY;
        double *cdata = leaf_dense_data(C);
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                cdata[i + (size_t)j * (size_t)m] += D[i + (size_t)j * (size_t)D_ld];
        return CHACAPK_HARITH_OK;
    }
    if (leaf_is_rk(C)) {
        int kt_c = leaf_rk_rank(C);
        double *C_full = (double*)malloc(sizeof(double) * (size_t)m * (size_t)n);
        if (!C_full) return CHACAPK_HARITH_ERR_NULL;
        /* C_full = U_c V_c^T */
        if (kt_c > 0) {
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                        m, n, kt_c, 1.0,
                        leaf_rk_U(C), m,
                        leaf_rk_V(C), n,
                        0.0, C_full, m);
        } else {
            memset(C_full, 0, sizeof(double) * (size_t)m * (size_t)n);
        }
        /* C_full += D (with possibly different leading dim) */
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                C_full[i + (size_t)j * (size_t)m] += D[i + (size_t)j * (size_t)D_ld];
        /* SVD-truncate to rk */
        double *U_new = NULL, *V_new = NULL; int k_new = 0;
        int kmin = (m < n) ? m : n;
        int rc = dense_to_rk_truncate(C_full, m, n, 1e-14, kmin,
                                       &U_new, &V_new, &k_new);
        free(C_full);
        if (rc != CHACAPK_HARITH_OK) { free(U_new); free(V_new); return rc; }
        free(C->leaf_mtx->a1); free(C->leaf_mtx->a2);
        C->leaf_mtx->a1 = U_new; C->leaf_mtx->a2 = V_new;
        C->leaf_mtx->kt = k_new;
        return CHACAPK_HARITH_OK;
    }
    /* Internal: recurse on children with sub-views of D. */
    int nr = C->nrsons, nc = C->ncsons;
    int row_base = C->row_cluster->nstrt;
    int col_base = C->col_cluster->nstrt;
    for (int j_s = 0; j_s < nc; j_s++) {
        for (int i_s = 0; i_s < nr; i_s++) {
            st_cHACApK_block_node_t *child = C->sons[i_s + j_s * nr];
            int cm = child->row_cluster->nsize;
            int cn = child->col_cluster->nsize;
            int row_off = child->row_cluster->nstrt - row_base;
            int col_off = child->col_cluster->nstrt - col_base;
            const double *D_sub = D + row_off + (size_t)col_off * (size_t)D_ld;
            int rc = add_dense_to_node(D_sub, D_ld, cm, cn, child);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
    }
    return CHACAPK_HARITH_OK;
}

/* Replace a node's content with a dense buffer D (overwrite, not add).
 * Used by htrsm mixed fallback to put the trsm result back into X's
 * structure.
 * - Dense leaf: memcpy.
 * - Rk leaf: SVD-truncate D.
 * - Internal: recurse with sub-views. */
static int set_node_from_dense(
    const double *D, int D_ld,
    int m, int n,
    st_cHACApK_block_node_t *C)
{
    if (!D || !C) return CHACAPK_HARITH_ERR_NULL;

    if (leaf_is_dense(C)) {
        if (leaf_rows(C) != m || leaf_cols(C) != n) return CHACAPK_HARITH_ERR_TOPOLOGY;
        double *cdata = leaf_dense_data(C);
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                cdata[i + (size_t)j * (size_t)m] = D[i + (size_t)j * (size_t)D_ld];
        return CHACAPK_HARITH_OK;
    }
    if (leaf_is_rk(C)) {
        /* Copy D into a contiguous (m x n) buffer for SVD if D_ld != m. */
        const double *D_use = D;
        double *D_packed = NULL;
        if (D_ld != m) {
            D_packed = (double*)malloc(sizeof(double) * (size_t)m * (size_t)n);
            if (!D_packed) return CHACAPK_HARITH_ERR_NULL;
            for (int j = 0; j < n; j++)
                memcpy(&D_packed[(size_t)j * (size_t)m],
                       &D[(size_t)j * (size_t)D_ld],
                       sizeof(double) * (size_t)m);
            D_use = D_packed;
        }
        double *U_new = NULL, *V_new = NULL; int k_new = 0;
        int kmin = (m < n) ? m : n;
        int rc = dense_to_rk_truncate(D_use, m, n, 1e-14, kmin,
                                       &U_new, &V_new, &k_new);
        free(D_packed);
        if (rc != CHACAPK_HARITH_OK) { free(U_new); free(V_new); return rc; }
        free(C->leaf_mtx->a1); free(C->leaf_mtx->a2);
        C->leaf_mtx->a1 = U_new; C->leaf_mtx->a2 = V_new;
        C->leaf_mtx->kt = k_new;
        return CHACAPK_HARITH_OK;
    }
    /* Internal: recurse on children with sub-views of D. */
    int nr = C->nrsons, nc = C->ncsons;
    int row_base = C->row_cluster->nstrt;
    int col_base = C->col_cluster->nstrt;
    for (int j_s = 0; j_s < nc; j_s++) {
        for (int i_s = 0; i_s < nr; i_s++) {
            st_cHACApK_block_node_t *child = C->sons[i_s + j_s * nr];
            int cm = child->row_cluster->nsize;
            int cn = child->col_cluster->nsize;
            int row_off = child->row_cluster->nstrt - row_base;
            int col_off = child->col_cluster->nstrt - col_base;
            const double *D_sub = D + row_off + (size_t)col_off * (size_t)D_ld;
            int rc = set_node_from_dense(D_sub, D_ld, cm, cn, child);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
    }
    return CHACAPK_HARITH_OK;
}


/* ---------- block-recursive H-matrix primitives ------------------ *
 *
 * Three primitives, each dispatching on operand leaf_kind:
 *
 *   h_addmul(alpha, A, B, C):       C += alpha * A * B
 *   htrsm_lln(L, X):             solve L * X' = X, L unit-lower (Doolittle)
 *   htrsm_run(U, X):             solve X' * U = X, U non-unit upper
 *
 * Dispatch table:
 *   - all dense leaves         -> direct BLAS call (base case)
 *   - all internal nodes       -> block-recursive descent
 *   - any rk leaf              -> CHACAPK_HARITH_ERR_LOWRANK_LEAF (Phase 3)
 *   - mixed leaf+internal      -> CHACAPK_HARITH_ERR_NEED_RECURSIVE (Phase 3)
 *
 * For build_deep_tree's uniform-split self-test the mixed case cannot occur:
 * a block-node's children are either all leaves (depth 0 subtree) or all
 * internal (depth >= 1 subtree). In real HACApK trees the mixed case can
 * appear at the admissibility boundary; Phase 3 handles it by densifying
 * the internal operand at the leaf size (still O(leaf^3), bounded).
 *
 * Block-recursive descent for h_addmul (sum-over-j):
 *
 *   C(i,k) += sum_{j} alpha * A(i,j) * B(j,k)
 *
 * For htrsm_lln on block-(s x s) L and block-(s x ncX) X:
 *
 *   for j in 0..ncX-1:                   // each column-block of X
 *     for i in 0..s-1:                   // each row-block of X
 *       for k in 0..i-1:                 // discharge above-row dependencies
 *         X(i,j) -= L(i,k) * X(k,j)      // h_addmul with alpha=-1
 *       solve L(i,i) * X(i,j) = X(i,j)   // recursive htrsm_lln
 *
 * htrsm_run is the mirror: walks j (column of X) left-to-right, accumulates
 * X(i,j) -= X(i,k) * U(k,j) for k < j, then trsm right-upper.
 */

/* Forward declarations: htrsm_* both call h_addmul; h_addmul is self-recursive. */
static int h_addmul(double alpha,
                 const st_cHACApK_block_node_t *A,
                 const st_cHACApK_block_node_t *B,
                 st_cHACApK_block_node_t *C);

static int htrsm_lln(const st_cHACApK_block_node_t *L,
                     st_cHACApK_block_node_t *X);

static int htrsm_run(const st_cHACApK_block_node_t *U,
                     st_cHACApK_block_node_t *X);


/* C += alpha * A * B  (block-recursive, leaf-leaf base case dispatches on rk/dense). */
static int h_addmul(double alpha,
                 const st_cHACApK_block_node_t *A,
                 const st_cHACApK_block_node_t *B,
                 st_cHACApK_block_node_t *C)
{
    if (!A || !B || !C) return CHACAPK_HARITH_ERR_NULL;

    /* Base case: all dense leaves -> direct dgemm. */
    if (leaf_is_dense(A) && leaf_is_dense(B) && leaf_is_dense(C)) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    leaf_rows(C), leaf_cols(C), leaf_cols(A),
                    alpha,
                    leaf_dense_data(A), leaf_rows(A),
                    leaf_dense_data(B), leaf_rows(B),
                    1.0, leaf_dense_data(C), leaf_rows(C));
        g_stats.n_dense_gemm++;
        return CHACAPK_HARITH_OK;
    }

    /* rk(A) * dense(B) -> dense(C):
     *   C += alpha * (U_a V_a^T) B  =  (alpha U_a) (V_a^T B)
     * Cost: O(kt_a * inner * nC) + O(mC * kt_a * nC). */
    if (leaf_is_rk(A) && leaf_is_dense(B) && leaf_is_dense(C)) {
        int kt = leaf_rk_rank(A);
        int mC = leaf_rows(C);
        int nC = leaf_cols(C);
        int inner = leaf_cols(A);  /* == leaf_rows(B) */
        double *W = (double*)malloc(sizeof(double) * (size_t)kt * (size_t)nC);
        if (!W) return CHACAPK_HARITH_ERR_NULL;
        /* W = V_a^T B  (kt x nC) */
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    kt, nC, inner, 1.0,
                    leaf_rk_V(A), inner,
                    leaf_dense_data(B), inner, 0.0,
                    W, kt);
        /* C += alpha * U_a * W */
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    mC, nC, kt, alpha,
                    leaf_rk_U(A), mC,
                    W, kt, 1.0,
                    leaf_dense_data(C), mC);
        free(W);
        g_stats.n_dense_gemm++;
        return CHACAPK_HARITH_OK;
    }

    /* dense(A) * rk(B) -> dense(C):
     *   C += alpha * A (U_b V_b^T)  =  (alpha A U_b) V_b^T */
    if (leaf_is_dense(A) && leaf_is_rk(B) && leaf_is_dense(C)) {
        int kt = leaf_rk_rank(B);
        int mC = leaf_rows(C);
        int nC = leaf_cols(C);
        int inner = leaf_cols(A);  /* == leaf_rows(B) */
        double *W = (double*)malloc(sizeof(double) * (size_t)mC * (size_t)kt);
        if (!W) return CHACAPK_HARITH_ERR_NULL;
        /* W = alpha * A * U_b  (mC x kt) */
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    mC, kt, inner, alpha,
                    leaf_dense_data(A), mC,
                    leaf_rk_U(B), inner, 0.0,
                    W, mC);
        /* C += W * V_b^T */
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    mC, nC, kt, 1.0,
                    W, mC,
                    leaf_rk_V(B), nC, 1.0,
                    leaf_dense_data(C), mC);
        free(W);
        g_stats.n_dense_gemm++;
        return CHACAPK_HARITH_OK;
    }

    /* rk(A) * rk(B) -> dense(C):
     *   C += alpha * U_a (V_a^T U_b) V_b^T */
    if (leaf_is_rk(A) && leaf_is_rk(B) && leaf_is_dense(C)) {
        int kA = leaf_rk_rank(A);
        int kB = leaf_rk_rank(B);
        int mC = leaf_rows(C);
        int nC = leaf_cols(C);
        int inner = leaf_cols(A);  /* == leaf_rows(B) */
        double *M = (double*)malloc(sizeof(double) * (size_t)kA * (size_t)kB);
        double *X = (double*)malloc(sizeof(double) * (size_t)mC * (size_t)kB);
        if (!M || !X) { free(M); free(X); return CHACAPK_HARITH_ERR_NULL; }
        /* M = V_a^T * U_b  (kA x kB) */
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    kA, kB, inner, 1.0,
                    leaf_rk_V(A), inner,
                    leaf_rk_U(B), inner, 0.0,
                    M, kA);
        /* X = alpha * U_a * M  (mC x kB) */
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    mC, kB, kA, alpha,
                    leaf_rk_U(A), mC,
                    M, kA, 0.0,
                    X, mC);
        /* C += X * V_b^T */
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    mC, nC, kB, 1.0,
                    X, mC,
                    leaf_rk_V(B), nC, 1.0,
                    leaf_dense_data(C), mC);
        free(M); free(X);
        g_stats.n_dense_gemm++;
        return CHACAPK_HARITH_OK;
    }

    /* ---- rk(C) cases: increment is rank-k, stack with C's factors then
     *      recompress via SVD truncation. ----
     *
     * Helper for the stack+recompress pattern: given an increment described
     * by (U_inc, V_inc) of size (m x k_inc, n x k_inc), update C's rk leaf
     * to represent  U_c_new V_c_new^T = U_c V_c^T + U_inc V_inc^T.
     *
     * Each rk(C) case below computes U_inc and V_inc differently, then calls
     * a common helper to perform the stack-and-recompress. */

    /* ---- rk(A) * rk(B) -> rk(C):
     *   increment = alpha U_a (V_a^T U_b) V_b^T
     *   U_inc = alpha * U_a * (V_a^T U_b)   (m x kt_B)
     *   V_inc = V_b                          (n x kt_B) */
    if (leaf_is_rk(A) && leaf_is_rk(B) && leaf_is_rk(C)) {
        int kA = leaf_rk_rank(A), kB = leaf_rk_rank(B), kC = leaf_rk_rank(C);
        int m = leaf_rows(C), n = leaf_cols(C);
        int inner = leaf_cols(A);
        int kw = kC + kB;
        /* M = V_a^T U_b  (kA x kB) */
        double *M = (double*)malloc(sizeof(double) * (size_t)kA * (size_t)kB);
        if (!M) return CHACAPK_HARITH_ERR_NULL;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    kA, kB, inner, 1.0,
                    leaf_rk_V(A), inner,
                    leaf_rk_U(B), inner, 0.0, M, kA);
        /* U_inc = alpha * U_a * M  (m x kB)  -- built directly inside U_widened */
        double *U_widened = (double*)malloc(sizeof(double) * (size_t)m * (size_t)kw);
        double *V_widened = (double*)malloc(sizeof(double) * (size_t)n * (size_t)kw);
        if (!U_widened || !V_widened) { free(M); free(U_widened); free(V_widened);
                                         return CHACAPK_HARITH_ERR_NULL; }
        memcpy(U_widened, leaf_rk_U(C), sizeof(double) * (size_t)m * (size_t)kC);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    m, kB, kA, alpha,
                    leaf_rk_U(A), m,
                    M, kA, 0.0,
                    U_widened + (size_t)m * (size_t)kC, m);
        free(M);
        memcpy(V_widened, leaf_rk_V(C), sizeof(double) * (size_t)n * (size_t)kC);
        memcpy(V_widened + (size_t)n * (size_t)kC, leaf_rk_V(B),
               sizeof(double) * (size_t)n * (size_t)kB);
        /* Recompress (tol = 1e-14, no rank cap for self-test fidelity). */
        double *U_new = NULL, *V_new = NULL; int k_new = 0;
        int rc = rkleaf_recompress(U_widened, V_widened, m, n, kw,
                                    1e-14, kw, &U_new, &V_new, &k_new);
        free(U_widened); free(V_widened);
        if (rc != CHACAPK_HARITH_OK) { free(U_new); free(V_new); return rc; }
        free(C->leaf_mtx->a1); free(C->leaf_mtx->a2);
        C->leaf_mtx->a1 = U_new; C->leaf_mtx->a2 = V_new;
        C->leaf_mtx->kt = k_new;
        g_stats.n_dense_gemm++;
        return CHACAPK_HARITH_OK;
    }

    /* ---- rk(A) * dense(B) -> rk(C):
     *   increment = alpha (U_a V_a^T) B = (alpha U_a) (B^T V_a)^T
     *   U_inc = alpha * U_a       (m x kA)
     *   V_inc = B^T V_a            (n x kA) */
    if (leaf_is_rk(A) && leaf_is_dense(B) && leaf_is_rk(C)) {
        int kA = leaf_rk_rank(A), kC = leaf_rk_rank(C);
        int m = leaf_rows(C), n = leaf_cols(C);
        int inner = leaf_cols(A);
        int kw = kC + kA;
        double *U_widened = (double*)malloc(sizeof(double) * (size_t)m * (size_t)kw);
        double *V_widened = (double*)malloc(sizeof(double) * (size_t)n * (size_t)kw);
        if (!U_widened || !V_widened) { free(U_widened); free(V_widened);
                                         return CHACAPK_HARITH_ERR_NULL; }
        memcpy(U_widened, leaf_rk_U(C), sizeof(double) * (size_t)m * (size_t)kC);
        /* U_inc = alpha * U_a placed into the right half of U_widened */
        double *U_inc_block = U_widened + (size_t)m * (size_t)kC;
        cblas_dcopy(m * kA, leaf_rk_U(A), 1, U_inc_block, 1);
        if (alpha != 1.0) cblas_dscal(m * kA, alpha, U_inc_block, 1);
        /* V_inc = B^T V_a  (n x kA) */
        memcpy(V_widened, leaf_rk_V(C), sizeof(double) * (size_t)n * (size_t)kC);
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    n, kA, inner, 1.0,
                    leaf_dense_data(B), inner,
                    leaf_rk_V(A), inner, 0.0,
                    V_widened + (size_t)n * (size_t)kC, n);
        /* Recompress. */
        double *U_new = NULL, *V_new = NULL; int k_new = 0;
        int rc = rkleaf_recompress(U_widened, V_widened, m, n, kw,
                                    1e-14, kw, &U_new, &V_new, &k_new);
        free(U_widened); free(V_widened);
        if (rc != CHACAPK_HARITH_OK) { free(U_new); free(V_new); return rc; }
        free(C->leaf_mtx->a1); free(C->leaf_mtx->a2);
        C->leaf_mtx->a1 = U_new; C->leaf_mtx->a2 = V_new;
        C->leaf_mtx->kt = k_new;
        g_stats.n_dense_gemm++;
        return CHACAPK_HARITH_OK;
    }

    /* ---- dense(A) * rk(B) -> rk(C):
     *   increment = alpha A (U_b V_b^T) = (alpha A U_b) V_b^T
     *   U_inc = alpha * A * U_b   (m x kB)
     *   V_inc = V_b               (n x kB) */
    if (leaf_is_dense(A) && leaf_is_rk(B) && leaf_is_rk(C)) {
        int kB = leaf_rk_rank(B), kC = leaf_rk_rank(C);
        int m = leaf_rows(C), n = leaf_cols(C);
        int inner = leaf_cols(A);
        int kw = kC + kB;
        double *U_widened = (double*)malloc(sizeof(double) * (size_t)m * (size_t)kw);
        double *V_widened = (double*)malloc(sizeof(double) * (size_t)n * (size_t)kw);
        if (!U_widened || !V_widened) { free(U_widened); free(V_widened);
                                         return CHACAPK_HARITH_ERR_NULL; }
        memcpy(U_widened, leaf_rk_U(C), sizeof(double) * (size_t)m * (size_t)kC);
        /* U_inc = alpha * A * U_b  (m x kB)  -- placed into right half of U_widened */
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    m, kB, inner, alpha,
                    leaf_dense_data(A), m,
                    leaf_rk_U(B), inner, 0.0,
                    U_widened + (size_t)m * (size_t)kC, m);
        memcpy(V_widened, leaf_rk_V(C), sizeof(double) * (size_t)n * (size_t)kC);
        memcpy(V_widened + (size_t)n * (size_t)kC, leaf_rk_V(B),
               sizeof(double) * (size_t)n * (size_t)kB);
        /* Recompress. */
        double *U_new = NULL, *V_new = NULL; int k_new = 0;
        int rc = rkleaf_recompress(U_widened, V_widened, m, n, kw,
                                    1e-14, kw, &U_new, &V_new, &k_new);
        free(U_widened); free(V_widened);
        if (rc != CHACAPK_HARITH_OK) { free(U_new); free(V_new); return rc; }
        free(C->leaf_mtx->a1); free(C->leaf_mtx->a2);
        C->leaf_mtx->a1 = U_new; C->leaf_mtx->a2 = V_new;
        C->leaf_mtx->kt = k_new;
        g_stats.n_dense_gemm++;
        return CHACAPK_HARITH_OK;
    }

    /* dense(A) * dense(B) -> rk(C) + all MIXED leaf+internal cases:
     * Phase 3.6 materialize-and-redo fallback.
     *
     *   1. Materialize A and B as flat dense buffers (handles all 3 cases:
     *      dense leaf -> memcpy, rk leaf -> U V^T, internal -> recursive).
     *   2. Compute D = alpha * A_dense * B_dense via cblas_dgemm.
     *   3. Add D to C via add_dense_to_node (handles dense leaf via in-place
     *      add, rk leaf via SVD recompression, internal via recursive
     *      sub-view distribution).
     *
     * Cost: O(m_A * n_A + n_A * n_B + m_A * n_A * n_B) per call.
     * For HACApK trees mixed cases occur at intermediate levels, bounded. */
    {
        int A_is_leaf = leaf_is_dense(A) || leaf_is_rk(A);
        int B_is_leaf = leaf_is_dense(B) || leaf_is_rk(B);
        int C_is_leaf = leaf_is_dense(C) || leaf_is_rk(C);
        int mixed = (!A_is_leaf || !B_is_leaf || !C_is_leaf) &&
                    (A_is_leaf  ||  B_is_leaf ||  C_is_leaf);
        int dense_dense_rk = leaf_is_dense(A) && leaf_is_dense(B) && leaf_is_rk(C);
        if (mixed || dense_dense_rk) {
            int m = C->row_cluster->nsize;
            int n = C->col_cluster->nsize;
            int inner = A->col_cluster->nsize;  /* == B->row_cluster->nsize */

            double *A_dense = materialize_node_as_dense(A);
            double *B_dense = materialize_node_as_dense(B);
            if (!A_dense || !B_dense) {
                free(A_dense); free(B_dense);
                return CHACAPK_HARITH_ERR_NULL;
            }
            double *D = (double*)malloc(sizeof(double) * (size_t)m * (size_t)n);
            if (!D) { free(A_dense); free(B_dense); return CHACAPK_HARITH_ERR_NULL; }
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                        m, n, inner, alpha,
                        A_dense, m,
                        B_dense, inner,
                        0.0, D, m);
            free(A_dense); free(B_dense);

            int rc = add_dense_to_node(D, m, m, n, C);
            free(D);
            g_stats.n_dense_gemm++;
            return rc;
        }
    }
    /* Should not reach here -- all leaf-leaf and mixed cases handled. */

    /* All internal: recurse on C's (i,k), summing over j.
     * Topology constraint: A and B share a common "inner" dimension. */
    int nrC = C->nrsons, ncC = C->ncsons;
    int nmA = A->ncsons;
    if (nrC != A->nrsons || ncC != B->ncsons || nmA != B->nrsons) {
        return CHACAPK_HARITH_ERR_TOPOLOGY;
    }
    for (int k = 0; k < ncC; k++) {
        for (int i = 0; i < nrC; i++) {
            for (int j = 0; j < nmA; j++) {
                int rc = h_addmul(alpha,
                               A->sons[i + j * nrC],
                               B->sons[j + k * nmA],
                               C->sons[i + k * nrC]);
                if (rc != CHACAPK_HARITH_OK) return rc;
            }
        }
    }
    return CHACAPK_HARITH_OK;
}


/* Solve L * X' = X  (X overwrites itself).
 * L is unit-lower (Doolittle: L_diag implicit unit, L_below stored).
 * X must have the same row partition as L. */
static int htrsm_lln(const st_cHACApK_block_node_t *L,
                     st_cHACApK_block_node_t *X)
{
    if (!L || !X) return CHACAPK_HARITH_ERR_NULL;

    /* L is the L factor of LU (unit-lower); never low-rank by construction. */
    if (leaf_is_rk(L)) {
        g_stats.n_lowrank_skip++;
        return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
    }

    if (leaf_is_dense(L) && leaf_is_dense(X)) {
        dense_tri_solve_left_lower(
            leaf_dense_data(L), leaf_rows(L),
            leaf_dense_data(X), leaf_rows(X), leaf_cols(X));
        return CHACAPK_HARITH_OK;
    }

    /* L * X = B with X rk (X = U_x V_x^T).  Solve L * Y = X for Y rk.
     *   Y = L^-1 X = L^-1 (U_x V_x^T) = (L^-1 U_x) V_x^T.
     * So U_y = L^-1 U_x (cheap dtrsm on the kt-column block U_x), V_y = V_x. */
    if (leaf_is_dense(L) && leaf_is_rk(X)) {
        int kt = leaf_rk_rank(X);
        int m  = leaf_rows(L);  /* == leaf_rows(X) */
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                    m, kt, 1.0,
                    leaf_dense_data(L), m,
                    leaf_rk_U(X), m);
        g_stats.n_dense_trsm++;
        return CHACAPK_HARITH_OK;
    }

    /* Phase 3.6 mixed fallback: materialize L and X, dtrsm, write back. */
    {
        int L_is_leaf = leaf_is_dense(L) || leaf_is_rk(L);  /* rk-L unreachable above */
        int X_is_leaf = leaf_is_dense(X) || leaf_is_rk(X);
        int mixed = (!L_is_leaf || !X_is_leaf) && (L_is_leaf || X_is_leaf);
        if (mixed) {
            int m = L->row_cluster->nsize;  /* == L->col_cluster->nsize (square) */
            int n = X->col_cluster->nsize;
            double *L_dense = materialize_node_as_dense(L);
            double *X_dense = materialize_node_as_dense(X);
            if (!L_dense || !X_dense) {
                free(L_dense); free(X_dense);
                return CHACAPK_HARITH_ERR_NULL;
            }
            /* L * Y = X (overwrite X_dense with Y). L is unit-lower (Doolittle). */
            cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                        m, n, 1.0, L_dense, m, X_dense, m);
            int rc = set_node_from_dense(X_dense, m, m, n, X);
            free(L_dense); free(X_dense);
            g_stats.n_dense_trsm++;
            return rc;
        }
    }
    /* Should not reach here -- all-internal handled below. */

    /* All internal. L is square (nrL == ncL = s); X has nrX = s, ncX arbitrary. */
    int s = L->nrsons;
    if (L->ncsons != s) return CHACAPK_HARITH_ERR_TOPOLOGY;
    if (X->nrsons != s) return CHACAPK_HARITH_ERR_TOPOLOGY;
    int ncX = X->ncsons;

    for (int j = 0; j < ncX; j++) {
        for (int i = 0; i < s; i++) {
            /* X(i,j) -= L(i,k) * X(k,j) for k = 0..i-1 */
            for (int k = 0; k < i; k++) {
                int rc = h_addmul(-1.0,
                               L->sons[i + k * s],
                               X->sons[k + j * s],
                               X->sons[i + j * s]);
                if (rc != CHACAPK_HARITH_OK) return rc;
            }
            /* solve L(i,i) * X(i,j) = X(i,j) */
            int rc = htrsm_lln(L->sons[i + i * s], X->sons[i + j * s]);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
    }
    return CHACAPK_HARITH_OK;
}


/* Solve X' * U = X  (X overwrites itself).
 * U is non-unit upper (in-place LU's upper part).
 * X must have the same column partition as U. */
static int htrsm_run(const st_cHACApK_block_node_t *U,
                     st_cHACApK_block_node_t *X)
{
    if (!U || !X) return CHACAPK_HARITH_ERR_NULL;

    /* U is the U factor of LU (non-unit upper); never low-rank by construction. */
    if (leaf_is_rk(U)) {
        g_stats.n_lowrank_skip++;
        return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
    }

    if (leaf_is_dense(U) && leaf_is_dense(X)) {
        dense_tri_solve_right_upper(
            leaf_dense_data(U), leaf_rows(U),
            leaf_dense_data(X), leaf_rows(X), leaf_cols(X));
        return CHACAPK_HARITH_OK;
    }

    /* X * U = B with X rk (X = U_x V_x^T).  Solve Y * U = X for Y rk.
     *   Y = X U^-1 = U_x (V_x^T U^-1) = U_x (U^-T V_x)^T.
     * So U_y = U_x, V_y = U^-T V_x. We solve U^T V_y = V_x in place. */
    if (leaf_is_dense(U) && leaf_is_rk(X)) {
        int kt = leaf_rk_rank(X);
        int n  = leaf_cols(U);  /* == leaf_cols(X) (rk's V_x has n rows) */
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
                    n, kt, 1.0,
                    leaf_dense_data(U), n,
                    leaf_rk_V(X), n);
        g_stats.n_dense_trsm++;
        return CHACAPK_HARITH_OK;
    }

    /* Phase 3.6 mixed fallback: materialize U and X, dtrsm right-upper, write back. */
    {
        int U_is_leaf = leaf_is_dense(U) || leaf_is_rk(U);  /* rk-U unreachable above */
        int X_is_leaf = leaf_is_dense(X) || leaf_is_rk(X);
        int mixed = (!U_is_leaf || !X_is_leaf) && (U_is_leaf || X_is_leaf);
        if (mixed) {
            int n_u = U->row_cluster->nsize;  /* U is square: ndl == ndt */
            int m_x = X->row_cluster->nsize;
            double *U_dense = materialize_node_as_dense(U);
            double *X_dense = materialize_node_as_dense(X);
            if (!U_dense || !X_dense) {
                free(U_dense); free(X_dense);
                return CHACAPK_HARITH_ERR_NULL;
            }
            /* X * U = X (overwrite). U non-unit upper. */
            cblas_dtrsm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit,
                        m_x, n_u, 1.0, U_dense, n_u, X_dense, m_x);
            int rc = set_node_from_dense(X_dense, m_x, m_x, n_u, X);
            free(U_dense); free(X_dense);
            g_stats.n_dense_trsm++;
            return rc;
        }
    }
    /* Should not reach here. */

    int s = U->nrsons;
    if (U->ncsons != s) return CHACAPK_HARITH_ERR_TOPOLOGY;
    if (X->ncsons != s) return CHACAPK_HARITH_ERR_TOPOLOGY;
    int nrX = X->nrsons;

    for (int j = 0; j < s; j++) {
        for (int i = 0; i < nrX; i++) {
            /* X(i,j) -= X(i,k) * U(k,j) for k = 0..j-1 */
            for (int k = 0; k < j; k++) {
                int rc = h_addmul(-1.0,
                               X->sons[i + k * nrX],
                               U->sons[k + j * s],
                               X->sons[i + j * nrX]);
                if (rc != CHACAPK_HARITH_OK) return rc;
            }
            int rc = htrsm_run(U->sons[j + j * s], X->sons[i + j * nrX]);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
    }
    return CHACAPK_HARITH_OK;
}


/* Recursive matvec: y -= node * x.
 * x and y are GLOBAL vectors (length N); node's row/col cluster nstrt offsets
 * select the active slices. Aliasing x == y is safe iff the node's row and
 * column clusters do not overlap (true for strictly off-diagonal blocks
 * during the forward/backward sweeps). */
static int hmatvec_subtract(
    const st_cHACApK_block_node_t *node,
    const double *x, double *y)
{
    if (!node) return CHACAPK_HARITH_ERR_NULL;

    if (leaf_is_dense(node)) {
        int m = leaf_rows(node), n = leaf_cols(node);
        int r0 = node->row_cluster->nstrt - 1;  /* 1-based -> 0-based */
        int c0 = node->col_cluster->nstrt - 1;
        cblas_dgemv(CblasColMajor, CblasNoTrans, m, n,
                    -1.0, leaf_dense_data(node), m,
                    &x[c0], 1, 1.0, &y[r0], 1);
        return CHACAPK_HARITH_OK;
    }

    /* rk leaf: y -= U_a (V_a^T x).  Two dgemvs through the kt-rank waist. */
    if (leaf_is_rk(node)) {
        int m = leaf_rows(node), n = leaf_cols(node);
        int kt = leaf_rk_rank(node);
        int r0 = node->row_cluster->nstrt - 1;
        int c0 = node->col_cluster->nstrt - 1;
        /* Stack allocation: kt is small (typically < 50 for ACA tol 1e-4). */
        double w_stack[256];
        double *w = (kt <= (int)(sizeof(w_stack)/sizeof(w_stack[0])))
                    ? w_stack
                    : (double*)malloc(sizeof(double) * (size_t)kt);
        if (!w) return CHACAPK_HARITH_ERR_NULL;
        /* w = V_a^T x  (V_a is n x kt, so transpose) */
        cblas_dgemv(CblasColMajor, CblasTrans, n, kt,
                    1.0, leaf_rk_V(node), n,
                    &x[c0], 1, 0.0, w, 1);
        /* y -= U_a * w */
        cblas_dgemv(CblasColMajor, CblasNoTrans, m, kt,
                    -1.0, leaf_rk_U(node), m,
                    w, 1, 1.0, &y[r0], 1);
        if (w != w_stack) free(w);
        return CHACAPK_HARITH_OK;
    }

    int nr = node->nrsons, nc = node->ncsons;
    for (int j = 0; j < nc; j++) {
        for (int i = 0; i < nr; i++) {
            int rc = hmatvec_subtract(node->sons[i + j * nr], x, y);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
    }
    return CHACAPK_HARITH_OK;
}


/* ---------- recursive H-LU --------------------------------------- */

static int hlu_rec(st_cHACApK_block_node_t *node)
{
    if (!node) return CHACAPK_HARITH_ERR_NULL;

    /* Leaf cases */
    if (leaf_is_rk(node)) {
        g_stats.n_lowrank_skip++;
        return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
    }
    if (leaf_is_dense(node)) {
        int m = leaf_rows(node), n = leaf_cols(node);
        if (m != n) {
            /* should not happen for a diagonal leaf */
            return CHACAPK_HARITH_ERR_TOPOLOGY;
        }
        return dense_lu_inplace(leaf_dense_data(node), m);
    }

    /* Internal node: square block structure required for LU */
    int s = node->nrsons;
    if (s <= 0 || node->ncsons != s) return CHACAPK_HARITH_ERR_TOPOLOGY;

    for (int i = 0; i < s; i++) {
        /* (1) factor diagonal block A_ii = L_ii U_ii (in place). */
        st_cHACApK_block_node_t *Aii = node->sons[i + i * s];
        int rc = hlu_rec(Aii);
        if (rc != CHACAPK_HARITH_OK) return rc;

        /* (2) for k > i: solve U_ik from L_ii * U_ik = A_ik (htrsm_lln). */
        for (int k = i + 1; k < s; k++) {
            int rc2 = htrsm_lln(Aii, node->sons[i + k * s]);
            if (rc2 != CHACAPK_HARITH_OK) return rc2;
        }
        /* (3) for j > i: solve L_ji from L_ji * U_ii = A_ji (htrsm_run). */
        for (int j = i + 1; j < s; j++) {
            int rc2 = htrsm_run(Aii, node->sons[j + i * s]);
            if (rc2 != CHACAPK_HARITH_OK) return rc2;
        }
        /* (4) trailing update: A_jk -= L_ji * U_ik for j,k > i (h_addmul). */
        for (int k = i + 1; k < s; k++) {
            for (int j = i + 1; j < s; j++) {
                int rc2 = h_addmul(-1.0,
                                node->sons[j + i * s],
                                node->sons[i + k * s],
                                node->sons[j + k * s]);
                if (rc2 != CHACAPK_HARITH_OK) return rc2;
            }
        }
    }
    return CHACAPK_HARITH_OK;
}

int cHACApK_hlu_decomp(st_cHACApK_block_node_t *root)
{
    stats_reset();
    clear_ipiv_registry();
    clock_t t0 = clock();
    int rc = hlu_rec(root);
    g_stats.t_decomp_sec = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
    return rc;
}

/* ---------- solve_vec (post-decomp) ------------------------------- *
 * Solve A x = b where A's factors are stored in `root` from a successful
 * cHACApK_hlu_decomp. Two passes:
 *   forward:  L y = b   (with L unit-lower stored in lower triangle)
 *   backward: U x = y   (with U non-unit upper)
 *
 * Implementation note: we operate on slices of x corresponding to each
 * leaf's row range. The dense-leaf LU stored its own pivot in the side
 * registry; we apply it during the leaf forward-solve.
 */

static int hlu_forward_rec(
    const st_cHACApK_block_node_t *node,
    double *b /* length n_global, but only the [row_cluster] slice touched */)
{
    if (!node) return CHACAPK_HARITH_ERR_NULL;

    if (leaf_is_rk(node)) return CHACAPK_HARITH_ERR_LOWRANK_LEAF;

    /* For a leaf: must be a DIAGONAL leaf to handle here -- the caller
     * (internal recursion) handles off-diagonal updates explicitly. */
    if (leaf_is_dense(node)) {
        int n = leaf_rows(node);
        int row0 = node->row_cluster->nstrt - 1;  /* HACApK is 1-based */
        double *A = leaf_dense_data(node);
        /* no-pivot LU: just forward-substitute with unit-lower */
        cblas_dtrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasUnit,
                    n, A, n, &b[row0], 1);
        return CHACAPK_HARITH_OK;
    }

    /* Internal node: i = 0..s-1
     *   forward(sons[i+i*s], b)            // solve L_ii y_i = b_i (recursive)
     *   for j > i:  b_j -= L_ji * y_i      // update lower b (hmatvec_subtract)
     *
     * Aliasing note: b is both x and y to hmatvec_subtract. Safe because the
     * off-diagonal block L_ji has row_cluster in the j-range and col_cluster
     * in the i-range, which are disjoint (j > i). */
    int s = node->nrsons;
    for (int i = 0; i < s; i++) {
        st_cHACApK_block_node_t *Aii = node->sons[i + i * s];
        int rc = hlu_forward_rec(Aii, b);
        if (rc != CHACAPK_HARITH_OK) return rc;
        for (int j = i + 1; j < s; j++) {
            int rc2 = hmatvec_subtract(node->sons[j + i * s], b, b);
            if (rc2 != CHACAPK_HARITH_OK) return rc2;
        }
    }
    return CHACAPK_HARITH_OK;
}

static int hlu_backward_rec(
    const st_cHACApK_block_node_t *node,
    double *x /* in/out */)
{
    if (!node) return CHACAPK_HARITH_ERR_NULL;
    if (leaf_is_rk(node)) return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
    if (leaf_is_dense(node)) {
        int n = leaf_rows(node);
        int row0 = node->row_cluster->nstrt - 1;
        cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                    n, leaf_dense_data(node), n, &x[row0], 1);
        return CHACAPK_HARITH_OK;
    }
    /* Internal node: i = s-1..0 (reverse order)
     *   for k > i:  x_i -= U_ik * x_k     // update upper x (hmatvec_subtract)
     *   backward(sons[i+i*s], x)          // solve U_ii x_i = x_i (recursive)
     *
     * Aliasing: U_ik has row_cluster in the i-range, col_cluster in the
     * k-range, disjoint (k > i). Safe. */
    int s = node->nrsons;
    for (int i = s - 1; i >= 0; i--) {
        st_cHACApK_block_node_t *Aii = node->sons[i + i * s];
        for (int k = i + 1; k < s; k++) {
            int rc2 = hmatvec_subtract(node->sons[i + k * s], x, x);
            if (rc2 != CHACAPK_HARITH_OK) return rc2;
        }
        int rc = hlu_backward_rec(Aii, x);
        if (rc != CHACAPK_HARITH_OK) return rc;
    }
    return CHACAPK_HARITH_OK;
}

int cHACApK_hlu_solve_vec(
    const st_cHACApK_block_node_t *root,
    const double *b,
    double *x,
    int n)
{
    if (!root || !b || !x || n <= 0) return CHACAPK_HARITH_ERR_NULL;
    if (b != x) memcpy(x, b, sizeof(double) * (size_t)n);
    int rc = hlu_forward_rec(root, x);
    if (rc != CHACAPK_HARITH_OK) return rc;
    return hlu_backward_rec(root, x);
}

/* ---------- self-test: 2x2 dense-leaf block-tree ------------------ *
 * Build a synthetic n_total x n_total matrix (n_total = 2 * n_per_block),
 * split into 4 dense leaves arranged as a 2x2 block grid. Run H-LU + solve
 * and compare against a reference dense LAPACK solve.
 *
 * Memory: ~ 2*n_total^2 * 8 bytes (one copy for H-LU, one for the reference).
 * For n_per_block=100, n_total=200, memory ~640 KB. Negligible.
 */
/* Helpers for the deep self-test. The tree is built bottom-up: each call
 * to build_deep_tree creates either a leaf (depth==0) or an internal node
 * with 4 children that recursively build their own subtrees. The cluster
 * tree mirrors the block tree exactly (each node points to a cluster with
 * the same nstrt/nsize). All cluster + leaf + block_node allocations are
 * registered in two parallel arrays for clean cleanup. */
typedef struct deep_state_t {
    /* allocation log for clean cleanup */
    st_cHACApK_cluster_t  **clt_log;   int n_clt;   int cap_clt;
    st_cHACApK_leafmtx_t  **lf_log;    int n_lf;    int cap_lf;
    st_cHACApK_block_node_t **bn_log;  int n_bn;    int cap_bn;
    /* source data: column-major NxN matrix, N = level0_size * 2^depth */
    const double *A_full;
    int N_global;
} deep_state_t;

static st_cHACApK_cluster_t *log_clt(deep_state_t *s, st_cHACApK_cluster_t *c)
{
    if (s->n_clt >= s->cap_clt) {
        s->cap_clt = s->cap_clt ? 2*s->cap_clt : 64;
        s->clt_log = (st_cHACApK_cluster_t**)realloc(s->clt_log, sizeof(*s->clt_log)*s->cap_clt);
    }
    s->clt_log[s->n_clt++] = c; return c;
}
static st_cHACApK_leafmtx_t *log_lf(deep_state_t *s, st_cHACApK_leafmtx_t *l)
{
    if (s->n_lf >= s->cap_lf) {
        s->cap_lf = s->cap_lf ? 2*s->cap_lf : 64;
        s->lf_log = (st_cHACApK_leafmtx_t**)realloc(s->lf_log, sizeof(*s->lf_log)*s->cap_lf);
    }
    s->lf_log[s->n_lf++] = l; return l;
}
static st_cHACApK_block_node_t *log_bn(deep_state_t *s, st_cHACApK_block_node_t *b)
{
    if (s->n_bn >= s->cap_bn) {
        s->cap_bn = s->cap_bn ? 2*s->cap_bn : 64;
        s->bn_log = (st_cHACApK_block_node_t**)realloc(s->bn_log, sizeof(*s->bn_log)*s->cap_bn);
    }
    s->bn_log[s->n_bn++] = b; return b;
}

/* Build the recursive block-tree.
 *   row_clt, col_clt   already-allocated row/col cluster for this node
 *   depth              0 = leaf with size = row_clt->nsize x col_clt->nsize
 *                      >0 = split into 2x2 children
 * Returns a block_node. */
static st_cHACApK_block_node_t *build_deep_tree(
    deep_state_t *s,
    st_cHACApK_cluster_t *row_clt, st_cHACApK_cluster_t *col_clt,
    int depth)
{
    st_cHACApK_block_node_t *bn =
        (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn));
    log_bn(s, bn);
    bn->row_cluster = row_clt;
    bn->col_cluster = col_clt;

    if (depth == 0) {
        int m = row_clt->nsize, n = col_clt->nsize;
        st_cHACApK_leafmtx_t *L = (st_cHACApK_leafmtx_t*)calloc(1, sizeof(*L));
        log_lf(s, L);
        L->ltmtx = 2; L->ndl = m; L->ndt = n;
        L->nstrtl = row_clt->nstrt; L->nstrtt = col_clt->nstrt;
        L->a1 = (double*)malloc(sizeof(double)*(size_t)m*(size_t)n);
        /* copy from A_full (column-major NxN, with offset (nstrt-1, nstrt-1)) */
        int r0 = row_clt->nstrt - 1, c0 = col_clt->nstrt - 1, N = s->N_global;
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                L->a1[i + j*m] = s->A_full[(r0+i) + (c0+j)*N];
        bn->leaf_mtx  = L;
        bn->leaf_kind = 2;
        return bn;
    }

    /* internal: split into 2 row x 2 col children */
    int half_r = row_clt->nsize / 2;
    int half_c = col_clt->nsize / 2;
    st_cHACApK_cluster_t *rc[2] = {
        log_clt(s, (st_cHACApK_cluster_t*)calloc(1, sizeof(**rc))),
        log_clt(s, (st_cHACApK_cluster_t*)calloc(1, sizeof(**rc)))
    };
    st_cHACApK_cluster_t *cc[2] = {
        log_clt(s, (st_cHACApK_cluster_t*)calloc(1, sizeof(**cc))),
        log_clt(s, (st_cHACApK_cluster_t*)calloc(1, sizeof(**cc)))
    };
    rc[0]->nstrt = row_clt->nstrt;          rc[0]->nsize = half_r;
    rc[1]->nstrt = row_clt->nstrt + half_r; rc[1]->nsize = row_clt->nsize - half_r;
    cc[0]->nstrt = col_clt->nstrt;          cc[0]->nsize = half_c;
    cc[1]->nstrt = col_clt->nstrt + half_c; cc[1]->nsize = col_clt->nsize - half_c;

    bn->nrsons = 2; bn->ncsons = 2;
    bn->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*bn->sons));
    for (int j = 0; j < 2; j++)
        for (int i = 0; i < 2; i++)
            bn->sons[i + j*2] = build_deep_tree(s, rc[i], cc[j], depth - 1);
    return bn;
}

static void deep_cleanup(deep_state_t *s)
{
    for (int i = 0; i < s->n_bn; i++) {
        if (s->bn_log[i]->sons) free(s->bn_log[i]->sons);
        free(s->bn_log[i]);
    }
    for (int i = 0; i < s->n_lf; i++) { free(s->lf_log[i]->a1); free(s->lf_log[i]); }
    for (int i = 0; i < s->n_clt; i++) free(s->clt_log[i]);
    free(s->bn_log); free(s->lf_log); free(s->clt_log);
}

double cHACApK_harith_self_test(int depth, int n_per_block)
{
    if (n_per_block <= 0 || depth < 0 || depth > 10) return -1.0;
    int nb = n_per_block;
    int N = (1 << depth) * nb;

    /* (1) generate a random diagonally-dominant matrix (so LU without
     * pivoting is stable). Diagonal magnitude ~ N (sum of row magnitudes). */
    double *A_full = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *A_ref  = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *b      = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_hlu  = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_ref  = (double*)malloc(sizeof(double) * (size_t)N);
    if (!A_full || !A_ref || !b || !x_hlu || !x_ref) return -2.0;

    /* xorshift PRNG for reproducibility */
    unsigned long seed = 1234567UL;
    #define RND() (seed = seed * 6364136223846793005UL + 1442695040888963407UL, \
                   (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5)
    for (int j = 0; j < N; j++) {
        double rowsum = 0.0;
        for (int i = 0; i < N; i++) {
            double v = RND();
            A_full[i + j * N] = v;
            if (i != j) rowsum += (v < 0.0 ? -v : v);
        }
        /* boost diagonal to dominate the row */
        A_full[j + j * N] = rowsum + (double)N;
    }
    for (int j = 0; j < N; j++) b[j] = RND() * 10.0;
    memcpy(A_ref, A_full, sizeof(double) * (size_t)N * (size_t)N);

    /* (2) reference: LAPACK dgesv */
    int *ref_ipiv = (int*)malloc(sizeof(int) * (size_t)N);
    memcpy(x_ref, b, sizeof(double) * (size_t)N);
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A_ref, N, ref_ipiv, x_ref, N);
    free(ref_ipiv);
    if (info != 0) {
        free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
        return -3.0;
    }

    /* (3) build recursive block-tree of depth `depth` via build_deep_tree */
    deep_state_t st;
    memset(&st, 0, sizeof(st));
    st.A_full = A_full;
    st.N_global = N;
    st_cHACApK_cluster_t *root_clt =
        log_clt(&st, (st_cHACApK_cluster_t*)calloc(1, sizeof(*root_clt)));
    root_clt->nstrt = 1; root_clt->nsize = N;
    st_cHACApK_block_node_t *root_node =
        build_deep_tree(&st, root_clt, root_clt, depth);

    /* (4) run H-LU + solve */
    int rc_dec = cHACApK_hlu_decomp(root_node);
    int rc_slv = (rc_dec == CHACAPK_HARITH_OK)
                  ? cHACApK_hlu_solve_vec(root_node, b, x_hlu, N)
                  : -99;

    /* (5) compute max relative error */
    double max_ref = 0.0, max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double r = (x_ref[i] < 0.0) ? -x_ref[i] : x_ref[i];
        if (r > max_ref) max_ref = r;
        double d = x_hlu[i] - x_ref[i]; if (d < 0.0) d = -d;
        if (d > max_err) max_err = d;
    }
    double rel = (max_ref > 0.0) ? (max_err / max_ref) : max_err;

    /* (6) cleanup */
    deep_cleanup(&st);
    clear_ipiv_registry();
    free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);

    if (rc_dec != CHACAPK_HARITH_OK) return -4.0 + (double)rc_dec * 0.001;
    if (rc_slv != CHACAPK_HARITH_OK) return -5.0 + (double)rc_slv * 0.001;
    return rel;
}


/* ---------- Phase 3.5 unit test: h_addmul rk-rk -> rk + recompression ---- *
 * Direct test of the rk(A) * rk(B) -> rk(C) path. Builds three rk leaves
 * with explicit random U/V factors of given ranks, computes the dense
 * ground truth M_truth = U_c V_c^T + alpha * (U_a V_a^T) (U_b V_b^T),
 * calls h_addmul(alpha, A, B, C), and compares C's new dense reconstruction
 * U_c_new V_c_new^T against M_truth.
 *
 * Returns max element-wise relative error. Expected ~ 1e-13 (rounding from
 * the QR + SVD recompression). */
double cHACApK_harith_self_test_addmul_rkrk(int m, int n, int inner,
                                              int kA, int kB, int kC)
{
    if (m <= 0 || n <= 0 || inner <= 0) return -1.0;
    if (kA <= 0 || kB <= 0 || kC <= 0) return -1.0;
    double alpha = -1.0;

    /* PRNG */
    unsigned long seed = 99887766UL;
    #define RND() (seed = seed * 6364136223846793005UL + 1442695040888963407UL, \
                   (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5)

    /* Allocate U/V factors. */
    double *U_a = (double*)malloc(sizeof(double)*(size_t)m*(size_t)kA);
    double *V_a = (double*)malloc(sizeof(double)*(size_t)inner*(size_t)kA);
    double *U_b = (double*)malloc(sizeof(double)*(size_t)inner*(size_t)kB);
    double *V_b = (double*)malloc(sizeof(double)*(size_t)n*(size_t)kB);
    double *U_c = (double*)malloc(sizeof(double)*(size_t)m*(size_t)kC);
    double *V_c = (double*)malloc(sizeof(double)*(size_t)n*(size_t)kC);
    if (!U_a || !V_a || !U_b || !V_b || !U_c || !V_c) {
        free(U_a); free(V_a); free(U_b); free(V_b); free(U_c); free(V_c);
        return -2.0;
    }
    for (int i = 0; i < m*kA; i++)     U_a[i] = RND();
    for (int i = 0; i < inner*kA; i++) V_a[i] = RND();
    for (int i = 0; i < inner*kB; i++) U_b[i] = RND();
    for (int i = 0; i < n*kB; i++)     V_b[i] = RND();
    for (int i = 0; i < m*kC; i++)     U_c[i] = RND();
    for (int i = 0; i < n*kC; i++)     V_c[i] = RND();

    /* Dense ground truth:
     *   A_dense = U_a V_a^T              (m x inner)
     *   B_dense = U_b V_b^T              (inner x n)
     *   M_truth = U_c V_c^T + alpha A_dense B_dense   (m x n) */
    double *A_dense = (double*)malloc(sizeof(double)*(size_t)m*(size_t)inner);
    double *B_dense = (double*)malloc(sizeof(double)*(size_t)inner*(size_t)n);
    double *M_truth = (double*)malloc(sizeof(double)*(size_t)m*(size_t)n);
    if (!A_dense || !B_dense || !M_truth) {
        free(U_a); free(V_a); free(U_b); free(V_b); free(U_c); free(V_c);
        free(A_dense); free(B_dense); free(M_truth); return -2.0;
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                m, inner, kA, 1.0, U_a, m, V_a, inner, 0.0, A_dense, m);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                inner, n, kB, 1.0, U_b, inner, V_b, n, 0.0, B_dense, inner);
    /* M_truth = U_c V_c^T */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                m, n, kC, 1.0, U_c, m, V_c, n, 0.0, M_truth, m);
    /* M_truth += alpha * A_dense * B_dense */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                m, n, inner, alpha, A_dense, m, B_dense, inner, 1.0, M_truth, m);
    free(A_dense); free(B_dense);

    /* Build the three leafmtx + block_node wrappers. */
    st_cHACApK_cluster_t *clt_m  = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_m));
    st_cHACApK_cluster_t *clt_in = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_in));
    st_cHACApK_cluster_t *clt_n  = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_n));
    clt_m->nstrt = 1;  clt_m->nsize = m;
    clt_in->nstrt = 1; clt_in->nsize = inner;
    clt_n->nstrt = 1;  clt_n->nsize = n;

    st_cHACApK_leafmtx_t lf_A = {0}, lf_B = {0}, lf_C = {0};
    lf_A.ltmtx = 1; lf_A.kt = kA; lf_A.ndl = m;     lf_A.ndt = inner;
    lf_A.a1 = (double*)malloc(sizeof(double)*(size_t)m*(size_t)kA);
    lf_A.a2 = (double*)malloc(sizeof(double)*(size_t)inner*(size_t)kA);
    memcpy(lf_A.a1, U_a, sizeof(double)*(size_t)m*(size_t)kA);
    memcpy(lf_A.a2, V_a, sizeof(double)*(size_t)inner*(size_t)kA);

    lf_B.ltmtx = 1; lf_B.kt = kB; lf_B.ndl = inner; lf_B.ndt = n;
    lf_B.a1 = (double*)malloc(sizeof(double)*(size_t)inner*(size_t)kB);
    lf_B.a2 = (double*)malloc(sizeof(double)*(size_t)n*(size_t)kB);
    memcpy(lf_B.a1, U_b, sizeof(double)*(size_t)inner*(size_t)kB);
    memcpy(lf_B.a2, V_b, sizeof(double)*(size_t)n*(size_t)kB);

    lf_C.ltmtx = 1; lf_C.kt = kC; lf_C.ndl = m;     lf_C.ndt = n;
    lf_C.a1 = (double*)malloc(sizeof(double)*(size_t)m*(size_t)kC);
    lf_C.a2 = (double*)malloc(sizeof(double)*(size_t)n*(size_t)kC);
    memcpy(lf_C.a1, U_c, sizeof(double)*(size_t)m*(size_t)kC);
    memcpy(lf_C.a2, V_c, sizeof(double)*(size_t)n*(size_t)kC);

    st_cHACApK_block_node_t bn_A = {0}, bn_B = {0}, bn_C = {0};
    bn_A.row_cluster = clt_m;  bn_A.col_cluster = clt_in; bn_A.leaf_mtx = &lf_A; bn_A.leaf_kind = 1;
    bn_B.row_cluster = clt_in; bn_B.col_cluster = clt_n;  bn_B.leaf_mtx = &lf_B; bn_B.leaf_kind = 1;
    bn_C.row_cluster = clt_m;  bn_C.col_cluster = clt_n;  bn_C.leaf_mtx = &lf_C; bn_C.leaf_kind = 1;

    /* Run the test. */
    int rc = h_addmul(alpha, &bn_A, &bn_B, &bn_C);
    if (rc != CHACAPK_HARITH_OK) {
        free(lf_A.a1); free(lf_A.a2);
        free(lf_B.a1); free(lf_B.a2);
        free(lf_C.a1); free(lf_C.a2);
        free(U_a); free(V_a); free(U_b); free(V_b); free(U_c); free(V_c);
        free(M_truth); free(clt_m); free(clt_in); free(clt_n);
        return -4.0 + rc * 0.001;
    }

    /* Reconstruct M_new from C's updated factors. */
    int kC_new = lf_C.kt;
    double *M_new = (double*)malloc(sizeof(double)*(size_t)m*(size_t)n);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                m, n, kC_new, 1.0, lf_C.a1, m, lf_C.a2, n, 0.0, M_new, m);

    /* Max relative error vs M_truth. */
    double max_truth = 0.0, max_err = 0.0;
    for (int i = 0; i < m*n; i++) {
        double t = M_truth[i]; if (t < 0.0) t = -t;
        if (t > max_truth) max_truth = t;
        double e = M_truth[i] - M_new[i]; if (e < 0.0) e = -e;
        if (e > max_err) max_err = e;
    }
    double rel = (max_truth > 0.0) ? (max_err / max_truth) : max_err;

    free(M_new); free(M_truth);
    free(lf_A.a1); free(lf_A.a2);
    free(lf_B.a1); free(lf_B.a2);
    free(lf_C.a1); free(lf_C.a2);
    free(U_a); free(V_a); free(U_b); free(V_b); free(U_c); free(V_c);
    free(clt_m); free(clt_in); free(clt_n);
    return rel;
}


/* ---------- rk-aware self-test (Phase 3 partial validation) ------------- *
 * Build a depth=1 (2x2) block-tree where:
 *   A00, A11 are DENSE leaves (random diagonally-dominant)
 *   A01, A10 are RK leaves of explicit rank rk_rank (A_ij = U_ij V_ij^T)
 *
 * Exercises:
 *   - dense LU on diagonal leaves            (existing)
 *   - htrsm_lln(L=dense, X=rk)               (Phase 3 partial)
 *   - htrsm_run(U=dense, X=rk)               (Phase 3 partial)
 *   - h_addmul(rk*rk -> dense)               (Phase 3 partial)
 *   - hmatvec_subtract on rk leaves          (Phase 3 partial)
 *
 * Depth >= 2 with rk off-diagonals would also exercise
 * rk(A)*rk(B)->rk(C) in the trailing update, which requires ACA
 * recompression (Phase 3.5) and is NOT yet implemented. */
double cHACApK_harith_self_test_rk(int n_per_block, int rk_rank)
{
    if (n_per_block <= 0 || rk_rank <= 0) return -1.0;
    int nb = n_per_block;
    int N = 2 * nb;
    int k = (rk_rank > nb) ? nb : rk_rank;

    /* Allocate workspace. */
    double *A_full = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *A_ref  = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *b      = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_hlu  = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_ref  = (double*)malloc(sizeof(double) * (size_t)N);
    double *U01    = (double*)malloc(sizeof(double) * (size_t)nb * (size_t)k);
    double *V01    = (double*)malloc(sizeof(double) * (size_t)nb * (size_t)k);
    double *U10    = (double*)malloc(sizeof(double) * (size_t)nb * (size_t)k);
    double *V10    = (double*)malloc(sizeof(double) * (size_t)nb * (size_t)k);
    if (!A_full || !A_ref || !b || !x_hlu || !x_ref ||
        !U01 || !V01 || !U10 || !V10) return -2.0;

    /* xorshift PRNG (reproducible). */
    unsigned long seed = 7654321UL;
    #define RND() (seed = seed * 6364136223846793005UL + 1442695040888963407UL, \
                   (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5)

    /* (1) build the rk factors for A01 and A10. */
    for (int j = 0; j < k; j++)
        for (int i = 0; i < nb; i++) {
            U01[i + j*nb] = RND();
            V01[i + j*nb] = RND();
            U10[i + j*nb] = RND();
            V10[i + j*nb] = RND();
        }

    /* (2) fill A_full.
     *   - A00, A11: random + boost diag.
     *   - A01: U01 @ V01^T.
     *   - A10: U10 @ V10^T.
     * Column-major A_full of size NxN. */
    /* A00 random */
    for (int j = 0; j < nb; j++)
        for (int i = 0; i < nb; i++)
            A_full[i + j*N] = RND();
    /* A11 random */
    for (int j = 0; j < nb; j++)
        for (int i = 0; i < nb; i++)
            A_full[(nb+i) + (nb+j)*N] = RND();
    /* A01 = U01 V01^T  (rows 0..nb-1, cols nb..2nb-1) */
    for (int j = 0; j < nb; j++) {
        for (int i = 0; i < nb; i++) {
            double s = 0.0;
            for (int p = 0; p < k; p++)
                s += U01[i + p*nb] * V01[j + p*nb];
            A_full[i + (nb+j)*N] = s;
        }
    }
    /* A10 = U10 V10^T  (rows nb..2nb-1, cols 0..nb-1) */
    for (int j = 0; j < nb; j++) {
        for (int i = 0; i < nb; i++) {
            double s = 0.0;
            for (int p = 0; p < k; p++)
                s += U10[i + p*nb] * V10[j + p*nb];
            A_full[(nb+i) + j*N] = s;
        }
    }
    /* Boost diagonal to dominate the row (now that off-diag rk is known). */
    for (int i = 0; i < N; i++) {
        double rowsum = 0.0;
        for (int j = 0; j < N; j++) {
            if (j == i) continue;
            double v = A_full[i + j*N];
            rowsum += (v < 0.0) ? -v : v;
        }
        A_full[i + i*N] = rowsum + (double)N;
    }
    /* RHS */
    for (int j = 0; j < N; j++) b[j] = RND() * 10.0;
    memcpy(A_ref, A_full, sizeof(double) * (size_t)N * (size_t)N);

    /* (3) reference: dgesv on dense A_full */
    int *ref_ipiv = (int*)malloc(sizeof(int) * (size_t)N);
    memcpy(x_ref, b, sizeof(double) * (size_t)N);
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A_ref, N, ref_ipiv, x_ref, N);
    free(ref_ipiv);
    if (info != 0) {
        free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
        free(U01); free(V01); free(U10); free(V10);
        return -3.0;
    }

    /* (4) build the tree manually:  2x2 root with dense diag + rk off-diag. */
    st_cHACApK_cluster_t *clt_root = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_root));
    st_cHACApK_cluster_t *clt0     = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt0));
    st_cHACApK_cluster_t *clt1     = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt1));
    clt_root->nstrt = 1;     clt_root->nsize = N;
    clt0->nstrt     = 1;     clt0->nsize     = nb;
    clt1->nstrt     = 1+nb;  clt1->nsize     = nb;

    st_cHACApK_block_node_t *root = (st_cHACApK_block_node_t*)calloc(1, sizeof(*root));
    root->row_cluster = clt_root;
    root->col_cluster = clt_root;
    root->nrsons = 2; root->ncsons = 2;
    root->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*root->sons));

    /* Leaf at (i_son, j_son), i,j in {0,1}. row_cluster = clt0 if i==0 else clt1, etc. */
    st_cHACApK_block_node_t *nodes[4];
    st_cHACApK_leafmtx_t    *leaves[4];
    st_cHACApK_cluster_t    *rclt[2] = {clt0, clt1};
    st_cHACApK_cluster_t    *cclt[2] = {clt0, clt1};
    for (int j_son = 0; j_son < 2; j_son++) {
        for (int i_son = 0; i_son < 2; i_son++) {
            int idx = i_son + j_son * 2;
            st_cHACApK_block_node_t *bn = (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn));
            bn->row_cluster = rclt[i_son];
            bn->col_cluster = cclt[j_son];
            st_cHACApK_leafmtx_t *lf = (st_cHACApK_leafmtx_t*)calloc(1, sizeof(*lf));
            lf->ndl = nb; lf->ndt = nb;
            lf->nstrtl = bn->row_cluster->nstrt;
            lf->nstrtt = bn->col_cluster->nstrt;
            if (i_son == j_son) {
                /* dense diagonal */
                lf->ltmtx = 2;
                lf->a1 = (double*)malloc(sizeof(double) * (size_t)nb * (size_t)nb);
                int r0 = bn->row_cluster->nstrt - 1;
                int c0 = bn->col_cluster->nstrt - 1;
                for (int jj = 0; jj < nb; jj++)
                    for (int ii = 0; ii < nb; ii++)
                        lf->a1[ii + jj*nb] = A_full[(r0+ii) + (c0+jj)*N];
                bn->leaf_kind = 2;
            } else {
                /* rk off-diagonal: a1 = U, a2 = V */
                lf->ltmtx = 1;
                lf->kt    = k;
                lf->a1 = (double*)malloc(sizeof(double) * (size_t)nb * (size_t)k);
                lf->a2 = (double*)malloc(sizeof(double) * (size_t)nb * (size_t)k);
                /* identify which rk factors to use */
                double *Usrc = (i_son == 0 && j_son == 1) ? U01 : U10;
                double *Vsrc = (i_son == 0 && j_son == 1) ? V01 : V10;
                memcpy(lf->a1, Usrc, sizeof(double) * (size_t)nb * (size_t)k);
                memcpy(lf->a2, Vsrc, sizeof(double) * (size_t)nb * (size_t)k);
                bn->leaf_kind = 1;
            }
            bn->leaf_mtx = lf;
            nodes[idx]   = bn;
            leaves[idx]  = lf;
            root->sons[idx] = bn;
        }
    }

    /* (5) run H-LU + solve */
    int rc_dec = cHACApK_hlu_decomp(root);
    int rc_slv = (rc_dec == CHACAPK_HARITH_OK)
                  ? cHACApK_hlu_solve_vec(root, b, x_hlu, N)
                  : -99;

    /* (6) max rel err */
    double max_ref = 0.0, max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double r = (x_ref[i] < 0.0) ? -x_ref[i] : x_ref[i];
        if (r > max_ref) max_ref = r;
        double d = x_hlu[i] - x_ref[i]; if (d < 0.0) d = -d;
        if (d > max_err) max_err = d;
    }
    double rel = (max_ref > 0.0) ? (max_err / max_ref) : max_err;

    /* (7) cleanup */
    for (int i = 0; i < 4; i++) {
        free(leaves[i]->a1);
        if (leaves[i]->a2) free(leaves[i]->a2);
        free(leaves[i]);
        free(nodes[i]);
    }
    free(root->sons); free(root);
    free(clt_root); free(clt0); free(clt1);
    free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
    free(U01); free(V01); free(U10); free(V10);

    if (rc_dec != CHACAPK_HARITH_OK) return -4.0 + (double)rc_dec * 0.001;
    if (rc_slv != CHACAPK_HARITH_OK) return -5.0 + (double)rc_slv * 0.001;
    return rel;
}


/* ---------- Phase 4 driver: H-LU on a real HACApK tree ----------- *
 *
 * The caller has already built a HACApK leafmtxp (in HACApK row-major
 * dense / V-first rk format) and computed y_orig = A * x_orig in the
 * user's original ordering. This driver:
 *
 *   1. Permutes x_orig, y_orig to HACApK internal ordering via control->lod.
 *   2. Converts the leafmtxp leaves to internal column-major format.
 *   3. Builds the block-tree view from leafmtxp + cluster_root.
 *   4. Runs cHACApK_hlu_decomp.
 *   5. Solves x_internal = A^-1 y_perm via cHACApK_hlu_solve_vec.
 *   6. Inverse-permutes x_internal back to original ordering.
 *   7. Compares against x_orig.
 *
 * The leafmtxp leaves now contain the LU factors (in internal layout).
 * The H-matrix is consumed -- can't be used for further MatVec. */
double cHACApK_hlu_run_on_hacapk(void *leafmtxp_void, void *control_void,
                                  const double *x_orig, const double *y_orig,
                                  int nffc)
{
    if (!leafmtxp_void || !control_void || !x_orig || !y_orig) return -1.0;
    if (nffc <= 0) return -1.0;
    st_cHACApK_leafmtxp_t *lp = (st_cHACApK_leafmtxp_t*)leafmtxp_void;
    st_cHACApK_lcontrol_t *lc = (st_cHACApK_lcontrol_t*)control_void;
    int nd = lp->nd;
    if (nd <= 0) return -1.0;
    if (!lp->st_clt_root) return -1.0;
    int *lod = lc->lod;
    if (!lod) return -1.0;

    /* (1) permute x_orig, y_orig to HACApK internal ordering: x_perm[i] = x_orig[lod[i+1]-1]. */
    double *x_perm = (double*)malloc(sizeof(double) * (size_t)nd);
    double *y_perm = (double*)malloc(sizeof(double) * (size_t)nd);
    double *x_solved = (double*)malloc(sizeof(double) * (size_t)nd);
    double *x_back = (double*)malloc(sizeof(double) * (size_t)nd);
    if (!x_perm || !y_perm || !x_solved || !x_back) {
        free(x_perm); free(y_perm); free(x_solved); free(x_back);
        return -2.0;
    }
    for (int i = 0; i < nd; i++) {
        int j = lod[i + 1] - 1;
        x_perm[i] = x_orig[j];
        y_perm[i] = y_orig[j];
    }

    /* (2) convert HACApK leaves to internal format. */
    cHACApK_convert_leafmtxp_to_internal(lp);

    /* (3) build block-tree view -- nffc bridges cluster (element units)
     * to leaf (DOF units). */
    st_cHACApK_block_node_t *root = cHACApK_build_block_tree_nffc(
        lp, lp->st_clt_root, lp->st_clt_root, nffc);
    if (!root) {
        free(x_perm); free(y_perm); free(x_solved); free(x_back);
        return -3.0;
    }

    /* (4) factor. */
    int rc_dec = cHACApK_hlu_decomp(root);
    /* (5) solve. */
    int rc_slv = CHACAPK_HARITH_OK;
    if (rc_dec == CHACAPK_HARITH_OK) {
        rc_slv = cHACApK_hlu_solve_vec(root, y_perm, x_solved, nd);
    }

    double rel = 0.0;
    if (rc_dec == CHACAPK_HARITH_OK && rc_slv == CHACAPK_HARITH_OK) {
        /* (6) inverse-permute x_solved back to original ordering. */
        for (int i = 0; i < nd; i++) {
            int j = lod[i + 1] - 1;
            x_back[j] = x_solved[i];
        }
        /* (7) compute max rel err. */
        double max_x = 0.0, max_err = 0.0;
        for (int i = 0; i < nd; i++) {
            double xa = (x_orig[i] < 0.0) ? -x_orig[i] : x_orig[i];
            if (xa > max_x) max_x = xa;
            double d = x_orig[i] - x_back[i]; if (d < 0.0) d = -d;
            if (d > max_err) max_err = d;
        }
        rel = (max_x > 0.0) ? (max_err / max_x) : max_err;
    }

    cHACApK_free_block_tree(root);
    free(x_perm); free(y_perm); free(x_solved); free(x_back);

    if (rc_dec != CHACAPK_HARITH_OK) return -4.0 + (double)rc_dec * 0.001;
    if (rc_slv != CHACAPK_HARITH_OK) return -5.0 + (double)rc_slv * 0.001;
    return rel;
}


/* ---------- depth=2 rk self-test (full Phase 1-3.5 integration) --------- *
 *
 * Builds a depth-2 block-tree (4x4 leaf grid) where:
 *   GLOBAL diagonal leaves (A00, A11, A22, A33) are DENSE
 *   12 GLOBAL off-diagonal leaves are RK of rank rk_rank
 *
 * Exercises through the recursive LU:
 *   Phase 1: dense LU on diagonal leaves
 *   Phase 2: recursive H-arith descent into internal nodes
 *   Phase 3 partial: htrsm with dense L and rk X
 *                    h_addmul rk*rk -> dense (when target is a diagonal leaf)
 *                    hmatvec_subtract on rk leaves
 *   Phase 3.5:       h_addmul rk*rk -> rk + recompression
 *                    (when target is an off-diagonal leaf inside a sub-block) */
double cHACApK_harith_self_test_rk_deep(int n_per_block, int rk_rank)
{
    if (n_per_block <= 0 || rk_rank <= 0) return -1.0;
    int nb = n_per_block;
    int N = 4 * nb;
    int k = (rk_rank > nb) ? nb : rk_rank;

    unsigned long seed = 333555777UL;
    #define RND2() (seed = seed * 6364136223846793005UL + 1442695040888963407UL, \
                    (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5)

    double *A_full = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *A_ref  = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *b      = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_hlu  = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_ref  = (double*)malloc(sizeof(double) * (size_t)N);
    if (!A_full || !A_ref || !b || !x_hlu || !x_ref) return -2.0;
    memset(A_full, 0, sizeof(double) * (size_t)N * (size_t)N);

    /* 12 rk factor pairs for off-diagonal leaves. */
    int n_offdiag = 12;
    double *U_off = (double*)malloc(sizeof(double) * (size_t)n_offdiag * (size_t)nb * (size_t)k);
    double *V_off = (double*)malloc(sizeof(double) * (size_t)n_offdiag * (size_t)nb * (size_t)k);
    int    off_ij[12][2];
    if (!U_off || !V_off) return -2.0;
    {
        int idx = 0;
        for (int j_grid = 0; j_grid < 4; j_grid++)
            for (int i_grid = 0; i_grid < 4; i_grid++)
                if (i_grid != j_grid) {
                    off_ij[idx][0] = i_grid;
                    off_ij[idx][1] = j_grid;
                    idx++;
                }
    }
    for (int idx = 0; idx < n_offdiag*nb*k; idx++) U_off[idx] = RND2();
    for (int idx = 0; idx < n_offdiag*nb*k; idx++) V_off[idx] = RND2();

    /* Fill A_full off-diag from rk factors. */
    for (int t = 0; t < n_offdiag; t++) {
        int i_grid = off_ij[t][0], j_grid = off_ij[t][1];
        int r0 = i_grid * nb, c0 = j_grid * nb;
        double *Ut = U_off + (size_t)t * (size_t)nb * (size_t)k;
        double *Vt = V_off + (size_t)t * (size_t)nb * (size_t)k;
        for (int jj = 0; jj < nb; jj++) {
            for (int ii = 0; ii < nb; ii++) {
                double s = 0.0;
                for (int p = 0; p < k; p++) s += Ut[ii + p*nb] * Vt[jj + p*nb];
                A_full[(r0+ii) + (c0+jj)*N] = s;
            }
        }
    }
    /* Fill A_full diagonal (random dense). */
    for (int g = 0; g < 4; g++) {
        int d0 = g * nb;
        for (int jj = 0; jj < nb; jj++)
            for (int ii = 0; ii < nb; ii++)
                A_full[(d0+ii) + (d0+jj)*N] = RND2();
    }
    /* Diag-dominate. */
    for (int i = 0; i < N; i++) {
        double rs = 0.0;
        for (int j = 0; j < N; j++) {
            if (j == i) continue;
            double v = A_full[i + j*N];
            rs += (v < 0.0) ? -v : v;
        }
        A_full[i + i*N] = rs + (double)N;
    }
    for (int j = 0; j < N; j++) b[j] = RND2() * 10.0;
    memcpy(A_ref, A_full, sizeof(double) * (size_t)N * (size_t)N);

    /* Reference dgesv. */
    int *ipiv_ref = (int*)malloc(sizeof(int) * (size_t)N);
    memcpy(x_ref, b, sizeof(double) * (size_t)N);
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A_ref, N, ipiv_ref, x_ref, N);
    free(ipiv_ref);
    if (info != 0) {
        free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
        free(U_off); free(V_off);
        return -3.0;
    }

    /* Cluster tree: root + 2 lvl-1 + 4 lvl-2. */
    st_cHACApK_cluster_t *clt_root = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_root));
    st_cHACApK_cluster_t *clt_lvl1[2];
    st_cHACApK_cluster_t *clt_lvl2[4];
    clt_root->nstrt = 1; clt_root->nsize = N;
    for (int g = 0; g < 2; g++) {
        clt_lvl1[g] = (st_cHACApK_cluster_t*)calloc(1, sizeof(**clt_lvl1));
        clt_lvl1[g]->nstrt = 1 + g*2*nb;
        clt_lvl1[g]->nsize = 2*nb;
    }
    for (int g = 0; g < 4; g++) {
        clt_lvl2[g] = (st_cHACApK_cluster_t*)calloc(1, sizeof(**clt_lvl2));
        clt_lvl2[g]->nstrt = 1 + g*nb;
        clt_lvl2[g]->nsize = nb;
    }

    /* Build 16 leaves + 4 lvl-1 internal + 1 root. */
    st_cHACApK_leafmtx_t    *all_lf[16]; int n_lf = 0;
    st_cHACApK_block_node_t *all_bn[21]; int n_bn = 0;
    st_cHACApK_block_node_t *lvl1_nodes[4];

    for (int j_root = 0; j_root < 2; j_root++) {
        for (int i_root = 0; i_root < 2; i_root++) {
            st_cHACApK_block_node_t *bn1 = (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn1));
            bn1->row_cluster = clt_lvl1[i_root];
            bn1->col_cluster = clt_lvl1[j_root];
            bn1->nrsons = 2; bn1->ncsons = 2;
            bn1->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*bn1->sons));
            all_bn[n_bn++] = bn1;

            for (int j_sub = 0; j_sub < 2; j_sub++) {
                for (int i_sub = 0; i_sub < 2; i_sub++) {
                    int i_grid = 2*i_root + i_sub;
                    int j_grid = 2*j_root + j_sub;
                    st_cHACApK_block_node_t *bn = (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn));
                    bn->row_cluster = clt_lvl2[i_grid];
                    bn->col_cluster = clt_lvl2[j_grid];
                    st_cHACApK_leafmtx_t *lf = (st_cHACApK_leafmtx_t*)calloc(1, sizeof(*lf));
                    lf->ndl = nb; lf->ndt = nb;
                    lf->nstrtl = bn->row_cluster->nstrt;
                    lf->nstrtt = bn->col_cluster->nstrt;
                    if (i_grid == j_grid) {
                        lf->ltmtx = 2;
                        lf->a1 = (double*)malloc(sizeof(double)*(size_t)nb*(size_t)nb);
                        int r0 = bn->row_cluster->nstrt - 1;
                        int c0 = bn->col_cluster->nstrt - 1;
                        for (int jj = 0; jj < nb; jj++)
                            for (int ii = 0; ii < nb; ii++)
                                lf->a1[ii + jj*nb] = A_full[(r0+ii) + (c0+jj)*N];
                        bn->leaf_kind = 2;
                    } else {
                        int slot = -1;
                        for (int t = 0; t < n_offdiag; t++) {
                            if (off_ij[t][0] == i_grid && off_ij[t][1] == j_grid) { slot = t; break; }
                        }
                        lf->ltmtx = 1; lf->kt = k;
                        lf->a1 = (double*)malloc(sizeof(double)*(size_t)nb*(size_t)k);
                        lf->a2 = (double*)malloc(sizeof(double)*(size_t)nb*(size_t)k);
                        memcpy(lf->a1, U_off + (size_t)slot*(size_t)nb*(size_t)k,
                               sizeof(double)*(size_t)nb*(size_t)k);
                        memcpy(lf->a2, V_off + (size_t)slot*(size_t)nb*(size_t)k,
                               sizeof(double)*(size_t)nb*(size_t)k);
                        bn->leaf_kind = 1;
                    }
                    bn->leaf_mtx = lf;
                    all_lf[n_lf++] = lf;
                    all_bn[n_bn++] = bn;
                    bn1->sons[i_sub + j_sub*2] = bn;
                }
            }
            lvl1_nodes[i_root + j_root*2] = bn1;
        }
    }

    st_cHACApK_block_node_t *root = (st_cHACApK_block_node_t*)calloc(1, sizeof(*root));
    root->row_cluster = clt_root;
    root->col_cluster = clt_root;
    root->nrsons = 2; root->ncsons = 2;
    root->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*root->sons));
    for (int idx = 0; idx < 4; idx++) root->sons[idx] = lvl1_nodes[idx];
    all_bn[n_bn++] = root;

    /* H-LU run. */
    int rc_dec = cHACApK_hlu_decomp(root);
    int rc_slv = (rc_dec == CHACAPK_HARITH_OK)
                  ? cHACApK_hlu_solve_vec(root, b, x_hlu, N)
                  : -99;

    /* Max rel err. */
    double max_ref = 0.0, max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double r = (x_ref[i] < 0.0) ? -x_ref[i] : x_ref[i];
        if (r > max_ref) max_ref = r;
        double d = x_hlu[i] - x_ref[i]; if (d < 0.0) d = -d;
        if (d > max_err) max_err = d;
    }
    double rel = (max_ref > 0.0) ? (max_err / max_ref) : max_err;

    /* Cleanup. */
    for (int i = 0; i < n_lf; i++) {
        free(all_lf[i]->a1);
        if (all_lf[i]->a2) free(all_lf[i]->a2);
        free(all_lf[i]);
    }
    for (int i = 0; i < n_bn; i++) {
        if (all_bn[i]->sons) free(all_bn[i]->sons);
        free(all_bn[i]);
    }
    free(clt_root);
    for (int i = 0; i < 2; i++) free(clt_lvl1[i]);
    for (int i = 0; i < 4; i++) free(clt_lvl2[i]);
    free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
    free(U_off); free(V_off);

    if (rc_dec != CHACAPK_HARITH_OK) return -4.0 + (double)rc_dec * 0.001;
    if (rc_slv != CHACAPK_HARITH_OK) return -5.0 + (double)rc_slv * 0.001;
    return rel;
}
