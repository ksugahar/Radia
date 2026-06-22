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

/* ---------- module-global recompression tolerance --------------- *
 *
 * Relative tolerance for rk-leaf SVD recompression during H-LU. The
 * KEY speed/accuracy knob: 1e-14 keeps machine-precision (ranks grow
 * toward full -> H-LU degenerates to dense + SVD overhead, O(N^3+));
 * 1e-4 (~ACA accuracy) keeps ranks LOW -> the intended O(N log^2 N)
 * scaling. Set via cHACApK_hlu_set_trunc_tol before a decomp.
 * Default 1e-14 preserves the high-accuracy synthetic self-tests. */
static double g_hlu_trunc_tol = 1e-14;

void cHACApK_hlu_set_trunc_tol(double tol)
{
    if (tol > 0.0 && tol < 1.0) g_hlu_trunc_tol = tol;
}
double cHACApK_hlu_get_trunc_tol(void) { return g_hlu_trunc_tol; }

/* ---------- block-level parallelism (ngcore TaskManager bridge) ---- *
 * Defined in cHACApK_harith_par.cpp. The H-LU leaf BLAS are too small for
 * MKL threading, so we parallelize the independent output blocks of the
 * recursive h_addmul descent. g_hlu_parallel toggles it (default on);
 * g_hlu_par_cutoff is the minimum C-block area (rows*cols) below which a
 * block runs serial to avoid tiny-task overhead. */
extern int  chacapk_max_threads(void);
extern void chacapk_par_region(void (*body)(void*), void *ctx);
extern void chacapk_par_for(int n, void (*body)(int, void*), void *ctx);

static int  g_hlu_parallel   = 1;
static long g_hlu_par_cutoff = 250000;  /* ~500x500 block */

void cHACApK_hlu_set_parallel(int on) { g_hlu_parallel = on ? 1 : 0; }
int  cHACApK_hlu_get_parallel(void)   { return g_hlu_parallel; }
void cHACApK_hlu_set_par_cutoff(long c){ if (c > 0) g_hlu_par_cutoff = c; }

/* ---------- accumulator (lazy recompression, Borm/Kriemann PP18) --- *
 * Instead of recompressing an rk-leaf on EVERY low-rank update (a QR-QR-SVD
 * per update), append the increment columns to the leaf's a1/a2 (the rk
 * representation U V^T stays mathematically EXACT as columns grow, so the
 * leaf can still be read mid-factorization), and recompress ONCE when the
 * accumulated rank exceeds g_hlu_accum_cap, plus a final flush pass after
 * decomp. This cuts the number of truncations (the accumulator idea), at the
 * cost of higher-rank intermediate ops. g_hlu_accum_cap = 0 disables it
 * (recompress every update = previous behavior) for A/B comparison. */
static int g_hlu_accum_cap = 64;
void cHACApK_hlu_set_accum_cap(int c) { g_hlu_accum_cap = (c < 0) ? 0 : c; }
int  cHACApK_hlu_get_accum_cap(void)  { return g_hlu_accum_cap; }

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

/* Flat accessor for pybind (avoids exposing the struct layout). */
extern long g_dbg_n_materialize, g_dbg_materialize_elems;
void cHACApK_hlu_get_timings(double *out_t_decomp, double *out_t_solve,
                              long *out_n_dense_lu, long *out_n_dense_gemm)
{
    if (out_t_decomp)   *out_t_decomp   = g_stats.t_decomp_sec;
    if (out_t_solve)    *out_t_solve    = g_stats.t_solve_sec;
    if (out_n_dense_lu) *out_n_dense_lu = g_stats.n_dense_lu;
    if (out_n_dense_gemm) *out_n_dense_gemm = g_stats.n_dense_gemm;
}

/* Profiling accessor for the mixed-case materialize fallback. */
void cHACApK_hlu_get_materialize_stats(long *out_n_calls, long *out_n_elems)
{
    if (out_n_calls) *out_n_calls = g_dbg_n_materialize;
    if (out_n_elems) *out_n_elems = g_dbg_materialize_elems;
}

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

    /* Economy QR ranks. When the stacked rank k_in exceeds the leaf
     * dimensions (common after rk-into-rk stacking: kw = kt_c + kt_inc
     * can exceed min(m,n)), the thin Q has only min(m,k_in) / min(n,k_in)
     * columns. Using k_in directly would make dorgqr fail with
     * "Parameter 2 incorrect" (n > m). */
    int ru = (m < k_in) ? m : k_in;   /* columns of Q_U / rows of R_U */
    int rv = (n < k_in) ? n : k_in;   /* columns of Q_V / rows of R_V */

    /* Workspace copies (LAPACK overwrites inputs). */
    double *U_work = (double*)malloc(sizeof(double) * (size_t)m * (size_t)k_in);
    double *V_work = (double*)malloc(sizeof(double) * (size_t)n * (size_t)k_in);
    double *R_U    = (double*)calloc((size_t)ru * (size_t)k_in, sizeof(double));
    double *R_V    = (double*)calloc((size_t)rv * (size_t)k_in, sizeof(double));
    double *tau    = (double*)malloc(sizeof(double) * (size_t)k_in);
    if (!U_work || !V_work || !R_U || !R_V || !tau) {
        free(U_work); free(V_work); free(R_U); free(R_V); free(tau);
        return CHACAPK_HARITH_ERR_NULL;
    }
    memcpy(U_work, U, sizeof(double) * (size_t)m * (size_t)k_in);
    memcpy(V_work, V, sizeof(double) * (size_t)n * (size_t)k_in);

    /* QR of U_work (m x k_in): ru = min(m,k_in) reflectors, R_U is ru x k_in
     * upper-trapezoidal in the top-left of U_work. */
    int info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m, k_in, U_work, m, tau);
    if (info != 0) goto lapack_err;
    for (int j = 0; j < k_in; j++)
        for (int i = 0; i <= j && i < ru; i++) R_U[i + j*ru] = U_work[i + j*m];
    info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, m, ru, ru, U_work, m, tau);
    if (info != 0) goto lapack_err;
    /* Now U_work[:, :ru] = Q_U (m x ru). */

    /* QR of V_work (n x k_in): rv = min(n,k_in) reflectors. */
    info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, n, k_in, V_work, n, tau);
    if (info != 0) goto lapack_err;
    for (int j = 0; j < k_in; j++)
        for (int i = 0; i <= j && i < rv; i++) R_V[i + j*rv] = V_work[i + j*n];
    info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, n, rv, rv, V_work, n, tau);
    if (info != 0) goto lapack_err;
    /* Now V_work[:, :rv] = Q_V (n x rv). */

    free(tau); tau = NULL;

    /* S = R_U R_V^T  (ru x rv). */
    double *S = (double*)malloc(sizeof(double) * (size_t)ru * (size_t)rv);
    if (!S) { free(U_work); free(V_work); free(R_U); free(R_V);
              return CHACAPK_HARITH_ERR_NULL; }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                ru, rv, k_in, 1.0,
                R_U, ru, R_V, rv, 0.0, S, ru);
    free(R_U); free(R_V);

    /* SVD: S = Us Sigma Vt.  p = min(ru, rv). */
    int p = (ru < rv) ? ru : rv;
    double *sigma  = (double*)malloc(sizeof(double) * (size_t)p);
    double *Us     = (double*)malloc(sizeof(double) * (size_t)ru * (size_t)p);
    double *Vt     = (double*)malloc(sizeof(double) * (size_t)p * (size_t)rv);
    double *superb = (double*)malloc(sizeof(double) * (size_t)(p > 1 ? p - 1 : 1));
    if (!sigma || !Us || !Vt || !superb) {
        free(sigma); free(Us); free(Vt); free(superb);
        free(S); free(U_work); free(V_work);
        return CHACAPK_HARITH_ERR_NULL;
    }
    info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'S', ru, rv, S, ru,
                          sigma, Us, ru, Vt, p, superb);
    free(superb); free(S);
    if (info != 0) { free(sigma); free(Us); free(Vt);
                     free(U_work); free(V_work);
                     return CHACAPK_HARITH_ERR_LAPACK; }

    /* Determine new rank (bounded by p = min(ru, rv) <= min(m, n)). */
    int k_new = 0;
    double sv_max = sigma[0];
    if (sv_max > 0.0) {
        double thresh = tol_rel * sv_max;
        for (int i = 0; i < p; i++) {
            if (sigma[i] > thresh) k_new++;
            else break;
        }
        if (k_new > k_max) k_new = k_max;
        if (k_new > p)     k_new = p;
        if (k_new < 1)     k_new = 1;
    } else {
        k_new = 1;  /* rank-1 to keep something */
    }

    /* U_new = Q_U Us[:,:k_new] diag(sigma[:k_new])  (m x k_new). */
    double *U_new = (double*)malloc(sizeof(double) * (size_t)m * (size_t)k_new);
    double *V_new = (double*)malloc(sizeof(double) * (size_t)n * (size_t)k_new);
    double *Us_scaled = (double*)malloc(sizeof(double) * (size_t)ru * (size_t)k_new);
    if (!U_new || !V_new || !Us_scaled) {
        free(U_new); free(V_new); free(Us_scaled);
        free(sigma); free(Us); free(Vt);
        free(U_work); free(V_work);
        return CHACAPK_HARITH_ERR_NULL;
    }
    for (int j = 0; j < k_new; j++) {
        double s = sigma[j];
        for (int i = 0; i < ru; i++)
            Us_scaled[i + j*ru] = Us[i + j*ru] * s;
    }
    /* U_new = Q_U (m x ru) * Us_scaled (ru x k_new). */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                m, k_new, ru, 1.0,
                U_work, m, Us_scaled, ru, 0.0, U_new, m);
    free(Us_scaled);

    /* V_new = Q_V Vs[:,:k_new], Vs = Vt^T (rv x p).  First k_new columns of
     * Vs = first k_new rows of Vt.  cblas Trans on Vt (ldb=p) accesses the
     * first k_new rows (B logically N x K = k_new x rv post-Trans). */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                n, k_new, rv, 1.0,
                V_work, n, Vt, p, 0.0, V_new, n);

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
/* DEBUG counters (profiling the mixed-case materialize fallback). */
long g_dbg_n_materialize = 0;
long g_dbg_materialize_elems = 0;
/* Split the materialize count by node kind: an INTERNAL-node densification is the
 * cubic-driving "materialize a whole subtree" cost (the thing the materialize-free
 * R(A*B) + htrsm-recursion fixes eliminate); a LEAF densification is the benign,
 * O(leaf^2) copy of a dense leaf's data for a leaf-level dgemm.  g_dbg_mat_internal
 * == 0 is the invariant that the near-cubic materialize is gone. */
long g_dbg_mat_internal = 0;
long g_dbg_mat_leaf = 0;
void cHACApK_hlu_get_materialize_split(long *out_internal, long *out_leaf)
{
    if (out_internal) *out_internal = g_dbg_mat_internal;
    if (out_leaf)     *out_leaf     = g_dbg_mat_leaf;
}

/* Breakdown of WHICH operand kinds trigger the materialize fallback, so the
 * optimization (rk-factored mixed multiply vs sub-view split) is targeted at
 * the dominant case rather than guessed. kind index: 0=internal, 1=rk leaf,
 * 2=dense leaf. addmul indexed [kindA*3 + kindB]; lln/run [kindL_or_U*3 + kindX]. */
long g_dbg_mixed_addmul[9] = {0};
long g_dbg_mixed_lln[9] = {0};
long g_dbg_mixed_run[9] = {0};

static inline int node_kind_idx(const st_cHACApK_block_node_t *n)
{
    if (leaf_is_rk(n)) return 1;
    if (leaf_is_dense(n)) return 2;
    return 0; /* internal */
}

void cHACApK_hlu_get_mixed_breakdown(long *out_addmul9, long *out_lln9, long *out_run9)
{
    if (out_addmul9) for (int i = 0; i < 9; i++) out_addmul9[i] = g_dbg_mixed_addmul[i];
    if (out_lln9)    for (int i = 0; i < 9; i++) out_lln9[i]    = g_dbg_mixed_lln[i];
    if (out_run9)    for (int i = 0; i < 9; i++) out_run9[i]    = g_dbg_mixed_run[i];
}

static double *materialize_node_as_dense(const st_cHACApK_block_node_t *node)
{
    if (!node) return NULL;
    int m = node->dof_nrows;
    int n = node->dof_ncols;
    g_dbg_n_materialize++;
    g_dbg_materialize_elems += (long)m * (long)n;
    if (node->nrsons > 0) g_dbg_mat_internal++; else g_dbg_mat_leaf++;
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
    int row_base = node->dof_row_start;
    int col_base = node->dof_col_start;
    for (int j_s = 0; j_s < nc; j_s++) {
        for (int i_s = 0; i_s < nr; i_s++) {
            const st_cHACApK_block_node_t *child = node->sons[i_s + j_s * nr];
            int cm = child->dof_nrows;
            int cn = child->dof_ncols;
            int row_off = child->dof_row_start - row_base;
            int col_off = child->dof_col_start - col_base;
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
        int rc = dense_to_rk_truncate(C_full, m, n, g_hlu_trunc_tol, kmin,
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
    int row_base = C->dof_row_start;
    int col_base = C->dof_col_start;
    for (int j_s = 0; j_s < nc; j_s++) {
        for (int i_s = 0; i_s < nr; i_s++) {
            st_cHACApK_block_node_t *child = C->sons[i_s + j_s * nr];
            int cm = child->dof_nrows;
            int cn = child->dof_ncols;
            int row_off = child->dof_row_start - row_base;
            int col_off = child->dof_col_start - col_base;
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
        int rc = dense_to_rk_truncate(D_use, m, n, g_hlu_trunc_tol, kmin,
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
    int row_base = C->dof_row_start;
    int col_base = C->dof_col_start;
    for (int j_s = 0; j_s < nc; j_s++) {
        for (int i_s = 0; i_s < nr; i_s++) {
            st_cHACApK_block_node_t *child = C->sons[i_s + j_s * nr];
            int cm = child->dof_nrows;
            int cn = child->dof_ncols;
            int row_off = child->dof_row_start - row_base;
            int col_off = child->dof_col_start - col_base;
            const double *D_sub = D + row_off + (size_t)col_off * (size_t)D_ld;
            int rc = set_node_from_dense(D_sub, D_ld, cm, cn, child);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
    }
    return CHACAPK_HARITH_OK;
}


/* ---------- matvec/matmat helpers for the materialize-free mixed multiply ----
 *
 * These apply a block-tree node (dense leaf / rk leaf / internal) to a dense
 * block of k columns WITHOUT ever densifying the node. They use LOCAL column-
 * major indexing: B is indexed [0, node->dof_ncols) x k, C is [0, node->dof_
 * nrows) x k (or transposed for the _trans variant). Internal nodes recurse,
 * slicing B/C by the child's DOF offset relative to the parent -- exactly the
 * pattern in materialize_node_as_dense, but the node stays compressed.
 *
 * This is the key to exploiting HACApK's compression: an rk-leaf operand in a
 * trailing update is handled as U(V^T B) via these helpers + add_lowrank_to_
 * node, instead of the Phase 3.6 materialize-and-redo fallback that densified
 * the whole internal operand (the dominant materialize cost at small leaf). */

/* C[dof_nrows x k] += alpha * node * B[dof_ncols x k]  (LOCAL indexing). */
static int node_matmat_local(const st_cHACApK_block_node_t *node, double alpha,
                             const double *B, int ldB, int k,
                             double *C, int ldC)
{
    if (!node || k <= 0) return CHACAPK_HARITH_OK;
    int m = node->dof_nrows, n = node->dof_ncols;

    if (leaf_is_dense(node)) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    m, k, n, alpha, leaf_dense_data(node), m,
                    B, ldB, 1.0, C, ldC);
        return CHACAPK_HARITH_OK;
    }
    if (leaf_is_rk(node)) {
        int kt = leaf_rk_rank(node);
        if (kt <= 0) return CHACAPK_HARITH_OK;
        double *T = (double*)malloc(sizeof(double) * (size_t)kt * (size_t)k);
        if (!T) return CHACAPK_HARITH_ERR_NULL;
        /* T = V^T B  (kt x k) */
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    kt, k, n, 1.0, leaf_rk_V(node), n, B, ldB, 0.0, T, kt);
        /* C += alpha U T  (m x k) */
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    m, k, kt, alpha, leaf_rk_U(node), m, T, kt, 1.0, C, ldC);
        free(T);
        return CHACAPK_HARITH_OK;
    }
    /* internal: recurse, slicing B by child col-offset, C by child row-offset. */
    int nr = node->nrsons, nc = node->ncsons;
    int rb = node->dof_row_start, cb = node->dof_col_start;
    for (int j_s = 0; j_s < nc; j_s++) {
        for (int i_s = 0; i_s < nr; i_s++) {
            const st_cHACApK_block_node_t *ch = node->sons[i_s + j_s * nr];
            int row_off = ch->dof_row_start - rb;
            int col_off = ch->dof_col_start - cb;
            int rc = node_matmat_local(ch, alpha, B + col_off, ldB, k,
                                       C + row_off, ldC);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
    }
    return CHACAPK_HARITH_OK;
}

/* C[dof_ncols x k] += alpha * node^T * B[dof_nrows x k]  (LOCAL indexing). */
static int node_matmat_trans_local(const st_cHACApK_block_node_t *node, double alpha,
                                   const double *B, int ldB, int k,
                                   double *C, int ldC)
{
    if (!node || k <= 0) return CHACAPK_HARITH_OK;
    int m = node->dof_nrows, n = node->dof_ncols;

    if (leaf_is_dense(node)) {
        /* C(n x k) += alpha * data^T(n x m) * B(m x k) */
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    n, k, m, alpha, leaf_dense_data(node), m,
                    B, ldB, 1.0, C, ldC);
        return CHACAPK_HARITH_OK;
    }
    if (leaf_is_rk(node)) {
        int kt = leaf_rk_rank(node);
        if (kt <= 0) return CHACAPK_HARITH_OK;
        double *T = (double*)malloc(sizeof(double) * (size_t)kt * (size_t)k);
        if (!T) return CHACAPK_HARITH_ERR_NULL;
        /* node^T = V U^T;  T = U^T B  (kt x k) */
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    kt, k, m, 1.0, leaf_rk_U(node), m, B, ldB, 0.0, T, kt);
        /* C += alpha V T  (n x k) */
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    n, k, kt, alpha, leaf_rk_V(node), n, T, kt, 1.0, C, ldC);
        free(T);
        return CHACAPK_HARITH_OK;
    }
    /* internal: child^T maps child row-space (B+row_off) to col-space (C+col_off). */
    int nr = node->nrsons, nc = node->ncsons;
    int rb = node->dof_row_start, cb = node->dof_col_start;
    for (int j_s = 0; j_s < nc; j_s++) {
        for (int i_s = 0; i_s < nr; i_s++) {
            const st_cHACApK_block_node_t *ch = node->sons[i_s + j_s * nr];
            int row_off = ch->dof_row_start - rb;
            int col_off = ch->dof_col_start - cb;
            int rc = node_matmat_trans_local(ch, alpha, B + row_off, ldB, k,
                                             C + col_off, ldC);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
    }
    return CHACAPK_HARITH_OK;
}

/* Accumulator append: grow an rk leaf's factors by the rank-kinc increment
 * (alpha*Uinc)(Vinc)^T WITHOUT recompressing. U V^T with appended columns is
 * the exact same matrix plus the increment, so the leaf stays valid to read
 * mid-factorization; recompression is deferred to rkleaf_flush_one. */
static int rkleaf_append(st_cHACApK_block_node_t *node, double alpha,
                         const double *Uinc, int ldU,
                         const double *Vinc, int ldV, int kinc)
{
    if (kinc <= 0) return CHACAPK_HARITH_OK;
    int m = node->dof_nrows, n = node->dof_ncols;
    int kc = node->leaf_mtx->kt;
    int kw = kc + kinc;
    double *a1n = (double*)malloc(sizeof(double) * (size_t)m * (size_t)kw);
    double *a2n = (double*)malloc(sizeof(double) * (size_t)n * (size_t)kw);
    if (!a1n || !a2n) { free(a1n); free(a2n); return CHACAPK_HARITH_ERR_NULL; }
    if (kc > 0) {
        memcpy(a1n, node->leaf_mtx->a1, sizeof(double) * (size_t)m * (size_t)kc);
        memcpy(a2n, node->leaf_mtx->a2, sizeof(double) * (size_t)n * (size_t)kc);
    }
    for (int j = 0; j < kinc; j++) {
        double *dU = a1n + (size_t)(kc + j) * (size_t)m;
        const double *sU = Uinc + (size_t)j * (size_t)ldU;
        for (int i = 0; i < m; i++) dU[i] = alpha * sU[i];
        memcpy(a2n + (size_t)(kc + j) * (size_t)n, Vinc + (size_t)j * (size_t)ldV,
               sizeof(double) * (size_t)n);
    }
    free(node->leaf_mtx->a1); free(node->leaf_mtx->a2);
    node->leaf_mtx->a1 = a1n; node->leaf_mtx->a2 = a2n; node->leaf_mtx->kt = kw;
    return CHACAPK_HARITH_OK;
}

/* Recompress an rk leaf to g_hlu_trunc_tol if forced, or if its accumulated
 * rank exceeds g_hlu_accum_cap. No-op otherwise (accumulator keeps growing).
 * g_hlu_accum_cap = 0 -> always recompress (accumulator disabled). */
static int rkleaf_flush_one(st_cHACApK_block_node_t *node, int force)
{
    if (!leaf_is_rk(node)) return CHACAPK_HARITH_OK;
    int kt = node->leaf_mtx->kt;
    if (kt <= 0) return CHACAPK_HARITH_OK;
    if (!force && g_hlu_accum_cap > 0 && kt <= g_hlu_accum_cap)
        return CHACAPK_HARITH_OK;
    int m = node->dof_nrows, n = node->dof_ncols;
    double *Un = NULL, *Vn = NULL; int kn = 0;
    int rc = rkleaf_recompress(node->leaf_mtx->a1, node->leaf_mtx->a2,
                               m, n, kt, g_hlu_trunc_tol, kt, &Un, &Vn, &kn);
    if (rc != CHACAPK_HARITH_OK) { free(Un); free(Vn); return rc; }
    free(node->leaf_mtx->a1); free(node->leaf_mtx->a2);
    node->leaf_mtx->a1 = Un; node->leaf_mtx->a2 = Vn; node->leaf_mtx->kt = kn;
    return CHACAPK_HARITH_OK;
}

/* Final pass: recompress every rk leaf (force) so the stored L/U factors are
 * compact for the solve. Called once after hlu_rec. */
static int flush_tree(st_cHACApK_block_node_t *node)
{
    if (!node) return CHACAPK_HARITH_OK;
    if (leaf_is_rk(node)) return rkleaf_flush_one(node, 1);
    if (node->sons) {
        int ns = node->nrsons * node->ncsons;
        for (int i = 0; i < ns; i++) {
            int rc = flush_tree(node->sons[i]);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
    }
    return CHACAPK_HARITH_OK;
}

/* C += alpha * U * W^T, where the increment is rank-k: U is (m x k, ld ldU),
 * W is (n x k, ld ldW). Mirrors add_dense_to_node but carries the low-rank
 * factors so an internal C never gets a full m x n densification.
 *   - dense leaf C: GEMM accumulate.
 *   - rk leaf C: append increment columns, lazy-recompress (accumulator).
 *   - internal C: recurse, slicing U by rows, W by rows (= product cols). */
static int add_lowrank_to_node(double alpha, const double *U, int ldU,
                               const double *W, int ldW, int m, int n, int k,
                               st_cHACApK_block_node_t *C)
{
    if (!C) return CHACAPK_HARITH_ERR_NULL;
    if (k <= 0) return CHACAPK_HARITH_OK;

    if (leaf_is_dense(C)) {
        if (leaf_rows(C) != m || leaf_cols(C) != n) return CHACAPK_HARITH_ERR_TOPOLOGY;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    m, n, k, alpha, U, ldU, W, ldW, 1.0, leaf_dense_data(C), m);
        return CHACAPK_HARITH_OK;
    }
    if (leaf_is_rk(C)) {
        if (leaf_rows(C) != m || leaf_cols(C) != n) return CHACAPK_HARITH_ERR_TOPOLOGY;
        /* Accumulator: append the increment columns (exact), recompress only
         * when the accumulated rank exceeds the cap (or cap=0 => always). */
        int rc = rkleaf_append(C, alpha, U, ldU, W, ldW, k);
        if (rc == CHACAPK_HARITH_OK) rc = rkleaf_flush_one(C, 0);
        return rc;
    }
    /* internal: recurse with row-slices of U and W. */
    int nr = C->nrsons, nc = C->ncsons;
    int rb = C->dof_row_start, cb = C->dof_col_start;
    for (int j_s = 0; j_s < nc; j_s++) {
        for (int i_s = 0; i_s < nr; i_s++) {
            st_cHACApK_block_node_t *ch = C->sons[i_s + j_s * nr];
            int cm = ch->dof_nrows, cn = ch->dof_ncols;
            int row_off = ch->dof_row_start - rb;
            int col_off = ch->dof_col_start - cb;
            int rc = add_lowrank_to_node(alpha, U + row_off, ldU,
                                         W + col_off, ldW, cm, cn, k, ch);
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

/* Context + body for parallelizing the all-internal h_addmul output blocks.
 * Each (i,k) output block C(i,k) is independent; the inner j-sum that
 * accumulates into one C(i,k) stays serial inside the body (no race on a
 * leaf), so distinct tasks touch distinct C leaves -> race-free. */
typedef struct {
    double alpha;
    const st_cHACApK_block_node_t *A;
    const st_cHACApK_block_node_t *B;
    st_cHACApK_block_node_t *C;
    int nrC, ncC, nmA;
    int *rc;   /* per-(i,k) return code, size nrC*ncC */
} addmul_ik_ctx_t;

static void addmul_ik_body(int idx, void *p)
{
    addmul_ik_ctx_t *c = (addmul_ik_ctx_t*)p;
    int i = idx % c->nrC;
    int k = idx / c->nrC;
    int rc = CHACAPK_HARITH_OK;
    for (int j = 0; j < c->nmA && rc == CHACAPK_HARITH_OK; j++) {
        rc = h_addmul(c->alpha,
                      c->A->sons[i + j * c->nrC],
                      c->B->sons[j + k * c->nmA],
                      c->C->sons[i + k * c->nrC]);
    }
    c->rc[idx] = rc;
}

/* ---------- materialize-FREE R(A*B) into a low-rank target ----------------- *
 *
 * The dominant near-cubic H-LU cost is the addmul case "A internal, B internal,
 * C an admissible (rk) leaf": the Phase-3.6 fallback densified the WHOLE internal
 * B subtree (materialize_node_as_dense), which is O(block^2) scratch and drives
 * the factor toward O(N^3) on deep trees.  Since C is admissible, A*B restricted
 * to C IS low-rank, so we approximate R(A*B) directly by randomized range-finding
 * (Halko-Martinsson-Tropp 2011, Alg 4.3 + 1 power iteration).  The action of A*B
 * and (A*B)^T on a dense sketch is applied THROUGH node_matmat_local/_trans_local
 * (compression-aware, never densifies A or B), so the cost is
 *   O((m+n+inner) * rr * rank(A,B))  -- sub-cubic.
 * The captured range is recompressed to g_hlu_trunc_tol and added to C via the
 * accumulator.  If the block is NOT low-rank within a rank ceiling (rare for an
 * admissible C), returns ADDMUL_RANK_SAT so the caller uses the EXACT materialize
 * path -- a correctness-preserving fallback (same result, slower), never a silent
 * wrong answer.  No static state (seed derived from C's DOF offsets) -> thread-safe
 * under the parallel h_addmul AND deterministic across runs/resumes. */
#define ADDMUL_RANK_SAT 1

static void rand_fill(double *buf, size_t cnt, unsigned long long seed)
{
    unsigned long long s = seed ? seed : 0x9E3779B97F4A7C15ULL;
    for (size_t i = 0; i < cnt; i++) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        buf[i] = (double)((s >> 33) & 0x7fffffff) / 2147483647.0 - 0.5;
    }
}

static int addmul_ii_to_rk(double alpha,
        const st_cHACApK_block_node_t *A, const st_cHACApK_block_node_t *B,
        st_cHACApK_block_node_t *C)
{
    const int m = C->dof_nrows, n = C->dof_ncols, inner = A->dof_ncols;
    const int mn = (m < n) ? m : n;
    if (mn <= 0 || inner <= 0) return CHACAPK_HARITH_OK;

    const int p = 8;                 /* oversampling */
    /* Rank ceiling: above ~mn/4 a low-rank sketch is no longer cheaper than the
     * exact dense path, so bail to materialize there (correct, and not worse for
     * that block).  But keep it generous (>= 256, up to a 1024 hard cap) so the
     * moderately-ranked admissible blocks that appear at larger N stay on the
     * materialize-FREE path instead of needlessly densifying (the residual
     * materialize seen at the old 256 ceiling). */
    int r_ceiling = mn / 4;
    if (r_ceiling < 256)  r_ceiling = 256;
    if (r_ceiling > 1024) r_ceiling = 1024;
    if (r_ceiling > mn)   r_ceiling = mn;
    /* Seed the target rank from C's current rank (the Schur update into an
     * admissible block has comparable rank), so well-ranked blocks resolve in
     * ONE sketch instead of doubling from a fixed 64. */
    int kc0 = leaf_rk_rank(C);
    int r = 2 * (kc0 > 0 ? kc0 : (g_hlu_accum_cap > 0 ? g_hlu_accum_cap : 16));
    if (r < 32) r = 32;
    if (r > mn) r = mn;

    for (;;) {
        int rr = r + p; if (rr > mn) rr = mn;
        const size_t snr = (size_t)inner * rr, smr = (size_t)m * rr, snn = (size_t)n * rr;

        double *Om = (double*)malloc(sizeof(double) * snn);
        double *S  = (double*)calloc(snr, sizeof(double));
        double *Y  = (double*)calloc(smr, sizeof(double));
        double *W  = (double*)calloc(snr, sizeof(double));
        double *S2 = (double*)calloc(snn, sizeof(double));
        double *tau= (double*)malloc(sizeof(double) * (size_t)rr);
        if (!Om||!S||!Y||!W||!S2||!tau) {
            free(Om);free(S);free(Y);free(W);free(S2);free(tau);
            return CHACAPK_HARITH_ERR_NULL;
        }
        unsigned long long seed =
            0x243F6A8885A308D3ULL
            ^ ((unsigned long long)C->dof_row_start * 2654435761ULL)
            ^ ((unsigned long long)C->dof_col_start * 40503ULL)
            ^ ((unsigned long long)rr * 2246822519ULL);
        rand_fill(Om, snn, seed);

        /* Y = (A*B) Om = A (B Om) */
        int rc = node_matmat_local(B, 1.0, Om, n, rr, S, inner);
        if (rc==CHACAPK_HARITH_OK) rc = node_matmat_local(A, 1.0, S, inner, rr, Y, m);
        /* one power iteration: Y <- (A*B)(A*B)^T Y */
        if (rc==CHACAPK_HARITH_OK) {
            memset(W, 0, sizeof(double)*snr); memset(S2, 0, sizeof(double)*snn);
            rc = node_matmat_trans_local(A, 1.0, Y, m, rr, W, inner);          /* W = A^T Y (inner x rr) */
            if (rc==CHACAPK_HARITH_OK) rc = node_matmat_trans_local(B, 1.0, W, inner, rr, S2, n); /* S2 = B^T W (n x rr) */
            if (rc==CHACAPK_HARITH_OK) {
                memset(S, 0, sizeof(double)*snr); memset(Y, 0, sizeof(double)*smr);
                rc = node_matmat_local(B, 1.0, S2, n, rr, S, inner);           /* S = B S2 */
                if (rc==CHACAPK_HARITH_OK) rc = node_matmat_local(A, 1.0, S, inner, rr, Y, m); /* Y = A S */
            }
        }
        if (rc != CHACAPK_HARITH_OK) { free(Om);free(S);free(Y);free(W);free(S2);free(tau); return rc; }

        /* QR(Y) -> orthonormal Q (m x rr) overwriting Y */
        int info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m, rr, Y, m, tau);
        if (info==0) info = LAPACKE_dorgqr(LAPACK_COL_MAJOR, m, rr, rr, Y, m, tau);
        if (info != 0) { free(Om);free(S);free(Y);free(W);free(S2);free(tau); return CHACAPK_HARITH_ERR_LAPACK; }

        /* V = (A*B)^T Q = B^T (A^T Q)  (n x rr), reusing W (inner x rr) and S2 (n x rr) */
        memset(W, 0, sizeof(double)*snr); memset(S2, 0, sizeof(double)*snn);
        rc = node_matmat_trans_local(A, 1.0, Y, m, rr, W, inner);
        if (rc==CHACAPK_HARITH_OK) rc = node_matmat_trans_local(B, 1.0, W, inner, rr, S2, n);
        if (rc != CHACAPK_HARITH_OK) { free(Om);free(S);free(Y);free(W);free(S2);free(tau); return rc; }

        /* A*B ~= Y * S2^T ; recompress to trunc tol */
        double *Un=NULL, *Vn=NULL; int kn=0;
        rc = rkleaf_recompress(Y, S2, m, n, rr, g_hlu_trunc_tol, rr, &Un, &Vn, &kn);
        free(Om);free(S);free(Y);free(W);free(S2);free(tau);
        if (rc != CHACAPK_HARITH_OK) { free(Un); free(Vn); return rc; }

        if (kn >= rr && rr < mn) {
            /* range not yet resolved at this rank */
            free(Un); free(Vn);
            if (rr >= r_ceiling) return ADDMUL_RANK_SAT;  /* not low-rank -> exact path */
            r = (rr * 2 > mn) ? mn : rr * 2;              /* grow + retry */
            continue;
        }
        rc = add_lowrank_to_node(alpha, Un, m, Vn, n, m, n, kn, C);
        free(Un); free(Vn);
        return rc;
    }
}


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
                                    g_hlu_trunc_tol, kw, &U_new, &V_new, &k_new);
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
                                    g_hlu_trunc_tol, kw, &U_new, &V_new, &k_new);
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
                                    g_hlu_trunc_tol, kw, &U_new, &V_new, &k_new);
        free(U_widened); free(V_widened);
        if (rc != CHACAPK_HARITH_OK) { free(U_new); free(V_new); return rc; }
        free(C->leaf_mtx->a1); free(C->leaf_mtx->a2);
        C->leaf_mtx->a1 = U_new; C->leaf_mtx->a2 = V_new;
        C->leaf_mtx->kt = k_new;
        g_stats.n_dense_gemm++;
        return CHACAPK_HARITH_OK;
    }

    /* ---- materialize-FREE mixed multiply when an operand is an rk leaf ----
     * These cover the dominant mixed cases at small leaf (rk x rk, internal x
     * rk, rk x internal -- ~73-93% of mixed addmul). Factor through the rk
     * waist so the OTHER (possibly internal) operand is only touched via
     * matvec, never densified:
     *   rk A = U_A V_A^T:  C += alpha U_A (V_A^T B) = add_lowrank(U_A, B^T V_A)
     *   rk B = U_B V_B^T:  C += alpha (A U_B) V_B^T = add_lowrank(A U_B, V_B) */
    if (leaf_is_rk(A)) {
        int m = C->dof_nrows, n = C->dof_ncols, inner = A->dof_ncols;
        int k = leaf_rk_rank(A);
        if (k <= 0) return CHACAPK_HARITH_OK;
        /* W (n x k) = B^T V_A  (V_A is inner x k = B's row space). */
        double *W = (double*)calloc((size_t)n * (size_t)k, sizeof(double));
        if (!W) return CHACAPK_HARITH_ERR_NULL;
        int rc = node_matmat_trans_local(B, 1.0, leaf_rk_V(A), inner, k, W, n);
        if (rc == CHACAPK_HARITH_OK)
            rc = add_lowrank_to_node(alpha, leaf_rk_U(A), m, W, n, m, n, k, C);
        free(W);
        g_stats.n_dense_gemm++;
        return rc;
    }
    if (leaf_is_rk(B)) {
        int m = C->dof_nrows, n = C->dof_ncols, inner = A->dof_ncols;
        int k = leaf_rk_rank(B);
        if (k <= 0) return CHACAPK_HARITH_OK;
        /* P (m x k) = A U_B  (U_B is inner x k = A's col space). */
        double *P = (double*)calloc((size_t)m * (size_t)k, sizeof(double));
        if (!P) return CHACAPK_HARITH_ERR_NULL;
        int rc = node_matmat_local(A, 1.0, leaf_rk_U(B), inner, k, P, m);
        if (rc == CHACAPK_HARITH_OK)
            rc = add_lowrank_to_node(alpha, P, m, leaf_rk_V(B), n, m, n, k, C);
        free(P);
        g_stats.n_dense_gemm++;
        return rc;
    }

    /* dense(A) * dense(B) -> rk(C) + remaining MIXED leaf+internal cases
     * (internal x internal, internal x dense, dense x internal, dense*dense->rk):
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
            /* Materialize-FREE fast path for the dominant near-cubic trigger:
             * A internal, B internal, C an admissible (rk) leaf.  Randomized
             * R(A*B) instead of densifying the whole internal B subtree.  On
             * ADDMUL_RANK_SAT (block not low-rank within ceiling -- rare) fall
             * through to the exact materialize path below (same result). */
            if (!leaf_is_dense(A) && !leaf_is_rk(A) &&
                !leaf_is_dense(B) && !leaf_is_rk(B) && leaf_is_rk(C)) {
                int rcq = addmul_ii_to_rk(alpha, A, B, C);
                if (rcq != ADDMUL_RANK_SAT) {
                    if (rcq == CHACAPK_HARITH_OK) g_stats.n_dense_gemm++;
                    return rcq;
                }
            }
            g_dbg_mixed_addmul[node_kind_idx(A)*3 + node_kind_idx(B)]++;
            int m = C->dof_nrows;
            int n = C->dof_ncols;
            int inner = A->dof_ncols;  /* == B->dof_nrows */
            int rc;

            if (!leaf_is_dense(A) && !leaf_is_rk(A)) {
                /* A internal (covers internal x internal, internal x dense):
                 * matmat THROUGH A so A's subtree compression is exploited
                 * (rk leaves handled as U(V^T...) ); materialize only B. */
                double *B_dense = materialize_node_as_dense(B);   /* inner x n */
                double *D = (double*)calloc((size_t)m * (size_t)n, sizeof(double));
                if (!B_dense || !D) { free(B_dense); free(D);
                                       return CHACAPK_HARITH_ERR_NULL; }
                rc = node_matmat_local(A, alpha, B_dense, inner, n, D, m);
                free(B_dense);
                if (rc == CHACAPK_HARITH_OK) rc = add_dense_to_node(D, m, m, n, C);
                free(D);
            } else if (!leaf_is_dense(B) && !leaf_is_rk(B)) {
                /* A dense leaf, B INTERNAL (densexinternal -- the MEASURED
                 * dominant materialize trigger, 2026-06-04). Compute
                 * D = A*B = (B^T A^T)^T via node_matmat_trans_local THROUGH B,
                 * so B's internal subtree is NEVER densified.  The old path
                 * materialize_node_as_dense(B) on an internal B was the bulk of
                 * the 2.0e9-elem / 1.21M-call materialize wall that made the
                 * factor scale ~N^5.  A is a dense leaf -> A^T is cheap; tmp =
                 * B^T A^T is (n x m) with m leaf-bounded (= A's row cluster). */
                const double *Ad = leaf_dense_data(A);            /* m x inner, ld m */
                double *At  = (double*)malloc(sizeof(double) * (size_t)inner * (size_t)m);
                double *tmp = (double*)calloc((size_t)n * (size_t)m, sizeof(double));
                double *D   = (double*)malloc(sizeof(double) * (size_t)m * (size_t)n);
                if (!Ad || !At || !tmp || !D) {
                    free(At); free(tmp); free(D);
                    return CHACAPK_HARITH_ERR_NULL;
                }
                /* At = A^T  (inner x m) */
                for (int jj = 0; jj < inner; jj++)
                    for (int ii = 0; ii < m; ii++)
                        At[jj + (size_t)ii * (size_t)inner] = Ad[ii + (size_t)jj * (size_t)m];
                /* tmp(n x m) = alpha * B^T * At  (B stays compressed) */
                rc = node_matmat_trans_local(B, alpha, At, inner, m, tmp, n);
                free(At);
                if (rc == CHACAPK_HARITH_OK) {
                    /* D(m x n) = tmp^T, then C += D */
                    for (int bcol = 0; bcol < n; bcol++)
                        for (int arow = 0; arow < m; arow++)
                            D[arow + (size_t)bcol * (size_t)m] =
                                tmp[bcol + (size_t)arow * (size_t)n];
                    rc = add_dense_to_node(D, m, m, n, C);
                }
                free(tmp); free(D);
            } else {
                /* A is a dense leaf, B a dense leaf (dense*dense->rk): both
                 * small, plain dgemm. */
                double *A_dense = materialize_node_as_dense(A);
                double *B_dense = materialize_node_as_dense(B);
                double *D = (double*)malloc(sizeof(double) * (size_t)m * (size_t)n);
                if (!A_dense || !B_dense || !D) {
                    free(A_dense); free(B_dense); free(D);
                    return CHACAPK_HARITH_ERR_NULL;
                }
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            m, n, inner, alpha, A_dense, m, B_dense, inner, 0.0, D, m);
                free(A_dense); free(B_dense);
                rc = add_dense_to_node(D, m, m, n, C);
                free(D);
            }
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

    /* Block-parallel over the independent output blocks (i,k) when C is large
     * enough to amortize task overhead. Distinct (i,k) -> distinct C leaves;
     * the inner j-sum stays serial inside addmul_ik_body -> race-free. ngcore
     * work-stealing handles the nested recursion. */
    long blk_area = (long)C->dof_nrows * (long)C->dof_ncols;
    int n_ik = nrC * ncC;
    if (g_hlu_parallel && n_ik > 1 && blk_area > g_hlu_par_cutoff) {
        int *rc = (int*)malloc(sizeof(int) * (size_t)n_ik);
        if (!rc) return CHACAPK_HARITH_ERR_NULL;
        for (int t = 0; t < n_ik; t++) rc[t] = CHACAPK_HARITH_OK;
        addmul_ik_ctx_t ctx = { alpha, A, B, C, nrC, ncC, nmA, rc };
        chacapk_par_for(n_ik, addmul_ik_body, &ctx);
        int rc_first = CHACAPK_HARITH_OK;
        for (int t = 0; t < n_ik; t++)
            if (rc[t] != CHACAPK_HARITH_OK) { rc_first = rc[t]; break; }
        free(rc);
        return rc_first;
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


/* Materialize-free forward solve  L * Z = Z  (overwrite), L unit-lower, applied
 * to a DENSE rhs block Z (L->dof_nrows x k, leading dim ldZ).  Recurses through
 * an internal L using node_matmat_local for the off-diagonal updates so L's
 * subtree is NEVER densified -- replaces the materialize_node_as_dense(L)
 * triangular-solve fallback that was ~half the H-LU materialize wall
 * (2026-06-04: lln = 966 of 1982 Me at C-type 6^3). */
static int htrsm_lln_dense_rhs(const st_cHACApK_block_node_t *L,
                               double *Z, int ldZ, int k)
{
    if (!L) return CHACAPK_HARITH_ERR_NULL;
    if (k <= 0) return CHACAPK_HARITH_OK;
    if (leaf_is_rk(L)) return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
    if (leaf_is_dense(L)) {
        int m = leaf_rows(L);
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                    m, k, 1.0, leaf_dense_data(L), m, Z, ldZ);
        return CHACAPK_HARITH_OK;
    }
    /* internal L: block forward substitution. row partition == col partition. */
    int s = L->nrsons;
    if (L->ncsons != s) return CHACAPK_HARITH_ERR_TOPOLOGY;
    int rb = L->dof_row_start, cb = L->dof_col_start;
    for (int i = 0; i < s; i++) {
        const st_cHACApK_block_node_t *Lii = L->sons[i + i * s];
        int ri = Lii->dof_row_start - rb;
        for (int kk = 0; kk < i; kk++) {
            const st_cHACApK_block_node_t *Lik = L->sons[i + kk * s];
            int ck = Lik->dof_col_start - cb;   /* local row of already-solved Z block kk */
            int rr = Lik->dof_row_start - rb;   /* local row of Z block i (== ri) */
            int rc = node_matmat_local(Lik, -1.0, Z + ck, ldZ, k, Z + rr, ldZ);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
        int rc = htrsm_lln_dense_rhs(Lii, Z + ri, ldZ, k);
        if (rc != CHACAPK_HARITH_OK) return rc;
    }
    return CHACAPK_HARITH_OK;
}

/* Materialize-free solve  U^T * V = V  (overwrite), U non-unit UPPER (so U^T is
 * lower), applied to a DENSE rhs block V (U->dof_nrows x k, ld ldV).  Recurses
 * through an internal U via node_matmat_trans_local.  This is the backward
 * (run) triangular solve building block: the rk case (Y U = U_x V_x^T) reduces
 * to U^T V_y = V_x, and the dense case (X U = X) reduces to U^T X^T = X^T. */
static int htrsm_lUTn_dense_rhs(const st_cHACApK_block_node_t *U,
                                double *V, int ldV, int k)
{
    if (!U) return CHACAPK_HARITH_ERR_NULL;
    if (k <= 0) return CHACAPK_HARITH_OK;
    if (leaf_is_rk(U)) return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
    if (leaf_is_dense(U)) {
        int n = leaf_rows(U);   /* square */
        cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit,
                    n, k, 1.0, leaf_dense_data(U), n, V, ldV);
        return CHACAPK_HARITH_OK;
    }
    /* internal U: forward substitution on U^T (lower). */
    int s = U->nrsons;
    if (U->ncsons != s) return CHACAPK_HARITH_ERR_TOPOLOGY;
    int rb = U->dof_row_start;
    for (int i = 0; i < s; i++) {
        const st_cHACApK_block_node_t *Uii = U->sons[i + i * s];
        int ri = Uii->dof_row_start - rb;
        for (int kk = 0; kk < i; kk++) {
            const st_cHACApK_block_node_t *Uki = U->sons[kk + i * s];  /* U[kk,i] */
            const st_cHACApK_block_node_t *Ukk = U->sons[kk + kk * s];
            int rkk = Ukk->dof_row_start - rb;   /* local row of V block kk */
            /* V_i -= U[kk,i]^T * V_kk */
            int rc = node_matmat_trans_local(Uki, -1.0, V + rkk, ldV, k, V + ri, ldV);
            if (rc != CHACAPK_HARITH_OK) return rc;
        }
        int rc = htrsm_lUTn_dense_rhs(Uii, V + ri, ldV, k);
        if (rc != CHACAPK_HARITH_OK) return rc;
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

    /* L INTERNAL x X leaf: solve materialize-FREE by recursing the forward
     * substitution through L (was: materialize the internal L, ~half the H-LU
     * materialize wall).  (L leaf x X internal cannot occur: a dense-leaf L has
     * leaf-level rows, so X -- same row partition -- is also a leaf.) */
    {
        int L_is_leaf = leaf_is_dense(L) || leaf_is_rk(L);  /* rk-L unreachable above */
        int X_is_leaf = leaf_is_dense(X) || leaf_is_rk(X);
        int mixed = (!L_is_leaf || !X_is_leaf) && (L_is_leaf || X_is_leaf);
        if (mixed) {
            g_dbg_mixed_lln[node_kind_idx(L)*3 + node_kind_idx(X)]++;
            int rc;
            if (!L_is_leaf) {
                if (leaf_is_rk(X))
                    /* L Y = U_x V_x^T -> Y = (L^-1 U_x) V_x^T; solve L U_x = U_x. */
                    rc = htrsm_lln_dense_rhs(L, leaf_rk_U(X), leaf_rows(X), leaf_rk_rank(X));
                else
                    rc = htrsm_lln_dense_rhs(L, leaf_dense_data(X), leaf_rows(X), leaf_cols(X));
            } else {
                /* L dense leaf, X internal.  X shares L's (terminal) ROW cluster,
                 * so X is subdivided by COLUMNS only (nrsons==1); recurse the solve
                 * over X's column sons -- materialize-FREE (was: densify L AND X, a
                 * dominant residual materialize at larger N).  L^-1 applies to each
                 * column son independently. */
                if (X->nrsons != 1) return CHACAPK_HARITH_ERR_TOPOLOGY;
                int sc = X->ncsons;
                rc = CHACAPK_HARITH_OK;
                for (int j = 0; j < sc && rc == CHACAPK_HARITH_OK; j++)
                    rc = htrsm_lln(L, X->sons[j]);   /* sons[0 + j*1] */
            }
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

    /* U INTERNAL x X leaf: solve materialize-FREE via the U^T forward solver
     * (was: materialize the internal U -- the other ~half of the materialize
     * wall).  rk X: Y U = U_x V_x^T -> U_y=U_x, solve U^T V_y = V_x.  dense X:
     * X U = X <=> U^T X^T = X^T (transpose, solve, transpose back). */
    {
        int U_is_leaf = leaf_is_dense(U) || leaf_is_rk(U);  /* rk-U unreachable above */
        int X_is_leaf = leaf_is_dense(X) || leaf_is_rk(X);
        int mixed = (!U_is_leaf || !X_is_leaf) && (U_is_leaf || X_is_leaf);
        if (mixed) {
            g_dbg_mixed_run[node_kind_idx(U)*3 + node_kind_idx(X)]++;
            int rc;
            if (!U_is_leaf) {
                if (leaf_is_rk(X)) {
                    rc = htrsm_lUTn_dense_rhs(U, leaf_rk_V(X), leaf_cols(X), leaf_rk_rank(X));
                } else {
                    int m = leaf_rows(X), n = leaf_cols(X);
                    double *Xt = (double*)malloc(sizeof(double) * (size_t)n * (size_t)m);
                    if (!Xt) return CHACAPK_HARITH_ERR_NULL;
                    double *Xd = leaf_dense_data(X);
                    for (int b = 0; b < n; b++)
                        for (int a = 0; a < m; a++)
                            Xt[b + (size_t)a * (size_t)n] = Xd[a + (size_t)b * (size_t)m];
                    rc = htrsm_lUTn_dense_rhs(U, Xt, n, m);
                    if (rc == CHACAPK_HARITH_OK)
                        for (int b = 0; b < n; b++)
                            for (int a = 0; a < m; a++)
                                Xd[a + (size_t)b * (size_t)m] = Xt[b + (size_t)a * (size_t)n];
                    free(Xt);
                }
            } else {
                /* U dense leaf, X internal.  X shares U's (terminal) COLUMN cluster,
                 * so X is subdivided by ROWS only (ncsons==1); recurse the solve over
                 * X's row sons -- materialize-FREE (was: densify U AND X).  Each row
                 * son X_i solves X_i U^-1 independently. */
                if (X->ncsons != 1) return CHACAPK_HARITH_ERR_TOPOLOGY;
                int sr = X->nrsons;
                rc = CHACAPK_HARITH_OK;
                for (int i = 0; i < sr && rc == CHACAPK_HARITH_OK; i++)
                    rc = htrsm_run(U, X->sons[i]);   /* sons[i + 0*nrsons] */
            }
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
        int r0 = node->dof_row_start;  /* 0-based DOF start */
        int c0 = node->dof_col_start;
        cblas_dgemv(CblasColMajor, CblasNoTrans, m, n,
                    -1.0, leaf_dense_data(node), m,
                    &x[c0], 1, 1.0, &y[r0], 1);
        return CHACAPK_HARITH_OK;
    }

    /* rk leaf: y -= U_a (V_a^T x).  Two dgemvs through the kt-rank waist. */
    if (leaf_is_rk(node)) {
        int m = leaf_rows(node), n = leaf_cols(node);
        int kt = leaf_rk_rank(node);
        int r0 = node->dof_row_start;
        int c0 = node->dof_col_start;
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

/* Context to run hlu_rec inside a TaskManager region (C-callable body). */
typedef struct { st_cHACApK_block_node_t *root; int rc; } hlu_decomp_ctx_t;
static void hlu_decomp_body(void *p)
{
    hlu_decomp_ctx_t *c = (hlu_decomp_ctx_t*)p;
    c->rc = hlu_rec(c->root);
}

int cHACApK_hlu_decomp(st_cHACApK_block_node_t *root)
{
    stats_reset();
    clear_ipiv_registry();
    g_dbg_n_materialize = 0;
    g_dbg_materialize_elems = 0;
    g_dbg_mat_internal = 0;
    g_dbg_mat_leaf = 0;
    for (int i = 0; i < 9; i++) { g_dbg_mixed_addmul[i] = 0; g_dbg_mixed_lln[i] = 0; g_dbg_mixed_run[i] = 0; }
    clock_t t0 = clock();
    /* Wrap in a TaskManager region so the block-parallel h_addmul (i,k) loop
     * actually uses the threadpool. Serial path if parallelism is disabled. */
    hlu_decomp_ctx_t ctx = { root, CHACAPK_HARITH_OK };
    if (g_hlu_parallel) chacapk_par_region(hlu_decomp_body, &ctx);
    else                hlu_decomp_body(&ctx);
    /* Accumulator final flush: compress all rk-leaf L/U factors that grew via
     * deferred recompression, so the solve uses compact factors. */
    if (ctx.rc == CHACAPK_HARITH_OK && g_hlu_accum_cap > 0)
        ctx.rc = flush_tree(root);
    g_stats.t_decomp_sec = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
    return ctx.rc;
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
        int row0 = node->dof_row_start;  /* 0-based DOF start (Phase 3.6 dof_* fields) */
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
        int row0 = node->dof_row_start;  /* 0-based DOF start (Phase 3.6 dof_* fields) */
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
    clock_t t0 = clock();
    int rc = hlu_forward_rec(root, x);
    if (rc == CHACAPK_HARITH_OK) rc = hlu_backward_rec(root, x);
    g_stats.t_solve_sec = (double)(clock() - t0) / (double)CLOCKS_PER_SEC;
    return rc;
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
    /* Synthetic test uses cluster nstrt/nsize as DOF (nffc=1). */
    bn->dof_nrows     = row_clt->nsize;
    bn->dof_ncols     = col_clt->nsize;
    bn->dof_row_start = row_clt->nstrt - 1;
    bn->dof_col_start = col_clt->nstrt - 1;

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

/* Build a "diagonal-refined" block tree: the H-matrix shape where the
 * near-field DIAGONAL is recursively refined into sub-trees while the
 * OFF-DIAGONAL far-field blocks are DENSE LEAVES.  This is the realistic
 * admissibility-driven H-matrix shape AND the documented validated subset
 * for H-LDL^T (off-diagonal blocks are dense; diagonal blocks may be deep).
 *
 *   is_diag = 1 : this node is on the block diagonal (row span == col span).
 *                 depth>0 -> split 2x2; the two diagonal children recurse
 *                 (is_diag=1), the two off-diagonal children are dense leaves.
 *   is_diag = 0 : off-diagonal node -> always a dense leaf.
 * depth==0 -> dense leaf regardless. */
static st_cHACApK_block_node_t *build_diag_refined_tree(
    deep_state_t *s,
    st_cHACApK_cluster_t *row_clt, st_cHACApK_cluster_t *col_clt,
    int depth, int is_diag)
{
    st_cHACApK_block_node_t *bn =
        (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn));
    log_bn(s, bn);
    bn->row_cluster = row_clt;
    bn->col_cluster = col_clt;
    bn->dof_nrows     = row_clt->nsize;
    bn->dof_ncols     = col_clt->nsize;
    bn->dof_row_start = row_clt->nstrt - 1;
    bn->dof_col_start = col_clt->nstrt - 1;

    if (depth == 0 || !is_diag) {
        int m = row_clt->nsize, n = col_clt->nsize;
        st_cHACApK_leafmtx_t *L = (st_cHACApK_leafmtx_t*)calloc(1, sizeof(*L));
        log_lf(s, L);
        L->ltmtx = 2; L->ndl = m; L->ndt = n;
        L->nstrtl = row_clt->nstrt; L->nstrtt = col_clt->nstrt;
        L->a1 = (double*)malloc(sizeof(double)*(size_t)m*(size_t)n);
        int r0 = row_clt->nstrt - 1, c0 = col_clt->nstrt - 1, N = s->N_global;
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                L->a1[i + (size_t)j*(size_t)m] = s->A_full[(r0+i) + (size_t)(c0+j)*(size_t)N];
        bn->leaf_mtx  = L;
        bn->leaf_kind = 2;
        return bn;
    }

    /* diagonal internal node: 2x2 split */
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
        for (int i = 0; i < 2; i++) {
            int child_is_diag = (i == j);   /* only (0,0) and (1,1) are diagonal */
            bn->sons[i + j*2] = build_diag_refined_tree(
                s, rc[i], cc[j], depth - 1, child_is_diag);
        }
    return bn;
}

static void deep_cleanup(deep_state_t *s)
{
    for (int i = 0; i < s->n_bn; i++) {
        if (s->bn_log[i]->sons) free(s->bn_log[i]->sons);
        free(s->bn_log[i]);
    }
    /* Free a1 always; a2 only for rk leaves (it is NULL for dense leaves, so
     * free(NULL) is harmless -- but a possibly-recompressed rk a2 must be
     * released to avoid a leak when a tree contains rk off-diagonal leaves). */
    for (int i = 0; i < s->n_lf; i++) { free(s->lf_log[i]->a1); free(s->lf_log[i]->a2); free(s->lf_log[i]); }
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
    bn_A.dof_nrows = m; bn_A.dof_ncols = inner;
    bn_B.row_cluster = clt_in; bn_B.col_cluster = clt_n;  bn_B.leaf_mtx = &lf_B; bn_B.leaf_kind = 1;
    bn_B.dof_nrows = inner; bn_B.dof_ncols = n;
    bn_C.row_cluster = clt_m;  bn_C.col_cluster = clt_n;  bn_C.leaf_mtx = &lf_C; bn_C.leaf_kind = 1;
    bn_C.dof_nrows = m; bn_C.dof_ncols = n;

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
    root->dof_nrows = N; root->dof_ncols = N;
    root->dof_row_start = 0; root->dof_col_start = 0;

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
            bn->dof_nrows = nb; bn->dof_ncols = nb;
            bn->dof_row_start = rclt[i_son]->nstrt - 1;
            bn->dof_col_start = cclt[j_son]->nstrt - 1;
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
int cHACApK_hlu_debug_materialize(void *leafmtxp_void, void *control_void,
                                    int nffc,
                                    double *A_perm_out,
                                    int *lod_out,
                                    int *nd_out)
{
    if (!leafmtxp_void || !control_void || !A_perm_out || !lod_out || !nd_out)
        return -1;
    if (nffc <= 0) return -1;
    st_cHACApK_leafmtxp_t *lp = (st_cHACApK_leafmtxp_t*)leafmtxp_void;
    st_cHACApK_lcontrol_t *lc = (st_cHACApK_lcontrol_t*)control_void;
    int nd = lp->nd;
    if (nd <= 0 || !lp->st_clt_root || !lc->lod) return -1;

    *nd_out = nd;
    /* Copy lod (1-based -> 0-based). */
    for (int i = 0; i < nd; i++) lod_out[i] = lc->lod[i + 1] - 1;

    /* Convert leaves to internal layout. */
    cHACApK_convert_leafmtxp_to_internal(lp);

    /* Build block-tree. */
    st_cHACApK_block_node_t *root = cHACApK_build_block_tree_nffc(
        lp, lp->st_clt_root, lp->st_clt_root, nffc);
    if (!root) return -3;

    /* Materialize root as dense (in permuted ordering, col-major). */
    double *dense = materialize_node_as_dense(root);
    if (!dense) { cHACApK_free_block_tree(root); return -2; }
    memcpy(A_perm_out, dense, sizeof(double) * (size_t)nd * (size_t)nd);
    free(dense);
    cHACApK_free_block_tree(root);
    return 0;
}

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


/* ---------- H-LU as a reusable PRECONDITIONER (factor once, apply many) ----
 * cHACApK_hlu_run_on_hacapk does convert+build+decomp+solve+free in ONE shot
 * (self-test). For a preconditioner we keep the factored tree and apply the
 * solve repeatedly (per GMRES iteration). These three wrappers expose that:
 *   factor : convert leaves + build block-tree + hlu_decomp -> opaque root
 *   apply  : permute(lod) + hlu_solve_vec + un-permute   (r,z in ORIGINAL order)
 *   free   : free the block-tree
 * Used by the HACApK-based H-LU preconditioner path (factor once, apply per GMRES iteration). */
void* cHACApK_hlu_factor_leafmtxp(void* leafmtxp_void, void* control_void, int nffc)
{
    (void)control_void;
    if (!leafmtxp_void || nffc <= 0) return NULL;
    st_cHACApK_leafmtxp_t *lp = (st_cHACApK_leafmtxp_t*)leafmtxp_void;
    if (lp->nd <= 0 || !lp->st_clt_root) return NULL;
    cHACApK_convert_leafmtxp_to_internal(lp);
    st_cHACApK_block_node_t *root = cHACApK_build_block_tree_nffc(
        lp, lp->st_clt_root, lp->st_clt_root, nffc);
    if (!root) return NULL;
    int rc = cHACApK_hlu_decomp(root);
    if (rc != CHACAPK_HARITH_OK) { cHACApK_free_block_tree(root); return NULL; }
    return (void*)root;
}

int cHACApK_hlu_apply(void* root_void, void* control_void,
                      const double* r, double* z, int nd)
{
    if (!root_void || !control_void || !r || !z || nd <= 0) return -1;
    st_cHACApK_block_node_t *root = (st_cHACApK_block_node_t*)root_void;
    st_cHACApK_lcontrol_t   *lc   = (st_cHACApK_lcontrol_t*)control_void;
    if (!lc->lod) return -1;
    double *rp = (double*)malloc(sizeof(double) * (size_t)nd);
    double *zp = (double*)malloc(sizeof(double) * (size_t)nd);
    if (!rp || !zp) { free(rp); free(zp); return -2; }
    for (int i = 0; i < nd; i++) rp[i] = r[lc->lod[i + 1] - 1];        /* permute */
    int rc = cHACApK_hlu_solve_vec(root, rp, zp, nd);
    if (rc == CHACAPK_HARITH_OK)
        for (int i = 0; i < nd; i++) z[lc->lod[i + 1] - 1] = zp[i];    /* un-permute */
    free(rp); free(zp);
    return (rc == CHACAPK_HARITH_OK) ? 0 : -3;
}

void cHACApK_hlu_free_factors(void* root_void)
{
    if (root_void) cHACApK_free_block_tree((st_cHACApK_block_node_t*)root_void);
}




/* ---------- Phase 4 debug: mixed-sibling synthetic test ---------------- *
 *
 * Mimics Radia's nx=3 leaf=10 tree shape (10 dense leaves, 3 internal,
 * depth 3) using uniform leaf sizes so we can isolate the "mixed sibling"
 * recursion bug from any non-uniform-size effects.
 *
 *   N = 4 * nb_small.  Root (2x2 internal) has children:
 *     sons[0,0] = TL: internal 2x2 with 4 leaves of size nb_small
 *     sons[1,0] = BL: LEAF of size 2*nb_small (full row range below TL,
 *                     full col range of TL)
 *     sons[0,1] = TR: LEAF of size 2*nb_small (transposed)
 *     sons[1,1] = BR: internal 2x2 with 4 leaves of size nb_small
 *
 * Total: 4 + 1 + 1 + 4 = 10 leaves, 3 internal -- exact Radia match.
 *
 * If our recursive H-LU fails this test, the mixed-sibling recursion has
 * a bug. If it passes, the real-Radia bug must come from non-uniform
 * leaf sizes (HACApK's element-count splits like 13 -> 6 + 7). */
double cHACApK_harith_self_test_mixed_sibling(int nb_small)
{
    if (nb_small <= 0 || nb_small > 64) return -1.0;
    int N = 4 * nb_small;
    int half = 2 * nb_small;  /* size of root-level sub-block */

    unsigned long seed = 0xC0FFEEUL;
    #define RND_MS() (seed = seed * 6364136223846793005UL + 1442695040888963407UL, \
                      (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5)

    /* Allocate */
    double *A_full = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *A_ref  = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *b      = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_hlu  = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_ref  = (double*)malloc(sizeof(double) * (size_t)N);
    if (!A_full || !A_ref || !b || !x_hlu || !x_ref) return -2.0;

    /* Random diag-dominant A. */
    for (int j = 0; j < N; j++) {
        double rs = 0.0;
        for (int i = 0; i < N; i++) {
            double v = RND_MS();
            A_full[i + (size_t)j*(size_t)N] = v;
            if (i != j) rs += (v < 0.0 ? -v : v);
        }
        A_full[j + (size_t)j*(size_t)N] = rs + (double)N;
    }
    for (int j = 0; j < N; j++) b[j] = RND_MS() * 10.0;
    memcpy(A_ref, A_full, sizeof(double) * (size_t)N * (size_t)N);

    /* Reference dgesv. */
    int *ipiv_ref = (int*)malloc(sizeof(int) * (size_t)N);
    memcpy(x_ref, b, sizeof(double) * (size_t)N);
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A_ref, N, ipiv_ref, x_ref, N);
    free(ipiv_ref);
    if (info != 0) {
        free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
        return -3.0;
    }

    /* Build clusters (synthetic: DOF units, nffc=1).
     * Cluster ranges:
     *   root: nstrt=1, nsize=N
     *   clt_TL: nstrt=1, nsize=2*nb
     *   clt_BR: nstrt=1+2*nb, nsize=2*nb
     *   clt_small[i] for i=0..3: nstrt = 1 + i*nb, nsize=nb */
    st_cHACApK_cluster_t *clt_root = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_root));
    st_cHACApK_cluster_t *clt_TL   = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL));
    st_cHACApK_cluster_t *clt_BR   = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_BR));
    st_cHACApK_cluster_t *clt_small[4];
    for (int i = 0; i < 4; i++) {
        clt_small[i] = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_small[i]));
        clt_small[i]->nstrt = 1 + i * nb_small;
        clt_small[i]->nsize = nb_small;
    }
    clt_root->nstrt = 1;       clt_root->nsize = N;
    clt_TL->nstrt   = 1;       clt_TL->nsize   = half;
    clt_BR->nstrt   = 1 + half;clt_BR->nsize   = half;

    /* Track allocations for cleanup. */
    st_cHACApK_block_node_t *all_bn[16]; int n_bn = 0;
    st_cHACApK_leafmtx_t    *all_lf[16]; int n_lf = 0;

    /* Helper macro for building a dense leaf node at (row_clt, col_clt). */
    #define BUILD_LEAF(out, rclt, cclt)                                            \
        do {                                                                       \
            st_cHACApK_block_node_t *bn = (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn)); \
            st_cHACApK_leafmtx_t    *lf = (st_cHACApK_leafmtx_t*)calloc(1, sizeof(*lf));    \
            bn->row_cluster = (rclt); bn->col_cluster = (cclt);                    \
            bn->leaf_mtx = lf; bn->leaf_kind = 2;                                  \
            bn->dof_nrows = (rclt)->nsize; bn->dof_ncols = (cclt)->nsize;          \
            bn->dof_row_start = (rclt)->nstrt - 1;                                 \
            bn->dof_col_start = (cclt)->nstrt - 1;                                 \
            lf->ltmtx = 2; lf->ndl = (rclt)->nsize; lf->ndt = (cclt)->nsize;       \
            lf->nstrtl = (rclt)->nstrt; lf->nstrtt = (cclt)->nstrt;                \
            lf->a1 = (double*)malloc(sizeof(double) * (size_t)lf->ndl * (size_t)lf->ndt); \
            for (int jj = 0; jj < lf->ndt; jj++)                                   \
                for (int ii = 0; ii < lf->ndl; ii++)                               \
                    lf->a1[ii + (size_t)jj*(size_t)lf->ndl] =                      \
                        A_full[((rclt)->nstrt - 1 + ii) +                          \
                               (size_t)((cclt)->nstrt - 1 + jj)*(size_t)N];        \
            all_lf[n_lf++] = lf; all_bn[n_bn++] = bn; (out) = bn;                  \
        } while (0)

    /* Build TL internal (2x2 with 4 small leaves). */
    st_cHACApK_block_node_t *TL = (st_cHACApK_block_node_t*)calloc(1, sizeof(*TL));
    TL->row_cluster = clt_TL; TL->col_cluster = clt_TL;
    TL->nrsons = 2; TL->ncsons = 2;
    TL->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*TL->sons));
    TL->dof_nrows = half; TL->dof_ncols = half;
    TL->dof_row_start = 0; TL->dof_col_start = 0;
    BUILD_LEAF(TL->sons[0 + 0*2], clt_small[0], clt_small[0]);
    BUILD_LEAF(TL->sons[1 + 0*2], clt_small[1], clt_small[0]);
    BUILD_LEAF(TL->sons[0 + 1*2], clt_small[0], clt_small[1]);
    BUILD_LEAF(TL->sons[1 + 1*2], clt_small[1], clt_small[1]);
    all_bn[n_bn++] = TL;

    /* Build BR internal (2x2 with 4 small leaves). */
    st_cHACApK_block_node_t *BR = (st_cHACApK_block_node_t*)calloc(1, sizeof(*BR));
    BR->row_cluster = clt_BR; BR->col_cluster = clt_BR;
    BR->nrsons = 2; BR->ncsons = 2;
    BR->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*BR->sons));
    BR->dof_nrows = half; BR->dof_ncols = half;
    BR->dof_row_start = half; BR->dof_col_start = half;
    BUILD_LEAF(BR->sons[0 + 0*2], clt_small[2], clt_small[2]);
    BUILD_LEAF(BR->sons[1 + 0*2], clt_small[3], clt_small[2]);
    BUILD_LEAF(BR->sons[0 + 1*2], clt_small[2], clt_small[3]);
    BUILD_LEAF(BR->sons[1 + 1*2], clt_small[3], clt_small[3]);
    all_bn[n_bn++] = BR;

    /* Build BL leaf and TR leaf (at root level). */
    st_cHACApK_block_node_t *BL = NULL, *TR = NULL;
    BUILD_LEAF(BL, clt_BR, clt_TL);  /* BL: rows = BR's range, cols = TL's range */
    BUILD_LEAF(TR, clt_TL, clt_BR);  /* TR: rows = TL's range, cols = BR's range */

    /* Build root. */
    st_cHACApK_block_node_t *root = (st_cHACApK_block_node_t*)calloc(1, sizeof(*root));
    root->row_cluster = clt_root; root->col_cluster = clt_root;
    root->nrsons = 2; root->ncsons = 2;
    root->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*root->sons));
    root->dof_nrows = N; root->dof_ncols = N;
    root->dof_row_start = 0; root->dof_col_start = 0;
    root->sons[0 + 0*2] = TL;
    root->sons[1 + 0*2] = BL;
    root->sons[0 + 1*2] = TR;
    root->sons[1 + 1*2] = BR;
    all_bn[n_bn++] = root;

    /* Run H-LU. */
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
    for (int i = 0; i < n_lf; i++) { free(all_lf[i]->a1); free(all_lf[i]); }
    for (int i = 0; i < n_bn; i++) {
        if (all_bn[i]->sons) free(all_bn[i]->sons);
        free(all_bn[i]);
    }
    free(clt_root); free(clt_TL); free(clt_BR);
    for (int i = 0; i < 4; i++) free(clt_small[i]);
    free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);

    #undef BUILD_LEAF
    #undef RND_MS

    if (rc_dec != CHACAPK_HARITH_OK) return -4.0 + (double)rc_dec * 0.001;
    if (rc_slv != CHACAPK_HARITH_OK) return -5.0 + (double)rc_slv * 0.001;
    return rel;
}


/* ---------- Phase 4 debug: EXACT Radia structure mimic ---------- *
 *
 * Mimics the EXACT tree shape + sizes produced by Radia nx=3 leaf=10:
 *   Root (162) INTERNAL 2x2
 *   ├─ TL (108) INTERNAL 2x2
 *   │  ├─ TL.TL (72) INTERNAL 2x2 with 4 sub-leaves (48,48), (24,48), (48,24), (24,24)
 *   │  ├─ TL.BL leaf (36 x 72)
 *   │  ├─ TL.TR leaf (72 x 36)
 *   │  └─ TL.BR leaf (36 x 36)
 *   ├─ BL leaf (54 x 108)
 *   ├─ TR leaf (108 x 54)
 *   └─ BR leaf (54 x 54)
 *
 * Total: 10 leaves, 3 internal, depth 3. SAME shape as Radia's tree
 * with the SAME size asymmetry. */
double cHACApK_harith_self_test_radia_exact(void)
{
    return cHACApK_harith_self_test_radia_exact_diag(2.0);
}

double cHACApK_harith_self_test_radia_exact_with_matrix(
    const double *A_full_in, const double *b_in)
{
    /* Same Radia-exact tree shape, but A_full + b come from caller. */
    int s_root_TL = 108, s_root_BR = 54;
    int s_TL_TL = 72, s_TL_BR = 36;
    int s_TLTL_TL = 48, s_TLTL_BR = 24;
    int N = s_root_TL + s_root_BR;

    double *A_full = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *A_ref  = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *b      = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_hlu  = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_ref  = (double*)malloc(sizeof(double) * (size_t)N);
    if (!A_full || !A_ref || !b || !x_hlu || !x_ref) return -2.0;

    memcpy(A_full, A_full_in, sizeof(double) * (size_t)N * (size_t)N);
    memcpy(A_ref,  A_full_in, sizeof(double) * (size_t)N * (size_t)N);
    memcpy(b, b_in, sizeof(double) * (size_t)N);

    int *ipiv_ref = (int*)malloc(sizeof(int) * (size_t)N);
    memcpy(x_ref, b, sizeof(double) * (size_t)N);
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A_ref, N, ipiv_ref, x_ref, N);
    free(ipiv_ref);
    if (info != 0) {
        free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
        return -3.0;
    }

    /* Cluster + tree construction same as radia_exact_diag (copy-paste). */
    st_cHACApK_cluster_t *clt_root     = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_root));
    st_cHACApK_cluster_t *clt_TL       = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL));
    st_cHACApK_cluster_t *clt_BR       = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_BR));
    st_cHACApK_cluster_t *clt_TL_TL    = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL_TL));
    st_cHACApK_cluster_t *clt_TL_BR    = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL_BR));
    st_cHACApK_cluster_t *clt_TLTL_TL  = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TLTL_TL));
    st_cHACApK_cluster_t *clt_TLTL_BR  = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TLTL_BR));
    clt_root->nstrt    = 1;                 clt_root->nsize    = N;
    clt_TL->nstrt      = 1;                 clt_TL->nsize      = s_root_TL;
    clt_BR->nstrt      = 1 + s_root_TL;     clt_BR->nsize      = s_root_BR;
    clt_TL_TL->nstrt   = 1;                 clt_TL_TL->nsize   = s_TL_TL;
    clt_TL_BR->nstrt   = 1 + s_TL_TL;       clt_TL_BR->nsize   = s_TL_BR;
    clt_TLTL_TL->nstrt = 1;                 clt_TLTL_TL->nsize = s_TLTL_TL;
    clt_TLTL_BR->nstrt = 1 + s_TLTL_TL;     clt_TLTL_BR->nsize = s_TLTL_BR;

    st_cHACApK_block_node_t *all_bn[24]; int n_bn = 0;
    st_cHACApK_leafmtx_t    *all_lf[16]; int n_lf = 0;

    #define BUILD_LEAF_WM(out, rclt, cclt)                                         \
        do {                                                                       \
            st_cHACApK_block_node_t *bn = (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn)); \
            st_cHACApK_leafmtx_t    *lf = (st_cHACApK_leafmtx_t*)calloc(1, sizeof(*lf));    \
            bn->row_cluster = (rclt); bn->col_cluster = (cclt);                    \
            bn->leaf_mtx = lf; bn->leaf_kind = 2;                                  \
            bn->dof_nrows = (rclt)->nsize; bn->dof_ncols = (cclt)->nsize;          \
            bn->dof_row_start = (rclt)->nstrt - 1;                                 \
            bn->dof_col_start = (cclt)->nstrt - 1;                                 \
            lf->ltmtx = 2; lf->ndl = (rclt)->nsize; lf->ndt = (cclt)->nsize;       \
            lf->nstrtl = (rclt)->nstrt; lf->nstrtt = (cclt)->nstrt;                \
            lf->a1 = (double*)malloc(sizeof(double) * (size_t)lf->ndl * (size_t)lf->ndt); \
            for (int jj = 0; jj < lf->ndt; jj++)                                   \
                for (int ii = 0; ii < lf->ndl; ii++)                               \
                    lf->a1[ii + (size_t)jj*(size_t)lf->ndl] =                      \
                        A_full[((rclt)->nstrt - 1 + ii) +                          \
                               (size_t)((cclt)->nstrt - 1 + jj)*(size_t)N];        \
            all_lf[n_lf++] = lf; all_bn[n_bn++] = bn; (out) = bn;                  \
        } while (0)

    st_cHACApK_block_node_t *TLTL = (st_cHACApK_block_node_t*)calloc(1, sizeof(*TLTL));
    TLTL->row_cluster = clt_TL_TL; TLTL->col_cluster = clt_TL_TL;
    TLTL->nrsons = 2; TLTL->ncsons = 2;
    TLTL->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*TLTL->sons));
    TLTL->dof_nrows = s_TL_TL; TLTL->dof_ncols = s_TL_TL;
    TLTL->dof_row_start = 0; TLTL->dof_col_start = 0;
    BUILD_LEAF_WM(TLTL->sons[0 + 0*2], clt_TLTL_TL, clt_TLTL_TL);
    BUILD_LEAF_WM(TLTL->sons[1 + 0*2], clt_TLTL_BR, clt_TLTL_TL);
    BUILD_LEAF_WM(TLTL->sons[0 + 1*2], clt_TLTL_TL, clt_TLTL_BR);
    BUILD_LEAF_WM(TLTL->sons[1 + 1*2], clt_TLTL_BR, clt_TLTL_BR);
    all_bn[n_bn++] = TLTL;

    st_cHACApK_block_node_t *TL = (st_cHACApK_block_node_t*)calloc(1, sizeof(*TL));
    TL->row_cluster = clt_TL; TL->col_cluster = clt_TL;
    TL->nrsons = 2; TL->ncsons = 2;
    TL->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*TL->sons));
    TL->dof_nrows = s_root_TL; TL->dof_ncols = s_root_TL;
    TL->dof_row_start = 0; TL->dof_col_start = 0;
    TL->sons[0 + 0*2] = TLTL;
    BUILD_LEAF_WM(TL->sons[1 + 0*2], clt_TL_BR, clt_TL_TL);
    BUILD_LEAF_WM(TL->sons[0 + 1*2], clt_TL_TL, clt_TL_BR);
    BUILD_LEAF_WM(TL->sons[1 + 1*2], clt_TL_BR, clt_TL_BR);
    all_bn[n_bn++] = TL;

    st_cHACApK_block_node_t *TR = NULL, *BL = NULL, *BR = NULL;
    BUILD_LEAF_WM(BL, clt_BR, clt_TL);
    BUILD_LEAF_WM(TR, clt_TL, clt_BR);
    BUILD_LEAF_WM(BR, clt_BR, clt_BR);

    st_cHACApK_block_node_t *root = (st_cHACApK_block_node_t*)calloc(1, sizeof(*root));
    root->row_cluster = clt_root; root->col_cluster = clt_root;
    root->nrsons = 2; root->ncsons = 2;
    root->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*root->sons));
    root->dof_nrows = N; root->dof_ncols = N;
    root->dof_row_start = 0; root->dof_col_start = 0;
    root->sons[0 + 0*2] = TL;
    root->sons[1 + 0*2] = BL;
    root->sons[0 + 1*2] = TR;
    root->sons[1 + 1*2] = BR;
    all_bn[n_bn++] = root;

    int rc_dec = cHACApK_hlu_decomp(root);
    int rc_slv = (rc_dec == CHACAPK_HARITH_OK)
                  ? cHACApK_hlu_solve_vec(root, b, x_hlu, N)
                  : -99;

    double max_ref = 0.0, max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double r = (x_ref[i] < 0.0) ? -x_ref[i] : x_ref[i];
        if (r > max_ref) max_ref = r;
        double d = x_hlu[i] - x_ref[i]; if (d < 0.0) d = -d;
        if (d > max_err) max_err = d;
    }
    double rel = (max_ref > 0.0) ? (max_err / max_ref) : max_err;

    for (int i = 0; i < n_lf; i++) { free(all_lf[i]->a1); free(all_lf[i]); }
    for (int i = 0; i < n_bn; i++) {
        if (all_bn[i]->sons) free(all_bn[i]->sons);
        free(all_bn[i]);
    }
    free(clt_root); free(clt_TL); free(clt_BR);
    free(clt_TL_TL); free(clt_TL_BR); free(clt_TLTL_TL); free(clt_TLTL_BR);
    free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);

    #undef BUILD_LEAF_WM

    if (rc_dec != CHACAPK_HARITH_OK) return -4.0 + (double)rc_dec * 0.001;
    if (rc_slv != CHACAPK_HARITH_OK) return -5.0 + (double)rc_slv * 0.001;
    return rel;
}

double cHACApK_harith_self_test_radia_exact_diag(double diag_boost)
{
    /* Match Radia exactly. diag_boost = additional diagonal magnitude over
     * row sum. Default 2.0 = mildly diag-dominant (similar to MSC).
     * Real Radia MSC has even weaker dominance. */
    int s_root_TL = 108, s_root_BR = 54;
    int s_TL_TL = 72, s_TL_BR = 36;
    int s_TLTL_TL = 48, s_TLTL_BR = 24;
    int N = s_root_TL + s_root_BR;  /* = 162 */

    unsigned long seed = 0xABCDEF42UL;
    #define RND_RX() (seed = seed * 6364136223846793005UL + 1442695040888963407UL, \
                      (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5)

    double *A_full = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *A_ref  = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *b      = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_hlu  = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_ref  = (double*)malloc(sizeof(double) * (size_t)N);
    if (!A_full || !A_ref || !b || !x_hlu || !x_ref) return -2.0;

    for (int j = 0; j < N; j++) {
        double rs = 0.0;
        for (int i = 0; i < N; i++) {
            double v = RND_RX();
            A_full[i + (size_t)j*(size_t)N] = v;
            if (i != j) rs += (v < 0.0 ? -v : v);
        }
        A_full[j + (size_t)j*(size_t)N] = rs + diag_boost;
    }
    for (int j = 0; j < N; j++) b[j] = RND_RX() * 10.0;
    memcpy(A_ref, A_full, sizeof(double) * (size_t)N * (size_t)N);

    int *ipiv_ref = (int*)malloc(sizeof(int) * (size_t)N);
    memcpy(x_ref, b, sizeof(double) * (size_t)N);
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A_ref, N, ipiv_ref, x_ref, N);
    free(ipiv_ref);
    if (info != 0) {
        free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
        return -3.0;
    }

    /* Clusters. */
    st_cHACApK_cluster_t *clt_root     = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_root));
    st_cHACApK_cluster_t *clt_TL       = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL));
    st_cHACApK_cluster_t *clt_BR       = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_BR));
    st_cHACApK_cluster_t *clt_TL_TL    = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL_TL));
    st_cHACApK_cluster_t *clt_TL_BR    = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL_BR));
    st_cHACApK_cluster_t *clt_TLTL_TL  = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TLTL_TL));
    st_cHACApK_cluster_t *clt_TLTL_BR  = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TLTL_BR));
    clt_root->nstrt    = 1;                 clt_root->nsize    = N;
    clt_TL->nstrt      = 1;                 clt_TL->nsize      = s_root_TL;
    clt_BR->nstrt      = 1 + s_root_TL;     clt_BR->nsize      = s_root_BR;
    clt_TL_TL->nstrt   = 1;                 clt_TL_TL->nsize   = s_TL_TL;
    clt_TL_BR->nstrt   = 1 + s_TL_TL;       clt_TL_BR->nsize   = s_TL_BR;
    clt_TLTL_TL->nstrt = 1;                 clt_TLTL_TL->nsize = s_TLTL_TL;
    clt_TLTL_BR->nstrt = 1 + s_TLTL_TL;     clt_TLTL_BR->nsize = s_TLTL_BR;

    st_cHACApK_block_node_t *all_bn[24]; int n_bn = 0;
    st_cHACApK_leafmtx_t    *all_lf[16]; int n_lf = 0;

    #define BUILD_LEAF_RX(out, rclt, cclt)                                         \
        do {                                                                       \
            st_cHACApK_block_node_t *bn = (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn)); \
            st_cHACApK_leafmtx_t    *lf = (st_cHACApK_leafmtx_t*)calloc(1, sizeof(*lf));    \
            bn->row_cluster = (rclt); bn->col_cluster = (cclt);                    \
            bn->leaf_mtx = lf; bn->leaf_kind = 2;                                  \
            bn->dof_nrows = (rclt)->nsize; bn->dof_ncols = (cclt)->nsize;          \
            bn->dof_row_start = (rclt)->nstrt - 1;                                 \
            bn->dof_col_start = (cclt)->nstrt - 1;                                 \
            lf->ltmtx = 2; lf->ndl = (rclt)->nsize; lf->ndt = (cclt)->nsize;       \
            lf->nstrtl = (rclt)->nstrt; lf->nstrtt = (cclt)->nstrt;                \
            lf->a1 = (double*)malloc(sizeof(double) * (size_t)lf->ndl * (size_t)lf->ndt); \
            for (int jj = 0; jj < lf->ndt; jj++)                                   \
                for (int ii = 0; ii < lf->ndl; ii++)                               \
                    lf->a1[ii + (size_t)jj*(size_t)lf->ndl] =                      \
                        A_full[((rclt)->nstrt - 1 + ii) +                          \
                               (size_t)((cclt)->nstrt - 1 + jj)*(size_t)N];        \
            all_lf[n_lf++] = lf; all_bn[n_bn++] = bn; (out) = bn;                  \
        } while (0)

    /* TL.TL (depth 2 internal, 72x72 split into 48 + 24). */
    st_cHACApK_block_node_t *TLTL = (st_cHACApK_block_node_t*)calloc(1, sizeof(*TLTL));
    TLTL->row_cluster = clt_TL_TL; TLTL->col_cluster = clt_TL_TL;
    TLTL->nrsons = 2; TLTL->ncsons = 2;
    TLTL->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*TLTL->sons));
    TLTL->dof_nrows = s_TL_TL; TLTL->dof_ncols = s_TL_TL;
    TLTL->dof_row_start = 0; TLTL->dof_col_start = 0;
    BUILD_LEAF_RX(TLTL->sons[0 + 0*2], clt_TLTL_TL, clt_TLTL_TL);  /* 48x48 */
    BUILD_LEAF_RX(TLTL->sons[1 + 0*2], clt_TLTL_BR, clt_TLTL_TL);  /* 24x48 */
    BUILD_LEAF_RX(TLTL->sons[0 + 1*2], clt_TLTL_TL, clt_TLTL_BR);  /* 48x24 */
    BUILD_LEAF_RX(TLTL->sons[1 + 1*2], clt_TLTL_BR, clt_TLTL_BR);  /* 24x24 */
    all_bn[n_bn++] = TLTL;

    /* TL (depth 1 internal, 108x108 split into 72 + 36). */
    st_cHACApK_block_node_t *TL = (st_cHACApK_block_node_t*)calloc(1, sizeof(*TL));
    TL->row_cluster = clt_TL; TL->col_cluster = clt_TL;
    TL->nrsons = 2; TL->ncsons = 2;
    TL->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*TL->sons));
    TL->dof_nrows = s_root_TL; TL->dof_ncols = s_root_TL;
    TL->dof_row_start = 0; TL->dof_col_start = 0;
    TL->sons[0 + 0*2] = TLTL;
    BUILD_LEAF_RX(TL->sons[1 + 0*2], clt_TL_BR, clt_TL_TL);  /* 36x72 */
    BUILD_LEAF_RX(TL->sons[0 + 1*2], clt_TL_TL, clt_TL_BR);  /* 72x36 */
    BUILD_LEAF_RX(TL->sons[1 + 1*2], clt_TL_BR, clt_TL_BR);  /* 36x36 */
    all_bn[n_bn++] = TL;

    /* Root-level leaves. */
    st_cHACApK_block_node_t *TR = NULL, *BL = NULL, *BR = NULL;
    BUILD_LEAF_RX(BL, clt_BR, clt_TL);  /* 54x108 */
    BUILD_LEAF_RX(TR, clt_TL, clt_BR);  /* 108x54 */
    BUILD_LEAF_RX(BR, clt_BR, clt_BR);  /* 54x54 */

    /* Root. */
    st_cHACApK_block_node_t *root = (st_cHACApK_block_node_t*)calloc(1, sizeof(*root));
    root->row_cluster = clt_root; root->col_cluster = clt_root;
    root->nrsons = 2; root->ncsons = 2;
    root->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*root->sons));
    root->dof_nrows = N; root->dof_ncols = N;
    root->dof_row_start = 0; root->dof_col_start = 0;
    root->sons[0 + 0*2] = TL;
    root->sons[1 + 0*2] = BL;
    root->sons[0 + 1*2] = TR;
    root->sons[1 + 1*2] = BR;
    all_bn[n_bn++] = root;

    int rc_dec = cHACApK_hlu_decomp(root);
    int rc_slv = (rc_dec == CHACAPK_HARITH_OK)
                  ? cHACApK_hlu_solve_vec(root, b, x_hlu, N)
                  : -99;

    double max_ref = 0.0, max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double r = (x_ref[i] < 0.0) ? -x_ref[i] : x_ref[i];
        if (r > max_ref) max_ref = r;
        double d = x_hlu[i] - x_ref[i]; if (d < 0.0) d = -d;
        if (d > max_err) max_err = d;
    }
    double rel = (max_ref > 0.0) ? (max_err / max_ref) : max_err;

    for (int i = 0; i < n_lf; i++) { free(all_lf[i]->a1); free(all_lf[i]); }
    for (int i = 0; i < n_bn; i++) {
        if (all_bn[i]->sons) free(all_bn[i]->sons);
        free(all_bn[i]);
    }
    free(clt_root); free(clt_TL); free(clt_BR);
    free(clt_TL_TL); free(clt_TL_BR); free(clt_TLTL_TL); free(clt_TLTL_BR);
    free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);

    #undef BUILD_LEAF_RX
    #undef RND_RX

    if (rc_dec != CHACAPK_HARITH_OK) return -4.0 + (double)rc_dec * 0.001;
    if (rc_slv != CHACAPK_HARITH_OK) return -5.0 + (double)rc_slv * 0.001;
    return rel;
}


/* ---------- Phase 4 debug: depth-3 asymmetric test ------------------- *
 *
 * Mimics Radia's exact tree shape: 10 leaves, 3 internal, depth 3.
 * Root: 2x2. TL: 2x2 internal (deep). TL.TL: 2x2 internal (deepest, 4 leaves).
 * TL.TR/TL.BL/TL.BR: leaves at depth 2. TR/BL/BR at root level: leaves at depth 1.
 *
 * N = 4 * nb_tiny (each TL.TL sub-leaf is nb_tiny x nb_tiny). */
double cHACApK_harith_self_test_depth3_asymmetric(int nb_tiny)
{
    if (nb_tiny <= 0 || nb_tiny > 32) return -1.0;
    int N = 4 * nb_tiny;       /* total */
    int q = nb_tiny;           /* deepest sub-leaf */
    int half = 2 * q;          /* TL.TL sub-block (=internal at depth 2) */
    /* TL has size 2*half = 4*q = N/?... wait, TL size = ? Let me think.
     * If TL.TL = half x half (= 2q x 2q), and TL has 2x2 children where TL.TL
     * is the (0,0) sub-block, the others (TL.TR, TL.BL, TL.BR) must have
     * matching sizes for TL to be valid: TL would be (2*half) x (2*half) = 4q x 4q.
     * Hmm, but then root TL is 4q x 4q, root BR is also some size, and N = TL.size + BR.size = 8q.
     *
     * Simpler choice: TL = 2*half, with 4 children all of size half (2q each).
     * - TL.TL: internal 2x2 (4 sub-leaves of size q each) -- so TL.TL = 2q x 2q.
     * - TL.TR/TL.BL/TL.BR: leaves of size 2q x 2q.
     * Then TL = 4q x 4q.
     *
     * Root: 4 children of size half_root each. If TL = 4q, root = 8q?
     * Or simpler: root has TL=4q x 4q, BR=4q x 4q, TR=4q x 4q, BL=4q x 4q.
     * Then N = 8q.
     */
    int qq = q;                 /* deepest sub-leaf (renamed for clarity) */
    int s2 = 2 * qq;            /* TL.TR/TL.BL/TL.BR leaf size (= half) */
    int s1 = 2 * s2;            /* TL itself = 2*s2 = 4*qq; root child size */
    /* root = 2*s1 = 8*qq. So N must be 8*qq. */
    N = 8 * qq;

    unsigned long seed = 0xBEEFCAFEUL;
    #define RND_D3() (seed = seed * 6364136223846793005UL + 1442695040888963407UL, \
                      (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5)

    double *A_full = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *A_ref  = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *b      = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_hlu  = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_ref  = (double*)malloc(sizeof(double) * (size_t)N);
    if (!A_full || !A_ref || !b || !x_hlu || !x_ref) return -2.0;

    for (int j = 0; j < N; j++) {
        double rs = 0.0;
        for (int i = 0; i < N; i++) {
            double v = RND_D3();
            A_full[i + (size_t)j*(size_t)N] = v;
            if (i != j) rs += (v < 0.0 ? -v : v);
        }
        A_full[j + (size_t)j*(size_t)N] = rs + (double)N;
    }
    for (int j = 0; j < N; j++) b[j] = RND_D3() * 10.0;
    memcpy(A_ref, A_full, sizeof(double) * (size_t)N * (size_t)N);

    int *ipiv_ref = (int*)malloc(sizeof(int) * (size_t)N);
    memcpy(x_ref, b, sizeof(double) * (size_t)N);
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A_ref, N, ipiv_ref, x_ref, N);
    free(ipiv_ref);
    if (info != 0) {
        free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
        return -3.0;
    }

    /* Cluster tree:
     *   root: nstrt=1, nsize=N (=8q)
     *   clt_root_TL: 1..s1 (TL = 4q)
     *   clt_root_BR: 1+s1..N (BR = 4q)
     *   Inside TL:
     *     clt_TL_TL: 1..s2 (= 2q) -- this is INTERNAL, splits further
     *     clt_TL_TR: 1+s2..s1
     *     clt_TL_BL: same as TL.TR
     *     clt_TL_BR: same as TL.TR
     *   Inside TL.TL: 4 sub-clusters of size q each.
     */
    st_cHACApK_cluster_t *clt_root   = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_root));
    st_cHACApK_cluster_t *clt_TL     = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL));   /* root sons[0,0] */
    st_cHACApK_cluster_t *clt_BR     = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_BR));   /* root sons[1,1] */
    st_cHACApK_cluster_t *clt_TLTL   = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TLTL)); /* TL sons[0,0], internal */
    st_cHACApK_cluster_t *clt_TLBR   = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TLBR));
    /* TL.TL sub-clusters (4 small): nstrt = 1, 1+q, 1+2q, ... */
    st_cHACApK_cluster_t *clt_q[4];
    for (int i = 0; i < 4; i++) {
        clt_q[i] = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_q[i]));
        clt_q[i]->nstrt = 1 + i * qq;
        clt_q[i]->nsize = qq;
    }
    clt_root->nstrt = 1;        clt_root->nsize = N;
    clt_TL->nstrt   = 1;        clt_TL->nsize   = s1;
    clt_BR->nstrt   = 1 + s1;   clt_BR->nsize   = s1;
    clt_TLTL->nstrt = 1;        clt_TLTL->nsize = s2;
    clt_TLBR->nstrt = 1 + s2;   clt_TLBR->nsize = s2;

    st_cHACApK_block_node_t *all_bn[24]; int n_bn = 0;
    st_cHACApK_leafmtx_t    *all_lf[16]; int n_lf = 0;

    #define BUILD_LEAF_D3(out, rclt, cclt)                                         \
        do {                                                                       \
            st_cHACApK_block_node_t *bn = (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn)); \
            st_cHACApK_leafmtx_t    *lf = (st_cHACApK_leafmtx_t*)calloc(1, sizeof(*lf));    \
            bn->row_cluster = (rclt); bn->col_cluster = (cclt);                    \
            bn->leaf_mtx = lf; bn->leaf_kind = 2;                                  \
            bn->dof_nrows = (rclt)->nsize; bn->dof_ncols = (cclt)->nsize;          \
            bn->dof_row_start = (rclt)->nstrt - 1;                                 \
            bn->dof_col_start = (cclt)->nstrt - 1;                                 \
            lf->ltmtx = 2; lf->ndl = (rclt)->nsize; lf->ndt = (cclt)->nsize;       \
            lf->nstrtl = (rclt)->nstrt; lf->nstrtt = (cclt)->nstrt;                \
            lf->a1 = (double*)malloc(sizeof(double) * (size_t)lf->ndl * (size_t)lf->ndt); \
            for (int jj = 0; jj < lf->ndt; jj++)                                   \
                for (int ii = 0; ii < lf->ndl; ii++)                               \
                    lf->a1[ii + (size_t)jj*(size_t)lf->ndl] =                      \
                        A_full[((rclt)->nstrt - 1 + ii) +                          \
                               (size_t)((cclt)->nstrt - 1 + jj)*(size_t)N];        \
            all_lf[n_lf++] = lf; all_bn[n_bn++] = bn; (out) = bn;                  \
        } while (0)

    /* TL.TL (depth 2 internal): 4 sub-leaves of size q. */
    st_cHACApK_block_node_t *TLTL = (st_cHACApK_block_node_t*)calloc(1, sizeof(*TLTL));
    TLTL->row_cluster = clt_TLTL; TLTL->col_cluster = clt_TLTL;
    TLTL->nrsons = 2; TLTL->ncsons = 2;
    TLTL->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*TLTL->sons));
    TLTL->dof_nrows = s2; TLTL->dof_ncols = s2;
    TLTL->dof_row_start = 0; TLTL->dof_col_start = 0;
    BUILD_LEAF_D3(TLTL->sons[0 + 0*2], clt_q[0], clt_q[0]);
    BUILD_LEAF_D3(TLTL->sons[1 + 0*2], clt_q[1], clt_q[0]);
    BUILD_LEAF_D3(TLTL->sons[0 + 1*2], clt_q[0], clt_q[1]);
    BUILD_LEAF_D3(TLTL->sons[1 + 1*2], clt_q[1], clt_q[1]);
    all_bn[n_bn++] = TLTL;

    /* TL (depth 1 internal): TL.TL = TLTL, others are leaves at depth 2. */
    st_cHACApK_block_node_t *TL = (st_cHACApK_block_node_t*)calloc(1, sizeof(*TL));
    TL->row_cluster = clt_TL; TL->col_cluster = clt_TL;
    TL->nrsons = 2; TL->ncsons = 2;
    TL->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*TL->sons));
    TL->dof_nrows = s1; TL->dof_ncols = s1;
    TL->dof_row_start = 0; TL->dof_col_start = 0;
    TL->sons[0 + 0*2] = TLTL;
    BUILD_LEAF_D3(TL->sons[1 + 0*2], clt_TLBR, clt_TLTL);  /* TL.BL: row range = TLBR, col range = TLTL */
    BUILD_LEAF_D3(TL->sons[0 + 1*2], clt_TLTL, clt_TLBR);  /* TL.TR */
    BUILD_LEAF_D3(TL->sons[1 + 1*2], clt_TLBR, clt_TLBR);  /* TL.BR */
    all_bn[n_bn++] = TL;

    /* TR / BL / BR (all leaves at root level). */
    st_cHACApK_block_node_t *TR = NULL, *BL = NULL, *BR = NULL;
    BUILD_LEAF_D3(BL, clt_BR, clt_TL);
    BUILD_LEAF_D3(TR, clt_TL, clt_BR);
    BUILD_LEAF_D3(BR, clt_BR, clt_BR);

    /* Root */
    st_cHACApK_block_node_t *root = (st_cHACApK_block_node_t*)calloc(1, sizeof(*root));
    root->row_cluster = clt_root; root->col_cluster = clt_root;
    root->nrsons = 2; root->ncsons = 2;
    root->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*root->sons));
    root->dof_nrows = N; root->dof_ncols = N;
    root->dof_row_start = 0; root->dof_col_start = 0;
    root->sons[0 + 0*2] = TL;
    root->sons[1 + 0*2] = BL;
    root->sons[0 + 1*2] = TR;
    root->sons[1 + 1*2] = BR;
    all_bn[n_bn++] = root;

    int rc_dec = cHACApK_hlu_decomp(root);
    int rc_slv = (rc_dec == CHACAPK_HARITH_OK)
                  ? cHACApK_hlu_solve_vec(root, b, x_hlu, N)
                  : -99;

    double max_ref = 0.0, max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double r = (x_ref[i] < 0.0) ? -x_ref[i] : x_ref[i];
        if (r > max_ref) max_ref = r;
        double d = x_hlu[i] - x_ref[i]; if (d < 0.0) d = -d;
        if (d > max_err) max_err = d;
    }
    double rel = (max_ref > 0.0) ? (max_err / max_ref) : max_err;

    for (int i = 0; i < n_lf; i++) { free(all_lf[i]->a1); free(all_lf[i]); }
    for (int i = 0; i < n_bn; i++) {
        if (all_bn[i]->sons) free(all_bn[i]->sons);
        free(all_bn[i]);
    }
    free(clt_root); free(clt_TL); free(clt_BR); free(clt_TLTL); free(clt_TLBR);
    for (int i = 0; i < 4; i++) free(clt_q[i]);
    free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);

    #undef BUILD_LEAF_D3
    #undef RND_D3

    if (rc_dec != CHACAPK_HARITH_OK) return -4.0 + (double)rc_dec * 0.001;
    if (rc_slv != CHACAPK_HARITH_OK) return -5.0 + (double)rc_slv * 0.001;
    return rel;
}


/* ---------- Phase 4 debug: conversion-path mixed-sibling test --------- *
 *
 * Mimics the EXACT data path used by cHACApK_hlu_run_on_hacapk on real
 * Radia trees: build dense leaves in HACApK row-major format, then call
 * the same transpose path used by cHACApK_convert_leafmtxp_to_internal,
 * then run H-LU.
 *
 * Walks the block-tree directly (not a HACApK leafmtxp) since the
 * synthetic test doesn't have one. */
static void transpose_dense_leaves_recursive(st_cHACApK_block_node_t *node)
{
    if (!node) return;
    if (leaf_is_dense(node)) {
        int ndl = node->leaf_mtx->ndl;
        int ndt = node->leaf_mtx->ndt;
        double *new_a1 = (double*)malloc(sizeof(double) * (size_t)ndl * (size_t)ndt);
        if (!new_a1) return;
        /* HACApK: a1[col + row*ndt] = M[row, col]  ->
         * internal: new_a1[row + col*ndl] = M[row, col] */
        for (int row = 0; row < ndl; row++)
            for (int col = 0; col < ndt; col++)
                new_a1[row + (size_t)col * (size_t)ndl] =
                    node->leaf_mtx->a1[col + (size_t)row * (size_t)ndt];
        free(node->leaf_mtx->a1);
        node->leaf_mtx->a1 = new_a1;
        return;
    }
    if (leaf_is_rk(node)) {
        /* rk: swap a1 (V) <-> a2 (U). */
        double *tmp = node->leaf_mtx->a1;
        node->leaf_mtx->a1 = node->leaf_mtx->a2;
        node->leaf_mtx->a2 = tmp;
        return;
    }
    int nsons = node->nrsons * node->ncsons;
    for (int i = 0; i < nsons; i++) transpose_dense_leaves_recursive(node->sons[i]);
}

double cHACApK_harith_self_test_mixed_sibling_via_conversion(int nb_small)
{
    if (nb_small <= 0 || nb_small > 64) return -1.0;
    int N = 4 * nb_small;
    int half = 2 * nb_small;

    unsigned long seed = 0xCAFEBABEUL;
    #define RND_VC() (seed = seed * 6364136223846793005UL + 1442695040888963407UL, \
                      (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5)

    double *A_full = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *A_ref  = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *b      = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_hlu  = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_ref  = (double*)malloc(sizeof(double) * (size_t)N);
    if (!A_full || !A_ref || !b || !x_hlu || !x_ref) return -2.0;

    for (int j = 0; j < N; j++) {
        double rs = 0.0;
        for (int i = 0; i < N; i++) {
            double v = RND_VC();
            A_full[i + (size_t)j*(size_t)N] = v;
            if (i != j) rs += (v < 0.0 ? -v : v);
        }
        A_full[j + (size_t)j*(size_t)N] = rs + (double)N;
    }
    for (int j = 0; j < N; j++) b[j] = RND_VC() * 10.0;
    memcpy(A_ref, A_full, sizeof(double) * (size_t)N * (size_t)N);

    int *ipiv_ref = (int*)malloc(sizeof(int) * (size_t)N);
    memcpy(x_ref, b, sizeof(double) * (size_t)N);
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A_ref, N, ipiv_ref, x_ref, N);
    free(ipiv_ref);
    if (info != 0) {
        free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
        return -3.0;
    }

    /* Build clusters. */
    st_cHACApK_cluster_t *clt_root = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_root));
    st_cHACApK_cluster_t *clt_TL   = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL));
    st_cHACApK_cluster_t *clt_BR   = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_BR));
    st_cHACApK_cluster_t *clt_small[4];
    for (int i = 0; i < 4; i++) {
        clt_small[i] = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_small[i]));
        clt_small[i]->nstrt = 1 + i * nb_small;
        clt_small[i]->nsize = nb_small;
    }
    clt_root->nstrt = 1;       clt_root->nsize = N;
    clt_TL->nstrt   = 1;       clt_TL->nsize   = half;
    clt_BR->nstrt   = 1 + half; clt_BR->nsize   = half;

    st_cHACApK_block_node_t *all_bn[16]; int n_bn = 0;
    st_cHACApK_leafmtx_t    *all_lf[16]; int n_lf = 0;

    /* Build leaf in HACApK row-major: a1[col + row*ndt] = M[row, col]. */
    #define BUILD_LEAF_HACAPK(out, rclt, cclt)                                     \
        do {                                                                       \
            st_cHACApK_block_node_t *bn = (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn)); \
            st_cHACApK_leafmtx_t    *lf = (st_cHACApK_leafmtx_t*)calloc(1, sizeof(*lf));    \
            bn->row_cluster = (rclt); bn->col_cluster = (cclt);                    \
            bn->leaf_mtx = lf; bn->leaf_kind = 2;                                  \
            bn->dof_nrows = (rclt)->nsize; bn->dof_ncols = (cclt)->nsize;          \
            bn->dof_row_start = (rclt)->nstrt - 1;                                 \
            bn->dof_col_start = (cclt)->nstrt - 1;                                 \
            lf->ltmtx = 2; lf->ndl = (rclt)->nsize; lf->ndt = (cclt)->nsize;       \
            lf->nstrtl = (rclt)->nstrt; lf->nstrtt = (cclt)->nstrt;                \
            lf->a1 = (double*)malloc(sizeof(double) * (size_t)lf->ndl * (size_t)lf->ndt); \
            for (int rr = 0; rr < lf->ndl; rr++)                                   \
                for (int cc = 0; cc < lf->ndt; cc++)                               \
                    lf->a1[cc + (size_t)rr*(size_t)lf->ndt] =                      \
                        A_full[((rclt)->nstrt - 1 + rr) +                          \
                               (size_t)((cclt)->nstrt - 1 + cc)*(size_t)N];        \
            all_lf[n_lf++] = lf; all_bn[n_bn++] = bn; (out) = bn;                  \
        } while (0)

    st_cHACApK_block_node_t *TL = (st_cHACApK_block_node_t*)calloc(1, sizeof(*TL));
    TL->row_cluster = clt_TL; TL->col_cluster = clt_TL;
    TL->nrsons = 2; TL->ncsons = 2;
    TL->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*TL->sons));
    TL->dof_nrows = half; TL->dof_ncols = half;
    TL->dof_row_start = 0; TL->dof_col_start = 0;
    BUILD_LEAF_HACAPK(TL->sons[0 + 0*2], clt_small[0], clt_small[0]);
    BUILD_LEAF_HACAPK(TL->sons[1 + 0*2], clt_small[1], clt_small[0]);
    BUILD_LEAF_HACAPK(TL->sons[0 + 1*2], clt_small[0], clt_small[1]);
    BUILD_LEAF_HACAPK(TL->sons[1 + 1*2], clt_small[1], clt_small[1]);
    all_bn[n_bn++] = TL;

    st_cHACApK_block_node_t *BR = (st_cHACApK_block_node_t*)calloc(1, sizeof(*BR));
    BR->row_cluster = clt_BR; BR->col_cluster = clt_BR;
    BR->nrsons = 2; BR->ncsons = 2;
    BR->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*BR->sons));
    BR->dof_nrows = half; BR->dof_ncols = half;
    BR->dof_row_start = half; BR->dof_col_start = half;
    BUILD_LEAF_HACAPK(BR->sons[0 + 0*2], clt_small[2], clt_small[2]);
    BUILD_LEAF_HACAPK(BR->sons[1 + 0*2], clt_small[3], clt_small[2]);
    BUILD_LEAF_HACAPK(BR->sons[0 + 1*2], clt_small[2], clt_small[3]);
    BUILD_LEAF_HACAPK(BR->sons[1 + 1*2], clt_small[3], clt_small[3]);
    all_bn[n_bn++] = BR;

    st_cHACApK_block_node_t *BL = NULL, *TR = NULL;
    BUILD_LEAF_HACAPK(BL, clt_BR, clt_TL);
    BUILD_LEAF_HACAPK(TR, clt_TL, clt_BR);

    st_cHACApK_block_node_t *root = (st_cHACApK_block_node_t*)calloc(1, sizeof(*root));
    root->row_cluster = clt_root; root->col_cluster = clt_root;
    root->nrsons = 2; root->ncsons = 2;
    root->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*root->sons));
    root->dof_nrows = N; root->dof_ncols = N;
    root->dof_row_start = 0; root->dof_col_start = 0;
    root->sons[0 + 0*2] = TL;
    root->sons[1 + 0*2] = BL;
    root->sons[0 + 1*2] = TR;
    root->sons[1 + 1*2] = BR;
    all_bn[n_bn++] = root;

    /* HACApK->internal: walk tree, transpose each dense leaf. */
    transpose_dense_leaves_recursive(root);

    /* Run H-LU. */
    int rc_dec = cHACApK_hlu_decomp(root);
    int rc_slv = (rc_dec == CHACAPK_HARITH_OK)
                  ? cHACApK_hlu_solve_vec(root, b, x_hlu, N)
                  : -99;

    double max_ref = 0.0, max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double r = (x_ref[i] < 0.0) ? -x_ref[i] : x_ref[i];
        if (r > max_ref) max_ref = r;
        double d = x_hlu[i] - x_ref[i]; if (d < 0.0) d = -d;
        if (d > max_err) max_err = d;
    }
    double rel = (max_ref > 0.0) ? (max_err / max_ref) : max_err;

    for (int i = 0; i < n_lf; i++) { free(all_lf[i]->a1); free(all_lf[i]); }
    for (int i = 0; i < n_bn; i++) {
        if (all_bn[i]->sons) free(all_bn[i]->sons);
        free(all_bn[i]);
    }
    free(clt_root); free(clt_TL); free(clt_BR);
    for (int i = 0; i < 4; i++) free(clt_small[i]);
    free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);

    #undef BUILD_LEAF_HACAPK
    #undef RND_VC

    if (rc_dec != CHACAPK_HARITH_OK) return -4.0 + (double)rc_dec * 0.001;
    if (rc_slv != CHACAPK_HARITH_OK) return -5.0 + (double)rc_slv * 0.001;
    return rel;
}


/* ---------- Phase 4 debug: non-uniform mixed-sibling test ------------- *
 *
 * Same shape as mixed_sibling but with asymmetric cluster splits.
 * Reproduces HACApK's element-count split pattern (e.g., 13 -> 6+7). */
double cHACApK_harith_self_test_mixed_sibling_nonuniform(
    int n1, int n2, int m1, int m3)
{
    if (n1 <= 0 || n2 <= 0 || m1 <= 0 || m3 <= 0) return -1.0;
    if (m1 >= n1 || m3 >= n2) return -1.0;
    int N = n1 + n2;

    unsigned long seed = 0xDEADBEEFUL;
    #define RND_NU() (seed = seed * 6364136223846793005UL + 1442695040888963407UL, \
                      (double)((seed >> 33) & 0x7fffffff) / 2147483647.0 - 0.5)

    double *A_full = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *A_ref  = (double*)malloc(sizeof(double) * (size_t)N * (size_t)N);
    double *b      = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_hlu  = (double*)malloc(sizeof(double) * (size_t)N);
    double *x_ref  = (double*)malloc(sizeof(double) * (size_t)N);
    if (!A_full || !A_ref || !b || !x_hlu || !x_ref) return -2.0;

    for (int j = 0; j < N; j++) {
        double rs = 0.0;
        for (int i = 0; i < N; i++) {
            double v = RND_NU();
            A_full[i + (size_t)j*(size_t)N] = v;
            if (i != j) rs += (v < 0.0 ? -v : v);
        }
        A_full[j + (size_t)j*(size_t)N] = rs + (double)N;
    }
    for (int j = 0; j < N; j++) b[j] = RND_NU() * 10.0;
    memcpy(A_ref, A_full, sizeof(double) * (size_t)N * (size_t)N);

    int *ipiv_ref = (int*)malloc(sizeof(int) * (size_t)N);
    memcpy(x_ref, b, sizeof(double) * (size_t)N);
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A_ref, N, ipiv_ref, x_ref, N);
    free(ipiv_ref);
    if (info != 0) {
        free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);
        return -3.0;
    }

    /* Build clusters. Each sub-cluster has different sizes. */
    st_cHACApK_cluster_t *clt_root = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_root));
    st_cHACApK_cluster_t *clt_TL   = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL));
    st_cHACApK_cluster_t *clt_BR   = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_BR));
    st_cHACApK_cluster_t *clt_TL_0 = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL_0));
    st_cHACApK_cluster_t *clt_TL_1 = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_TL_1));
    st_cHACApK_cluster_t *clt_BR_0 = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_BR_0));
    st_cHACApK_cluster_t *clt_BR_1 = (st_cHACApK_cluster_t*)calloc(1, sizeof(*clt_BR_1));

    clt_root->nstrt = 1;          clt_root->nsize = N;
    clt_TL->nstrt   = 1;          clt_TL->nsize   = n1;
    clt_BR->nstrt   = 1 + n1;     clt_BR->nsize   = n2;
    clt_TL_0->nstrt = 1;          clt_TL_0->nsize = m1;
    clt_TL_1->nstrt = 1 + m1;     clt_TL_1->nsize = n1 - m1;
    clt_BR_0->nstrt = 1 + n1;     clt_BR_0->nsize = m3;
    clt_BR_1->nstrt = 1 + n1 + m3; clt_BR_1->nsize = n2 - m3;

    st_cHACApK_block_node_t *all_bn[16]; int n_bn = 0;
    st_cHACApK_leafmtx_t    *all_lf[16]; int n_lf = 0;

    #define BUILD_LEAF_NU(out, rclt, cclt)                                         \
        do {                                                                       \
            st_cHACApK_block_node_t *bn = (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn)); \
            st_cHACApK_leafmtx_t    *lf = (st_cHACApK_leafmtx_t*)calloc(1, sizeof(*lf));    \
            bn->row_cluster = (rclt); bn->col_cluster = (cclt);                    \
            bn->leaf_mtx = lf; bn->leaf_kind = 2;                                  \
            bn->dof_nrows = (rclt)->nsize; bn->dof_ncols = (cclt)->nsize;          \
            bn->dof_row_start = (rclt)->nstrt - 1;                                 \
            bn->dof_col_start = (cclt)->nstrt - 1;                                 \
            lf->ltmtx = 2; lf->ndl = (rclt)->nsize; lf->ndt = (cclt)->nsize;       \
            lf->nstrtl = (rclt)->nstrt; lf->nstrtt = (cclt)->nstrt;                \
            lf->a1 = (double*)malloc(sizeof(double) * (size_t)lf->ndl * (size_t)lf->ndt); \
            for (int jj = 0; jj < lf->ndt; jj++)                                   \
                for (int ii = 0; ii < lf->ndl; ii++)                               \
                    lf->a1[ii + (size_t)jj*(size_t)lf->ndl] =                      \
                        A_full[((rclt)->nstrt - 1 + ii) +                          \
                               (size_t)((cclt)->nstrt - 1 + jj)*(size_t)N];        \
            all_lf[n_lf++] = lf; all_bn[n_bn++] = bn; (out) = bn;                  \
        } while (0)

    /* TL internal */
    st_cHACApK_block_node_t *TL = (st_cHACApK_block_node_t*)calloc(1, sizeof(*TL));
    TL->row_cluster = clt_TL; TL->col_cluster = clt_TL;
    TL->nrsons = 2; TL->ncsons = 2;
    TL->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*TL->sons));
    TL->dof_nrows = n1; TL->dof_ncols = n1;
    TL->dof_row_start = 0; TL->dof_col_start = 0;
    BUILD_LEAF_NU(TL->sons[0 + 0*2], clt_TL_0, clt_TL_0);
    BUILD_LEAF_NU(TL->sons[1 + 0*2], clt_TL_1, clt_TL_0);
    BUILD_LEAF_NU(TL->sons[0 + 1*2], clt_TL_0, clt_TL_1);
    BUILD_LEAF_NU(TL->sons[1 + 1*2], clt_TL_1, clt_TL_1);
    all_bn[n_bn++] = TL;

    /* BR internal */
    st_cHACApK_block_node_t *BR = (st_cHACApK_block_node_t*)calloc(1, sizeof(*BR));
    BR->row_cluster = clt_BR; BR->col_cluster = clt_BR;
    BR->nrsons = 2; BR->ncsons = 2;
    BR->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*BR->sons));
    BR->dof_nrows = n2; BR->dof_ncols = n2;
    BR->dof_row_start = n1; BR->dof_col_start = n1;
    BUILD_LEAF_NU(BR->sons[0 + 0*2], clt_BR_0, clt_BR_0);
    BUILD_LEAF_NU(BR->sons[1 + 0*2], clt_BR_1, clt_BR_0);
    BUILD_LEAF_NU(BR->sons[0 + 1*2], clt_BR_0, clt_BR_1);
    BUILD_LEAF_NU(BR->sons[1 + 1*2], clt_BR_1, clt_BR_1);
    all_bn[n_bn++] = BR;

    /* BL leaf (n2 x n1) and TR leaf (n1 x n2). */
    st_cHACApK_block_node_t *BL = NULL, *TR = NULL;
    BUILD_LEAF_NU(BL, clt_BR, clt_TL);
    BUILD_LEAF_NU(TR, clt_TL, clt_BR);

    /* Root */
    st_cHACApK_block_node_t *root = (st_cHACApK_block_node_t*)calloc(1, sizeof(*root));
    root->row_cluster = clt_root; root->col_cluster = clt_root;
    root->nrsons = 2; root->ncsons = 2;
    root->sons = (st_cHACApK_block_node_t**)calloc(4, sizeof(*root->sons));
    root->dof_nrows = N; root->dof_ncols = N;
    root->dof_row_start = 0; root->dof_col_start = 0;
    root->sons[0 + 0*2] = TL;
    root->sons[1 + 0*2] = BL;
    root->sons[0 + 1*2] = TR;
    root->sons[1 + 1*2] = BR;
    all_bn[n_bn++] = root;

    int rc_dec = cHACApK_hlu_decomp(root);
    int rc_slv = (rc_dec == CHACAPK_HARITH_OK)
                  ? cHACApK_hlu_solve_vec(root, b, x_hlu, N)
                  : -99;

    double max_ref = 0.0, max_err = 0.0;
    for (int i = 0; i < N; i++) {
        double r = (x_ref[i] < 0.0) ? -x_ref[i] : x_ref[i];
        if (r > max_ref) max_ref = r;
        double d = x_hlu[i] - x_ref[i]; if (d < 0.0) d = -d;
        if (d > max_err) max_err = d;
    }
    double rel = (max_ref > 0.0) ? (max_err / max_ref) : max_err;

    for (int i = 0; i < n_lf; i++) { free(all_lf[i]->a1); free(all_lf[i]); }
    for (int i = 0; i < n_bn; i++) {
        if (all_bn[i]->sons) free(all_bn[i]->sons);
        free(all_bn[i]);
    }
    free(clt_root); free(clt_TL); free(clt_BR);
    free(clt_TL_0); free(clt_TL_1); free(clt_BR_0); free(clt_BR_1);
    free(A_full); free(A_ref); free(b); free(x_hlu); free(x_ref);

    #undef BUILD_LEAF_NU
    #undef RND_NU

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
            bn1->dof_nrows = 2*nb; bn1->dof_ncols = 2*nb;
            bn1->dof_row_start = clt_lvl1[i_root]->nstrt - 1;
            bn1->dof_col_start = clt_lvl1[j_root]->nstrt - 1;
            all_bn[n_bn++] = bn1;

            for (int j_sub = 0; j_sub < 2; j_sub++) {
                for (int i_sub = 0; i_sub < 2; i_sub++) {
                    int i_grid = 2*i_root + i_sub;
                    int j_grid = 2*j_root + j_sub;
                    st_cHACApK_block_node_t *bn = (st_cHACApK_block_node_t*)calloc(1, sizeof(*bn));
                    bn->row_cluster = clt_lvl2[i_grid];
                    bn->col_cluster = clt_lvl2[j_grid];
                    bn->dof_nrows = nb; bn->dof_ncols = nb;
                    bn->dof_row_start = clt_lvl2[i_grid]->nstrt - 1;
                    bn->dof_col_start = clt_lvl2[j_grid]->nstrt - 1;
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
    root->dof_nrows = N; root->dof_ncols = N;
    root->dof_row_start = 0; root->dof_col_start = 0;
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
