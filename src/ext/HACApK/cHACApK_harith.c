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

/* C -= A * B,  all dense column-major: A is mA x k, B is k x nB, C is mA x nB. */
static void dense_gemm_subtract(
    const double *A, int mA,
    const double *B, int kB, int nB,
    double *C, int ldC)
{
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                mA, nB, kB, -1.0, A, mA, B, kB, 1.0, C, ldC);
    g_stats.n_dense_gemm++;
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
        /* (1) factor diagonal block A_ii = L_ii U_ii */
        st_cHACApK_block_node_t *Aii = node->sons[i + i * s];
        int rc = hlu_rec(Aii);
        if (rc != CHACAPK_HARITH_OK) return rc;

        /* (2) for k > i: solve U_ik from L_ii U_ik = A_ik   (left-lower-unit) */
        for (int k = i + 1; k < s; k++) {
            st_cHACApK_block_node_t *Aik = node->sons[i + k * s];
            if (!leaf_is_dense(Aii) || !leaf_is_dense(Aik))
                return CHACAPK_HARITH_ERR_NEED_RECURSIVE; /* Phase 2: recursive trsm */
            dense_tri_solve_left_lower(
                leaf_dense_data(Aii), leaf_rows(Aii),
                leaf_dense_data(Aik), leaf_rows(Aik), leaf_cols(Aik));
        }
        /* (3) for j > i: solve L_ji from L_ji U_ii = A_ji   (right-upper-nonunit) */
        for (int j = i + 1; j < s; j++) {
            st_cHACApK_block_node_t *Aji = node->sons[j + i * s];
            if (!leaf_is_dense(Aii) || !leaf_is_dense(Aji))
                return CHACAPK_HARITH_ERR_NEED_RECURSIVE; /* Phase 2: recursive trsm */
            dense_tri_solve_right_upper(
                leaf_dense_data(Aii), leaf_rows(Aii),
                leaf_dense_data(Aji), leaf_rows(Aji), leaf_cols(Aji));
        }
        /* (4) trailing update: A_jk -= L_ji * U_ik  for j>i, k>i */
        for (int k = i + 1; k < s; k++) {
            for (int j = i + 1; j < s; j++) {
                st_cHACApK_block_node_t *Ajk = node->sons[j + k * s];
                st_cHACApK_block_node_t *Lji = node->sons[j + i * s];
                st_cHACApK_block_node_t *Uik = node->sons[i + k * s];
                if (!leaf_is_dense(Ajk) || !leaf_is_dense(Lji) || !leaf_is_dense(Uik))
                    return CHACAPK_HARITH_ERR_NEED_RECURSIVE; /* Phase 2: recursive gemm */
                dense_gemm_subtract(
                    leaf_dense_data(Lji), leaf_rows(Lji),
                    leaf_dense_data(Uik), leaf_rows(Uik), leaf_cols(Uik),
                    leaf_dense_data(Ajk), leaf_rows(Ajk));
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
     *   forward(sons[i+i*s], b)            // solve L_ii y_i = b_i
     *   for j > i:  b_j -= L_ji * y_i      // update remaining */
    int s = node->nrsons;
    for (int i = 0; i < s; i++) {
        st_cHACApK_block_node_t *Aii = node->sons[i + i * s];
        int rc = hlu_forward_rec(Aii, b);
        if (rc != CHACAPK_HARITH_OK) return rc;
        int row_i = Aii->row_cluster->nstrt - 1;
        int n_i = Aii->row_cluster->nsize;
        for (int j = i + 1; j < s; j++) {
            st_cHACApK_block_node_t *Lji = node->sons[j + i * s];
            if (!leaf_is_dense(Lji)) return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
            int row_j = Lji->row_cluster->nstrt - 1;
            int n_j = leaf_rows(Lji);
            /* b_j -= L_ji * y_i */
            cblas_dgemv(CblasColMajor, CblasNoTrans, n_j, n_i,
                        -1.0, leaf_dense_data(Lji), n_j,
                        &b[row_i], 1, 1.0, &b[row_j], 1);
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
    int s = node->nrsons;
    for (int i = s - 1; i >= 0; i--) {
        st_cHACApK_block_node_t *Aii = node->sons[i + i * s];
        int row_i = Aii->row_cluster->nstrt - 1;
        int n_i = Aii->row_cluster->nsize;
        /* update with upper off-diagonals: x_i -= U_ik * x_k, k > i */
        for (int k = i + 1; k < s; k++) {
            st_cHACApK_block_node_t *Uik = node->sons[i + k * s];
            if (!leaf_is_dense(Uik)) return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
            int col_k = Uik->col_cluster->nstrt - 1;
            int n_k = Uik->col_cluster->nsize;
            cblas_dgemv(CblasColMajor, CblasNoTrans, n_i, n_k,
                        -1.0, leaf_dense_data(Uik), n_i,
                        &x[col_k], 1, 1.0, &x[row_i], 1);
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
