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


/* C += alpha * A * B  (block-recursive, dense-leaf base case). */
static int h_addmul(double alpha,
                 const st_cHACApK_block_node_t *A,
                 const st_cHACApK_block_node_t *B,
                 st_cHACApK_block_node_t *C)
{
    if (!A || !B || !C) return CHACAPK_HARITH_ERR_NULL;

    if (leaf_is_rk(A) || leaf_is_rk(B) || leaf_is_rk(C)) {
        g_stats.n_lowrank_skip++;
        return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
    }

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

    /* Mixed leaf+internal: reserved for Phase 3 (densify the internal one). */
    if (leaf_is_dense(A) || leaf_is_dense(B) || leaf_is_dense(C)) {
        return CHACAPK_HARITH_ERR_NEED_RECURSIVE;
    }

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

    if (leaf_is_rk(L) || leaf_is_rk(X)) {
        g_stats.n_lowrank_skip++;
        return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
    }

    if (leaf_is_dense(L) && leaf_is_dense(X)) {
        dense_tri_solve_left_lower(
            leaf_dense_data(L), leaf_rows(L),
            leaf_dense_data(X), leaf_rows(X), leaf_cols(X));
        return CHACAPK_HARITH_OK;
    }

    if (leaf_is_dense(L) || leaf_is_dense(X)) {
        return CHACAPK_HARITH_ERR_NEED_RECURSIVE;
    }

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

    if (leaf_is_rk(U) || leaf_is_rk(X)) {
        g_stats.n_lowrank_skip++;
        return CHACAPK_HARITH_ERR_LOWRANK_LEAF;
    }

    if (leaf_is_dense(U) && leaf_is_dense(X)) {
        dense_tri_solve_right_upper(
            leaf_dense_data(U), leaf_rows(U),
            leaf_dense_data(X), leaf_rows(X), leaf_cols(X));
        return CHACAPK_HARITH_OK;
    }

    if (leaf_is_dense(U) || leaf_is_dense(X)) {
        return CHACAPK_HARITH_ERR_NEED_RECURSIVE;
    }

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
    if (leaf_is_rk(node)) return CHACAPK_HARITH_ERR_LOWRANK_LEAF;

    if (leaf_is_dense(node)) {
        int m = leaf_rows(node), n = leaf_cols(node);
        int r0 = node->row_cluster->nstrt - 1;  /* 1-based -> 0-based */
        int c0 = node->col_cluster->nstrt - 1;
        cblas_dgemv(CblasColMajor, CblasNoTrans, m, n,
                    -1.0, leaf_dense_data(node), m,
                    &x[c0], 1, 1.0, &y[r0], 1);
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
