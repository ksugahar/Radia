/* cHACApK_harith.h
 *
 * H-matrix arithmetic for HACApK -- port of H2Lib's harith.c on top of
 * the block-tree view (cHACApK_block_tree.h).
 *
 * This is the operational core of H-LU: H-triangular solve + H-matrix
 * multiply-add. Both functions are block-recursive over the block tree
 * built by cHACApK_build_block_tree.
 *
 * SCOPE OF THIS REVISION (Phase 1/2 minimal):
 *   - dense leaf x dense leaf  (LAPACK dtrsm + dgemm)
 *   - the recursive block driver (no low-rank leaves yet)
 *
 * NOT YET SUPPORTED (will be added in Phase 2 proper):
 *   - low-rank leaves (rkmatrix UV^T)
 *   - rank growth recompression via ACA / SVD truncation
 *
 * Existing HACApK builds may have low-rank leaves; functions in this
 * revision detect them and return a non-zero error code so the caller
 * can decide to fall back.
 */

#ifndef CHACAPK_HARITH_H_INCLUDED
#define CHACAPK_HARITH_H_INCLUDED

#include "cHACApK_base.h"
#include "cHACApK_block_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Error codes. */
enum {
  CHACAPK_HARITH_OK = 0,
  CHACAPK_HARITH_ERR_LOWRANK_LEAF = -1,  /* low-rank leaf not yet handled (Phase 2) */
  CHACAPK_HARITH_ERR_LAPACK       = -2,  /* LAPACK info != 0 */
  CHACAPK_HARITH_ERR_TOPOLOGY     = -3,  /* block-tree sons-count mismatch */
  CHACAPK_HARITH_ERR_NULL         = -4,
  CHACAPK_HARITH_ERR_NEED_RECURSIVE = -5 /* off-diag op needs recursive H-arith (Phase 2) */
};

/* H-LU FACTORIZATION (in-place).
 *
 * Computes the LU factorization A = L U on the H-matrix represented by
 * root. The factors are stored in-place over the same leaf data
 * (lower-triangle of dense diagonal blocks holds L, upper holds U).
 * Off-diagonal dense blocks store L_ij (below diag) or U_ij (above).
 *
 * Returns CHACAPK_HARITH_OK on success. The dense-leaf LAPACK uses
 * LAPACKE_dgetrf WITHOUT pivoting (we ignore the pivot vector in this
 * revision; high-mu_r non-symmetric matrices may need pivoting in
 * Phase 3 -- the no-pivoting choice matches H2Lib's lrdecomp_hmatrix
 * which also assumes a stable factorization). */
int cHACApK_hlu_decomp(st_cHACApK_block_node_t *root);

/* Apply the H-LU factorization to a dense vector: solve A x = b where
 * the factors are in `root` (after cHACApK_hlu_decomp).
 *
 * b and x are dense vectors of length n (must match the cluster's
 * nsize). For an in-place solve, pass x == b.
 *
 * Returns CHACAPK_HARITH_OK on success. */
int cHACApK_hlu_solve_vec(
    const st_cHACApK_block_node_t *root,
    const double *b,
    double *x,
    int n);

/* Diagnostic: count operations during the most recent decomp. */
typedef struct cHACApK_hlu_stats_t {
  long n_dense_lu;       /* leaf dgetrf calls                */
  long n_dense_trsm;     /* dtrsm calls on off-diag dense    */
  long n_dense_gemm;     /* dgemm calls in trailing update   */
  long n_lowrank_skip;   /* low-rank leaves we skipped       */
  double t_decomp_sec;   /* wall time of the most recent decomp */
} cHACApK_hlu_stats_t;

const cHACApK_hlu_stats_t *cHACApK_hlu_last_stats(void);

/* Self-test: construct a synthetic dense matrix recursively split into a
 * (2^depth) x (2^depth) grid of dense leaves, run cHACApK_hlu_decomp +
 * cHACApK_hlu_solve_vec, compare against a reference LAPACKE_dgesv solve.
 * Returns the maximum componentwise error |x_hlu - x_ref| / max|x_ref|.
 * Should be < 1e-10 for a well-conditioned diag-dominant matrix.
 *
 * depth = 1: 2x2 block-tree (4 leaves) -- shallow recursion sanity
 * depth = 2: root has 4 sons, each with 4 sons (16 leaves) -- recursion
 * depth = 3: 64 leaves (~ realistic HACApK leaf count at moderate scale) */
double cHACApK_harith_self_test(int depth, int n_per_block);

/* Rk-aware self-test (Phase 3 partial validation, depth=1 only).
 *
 * Builds a 2x2 block-tree with DENSE diagonal leaves (random
 * diagonally-dominant) and explicit-rank RK off-diagonal leaves
 * (constructed as U_ij V_ij^T from random U_ij, V_ij of rank rk_rank).
 *
 * Exercises the Phase 3 partial paths:
 *   - htrsm_lln(L=dense, X=rk),  htrsm_run(U=dense, X=rk)
 *   - h_addmul(rk*rk -> dense) on the trailing update
 *   - hmatvec_subtract on rk leaves
 *
 * Depth >= 2 with rk off-diagonals would also need rk(A)*rk(B) -> rk(C)
 * in trailing updates, which requires ACA recompression (Phase 3.5).
 *
 * Returns max relative error vs LAPACKE_dgesv (should be ~ machine
 * precision for a well-conditioned diag-dominant matrix). */
double cHACApK_harith_self_test_rk(int n_per_block, int rk_rank);

#ifdef __cplusplus
}
#endif

#endif /* CHACAPK_HARITH_H_INCLUDED */
