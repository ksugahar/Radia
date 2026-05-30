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

/* Phase 3.5 unit test: h_addmul rk(A) * rk(B) -> rk(C) with recompression.
 * Builds 3 random rk leaves, computes the dense ground truth + alpha A B,
 * then calls h_addmul and verifies the result. Returns max element-wise
 * relative error. */
double cHACApK_harith_self_test_addmul_rkrk(int m, int n, int inner,
                                              int kA, int kB, int kC);

/* Phase 3.5 integration test: depth=2 H-LU with rk off-diagonal leaves.
 * Builds 4x4 leaf grid with dense diagonal + 12 rk off-diagonal leaves of
 * rank rk_rank, runs hlu_decomp + solve, compares to LAPACK dgesv.
 * Exercises ALL Phase 1-3.5 paths (the rk-rk -> rk via root-level trailing
 * update inside off-diagonal sub-blocks). */
double cHACApK_harith_self_test_rk_deep(int n_per_block, int rk_rank);

/* --------- Phase 4: storage-convention conversion + driver --------- *
 *
 * HACApK stores dense leaves ROW-MAJOR (a1[col + row*ndt] = M[row, col])
 * and rk leaves with a1=V, a2=U.  The H-LU code in this file expects
 * dense COLUMN-MAJOR (a1[row + col*ndl] = M[row, col]) and rk with
 * a1=U, a2=V.  Use the two helpers below to flip between conventions
 * IN PLACE on an existing leafmtxp.  Round-trip is identity.
 *
 * After convert_to_internal + cHACApK_hlu_decomp, the leaves contain
 * the in-place LU factors in internal format.  convert_to_hacapk WILL
 * transpose those factors back to HACApK row-major (so subsequent
 * HACApK MatVec would compute factor * x instead of A * x).  Caller's
 * choice. */
void cHACApK_convert_leafmtxp_to_internal(struct st_cHACApK_leafmtxp_t *lp);
void cHACApK_convert_leafmtxp_to_hacapk (struct st_cHACApK_leafmtxp_t *lp);

/* Phase 4 driver: run cHACApK_hlu_decomp + solve on a real HACApK tree.
 *
 * Inputs:
 *   leafmtxp_void    : st_cHACApK_leafmtxp_t* (in HACApK format, will be
 *                      converted to internal format and LU-factored).
 *                      The cluster-tree root must be preserved at
 *                      leafmtxp->st_clt_root (true after recent build
 *                      wrappers; see Phase 4 ground work commit).
 *   control_void     : st_cHACApK_lcontrol_t* (for lod permutation).
 *   x_orig, y_orig   : Length-N vectors in ORIGINAL (user) ordering.
 *                      Caller must precompute y_orig = A * x_orig via
 *                      HACApK_matvec_wrapper before calling this function.
 *
 * Returns max relative error |x_solved - x_orig| / max|x_orig| in the
 * original ordering. The leafmtxp is left in internal format with LU
 * factors -- use cHACApK_convert_leafmtxp_to_hacapk to restore if
 * needed, but the LU OVERWROTE the original A entries so MatVec on
 * the restored leafmtxp would be meaningless. Re-build the H-matrix
 * to recover A.
 *
 * On internal error (allocation, block-tree build, LU failure), returns
 * a negative sentinel:
 *   -1.0     : null pointer
 *   -2.0     : alloc failure
 *   -3.0     : block-tree build failure
 *   -4.0 + rc * 0.001 : hlu_decomp returned rc (e.g., -4.005 = NEED_RECURSIVE)
 *   -5.0 + rc * 0.001 : hlu_solve_vec returned rc */
double cHACApK_hlu_run_on_hacapk(void *leafmtxp_void, void *control_void,
                                  const double *x_orig, const double *y_orig,
                                  int nffc);

#ifdef __cplusplus
}
#endif

#endif /* CHACAPK_HARITH_H_INCLUDED */
