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
  double t_solve_sec;    /* wall time of the most recent solve  */
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

/* Phase 4 debug test: mimic Radia's mixed-sibling tree structure.
 *
 * Builds a depth=3 tree where the root has 2x2 children but two of them
 * are internal (TL and BR, each containing 2x2 leaves) and two are
 * leaves (TR and BL). This matches the structure HACApK produces for
 * Radia magnetostatic H-matrix blocks at small leaf_size (e.g., nx=3 hex cube with
 * leaf_size=10 elements: 10 leaves, 3 internal nodes, depth 3).
 *
 * Both root-level children are square sub-blocks (2*nb_small x 2*nb_small),
 * so the LU partition is well-defined. All leaves are dense, no rk.
 *
 * Returns max relative error vs LAPACKE_dgesv on the same matrix.
 * If this test PASSES while real Radia trees fail, the bug is in
 * non-uniform leaf sizes (HACApK splits elements asymmetrically:
 * 13 -> 6+7). If this test FAILS, the bug is in the mixed-sibling
 * recursion itself. */
double cHACApK_harith_self_test_mixed_sibling(int nb_small);

/* Phase 4 debug: non-uniform-size variant of mixed_sibling.
 * Root splits into (n1, n2) where n1 != n2 (asymmetric, mimics HACApK's
 * element-count splits like 13 -> 6+7). TL = n1 x n1 internal split into
 * (m1, n1-m1) sub-leaves. BR = n2 x n2 internal split into (m3, n2-m3).
 * TR / BL leaves are rectangular (n1 x n2 / n2 x n1).
 *
 * If uniform mixed_sibling passes but this fails, the bug is in
 * non-uniform leaf sub-views (Phase 3.6 mixed materialize/distribute
 * with sibling leaves of different sizes). */
double cHACApK_harith_self_test_mixed_sibling_nonuniform(
    int n1, int n2, int m1, int m3);

/* Phase 4 debug: depth=3 asymmetric tree (mimics Radia's nx=3 leaf=10
 * tree shape exactly: 10 leaves, 3 internal nodes, max depth 3).
 *
 *   Root (2x2 internal)
 *     TL (2x2 internal)
 *       TL.TL (2x2 internal) -- 4 small leaves at depth 3
 *       TL.TR, TL.BL, TL.BR -- 3 leaves at depth 2
 *     TR, BL, BR -- 3 leaves at depth 1
 *
 * If mixed_sibling (depth 2) passes but this fails, the bug is in the
 * deeper recursion with mixed leaf+internal at multiple levels. */
double cHACApK_harith_self_test_depth3_asymmetric(int nb_tiny);

/* Phase 4 debug: EXACT mimic of Radia nx=3 leaf=10 tree (sizes 108/54
 * at root, 72/36 at TL, 48/24 at TL.TL). Hardcoded. */
double cHACApK_harith_self_test_radia_exact(void);

/* Phase 4 debug: same as radia_exact but with adjustable diag_boost
 * (default 2.0 = mildly diagonally dominant; lower values approach weakly
 * matrix's weaker dominance to test no-pivot LU stability). */
double cHACApK_harith_self_test_radia_exact_diag(double diag_boost);

/* Phase 4 debug: same Radia-exact tree shape but with EXTERNAL matrix.
 * A_full is 162x162 column-major; b is the RHS (caller provides).
 * Returns max rel err vs LAPACKE_dgesv on the same matrix. */
double cHACApK_harith_self_test_radia_exact_with_matrix(
    const double *A_full, const double *b);

/* Phase 4 debug: mixed_sibling test with HACApK row-major leaves, then
 * convert to internal format before H-LU.  Mimics the EXACT path used
 * by cHACApK_hlu_run_on_hacapk on real Radia trees, but with synthetic
 * uniform sizes.  If this fails while plain mixed_sibling passes,
 * the bug is in cHACApK_convert_leafmtxp_to_internal or its
 * interaction with H-LU. */
double cHACApK_harith_self_test_mixed_sibling_via_conversion(int nb_small);


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

/* H-LU as a reusable PRECONDITIONER: factor a leafmtxp once (convert + build
 * block-tree + hlu_decomp -> opaque block-tree root), apply the solve many
 * times (permute via control->lod + hlu_solve_vec + un-permute; r,z in ORIGINAL
 * ordering), then free. Used to H-LU-precondition the A_SS = S^T A S star block. */
void* cHACApK_hlu_factor_leafmtxp(void* leafmtxp_void, void* control_void, int nffc);
int   cHACApK_hlu_apply(void* root_void, void* control_void,
                        const double* r, double* z, int nd);
void  cHACApK_hlu_free_factors(void* root_void);

/* Phase 4 debug: materialize the post-convert tree as a dense matrix.
 * Returns the matrix in PERMUTED ordering (caller can apply lod to
 * compare with HMatrixDensify which gives original ordering).
 *
 * Inputs:
 *   leafmtxp_void, control_void: HACApK structures (will be converted)
 *   nffc: uniform DOF per element
 *
 * Outputs:
 *   A_perm_out: caller-allocated [nd * nd] double buffer (column-major,
 *               stores A_perm = P A_orig P^T in permuted ordering)
 *   lod_out:    caller-allocated [nd] int buffer (0-based, lod_out[i] =
 *               original index of permuted position i)
 *   nd_out:     written with the DOF count
 *
 * Returns 0 on success, negative on error.
 * Note: the leafmtxp is left in INTERNAL layout after this call.
 *       Cluster tree on leafmtxp is preserved. */
int cHACApK_hlu_debug_materialize(void *leafmtxp_void, void *control_void,
                                    int nffc,
                                    double *A_perm_out,
                                    int *lod_out,
                                    int *nd_out);

#ifdef __cplusplus
}
#endif

#endif /* CHACAPK_HARITH_H_INCLUDED */
