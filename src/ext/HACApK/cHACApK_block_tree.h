/* cHACApK_block_tree.h
 *
 * Block-tree VIEW on top of HACApK's flat leaf array.
 *
 * Background. HACApK stores the H-matrix as a flat array of leaves
 * (st_leafmtxp.st_lf[]), each annotated with its (row_start, row_size,
 * col_start, col_size) position in the global matrix. The HACApK MatVec
 * iterates over this flat array.
 *
 * H-LU (port from H2Lib::lrdecomp_hmatrix and triangularinvmul_hmatrix)
 * needs a RECURSIVE block-tree: each node either (a) is a leaf pointing
 * at a HACApK leafmtx, or (b) has 4 children (2x2 block split) whose
 * row/col clusters are the children of this node's row/col cluster.
 *
 * The block-tree built here is a READ-ONLY VIEW. It owns no leaf data;
 * it only owns the tree topology (struct st_cHACApK_block_node nodes
 * and the sons-array pointers). All leaf data continues to live in the
 * HACApK leaf array, untouched. HACApK MatVec is unaffected.
 *
 * Lifecycle. Build once with cHACApK_build_block_tree (input:
 * st_leafmtxp + root cluster), free with cHACApK_free_block_tree.
 * The cluster tree (st_clt) MUST outlive the block-tree.
 *
 * Phase 0.7 of the H-LU project (Tasks #30 onward). See
 * docs/research/h-lu/H2LIB_TO_HACAPK_MAPPING.ipynb.
 */

#ifndef CHACAPK_BLOCK_TREE_H_INCLUDED
#define CHACAPK_BLOCK_TREE_H_INCLUDED

#include "cHACApK_base.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Node of the block-tree.
 *
 * Either:
 *   - LEAF: leaf_mtx != NULL, sons == NULL. nrsons = ncsons = 0.
 *     Points at one of st_leafmtxp.st_lf[].
 *   - INTERNAL: leaf_mtx == NULL, sons != NULL.
 *     nrsons = number of row clusters' children, ncsons = column's.
 *     sons[i + j*nrsons] is the child for (row_cluster->pc_sons[i+1],
 *     col_cluster->pc_sons[j+1])  (HACApK indexes children 1-based).
 *
 * row_cluster, col_cluster point into the existing cluster tree; the
 * block-tree does NOT own them.
 *
 * For convenient algorithmic dispatch in the H-LU code, we also cache
 * leaf_mtx->ltmtx (1 = low-rank, 2 = dense) as `leaf_kind` (0 for
 * non-leaf nodes). This avoids a NULL pointer chase in inner loops.
 */
typedef struct st_cHACApK_block_node_t {
  st_cHACApK_cluster_t *row_cluster;
  st_cHACApK_cluster_t *col_cluster;
  st_cHACApK_leafmtx_t *leaf_mtx;       /* NULL if internal node */
  int                   leaf_kind;       /* 0 internal, 1 rk leaf, 2 dense leaf */
  int                   nrsons;          /* 0 if leaf */
  int                   ncsons;          /* 0 if leaf */
  struct st_cHACApK_block_node_t **sons; /* size nrsons*ncsons, col-major */
  /* Phase 3.6 / Phase 4: DOF-grain dimensions stored alongside the
   * element-grain cluster pointers. For synthetic tests (cHACApK_build_
   * block_tree, nffc=1) these equal cluster->nsize / cluster->nstrt-1.
   * For real HACApK trees (cHACApK_build_block_tree_nffc, nffc > 1)
   * these are cluster->nsize * nffc / (cluster->nstrt - 1) * nffc.
   * Used by materialize/add/set helpers and hmatvec_subtract so they
   * don't have to thread nffc through every recursion. */
  int                   dof_nrows;       /* DOF count of row_cluster */
  int                   dof_ncols;       /* DOF count of col_cluster */
  int                   dof_row_start;   /* 0-based DOF start of row */
  int                   dof_col_start;   /* 0-based DOF start of col */
} st_cHACApK_block_node_t;

/* Build the block-tree view for the H-matrix described by leafmtxp,
 * with root_row_cluster and root_col_cluster as the row/column index
 * spans (typically both equal to the root cluster for a square matrix).
 *
 * Memory: O(#nodes) = O(#leaves + #internal-nodes).
 *
 * Returns NULL on allocation failure or topology mismatch (the latter
 * indicates a HACApK build that does not align with the cluster tree --
 * very rare; means the leaf array partition is incomplete). */
st_cHACApK_block_node_t *cHACApK_build_block_tree(
    st_cHACApK_leafmtxp_t *leafmtxp,
    st_cHACApK_cluster_t  *root_row_cluster,
    st_cHACApK_cluster_t  *root_col_cluster);

/* Phase 4 variant: cluster ranges in ELEMENT units, leaf positions in
 * DOF units. nffc = uniform DOF per element (3 tet / 6 hex MSC / etc.).
 *
 * HACApK's cHACApK_generate_cbitree builds clusters in element units
 * (nstrt, nsize count elements). The leafmtxp leaves are filled with
 * nstrtl/ndl in DOF units (cumulative DOF position of permuted
 * elements). This variant scales the cluster ranges by nffc to bridge
 * the two grains.
 *
 * Restriction: requires UNIFORM DOF per element. For mixed-element
 * varDOF cases (tet+hex+wedge), use the cumdof-aware variant (TODO).
 *
 * Pass nffc = 1 to recover the cHACApK_build_block_tree behavior. */
st_cHACApK_block_node_t *cHACApK_build_block_tree_nffc(
    st_cHACApK_leafmtxp_t *leafmtxp,
    st_cHACApK_cluster_t  *root_row_cluster,
    st_cHACApK_cluster_t  *root_col_cluster,
    int                    nffc);

/* Recursively free the block-tree (the leaves' underlying st_lf[] data
 * is NOT freed -- it belongs to HACApK). */
void cHACApK_free_block_tree(st_cHACApK_block_node_t *root);

/* Diagnostic: count leaves + internal nodes. */
void cHACApK_block_tree_stats(
    const st_cHACApK_block_node_t *root,
    int *out_n_leaves_rk,
    int *out_n_leaves_dense,
    int *out_n_internal,
    int *out_max_depth);

#ifdef __cplusplus
}
#endif

#endif /* CHACAPK_BLOCK_TREE_H_INCLUDED */
