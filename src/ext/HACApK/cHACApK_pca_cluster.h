/* cHACApK_pca_cluster.h
 *
 * PCA-based cluster split for HACApK. See cHACApK_pca_cluster.c for the
 * algorithm; this header exposes the strategy switch + the split helper.
 *
 * The split helper is intended to be called from inside
 * cHACApK_generate_cbitree as a drop-in alternative to the current
 * "midpoint of longest bbox axis" split, gated by the strategy flag.
 */

#ifndef CHACAPK_PCA_CLUSTER_H_INCLUDED
#define CHACAPK_PCA_CLUSTER_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

/* Cluster strategy enum. Currently:
 *   BBOX  : axis-aligned bbox split at midpoint of longest extent (default,
 *           the historical HACApK behaviour, what all existing tests rely on)
 *   PCA   : split perpendicular to the dominant covariance eigenvector at
 *           the center of mass (handles flat / elongated geometries cleanly) */
enum cHACApK_cluster_strategy {
  CHACAPK_CLUSTER_BBOX = 0,
  CHACAPK_CLUSTER_PCA  = 1
};

/* Setter / getter for the strategy. Global because the cluster build is
 * already a process-global operation in HACApK; matching scope keeps the
 * API simple. Default: BBOX (backward compatible). */
void cHACApK_set_cluster_strategy(int strategy);
int  cHACApK_get_cluster_strategy(void);

/* PCA split (see .c file for full semantics).
 * Returns:
 *   n_left in [1, cur_nd-1]  : number of indices on the left side; left =
 *                              cur_lod[1..n_left], right = cur_lod[n_left+1..cur_nd]
 *   -1                       : failure -> caller MUST fall back to bbox split */
int cHACApK_pca_split(
    double **zgmid_t,
    int *cur_lod,
    int cur_nd,
    int ndim);

#ifdef __cplusplus
}
#endif

#endif /* CHACAPK_PCA_CLUSTER_H_INCLUDED */
