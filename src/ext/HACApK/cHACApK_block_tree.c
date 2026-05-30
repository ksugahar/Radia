/* cHACApK_block_tree.c -- thin block-tree view on HACApK's flat leaves.
 *
 * Builder strategy:
 *   1. Index the flat st_lf[] array by (nstrtl, nstrtt) -> leaf pointer
 *      via a sorted lookup (small log-N cost per query, well-behaved on
 *      the lattice-like H-matrix HACApK produces).
 *   2. Recursively walk the cluster-pair tree from (root_rc, root_cc).
 *      At each visited (rc, cc):
 *        - probe the index for a leaf at (rc->nstrt, cc->nstrt) with
 *          matching size (rc->nsize, cc->nsize)
 *        - if hit -> leaf block-tree node (cache ltmtx as leaf_kind)
 *        - if miss -> internal node, recurse on rc->pc_sons x cc->pc_sons
 *          (HACApK indexes children 1..nnson, the [0] slot is unused).
 *
 * Memory: each node ~ 56 bytes + sons-array. For a balanced binary
 * cluster tree of N elements with leaf_size=64, total nodes ~ 2*N/64
 * = N/32. At N=165000 (C-type 20^3) that's ~5k nodes ~ 280 KB. Fine.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cHACApK_block_tree.h"

/* ---------------------------------------------------------------- *
 * Internal helpers                                                  *
 * ---------------------------------------------------------------- */

/* qsort comparator: order leaves by (nstrtl, nstrtt). */
static int leaf_key_cmp(const void *aa, const void *bb)
{
    const st_cHACApK_leafmtx_t *a = *(const st_cHACApK_leafmtx_t * const *)aa;
    const st_cHACApK_leafmtx_t *b = *(const st_cHACApK_leafmtx_t * const *)bb;
    if (a->nstrtl != b->nstrtl) return (a->nstrtl < b->nstrtl) ? -1 : 1;
    if (a->nstrtt != b->nstrtt) return (a->nstrtt < b->nstrtt) ? -1 : 1;
    return 0;
}

/* Binary search the sorted leaf-pointer array for a leaf with exact
 * (nstrtl, nstrtt) match. Returns NULL on miss. Caller verifies sizes. */
static st_cHACApK_leafmtx_t *find_leaf_sorted(
    st_cHACApK_leafmtx_t **sorted, int n_sorted,
    int nstrtl, int nstrtt)
{
    int lo = 0, hi = n_sorted - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        const st_cHACApK_leafmtx_t *m = sorted[mid];
        if (m->nstrtl < nstrtl || (m->nstrtl == nstrtl && m->nstrtt < nstrtt)) {
            lo = mid + 1;
        } else if (m->nstrtl > nstrtl || (m->nstrtl == nstrtl && m->nstrtt > nstrtt)) {
            hi = mid - 1;
        } else {
            return sorted[mid];
        }
    }
    return NULL;
}

/* Allocate + zero-init a node. */
static st_cHACApK_block_node_t *new_node(
    st_cHACApK_cluster_t *rc, st_cHACApK_cluster_t *cc)
{
    st_cHACApK_block_node_t *n = (st_cHACApK_block_node_t *)
        calloc(1, sizeof(st_cHACApK_block_node_t));
    if (!n) return NULL;
    n->row_cluster = rc;
    n->col_cluster = cc;
    return n;
}

/* Recursive builder.  nffc = uniform DOF per element (or 1 if cluster
 * tree already in DOF units, as in the synthetic tests). */
static st_cHACApK_block_node_t *build_rec(
    st_cHACApK_leafmtx_t **sorted, int n_sorted,
    st_cHACApK_cluster_t *rc, st_cHACApK_cluster_t *cc,
    int nffc)
{
    if (!rc || !cc) return NULL;
    st_cHACApK_block_node_t *node = new_node(rc, cc);
    if (!node) return NULL;

    /* Convert element-grain cluster positions to DOF grain via uniform
     * nffc. nffc == 1 makes this a no-op. */
    int rc_dof_start = (rc->nstrt - 1) * nffc + 1;  /* 1-based */
    int cc_dof_start = (cc->nstrt - 1) * nffc + 1;
    int rc_dof_size  = rc->nsize * nffc;
    int cc_dof_size  = cc->nsize * nffc;

    /* Probe for a leaf at this exact DOF-position (rc, cc). */
    st_cHACApK_leafmtx_t *leaf =
        find_leaf_sorted(sorted, n_sorted, rc_dof_start, cc_dof_start);
    if (leaf && leaf->ndl == rc_dof_size && leaf->ndt == cc_dof_size) {
        node->leaf_mtx  = leaf;
        node->leaf_kind = leaf->ltmtx;   /* 1 = rk, 2 = dense */
        return node;
    }

    /* Internal node -- must have children to descend. */
    if (rc->nnson <= 0 || cc->nnson <= 0) {
        /* No children but no matching leaf either -> topology mismatch. */
        free(node);
        fprintf(stderr,
                "[cHACApK_block_tree] topology mismatch at "
                "elem (nstrt=%d size=%d) x (nstrt=%d size=%d) "
                "= dof (nstrt=%d size=%d) x (nstrt=%d size=%d): "
                "no leaf, no children, nffc=%d. Aborting build.\n",
                rc->nstrt, rc->nsize, cc->nstrt, cc->nsize,
                rc_dof_start, rc_dof_size, cc_dof_start, cc_dof_size,
                nffc);
        return NULL;
    }

    node->nrsons = rc->nnson;
    node->ncsons = cc->nnson;
    node->sons = (st_cHACApK_block_node_t **)
        calloc((size_t)node->nrsons * (size_t)node->ncsons,
               sizeof(st_cHACApK_block_node_t *));
    if (!node->sons) { free(node); return NULL; }

    for (int j = 0; j < node->ncsons; j++) {
        /* HACApK uses 1-based child indexing: pc_sons[1..nnson]. */
        st_cHACApK_cluster_t *cc_son = cc->pc_sons[j + 1];
        for (int i = 0; i < node->nrsons; i++) {
            st_cHACApK_cluster_t *rc_son = rc->pc_sons[i + 1];
            st_cHACApK_block_node_t *child =
                build_rec(sorted, n_sorted, rc_son, cc_son, nffc);
            if (!child) {
                for (int k = 0; k < i + j * node->nrsons; k++) {
                    cHACApK_free_block_tree(node->sons[k]);
                }
                free(node->sons);
                free(node);
                return NULL;
            }
            node->sons[i + j * node->nrsons] = child;
        }
    }
    return node;
}

/* ---------------------------------------------------------------- *
 * Public API                                                        *
 * ---------------------------------------------------------------- */

st_cHACApK_block_node_t *cHACApK_build_block_tree(
    st_cHACApK_leafmtxp_t *leafmtxp,
    st_cHACApK_cluster_t  *root_row_cluster,
    st_cHACApK_cluster_t  *root_col_cluster)
{
    return cHACApK_build_block_tree_nffc(leafmtxp, root_row_cluster,
                                          root_col_cluster, 1);
}

st_cHACApK_block_node_t *cHACApK_build_block_tree_nffc(
    st_cHACApK_leafmtxp_t *leafmtxp,
    st_cHACApK_cluster_t  *root_row_cluster,
    st_cHACApK_cluster_t  *root_col_cluster,
    int                    nffc)
{
    if (!leafmtxp || !root_row_cluster || !root_col_cluster) return NULL;
    if (leafmtxp->nlf <= 0 || !leafmtxp->st_lf) return NULL;
    if (nffc <= 0) return NULL;

    /* Build a sorted index of leaf pointers for O(log n) lookup. HACApK
     * st_lf[] is 1-based: slot 0 unused, slots 1..nlf hold the leaves. */
    st_cHACApK_leafmtx_t **sorted = (st_cHACApK_leafmtx_t **)
        malloc(sizeof(st_cHACApK_leafmtx_t *) * (size_t)leafmtxp->nlf);
    if (!sorted) return NULL;
    int n_sorted = 0;
    for (int i = 1; i <= leafmtxp->nlf; i++) {
        if (leafmtxp->st_lf[i]) sorted[n_sorted++] = leafmtxp->st_lf[i];
    }
    qsort(sorted, (size_t)n_sorted, sizeof(st_cHACApK_leafmtx_t *), leaf_key_cmp);

    st_cHACApK_block_node_t *root =
        build_rec(sorted, n_sorted, root_row_cluster, root_col_cluster, nffc);

    free(sorted);
    return root;
}

void cHACApK_free_block_tree(st_cHACApK_block_node_t *root)
{
    if (!root) return;
    if (root->sons) {
        int nsons = root->nrsons * root->ncsons;
        for (int i = 0; i < nsons; i++) cHACApK_free_block_tree(root->sons[i]);
        free(root->sons);
    }
    free(root);
}

static void block_tree_stats_rec(
    const st_cHACApK_block_node_t *n, int depth,
    int *nrk, int *ndense, int *nint, int *maxd)
{
    if (!n) return;
    if (depth > *maxd) *maxd = depth;
    if (n->leaf_mtx) {
        if (n->leaf_kind == 1) (*nrk)++;
        else (*ndense)++;
    } else {
        (*nint)++;
        int nsons = n->nrsons * n->ncsons;
        for (int i = 0; i < nsons; i++)
            block_tree_stats_rec(n->sons[i], depth + 1, nrk, ndense, nint, maxd);
    }
}

void cHACApK_block_tree_stats(
    const st_cHACApK_block_node_t *root,
    int *out_n_leaves_rk,
    int *out_n_leaves_dense,
    int *out_n_internal,
    int *out_max_depth)
{
    int nrk = 0, nd = 0, ni = 0, md = 0;
    block_tree_stats_rec(root, 0, &nrk, &nd, &ni, &md);
    if (out_n_leaves_rk)    *out_n_leaves_rk    = nrk;
    if (out_n_leaves_dense) *out_n_leaves_dense = nd;
    if (out_n_internal)     *out_n_internal     = ni;
    if (out_max_depth)      *out_max_depth      = md;
}
