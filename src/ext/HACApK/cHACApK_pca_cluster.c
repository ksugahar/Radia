/* cHACApK_pca_cluster.c
 *
 * Principal-component analysis (PCA) cluster split helper for HACApK.
 *
 * Replaces the default axis-aligned "split at midpoint of longest bbox axis"
 * with "split at center of mass perpendicular to the dominant covariance
 * eigenvector". This recovers tight, compact bounding boxes on flat /
 * elongated / curved meshes (thin plates, shielding sheets, beam-pipe liners,
 * flat coils), where axis-aligned bbox splits over-allocate bbox volume and
 * the H-matrix admissibility condition (eta * diam(B1) * diam(B2) > dist(B1,B2))
 * fails too often, degrading ACA compression and inflating storage.
 *
 * The covariance eigendecomposition uses LAPACKE_dsyev (MKL via mkl.h, already
 * linked by Radia/HACApK). For typical dim=3 meshes this is a 3x3 dsyev call,
 * negligible cost (~us per split). Falls back to the original bbox split if
 * dsyev fails or the covariance is rank-deficient (e.g. all points collinear).
 *
 * Algorithm (Hackbusch, "Hierarchical Matrices", Algorithm 1.16):
 *   1. compute center of mass of the cluster's points
 *   2. assemble the 3x3 covariance matrix C = sum_i (x_i - cm)(x_i - cm)^T
 *   3. eigendecompose: C = Q diag(lambda) Q^T  (lambda ascending)
 *   4. split normal v = eigenvector of LARGEST lambda (= Q[:, dim-1])
 *   5. partition indices by sign of (x_i - cm) . v
 *
 * This file is intentionally self-contained: the only dependencies are
 * <math.h>, <stdlib.h>, <mkl.h> (for LAPACKE_dsyev), and the cHACApK basic
 * types via cHACApK_base.h. The dispatch from cHACApK_generate_cbitree is
 * added as a separate edit so this file can be compiled and unit-tested
 * standalone.
 *
 * Reference: H2Lib/Library/cluster.c::build_pca_cluster (LGPL).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl.h>

#include "cHACApK_base.h"
#include "cHACApK_pca_cluster.h"

/* Static global flag controlling whether PCA split is used. Default OFF for
 * backward-compatible behaviour. Switched on via cHACApK_set_cluster_strategy.
 * Static-global avoids touching the public struct layouts for now (Option A
 * port plan; see docs/research/h-lu/H2LIB_TO_HACAPK_MAPPING.ipynb). */
static int g_cluster_strategy = CHACAPK_CLUSTER_BBOX;

void cHACApK_set_cluster_strategy(int strategy)
{
  g_cluster_strategy = strategy;
}

int cHACApK_get_cluster_strategy(void)
{
  return g_cluster_strategy;
}

/* PCA-based split of a cluster.
 *
 * Inputs:
 *   zgmid_t  [ndim+1][nofc+1] : geometry midpoints, 1-based indexing as the
 *                                rest of HACApK uses
 *   cur_lod  [cur_nd+1]       : 1-based index list into zgmid_t; this routine
 *                                permutes it in place so the first <retval>
 *                                entries are "left" of the split plane and
 *                                the remaining are "right"
 *   cur_nd                    : number of indices in cur_lod
 *   ndim                      : geometric dimension (1, 2, or 3)
 *
 * Returns:
 *   nl in [1, cur_nd-1] (1-based)   : number of indices on the "left" side
 *                                      of the split plane (= signed projection
 *                                      onto principal axis <= 0)
 *   -1 on failure (covariance singular, dsyev info != 0, or split degenerate
 *                   = all points on the same side); caller MUST fall back to
 *                   the bbox split path.
 *
 * Side effect: cur_lod is permuted (in place) so cur_lod[1..nl] are "left"
 * and cur_lod[nl+1..cur_nd] are "right".
 */
int cHACApK_pca_split(
    double **zgmid_t,
    int *cur_lod,
    int cur_nd,
    int ndim)
{
  if (cur_nd < 2 || ndim < 1 || ndim > 3) return -1;

  /* 1. center of mass (use 0-based local arrays for clarity) */
  double cm[3] = {0.0, 0.0, 0.0};
  for (int i = 1; i <= cur_nd; i++) {
    int idx = cur_lod[i];
    for (int j = 0; j < ndim; j++) cm[j] += zgmid_t[j+1][idx];
  }
  double inv_n = 1.0 / (double)cur_nd;
  for (int j = 0; j < ndim; j++) cm[j] *= inv_n;

  /* 2. covariance matrix, column-major 3x3 (always size 3 for simplicity;
   * the unused tail entries stay 0). LAPACKE_dsyev expects column-major
   * symmetric matrix; we fill the upper triangle. */
  double C[9] = {0.0};
  const int ldC = (ndim > 0) ? ndim : 1;
  for (int i = 1; i <= cur_nd; i++) {
    int idx = cur_lod[i];
    double y[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < ndim; j++) y[j] = zgmid_t[j+1][idx] - cm[j];
    for (int j = 0; j < ndim; j++) {
      for (int k = j; k < ndim; k++) {
        C[j + k*ldC] += y[j] * y[k];   /* upper triangle, column-major */
      }
    }
  }

  /* 3. eigendecomposition. dsyev: 'V' = compute vectors, 'U' = upper triangle
   * stored. After return, C contains eigenvectors as columns, w contains
   * eigenvalues in ASCENDING order. */
  double w[3] = {0.0, 0.0, 0.0};
  int info = (int)LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', ndim, C, ldC, w);
  if (info != 0) return -1;

  /* 4. largest eigenvalue is the last; corresponding eigenvector is column
   * (ndim-1) of the overwritten C. Guard against rank-deficient case
   * (smallest eigenvalue ~ 0 is fine; we just need the largest to be > 0). */
  if (w[ndim-1] <= 0.0) return -1;
  double v[3] = {0.0, 0.0, 0.0};
  for (int j = 0; j < ndim; j++) v[j] = C[j + (ndim-1)*ldC];

  /* 5. partition by sign of (x_i - cm) . v. We use a 2-pointer in-place
   * swap (same logic as the existing bbox path in cHACApK_generate_cbitree).
   * Left = sign <= 0, right = sign > 0. */
  int nl = 1, nr = cur_nd;
  while (nl < nr) {
    while (nl < cur_nd) {
      int idx = cur_lod[nl];
      double s = 0.0;
      for (int j = 0; j < ndim; j++) s += (zgmid_t[j+1][idx] - cm[j]) * v[j];
      if (s <= 0.0) nl++;
      else break;
    }
    while (nr >= 1) {
      int idx = cur_lod[nr];
      double s = 0.0;
      for (int j = 0; j < ndim; j++) s += (zgmid_t[j+1][idx] - cm[j]) * v[j];
      if (s > 0.0) nr--;
      else break;
    }
    if (nl < nr) {
      int tmp = cur_lod[nl];
      cur_lod[nl] = cur_lod[nr];
      cur_lod[nr] = tmp;
    }
  }

  /* nl now is the 1-based index of the FIRST right-side entry; the left side
   * count is (nl - 1). Return as the C-style count (1-based, see header). */
  int n_left = nl - 1;   /* points on left side (inclusive of cur_lod[1..n_left]) */

  /* degenerate split (all on one side): caller falls back to bbox path */
  if (n_left <= 0 || n_left >= cur_nd) return -1;

  return n_left;
}
