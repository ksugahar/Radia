/**
 * @file hacapk_openmp.c
 * @brief OpenMP-based H-matrix implementation
 *
 * @date 2025-11-07
 */

#include "hacapk_openmp.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

// ============================================================================
// Bounding Box Functions
// ============================================================================

void hmatrix_compute_bbox(
	const Point3D *points,
	const int *indices,
	int start,
	int size,
	BoundingBox *bbox
) {
	// Initialize with first point
	int idx = indices[start];
	bbox->min[0] = bbox->max[0] = points[idx].x;
	bbox->min[1] = bbox->max[1] = points[idx].y;
	bbox->min[2] = bbox->max[2] = points[idx].z;

	// Expand to include all points
	for (int i = 1; i < size; i++) {
		idx = indices[start + i];
		const Point3D *p = &points[idx];

		if (p->x < bbox->min[0]) bbox->min[0] = p->x;
		if (p->x > bbox->max[0]) bbox->max[0] = p->x;
		if (p->y < bbox->min[1]) bbox->min[1] = p->y;
		if (p->y > bbox->max[1]) bbox->max[1] = p->y;
		if (p->z < bbox->min[2]) bbox->min[2] = p->z;
		if (p->z > bbox->max[2]) bbox->max[2] = p->z;
	}
}

double bbox_distance(const BoundingBox *b1, const BoundingBox *b2) {
	double dist_sq = 0.0;

	for (int i = 0; i < 3; i++) {
		double d;
		if (b1->max[i] < b2->min[i]) {
			d = b2->min[i] - b1->max[i];
		} else if (b2->max[i] < b1->min[i]) {
			d = b1->min[i] - b2->max[i];
		} else {
			d = 0.0;  // Boxes overlap in this dimension
		}
		dist_sq += d * d;
	}

	return sqrt(dist_sq);
}

int hmatrix_is_admissible(
	const BoundingBox *bbox1,
	const BoundingBox *bbox2,
	double eta
) {
	// Admissibility condition: distance >= eta * min(diam1, diam2)
	double dist = bbox_distance(bbox1, bbox2);
	double diam1 = bbox_diameter(bbox1);
	double diam2 = bbox_diameter(bbox2);
	double min_diam = (diam1 < diam2) ? diam1 : diam2;

	return (dist >= eta * min_diam);
}

// ============================================================================
// Cluster Tree Functions
// ============================================================================

Cluster* hmatrix_create_cluster(
	const Point3D *points,
	int *indices,
	int start,
	int size,
	int depth,
	int max_leaf_size
) {
	Cluster *cluster = (Cluster*)malloc(sizeof(Cluster));
	cluster->start_idx = start;
	cluster->size = size;
	cluster->depth = depth;
	cluster->child[0] = NULL;
	cluster->child[1] = NULL;

	// Compute bounding box
	hmatrix_compute_bbox(points, indices, start, size, &cluster->bbox);

	// If small enough, make it a leaf
	if (size <= max_leaf_size) {
		return cluster;
	}

	// Find longest dimension
	double dx = cluster->bbox.max[0] - cluster->bbox.min[0];
	double dy = cluster->bbox.max[1] - cluster->bbox.min[1];
	double dz = cluster->bbox.max[2] - cluster->bbox.min[2];

	int split_dim = 0;
	if (dy > dx && dy > dz) split_dim = 1;
	else if (dz > dx && dz > dy) split_dim = 2;

	// Find split point (median)
	double split_val = (cluster->bbox.min[split_dim] + cluster->bbox.max[split_dim]) / 2.0;

	// Partition indices around split
	int left = start;
	int right = start + size - 1;

	while (left <= right) {
		const Point3D *p = &points[indices[left]];
		double val = (split_dim == 0) ? p->x : (split_dim == 1) ? p->y : p->z;

		if (val < split_val) {
			left++;
		} else {
			// Swap
			int temp = indices[left];
			indices[left] = indices[right];
			indices[right] = temp;
			right--;
		}
	}

	int left_size = left - start;
	int right_size = size - left_size;

	// Ensure both children have at least one element
	if (left_size == 0 || right_size == 0) {
		left_size = size / 2;
		right_size = size - left_size;
		left = start + left_size;
	}

	// Recursively create children
	cluster->child[0] = hmatrix_create_cluster(points, indices, start, left_size, depth + 1, max_leaf_size);
	cluster->child[1] = hmatrix_create_cluster(points, indices, left, right_size, depth + 1, max_leaf_size);

	return cluster;
}

void hmatrix_free_cluster(Cluster *cluster) {
	if (!cluster) return;

	if (cluster->child[0]) hmatrix_free_cluster(cluster->child[0]);
	if (cluster->child[1]) hmatrix_free_cluster(cluster->child[1]);

	free(cluster);
}

// ============================================================================
// ACA (Adaptive Cross Approximation)
// ============================================================================

LowRankBlock* hmatrix_aca(
	int m, int n,
	double (*kernel_func)(int i, int j, void *data),
	void *kernel_data,
	double epsilon
) {
	// Allocate low-rank block
	LowRankBlock *block = (LowRankBlock*)malloc(sizeof(LowRankBlock));
	block->m = m;
	block->n = n;
	block->rank = 0;

	int max_rank = (m < n) ? m : n;
	if (max_rank > 50) max_rank = 50;  // Limit rank for memory

	// Allocate temporary storage
	block->U = (double*)malloc(m * max_rank * sizeof(double));
	block->V = (double*)malloc(n * max_rank * sizeof(double));

	// ACA algorithm
	double *row = (double*)malloc(n * sizeof(double));
	double *col = (double*)malloc(m * sizeof(double));
	int *used_rows = (int*)calloc(m, sizeof(int));
	int *used_cols = (int*)calloc(n, sizeof(int));

	double norm_A = 0.0;
	double norm_R = 0.0;

	for (int k = 0; k < max_rank; k++) {
		// Find pivot row (not yet used, with maximum residual)
		int pivot_i = -1;
		double max_val = 0.0;

		for (int i = 0; i < m; i++) {
			if (used_rows[i]) continue;

			// Compute residual row
			for (int j = 0; j < n; j++) {
				double val = kernel_func(i, j, kernel_data);

				// Subtract previous approximation
				for (int r = 0; r < k; r++) {
					val -= block->U[i * max_rank + r] * block->V[j * max_rank + r];
				}

				row[j] = val;
			}

			// Find maximum in this row
			for (int j = 0; j < n; j++) {
				if (fabs(row[j]) > max_val) {
					max_val = fabs(row[j]);
					pivot_i = i;
				}
			}
		}

		if (pivot_i == -1 || max_val < epsilon) break;

		// Compute full residual row for pivot
		for (int j = 0; j < n; j++) {
			double val = kernel_func(pivot_i, j, kernel_data);

			for (int r = 0; r < k; r++) {
				val -= block->U[pivot_i * max_rank + r] * block->V[j * max_rank + r];
			}

			row[j] = val;
		}

		// Find pivot column
		int pivot_j = 0;
		max_val = fabs(row[0]);
		for (int j = 1; j < n; j++) {
			if (fabs(row[j]) > max_val) {
				max_val = fabs(row[j]);
				pivot_j = j;
			}
		}

		if (max_val < epsilon) break;

		// Compute residual column
		for (int i = 0; i < m; i++) {
			double val = kernel_func(i, pivot_j, kernel_data);

			for (int r = 0; r < k; r++) {
				val -= block->U[i * max_rank + r] * block->V[pivot_j * max_rank + r];
			}

			col[i] = val / row[pivot_j];  // Normalize
		}

		// Store rank-1 update
		for (int i = 0; i < m; i++) {
			block->U[i * max_rank + k] = col[i];
		}
		for (int j = 0; j < n; j++) {
			block->V[j * max_rank + k] = row[j];
		}

		used_rows[pivot_i] = 1;
		used_cols[pivot_j] = 1;

		// Update norms
		norm_R += max_val * max_val;
		norm_A += norm_R;

		// Check convergence
		if (norm_R < epsilon * epsilon * norm_A) {
			block->rank = k + 1;
			break;
		}

		block->rank = k + 1;
	}

	// Reallocate to actual rank
	if (block->rank < max_rank) {
		block->U = (double*)realloc(block->U, m * block->rank * sizeof(double));
		block->V = (double*)realloc(block->V, n * block->rank * sizeof(double));
	}

	free(row);
	free(col);
	free(used_rows);
	free(used_cols);

	return block;
}

void hmatrix_free_lowrank(LowRankBlock *block) {
	if (!block) return;
	free(block->U);
	free(block->V);
	free(block);
}

// ============================================================================
// H-Matrix Construction (Simplified - single block for now)
// ============================================================================

HMatrix* hmatrix_create(
	const Point3D *source_points,
	int n_source,
	const Point3D *target_points,
	int n_target,
	double eta,
	double epsilon,
	int max_leaf_size
) {
	HMatrix *hmat = (HMatrix*)malloc(sizeof(HMatrix));
	hmat->n_source = n_source;
	hmat->n_target = n_target;
	hmat->eta = eta;
	hmat->epsilon = epsilon;
	hmat->blocks = NULL;
	hmat->n_blocks = 0;

	// Create index arrays
	int *source_indices = (int*)malloc(n_source * sizeof(int));
	int *target_indices = (int*)malloc(n_target * sizeof(int));

	for (int i = 0; i < n_source; i++) source_indices[i] = i;
	for (int i = 0; i < n_target; i++) target_indices[i] = i;

	// Build cluster trees
	hmat->source_tree = hmatrix_create_cluster(source_points, source_indices, 0, n_source, 0, max_leaf_size);
	hmat->target_tree = hmatrix_create_cluster(target_points, target_indices, 0, n_target, 0, max_leaf_size);

	free(source_indices);
	free(target_indices);

	return hmat;
}

void hmatrix_free(HMatrix *hmat) {
	if (!hmat) return;

	hmatrix_free_cluster(hmat->source_tree);
	hmatrix_free_cluster(hmat->target_tree);

	for (int i = 0; i < hmat->n_blocks; i++) {
		hmatrix_free_lowrank(hmat->blocks[i]);
	}
	free(hmat->blocks);

	free(hmat);
}

// ============================================================================
// Matrix-Vector Multiplication (OpenMP parallelized)
// ============================================================================

void hmatrix_matvec(
	const HMatrix *hmat,
	const double *x,
	double *y
) {
	int b;

	// Initialize output
	memset(y, 0, hmat->n_target * sizeof(double));

	// For each low-rank block: y += U * (V^T * x)
	// Note: OpenMP 2.0 (MSVC) requires loop variable to be declared before pragma
	#pragma omp parallel for private(b)
	for (b = 0; b < hmat->n_blocks; b++) {
		const LowRankBlock *block = hmat->blocks[b];
		double *temp;
		int r, i, j;

		// Compute temp = V^T * x
		temp = (double*)malloc(block->rank * sizeof(double));

		for (r = 0; r < block->rank; r++) {
			temp[r] = 0.0;
			for (j = 0; j < block->n; j++) {
				temp[r] += block->V[j * block->rank + r] * x[j];
			}
		}

		// Compute y += U * temp
		#pragma omp critical
		{
			for (i = 0; i < block->m; i++) {
				for (r = 0; r < block->rank; r++) {
					y[i] += block->U[i * block->rank + r] * temp[r];
				}
			}
		}

		free(temp);
	}
}
