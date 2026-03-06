/**
 * @file hacapk_openmp.h
 * @brief OpenMP-based H-matrix implementation for Radia
 *
 * This is a simplified H-matrix implementation using OpenMP instead of MPI.
 * Designed specifically for Radia's magnetic field calculations.
 *
 * Key differences from HACApK:
 * - OpenMP (shared memory) instead of MPI (distributed memory)
 * - Simplified for single-node parallelization
 * - Integrated with Radia's existing OpenMP infrastructure
 *
 * @date 2025-11-07
 */

#ifndef HACAPK_OPENMP_H
#define HACAPK_OPENMP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <omp.h>

// ============================================================================
// Data Structures
// ============================================================================

/**
 * 3D Point structure
 */
typedef struct {
	double x, y, z;
} Point3D;

/**
 * Bounding box for spatial clustering
 */
typedef struct {
	double min[3];  // Minimum coordinates (x, y, z)
	double max[3];  // Maximum coordinates (x, y, z)
} BoundingBox;

/**
 * Cluster node in hierarchical tree
 */
typedef struct Cluster {
	int start_idx;          // Starting index in point array
	int size;               // Number of points in this cluster
	int depth;              // Depth in tree
	BoundingBox bbox;       // Bounding box
	struct Cluster *child[2]; // Binary tree: left and right children (NULL if leaf)
} Cluster;

/**
 * Low-rank matrix block (ACA approximation)
 * Represents: M ≈ U * V^T
 */
typedef struct {
	double *U;     // Left factor matrix (m × rank)
	double *V;     // Right factor matrix (n × rank)
	int m, n;      // Dimensions
	int rank;      // Approximation rank
} LowRankBlock;

/**
 * H-matrix structure
 */
typedef struct {
	int n_source;          // Number of source points
	int n_target;          // Number of target points
	Cluster *source_tree;  // Source point cluster tree
	Cluster *target_tree;  // Target point cluster tree
	LowRankBlock **blocks; // Low-rank blocks
	int n_blocks;          // Number of blocks
	double eta;            // Admissibility parameter
	double epsilon;        // ACA tolerance
} HMatrix;

// ============================================================================
// Cluster Tree Functions
// ============================================================================

/**
 * Create a cluster from points
 */
Cluster* hmatrix_create_cluster(
	const Point3D *points,
	int *indices,
	int start,
	int size,
	int depth,
	int max_leaf_size
);

/**
 * Free cluster tree
 */
void hmatrix_free_cluster(Cluster *cluster);

/**
 * Compute bounding box for a set of points
 */
void hmatrix_compute_bbox(
	const Point3D *points,
	const int *indices,
	int start,
	int size,
	BoundingBox *bbox
);

/**
 * Check if two clusters are admissible (far-field)
 * Returns 1 if admissible, 0 otherwise
 */
int hmatrix_is_admissible(
	const BoundingBox *bbox1,
	const BoundingBox *bbox2,
	double eta
);

// ============================================================================
// ACA (Adaptive Cross Approximation) Functions
// ============================================================================

/**
 * Compute low-rank approximation using ACA
 * kernel_func: Function to compute matrix element (i, j)
 */
LowRankBlock* hmatrix_aca(
	int m, int n,
	double (*kernel_func)(int i, int j, void *data),
	void *kernel_data,
	double epsilon
);

/**
 * Free low-rank block
 */
void hmatrix_free_lowrank(LowRankBlock *block);

// ============================================================================
// H-Matrix Construction
// ============================================================================

/**
 * Create H-matrix structure
 */
HMatrix* hmatrix_create(
	const Point3D *source_points,
	int n_source,
	const Point3D *target_points,
	int n_target,
	double eta,        // Admissibility parameter (typically 2.0)
	double epsilon,    // ACA tolerance (typically 1e-6)
	int max_leaf_size  // Maximum points per leaf (typically 50)
);

/**
 * Free H-matrix
 */
void hmatrix_free(HMatrix *hmat);

// ============================================================================
// Matrix-Vector Multiplication (OpenMP parallelized)
// ============================================================================

/**
 * Compute y = H * x using H-matrix
 * Uses OpenMP for parallelization
 */
void hmatrix_matvec(
	const HMatrix *hmat,
	const double *x,
	double *y
);

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Get distance between two points
 */
static inline double point_distance(const Point3D *p1, const Point3D *p2) {
	double dx = p1->x - p2->x;
	double dy = p1->y - p2->y;
	double dz = p1->z - p2->z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}

/**
 * Get distance between two bounding boxes (minimum distance)
 */
double bbox_distance(const BoundingBox *b1, const BoundingBox *b2);

/**
 * Get diameter of bounding box
 */
static inline double bbox_diameter(const BoundingBox *bbox) {
	double dx = bbox->max[0] - bbox->min[0];
	double dy = bbox->max[1] - bbox->min[1];
	double dz = bbox->max[2] - bbox->min[2];
	return sqrt(dx*dx + dy*dy + dz*dz);
}

// ============================================================================
// OpenMP Configuration
// ============================================================================

/**
 * Set number of OpenMP threads
 */
static inline void hmatrix_set_threads(int n_threads) {
	omp_set_num_threads(n_threads);
}

/**
 * Get number of OpenMP threads
 */
static inline int hmatrix_get_threads(void) {
	return omp_get_max_threads();
}

#ifdef __cplusplus
}
#endif

#endif // HACAPK_OPENMP_H
