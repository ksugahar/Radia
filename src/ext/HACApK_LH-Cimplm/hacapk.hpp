/**
 * @file hacapk.hpp
 * @brief HACApK C++ implementation with OpenMP (Fortran-free version)
 *
 * This is a C++17 reimplementation of HACApK (Hierarchical Adaptive Cross
 * Approximation on Kernel matrices) originally written in Fortran + MPI.
 *
 * Key changes from original HACApK:
 * - C++17 instead of Fortran
 * - OpenMP instead of MPI for parallelization
 * - Modern C++ features (smart pointers, STL containers)
 * - Integrated with Radia's build system
 *
 * Original HACApK:
 * - Copyright (c) 2015 Akihiro Ida and Takeshi Iwashita
 * - License: MIT License
 * - ppOpen-HPC project
 *
 * @date 2025-11-07
 */

#ifndef HACAPK_HPP
#define HACAPK_HPP

#include <vector>
#include <memory>
#include <functional>
#include <cstdint>
#include <omp.h>

namespace hacapk {

// ============================================================================
// Data Structures
// ============================================================================

/**
 * 3D Point structure
 */
struct Point3D {
	double x, y, z;

	Point3D() : x(0), y(0), z(0) {}
	Point3D(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

/**
 * Bounding box for spatial clustering
 */
struct BoundingBox {
	std::vector<double> min;  // Minimum coordinates
	std::vector<double> max;  // Maximum coordinates

	BoundingBox(int ndim = 3) : min(ndim, 0.0), max(ndim, 0.0) {}

	double width() const;
	double diameter() const;
};

/**
 * Cluster node in hierarchical tree
 */
class Cluster {
public:
	int ndim;           // Dimension
	int nstrt;          // Starting index
	int nsize;          // Number of points
	int ndpth;          // Depth in tree
	int nnson;          // Number of sons (children)
	int nmbr;           // Member number
	int ndscd;          // Number of descendants

	BoundingBox bbox;   // Bounding box
	double zwdth;       // Width

	std::vector<std::shared_ptr<Cluster>> sons;  // Children

	Cluster(int ndim_ = 3);
	~Cluster() = default;

	bool is_leaf() const { return sons.empty(); }
};

/**
 * Low-rank matrix block (ACA approximation)
 * Represents: M approx U * V^T
 */
struct LowRankBlock {
	int ltmtx;          // Type: 1=low-rank, 2=full, 3=H-matrix, 4=block
	int kt;             // Rank for low-rank blocks

	int nstrtl, ndl;    // Row start and size
	int nstrtt, ndt;    // Column start and size

	std::vector<double> a1;  // U matrix (ndl * kt)
	std::vector<double> a2;  // V matrix (ndt * kt)

	// For hierarchical blocks
	std::vector<LowRankBlock> sublocks;

	LowRankBlock();
	~LowRankBlock() = default;

	bool is_lowrank() const { return ltmtx == 1; }
	bool is_full() const { return ltmtx == 2; }
	bool is_hierarchical() const { return ltmtx == 3; }
	bool is_block() const { return ltmtx == 4; }

	size_t memory_usage() const;
};

/**
 * H-matrix structure (leaf matrix parameters)
 */
class HMatrix {
public:
	int nd;             // Total number of unknowns
	int nlf;            // Number of leaf blocks in this process
	int nlfkt;          // Number of low-rank blocks
	int ktmax;          // Maximum rank

	std::vector<LowRankBlock> blocks;  // Leaf blocks

	// Block structure
	std::vector<int> lbstrtl, lbstrtt;  // Block start indices (row/col)
	std::vector<int> lbndl, lbndt;      // Block sizes (row/col)

	HMatrix();
	~HMatrix() = default;

	size_t memory_usage() const;
	double compression_ratio() const;
};

/**
 * Control parameters
 */
struct ControlParams {
	std::vector<double> param;  // Numerical parameters
	std::vector<double> time;   // Timing data

	int nthr;           // Number of OpenMP threads
	int print_level;    // Print level (0=errors, 1=std, 2=debug)

	// Clustering parameters
	double leaf_size;           // Leaf size (param[21])
	double max_leaf_size_ratio; // Max leaf size ratio (param[22])

	// H-matrix parameters
	double eta;                 // Distance parameter (param[51])
	double eps_aca;             // ACA tolerance (param[63])
	int aca_type;               // ACA type: 1=ACA, 2=ACA+ (param[60])

	ControlParams();
	~ControlParams() = default;
};

// ============================================================================
// Kernel Function Type
// ============================================================================

/**
 * Kernel function signature
 * @param i Row index
 * @param j Column index
 * @param user_data User-provided data pointer
 * @return Matrix element K(i,j)
 */
using KernelFunction = std::function<double(int, int, void*)>;

// ============================================================================
// Cluster Tree Functions
// ============================================================================

/**
 * Generate cluster tree from points
 */
std::shared_ptr<Cluster> generate_cluster(
	const std::vector<Point3D>& points,
	std::vector<int>& indices,
	int start,
	int size,
	int depth,
	const ControlParams& params
);

/**
 * Compute bounding box for cluster
 */
void compute_bounding_box(
	Cluster& cluster,
	const std::vector<Point3D>& points,
	const std::vector<int>& indices
);

/**
 * Check admissibility (far-field condition)
 */
bool is_admissible(
	const BoundingBox& box1,
	const BoundingBox& box2,
	double eta
);

// ============================================================================
// H-Matrix Construction
// ============================================================================

/**
 * Build H-matrix from kernel function
 */
std::unique_ptr<HMatrix> build_hmatrix(
	const std::vector<Point3D>& source_points,
	const std::vector<Point3D>& target_points,
	KernelFunction kernel,
	void* kernel_data,
	const ControlParams& params
);

/**
 * Generate leaf blocks recursively
 */
void generate_leaf_blocks(
	HMatrix& hmat,
	const std::shared_ptr<Cluster>& source_cluster,
	const std::shared_ptr<Cluster>& target_cluster,
	KernelFunction kernel,
	void* kernel_data,
	const ControlParams& params
);

// ============================================================================
// ACA (Adaptive Cross Approximation)
// ============================================================================

/**
 * ACA algorithm for low-rank approximation
 */
void aca_approximation(
	LowRankBlock& block,
	KernelFunction kernel,
	void* kernel_data,
	double eps
);

/**
 * ACA+ algorithm (improved version)
 */
void aca_plus_approximation(
	LowRankBlock& block,
	KernelFunction kernel,
	void* kernel_data,
	double eps
);

// ============================================================================
// Matrix-Vector Multiplication
// ============================================================================

/**
 * H-matrix matrix-vector product: y = H * x
 * OpenMP parallelized
 */
void hmatrix_matvec(
	const HMatrix& hmat,
	const std::vector<double>& x,
	std::vector<double>& y
);

/**
 * Low-rank block matrix-vector product: y += U * (V^T * x)
 */
void lowrank_matvec(
	const LowRankBlock& block,
	const std::vector<double>& x,
	std::vector<double>& y
);

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Distance between two points
 */
inline double point_distance(const Point3D& p1, const Point3D& p2) {
	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	double dz = p1.z - p2.z;
	return std::sqrt(dx*dx + dy*dy + dz*dz);
}

/**
 * Distance between two bounding boxes
 */
double bbox_distance(const BoundingBox& b1, const BoundingBox& b2);

/**
 * Set number of OpenMP threads
 */
inline void set_num_threads(int n) {
	omp_set_num_threads(n);
}

/**
 * Get number of OpenMP threads
 */
inline int get_num_threads() {
	return omp_get_max_threads();
}

} // namespace hacapk

#endif // HACAPK_HPP
