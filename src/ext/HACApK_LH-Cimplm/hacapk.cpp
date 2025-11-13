/*
!=====================================================================*
!                                                                     *
!   Software Name : HACApK                                            *
!         Version : 1.3.0                                             *
!                                                                     *
!   License                                                           *
!     This file is part of HACApK.                                    *
!     HACApK is a free software, you can use it under the terms       *
!     of The MIT License (MIT). See LICENSE file and User's guide     *
!     for more details.                                               *
!                                                                     *
!   ppOpen-HPC project:                                               *
!     Open Source Infrastructure for Development and Execution of     *
!     Large-Scale Scientific Applications on Post-Peta-Scale          *
!     Supercomputers with Automatic Tuning (AT).                      *
!                                                                     *
!   Sponsorship:                                                      *
!     Japan Science and Technology Agency (JST), Basic Research       *
!     Programs: CREST, Development of System Software Technologies    *
!     for post-Peta Scale High Performance Computing.                 *
!                                                                     *
!   Copyright (c) 2015 <Akihiro Ida and Takeshi Iwashita>             *
!                                                                     *
!=====================================================================*
 * @file hacapk.cpp
 * @brief HACApK C++ implementation
 *
 * @date 2025-11-07
*/

#include "hacapk.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

namespace hacapk {

// ============================================================================
// BoundingBox Implementation
// ============================================================================

double BoundingBox::width() const {
	double max_width = 0.0;
	for (size_t i = 0; i < min.size(); i++) {
		double w = max[i] - min[i];
		if (w > max_width) max_width = w;
	}
	return max_width;
}

double BoundingBox::diameter() const {
	double sum_sq = 0.0;
	for (size_t i = 0; i < min.size(); i++) {
		double d = max[i] - min[i];
		sum_sq += d * d;
	}
	return std::sqrt(sum_sq);
}

double bbox_distance(const BoundingBox& b1, const BoundingBox& b2) {
	double dist_sq = 0.0;

	for (size_t i = 0; i < b1.min.size(); i++) {
		double d = 0.0;
		if (b1.max[i] < b2.min[i]) {
			d = b2.min[i] - b1.max[i];
		} else if (b2.max[i] < b1.min[i]) {
			d = b1.min[i] - b2.max[i];
		}
		// Otherwise boxes overlap in this dimension
		dist_sq += d * d;
	}

	return std::sqrt(dist_sq);
}

// ============================================================================
// Cluster Implementation
// ============================================================================

Cluster::Cluster(int ndim_)
	: ndim(ndim_)
	, nstrt(0)
	, nsize(0)
	, ndpth(0)
	, nnson(0)
	, nmbr(0)
	, ndscd(0)
	, bbox(ndim_)
	, zwdth(0.0)
{
}

void compute_bounding_box(
	Cluster& cluster,
	const std::vector<Point3D>& points,
	const std::vector<int>& indices
) {
	if (cluster.nsize == 0) return;

	// Initialize with first point
	int idx = indices[cluster.nstrt];
	cluster.bbox.min[0] = cluster.bbox.max[0] = points[idx].x;
	cluster.bbox.min[1] = cluster.bbox.max[1] = points[idx].y;
	cluster.bbox.min[2] = cluster.bbox.max[2] = points[idx].z;

	// Expand to include all points
	for (int i = 1; i < cluster.nsize; i++) {
		idx = indices[cluster.nstrt + i];
		const Point3D& p = points[idx];

		cluster.bbox.min[0] = std::min(cluster.bbox.min[0], p.x);
		cluster.bbox.max[0] = std::max(cluster.bbox.max[0], p.x);
		cluster.bbox.min[1] = std::min(cluster.bbox.min[1], p.y);
		cluster.bbox.max[1] = std::max(cluster.bbox.max[1], p.y);
		cluster.bbox.min[2] = std::min(cluster.bbox.min[2], p.z);
		cluster.bbox.max[2] = std::max(cluster.bbox.max[2], p.z);
	}

	cluster.zwdth = cluster.bbox.width();
}

bool is_admissible(
	const BoundingBox& box1,
	const BoundingBox& box2,
	double eta
) {
	double dist = bbox_distance(box1, box2);
	double diam1 = box1.diameter();
	double diam2 = box2.diameter();
	double min_diam = std::min(diam1, diam2);

	return (dist >= eta * min_diam);
}

// ============================================================================
// Cluster Tree Generation
// ============================================================================

std::shared_ptr<Cluster> generate_cluster(
	const std::vector<Point3D>& points,
	std::vector<int>& indices,
	int start,
	int size,
	int depth,
	const ControlParams& params
) {
	auto cluster = std::make_shared<Cluster>(3);
	cluster->nstrt = start;
	cluster->nsize = size;
	cluster->ndpth = depth;

	// Compute bounding box
	compute_bounding_box(*cluster, points, indices);

	// Check if this should be a leaf
	if (size <= params.leaf_size) {
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
		const Point3D& p = points[indices[left]];
		double val = (split_dim == 0) ? p.x : (split_dim == 1) ? p.y : p.z;

		if (val < split_val) {
			left++;
		} else {
			std::swap(indices[left], indices[right]);
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

	// Recursively create children (binary tree)
	cluster->nnson = 2;
	cluster->sons.resize(2);
	cluster->sons[0] = generate_cluster(points, indices, start, left_size, depth + 1, params);
	cluster->sons[1] = generate_cluster(points, indices, left, right_size, depth + 1, params);

	cluster->ndscd = left_size + right_size;

	return cluster;
}

// ============================================================================
// LowRankBlock Implementation
// ============================================================================

LowRankBlock::LowRankBlock()
	: ltmtx(0)
	, kt(0)
	, nstrtl(0)
	, ndl(0)
	, nstrtt(0)
	, ndt(0)
{
}

size_t LowRankBlock::memory_usage() const {
	size_t mem = 0;
	if (is_lowrank()) {
		mem = a1.size() * sizeof(double) + a2.size() * sizeof(double);
	} else if (is_full()) {
		mem = a1.size() * sizeof(double);
	} else if (is_hierarchical() || is_block()) {
		for (const auto& sub : sublocks) {
			mem += sub.memory_usage();
		}
	}
	return mem;
}

// ============================================================================
// HMatrix Implementation
// ============================================================================

ControlParams::ControlParams()
	: param(100, 0.0)
	, time(10, 0.0)
	, nthr(omp_get_max_threads())
	, print_level(1)
	, leaf_size(15)
	, max_leaf_size_ratio(1.0)
	, eta(2.0)
	, eps_aca(1e-6)
	, aca_type(2)
{
	param[1] = 1.0;    // Print level
	param[21] = 15.0;  // Leaf size
	param[22] = 1.0;   // Max leaf size ratio
	param[51] = 2.0;   // Eta parameter
	param[60] = 2.0;   // ACA type (2=ACA+)
	param[63] = 1e-6;  // ACA tolerance
}

HMatrix::HMatrix()
	: nd(0)
	, nlf(0)
	, nlfkt(0)
	, ktmax(0)
{
}

size_t HMatrix::memory_usage() const {
	size_t mem = 0;
	for (const auto& block : blocks) {
		mem += block.memory_usage();
	}
	return mem;
}

double HMatrix::compression_ratio() const {
	if (nd == 0) return 0.0;

	size_t hmatrix_mem = memory_usage();
	size_t full_mem = static_cast<size_t>(nd) * static_cast<size_t>(nd) * sizeof(double);

	return static_cast<double>(full_mem) / static_cast<double>(hmatrix_mem);
}

// ============================================================================
// ACA Implementation
// ============================================================================

void aca_approximation(
	LowRankBlock& block,
	KernelFunction kernel,
	void* kernel_data,
	double eps
) {
	int m = block.ndl;
	int n = block.ndt;

	// Ensure kernel function is valid
	if (!kernel) {
		block.kt = 0;
		block.ltmtx = 2;  // Mark as full matrix (not approximated)
		return;
	}

	int max_rank = std::min(m, n);
	if (max_rank > 50) max_rank = 50;  // Limit rank
	if (max_rank == 0) {
		block.kt = 0;
		block.ltmtx = 2;
		return;
	}

	// Allocate temporary storage
	std::vector<double> U_temp(m * max_rank, 0.0);
	std::vector<double> V_temp(n * max_rank, 0.0);

	std::vector<double> row(n, 0.0);
	std::vector<double> col(m, 0.0);
	std::vector<int> used_rows(m, 0);
	std::vector<int> used_cols(n, 0);

	double norm_A = 0.0;
	double norm_R = 0.0;
	int rank = 0;

	// Debug: Print initial state
	// std::cout << "[DEBUG] ACA: m=" << m << ", n=" << n << ", max_rank=" << max_rank << ", eps=" << eps << std::endl;

	for (int k = 0; k < max_rank; k++) {
		// Find pivot row
		int pivot_i = 0;  // Start with first unused row
		double max_val = 0.0;

		// Find pivot row with maximum residual
		for (int i = 0; i < m; i++) {
			if (used_rows[i]) continue;

			// Compute residual row
			for (int j = 0; j < n; j++) {
				double val = kernel(block.nstrtl + i, block.nstrtt + j, kernel_data);

				// Subtract previous approximation
				for (int r = 0; r < k; r++) {
					val -= U_temp[i * max_rank + r] * V_temp[j * max_rank + r];
				}

				row[j] = val;
			}

			// Find maximum in this row
			double row_max = 0.0;
			for (int j = 0; j < n; j++) {
				row_max = std::max(row_max, std::abs(row[j]));
			}

			if (row_max > max_val) {
				max_val = row_max;
				pivot_i = i;
			}
		}

		// std::cout << "[DEBUG] ACA iter " << k << ": pivot_i=" << pivot_i << ", max_val=" << max_val << ", eps=" << eps << std::endl;

		if (max_val < eps) {
			// std::cout << "[DEBUG] ACA stopped: max_val < eps" << std::endl;
			break;
		}

		// Compute full residual row for pivot
		for (int j = 0; j < n; j++) {
			double val = kernel(block.nstrtl + pivot_i, block.nstrtt + j, kernel_data);

			for (int r = 0; r < k; r++) {
				val -= U_temp[pivot_i * max_rank + r] * V_temp[j * max_rank + r];
			}

			row[j] = val;
		}

		// Find pivot column
		int pivot_j = 0;
		max_val = std::abs(row[0]);
		for (int j = 1; j < n; j++) {
			if (std::abs(row[j]) > max_val) {
				max_val = std::abs(row[j]);
				pivot_j = j;
			}
		}

		if (max_val < eps) break;

		// Compute residual column
		for (int i = 0; i < m; i++) {
			double val = kernel(block.nstrtl + i, block.nstrtt + pivot_j, kernel_data);

			for (int r = 0; r < k; r++) {
				val -= U_temp[i * max_rank + r] * V_temp[pivot_j * max_rank + r];
			}

			col[i] = val / row[pivot_j];  // Normalize
		}

		// Store rank-1 update
		for (int i = 0; i < m; i++) {
			U_temp[i * max_rank + k] = col[i];
		}
		for (int j = 0; j < n; j++) {
			V_temp[j * max_rank + k] = row[j];
		}

		used_rows[pivot_i] = 1;
		used_cols[pivot_j] = 1;

		rank = k + 1;

		// Update norms
		norm_R += max_val * max_val;
		norm_A += norm_R;

		// Check convergence
		if (norm_R < eps * eps * norm_A) {
			// std::cout << "[DEBUG] ACA converged: norm_R=" << norm_R << ", eps^2*norm_A=" << (eps*eps*norm_A) << std::endl;
			break;
		}
	}

	// std::cout << "[DEBUG] ACA finished: rank=" << rank << std::endl;

	// Store final result
	block.kt = rank;
	block.ltmtx = 1;  // Low-rank

	block.a1.resize(m * rank);
	block.a2.resize(n * rank);

	for (int i = 0; i < m; i++) {
		for (int r = 0; r < rank; r++) {
			block.a1[i * rank + r] = U_temp[i * max_rank + r];
		}
	}

	for (int j = 0; j < n; j++) {
		for (int r = 0; r < rank; r++) {
			block.a2[j * rank + r] = V_temp[j * max_rank + r];
		}
	}
}

void aca_plus_approximation(
	LowRankBlock& block,
	KernelFunction kernel,
	void* kernel_data,
	double eps
) {
	// For now, use standard ACA
	// ACA+ implementation can be added later
	aca_approximation(block, kernel, kernel_data, eps);
}

// ============================================================================
// H-Matrix Construction
// ============================================================================

void generate_leaf_blocks(
	HMatrix& hmat,
	const std::shared_ptr<Cluster>& source_cluster,
	const std::shared_ptr<Cluster>& target_cluster,
	KernelFunction kernel,
	void* kernel_data,
	const ControlParams& params
) {
	// Check admissibility
	bool admissible = is_admissible(source_cluster->bbox, target_cluster->bbox, params.eta);

	// Both are leaves
	if (source_cluster->is_leaf() && target_cluster->is_leaf()) {
		LowRankBlock block;
		block.nstrtl = source_cluster->nstrt;
		block.ndl = source_cluster->nsize;
		block.nstrtt = target_cluster->nstrt;
		block.ndt = target_cluster->nsize;

		if (admissible) {
			// Use ACA for low-rank approximation
			if (params.aca_type == 2) {
				aca_plus_approximation(block, kernel, kernel_data, params.eps_aca);
			} else {
				aca_approximation(block, kernel, kernel_data, params.eps_aca);
			}
		} else {
			// Store as full matrix
			block.ltmtx = 2;
			block.kt = 0;
			int m = block.ndl;
			int n = block.ndt;
			block.a1.resize(m * n);

			// Compute full matrix using kernel function
			// Use same indexing as ACA
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					block.a1[i * n + j] = kernel(block.nstrtl + i, block.nstrtt + j, kernel_data);
				}
			}
		}

		hmat.blocks.push_back(block);
		hmat.nlf++;
		if (block.is_lowrank()) {
			hmat.nlfkt++;
			hmat.ktmax = std::max(hmat.ktmax, block.kt);
		}

		return;
	}

	// Recursively split
	if (!source_cluster->is_leaf() && !target_cluster->is_leaf()) {
		for (const auto& src_son : source_cluster->sons) {
			for (const auto& tgt_son : target_cluster->sons) {
				generate_leaf_blocks(hmat, src_son, tgt_son, kernel, kernel_data, params);
			}
		}
	} else if (!source_cluster->is_leaf()) {
		for (const auto& src_son : source_cluster->sons) {
			generate_leaf_blocks(hmat, src_son, target_cluster, kernel, kernel_data, params);
		}
	} else {
		for (const auto& tgt_son : target_cluster->sons) {
			generate_leaf_blocks(hmat, source_cluster, tgt_son, kernel, kernel_data, params);
		}
	}
}

std::unique_ptr<HMatrix> build_hmatrix(
	const std::vector<Point3D>& source_points,
	const std::vector<Point3D>& target_points,
	KernelFunction kernel,
	void* kernel_data,
	const ControlParams& params
) {
	auto hmat = std::make_unique<HMatrix>();
	hmat->nd = source_points.size();

	// Create index arrays
	std::vector<int> source_indices(source_points.size());
	std::vector<int> target_indices(target_points.size());

	for (size_t i = 0; i < source_points.size(); i++) source_indices[i] = i;
	for (size_t i = 0; i < target_points.size(); i++) target_indices[i] = i;

	// Build cluster trees
	auto source_tree = generate_cluster(source_points, source_indices, 0, source_points.size(), 0, params);
	auto target_tree = generate_cluster(target_points, target_indices, 0, target_points.size(), 0, params);

	// Generate leaf blocks
	generate_leaf_blocks(*hmat, source_tree, target_tree, kernel, kernel_data, params);

	return hmat;
}

// ============================================================================
// Matrix-Vector Multiplication
// ============================================================================

void lowrank_matvec(
	const LowRankBlock& block,
	const std::vector<double>& x,
	std::vector<double>& y
) {
	if (!block.is_lowrank()) return;

	int m = block.ndl;
	int n = block.ndt;
	int k = block.kt;

	// Compute temp = V^T * x
	std::vector<double> temp(k, 0.0);

	for (int r = 0; r < k; r++) {
		for (int j = 0; j < n; j++) {
			temp[r] += block.a2[j * k + r] * x[block.nstrtt + j];
		}
	}

	// Compute y += U * temp
	for (int i = 0; i < m; i++) {
		for (int r = 0; r < k; r++) {
			y[block.nstrtl + i] += block.a1[i * k + r] * temp[r];
		}
	}
}

void hmatrix_matvec(
	const HMatrix& hmat,
	const std::vector<double>& x,
	std::vector<double>& y
) {
	// Initialize output
	std::fill(y.begin(), y.end(), 0.0);

	// OpenMP parallelized loop over blocks
	int nblocks = hmat.blocks.size();

	#pragma omp parallel for
	for (int b = 0; b < nblocks; b++) {
		const LowRankBlock& block = hmat.blocks[b];

		if (block.is_lowrank()) {
			// Need temporary vector for thread safety
			std::vector<double> y_local(y.size(), 0.0);
			lowrank_matvec(block, x, y_local);

			// Add to global y (critical section)
			#pragma omp critical
			{
				for (size_t i = 0; i < y.size(); i++) {
					y[i] += y_local[i];
				}
			}
		}
		else if (block.is_full()) {
			// Full matrix multiplication: y = A * x
			// A[i,j] stored in row-major format
			int m = block.ndl;
			int n = block.ndt;

			std::vector<double> y_local(y.size(), 0.0);

			for (int i = 0; i < m; i++) {
				double sum = 0.0;
				for (int j = 0; j < n; j++) {
					sum += block.a1[i * n + j] * x[block.nstrtt + j];
				}
				y_local[block.nstrtl + i] = sum;
			}

			// Add to global y (critical section)
			#pragma omp critical
			{
				for (size_t i = 0; i < y.size(); i++) {
					y[i] += y_local[i];
				}
			}
		}
	}
}

} // namespace hacapk
