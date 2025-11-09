/*-------------------------------------------------------------------------
*
* File name:      radintrc_hmat.h
*
* Project:        RADIA
*
* Description:    H-matrix acceleration for relaxation solver
*                 Provides fast interaction matrix construction and
*                 matrix-vector multiplication using HACApK library
*
* Author(s):      Implemented with Claude Code
*
* First release:  2025-11-07
*
* Copyright (C):  2025
*
*-------------------------------------------------------------------------*/

#ifndef __RADINTRC_HMAT_H
#define __RADINTRC_HMAT_H

#include "rad_interaction.h"
#include "rad_geometry_3d.h"
#include "rad_group.h"
#include "../ext/HACApK_LH-Cimplm/hacapk.hpp"  // HACApK library
#include <vector>
#include <memory>

//-------------------------------------------------------------------------
// Forward declarations
//-------------------------------------------------------------------------

class radTInteraction;

//-------------------------------------------------------------------------
// H-Matrix configuration for relaxation solver
//-------------------------------------------------------------------------

struct radTHMatrixSolverConfig
{
	double eps;              // ACA tolerance (default: 1e-6)
	int max_rank;            // Maximum rank for low-rank blocks (default: 50)
	int min_cluster_size;    // Minimum cluster size (default: 10)
	bool use_openmp;         // Enable OpenMP parallelization (default: true)
	int num_threads;         // Number of OpenMP threads (0 = auto-detect)

	radTHMatrixSolverConfig()
	{
		eps = 1e-6;
		max_rank = 50;
		min_cluster_size = 10;
		use_openmp = true;
		num_threads = 0;  // Auto-detect
	}

	radTHMatrixSolverConfig(double e, int mr, int mcs, bool omp, int nt)
		: eps(e), max_rank(mr), min_cluster_size(mcs), use_openmp(omp), num_threads(nt)
	{
	}
};

//-------------------------------------------------------------------------
// H-Matrix-accelerated interaction matrix
//
// Purpose: Replaces dense NxN interaction matrix with H-matrix
//          representation, providing O(N log N) operations instead of O(N^2)
//
// Usage:
//   1. Create: radTHMatrixInteraction hmat(interaction_ptr, config);
//   2. Build: hmat.BuildHMatrix();
//   3. Use: hmat.MatVec(M_in, H_out);
//-------------------------------------------------------------------------

class radTHMatrixInteraction
{
public:
	radTInteraction* intrct_ptr;     // Parent interaction object
	radTHMatrixSolverConfig config;  // H-matrix configuration

	// Element data
	int n_elem;                      // Number of relaxation elements
	double* elem_coords;             // Element coordinates [3*n_elem]
	radTg3dRelax** elem_ptrs;        // Pointers to element objects

	// HACApK data structures
	std::vector<hacapk::Point3D> points;  // Element center points
	hacapk::ControlParams hacapk_params;  // HACApK parameters

	// H-matrix storage for 3x3 tensor interaction matrix
	// We store 9 scalar H-matrices, one for each tensor component
	std::unique_ptr<hacapk::HMatrix> hmat[9];  // [0]=M[0][0], [1]=M[0][1], ..., [8]=M[2][2]

	// Cached symmetry transformations for each element (for performance)
	std::vector<std::vector<radTrans*>> cached_trans_vect;  // [j] = list of transformations for element j

	bool is_built;                   // H-matrix built flag

	// Statistics
	size_t memory_used;              // Estimated memory usage (bytes)
	double compression_ratio;        // Compression ratio vs dense
	double construction_time;        // Construction time (seconds)

public:
	radTHMatrixInteraction(radTInteraction* intrct, const radTHMatrixSolverConfig& cfg);
	~radTHMatrixInteraction();

	// Build H-matrix from interaction data
	int BuildHMatrix();

	// H-matrix-vector multiplication
	// Computes: H_field = InteractMatrix * M_vector
	// Input: M_in[n_elem] - magnetization vectors
	// Output: H_out[n_elem] - field vectors (without external field)
	void MatVec(const TVector3d* M_in, TVector3d* H_out);

	// Memory and statistics
	void PrintStatistics();
	size_t EstimateMemoryUsage() const;
	double GetCompressionRatio() const;

private:
	// Extract element coordinates from radTInteraction
	void ExtractElementData();

	// Kernel function for interaction matrix computation
	// Computes the 3x3 interaction matrix between elements i and j
	void ComputeInteractionKernel(int i, int j, TMatrix3df& result);

	// Kernel function wrapper data (for HACApK callback)
	struct KernelData
	{
		radTHMatrixInteraction* hmat_ptr;  // Pointer to this object
		int tensor_row;                     // Tensor component row (0, 1, 2)
		int tensor_col;                     // Tensor component col (0, 1, 2)
	};

	// Static kernel function for HACApK (extracts tensor component)
	static double KernelFunction(int i, int j, void* user_data);
};

#endif
