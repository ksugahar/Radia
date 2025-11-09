/*-------------------------------------------------------------------------
*
* File name:      radintrc_hmat.cpp
*
* Project:        RADIA
*
* Description:    H-matrix acceleration for relaxation solver
*                 Implementation of radTHMatrixInteraction class
*
* Author(s):      Implemented with Claude Code
*
* First release:  2025-11-07
*
* Copyright (C):  2025
*
*-------------------------------------------------------------------------*/

#include "rad_intrc_hmat.h"
#include "rad_geometry_3d.h"
#include "rad_group.h"
#include "rad_transform_def.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <chrono>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

//-------------------------------------------------------------------------
// Constructor
//-------------------------------------------------------------------------

radTHMatrixInteraction::radTHMatrixInteraction(radTInteraction* intrct, const radTHMatrixSolverConfig& cfg)
	: intrct_ptr(intrct), config(cfg)
{
	n_elem = 0;
	elem_coords = nullptr;
	elem_ptrs = nullptr;
	is_built = false;
	memory_used = 0;
	compression_ratio = 0.0;
	construction_time = 0.0;

	// Extract element data from interaction object
	ExtractElementData();
}

//-------------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------------

radTHMatrixInteraction::~radTHMatrixInteraction()
{
	// Clean up allocated memory
	if(elem_coords != nullptr)
	{
		delete[] elem_coords;
		elem_coords = nullptr;
	}

	if(elem_ptrs != nullptr)
	{
		delete[] elem_ptrs;
		elem_ptrs = nullptr;
	}

	// Note: H-matrix data structures will be cleaned up here
	// when HACApK integration is complete
}

//-------------------------------------------------------------------------
// Extract element coordinates and pointers from radTInteraction
//-------------------------------------------------------------------------

void radTHMatrixInteraction::ExtractElementData()
{
	if(intrct_ptr == nullptr)
	{
		throw std::runtime_error("H-matrix solver: Invalid interaction pointer");
	}

	n_elem = intrct_ptr->AmOfMainElem;

	if(n_elem <= 0)
	{
		throw std::runtime_error("H-matrix solver: No relaxation elements found");
	}

	// Allocate coordinate array
	elem_coords = new double[3 * n_elem];
	elem_ptrs = new radTg3dRelax*[n_elem];

	// Prepare HACApK point array
	points.reserve(n_elem);

	// Extract element centers and pointers
	for(int i = 0; i < n_elem; i++)
	{
		radTg3dRelax* elem = intrct_ptr->g3dRelaxPtrVect[i];
		elem_ptrs[i] = elem;

		// Get element center point
		TVector3d center = elem->ReturnCentrPoint();

		// Apply transformation if needed
		radTrans* trans = intrct_ptr->MainTransPtrArray[i];
		if(trans != nullptr)
		{
			center = trans->TrPoint(center);
		}

		// Store coordinates in legacy array
		elem_coords[3*i + 0] = center.x;
		elem_coords[3*i + 1] = center.y;
		elem_coords[3*i + 2] = center.z;

		// Store in HACApK Point3D array
		points.emplace_back(center.x, center.y, center.z);
	}

	// Initialize HACApK parameters
	hacapk_params.eps_aca = config.eps;
	hacapk_params.leaf_size = config.min_cluster_size;
	hacapk_params.eta = 0.8;  // Admissibility parameter
	hacapk_params.aca_type = 2;  // Use ACA+ (improved version)
	hacapk_params.nthr = config.num_threads;
	if(hacapk_params.nthr <= 0)
	{
#ifdef _OPENMP
		hacapk_params.nthr = omp_get_max_threads();
#else
		hacapk_params.nthr = 1;
#endif
	}
	hacapk_params.print_level = 1;  // Standard output

	// Cache symmetry transformations for each element (for performance)
	// This avoids calling FillInTransPtrVectForElem() millions of times
	std::cout << "Caching symmetry transformations for " << n_elem << " elements..." << std::endl;
	cached_trans_vect.resize(n_elem);
	for(int j = 0; j < n_elem; j++)
	{
		intrct_ptr->FillInTransPtrVectForElem(j, 'I');
		cached_trans_vect[j] = intrct_ptr->TransPtrVect;  // Copy the vector
		intrct_ptr->EmptyTransPtrVect();
	}
	std::cout << "Cached " << n_elem << " transformation lists" << std::endl;
}

//-------------------------------------------------------------------------
// Build H-matrix structure
// Phase 2: Full implementation
//-------------------------------------------------------------------------

int radTHMatrixInteraction::BuildHMatrix()
{
	auto t_start = std::chrono::high_resolution_clock::now();

	try
	{
		std::cout << "\n========================================" << std::endl;
		std::cout << "Building True H-Matrix for Relaxation Solver" << std::endl;
		std::cout << "========================================" << std::endl;
		std::cout << "Number of elements: " << n_elem << std::endl;
		std::cout << "ACA tolerance: " << hacapk_params.eps_aca << std::endl;
		std::cout << "Admissibility param: " << hacapk_params.eta << std::endl;
		std::cout << "Min cluster size: " << hacapk_params.leaf_size << std::endl;
		std::cout << "OpenMP threads: " << hacapk_params.nthr << std::endl;

		// Build 9 H-matrices for the 3x3 tensor interaction matrix
		// Each H-matrix corresponds to one tensor component M[row][col]
		std::cout << "\nBuilding 9 H-matrices (3x3 tensor components)";
		if(config.use_openmp && n_elem > 100)
		{
			std::cout << " in parallel..." << std::endl;
		}
		else
		{
			std::cout << " sequentially..." << std::endl;
		}

		// Thread-safe memory accumulation (each component stores its own memory usage)
		std::vector<size_t> memory_per_component(9);
		memory_used = 0;

		// Parallel construction of 9 H-matrices
		// Use dynamic scheduling for load balancing (different components may have different ranks)
		// Only parallelize for large problems (n_elem > 100) to avoid OpenMP overhead
		#pragma omp parallel for schedule(dynamic) if(config.use_openmp && n_elem > 100)
		for(int idx = 0; idx < 9; idx++)
		{
			int row = idx / 3;
			int col = idx % 3;

			// Thread-safe output (critical section)
			#pragma omp critical
			{
				std::cout << "  Component [" << row << "][" << col << "]... " << std::flush;
			}

			// Prepare kernel data
			KernelData kdata;
			kdata.hmat_ptr = this;
			kdata.tensor_row = row;
			kdata.tensor_col = col;

			// Build H-matrix for this tensor component
			// Each thread builds its own H-matrix independently
			hmat[idx] = hacapk::build_hmatrix(
				points,             // Source points
				points,             // Target points (same for self-interaction)
				KernelFunction,     // Kernel function callback
				&kdata,             // User data
				hacapk_params       // Control parameters
			);

			if(!hmat[idx])
			{
				#pragma omp critical
				{
					std::cerr << "Failed to build H-matrix for component ["
					          << row << "][" << col << "]" << std::endl;
				}
				throw std::runtime_error("Failed to build H-matrix for component ["
				                         + std::to_string(row) + "][" + std::to_string(col) + "]");
			}

			// Store memory usage (thread-safe: each thread writes to different index)
			memory_per_component[idx] = hmat[idx]->memory_usage();

			// Thread-safe output (critical section)
			#pragma omp critical
			{
				std::cout << "rank=" << hmat[idx]->ktmax
				          << ", blocks=" << hmat[idx]->nlf
				          << ", memory=" << (hmat[idx]->memory_usage() / 1024) << " KB" << std::endl;
			}
		}

		// Sum up memory usage from all components (after parallel region)
		for(int idx = 0; idx < 9; idx++)
		{
			memory_used += memory_per_component[idx];
		}

		is_built = true;

		// Calculate compression ratio
		size_t dense_memory = (size_t)n_elem * (size_t)n_elem * 9 * sizeof(double);
		compression_ratio = (double)memory_used / (double)dense_memory;

		auto t_end = std::chrono::high_resolution_clock::now();
		construction_time = std::chrono::duration<double>(t_end - t_start).count();

		std::cout << "\n========================================" << std::endl;
		std::cout << "H-Matrix Construction Complete" << std::endl;
		std::cout << "========================================" << std::endl;
		std::cout << "Construction time: " << construction_time << " s" << std::endl;
		std::cout << "Total memory used: " << (memory_used / 1024 / 1024) << " MB" << std::endl;
		std::cout << "Dense would be: " << (dense_memory / 1024 / 1024) << " MB" << std::endl;
		std::cout << "Compression ratio: " << (compression_ratio * 100) << "%" << std::endl;
		std::cout << "========================================\n" << std::endl;

		return 1; // Success
	}
	catch(const std::exception& e)
	{
		std::cerr << "H-matrix construction failed: " << e.what() << std::endl;
		is_built = false;
		return 0; // Failure
	}
}

//-------------------------------------------------------------------------
// Compute interaction kernel between elements i and j
// Phase 2: Full implementation with symmetry handling
// This matches the logic in radTInteraction::SetupInteractMatrix()
//-------------------------------------------------------------------------

void radTHMatrixInteraction::ComputeInteractionKernel(int i, int j, TMatrix3df& result)
{
	// Setup field computation
	radTFieldKey FieldKeyInteract;
	FieldKeyInteract.B_ = FieldKeyInteract.H_ = FieldKeyInteract.PreRelax_ = 1;
	TVector3d ZeroVect(0., 0., 0.);

	// Get source element (ColNo in SetupInteractMatrix)
	radTg3dRelax* g3dRelaxPtrColNo = elem_ptrs[j];

	// Get observation point at element i (StrNo in SetupInteractMatrix)
	TVector3d InitObsPoiVect = intrct_ptr->MainTransPtrArray[i]->TrPoint(elem_ptrs[i]->ReturnCentrPoint());

	// Use cached symmetry transformations (avoid expensive FillInTransPtrVectForElem call)
	const std::vector<radTrans*>& trans_vect = cached_trans_vect[j];

	TMatrix3d SubMatrix(ZeroVect, ZeroVect, ZeroVect), BufSubMatrix;

	// Loop over all symmetry transformations for element j
	for(unsigned i_trans=0; i_trans < trans_vect.size(); i_trans++)
	{
		radTrans* TransPtr = trans_vect[i_trans];

		TVector3d ObsPoiVect = TransPtr->TrPoint_inv(InitObsPoiVect);

		radTField Field(FieldKeyInteract, intrct_ptr->CompCriterium, ObsPoiVect,
		                ZeroVect, ZeroVect, ZeroVect, ZeroVect, 0.);

		g3dRelaxPtrColNo->B_comp(&Field);

		BufSubMatrix.Str0 = Field.B;
		BufSubMatrix.Str1 = Field.H;
		BufSubMatrix.Str2 = Field.A;

		TransPtr->TrMatrix(BufSubMatrix);
		SubMatrix += BufSubMatrix;
	}

	intrct_ptr->MainTransPtrArray[i]->TrMatrix_inv(SubMatrix);

	// Convert to float matrix
	result = SubMatrix;
}

//-------------------------------------------------------------------------
// Matrix-vector multiplication: H = InteractMatrix * M
// True H-matrix implementation with O(N log N) complexity
//-------------------------------------------------------------------------

void radTHMatrixInteraction::MatVec(const TVector3d* M_in, TVector3d* H_out)
{
	if(!is_built)
	{
		throw std::runtime_error("H-matrix solver: H-matrix not built yet");
	}

	// Convert input magnetization vectors to scalar arrays for HACApK
	// M_in[i] = (Mx, My, Mz) -> M_x[i], M_y[i], M_z[i]
	std::vector<double> M_x(n_elem), M_y(n_elem), M_z(n_elem);
	for(int i = 0; i < n_elem; i++)
	{
		M_x[i] = M_in[i].x;
		M_y[i] = M_in[i].y;
		M_z[i] = M_in[i].z;
	}

	// Allocate result arrays for 9 H-matrix-vector products
	// Result[row][col] = hmat[row*3+col] * M_col
	std::vector<double> result[3][3];
	for(int row = 0; row < 3; row++)
	{
		for(int col = 0; col < 3; col++)
		{
			result[row][col].resize(n_elem, 0.0);
		}
	}

	// Perform 9 H-matrix-vector multiplications using HACApK
	// InteractMatrix is 3x3 tensor: H[row] = sum_col (M[row][col] * M_in[col])
	for(int row = 0; row < 3; row++)
	{
		// H[row].x = M[row][0]*M_x + M[row][1]*M_y + M[row][2]*M_z
		int idx_0 = row * 3 + 0;  // M[row][0]
		int idx_1 = row * 3 + 1;  // M[row][1]
		int idx_2 = row * 3 + 2;  // M[row][2]

		// result[row][0] = hmat[idx_0] * M_x
		hacapk::hmatrix_matvec(*hmat[idx_0], M_x, result[row][0]);

		// result[row][1] = hmat[idx_1] * M_y
		hacapk::hmatrix_matvec(*hmat[idx_1], M_y, result[row][1]);

		// result[row][2] = hmat[idx_2] * M_z
		hacapk::hmatrix_matvec(*hmat[idx_2], M_z, result[row][2]);
	}

	// Combine results into output field vectors
	// H_out[i] = (H_x[i], H_y[i], H_z[i])
	// H_x[i] = result[0][0][i] + result[0][1][i] + result[0][2][i]
	// H_y[i] = result[1][0][i] + result[1][1][i] + result[1][2][i]
	// H_z[i] = result[2][0][i] + result[2][1][i] + result[2][2][i]
	for(int i = 0; i < n_elem; i++)
	{
		H_out[i].x = result[0][0][i] + result[0][1][i] + result[0][2][i];
		H_out[i].y = result[1][0][i] + result[1][1][i] + result[1][2][i];
		H_out[i].z = result[2][0][i] + result[2][1][i] + result[2][2][i];
	}
}

//-------------------------------------------------------------------------
// Print statistics
//-------------------------------------------------------------------------

void radTHMatrixInteraction::PrintStatistics()
{
	if(!is_built)
	{
		std::cout << "H-matrix solver: Not built yet" << std::endl;
		return;
	}

	std::cout << "\n";
	std::cout << "========================================" << std::endl;
	std::cout << "H-Matrix Solver Statistics" << std::endl;
	std::cout << "========================================" << std::endl;
	std::cout << "Number of elements: " << n_elem << std::endl;
	std::cout << "Construction time: " << construction_time << " s" << std::endl;

	size_t dense_memory = (size_t)n_elem * (size_t)n_elem * 9 * sizeof(float);
	std::cout << "Dense matrix memory: " << dense_memory / (1024*1024) << " MB" << std::endl;
	std::cout << "H-matrix memory: " << memory_used / (1024*1024) << " MB" << std::endl;
	std::cout << "Compression ratio: " << compression_ratio * 100 << "%" << std::endl;

	std::cout << "\nConfiguration:" << std::endl;
	std::cout << "  eps = " << config.eps << std::endl;
	std::cout << "  max_rank = " << config.max_rank << std::endl;
	std::cout << "  min_cluster_size = " << config.min_cluster_size << std::endl;
	std::cout << "  use_openmp = " << (config.use_openmp ? "yes" : "no") << std::endl;
	std::cout << "  num_threads = " << config.num_threads << std::endl;
	std::cout << "========================================" << std::endl;
}

//-------------------------------------------------------------------------
// Estimate memory usage
//-------------------------------------------------------------------------

size_t radTHMatrixInteraction::EstimateMemoryUsage() const
{
	return memory_used;
}

//-------------------------------------------------------------------------
// Get compression ratio
//-------------------------------------------------------------------------

double radTHMatrixInteraction::GetCompressionRatio() const
{
	return compression_ratio;
}

//-------------------------------------------------------------------------
// Kernel function wrapper for HACApK
// Extracts a specific tensor component from the 3x3 interaction matrix
//-------------------------------------------------------------------------

double radTHMatrixInteraction::KernelFunction(int i, int j, void* user_data)
{
	// Extract kernel data
	KernelData* kdata = static_cast<KernelData*>(user_data);
	radTHMatrixInteraction* hmat = kdata->hmat_ptr;
	int row = kdata->tensor_row;
	int col = kdata->tensor_col;

	// Get source element j
	radTg3dRelax* elem_j = hmat->elem_ptrs[j];

	// Save original magnetization
	TVector3d originalMagn = elem_j->Magn;

	// Set unit magnetization in the col direction
	// This is critical because B_comp() uses the element's current magnetization
	TVector3d unitMagn(0., 0., 0.);
	switch(col)
	{
	case 0: unitMagn.x = 1.0; break;
	case 1: unitMagn.y = 1.0; break;
	case 2: unitMagn.z = 1.0; break;
	}

	// Result variable (must be declared outside critical section)
	double result = 0.0;

	// Temporarily set unit magnetization
	// CRITICAL SECTION: This modification must be thread-safe
#ifdef _OPENMP
	#pragma omp critical(kernel_magn_access)
	{
#endif
		elem_j->Magn = unitMagn;

		// Compute full 3x3 interaction kernel with unit magnetization
		TMatrix3df kernel;
		hmat->ComputeInteractionKernel(i, j, kernel);

		// Restore original magnetization
		elem_j->Magn = originalMagn;

		// Extract the requested tensor component
		// TMatrix3df is stored as: Str0 = row 0, Str1 = row 1, Str2 = row 2
		switch(row)
		{
		case 0:
			switch(col)
			{
			case 0: result = kernel.Str0.x; break;
			case 1: result = kernel.Str0.y; break;
			case 2: result = kernel.Str0.z; break;
			}
			break;
		case 1:
			switch(col)
			{
			case 0: result = kernel.Str1.x; break;
			case 1: result = kernel.Str1.y; break;
			case 2: result = kernel.Str1.z; break;
			}
			break;
		case 2:
			switch(col)
			{
			case 0: result = kernel.Str2.x; break;
			case 1: result = kernel.Str2.y; break;
			case 2: result = kernel.Str2.z; break;
			}
			break;
		}

#ifdef _OPENMP
	}  // End critical section
#endif

	return result;
}
