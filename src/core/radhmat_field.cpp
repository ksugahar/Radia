/*-------------------------------------------------------------------------
*
* File name:      radhmat_field.cpp
*
* Project:        RADIA
*
* Description:    H-matrix field evaluator implementation
*                 Fast batch field evaluation using HACApK
*
* Author(s):      Radia Development Team
*
* First release:  2025
*
-------------------------------------------------------------------------*/

#include "radhmat.h"
#include "radgroup.h"
#include "radexcep.h"
#include "radrec.h"
#include "radsend.h"
#include "radg3d.h"
#include <iostream>
#include <chrono>
#include <cmath>
#include <functional>

#ifdef _OPENMP
#include <omp.h>
#endif

// HACApK includes
#include "hacapk.hpp"

//-------------------------------------------------------------------------
// radTHMatrixFieldEvaluator Implementation
//-------------------------------------------------------------------------

radTHMatrixFieldEvaluator::radTHMatrixFieldEvaluator(const radTHMatrixConfig& cfg)
	: config(cfg)
	, is_built(false)
	, num_sources(0)
	, hmatrix_data(nullptr)
	, source_cluster_tree(nullptr)
	, target_cluster_tree(nullptr)
	, geometry_hash(0)
	, memory_usage(0)
	, build_time(0.0)
	, last_eval_time(0.0)
	, num_evaluations(0)
{
	std::cout << "[HMatrix Field] Evaluator created" << std::endl;
}

//-------------------------------------------------------------------------

radTHMatrixFieldEvaluator::~radTHMatrixFieldEvaluator()
{
	Clear();
	std::cout << "[HMatrix Field] Evaluator destroyed" << std::endl;
}

//-------------------------------------------------------------------------

void radTHMatrixFieldEvaluator::Clear()
{
	// Free HACApK data structures
	if(hmatrix_data) {
		// TODO: Delete HACApK H-matrix
		hmatrix_data = nullptr;
	}

	if(source_cluster_tree) {
		// TODO: Delete source cluster tree
		source_cluster_tree = nullptr;
	}

	if(target_cluster_tree) {
		// TODO: Delete target cluster tree
		target_cluster_tree = nullptr;
	}

	source_positions.clear();
	source_moments.clear();

	is_built = false;
	num_sources = 0;
	geometry_hash = 0;
	memory_usage = 0;
	num_evaluations = 0;

	std::cout << "[HMatrix Field] Cleared" << std::endl;
}

//-------------------------------------------------------------------------

size_t radTHMatrixFieldEvaluator::ComputeGeometryHash(radTGroup* source_group)
{
	if(!source_group) return 0;

	// Simple hash based on number of elements and first few positions
	std::hash<double> hasher;
	size_t hash = 0;

	// Hash number of elements
	hash ^= std::hash<int>{}(source_group->GroupMapOfHandlers.size());

	// Hash first few element positions (for quick check)
	int count = 0;
	for(auto& elem_pair : source_group->GroupMapOfHandlers) {
		radTg3d* g3d = (radTg3d*)(elem_pair.second.rep);
		if(!g3d) continue;

		TVector3d center = g3d->CentrPoint;
		hash ^= hasher(center.x) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		hash ^= hasher(center.y) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		hash ^= hasher(center.z) + 0x9e3779b9 + (hash << 6) + (hash >> 2);

		if(++count >= 10) break;  // Only hash first 10 elements for speed
	}

	return hash;
}

//-------------------------------------------------------------------------

bool radTHMatrixFieldEvaluator::IsValid(radTGroup* source_group)
{
	if(!is_built) return false;
	if(!source_group) return false;

	size_t current_hash = ComputeGeometryHash(source_group);
	return (current_hash == geometry_hash);
}

//-------------------------------------------------------------------------

int radTHMatrixFieldEvaluator::ExtractSourceGeometry(radTGroup* source_group)
{
	if(!source_group) {
		std::cerr << "[HMatrix Field] NULL source group" << std::endl;
		return 0;
	}

	num_sources = static_cast<int>(source_group->GroupMapOfHandlers.size());

	if(num_sources < 1) {
		std::cerr << "[HMatrix Field] No source elements" << std::endl;
		return 0;
	}

	std::cout << "[HMatrix Field] Extracting geometry from " << num_sources << " elements" << std::endl;

	// Reserve space
	source_positions.reserve(num_sources * 3);
	source_moments.reserve(num_sources * 3);

	// Extract center points and magnetic moments
	int extracted_count = 0;

	for(auto& elem_pair : source_group->GroupMapOfHandlers) {
		radTg3d* g3d = (radTg3d*)(elem_pair.second.rep);
		if(!g3d) continue;

		// Try to cast to radTg3dRelax (magnetized volume)
		radTg3dRelax* relaxable = dynamic_cast<radTg3dRelax*>(g3d);
		if(!relaxable) continue;  // Skip non-relaxable elements

		// Get magnetization
		TVector3d M = relaxable->Magn;

		// Center point (mm)
		TVector3d center = g3d->CentrPoint;
		source_positions.push_back(center.x);
		source_positions.push_back(center.y);
		source_positions.push_back(center.z);

		// Get volume (mm^3)
		double volume = relaxable->Volume();

		// Magnetic moment = M * V (A*m^2)
		// Note: Volume() returns mm^3, M is A/m
		// Need to convert: volume_m3 = volume_mm3 * 1e-9
		double vol_m3 = volume * 1e-9;  // mm^3 to m^3
		TVector3d moment = M * vol_m3;

		source_moments.push_back(moment.x);
		source_moments.push_back(moment.y);
		source_moments.push_back(moment.z);

		extracted_count++;
	}

	std::cout << "[HMatrix Field] Extracted " << source_positions.size() / 3 << " source points" << std::endl;

	return 1;
}

//-------------------------------------------------------------------------

int radTHMatrixFieldEvaluator::Build(radTGroup* source_group)
{
	if(!source_group) {
		std::cerr << "[HMatrix Field] NULL source group" << std::endl;
		return 0;
	}

	auto start_time = std::chrono::high_resolution_clock::now();

	// Clear previous data
	Clear();

	// Extract geometry
	if(!ExtractSourceGeometry(source_group)) {
		std::cerr << "[HMatrix Field] Failed to extract geometry" << std::endl;
		return 0;
	}

	// Check problem size
	if(num_sources < config.min_cluster_size) {
		std::cout << "[HMatrix Field] Problem too small (N=" << num_sources
		          << "), use direct calculation instead" << std::endl;
		return 0;  // Will fallback to direct calculation
	}

	// Compute geometry hash for caching
	geometry_hash = ComputeGeometryHash(source_group);

	// TODO: Build source cluster tree
	// source_cluster_tree = BuildClusterTree(source_positions, config.min_cluster_size);

	is_built = true;

	auto end_time = std::chrono::high_resolution_clock::now();
	build_time = std::chrono::duration<double>(end_time - start_time).count();

	std::cout << "[HMatrix Field] Build completed in " << build_time << " seconds" << std::endl;
	std::cout << "[HMatrix Field] Memory usage: " << memory_usage / 1024.0 / 1024.0 << " MB" << std::endl;

	return 1;
}

//-------------------------------------------------------------------------

int radTHMatrixFieldEvaluator::EvaluateFieldDirect(
	const std::vector<TVector3d>& obs_points,
	std::vector<TVector3d>& field_out,
	char field_type)
{
	int M = static_cast<int>(obs_points.size());
	int N = num_sources;

	if(M == 0 || N == 0) {
		std::cerr << "[HMatrix Field] No points to evaluate" << std::endl;
		return 0;
	}

	std::cout << "[HMatrix Field] Direct calculation: M=" << M << ", N=" << N << std::endl;

	field_out.resize(M, TVector3d(0, 0, 0));

	const double PI = 3.14159265358979;
	const double MU0 = 4.0 * PI * 1e-7;  // H/m
	const double C = 1.0 / (4.0 * PI);   // For H-field

	auto start_time = std::chrono::high_resolution_clock::now();

	// Direct summation: H(r) = Σ [3(m·r̂)r̂ - m] / (4π|r|³)
	#pragma omp parallel for if(M > 100)
	for(int i = 0; i < M; i++) {
		TVector3d H_total(0, 0, 0);
		const TVector3d& obs = obs_points[i];

		for(int j = 0; j < N; j++) {
			// Source position (mm)
			double sx = source_positions[j*3 + 0];
			double sy = source_positions[j*3 + 1];
			double sz = source_positions[j*3 + 2];

			// Magnetic moment (A*m^2)
			double mx = source_moments[j*3 + 0];
			double my = source_moments[j*3 + 1];
			double mz = source_moments[j*3 + 2];

			// Vector from source to observation point (mm)
			double rx = obs.x - sx;
			double ry = obs.y - sy;
			double rz = obs.z - sz;

			// Distance (convert to meters)
			double r_mm = std::sqrt(rx*rx + ry*ry + rz*rz);
			if(r_mm < 1e-6) continue;  // Skip if too close

			double r_m = r_mm * 1e-3;  // mm to m
			double r3 = r_m * r_m * r_m;
			double r5 = r3 * r_m * r_m;

			// Unit vector
			double rx_m = rx * 1e-3;
			double ry_m = ry * 1e-3;
			double rz_m = rz * 1e-3;

			// m · r
			double m_dot_r = mx*rx_m + my*ry_m + mz*rz_m;

			// H = (1/4π) * [3(m·r̂)r̂ - m] / r³
			// = (1/4π) * [3(m·r)r/r⁵ - m/r³]
			double coeff1 = C * 3.0 * m_dot_r / r5;
			double coeff2 = C / r3;

			H_total.x += coeff1 * rx_m - coeff2 * mx;
			H_total.y += coeff1 * ry_m - coeff2 * my;
			H_total.z += coeff1 * rz_m - coeff2 * mz;
		}

		field_out[i] = H_total;
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	last_eval_time = std::chrono::duration<double>(end_time - start_time).count();

	std::cout << "[HMatrix Field] Direct evaluation completed in " << last_eval_time << " seconds" << std::endl;

	// Convert to requested field type
	if(field_type == 'b') {
		// B = μ₀ * H
		for(auto& field : field_out) {
			field.x *= MU0;
			field.y *= MU0;
			field.z *= MU0;
		}
	}

	return 1;
}

//-------------------------------------------------------------------------

int radTHMatrixFieldEvaluator::BuildTargetClusterTree(const std::vector<TVector3d>& obs_points)
{
	// TODO: Implement target cluster tree construction using HACApK
	std::cout << "[HMatrix Field] Building target cluster tree for " << obs_points.size() << " points" << std::endl;

	// For now, return success (will use direct calculation)
	return 1;
}

//-------------------------------------------------------------------------

int radTHMatrixFieldEvaluator::BuildFieldHMatrix()
{
	// TODO: Implement H-matrix construction for field evaluation
	std::cout << "[HMatrix Field] Building field H-matrix" << std::endl;

	// For now, return 0 to trigger direct calculation
	return 0;
}

//-------------------------------------------------------------------------

double radTHMatrixFieldEvaluator::FieldKernel(int i, int j, void* kernel_data)
{
	// TODO: Implement magnetic field kernel
	// For now, return 0
	return 0.0;
}

//-------------------------------------------------------------------------

int radTHMatrixFieldEvaluator::EvaluateField(
	const std::vector<TVector3d>& obs_points,
	std::vector<TVector3d>& field_out,
	char field_type)
{
	int M = static_cast<int>(obs_points.size());

	if(M == 0) {
		std::cerr << "[HMatrix Field] No observation points" << std::endl;
		return 0;
	}

	if(!is_built) {
		std::cerr << "[HMatrix Field] H-matrix not built, using direct calculation" << std::endl;
		return EvaluateFieldDirect(obs_points, field_out, field_type);
	}

	// Check if H-matrix acceleration is worthwhile
	if(M < 100 || num_sources < 100) {
		std::cout << "[HMatrix Field] Problem too small for H-matrix, using direct" << std::endl;
		return EvaluateFieldDirect(obs_points, field_out, field_type);
	}

	num_evaluations++;

	auto start_time = std::chrono::high_resolution_clock::now();

	// TODO: Implement H-matrix accelerated evaluation
	// For now, fallback to direct calculation
	std::cout << "[HMatrix Field] H-matrix evaluation not yet implemented, using direct" << std::endl;
	return EvaluateFieldDirect(obs_points, field_out, field_type);

	// Future implementation:
	// 1. BuildTargetClusterTree(obs_points)
	// 2. BuildFieldHMatrix()
	// 3. H-matrix vector multiplication

	auto end_time = std::chrono::high_resolution_clock::now();
	last_eval_time = std::chrono::duration<double>(end_time - start_time).count();

	std::cout << "[HMatrix Field] Evaluation completed in " << last_eval_time << " seconds" << std::endl;
	std::cout << "[HMatrix Field] Reused " << num_evaluations << " times" << std::endl;

	return 1;
}

//-------------------------------------------------------------------------
