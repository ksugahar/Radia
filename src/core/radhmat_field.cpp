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
#include "rad_rectangular_block.h"
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
	, hmatrix_data_x(nullptr)
	, hmatrix_data_y(nullptr)
	, hmatrix_data_z(nullptr)
	, source_cluster_tree(nullptr)
	, target_cluster_tree(nullptr)
	, current_component(0)
	, geometry_hash(0)
	, target_hash(0)
	, memory_usage(0)
	, build_time(0.0)
	, last_eval_time(0.0)
	, num_evaluations(0)
{
	std::cout << "[HMatrix Field] Evaluator created (3D vector field)" << std::endl;
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
	// Free HACApK data structures (3 H-matrices for vector field)
	if(hmatrix_data_x) {
		delete static_cast<hacapk::HMatrix*>(hmatrix_data_x);
		hmatrix_data_x = nullptr;
	}
	if(hmatrix_data_y) {
		delete static_cast<hacapk::HMatrix*>(hmatrix_data_y);
		hmatrix_data_y = nullptr;
	}
	if(hmatrix_data_z) {
		delete static_cast<hacapk::HMatrix*>(hmatrix_data_z);
		hmatrix_data_z = nullptr;
	}

	if(source_cluster_tree) {
		delete static_cast<std::shared_ptr<hacapk::Cluster>*>(source_cluster_tree);
		source_cluster_tree = nullptr;
	}

	if(target_cluster_tree) {
		delete static_cast<std::shared_ptr<hacapk::Cluster>*>(target_cluster_tree);
		target_cluster_tree = nullptr;
	}

	source_positions.clear();
	source_moments.clear();
	target_points.clear();

	is_built = false;
	num_sources = 0;
	geometry_hash = 0;
	target_hash = 0;
	memory_usage = 0;
	num_evaluations = 0;

	std::cout << "[HMatrix Field] Cleared (3 H-matrices freed)" << std::endl;
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

size_t ComputeTargetHash(const std::vector<TVector3d>& obs_points)
{
	// Hash observation points for cache validation
	std::hash<double> hasher;
	size_t hash = 0;

	// Hash number of points
	hash ^= std::hash<int>{}(obs_points.size());

	// Hash all points (or first 100 for speed)
	int count = 0;
	int max_points = std::min(100, (int)obs_points.size());
	for(int i = 0; i < max_points; i++) {
		const auto& pt = obs_points[i];
		hash ^= hasher(pt.x) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		hash ^= hasher(pt.y) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		hash ^= hasher(pt.z) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
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

	int num_top_level = static_cast<int>(source_group->GroupMapOfHandlers.size());

	if(num_top_level < 1) {
		std::cerr << "[HMatrix Field] No source elements" << std::endl;
		return 0;
	}

	std::cout << "[HMatrix Field] Extracting geometry from " << num_top_level << " top-level elements" << std::endl;

	// Clear existing data
	source_positions.clear();
	source_moments.clear();

	// Recursively extract all sub-elements
	for(auto& elem_pair : source_group->GroupMapOfHandlers) {
		radTg3d* g3d = (radTg3d*)(elem_pair.second.rep);
		if(!g3d) continue;

		// Recursively extract leaf elements (handles subdivisions)
		ExtractLeafElements(g3d, 0);
	}

	num_sources = static_cast<int>(source_positions.size() / 3);

	std::cout << "[HMatrix Field] Extracted " << num_sources << " source points (including sub-elements)" << std::endl;

	if(num_sources < 1) {
		std::cerr << "[HMatrix Field] No valid source points extracted" << std::endl;
		return 0;
	}

	return 1;
}

//-------------------------------------------------------------------------

// Helper: recursively extract leaf elements (sub-elements)
void radTHMatrixFieldEvaluator::ExtractLeafElements(radTg3d* g3d, int depth)
{
	if(!g3d) return;

	// Try to cast to radTGroup (container with sub-elements)
	radTGroup* group = dynamic_cast<radTGroup*>(g3d);

	if(group && group->GroupMapOfHandlers.size() > 0) {
		// Container - recursively extract sub-elements
		for(auto& elem_pair : group->GroupMapOfHandlers) {
			radTg3d* sub_elem = (radTg3d*)(elem_pair.second.rep);
			ExtractLeafElements(sub_elem, depth + 1);
		}
	}
	else {
		// Leaf element - extract it
		radTg3dRelax* relaxable = dynamic_cast<radTg3dRelax*>(g3d);
		if(!relaxable) return;

		TVector3d M = relaxable->Magn;
		TVector3d center = g3d->CentrPoint;
		double volume = relaxable->Volume();
		double vol_m3 = volume * 1e-9;
		TVector3d moment = M * vol_m3;

		source_positions.push_back(center.x);
		source_positions.push_back(center.y);
		source_positions.push_back(center.z);
		source_moments.push_back(moment.x);
		source_moments.push_back(moment.y);
		source_moments.push_back(moment.z);
	}
}

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

	// Build source cluster tree
	std::vector<hacapk::Point3D> source_hacapk;
	source_hacapk.reserve(num_sources);

	for(int i = 0; i < num_sources; i++) {
		source_hacapk.emplace_back(
			source_positions[i*3 + 0],
			source_positions[i*3 + 1],
			source_positions[i*3 + 2]
		);
	}

	hacapk::ControlParams params;
	params.leaf_size = config.min_cluster_size;
	params.eta = 2.0;  // Admissibility parameter

	std::vector<int> indices(num_sources);
	for(int i = 0; i < num_sources; i++) indices[i] = i;

	auto cluster = hacapk::generate_cluster(
		source_hacapk,
		indices,
		0,              // start
		num_sources,    // size
		0,              // depth
		params
	);

	if(!cluster) {
		std::cerr << "[HMatrix Field] Failed to build source cluster tree" << std::endl;
		return 0;
	}

	// Store cluster tree (convert shared_ptr to void*)
	source_cluster_tree = new std::shared_ptr<hacapk::Cluster>(cluster);

	std::cout << "[HMatrix Field] Source cluster tree built successfully" << std::endl;

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

	// Direct summation: H(r) = sum [3(m*r_hatr_hat- m] / (4pi|r|^3)
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

			// m * r
			double m_dot_r = mx*rx_m + my*ry_m + mz*rz_m;

			// H = (1/4pi) * [3(m*r_hatr_hat- m] / r^3
			// = (1/4pi) * [3(m*r)r/r^5 - m/r^3]
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
		// B = mu0 * H
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
	int M = static_cast<int>(obs_points.size());

	std::cout << "[HMatrix Field] Building target cluster tree for " << M << " points" << std::endl;

	// Convert observation points to HACApK format
	target_points.clear();
	target_points.reserve(M * 3);

	std::vector<hacapk::Point3D> target_hacapk;
	target_hacapk.reserve(M);

	for(int i = 0; i < M; i++) {
		const auto& pt = obs_points[i];
		target_points.push_back(pt.x);
		target_points.push_back(pt.y);
		target_points.push_back(pt.z);

		target_hacapk.emplace_back(pt.x, pt.y, pt.z);
	}

	// Build cluster tree using HACApK
	hacapk::ControlParams params;
	params.leaf_size = config.min_cluster_size;
	params.eta = 2.0;  // Admissibility parameter

	std::vector<int> indices(M);
	for(int i = 0; i < M; i++) indices[i] = i;

	auto cluster = hacapk::generate_cluster(
		target_hacapk,
		indices,
		0,      // start
		M,      // size
		0,      // depth
		params
	);

	if(!cluster) {
		std::cerr << "[HMatrix Field] Failed to build target cluster tree" << std::endl;
		return 0;
	}

	// Store cluster tree (convert shared_ptr to void*)
	target_cluster_tree = new std::shared_ptr<hacapk::Cluster>(cluster);

	std::cout << "[HMatrix Field] Target cluster tree built successfully" << std::endl;

	return 1;
}

//-------------------------------------------------------------------------

int radTHMatrixFieldEvaluator::BuildFieldHMatrix()
{
	std::cout << "[HMatrix Field] Building 3 H-matrices for vector field (Hx, Hy, Hz)" << std::endl;

	if(!source_cluster_tree || !target_cluster_tree) {
		std::cerr << "[HMatrix Field] Cluster trees not built" << std::endl;
		return 0;
	}

	// Prepare HACApK points
	int M = static_cast<int>(target_points.size() / 3);
	int N = num_sources;

	std::vector<hacapk::Point3D> source_hacapk(N);
	std::vector<hacapk::Point3D> target_hacapk(M);

	for(int i = 0; i < N; i++) {
		source_hacapk[i] = hacapk::Point3D(
			source_positions[i*3 + 0],
			source_positions[i*3 + 1],
			source_positions[i*3 + 2]
		);
	}

	for(int i = 0; i < M; i++) {
		target_hacapk[i] = hacapk::Point3D(
			target_points[i*3 + 0],
			target_points[i*3 + 1],
			target_points[i*3 + 2]
		);
	}

	// Setup kernel parameters
	hacapk::ControlParams params;
	params.eps_aca = config.eps;
	params.leaf_size = config.min_cluster_size;
	params.eta = 2.0;  // Admissibility parameter

	// Kernel function wrapper
	auto kernel_func = [](int i, int j, void* data) -> double {
		return radTHMatrixFieldEvaluator::FieldKernel(i, j, data);
	};

	auto start_time = std::chrono::high_resolution_clock::now();

	memory_usage = 0;

	// Build 3 H-matrices, one for each component
	const char* comp_names[] = {"Hx", "Hy", "Hz"};
	void** hmatrix_ptrs[] = {&hmatrix_data_x, &hmatrix_data_y, &hmatrix_data_z};

	for(int comp = 0; comp < 3; comp++) {
		std::cout << "[HMatrix Field] Building H-matrix for " << comp_names[comp] << "..." << std::endl;

		// Set current component for kernel
		current_component = comp;

		// Build H-matrix
		auto hmatrix = hacapk::build_hmatrix(
			source_hacapk,
			target_hacapk,
			kernel_func,
			this,  // kernel_data = this evaluator
			params
		);

		if(!hmatrix) {
			std::cerr << "[HMatrix Field] Failed to build H-matrix for " << comp_names[comp] << std::endl;
			return 0;
		}

		// Store H-matrix
		*hmatrix_ptrs[comp] = hmatrix.release();

		// Update memory usage
		auto* hmat = static_cast<hacapk::HMatrix*>(*hmatrix_ptrs[comp]);
		size_t comp_memory = 0;
		for(const auto& block : hmat->blocks) {
			comp_memory += block.memory_usage();
		}
		memory_usage += comp_memory;

		std::cout << "[HMatrix Field]   " << comp_names[comp] << ": "
		          << hmat->blocks.size() << " blocks, "
		          << comp_memory / 1024.0 / 1024.0 << " MB" << std::endl;
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	double construction_time = std::chrono::duration<double>(end_time - start_time).count();

	std::cout << "[HMatrix Field] All 3 H-matrices built in " << construction_time << " seconds" << std::endl;
	std::cout << "[HMatrix Field] Matrix size: " << M << " x " << N << std::endl;
	std::cout << "[HMatrix Field] Total memory: " << memory_usage / 1024.0 / 1024.0 << " MB" << std::endl;

	return 1;
}

//-------------------------------------------------------------------------

double radTHMatrixFieldEvaluator::FieldKernel(int i, int j, void* kernel_data)
{
	// Magnetic dipole kernel for H-field calculation
	// H(r) = (1/4pi) * [3(m*r_hatr_hat- m] / |r|^3
	//
	// For H-matrix, we return scalar kernel value for each component
	// The actual vector field is assembled from 9 H-matrices (3x3 tensor)

	auto* eval = static_cast<radTHMatrixFieldEvaluator*>(kernel_data);

	// Target (observation) point position (mm)
	double tx = eval->target_points[i*3 + 0];
	double ty = eval->target_points[i*3 + 1];
	double tz = eval->target_points[i*3 + 2];

	// Source position (mm)
	double sx = eval->source_positions[j*3 + 0];
	double sy = eval->source_positions[j*3 + 1];
	double sz = eval->source_positions[j*3 + 2];

	// Magnetic moment (A*m^2)
	double mx = eval->source_moments[j*3 + 0];
	double my = eval->source_moments[j*3 + 1];
	double mz = eval->source_moments[j*3 + 2];

	// Vector from source to observation point (mm)
	double rx = tx - sx;
	double ry = ty - sy;
	double rz = tz - sz;

	// Distance (convert to meters for SI units)
	double r_mm = std::sqrt(rx*rx + ry*ry + rz*rz);
	if(r_mm < 1e-6) return 0.0;  // Singularity avoidance

	double r_m = r_mm * 1e-3;  // mm to m
	double r3 = r_m * r_m * r_m;
	double r5 = r3 * r_m * r_m;

	// Unit vector (dimensionless)
	double rx_m = rx * 1e-3;
	double ry_m = ry * 1e-3;
	double rz_m = rz * 1e-3;

	// m * r (A*m^3)
	double m_dot_r = mx*rx_m + my*ry_m + mz*rz_m;

	// Compute H = (1/4pi) * [3(m*r_hat)r_hat - m] / r^3
	const double PI = 3.14159265358979;
	const double C = 1.0 / (4.0 * PI);

	// H-field components (A/m)
	double Hx = C * (3.0 * m_dot_r * rx_m / r5 - mx / r3);
	double Hy = C * (3.0 * m_dot_r * ry_m / r5 - my / r3);
	double Hz = C * (3.0 * m_dot_r * rz_m / r5 - mz / r3);

	// Return component specified by current_component (0=x, 1=y, 2=z)
	switch(eval->current_component) {
		case 0: return Hx;
		case 1: return Hy;
		case 2: return Hz;
		default: return 0.0;
	}
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

	// Compute hash of observation points
	size_t current_target_hash = ComputeTargetHash(obs_points);

	// Check if we can reuse existing H-matrix
	bool need_rebuild = (current_target_hash != target_hash);

	if(need_rebuild) {
		std::cout << "[HMatrix Field] Observation points changed, rebuilding H-matrix..." << std::endl;

		// Step 1: Build target cluster tree for observation points
		if(!BuildTargetClusterTree(obs_points)) {
			std::cerr << "[HMatrix Field] Failed to build target cluster tree, falling back to direct" << std::endl;
			return EvaluateFieldDirect(obs_points, field_out, field_type);
		}

		// Step 2: Build H-matrix
		if(!BuildFieldHMatrix()) {
			std::cerr << "[HMatrix Field] Failed to build H-matrix, falling back to direct" << std::endl;
			return EvaluateFieldDirect(obs_points, field_out, field_type);
		}

		// Update cache
		target_hash = current_target_hash;
	}
	else {
		std::cout << "[HMatrix Field] Reusing cached H-matrix (same observation points)" << std::endl;
	}

	// Step 3: H-matrix vector multiplication for all 3 components
	std::cout << "[HMatrix Field] Computing field with 3 H-matrices" << std::endl;

	// Input vector: all ones (FieldKernel already includes magnetic moments)
	std::vector<double> ones(num_sources, 1.0);

	// Output vectors
	std::vector<double> Hx(M), Hy(M), Hz(M);

	// H-matrix vector multiplication for each component
	auto* hmat_x = static_cast<hacapk::HMatrix*>(hmatrix_data_x);
	auto* hmat_y = static_cast<hacapk::HMatrix*>(hmatrix_data_y);
	auto* hmat_z = static_cast<hacapk::HMatrix*>(hmatrix_data_z);

	hacapk::hmatrix_matvec(*hmat_x, ones, Hx);
	hacapk::hmatrix_matvec(*hmat_y, ones, Hy);
	hacapk::hmatrix_matvec(*hmat_z, ones, Hz);

	// Assemble field output - full 3D vector field
	field_out.resize(M, TVector3d(0, 0, 0));
	for(int i = 0; i < M; i++) {
		field_out[i].x = Hx[i];
		field_out[i].y = Hy[i];
		field_out[i].z = Hz[i];
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	last_eval_time = std::chrono::duration<double>(end_time - start_time).count();

	// Convert to requested field type
	const double MU0 = 4.0 * 3.14159265358979 * 1e-7;  // H/m
	if(field_type == 'b') {
		// B = mu0 * H
		for(auto& field : field_out) {
			field.x *= MU0;
			field.y *= MU0;
			field.z *= MU0;
		}
	}

	std::cout << "[HMatrix Field] H-matrix evaluation completed in " << last_eval_time << " seconds" << std::endl;
	std::cout << "[HMatrix Field] Speedup vs direct: "
	          << (M * num_sources * 1e-9) / last_eval_time << " GFlop/s equivalent" << std::endl;
	std::cout << "[HMatrix Field] Reused " << num_evaluations << " times" << std::endl;

	return 1;
}

//-------------------------------------------------------------------------
