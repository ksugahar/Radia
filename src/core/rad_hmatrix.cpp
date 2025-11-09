/*-------------------------------------------------------------------------
*
* File name:      radhmat.cpp
*
* Project:        RADIA
*
* Description:    H-matrix (Hierarchical Matrix) based field source
*                 using HACApK for fast field computation
*
* Author(s):      Radia Development Team
*
* First release:  2025
*
-------------------------------------------------------------------------*/

#include "rad_hmatrix.h"
#include "rad_group.h"
#include "rad_exception.h"
#include "rad_rectangular_block.h"
#include <iostream>
#include <chrono>
#include <cmath>

// HACApK includes
#include "hacapk.hpp"

// Radia application includes
#include "rad_application.h"
#include "rad_serialization.h"

//-------------------------------------------------------------------------
// Helper structures for kernel function
//-------------------------------------------------------------------------

struct KernelData {
	std::vector<hacapk::Point3D>* source_points;
	std::vector<hacapk::Point3D>* target_points;
	std::vector<TVector3d>* magnetic_moments;  // Magnetic moments (A*m^2)
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

radTHMatrixFieldSource::radTHMatrixFieldSource(
	radTGroup* group,
	const radTHMatrixConfig& cfg)
	: config(cfg)
	, is_built(false)
	, hmatrix_data(nullptr)
	, cluster_tree_data(nullptr)
	, num_elements(0)
	, hmatrix_memory(0)
	, construction_time(0.0)
{
	if(!group) {
		throw radTException("radTHMatrixFieldSource: NULL group pointer");
	}

	// Copy source elements from group
	source_elements = group->GroupMapOfHandlers;
	num_elements = static_cast<int>(source_elements.size());

	std::cout << "[HMatrix] Created field source with " << num_elements << " elements" << std::endl;

	// Extract geometry for future H-matrix construction
	if(!ExtractGeometry()) {
		std::cerr << "[HMatrix] Warning: Geometry extraction failed" << std::endl;
	}
}

//-------------------------------------------------------------------------

radTHMatrixFieldSource::radTHMatrixFieldSource(
	CAuxBinStrVect& inStr,
	map<int, int>& mKeysOldNew,
	radTmhg& gMapOfHandlers)
	: config()
	, is_built(false)
	, hmatrix_data(nullptr)
	, cluster_tree_data(nullptr)
	, num_elements(0)
	, hmatrix_memory(0)
	, construction_time(0.0)
{
	// Binary deserialization not fully supported yet
	DumpBinParse_g3d(inStr, mKeysOldNew, gMapOfHandlers);

	// Note: H-matrix needs to be rebuilt after deserialization
	std::cerr << "[HMatrix] Warning: Deserialized H-matrix needs rebuild" << std::endl;
}

//-------------------------------------------------------------------------

radTHMatrixFieldSource::radTHMatrixFieldSource(const radTHMatrixFieldSource& src)
	: radTg3d(src)
	, source_elements(src.source_elements)
	, config(src.config)
	, is_built(false)  // Do not copy H-matrix data (rebuild required)
	, hmatrix_data(nullptr)
	, cluster_tree_data(nullptr)
	, element_positions(src.element_positions)
	, element_moments(src.element_moments)
	, num_elements(src.num_elements)
	, hmatrix_memory(0)
	, construction_time(0.0)
{
	// Note: H-matrix needs to be rebuilt for copied object
}

//-------------------------------------------------------------------------

radTHMatrixFieldSource::~radTHMatrixFieldSource()
{
	// Clean up H-matrix data structures
	if(hmatrix_data) {
		delete static_cast<hacapk::HMatrix*>(hmatrix_data);
		hmatrix_data = nullptr;
	}

	// Cluster tree is owned by shared_ptr, will be cleaned up automatically
	cluster_tree_data = nullptr;
}

//-------------------------------------------------------------------------

int radTHMatrixFieldSource::ExtractGeometry()
{
	if(source_elements.empty()) {
		std::cerr << "[HMatrix] No source elements to extract" << std::endl;
		return 0;
	}

	// Reserve space
	element_positions.reserve(num_elements * 3);
	element_moments.reserve(num_elements * 3);

	int extracted_count = 0;

	// Extract center positions and magnetic properties
	for(const auto& pair : source_elements) {
		radTg3d* elem = static_cast<radTg3d*>(pair.second.rep);
		if(!elem) continue;

		// Get element center and magnetization
		TVector3d center(0, 0, 0);
		TVector3d M(0, 0, 0);
		double volume = 0.0;

		// Try to cast to known types and extract geometry
		// Type 1: radTg3dRelax (magnetized volume)
		radTg3dRelax* relaxable = dynamic_cast<radTg3dRelax*>(elem);
		if(relaxable) {
			// Get magnetization
			M = relaxable->Magn;

			// Get actual center from radTg3d::CentrPoint
			center = elem->CentrPoint;

			// Get volume
			volume = relaxable->Volume();

			// Store data
			element_positions.push_back(center.x);
			element_positions.push_back(center.y);
			element_positions.push_back(center.z);

			// Magnetic moment = M * Volume (A*m^2)
			TVector3d moment = M * volume;
			element_moments.push_back(moment.x);
			element_moments.push_back(moment.y);
			element_moments.push_back(moment.z);

			extracted_count++;
		}
	}

	if(extracted_count == 0) {
		std::cerr << "[HMatrix] Warning: No geometry could be extracted" << std::endl;
		return 0;
	}

	std::cout << "[HMatrix] Extracted geometry from " << extracted_count << " elements" << std::endl;

	// Print statistics for verification
	if(extracted_count > 0) {
		// Calculate bounding box
		double xmin = element_positions[0], xmax = element_positions[0];
		double ymin = element_positions[1], ymax = element_positions[1];
		double zmin = element_positions[2], zmax = element_positions[2];

		for(int i = 0; i < extracted_count; i++) {
			double x = element_positions[3*i];
			double y = element_positions[3*i+1];
			double z = element_positions[3*i+2];

			if(x < xmin) xmin = x; if(x > xmax) xmax = x;
			if(y < ymin) ymin = y; if(y > ymax) ymax = y;
			if(z < zmin) zmin = z; if(z > zmax) zmax = z;
		}

		std::cout << "[HMatrix] Bounding box: "
		          << "[" << xmin << "," << xmax << "] x "
		          << "[" << ymin << "," << ymax << "] x "
		          << "[" << zmin << "," << zmax << "] mm" << std::endl;

		// Calculate average magnetization magnitude
		double total_M = 0.0;
		for(int i = 0; i < extracted_count; i++) {
			double mx = element_moments[3*i];
			double my = element_moments[3*i+1];
			double mz = element_moments[3*i+2];
			double M_mag = std::sqrt(mx*mx + my*my + mz*mz);
			total_M += M_mag;
		}
		double avg_M = total_M / extracted_count;

		std::cout << "[HMatrix] Average magnetic moment: " << avg_M << " A*m^2" << std::endl;
	}

	return 1;
}

//-------------------------------------------------------------------------

int radTHMatrixFieldSource::BuildHMatrix()
{
	if(is_built) {
		std::cout << "[HMatrix] Already built, skipping" << std::endl;
		return 1;
	}

	if(num_elements == 0) {
		std::cerr << "[HMatrix] Error: No source elements" << std::endl;
		return 0;
	}

	std::cout << "[HMatrix] Building H-matrix for " << num_elements << " elements..." << std::endl;
	std::cout << "[HMatrix] Configuration:" << std::endl;
	std::cout << "  eps = " << config.eps << std::endl;
	std::cout << "  max_rank = " << config.max_rank << std::endl;
	std::cout << "  min_cluster_size = " << config.min_cluster_size << std::endl;
	std::cout << "  use_openmp = " << (config.use_openmp ? "true" : "false") << std::endl;

	auto start_time = std::chrono::high_resolution_clock::now();

	try {
		// Convert element positions to hacapk::Point3D
		std::vector<hacapk::Point3D> points;
		points.reserve(num_elements);

		for(int i = 0; i < num_elements; i++) {
			points.emplace_back(
				element_positions[3*i],
				element_positions[3*i+1],
				element_positions[3*i+2]
			);
		}

		// Convert element moments to TVector3d
		std::vector<TVector3d> moments;
		moments.reserve(num_elements);

		for(int i = 0; i < num_elements; i++) {
			moments.emplace_back(
				element_moments[3*i],
				element_moments[3*i+1],
				element_moments[3*i+2]
			);
		}

		// Setup HACApK control parameters
		hacapk::ControlParams params;
		params.eps_aca = config.eps;
		params.leaf_size = static_cast<double>(config.min_cluster_size);
		params.aca_type = 1;  // Standard ACA
		params.eta = 2.0;     // Distance parameter for admissibility
		params.print_level = 1;

		if(config.use_openmp) {
			if(config.num_threads > 0) {
				hacapk::set_num_threads(config.num_threads);
			}
			std::cout << "[HMatrix] Using " << hacapk::get_num_threads() << " OpenMP threads" << std::endl;
		}

		// Prepare kernel data
		KernelData kdata;
		kdata.source_points = &points;
		kdata.target_points = &points;  // Self-interaction for now
		kdata.magnetic_moments = &moments;

		// Define kernel function (Biot-Savart law)
		auto kernel = [](int i, int j, void* user_data) -> double {
			KernelData* kd = static_cast<KernelData*>(user_data);

			const hacapk::Point3D& pi = (*kd->target_points)[i];
			const hacapk::Point3D& pj = (*kd->source_points)[j];
			const TVector3d& mj = (*kd->magnetic_moments)[j];

			// Distance vector (mm)
			double dx = pi.x - pj.x;
			double dy = pi.y - pj.y;
			double dz = pi.z - pj.z;
			double r = std::sqrt(dx*dx + dy*dy + dz*dz);

			if(r < 1e-10) {
				// Self-interaction or very close points
				return 0.0;
			}

			// Biot-Savart law: B = (mu0/4*pi) * (m x r) / r^3
			// For simplicity, return scalar influence (magnitude)
			// Full implementation should handle vector field
			const double mu0_over_4pi = 1e-7;  // T*m/A
			double r3 = r*r*r;

			// Cross product: m x r
			// For now, approximate with m*r term (simplified)
			double influence = mu0_over_4pi / r3;

			return influence;
		};

		// Build H-matrix
		std::cout << "[HMatrix] Calling hacapk::build_hmatrix..." << std::endl;
		std::unique_ptr<hacapk::HMatrix> hmat = hacapk::build_hmatrix(
			points, points, kernel, &kdata, params
		);

		if(!hmat) {
			std::cerr << "[HMatrix] Error: H-matrix construction failed" << std::endl;
			return 0;
		}

		// Store H-matrix
		hmatrix_data = hmat.release();
		hmatrix_memory = static_cast<hacapk::HMatrix*>(hmatrix_data)->memory_usage();

		is_built = true;

		std::cout << "[HMatrix] H-matrix construction successful" << std::endl;
		std::cout << "  Number of blocks: " << static_cast<hacapk::HMatrix*>(hmatrix_data)->nlf << std::endl;
		std::cout << "  Low-rank blocks: " << static_cast<hacapk::HMatrix*>(hmatrix_data)->nlfkt << std::endl;
		std::cout << "  Max rank: " << static_cast<hacapk::HMatrix*>(hmatrix_data)->ktmax << std::endl;
		std::cout << "  Memory usage: " << hmatrix_memory / 1024.0 / 1024.0 << " MB" << std::endl;
		std::cout << "  Compression ratio: " << static_cast<hacapk::HMatrix*>(hmatrix_data)->compression_ratio() << std::endl;

	} catch(const std::exception& e) {
		std::cerr << "[HMatrix] Exception during construction: " << e.what() << std::endl;
		return 0;
	} catch(...) {
		std::cerr << "[HMatrix] Unknown exception during construction" << std::endl;
		return 0;
	}

	auto end_time = std::chrono::high_resolution_clock::now();
	construction_time = std::chrono::duration<double>(end_time - start_time).count();

	std::cout << "[HMatrix] Construction completed in " << construction_time << " seconds" << std::endl;

	return 1;
}

//-------------------------------------------------------------------------

void radTHMatrixFieldSource::B_comp(radTField* FieldPtr)
{
	if(!FieldPtr) return;

	// If H-matrix is not available, use direct calculation
	if(!is_built || hmatrix_data == nullptr) {
		B_comp_direct(FieldPtr);
		return;
	}

	// Phase 3: Fast field evaluation for arbitrary points
	try {
		// Evaluation point (mm)
		TVector3d P = FieldPtr->P;

		// For arbitrary points, we compute direct contribution from each element
		// This is still O(N), but we can use OpenMP for parallelization
		// Note: H-matrix is most beneficial for batch evaluation or
		// field calculations at element centers (self-consistent fields)

		#ifdef _OPENMP
		if(config.use_openmp) {
			B_comp_direct_openmp(FieldPtr);
		} else {
			B_comp_direct(FieldPtr);
		}
		#else
		B_comp_direct(FieldPtr);
		#endif

	} catch(...) {
		// Fallback to direct calculation on error
		B_comp_direct(FieldPtr);
	}
}

//-------------------------------------------------------------------------

void radTHMatrixFieldSource::B_intComp(radTField* FieldPtr)
{
	if(!FieldPtr) return;

	// For field integrals, use direct calculation
	// H-matrix is primarily for point evaluations
	B_comp_direct(FieldPtr);
}

//-------------------------------------------------------------------------

void radTHMatrixFieldSource::B_comp_batch(std::vector<radTField*>& fields)
{
	if(fields.empty()) return;

	if(!is_built || hmatrix_data == nullptr) {
		// Use direct calculation for each point
		#ifdef _OPENMP
		if(config.use_openmp) {
			// Parallel batch evaluation
			int n_fields = static_cast<int>(fields.size());

			#pragma omp parallel for schedule(dynamic)
			for(int i = 0; i < n_fields; i++) {
				if(fields[i]) {
					B_comp_direct(fields[i]);
				}
			}
		} else {
			// Serial batch evaluation
			for(radTField* field : fields) {
				if(field) {
					B_comp_direct(field);
				}
			}
		}
		#else
		// Serial batch evaluation (no OpenMP)
		for(radTField* field : fields) {
			if(field) {
				B_comp_direct(field);
			}
		}
		#endif
		return;
	}

	// Phase 3: Optimized batch evaluation with OpenMP
	try {
		#ifdef _OPENMP
		if(config.use_openmp) {
			// Parallel batch evaluation
			int n_fields = static_cast<int>(fields.size());

			#pragma omp parallel for schedule(dynamic)
			for(int i = 0; i < n_fields; i++) {
				if(fields[i]) {
					B_comp_direct(fields[i]);
				}
			}
		} else {
			// Serial batch evaluation
			for(radTField* field : fields) {
				if(field) {
					B_comp(field);
				}
			}
		}
		#else
		// Serial batch evaluation (no OpenMP)
		for(radTField* field : fields) {
			if(field) {
				B_comp(field);
			}
		}
		#endif

	} catch(...) {
		// Fallback to direct calculation
		for(radTField* field : fields) {
			if(field) {
				B_comp_direct(field);
			}
		}
	}
}

//-------------------------------------------------------------------------

void radTHMatrixFieldSource::B_comp_direct(radTField* FieldPtr)
{
	if(!FieldPtr) return;

	// Direct calculation: sum contributions from all source elements
	for(const auto& pair : source_elements) {
		radTg3d* elem = static_cast<radTg3d*>(pair.second.rep);
		if(elem) {
			elem->B_genComp(FieldPtr);
		}
	}
}

//-------------------------------------------------------------------------

void radTHMatrixFieldSource::B_comp_direct_openmp(radTField* FieldPtr)
{
	if(!FieldPtr) return;

	#ifdef _OPENMP
	// OpenMP-parallelized direct calculation
	// Each thread computes contributions from a subset of elements

	// Collect element pointers into vector for indexed access
	std::vector<radTg3d*> elements;
	elements.reserve(source_elements.size());

	for(const auto& pair : source_elements) {
		radTg3d* elem = static_cast<radTg3d*>(pair.second.rep);
		if(elem) {
			elements.push_back(elem);
		}
	}

	int n_elem = static_cast<int>(elements.size());

	// Each thread accumulates its contribution separately
	// Note: Cannot directly update FieldPtr->B in parallel due to race condition
	// Instead, we compute contributions and sum them at the end

	TVector3d B_total(0, 0, 0);
	TVector3d H_total(0, 0, 0);
	TVector3d A_total(0, 0, 0);

	// Store original field values
	TVector3d B_orig = FieldPtr->B;
	TVector3d H_orig = FieldPtr->H;
	TVector3d A_orig = FieldPtr->A;

	#pragma omp parallel
	{
		// Each thread has its own field accumulator
		radTField thread_field = *FieldPtr;
		thread_field.B = TVector3d(0, 0, 0);
		thread_field.H = TVector3d(0, 0, 0);
		thread_field.A = TVector3d(0, 0, 0);

		#pragma omp for nowait
		for(int i = 0; i < n_elem; i++) {
			elements[i]->B_genComp(&thread_field);
		}

		// Combine results using critical section
		#pragma omp critical
		{
			B_total += thread_field.B;
			H_total += thread_field.H;
			A_total += thread_field.A;
		}
	}

	// Update field with accumulated values
	FieldPtr->B = B_orig + B_total;
	FieldPtr->H = H_orig + H_total;
	FieldPtr->A = A_orig + A_total;

	#else
	// OpenMP not available, fall back to serial
	B_comp_direct(FieldPtr);
	#endif
}

//-------------------------------------------------------------------------

double radTHMatrixFieldSource::KernelFunction(int i, int j, void* kernel_data)
{
	// This is now handled by the lambda in BuildHMatrix()
	return 0.0;
}

//-------------------------------------------------------------------------

void radTHMatrixFieldSource::Dump(std::ostream& o, int ShortSign)
{
	radTg3d::Dump(o);

	o << "H-matrix field source (HACApK-based acceleration)";
	if(ShortSign == 1) return;

	o << std::endl;
	o << "   Number of source elements: " << num_elements << std::endl;
	o << "   H-matrix status: " << (is_built ? "built" : "not built") << std::endl;

	if(is_built && hmatrix_data) {
		hacapk::HMatrix* hmat = static_cast<hacapk::HMatrix*>(hmatrix_data);
		o << "   Construction time: " << construction_time << " seconds" << std::endl;
		o << "   Memory usage: " << hmatrix_memory / 1024.0 / 1024.0 << " MB" << std::endl;
		o << "   Number of blocks: " << hmat->nlf << std::endl;
		o << "   Low-rank blocks: " << hmat->nlfkt << std::endl;
		o << "   Max rank: " << hmat->ktmax << std::endl;
		o << "   Compression ratio: " << hmat->compression_ratio() << std::endl;
	}

	o << "   Configuration:" << std::endl;
	o << "     eps = " << config.eps << std::endl;
	o << "     max_rank = " << config.max_rank << std::endl;
	o << "     min_cluster_size = " << config.min_cluster_size << std::endl;
	o << "     use_openmp = " << (config.use_openmp ? "yes" : "no") << std::endl;

	o << "   Memory occupied: " << SizeOfThis() << " bytes" << std::endl;
}

//-------------------------------------------------------------------------

void radTHMatrixFieldSource::DumpBin(
	CAuxBinStrVect& oStr,
	vector<int>& vElemKeysOut,
	map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers,
	int& gUniqueMapKey,
	int elemKey)
{
	// Binary serialization not fully supported
	// H-matrix data cannot be easily serialized
	std::cerr << "[HMatrix] Warning: Binary serialization not supported" << std::endl;
}

//-------------------------------------------------------------------------

radTg3dGraphPresent* radTHMatrixFieldSource::CreateGraphPresent()
{
	// No direct geometry visualization for H-matrix
	// Could visualize source elements if needed
	return nullptr;
}

//-------------------------------------------------------------------------

int radTHMatrixFieldSource::DuplicateItself(
	radThg& hg,
	radTApplication* radPtr,
	char PutNewStuffIntoGenCont)
{
	radTHMatrixFieldSource* pNew = nullptr;
	try {
		pNew = new radTHMatrixFieldSource(*this);
		if(!pNew) {
			radTSend Send;
			Send.ErrorMessage("Radia::Error900");
			return 0;
		}

		return FinishDuplication(pNew, hg);
	}
	catch(...) {
		if(pNew) delete pNew;  // Clean up if exception in FinishDuplication
		radTSend Send;
		Send.ErrorMessage("Radia::Error900");
		return 0;
	}
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// Global functions for Python/C API
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

extern radTApplication* g_pRadApp; // Defined in radappl.cpp

void CreateHMatrixFieldSource(int grpKey, double eps, int max_rank, int min_cluster_size, int use_openmp, int num_threads)
{
	radTHMatrixFieldSource* pHMat = nullptr;
	try {
		// Get group object
		radThg hGroup;
		if(!g_pRadApp->RetrieveObject(grpKey, hGroup)) {
			throw std::runtime_error("H-matrix: Invalid group key");
		}

		radTGroup* pGroup = dynamic_cast<radTGroup*>(static_cast<radTg*>(hGroup.rep));
		if(!pGroup) {
			throw std::runtime_error("H-matrix: Object is not a group");
		}

		// Create H-matrix configuration
		radTHMatrixConfig config;
		config.eps = eps;
		config.max_rank = max_rank;
		config.min_cluster_size = min_cluster_size;
		config.use_openmp = (use_openmp != 0);
		config.num_threads = num_threads;

		// Create H-matrix field source
		pHMat = new radTHMatrixFieldSource(pGroup, config);
		if(!pHMat) {
			throw std::runtime_error("H-matrix: Failed to create H-matrix field source");
		}

		// Register in application
		radThg hHMat;
		hHMat.rep = pHMat;
		pHMat = nullptr;  // Ownership transferred to handle
		int newKey = g_pRadApp->AddElementToContainer(hHMat);

		// Output key
		g_pRadApp->OutInt(newKey);

	} catch(const radTException& ex) {
		if(pHMat) delete pHMat;  // Clean up if exception before handle ownership transfer
		radTSend Send;
		Send.ErrorMessage(ex.what());
		g_pRadApp->OutInt(0);
	} catch(const std::exception& ex) {
		if(pHMat) delete pHMat;  // Clean up if exception before handle ownership transfer
		radTSend Send;
		Send.ErrorMessage(ex.what());
		g_pRadApp->OutInt(0);
	} catch(...) {
		if(pHMat) delete pHMat;  // Clean up if exception before handle ownership transfer
		radTSend Send;
		Send.ErrorMessage("H-matrix: Unknown error");
		g_pRadApp->OutInt(0);
	}
}

//-------------------------------------------------------------------------

void BuildHMatrixFieldSource(int hmatKey)
{
	try {
		// Get H-matrix object
		radThg hHMat;
		if(!g_pRadApp->RetrieveObject(hmatKey, hHMat)) {
			throw std::runtime_error("H-matrix: Invalid H-matrix key");
		}

		radTHMatrixFieldSource* pHMat = dynamic_cast<radTHMatrixFieldSource*>(static_cast<radTg*>(hHMat.rep));
		if(!pHMat) {
			throw std::runtime_error("H-matrix: Object is not an H-matrix field source");
		}

		// Build H-matrix
		int result = pHMat->BuildHMatrix();
		if(result != 1) {
			throw std::runtime_error("H-matrix: Build failed");
		}

	} catch(const radTException& ex) {
		radTSend Send;
		Send.ErrorMessage(ex.what());
	} catch(const std::exception& ex) {
		radTSend Send;
		Send.ErrorMessage(ex.what());
	} catch(...) {
		radTSend Send;
		Send.ErrorMessage("H-matrix: Unknown error during build");
	}
}

//-------------------------------------------------------------------------
