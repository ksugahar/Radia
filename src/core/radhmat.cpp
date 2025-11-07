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

#include "radhmat.h"
#include "radgroup.h"
#include "radexcep.h"
#include <iostream>
#include <chrono>
#include <cmath>

// HACApK includes (will be used in Phase 2)
// #include "hacapk.hpp"

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
	// Phase 2: Add proper cleanup of HACApK objects
	// if (hmatrix_data) {
	//     delete static_cast<hacapk::HMatrix*>(hmatrix_data);
	// }
	// if (cluster_tree_data) {
	//     delete static_cast<hacapk::Cluster*>(cluster_tree_data);
	// }

	hmatrix_data = nullptr;
	cluster_tree_data = nullptr;
}

//-------------------------------------------------------------------------

int radTHMatrixFieldSource::ExtractGeometry()
{
	if(source_elements.empty()) {
		return 0;
	}

	// Reserve space
	element_positions.reserve(num_elements * 3);
	element_moments.reserve(num_elements * 3);

	// Extract center positions and magnetic properties
	for(const auto& pair : source_elements) {
		radTg3d* elem = static_cast<radTg3d*>(pair.second.rep);

		// Get element center (simplified - needs proper implementation)
		TVector3d center(0, 0, 0);

		// Try to get magnetization
		// This is a placeholder - proper implementation in Phase 2
		TVector3d M(0, 0, 0);

		// Store data
		element_positions.push_back(center.x);
		element_positions.push_back(center.y);
		element_positions.push_back(center.z);

		element_moments.push_back(M.x);
		element_moments.push_back(M.y);
		element_moments.push_back(M.z);
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

	// Phase 2: Implement H-matrix construction using HACApK
	// 1. Build cluster tree
	// 2. Assemble H-matrix with ACA
	// 3. Store in hmatrix_data

	// Placeholder for Phase 1
	std::cout << "[HMatrix] Phase 1: H-matrix construction not yet implemented" << std::endl;
	std::cout << "[HMatrix] Falling back to direct calculation" << std::endl;

	auto end_time = std::chrono::high_resolution_clock::now();
	construction_time = std::chrono::duration<double>(end_time - start_time).count();

	// Mark as "built" even though we're using direct calculation
	// This allows testing the rest of the infrastructure
	is_built = true;

	std::cout << "[HMatrix] Construction completed in " << construction_time << " seconds" << std::endl;

	return 1;
}

//-------------------------------------------------------------------------

void radTHMatrixFieldSource::B_comp(radTField* FieldPtr)
{
	if(!FieldPtr) return;

	// Phase 1: Use direct calculation (fallback)
	// Phase 2: Use H-matrix if available
	if(!is_built || hmatrix_data == nullptr) {
		B_comp_direct(FieldPtr);
		return;
	}

	// Phase 2: Implement H-matrix field evaluation
	// TVector3d P = FieldPtr->P;
	// TVector3d B = hmatrix_eval_field(P);
	// if(FieldPtr->FieldKey.B_) FieldPtr->B += B;
	// if(FieldPtr->FieldKey.H_) FieldPtr->H += B;

	// For now, use direct calculation
	B_comp_direct(FieldPtr);
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

	// Phase 2: Implement efficient batch evaluation with H-matrix
	// For now, use direct calculation for each point
	for(radTField* field : fields) {
		if(field) {
			B_comp(field);
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

double radTHMatrixFieldSource::KernelFunction(int i, int j, void* kernel_data)
{
	// Phase 2: Implement kernel function for magnetic field
	// This should compute the field influence from element j to element i
	// using Biot-Savart law or dipole approximation

	// Placeholder
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

	if(is_built) {
		o << "   Construction time: " << construction_time << " seconds" << std::endl;
		o << "   Memory usage: " << hmatrix_memory << " bytes" << std::endl;
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
	radTHMatrixFieldSource* pNew = new radTHMatrixFieldSource(*this);
	if(!pNew) {
		radTSend Send;
		Send.ErrorMessage("Radia::Error900");
		return 0;
	}

	return FinishDuplication(pNew, hg);
}

//-------------------------------------------------------------------------
