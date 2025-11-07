/*-------------------------------------------------------------------------
*
* File name:      radhmat.h
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

#ifndef __RADHMAT_H
#define __RADHMAT_H

#include "radg3d.h"
#include "radgroup.h"
#include <vector>
#include <memory>

// Forward declarations to avoid including HACApK in header
namespace hacapk {
	class HMatrix;
	struct Point3D;
	class Cluster;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

/**
 * Configuration parameters for H-matrix construction
 */
struct radTHMatrixConfig {
	double eps;             // ACA tolerance (default: 1e-6)
	int max_rank;           // Maximum rank for low-rank blocks (default: 50)
	int min_cluster_size;   // Minimum cluster size (default: 10)
	bool use_openmp;        // Enable OpenMP parallelization (default: true)
	int num_threads;        // Number of OpenMP threads (default: 0 = auto)

	radTHMatrixConfig()
		: eps(1e-6)
		, max_rank(50)
		, min_cluster_size(10)
		, use_openmp(true)
		, num_threads(0)
	{}
};

//-------------------------------------------------------------------------

/**
 * H-matrix based field source for fast field computation
 *
 * This class uses HACApK (Hierarchical Adaptive Cross Approximation)
 * to accelerate magnetic field calculations for systems with many elements.
 *
 * Complexity:
 * - Construction: O(N log N)
 * - Field evaluation: O(log N) per point
 * - Batch evaluation: O((N+M) log N) for M points
 *
 * Recommended for N > 100 magnetic elements.
 */
class radTHMatrixFieldSource : public radTg3d {
public:
	radTmhg source_elements;           // Original magnetic elements
	radTHMatrixConfig config;          // Configuration parameters
	bool is_built;                     // H-matrix construction status

	// Internal H-matrix data structures (opaque pointers)
	void* hmatrix_data;                // HACApK H-matrix (void* to avoid header dependency)
	void* cluster_tree_data;           // Cluster tree
	std::vector<double> element_positions;  // Element center positions (flattened x,y,z)
	std::vector<double> element_moments;    // Magnetic moments or currents

	// Statistics
	int num_elements;                  // Number of source elements
	size_t hmatrix_memory;             // Memory used by H-matrix (bytes)
	double construction_time;          // H-matrix construction time (seconds)

public:
	// Constructor from radTGroup
	radTHMatrixFieldSource(radTGroup* group, const radTHMatrixConfig& cfg = radTHMatrixConfig());

	// Constructor for binary deserialization (not fully supported yet)
	radTHMatrixFieldSource(CAuxBinStrVect& inStr, map<int, int>& mKeysOldNew, radTmhg& gMapOfHandlers);

	// Copy constructor
	radTHMatrixFieldSource(const radTHMatrixFieldSource& src);

	// Destructor
	virtual ~radTHMatrixFieldSource();

	// Type identifier
	int Type_g3d() override { return 100; }  // New type for H-matrix

	/**
	 * Build H-matrix from source elements
	 *
	 * This performs:
	 * 1. Geometry extraction from magnetic elements
	 * 2. Cluster tree construction
	 * 3. H-matrix assembly with ACA
	 *
	 * @return 1 on success, 0 on failure
	 */
	int BuildHMatrix();

	/**
	 * Compute magnetic field at a single point
	 *
	 * If H-matrix is not built, falls back to direct calculation.
	 *
	 * @param FieldPtr Pointer to field structure containing evaluation point
	 */
	void B_comp(radTField* FieldPtr) override;

	/**
	 * Compute field integral
	 *
	 * @param FieldPtr Pointer to field structure
	 */
	void B_intComp(radTField* FieldPtr) override;

	/**
	 * Batch field computation (efficient for multiple points)
	 *
	 * This is more efficient than calling B_comp repeatedly.
	 *
	 * @param fields Vector of field structures to evaluate
	 */
	void B_comp_batch(std::vector<radTField*>& fields);

	/**
	 * Get H-matrix construction status
	 *
	 * @return true if H-matrix is built and ready
	 */
	bool IsBuilt() const { return is_built; }

	/**
	 * Get number of source elements
	 */
	int GetNumElements() const { return num_elements; }

	/**
	 * Get H-matrix memory usage in bytes
	 */
	size_t GetMemoryUsage() const { return hmatrix_memory; }

	/**
	 * Get H-matrix construction time in seconds
	 */
	double GetConstructionTime() const { return construction_time; }

	// Standard radTg3d interface
	void Dump(std::ostream& o, int ShortSign = 0) override;
	void DumpBin(CAuxBinStrVect& oStr, vector<int>& vElemKeysOut,
	             map<int, radTHandle<radTg>, less<int> >& gMapOfHandlers,
	             int& gUniqueMapKey, int elemKey) override;
	radTg3dGraphPresent* CreateGraphPresent() override;
	int DuplicateItself(radThg& hg, radTApplication*, char) override;
	int SizeOfThis() override { return sizeof(radTHMatrixFieldSource); }

private:
	/**
	 * Extract geometry and magnetic properties from source elements
	 *
	 * @return 1 on success, 0 on failure
	 */
	int ExtractGeometry();

	/**
	 * Kernel function for H-matrix: magnetic field influence
	 *
	 * Computes the magnetic field at element i due to element j.
	 *
	 * @param i Target element index
	 * @param j Source element index
	 * @param kernel_data User data pointer (unused)
	 * @return Field influence value
	 */
	static double KernelFunction(int i, int j, void* kernel_data);

	/**
	 * Direct field computation (fallback)
	 *
	 * Used when H-matrix is not available or for debugging.
	 *
	 * @param FieldPtr Pointer to field structure
	 */
	void B_comp_direct(radTField* FieldPtr);
};

//-------------------------------------------------------------------------

#endif
