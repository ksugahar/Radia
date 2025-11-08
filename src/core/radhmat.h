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

	/**
	 * OpenMP-parallelized direct field computation
	 *
	 * Accelerates field calculation using OpenMP parallelization.
	 * Each thread computes contributions from a subset of elements.
	 *
	 * @param FieldPtr Pointer to field structure
	 */
	void B_comp_direct_openmp(radTField* FieldPtr);
};

//-------------------------------------------------------------------------

/**
 * H-matrix based field evaluator for arbitrary observation points
 *
 * This class accelerates rad.Fld() batch evaluation using H-matrix
 * approximation. Achieves O((M+N)log(M+N)) complexity instead of O(M×N).
 *
 * Key features:
 * - Caching for non-linear iterations
 * - Automatic fallback for small problems (N<100)
 * - Batch evaluation only (single points use direct calculation)
 *
 * Recommended for:
 * - N ≥ 100 source elements
 * - M ≥ 100 observation points
 * - Non-linear iteration problems
 */
class radTHMatrixFieldEvaluator {
public:
	// Configuration
	radTHMatrixConfig config;
	bool is_built;

	// Source geometry (from magnetic elements)
	int num_sources;
	std::vector<double> source_positions;   // Flattened [x,y,z,x,y,z,...]
	std::vector<double> source_moments;     // Magnetic moments [Mx,My,Mz,...]

	// H-matrix data structures (opaque pointers)
	void* hmatrix_data;           // HACApK H-matrix
	void* source_cluster_tree;    // Source cluster tree
	void* target_cluster_tree;    // Observation point cluster tree

	// Cache key for invalidation
	size_t geometry_hash;

	// Statistics
	size_t memory_usage;          // Total memory (bytes)
	double build_time;            // Construction time (seconds)
	double last_eval_time;        // Last evaluation time (seconds)
	int num_evaluations;          // Number of times reused

public:
	/**
	 * Constructor
	 */
	radTHMatrixFieldEvaluator(const radTHMatrixConfig& cfg = radTHMatrixConfig());

	/**
	 * Destructor
	 */
	~radTHMatrixFieldEvaluator();

	/**
	 * Build H-matrix from source elements
	 *
	 * Extracts geometry from radTGroup and constructs cluster trees
	 * and H-matrix for field evaluation.
	 *
	 * @param source_group Group containing magnetic elements
	 * @return 1 on success, 0 on failure
	 */
	int Build(radTGroup* source_group);

	/**
	 * Evaluate field at multiple observation points (batch)
	 *
	 * This is the main interface for H-matrix accelerated field evaluation.
	 *
	 * Complexity: O((M+N)log(M+N)) instead of O(M×N)
	 *
	 * @param obs_points Vector of observation points (3D coordinates)
	 * @param field_out Vector to store output field values (3D vectors)
	 * @param field_type Field type: 'h'=H-field, 'b'=B-field, 'a'=vector potential
	 * @return 1 on success, 0 on failure
	 */
	int EvaluateField(
		const std::vector<TVector3d>& obs_points,
		std::vector<TVector3d>& field_out,
		char field_type = 'h'
	);

	/**
	 * Check if H-matrix is valid for current geometry
	 *
	 * Compares stored geometry hash with current geometry.
	 *
	 * @param source_group Current source geometry
	 * @return true if cache is valid
	 */
	bool IsValid(radTGroup* source_group);

	/**
	 * Clear H-matrix data and free memory
	 */
	void Clear();

	/**
	 * Get memory usage
	 */
	size_t GetMemoryUsage() const { return memory_usage; }

	/**
	 * Get build time
	 */
	double GetBuildTime() const { return build_time; }

	/**
	 * Get number of times H-matrix was reused
	 */
	int GetNumEvaluations() const { return num_evaluations; }

private:
	/**
	 * Extract geometry from source group
	 *
	 * @param source_group Group containing elements
	 * @return 1 on success
	 */
	int ExtractSourceGeometry(radTGroup* source_group);

	/**
	 * Build cluster tree for observation points
	 *
	 * @param obs_points Vector of observation points
	 * @return 1 on success
	 */
	int BuildTargetClusterTree(const std::vector<TVector3d>& obs_points);

	/**
	 * Build H-matrix for field evaluation
	 *
	 * Constructs H-matrix representing field influence from sources to targets.
	 *
	 * @return 1 on success
	 */
	int BuildFieldHMatrix();

	/**
	 * Magnetic field kernel function
	 *
	 * Computes field at target i due to source j.
	 * Kernel: H(r) = (3(m·r̂)r̂ - m) / (4π|r|³)
	 *
	 * @param i Target point index
	 * @param j Source element index
	 * @param kernel_data Pointer to kernel data structure
	 * @return Kernel value
	 */
	static double FieldKernel(int i, int j, void* kernel_data);

	/**
	 * Compute geometry hash for cache validation
	 *
	 * @param source_group Source geometry
	 * @return Hash value
	 */
	size_t ComputeGeometryHash(radTGroup* source_group);

	/**
	 * Fallback to direct calculation
	 *
	 * Used when H-matrix construction fails or for debugging.
	 *
	 * @param obs_points Observation points
	 * @param field_out Output field values
	 * @param field_type Field type
	 * @return 1 on success
	 */
	int EvaluateFieldDirect(
		const std::vector<TVector3d>& obs_points,
		std::vector<TVector3d>& field_out,
		char field_type
	);
};

//-------------------------------------------------------------------------

#endif
