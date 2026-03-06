/**
 * @file test_hmatrix_basic.cpp
 * @brief Basic H-matrix functionality test with OpenMP
 *
 * Tests:
 * - Cluster tree generation
 * - Bounding box computation
 * - Admissibility check
 * - ACA low-rank approximation
 *
 * @date 2025-11-07
 */

#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>

extern "C" {
#include "hacapk_openmp.h"
}

using namespace std;

// Test result tracking
struct TestResult {
	int passed = 0;
	int failed = 0;

	void report(const char* test_name, bool success) {
		if (success) {
			cout << "[PASS] " << test_name << endl;
			passed++;
		} else {
			cout << "[FAIL] " << test_name << endl;
			failed++;
		}
	}

	void summary() {
		cout << "\n" << string(70, '=') << endl;
		cout << "Test Summary" << endl;
		cout << string(70, '=') << endl;
		cout << "Passed: " << passed << endl;
		cout << "Failed: " << failed << endl;
		cout << "Total:  " << (passed + failed) << endl;
		cout << string(70, '=') << endl;
	}
};

/**
 * Test 1: OpenMP configuration
 */
bool test_openmp_config() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 1: OpenMP Configuration" << endl;
	cout << string(70, '-') << endl;

	int n_threads = hmatrix_get_threads();
	cout << "  Default OpenMP threads: " << n_threads << endl;

	// Try setting threads
	hmatrix_set_threads(4);
	int new_threads = hmatrix_get_threads();
	cout << "  After setting to 4: " << new_threads << endl;

	return (n_threads > 0 && new_threads > 0);
}

/**
 * Test 2: Bounding box computation
 */
bool test_bounding_box() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 2: Bounding Box Computation" << endl;
	cout << string(70, '-') << endl;

	// Create test points
	Point3D points[] = {
		{0.0, 0.0, 0.0},
		{1.0, 0.0, 0.0},
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
		{0.5, 0.5, 0.5}
	};
	int n_points = 5;

	int indices[] = {0, 1, 2, 3, 4};

	BoundingBox bbox;
	hmatrix_compute_bbox(points, indices, 0, n_points, &bbox);

	cout << "  Bounding box:" << endl;
	cout << "    X: [" << bbox.min[0] << ", " << bbox.max[0] << "]" << endl;
	cout << "    Y: [" << bbox.min[1] << ", " << bbox.max[1] << "]" << endl;
	cout << "    Z: [" << bbox.min[2] << ", " << bbox.max[2] << "]" << endl;

	double diam = bbox_diameter(&bbox);
	cout << "    Diameter: " << diam << endl;

	// Check correctness
	bool correct = (
		bbox.min[0] == 0.0 && bbox.max[0] == 1.0 &&
		bbox.min[1] == 0.0 && bbox.max[1] == 1.0 &&
		bbox.min[2] == 0.0 && bbox.max[2] == 1.0
	);

	return correct;
}

/**
 * Test 3: Cluster tree generation
 */
bool test_cluster_tree() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 3: Cluster Tree Generation" << endl;
	cout << string(70, '-') << endl;

	// Create grid of points
	int n = 8;
	vector<Point3D> points;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				Point3D p;
				p.x = i * 10.0;
				p.y = j * 10.0;
				p.z = k * 10.0;
				points.push_back(p);
			}
		}
	}

	int n_points = points.size();
	cout << "  Created " << n_points << " points in grid" << endl;

	// Create index array
	vector<int> indices(n_points);
	for (int i = 0; i < n_points; i++) {
		indices[i] = i;
	}

	// Build cluster tree
	int max_leaf_size = 32;
	Cluster *cluster = hmatrix_create_cluster(
		points.data(), indices.data(), 0, n_points, 0, max_leaf_size
	);

	bool success = (cluster != nullptr);

	if (success) {
		cout << "  Cluster tree created successfully" << endl;
		cout << "    Root cluster size: " << cluster->size << endl;
		cout << "    Has children: " << (cluster->child[0] != nullptr ? "Yes" : "No") << endl;

		// Free cluster
		hmatrix_free_cluster(cluster);
	} else {
		cout << "  [ERROR] Failed to create cluster tree" << endl;
	}

	return success;
}

/**
 * Test 4: Admissibility check
 */
bool test_admissibility() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 4: Admissibility Check" << endl;
	cout << string(70, '-') << endl;

	// Create two well-separated bounding boxes
	BoundingBox b1, b2;

	b1.min[0] = 0.0; b1.max[0] = 1.0;
	b1.min[1] = 0.0; b1.max[1] = 1.0;
	b1.min[2] = 0.0; b1.max[2] = 1.0;

	b2.min[0] = 10.0; b2.max[0] = 11.0;
	b2.min[1] = 10.0; b2.max[1] = 11.0;
	b2.min[2] = 10.0; b2.max[2] = 11.0;

	double eta = 2.0;
	int admissible = hmatrix_is_admissible(&b1, &b2, eta);

	cout << "  Box 1: [0,1]³" << endl;
	cout << "  Box 2: [10,11]³" << endl;
	cout << "  Distance: " << bbox_distance(&b1, &b2) << endl;
	cout << "  Eta parameter: " << eta << endl;
	cout << "  Admissible: " << (admissible ? "Yes" : "No") << endl;

	// Should be admissible (far apart)
	return (admissible == 1);
}

/**
 * Test 5: Simple ACA approximation
 */
bool test_aca() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 5: ACA Low-Rank Approximation" << endl;
	cout << string(70, '-') << endl;

	// Simple test kernel: K(i,j) = 1/(1 + |i-j|)
	// This should have low rank
	struct KernelData {
		// Empty for this simple test
	};

	KernelData data;

	auto kernel = [](int i, int j, void *user_data) -> double {
		(void)user_data;
		double dist = abs(i - j);
		return 1.0 / (1.0 + dist);
	};

	int m = 20, n = 20;
	double epsilon = 1e-3;

	cout << "  Matrix size: " << m << " × " << n << endl;
	cout << "  Tolerance: " << epsilon << endl;

	LowRankBlock *block = hmatrix_aca(m, n, kernel, &data, epsilon);

	bool success = (block != nullptr && block->rank > 0 && block->rank < m);

	if (success) {
		cout << "  ACA completed successfully" << endl;
		cout << "    Achieved rank: " << block->rank << " (max " << m << ")" << endl;
		cout << "    Compression ratio: " << (double)(m*n) / (double)(block->rank * (m+n)) << endl;

		hmatrix_free_lowrank(block);
	} else {
		cout << "  [ERROR] ACA failed" << endl;
	}

	return success;
}

/**
 * Test 6: H-matrix creation
 */
bool test_hmatrix_creation() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 6: H-Matrix Creation" << endl;
	cout << string(70, '-') << endl;

	// Create source and target points
	int n_source = 64;
	int n_target = 64;

	vector<Point3D> source_points(n_source);
	vector<Point3D> target_points(n_target);

	// Create grid of source points
	for (int i = 0; i < n_source; i++) {
		source_points[i].x = (i % 8) * 10.0;
		source_points[i].y = ((i / 8) % 8) * 10.0;
		source_points[i].z = (i / 64) * 10.0;
	}

	// Create grid of target points (offset)
	for (int i = 0; i < n_target; i++) {
		target_points[i].x = (i % 8) * 10.0 + 5.0;
		target_points[i].y = ((i / 8) % 8) * 10.0 + 5.0;
		target_points[i].z = (i / 64) * 10.0 + 5.0;
	}

	cout << "  Source points: " << n_source << endl;
	cout << "  Target points: " << n_target << endl;

	HMatrix *hmat = hmatrix_create(
		source_points.data(), n_source,
		target_points.data(), n_target,
		2.0,   // eta
		1e-6,  // epsilon
		16     // max_leaf_size
	);

	bool success = (hmat != nullptr);

	if (success) {
		cout << "  H-matrix created successfully" << endl;
		cout << "    Source tree depth: " << hmat->source_tree->depth << endl;
		cout << "    Target tree depth: " << hmat->target_tree->depth << endl;

		hmatrix_free(hmat);
	} else {
		cout << "  [ERROR] H-matrix creation failed" << endl;
	}

	return success;
}

/**
 * Main test driver
 */
int main(int argc, char** argv) {
	cout << string(70, '=') << endl;
	cout << "H-Matrix OpenMP Basic Functionality Tests" << endl;
	cout << string(70, '=') << endl;
	cout << "\nThese tests verify OpenMP-based H-matrix operations:" << endl;
	cout << "  - OpenMP configuration" << endl;
	cout << "  - Bounding box computation" << endl;
	cout << "  - Cluster tree generation" << endl;
	cout << "  - Admissibility conditions" << endl;
	cout << "  - ACA low-rank approximation" << endl;
	cout << "  - H-matrix structure creation" << endl;

	TestResult results;

	// Run tests
	results.report("OpenMP Configuration", test_openmp_config());
	results.report("Bounding Box Computation", test_bounding_box());
	results.report("Cluster Tree Generation", test_cluster_tree());
	results.report("Admissibility Check", test_admissibility());
	results.report("ACA Low-Rank Approximation", test_aca());
	results.report("H-Matrix Creation", test_hmatrix_creation());

	// Print summary
	results.summary();

	return (results.failed == 0) ? 0 : 1;
}
