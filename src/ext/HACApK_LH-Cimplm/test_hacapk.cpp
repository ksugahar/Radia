/**
 * @file test_hacapk.cpp
 * @brief Test program for HACApK C++ implementation
 *
 * @date 2025-11-07
 */

#define _USE_MATH_DEFINES
#include "hacapk.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>

using namespace hacapk;
using namespace std;

// ============================================================================
// Test Utilities
// ============================================================================

struct TestResults {
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

// ============================================================================
// Test Kernels
// ============================================================================

/**
 * Simple 1D kernel for testing: K(i,j) = 1/(1 + |i-j|)
 */
double kernel_1d(int i, int j, void* data) {
	(void)data;
	double dist = abs(i - j);
	return 1.0 / (1.0 + dist);
}

/**
 * 3D Laplace kernel: K(i,j) = 1/||x_i - x_j||
 */
struct Laplace3DData {
	const vector<Point3D>* points;
};

double kernel_laplace_3d(int i, int j, void* data) {
	auto* d = static_cast<Laplace3DData*>(data);
	const Point3D& pi = (*d->points)[i];
	const Point3D& pj = (*d->points)[j];

	double dist = point_distance(pi, pj);
	if (dist < 1e-10) return 0.0;  // Avoid division by zero

	return 1.0 / dist;
}

// ============================================================================
// Test 1: BoundingBox
// ============================================================================

bool test_bounding_box() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 1: BoundingBox" << endl;
	cout << string(70, '-') << endl;

	BoundingBox bbox(3);
	bbox.min[0] = 0.0; bbox.max[0] = 1.0;
	bbox.min[1] = 0.0; bbox.max[1] = 2.0;
	bbox.min[2] = 0.0; bbox.max[2] = 3.0;

	double width = bbox.width();
	double diam = bbox.diameter();

	cout << "  Width: " << width << " (expected 3.0)" << endl;
	cout << "  Diameter: " << diam << " (expected " << sqrt(14.0) << ")" << endl;

	bool correct = (abs(width - 3.0) < 1e-10) && (abs(diam - sqrt(14.0)) < 1e-10);

	return correct;
}

// ============================================================================
// Test 2: Cluster Generation
// ============================================================================

bool test_cluster_generation() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 2: Cluster Generation" << endl;
	cout << string(70, '-') << endl;

	// Create grid of points
	int n = 8;
	vector<Point3D> points;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				points.emplace_back(i * 10.0, j * 10.0, k * 10.0);
			}
		}
	}

	cout << "  Created " << points.size() << " points in grid" << endl;

	// Create cluster
	vector<int> indices(points.size());
	for (size_t i = 0; i < points.size(); i++) indices[i] = i;

	ControlParams params;
	params.leaf_size = 32;

	auto cluster = generate_cluster(points, indices, 0, points.size(), 0, params);

	bool success = (cluster != nullptr);

	if (success) {
		cout << "  Cluster created successfully" << endl;
		cout << "    Size: " << cluster->nsize << endl;
		cout << "    Depth: " << cluster->ndpth << endl;
		cout << "    Is leaf: " << (cluster->is_leaf() ? "Yes" : "No") << endl;
		cout << "    Num children: " << cluster->sons.size() << endl;
	} else {
		cout << "  [ERROR] Failed to create cluster" << endl;
	}

	return success;
}

// ============================================================================
// Test 3: Admissibility Check
// ============================================================================

bool test_admissibility() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 3: Admissibility Check" << endl;
	cout << string(70, '-') << endl;

	// Two well-separated boxes
	BoundingBox b1(3), b2(3);

	b1.min[0] = 0.0; b1.max[0] = 1.0;
	b1.min[1] = 0.0; b1.max[1] = 1.0;
	b1.min[2] = 0.0; b1.max[2] = 1.0;

	b2.min[0] = 10.0; b2.max[0] = 11.0;
	b2.min[1] = 10.0; b2.max[1] = 11.0;
	b2.min[2] = 10.0; b2.max[2] = 11.0;

	double eta = 2.0;
	bool admissible = is_admissible(b1, b2, eta);

	cout << "  Box 1: [0,1]³" << endl;
	cout << "  Box 2: [10,11]³" << endl;
	cout << "  Distance: " << bbox_distance(b1, b2) << endl;
	cout << "  Eta parameter: " << eta << endl;
	cout << "  Admissible: " << (admissible ? "Yes" : "No") << endl;

	// Should be admissible (far apart)
	return admissible;
}

// ============================================================================
// Test 4: ACA Approximation
// ============================================================================

bool test_aca() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 4: ACA Low-Rank Approximation" << endl;
	cout << string(70, '-') << endl;

	LowRankBlock block;
	block.nstrtl = 0;
	block.nstrtt = 0;
	block.ndl = 20;
	block.ndt = 20;

	double epsilon = 1e-3;

	cout << "  Matrix size: " << block.ndl << " × " << block.ndt << endl;
	cout << "  Tolerance: " << epsilon << endl;

	aca_approximation(block, kernel_1d, nullptr, epsilon);

	cout << "  After ACA:" << endl;
	cout << "    ltmtx: " << block.ltmtx << " (1=lowrank, 2=full)" << endl;
	cout << "    kt (rank): " << block.kt << endl;
	cout << "    is_lowrank(): " << (block.is_lowrank() ? "true" : "false") << endl;

	bool success = (block.is_lowrank() && block.kt > 0 && block.kt < block.ndl);

	if (success) {
		cout << "  ACA completed successfully" << endl;
		cout << "    Achieved rank: " << block.kt << " (max " << block.ndl << ")" << endl;
		double comp_ratio = (double)(block.ndl * block.ndt) / (double)(block.kt * (block.ndl + block.ndt));
		cout << "    Compression ratio: " << comp_ratio << endl;

		// Verify approximation quality
		double max_error = 0.0;
		for (int i = 0; i < block.ndl; i++) {
			for (int j = 0; j < block.ndt; j++) {
				double exact = kernel_1d(i, j, nullptr);

				double approx = 0.0;
				for (int k = 0; k < block.kt; k++) {
					approx += block.a1[i * block.kt + k] * block.a2[j * block.kt + k];
				}

				double error = abs(exact - approx);
				max_error = max(max_error, error);
			}
		}

		cout << "    Maximum error: " << max_error << endl;
		success = (max_error < 0.01);  // Should be accurate
	} else {
		cout << "  [ERROR] ACA failed" << endl;
	}

	return success;
}

// ============================================================================
// Test 5: H-Matrix Construction
// ============================================================================

bool test_hmatrix_construction() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 5: H-Matrix Construction" << endl;
	cout << string(70, '-') << endl;

	// Create grid of points
	int n = 64;
	vector<Point3D> points;

	for (int i = 0; i < n; i++) {
		double x = (i % 8) * 10.0;
		double y = ((i / 8) % 8) * 10.0;
		double z = (i / 64) * 10.0;
		points.emplace_back(x, y, z);
	}

	cout << "  Number of points: " << points.size() << endl;

	ControlParams params;
	params.leaf_size = 16;
	params.eta = 2.0;
	params.eps_aca = 1e-3;

	Laplace3DData kernel_data;
	kernel_data.points = &points;

	auto t_start = chrono::high_resolution_clock::now();

	auto hmat = build_hmatrix(points, points, kernel_laplace_3d, &kernel_data, params);

	auto t_end = chrono::high_resolution_clock::now();
	double time_ms = chrono::duration<double, milli>(t_end - t_start).count();

	bool success = (hmat != nullptr && hmat->nlf > 0);

	if (success) {
		cout << "  H-matrix created successfully" << endl;
		cout << "    Number of leaf blocks: " << hmat->nlf << endl;
		cout << "    Number of low-rank blocks: " << hmat->nlfkt << endl;
		cout << "    Maximum rank: " << hmat->ktmax << endl;
		cout << "    Memory usage: " << hmat->memory_usage() / 1024.0 << " KB" << endl;
		cout << "    Compression ratio: " << hmat->compression_ratio() << endl;
		cout << "    Construction time: " << time_ms << " ms" << endl;
	} else {
		cout << "  [ERROR] H-matrix construction failed" << endl;
	}

	return success;
}

// ============================================================================
// Test 6: Matrix-Vector Multiplication
// ============================================================================

bool test_matvec() {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 6: H-Matrix Matrix-Vector Multiplication" << endl;
	cout << string(70, '-') << endl;

	// Create small grid
	int n = 32;
	vector<Point3D> points;

	for (int i = 0; i < n; i++) {
		double x = (i % 4) * 10.0;
		double y = ((i / 4) % 4) * 10.0;
		double z = (i / 16) * 10.0;
		points.emplace_back(x, y, z);
	}

	ControlParams params;
	params.leaf_size = 8;
	params.eta = 2.0;
	params.eps_aca = 1e-3;

	Laplace3DData kernel_data;
	kernel_data.points = &points;

	auto hmat = build_hmatrix(points, points, kernel_laplace_3d, &kernel_data, params);

	// Create test vector
	vector<double> x(n, 1.0);  // All ones
	vector<double> y(n, 0.0);

	cout << "  Vector size: " << n << endl;
	cout << "  Number of threads: " << get_num_threads() << endl;

	auto t_start = chrono::high_resolution_clock::now();

	hmatrix_matvec(*hmat, x, y);

	auto t_end = chrono::high_resolution_clock::now();
	double time_ms = chrono::duration<double, milli>(t_end - t_start).count();

	// Check result is non-zero
	double norm_y = 0.0;
	for (double val : y) {
		norm_y += val * val;
	}
	norm_y = sqrt(norm_y);

	cout << "  ||y||_2 = " << norm_y << endl;
	cout << "  Time: " << time_ms << " ms" << endl;

	bool success = (norm_y > 1e-6);  // Should be non-zero

	return success;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
	(void)argc;
	(void)argv;

	cout << string(70, '=') << endl;
	cout << "HACApK C++ Implementation Tests" << endl;
	cout << string(70, '=') << endl;
	cout << "\nC++ version of HACApK (Fortran-free)" << endl;
	cout << "OpenMP threads: " << get_num_threads() << endl;

	TestResults results;

	// Run tests
	results.report("BoundingBox", test_bounding_box());
	results.report("Cluster Generation", test_cluster_generation());
	results.report("Admissibility Check", test_admissibility());
	results.report("ACA Low-Rank Approximation", test_aca());
	results.report("H-Matrix Construction", test_hmatrix_construction());
	results.report("Matrix-Vector Multiplication", test_matvec());

	// Print summary
	results.summary();

	return (results.failed == 0) ? 0 : 1;
}
