/**
 * @file test_hmatrix_radia.cpp
 * @brief Radia integration test with OpenMP H-matrix acceleration
 *
 * Demonstrates how H-matrix can accelerate magnetic field calculations
 * for large coil systems using OpenMP parallelization.
 *
 * @date 2025-11-07
 */

#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

extern "C" {
#include "hacapk_openmp.h"
}

using namespace std;

// Physical constants
const double MU0 = 4.0 * M_PI * 1e-7;  // Permeability of free space [H/m]
const double MU0_4PI = MU0 / (4.0 * M_PI);

/**
 * Current element (simplified Radia-style)
 */
struct CurrentElement {
	Point3D position;    // Center position [mm]
	Point3D direction;   // Current direction (unit vector)
	double length;       // Length [mm]
	double current;      // Current [A]
};

/**
 * Create a circular coil with N segments
 */
vector<CurrentElement> create_circular_coil(
	double radius,      // [mm]
	int n_segments,
	double current      // [A]
) {
	vector<CurrentElement> elements;
	elements.reserve(n_segments);

	double d_theta = 2.0 * M_PI / n_segments;

	for (int i = 0; i < n_segments; i++) {
		double theta = i * d_theta;

		CurrentElement elem;

		// Center of segment
		elem.position.x = radius * cos(theta + d_theta/2);
		elem.position.y = radius * sin(theta + d_theta/2);
		elem.position.z = 0.0;

		// Direction (tangent to circle)
		elem.direction.x = -sin(theta + d_theta/2);
		elem.direction.y = cos(theta + d_theta/2);
		elem.direction.z = 0.0;

		// Segment length
		elem.length = radius * d_theta;
		elem.current = current;

		elements.push_back(elem);
	}

	return elements;
}

/**
 * Direct magnetic field calculation (Biot-Savart law)
 * B = (μ₀/4π) * I * (dl × R) / |R|³
 */
Point3D compute_field_direct(
	const vector<CurrentElement>& elements,
	const Point3D& field_point
) {
	Point3D B_total = {0, 0, 0};

	for (const auto& elem : elements) {
		// Vector from element to field point
		Point3D R;
		R.x = field_point.x - elem.position.x;
		R.y = field_point.y - elem.position.y;
		R.z = field_point.z - elem.position.z;

		double r = sqrt(R.x*R.x + R.y*R.y + R.z*R.z) + 1e-10;  // Avoid division by zero
		double r3 = r * r * r;

		// dl × R (cross product)
		double dl_cross_R_x = elem.direction.y * R.z - elem.direction.z * R.y;
		double dl_cross_R_y = elem.direction.z * R.x - elem.direction.x * R.z;
		double dl_cross_R_z = elem.direction.x * R.y - elem.direction.y * R.x;

		// Biot-Savart contribution
		double factor = MU0_4PI * elem.current * elem.length / r3;

		B_total.x += factor * dl_cross_R_x;
		B_total.y += factor * dl_cross_R_y;
		B_total.z += factor * dl_cross_R_z;
	}

	return B_total;
}

/**
 * OpenMP-parallelized direct field calculation
 */
Point3D compute_field_direct_openmp(
	const vector<CurrentElement>& elements,
	const Point3D& field_point
) {
	// OpenMP 2.0 doesn't support struct member reduction
	// Use separate variables instead
	double Bx_sum = 0.0, By_sum = 0.0, Bz_sum = 0.0;
	int n_elements = (int)elements.size();

	#pragma omp parallel for reduction(+:Bx_sum,By_sum,Bz_sum)
	for (int idx = 0; idx < n_elements; idx++) {
		const CurrentElement& elem = elements[idx];

		// Vector from element to field point
		Point3D R;
		R.x = field_point.x - elem.position.x;
		R.y = field_point.y - elem.position.y;
		R.z = field_point.z - elem.position.z;

		double r = sqrt(R.x*R.x + R.y*R.y + R.z*R.z) + 1e-10;
		double r3 = r * r * r;

		// dl × R (cross product)
		double dl_cross_R_x = elem.direction.y * R.z - elem.direction.z * R.y;
		double dl_cross_R_y = elem.direction.z * R.x - elem.direction.x * R.z;
		double dl_cross_R_z = elem.direction.x * R.y - elem.direction.y * R.x;

		// Biot-Savart contribution
		double factor = MU0_4PI * elem.current * elem.length / r3;

		Bx_sum += factor * dl_cross_R_x;
		By_sum += factor * dl_cross_R_y;
		Bz_sum += factor * dl_cross_R_z;
	}

	Point3D B_total;
	B_total.x = Bx_sum;
	B_total.y = By_sum;
	B_total.z = Bz_sum;
	return B_total;
}

/**
 * Test: Compare direct vs OpenMP parallelized calculation
 */
void test_field_calculation() {
	cout << "\n" << string(70, '=') << endl;
	cout << "Test: Magnetic Field Calculation with OpenMP" << endl;
	cout << string(70, '=') << endl;

	// Create a circular coil
	double radius = 100.0;  // mm
	int n_segments = 256;   // More segments for performance test
	double current = 1000.0;  // A

	cout << "\nCreating circular coil:" << endl;
	cout << "  Radius:     " << radius << " mm" << endl;
	cout << "  Segments:   " << n_segments << endl;
	cout << "  Current:    " << current << " A" << endl;

	auto elements = create_circular_coil(radius, n_segments, current);

	// Field point (on axis, 50mm above coil)
	Point3D field_point = {0, 0, 50};

	cout << "\nField point: (" << field_point.x << ", "
	     << field_point.y << ", " << field_point.z << ") mm" << endl;

	// Serial calculation
	cout << "\n" << string(70, '-') << endl;
	cout << "Serial Calculation" << endl;
	cout << string(70, '-') << endl;

	auto t_start = chrono::high_resolution_clock::now();
	Point3D B_serial = compute_field_direct(elements, field_point);
	auto t_end = chrono::high_resolution_clock::now();
	double time_serial = chrono::duration<double, milli>(t_end - t_start).count();

	cout << "  Bx = " << B_serial.x * 1000 << " mT" << endl;
	cout << "  By = " << B_serial.y * 1000 << " mT" << endl;
	cout << "  Bz = " << B_serial.z * 1000 << " mT" << endl;
	double B_mag = sqrt(B_serial.x*B_serial.x + B_serial.y*B_serial.y + B_serial.z*B_serial.z);
	cout << "  |B| = " << B_mag * 1000 << " mT" << endl;
	cout << "  Time: " << time_serial << " ms" << endl;

	// OpenMP parallel calculation
	cout << "\n" << string(70, '-') << endl;
	cout << "OpenMP Parallel Calculation" << endl;
	cout << string(70, '-') << endl;

	int n_threads = hmatrix_get_threads();
	cout << "  Using " << n_threads << " OpenMP threads" << endl;

	t_start = chrono::high_resolution_clock::now();
	Point3D B_openmp = compute_field_direct_openmp(elements, field_point);
	t_end = chrono::high_resolution_clock::now();
	double time_openmp = chrono::duration<double, milli>(t_end - t_start).count();

	cout << "  Bx = " << B_openmp.x * 1000 << " mT" << endl;
	cout << "  By = " << B_openmp.y * 1000 << " mT" << endl;
	cout << "  Bz = " << B_openmp.z * 1000 << " mT" << endl;
	B_mag = sqrt(B_openmp.x*B_openmp.x + B_openmp.y*B_openmp.y + B_openmp.z*B_openmp.z);
	cout << "  |B| = " << B_mag * 1000 << " mT" << endl;
	cout << "  Time: " << time_openmp << " ms" << endl;

	double speedup = time_serial / time_openmp;
	cout << "\n  Speedup: " << speedup << "x" << endl;

	// Verify correctness
	double error = sqrt(
		pow(B_serial.x - B_openmp.x, 2) +
		pow(B_serial.y - B_openmp.y, 2) +
		pow(B_serial.z - B_openmp.z, 2)
	);

	cout << "  Error: " << error * 1e6 << " µT" << endl;

	if (error < 1e-10) {
		cout << "  [PASS] Results match!" << endl;
	} else {
		cout << "  [FAIL] Results differ!" << endl;
	}
}

/**
 * Test: H-matrix structure for coil elements
 */
void test_coil_clustering() {
	cout << "\n" << string(70, '=') << endl;
	cout << "Test: Hierarchical Clustering of Coil Elements" << endl;
	cout << string(70, '=') << endl;

	// Create a multi-turn coil
	double radius = 100.0;
	int turns = 4;
	int segments_per_turn = 32;
	double spacing = 10.0;  // mm between turns
	double current = 500.0;

	cout << "\nCreating multi-turn coil:" << endl;
	cout << "  Turns:               " << turns << endl;
	cout << "  Segments per turn:   " << segments_per_turn << endl;
	cout << "  Total elements:      " << turns * segments_per_turn << endl;

	vector<CurrentElement> all_elements;
	for (int turn = 0; turn < turns; turn++) {
		auto turn_elements = create_circular_coil(radius, segments_per_turn, current);

		// Offset each turn in z-direction
		for (auto& elem : turn_elements) {
			elem.position.z = turn * spacing;
		}

		all_elements.insert(all_elements.end(), turn_elements.begin(), turn_elements.end());
	}

	// Convert to Point3D array
	vector<Point3D> points(all_elements.size());
	for (size_t i = 0; i < all_elements.size(); i++) {
		points[i] = all_elements[i].position;
	}

	cout << "\nBuilding H-matrix structure..." << endl;

	// Create H-matrix (source = target for field calculation within coil)
	HMatrix *hmat = hmatrix_create(
		points.data(), points.size(),
		points.data(), points.size(),
		2.0,   // eta
		1e-6,  // epsilon
		8      // max_leaf_size
	);

	if (hmat) {
		cout << "  [SUCCESS] H-matrix structure created" << endl;
		cout << "    Source points: " << hmat->n_source << endl;
		cout << "    Target points: " << hmat->n_target << endl;
		cout << "    Eta parameter: " << hmat->eta << endl;
		cout << "    ACA tolerance: " << hmat->epsilon << endl;
		cout << "\n  This structure enables O(N log N) field calculations" << endl;

		hmatrix_free(hmat);
	} else {
		cout << "  [ERROR] Failed to create H-matrix structure" << endl;
	}
}

/**
 * Main test driver
 */
int main(int argc, char** argv) {
	cout << string(70, '=') << endl;
	cout << "H-Matrix + Radia Integration Tests (OpenMP)" << endl;
	cout << string(70, '=') << endl;
	cout << "\nThese tests demonstrate OpenMP-accelerated magnetic field" << endl;
	cout << "calculations using hierarchical matrices." << endl;

	// Run tests
	test_field_calculation();
	test_coil_clustering();

	cout << "\n" << string(70, '=') << endl;
	cout << "Summary" << endl;
	cout << string(70, '=') << endl;
	cout << "\nOpenMP Parallel Benefits:" << endl;
	cout << "  - Shared memory (no MPI complexity)" << endl;
	cout << "  - Easy to integrate with existing Radia code" << endl;
	cout << "  - Good speedup on multi-core systems" << endl;
	cout << "\nH-Matrix Benefits:" << endl;
	cout << "  - O(N log N) instead of O(N²) complexity" << endl;
	cout << "  - Reduced memory usage" << endl;
	cout << "  - 10-100x speedup for large systems (N > 1000)" << endl;
	cout << "\nNext Steps:" << endl;
	cout << "  1. Implement full H-matrix matrix-vector multiplication" << endl;
	cout << "  2. Integrate with Radia's field calculation routines" << endl;
	cout << "  3. Benchmark on large coil systems" << endl;
	cout << string(70, '=') << endl;

	return 0;
}
