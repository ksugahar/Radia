/**
 * @file test_radhmat.cpp
 * @brief Test program for radTHMatrixFieldSource (Phase 2 validation)
 *
 * This test validates the H-matrix integration into Radia:
 * - Geometry extraction from radTGroup
 * - H-matrix construction with HACApK
 * - Field calculation accuracy vs direct calculation
 * - Performance comparison
 *
 * @date 2025-11-07
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <cmath>

// Radia core headers
#include "radentry.h"
#include "radgroup.h"
#include "radrec.h"
#include "radmater.h"
#include "radhmat.h"
#include "radsend.h"

using namespace std;

//-------------------------------------------------------------------------
// Test utilities
//-------------------------------------------------------------------------

struct TestResult {
	int passed = 0;
	int failed = 0;

	void report(const char* test_name, bool success) {
		if(success) {
			cout << "  [PASS] " << test_name << endl;
			passed++;
		} else {
			cout << "  [FAIL] " << test_name << endl;
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

//-------------------------------------------------------------------------
// Test 1: Create simple magnetic system
//-------------------------------------------------------------------------

bool test_create_magnetic_system(radTGroup*& group) {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 1: Create Simple Magnetic System" << endl;
	cout << string(70, '-') << endl;

	try {
		// Create a simple array of rectangular magnets
		int n_magnets = 10;
		double spacing = 50.0;  // mm

		cout << "  Creating " << n_magnets << " rectangular magnets..." << endl;

		group = new radTGroup();

		for(int i = 0; i < n_magnets; i++) {
			// Center position
			TVector3d center(i * spacing, 0, 0);

			// Dimensions (mm)
			TVector3d size(20.0, 20.0, 20.0);

			// Magnetization (T)
			TVector3d M(0, 0, 1.0);

			// Create rectangular magnet
			radTRecMag* mag = new radTRecMag(center, size, M);

			// Add to group
			radThg hMag;
			hMag.rep = mag;
			group->AddElement(radHandle(), hMag);
		}

		cout << "  Created group with " << group->GroupMapOfHandlers.size() << " elements" << endl;
		return true;

	} catch(const exception& e) {
		cerr << "  [ERROR] Exception: " << e.what() << endl;
		return false;
	} catch(...) {
		cerr << "  [ERROR] Unknown exception" << endl;
		return false;
	}
}

//-------------------------------------------------------------------------
// Test 2: Create H-matrix field source
//-------------------------------------------------------------------------

bool test_create_hmatrix_source(radTGroup* group, radTHMatrixFieldSource*& hmat_source) {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 2: Create H-Matrix Field Source" << endl;
	cout << string(70, '-') << endl;

	if(!group) {
		cerr << "  [ERROR] NULL group pointer" << endl;
		return false;
	}

	try {
		// Configure H-matrix
		radTHMatrixConfig config;
		config.eps = 1e-6;
		config.max_rank = 50;
		config.min_cluster_size = 4;
		config.use_openmp = true;
		config.num_threads = 4;

		cout << "  Configuration:" << endl;
		cout << "    eps = " << config.eps << endl;
		cout << "    max_rank = " << config.max_rank << endl;
		cout << "    min_cluster_size = " << config.min_cluster_size << endl;
		cout << "    use_openmp = " << (config.use_openmp ? "true" : "false") << endl;

		// Create H-matrix field source
		hmat_source = new radTHMatrixFieldSource(group, config);

		cout << "  H-matrix field source created" << endl;
		cout << "  Number of elements: " << hmat_source->GetNumElements() << endl;

		return true;

	} catch(const exception& e) {
		cerr << "  [ERROR] Exception: " << e.what() << endl;
		return false;
	} catch(...) {
		cerr << "  [ERROR] Unknown exception" << endl;
		return false;
	}
}

//-------------------------------------------------------------------------
// Test 3: Build H-matrix
//-------------------------------------------------------------------------

bool test_build_hmatrix(radTHMatrixFieldSource* hmat_source) {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 3: Build H-Matrix" << endl;
	cout << string(70, '-') << endl;

	if(!hmat_source) {
		cerr << "  [ERROR] NULL H-matrix source pointer" << endl;
		return false;
	}

	try {
		auto t_start = chrono::high_resolution_clock::now();

		int result = hmat_source->BuildHMatrix();

		auto t_end = chrono::high_resolution_clock::now();
		double build_time = chrono::duration<double>(t_end - t_start).count();

		if(result != 1) {
			cerr << "  [ERROR] BuildHMatrix() returned " << result << endl;
			return false;
		}

		cout << "  H-matrix built successfully" << endl;
		cout << "  Build time: " << build_time << " seconds" << endl;
		cout << "  Memory usage: " << hmat_source->GetMemoryUsage() / 1024.0 / 1024.0 << " MB" << endl;
		cout << "  Is built: " << (hmat_source->IsBuilt() ? "Yes" : "No") << endl;

		return hmat_source->IsBuilt();

	} catch(const exception& e) {
		cerr << "  [ERROR] Exception: " << e.what() << endl;
		return false;
	} catch(...) {
		cerr << "  [ERROR] Unknown exception" << endl;
		return false;
	}
}

//-------------------------------------------------------------------------
// Test 4: Field calculation comparison
//-------------------------------------------------------------------------

bool test_field_calculation(radTGroup* group, radTHMatrixFieldSource* hmat_source) {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 4: Field Calculation Accuracy" << endl;
	cout << string(70, '-') << endl;

	if(!group || !hmat_source) {
		cerr << "  [ERROR] NULL pointer" << endl;
		return false;
	}

	try {
		// Test points (mm)
		vector<TVector3d> test_points = {
			TVector3d(0, 0, 100),      // Above first magnet
			TVector3d(250, 0, 100),    // Above middle
			TVector3d(500, 0, 100),    // Above last magnet
			TVector3d(250, 50, 0),     // Side of middle
			TVector3d(250, 0, 0)       // Center height
		};

		cout << "  Testing " << test_points.size() << " field evaluation points..." << endl;
		cout << fixed << setprecision(6);

		bool all_pass = true;

		for(size_t i = 0; i < test_points.size(); i++) {
			// Direct calculation (reference)
			radTField field_direct;
			field_direct.P = test_points[i];
			field_direct.FieldKey.B_ = 1;
			field_direct.B = TVector3d(0, 0, 0);

			group->B_comp(&field_direct);

			// H-matrix calculation
			radTField field_hmat;
			field_hmat.P = test_points[i];
			field_hmat.FieldKey.B_ = 1;
			field_hmat.B = TVector3d(0, 0, 0);

			hmat_source->B_comp(&field_hmat);

			// Compare
			TVector3d diff = field_direct.B - field_hmat.B;
			double error = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);

			cout << "\n  Point " << (i+1) << ": ("
			     << test_points[i].x << ", "
			     << test_points[i].y << ", "
			     << test_points[i].z << ") mm" << endl;
			cout << "    Direct:  B = ("
			     << field_direct.B.x << ", "
			     << field_direct.B.y << ", "
			     << field_direct.B.z << ") T" << endl;
			cout << "    H-matrix: B = ("
			     << field_hmat.B.x << ", "
			     << field_hmat.B.y << ", "
			     << field_hmat.B.z << ") T" << endl;
			cout << "    Error: " << error * 1e6 << " µT" << endl;

			if(error > 1e-6) {
				cout << "    [WARNING] Error exceeds tolerance" << endl;
				all_pass = false;
			} else {
				cout << "    [OK]" << endl;
			}
		}

		return all_pass;

	} catch(const exception& e) {
		cerr << "  [ERROR] Exception: " << e.what() << endl;
		return false;
	} catch(...) {
		cerr << "  [ERROR] Unknown exception" << endl;
		return false;
	}
}

//-------------------------------------------------------------------------
// Test 5: Performance comparison
//-------------------------------------------------------------------------

bool test_performance(radTGroup* group, radTHMatrixFieldSource* hmat_source) {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 5: Performance Comparison" << endl;
	cout << string(70, '-') << endl;

	if(!group || !hmat_source) {
		cerr << "  [ERROR] NULL pointer" << endl;
		return false;
	}

	try {
		// Create grid of test points
		int n_points = 1000;
		vector<TVector3d> test_points;
		test_points.reserve(n_points);

		for(int i = 0; i < n_points; i++) {
			double x = (i % 10) * 50.0;
			double y = ((i / 10) % 10) * 50.0;
			double z = (i / 100) * 50.0;
			test_points.push_back(TVector3d(x, y, z));
		}

		cout << "  Evaluating field at " << n_points << " points..." << endl;

		// Direct calculation timing
		auto t_start = chrono::high_resolution_clock::now();

		for(const auto& point : test_points) {
			radTField field;
			field.P = point;
			field.FieldKey.B_ = 1;
			field.B = TVector3d(0, 0, 0);
			group->B_comp(&field);
		}

		auto t_end = chrono::high_resolution_clock::now();
		double time_direct = chrono::duration<double>(t_end - t_start).count();

		cout << "  Direct calculation: " << time_direct << " seconds" << endl;

		// H-matrix calculation timing
		t_start = chrono::high_resolution_clock::now();

		for(const auto& point : test_points) {
			radTField field;
			field.P = point;
			field.FieldKey.B_ = 1;
			field.B = TVector3d(0, 0, 0);
			hmat_source->B_comp(&field);
		}

		t_end = chrono::high_resolution_clock::now();
		double time_hmat = chrono::duration<double>(t_end - t_start).count();

		cout << "  H-matrix calculation: " << time_hmat << " seconds" << endl;

		double speedup = time_direct / time_hmat;
		cout << "  Speedup: " << speedup << "x" << endl;

		// For small N (10 elements), direct might be faster
		// Just check that both methods work
		return true;

	} catch(const exception& e) {
		cerr << "  [ERROR] Exception: " << e.what() << endl;
		return false;
	} catch(...) {
		cerr << "  [ERROR] Unknown exception" << endl;
		return false;
	}
}

//-------------------------------------------------------------------------
// Test 6: Dump H-matrix information
//-------------------------------------------------------------------------

bool test_dump_info(radTHMatrixFieldSource* hmat_source) {
	cout << "\n" << string(70, '-') << endl;
	cout << "Test 6: H-Matrix Information Dump" << endl;
	cout << string(70, '-') << endl;

	if(!hmat_source) {
		cerr << "  [ERROR] NULL H-matrix source pointer" << endl;
		return false;
	}

	try {
		cout << "\n";
		hmat_source->Dump(cout, 0);
		cout << endl;

		return true;

	} catch(const exception& e) {
		cerr << "  [ERROR] Exception: " << e.what() << endl;
		return false;
	} catch(...) {
		cerr << "  [ERROR] Unknown exception" << endl;
		return false;
	}
}

//-------------------------------------------------------------------------
// Main test driver
//-------------------------------------------------------------------------

int main(int argc, char** argv) {
	cout << string(70, '=') << endl;
	cout << "radTHMatrixFieldSource Test Program (Phase 2 Validation)" << endl;
	cout << string(70, '=') << endl;
	cout << "\nThis test validates the H-matrix integration:" << endl;
	cout << "  - Geometry extraction from radTGroup" << endl;
	cout << "  - H-matrix construction with HACApK" << endl;
	cout << "  - Field calculation accuracy" << endl;
	cout << "  - Performance comparison" << endl;

	TestResult results;

	// Test objects
	radTGroup* group = nullptr;
	radTHMatrixFieldSource* hmat_source = nullptr;

	// Run tests
	results.report("Create Magnetic System", test_create_magnetic_system(group));
	results.report("Create H-Matrix Source", test_create_hmatrix_source(group, hmat_source));
	results.report("Build H-Matrix", test_build_hmatrix(hmat_source));
	results.report("Field Calculation Accuracy", test_field_calculation(group, hmat_source));
	results.report("Performance Comparison", test_performance(group, hmat_source));
	results.report("Information Dump", test_dump_info(hmat_source));

	// Cleanup
	if(hmat_source) {
		delete hmat_source;
	}
	if(group) {
		delete group;
	}

	// Print summary
	results.summary();

	return (results.failed == 0) ? 0 : 1;
}
