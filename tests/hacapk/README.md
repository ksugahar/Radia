# HACApK Integration Tests

This directory contains tests for integrating HACApK (Hierarchical Adaptive Cross Approximation with post-K) library with Radia.

## About HACApK

HACApK is a library for hierarchical matrices (H-matrices) specifically designed for boundary element methods (BEM). It uses Adaptive Cross Approximation (ACA) algorithm to efficiently represent dense matrices that arise in BEM calculations.

**Key Features:**
- Fast matrix-vector multiplication: O(N log N) instead of O(N²)
- Memory efficient: O(N log N) instead of O(N²)
- Suitable for large-scale BEM problems
- MPI parallel support

## Directory Structure

```
tests/hacapk/
├── README.md                    # This file
├── test_hacapk_basic.cpp        # Basic HACApK functionality test
├── test_hacapk_cluster.cpp      # Cluster tree generation test
├── test_hacapk_radia.cpp        # Radia + HACApK integration test
└── CMakeLists.txt               # CMake build configuration
```

## Building Tests

### Prerequisites

- GCC or Intel C/C++ compiler with OpenMP support
- MPI library (OpenMPI or Intel MPI)
- HACApK library (located in `../../2025_10_31_HaCapK/HACApK_LH-Cimplm`)

### Build Instructions

```bash
cd tests/hacapk
mkdir build
cd build
cmake ..
make
```

## Running Tests

### Basic HACApK Test

```bash
./test_hacapk_basic
```

Tests basic HACApK functionality:
- Cluster tree generation
- Leaf matrix creation
- Memory allocation/deallocation

### Cluster Tree Test

```bash
./test_hacapk_cluster
```

Tests hierarchical clustering algorithm:
- Bounding box calculation
- Binary tree generation
- Cluster organization

### Radia Integration Test

```bash
./test_hacapk_radia
```

Tests integration with Radia:
- Magnetic field calculation with H-matrix acceleration
- Comparison with direct Radia calculation
- Performance benchmarking

## Test Results

Expected output:
- All tests should report `[PASS]` status
- Performance tests show speedup compared to direct calculation
- Memory usage should be significantly reduced for large problems

## Integration Plan

1. **Phase 1:** Basic HACApK integration (current)
   - Standalone HACApK tests
   - C interface wrapping

2. **Phase 2:** Radia magnetic field acceleration
   - Replace direct field summation with H-matrix
   - Implement kernel function for magnetic field

3. **Phase 3:** Optimization
   - Performance tuning
   - Memory optimization
   - Parallel execution with OpenMP/MPI

## References

- HACApK paper: Ida et al., "Implementation of Hierarchical Matrices for BEM"
- Radia documentation: https://github.com/ochubar/Radia
- H-matrix theory: Hackbusch, "Hierarchical Matrices"

---
**Created:** 2025-11-07
**Author:** Radia Development Team
