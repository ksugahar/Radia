# HACApK Tests - Quick Start Guide

This guide will help you build and run the HACApK integration tests.

## Prerequisites

### Windows (MSVC)

1. **Visual Studio 2019 or later** with C++ development tools
2. **CMake 3.15+** - Download from https://cmake.org/download/
3. **OpenMP support** - Included with Visual Studio

### Windows (MinGW/GCC)

1. **MinGW-w64** with GCC
2. **CMake 3.15+**
3. **OpenMP** - Usually included with GCC

## Quick Build (Windows)

### Method 1: Using build script (Recommended)

```cmd
cd S:\Radia\01_GitHub\tests\hacapk
build.cmd
```

This will:
- Configure CMake with Visual Studio generator
- Build in Release mode
- Create executables in `build/Release/`

### Method 2: Manual CMake

```cmd
cd S:\Radia\01_GitHub\tests\hacapk
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022" -A x64
cmake --build . --config Release
```

## Running Tests

### Test 1: Basic HACApK Functionality

```cmd
cd build\Release
test_hacapk_basic.exe
```

**Expected output:**
```
======================================================================
HACApK Basic Functionality Tests
======================================================================

----------------------------------------------------------------------
Test 1: Cluster Generation
----------------------------------------------------------------------
  Cluster created successfully
    ndim:   3
    nstrt:  0
    nsize:  64
    nnson:  2
    nmbr:   1
[PASS] Cluster Generation

----------------------------------------------------------------------
Test 2: Bounding Box Calculation
----------------------------------------------------------------------
  Bounding box calculated successfully
    X: [0, 30]
    Y: [0, 30]
    Z: [0, 30]
    Width: 51.9615
[PASS] Bounding Box Calculation

----------------------------------------------------------------------
Test 3: Binary Tree Generation
----------------------------------------------------------------------
  Binary tree generated successfully
    Depth:    4
    Clusters: 15
[PASS] Binary Tree Generation

======================================================================
Test Summary
======================================================================
Passed: 3
Failed: 0
Total:  3
======================================================================
```

### Test 2: Radia Integration Concept

```cmd
cd build\Release
test_hacapk_radia_concept.exe
```

**Expected output:**
```
======================================================================
HACApK + Radia Integration Tests (Conceptual)
======================================================================

======================================================================
Test: Magnetic Field Calculation
======================================================================

Creating circular coil:
  Radius:     100 mm
  Segments:   64
  Current:    1000 A

Field point: (0, 0, 50) mm

----------------------------------------------------------------------
Direct Calculation (Biot-Savart)
----------------------------------------------------------------------
  Bx = 0.000xxx mT
  By = 0.000xxx mT
  Bz = 2.xxxxx mT
  |B| = 2.xxxxx mT
  Time: 0.xxx ms

----------------------------------------------------------------------
Expected Performance for Large Problems
----------------------------------------------------------------------
  Current problem size: 64 elements
  For N=10,000 elements and M=1,000 field points:
    Direct:   O(N*M) = 10,000,000 operations
    H-matrix: O((N+M)*log(N+M)) ≈ 100,000 operations
    Expected speedup: ~100x
```

## Running All Tests with CTest

```cmd
cd build
ctest -C Release --verbose
```

This runs all registered tests automatically.

## Troubleshooting

### Error: "HACApK directory not found"

**Solution:** Check that the HACApK library is in the correct location:
```
S:\Radia\2025_10_31_HaCapK\HACApK_LH-Cimplm\
```

If it's in a different location, edit `CMakeLists.txt` line 20:
```cmake
set(HACAPK_DIR "path/to/your/HACApK/directory")
```

### Error: "OpenMP not found"

**Solution:**
- **Visual Studio:** Make sure C++ development tools are installed
- **MinGW:** Install GCC with OpenMP support (`-fopenmp` flag)

### Build warnings about unused variables

These are expected and can be ignored. They come from the HACApK C code.

### Link errors

If you get linker errors:
1. Make sure you're building in Release mode (not Debug)
2. Check that OpenMP is properly configured
3. Try cleaning and rebuilding: `build.cmd clean`

## Understanding the Tests

### test_hacapk_basic.cpp

Tests fundamental HACApK operations:
- **Cluster generation:** Creates hierarchical tree structure
- **Bounding box:** Calculates spatial bounds for clusters
- **Binary tree:** Builds adaptive tree based on point distribution

These are building blocks for H-matrix acceleration.

### test_hacapk_radia_concept.cpp

Demonstrates how HACApK would integrate with Radia:
- **Direct field calculation:** Standard Biot-Savart law (slow)
- **H-matrix concept:** How to accelerate with hierarchical matrices
- **Performance analysis:** Expected speedup for large problems

This is a **proof-of-concept**. Full implementation requires:
1. Kernel function for magnetic field interactions
2. ACA algorithm for matrix compression
3. Integration with Radia's field calculation pipeline

## Next Steps

1. ✅ **Run basic tests** - Verify HACApK library works
2. ✅ **Review concept test** - Understand integration approach
3. ⏳ **Implement kernel function** - Magnetic field calculation
4. ⏳ **Add ACA compression** - Use HACApK's adaptive algorithm
5. ⏳ **Benchmark performance** - Test on large coil systems
6. ⏳ **Integrate with Radia** - Replace direct field summation

## References

- **HACApK documentation:** See PDF in `S:\Radia\2025_10_31_HaCapK\`
- **Radia source:** `S:\Radia\01_GitHub\src\core\`
- **H-matrix theory:** Hackbusch, "Hierarchical Matrices" (book)

## Support

For issues or questions:
1. Check this QUICKSTART.md
2. Review README.md in this directory
3. Consult HACApK documentation
4. Review Radia source code

---
**Last Updated:** 2025-11-07
