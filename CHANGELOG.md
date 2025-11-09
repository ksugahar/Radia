# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

## [1.0.10] - 2025-11-10

### Fixed
- **H-matrix Implementation**
  - Implemented full matrix block storage in HACApK library
  - Implemented full matrix-vector multiplication
  - Fixed kernel function to use accurate B_comp() instead of dipole approximation
  - Achieved <1% accuracy for N=100 test (0.0119% max error)

### Changed
- **Test Suite**
  - Fixed all test return values to use assert instead of return (pytest compliance)
  - Fixed rad.Solve() return value checks (returns list, not int)
  - Added material to relaxation performance test
  - Adjusted transformation inversion test tolerance for numerical precision
  - All 76 tests now passing (100%)

### Removed
- Removed temporary development and debug files
- Removed old dipole approximation benchmark
- Cleaned up plan and status documentation files

### Added (2025-11-10)

- **H-matrix Benchmarks**
  - Added H-matrix field evaluation benchmark (examples/solver_benchmarks/)
  - Demonstrates <1% accuracy with 100+ magnetic elements

## [1.0.9] - 2025-11-02

### Changed (2025-11-02)

- **Package Name Change**
  - Renamed PyPI package from `radia` to `radia-ngsolve`
  - Version reset to 1.0.0 for new package
  - Updated all documentation with new package name
  - Installation: `pip install radia-ngsolve`

- **Build Scripts Migration**
  - Migrated from .cmd to PowerShell (.ps1) scripts
  - Unified build and upload into single `Publish_to_PyPI.ps1` script
  - Improved error handling and colored output
  - Updated documentation to reference new scripts

### Added (2025-11-02)

- **PyPI Distribution**
  - Published to PyPI as `radia-ngsolve`
  - Includes pre-built binaries for Windows Python 3.12
  - Complete LGPL-2.1 + Original RADIA BSD-style license
  - English documentation for international users

### Added (2025-11-01)

- **rad_ngsolve Unified Interface**
  - New unified `RadiaField` class supporting all field types
  - Field type selection: 'b' (flux density), 'h' (field), 'a' (vector potential), 'm' (magnetization)
  - Removed legacy interfaces (RadBfield, RadHfield, RadAfield) for cleaner API
  - Updated all examples and tests to use new interface

- **VTK Export Improvements**
  - Automatic mm → m unit conversion in `radia_vtk_export.py`
  - Consistent units across Radia (mm) and visualization tools (m)

- **Project Documentation**
  - Created `claude.md` with coding standards and project guidelines
  - Updated `.gitignore` to preserve `.pvsm` files and small example `.vtk` files
  - Consolidated build documentation in `README_BUILD.md`

- **Examples Cleanup**
  - Created two NGSolve integration example directories:
  - `examples/Radia_to_NGSolve_CoefficientFunction/` - Use Radia fields in NGSolve
  - `examples/NGSolve_CoefficientFunction_to_Radia_BackgroundField/` - Use background fields in Radia
  - Added `demo_field_types.py` demonstrating all field types
  - Removed obsolete documentation and test files
  - Reduced directory size from 11MB to 52KB

### Changed (2025-11-01)

- **Coding Standards**
  - Standardized on TAB characters for indentation (not 4 spaces)
  - Updated all Python and C++ files to follow new standards
  - Documented standards in `claude.md`

- **rad_ngsolve API Simplification**
  - Single unified interface: `RadiaField(obj, field_type)`
  - Removed backward compatibility layer
  - Cleaner, more maintainable codebase

### Fixed (2025-11-01)

- **Unit Conversion**
  - Fixed VTK export to properly convert mm → m
  - Consistent coordinate systems across Radia/NGSolve integration

### Added (2025-10-30)

- **Test Suite Reorganization**
  - Created `tests/` directory with standard Python project structure
  - Added comprehensive `tests/README.md` with detailed testing guide
  - Added `pytest.ini` for pytest configuration
  - Added `pyproject.toml` for modern Python project metadata
  - Added `.gitignore` with comprehensive ignore patterns
  - Moved all test files to `tests/` directory
  - Moved all benchmarks to `tests/benchmarks/` directory

- **Documentation**
  - Added `SECURITY_FIXES.md` documenting all security vulnerabilities and fixes
  - Added `docs/scripts/README.md` for development utility scripts
  - Updated main `README.md` with new test paths

### Fixed (2025-10-30)

- **Critical Security Vulnerabilities**
  - Fixed buffer overflow in `CombErStr` function (src/python/radpy.cpp:29-49)
  - Fixed array bounds overflow in `CopyPyStringToC` (src/python/pyparse.h:604)
  - Removed 43 unnecessary `Py_XINCREF` calls causing memory leaks
  - Fixed test suite material database issue (using valid 'Steel37' instead of invalid 'Iron')

### Changed (2025-10-30)

- **Test Organization**
  - Reorganized test files from project root to `tests/` directory
  - Updated import paths in all test files for new location
  - Test results improved from 5/7 (71.4%) to 7/7 (100%) passing

- **Build Artifacts**
  - Rebuilt `dist/radia.pyd` with all security fixes applied
  - Module size: 1.86 MB
  - Build: MSVC 19.44, Release mode, OpenMP enabled

## [4.32] - 2025-10-29

### Added
- OpenMP 2.0 parallelization for field computation
- PyVista viewer support for 3D visualization
- Comprehensive benchmark suite
- Performance reports and documentation

### Changed
- Migrated to Python 3.12 only (dropped Python 2.7, 3.6-3.11)
- Modernized build system with CMake
- Converted all indentation to tabs
- Removed legacy Igor Pro, Mathematica, GLUT, and MPI support

### Performance
- 2.7x speedup on 8-core systems for complex geometries
- OpenMP parallel field computation
- Optimized for Python 3.12

## Security Fixes Summary

### Critical (2025-10-30)
- **CVE-BUFFER-001**: Buffer overflow in error string concatenation
  - **Impact**: Potential arbitrary code execution
  - **Status**: ✅ Fixed

### High (2025-10-30)
- **CVE-BOUNDS-001**: Off-by-one array bounds error
  - **Impact**: Stack corruption when copying Python strings
  - **Status**: ✅ Fixed

### Medium (2025-10-30)
- **CVE-MEMORY-001**: Reference counting memory leaks (43 locations)
  - **Impact**: Memory leaks in long-running Python scripts
  - **Status**: ✅ Fixed

## Test Results

### Current (2025-11-01)
```
rad_ngsolve tests: 4/4 passed (100%)
radia core tests: 7/7 passed (100%)
```

## Links

- [Coding Standards](claude.md) (not in repository)
- [Security Fixes Documentation](SECURITY_FIXES.md)
- [Testing Guide](tests/README.md)
- [Build Instructions](README_BUILD.md)

---

**Maintained by**: Radia Development Team
**Python Version**: 3.12
**Last Updated**: 2025-11-01
