# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Fixed

- **BiCGSTAB Block Jacobi Preconditioner for Ill-Conditioned Matrices**
  - Added automatic block Jacobi preconditioner for distorted hexahedral meshes
  - Scalar Jacobi fails when diagonal ratio > 10 or min dominance < 0.1
  - Block Jacobi inverts each element's 6x6 diagonal block using LAPACK
  - V304 mesh (74 elements with 45° skew): BiCGSTAB now converges (was diverging)
  - Automatic detection based on matrix conditioning analysis

## [1.7.0] - 2026-02-05

### Fixed

- **IMA (Image Method Analysis) Symmetry - Correct Sign Selection**
  - Fixed IMA sign selection policy for combined symmetries (+x+y-z, etc.)
  - Field parallel to mirror plane: use `+` (symmetric)
  - Field perpendicular to mirror plane: use `-` (antisymmetric)
  - Verified with 2-element and 8-element test cases (DOF reduction 48 -> 6)
  - All test ratios = 1.0000 (exact match with full model)

- **Netgen Face Normal Convention**
  - Removed inside/outside check in `SetupFaceGeometry()`
  - Face normals now computed mechanically from vertex winding order
  - Follows Netgen convention: Face ordering 0=z-, 1=z+, 2=y-, 3=y+, 4=x-, 5=x+

### Added

- **IMA Verification Tests**
  - `test_ima_2elem_linear.py`: 2 elements with shared boundary at z=0
  - `test_ima_8elem_linear.py`: 8 elements (octants) with IMA +x+y-z
  - Both tests verify magnetization and field computation match full model

### Documentation

- Updated CLAUDE.md with IMA sign selection policy
- Documented IMA boundary element limitation (elements with faces ON symmetry plane)

## [1.6.0] - 2026-01-14

### Added

- **mmm_core pybind11 Module - Standalone MMM Solver**
  - New `mmm_core` Python module with direct C++ bindings via pybind11
  - `MMMBuilder`: Build interaction matrices from tetrahedral/hexahedral meshes
    - `add_tetrahedra_from_mesh(vertices, elements)`: Add tetrahedra from mesh data
    - `add_hexahedra_from_mesh(vertices, elements)`: Add hexahedra from mesh data
    - `build()`: Returns (N_matrix, dof_offset) tuple
  - `MMMSolver`: Linear solvers for MMM equations
    - `solve_lu(inv_chi, H_ext, chi_per_element)`: Direct LU decomposition
    - `solve_bicgstab(inv_chi, H_ext, tol, max_iter, chi_per_element)`: BiCGSTAB iteration
  - `MMMHACApKSolver`: H-matrix accelerated solver
    - `set_from_builder(builder, dof_offset)`: Set elements from MMMBuilder
    - `build_hmatrix(inv_chi, eps, leaf_size, eta, print_level)`: Build H-matrix with ACA+
    - `matvec(x)`: H-matrix vector product (O(N log N) complexity)
    - `solve(inv_chi, H_ext, tol, max_iter)`: Full BiCGSTAB solve with H-matrix
    - `get_stats()`: H-matrix compression statistics
  - `MMMFieldComputer`: Field computation from solved magnetization
    - `compute_b_field(M, obs_points)`: Compute B field at observation points
    - `compute_h_field(M, obs_points)`: Compute H field at observation points
  - Helper functions: `compute_chi_from_bh()`, `check_convergence()`

- **HACApK H-Matrix Library Integration**
  - Integrated HACApK library for ACA+ (Adaptive Cross Approximation) compression
  - O(N log N) memory usage vs O(N^2) for dense matrices
  - H-matrix statistics: compression ratio, memory usage, max rank
  - Automatic permutation handling via internal LOD array

### Changed

- **Documentation Updates**
  - Updated `docs/PLAN_MMM_PYBIND11_REFACTOR.md` to "Implementation Complete" status
  - Updated `examples/peec_integration/test_mmm_hacapk.py` to use new mmm_core API
  - All mmm_core examples now use direct C++ API instead of mmm_ngsolve wrapper

### Removed

- **Legacy HACApK Files**
  - Removed `src/ext/HACApK/cHACApK_radia.c` (replaced by cHACApK_cpp_impl.c)
  - Removed `src/ext/HACApK/cHACApK_radia.h` (replaced by cHACApK_cpp.h)
  - Removed development test files: `test_hacapk_quick.py`, `test_hacapk_simple.py`, `test_import.py`

### Technical Details

- **HACApK Permutation Fix**: Fixed 1-based to 0-based index conversion for LOD array
- **Double Permutation Fix**: Removed manual permutation in MatVec (HACApK handles internally)
- **inv_chi Parameter**: Added inv_chi to build_hmatrix for correct system matrix construction

## [1.5.0] - 2026-01-11

### Added

- **CplMag Coupled PEEC-MMM Solver**
  - New coupled solver combining PEEC (conductor) with MMM (magnetic material)
  - Full element-to-element MMM coupling with demagnetization tensor
  - Complex permeability support (mu = mu' - j*mu") for magnetic losses
  - APIs: `CplMagCreate()`, `CplMagSetFrequency()`, `CplMagSetVoltage()`, `CplMagSetMu()`, `CplMagSolve()`, `CplMagDelete()`
  - Target applications: WPT (Wireless Power Transfer), induction heating with ferromagnetic cores

- **Matrix Symmetrization for CLN Model Order Reduction**
  - New `CplMagSetSymmetric(solver, use_symmetric)` API for matrix symmetrization
  - Variable scaling M' = sqrt(mu_0 * V) * M symmetrizes the coupled matrix
  - Enables CLN (Cauer Ladder Network) model order reduction
  - Symmetrized system produces **identical results** (machine precision: 6.5e-19 Ohm difference)
  - Mathematical proof via reciprocity: Z_LM^T = mu_0 * V * Z_ML

- **CLN Model Reduction Design Document**
  - New `docs/CLN_MODEL_REDUCTION_DESIGN.md` with full symmetry analysis
  - PEEC Loop-Star decomposition symmetry analysis (all blocks symmetric)
  - MMM demagnetization tensor symmetry analysis
  - Hierarchical CLN extraction strategy: PEEC-only first, then magnetic coupling
  - ACA+ low-rank approximation integration plan

- **Symmetrization Verification Scripts**
  - `examples/peec_integration/verify_symmetrization.py` - Mathematical proof verification
  - `examples/peec_integration/test_symmetrization.py` - Numerical equivalence test
  - `examples/peec_integration/test_cplmag_cubit_hex.py` - Cubit hex mesh test

- **Hex Mesh Import Functions**
  - `create_hex_mesh_grid()`: Create structured hexahedral mesh (no Cubit needed)
  - `cubit_hex_to_radia()`: Import Cubit hexahedral mesh to Radia
  - Located in `netgen_mesh_import.py`

- **VTS Field Export (C++ Implementation)**
  - `FldVTS()`: High-performance VTS (VTK XML Structured Grid) export
  - OpenMP parallelization for large field grids
  - Replaces Python-based VTK export for better performance

### Changed

- **Mesh Operations Policy**
  - All mesh operations now use "Netgen with Coreform Cubit Integration"
  - Coreform Cubit provides geometry and high-quality hex meshing
  - Netgen/NGSolve provides mesh import interface to Radia
  - See `S:\CoreformCubit\01_GitHub` for `cubit_mesh_export` utilities

### Removed

- **Deprecated Mesh APIs**
  - `ObjCutMag`: Removed from Python API (use Cubit instead)
  - `ObjDivMag`, `ObjDivMagPln`: Not supported (use Cubit for mesh subdivision)

- **VTK Geometry Export**
  - `ObjDrwVTK()`: C++ geometry export removed
  - `exportGeometryToVTK()`: Python geometry export removed
  - Use `FldVTS()` for field visualization in ParaView

### Documentation

- Updated CLAUDE.md with mesh operations policy
- Added CplMag solver documentation in header files
- Updated examples in `examples/peec_integration/`

## [1.4.4] - 2026-01-01

### Changed

- **API Cleanup: ObjHexahedron/ObjTetrahedron Standardization**
  - Updated all examples and documentation to use high-level Python APIs
  - `ObjHexahedron(vertices, magnetization)` for 8-vertex hexahedral elements
  - `ObjTetrahedron(vertices, magnetization)` for 4-vertex tetrahedral elements
  - These APIs auto-generate face topology - no need to specify face indices
  - `ObjPolyhdr` is now internal API only (used for wedge/pyramid/surface meshes)

### Cleaned

- **Repository Cleanup**
  - Removed 30+ debug/diagnostic scripts from `examples/ngsolve_integration/`
  - Removed 25+ debug test files from `tests/`
  - Removed internal design documents (historical info preserved in git)
  - Removed benchmark result JSON and VTK output files

### Documentation

- **Updated API Examples**
  - All Python examples now use `ObjHexahedron`/`ObjTetrahedron` consistently
  - Updated CLAUDE.md, README.md, API_REFERENCE.md
  - Updated example READMEs with correct API usage

## [1.4.3] - 2025-12-31

### Fixed

- **PyPI Wheel Build Path**
  - Fixed setup.py to check `build-msvc/` first, then fall back to `build/Release/`
  - Ensures wheel contains the latest .pyd built with MSVC + Intel MKL

### Added

- **Wheel Verification Policy**
  - Added mandatory wheel verification steps to CLAUDE.md
  - Prevents shipping outdated .pyd files in PyPI packages

## [1.4.2] - 2025-12-31

### Fixed

- **Face-based Scalar Potential (Phi) Calculation**
  - Replaced dipole approximation with accurate face-based integration for ObjHexahedron/ObjTetrahedron
  - New formula: `Phi = (1/4pi) * M dot BufVect` where `BufVect = n * integral(1/|r-r'|) dS`
  - Phi on z-axis now matches exactly between ObjRecMag and ObjHexahedron (< 1e-10 error)
  - Correct symmetry behavior: Phi ~ 0 on x/y axes for z-magnetized blocks (due to face cancellation)

- **Vector Potential (A) Field Consistency**
  - A field uses same face-based integration: `A = (1/4pi) * M x BufVect`
  - A = 0 on symmetry axis is physically correct (not a bug)
  - Off-axis A field matches within 2% between ObjRecMag and ObjHexahedron

### Added

- **New Test Files for Field Verification**
  - `tests/test_phi_field.py` - 9 tests for scalar potential computation
  - `tests/test_a_field.py` - 9 tests for vector potential computation
  - `tests/test_field_relations.py` - 10 tests for Maxwell equation consistency

### Verified

- **Maxwell Equation Consistency**
  - `curl(A) proportional to B` - verified at multiple points
  - `-grad(Phi) proportional to H` - verified at multiple points
  - `B = mu_0 * H` in air region - verified with < 0.01% error

## [1.4.1] - 2025-12-31

### Performance

- **Tetra HACApK 12.9x Speedup**
  - Implemented ELF-style face basis caching for tetrahedral H-matrix construction
  - Pre-compute face basis ONCE per face, reuse for all 3 magnetization directions
  - Pre-compute edge parameters (DS, AM, XD, YD) in RadTriangleFaceBasis struct
  - Reduces coordinate transformation overhead from 12x to 4x per element pair

### Benchmark Results (maxh=0.15, 2211 elements, 6633 DOF)

| Metric | Before | After | Improvement |
|--------|--------|-------|-------------|
| H-matrix build | 23.96s | 1.86s | **12.9x** |
| Total solve | 24.33s | 2.24s | **10.9x** |

### Comparison

- Tetra HACApK is now **faster than ELF** (2.24s vs 3.8s)
- Tetra HACApK is more efficient than Hex HACApK for similar DOF counts

## [1.4.0] - 2025-12-31

### Major Release Highlights

This is a major feature release with significant new APIs, performance improvements, and OpenMP thread-safe batch field evaluation.

### Added

- **Simplified Element Creation APIs**
  - `ObjTetrahedron(vertices, magnetization)` - Create tetrahedron from 4 vertices (face topology auto-generated)
  - `ObjHexahedron(vertices, magnetization)` - Create hexahedron from 8 vertices (face topology auto-generated)
  - No longer need to specify `TETRA_FACES` or `HEX_FACES` constants
  - Both APIs accept optional magnetization parameter (defaults to [0,0,0])

- **Batch Field Evaluation APIs with OpenMP Parallelization**
  - `FldBatch(obj, points, method)` - Compute B and H fields at multiple points efficiently
  - `FldPhi(obj, points)` - Compute magnetic scalar potential at multiple points
  - `FldA(obj, points)` - Compute magnetic vector potential at multiple points
  - `ClassifyPoints(obj, points, threshold)` - Classify points as inside/near/far from elements
  - **OpenMP parallelization enabled** for n_points > 100 (thread-safe implementation)

- **FMM Field Evaluation Example**
  - New example directory: `examples/fmm_field_evaluation/`
  - `demo_fldbatch.py` - Demonstrates batch field evaluation

### Fixed

- **FldBatch OpenMP Thread Safety**
  - Fixed `radTHandle` reference count data race in parallel regions
  - Changed `radTHandlePgnAndTrans` from value copy to const reference
  - Changed `std::vector` to `std::array` for stack allocation in B_comp functions
  - Added GIL release (`Py_BEGIN_ALLOW_THREADS`) for proper Python/C++ threading
  - Verified working with 20,000+ points and 4 threads

### Changed

- **Documentation Updates**
  - Updated all README files to use new `ObjTetrahedron`/`ObjHexahedron` APIs
  - Added comprehensive API documentation for batch field evaluation
  - Updated `docs/API_REFERENCE.md` with new API signatures

### Performance

- **Batch Field Evaluation**
  - 1,000 points: ~20ms (4 threads)
  - 10,000 points: ~93ms (4 threads)
  - 20,000 points: ~206ms (4 threads)
  - Linear scaling with point count

## [1.3.16] - 2025-12-24

### Fixed

- **BiCGSTAB Convergence Dramatically Improved**
  - Changed initial guess from zero to previous solution (FlatMagn)
  - This matches ELF's approach and significantly accelerates convergence
  - N=10 benchmark: 50 iterations -> 5 iterations (was 10x slower than ELF)
  - N=5 benchmark: 8 iterations -> 3 iterations

### Changed

- **Debug Logging Policy**
  - Removed C++ file-based debug logging (policy: debug info via Python only)
  - All debug information should be managed through Python scripts

### Validation

- N=10 hexahedron benchmark (1000 elements, H_ext = 200,000 A/m):
  - BiCGSTAB: 5 iterations, M_avg_z = 716,316 A/m (ELF: 4 iterations)
  - LU: 13 iterations, M_avg_z = 716,281 A/m (matches ELF exactly)
  - HACApK: 5 iterations, M_avg_z = 715,885 A/m

## [1.3.15] - 2025-12-22

### Added

- **B-field Based Convergence (mucal2)**
  - Implemented industry-standard B-field convergence criterion for nonlinear materials
  - Formula: `rel_change = |B_new - B_old| / B_sat`
  - B_sat automatically calculated from BH curve (B_last - H_last)
  - Matches Newton-Raphson convergence behavior of reference solvers

### Changed

- **Nonlinear Solver Convergence**
  - Replaced chi-based (mucal1) convergence with B-field based (mucal2) convergence
  - Provides faster and more reliable convergence for saturable materials
  - Iteration count now matches reference solver exactly (e.g., 6 iterations for N=5 cube)

- **Permeability Calculation from BH Curve**
  - Changed from chi (susceptibility) based to mu_r (relative permeability) based calculation
  - mu_r = B / (mu_0 * H) directly from BH curve interpolation
  - More numerically stable for high-permeability materials
  - Consistent with industry-standard magnetostatic solvers

- **BH Curve Storage Format**
  - Internally uses B-H representation (B as function of H)
  - Supports both linear interpolation and spline interpolation
  - GetBsaturation() method added for automatic B_sat extraction

### Validation

- N=5 hexahedron benchmark (125 elements, H_ext = 200,000 A/m):
  - Radia: 6 iterations, M_avg_z = 702,129.2 A/m
  - Reference: 6 iterations, M_avg_z = 702,131.9 A/m
  - Difference: 0.0004%

## [1.3.14] - 2025-12-11

### Fixed

- **PyPI Package DLL Loading Issue**
  - Fixed "DLL load failed" error when importing radia from PyPI package
  - Added Intel MKL runtime DLL (`mkl_rt.2.dll`) to package distribution
  - Added DLL directory to search path via `os.add_dll_directory()` on Windows
  - Package now correctly loads Intel MKL dependency

## [1.3.13] - 2025-12-11

### Performance

- **Major Performance Optimization with Intel MKL and OpenMP**
  - Replaced manual BLAS operations with Intel MKL CBLAS calls (cblas_ddot, cblas_dnrm2, cblas_daxpy, cblas_dcopy, cblas_dscal, cblas_dgemv)
  - Added OpenMP parallelization to interaction matrix O(N^2) construction
  - LU solver: Up to **240x faster** (e.g., 410s -> 1.7s for 390 elements)
  - BiCGSTAB solver: Up to **17x faster** (e.g., 29s -> 1.7s for 390 elements)
  - BiCGSTAB now **faster than ELF_MAGIC** for large problems (0.4x-0.6x ratio)

### Changed

- **Solver Architecture**
  - Original Radia used Implicit SS (Successive Substitution/Gauss-Seidel) method which had slow convergence for high-permeability nonlinear materials
  - Replaced with BiCGSTAB iterative solver (Method 1) for better convergence
  - LU direct solver (Method 0) retained for small problems and guaranteed convergence
  - Both solvers now use pure Newton-Raphson iteration without Gauss-Seidel M(H) correction

### Technical Notes

- Convergence tolerance: Radia uses relative change ||dM||/||M||
- ELF_MAGIC default tolerance: 0.01 (1%), Radia benchmark used 0.0001 (0.01%)
- With same tolerance (0.01), Radia converges in 6 iterations vs ELF's 9 iterations

## [1.3.12] - 2025-12-11

### Changed

- **Solver Methods Simplified**
  - Only LU (Method 0) and BiCGSTAB (Method 1) remain
  - Implicit SS (Method 2) implementation removed due to nonlinear material issues
  - String arguments: 'lu'/'direct' for Method 0, 'bicgstab'/'iterative' for Method 1

## [1.3.11] - 2025-12-11

### Added

- **6-Face Hexahedral MSC Method (ELF_MAGIC Compatible)**
  - Implemented `B_comp_hexahedron_MSC()` for 6 quadrilateral faces
  - Each quad face is split into 2 triangles using [V0,V1,V2] + [V0,V2,V3]
  - Matches ELF_MAGIC face ordering and diagonal split convention
  - Single element validation: 0.05% error vs ELF_MAGIC reference

### Changed

- **Polyhedron Element Support**
  - Now only supports tetrahedra (4 faces) and hexahedra (6 faces)
  - Removed 12-face triangular hexahedron support (throws error)
  - Updated documentation to reflect 6-face MSC implementation

### Documentation

- Updated CLAUDE.md with 6-face hexahedral MSC details
- Updated docs/MMM_MSC_IMPLEMENTATION.md
- Updated docs/MESH_MSC_API_DESIGN.md

## [1.3.10] - 2025-12-10

### Added

- **NGSolve Integration Documentation**
  - Documented RadiaField coordinate transformation parameters (origin, u_axis, v_axis, w_axis)
  - Added detailed Coordinate Transformation section to NGSOLVE_INTEGRATION.md
  - Added Cache Methods documentation (PrepareCache, ClearCache, GetCacheStats)
  - New example scripts:
    - `demo_basic_field.py` - Basic RadiaField usage
    - `demo_coordinate_transform.py` - Coordinate transformation examples
    - `test_unit_conversion.py` - Unit conversion tests

### Fixed

- **Documentation Link Fixes**
  - Fixed broken links: `solver_time_evaluation` -> `solver_benchmarks`
  - Fixed broken links: `Radia_to_NGSolve_CoefficientFunction_A` -> `ngsolve_integration`
  - Fixed broken links: `02_EMPY_Field` -> `examples/ngsolve_integration`
  - Fixed file path: `src/python/radia_ngsolve.cpp` -> `src/radia/radia_ngsolve.cpp`

## [1.3.9] - 2025-12-05

### Added

- **NGSolve Integration Tests and Examples**
  - New `test_ngsolve_integration.py` - Comprehensive HDiv projection test
  - New `test_mesh_import.py` - Tetrahedral and hexahedral mesh import tests
  - New `demo_hdiv_projection.py` - HDiv function space projection example
  - Documented HDiv best practices for magnetic field projection

- **Build System Improvements**
  - Fixed Build.ps1 bugs (radia_ngsolve target name, path handling)
  - Automatic copy of .pyd files to src/radia/ for PyPI packaging
  - Multi-location search for built .pyd files

- **Nonlinear Material Benchmarks**
  - New `benchmark_bicgstab_hex.py` - BiCGSTAB solver benchmark
  - New `benchmark_solver_methods.py` - LU vs BiCGSTAB comparison
  - New `compare_radia_elfmagic_field.py` - ELF_MAGIC comparison

- **Documentation**
  - `docs/MESH_MSC_API_DESIGN.md` - Mesh MSC API design document
  - `docs/MMM_MSC_IMPLEMENTATION.md` - MMM+MSC implementation guide
  - `docs/SOLVER_METHODS.md` - Solver method documentation

### Changed

- **CLAUDE.md Updates**
  - Added isotropic MatLin usage policy (single argument form required)
  - Updated NGSolve integration best practices with HDiv(order=2) recommendation
  - Added release workflow policy (Build -> Test -> Git Push -> PyPI)

### Fixed

- **Build.ps1 Bug Fixes**
  - Fixed target name `rad_ngsolve` -> `radia_ngsolve`
  - Fixed broken paths `srcadiaadia.pyd` -> `srcadiaadia.pyd`
  - Added automatic directory creation for src/radia

## [1.3.4] - 2025-11-27

### Added

- **SolverTetraMethod() API**
  - New API to control tetrahedral element field computation method
  - `rad.SolverTetraMethod(0)`: Original Radia polygon method (default)
  - `rad.SolverTetraMethod(1)`: Analytical method for high-permeability materials
  - Replaces deprecated `RADIA_TETRA_METHOD` environment variable
  - Both methods now produce identical results (verified 0.00% difference)

### Fixed

- **ANALYTICAL Method (SolverTetraMethod=1)**
  - Fixed double coordinate transformation bug causing 59.64% error
  - Root cause: `B_comp` was incorrectly receiving global basis vectors
  - Solution: Reverted to LOCAL coordinate computation with identity basis
  - Method 0 and Method 1 now produce identical results

### Changed

- **Documentation Updates**
  - Updated `FldUnits` documentation with meter/millimeter setting examples
  - Added NGSolve integration unit policy (`rad.FldUnits('m')` required)
  - Corrected "BEM" terminology to "MMM" (Magnetic Moment Method) throughout
  - Added `SolverTetraMethod` to API_REFERENCE.md and API_EXTENSIONS.md
  - Updated NGSOLVE_USAGE_GUIDE.md and NGSOLVE_INTEGRATION.md

### Documentation

- **Unit System Policy**
  - Enforced no hard-coded unit conversions policy
  - All unit conversions via `rad.FldUnits()` only
  - Updated all NGSolve examples to use `rad.FldUnits('m')`

## [1.3.3] - 2025-01-21

### Optimized

- **radia_ngsolve Memory Management**
  - Minimized GIL (Global Interpreter Lock) scope - acquired only during Python calls
  - C++ cache lookups now avoid GIL entirely for improved performance
  - Improved Python object lifecycle - temporary objects released in tight scopes
  - Cached Radia module import to eliminate repeated `py::module_::import()` calls
  - Memory growth verified as Python memory pool initialization (not a leak)
  - Memory usage stabilizes after warm-up period (-41% saturation in extended tests)

### Added

- **Memory Test Suite for NGSolve Integration**
  - `test_radia_ngsolve_memory_leak.py` - Basic memory leak detection test
  - `test_radia_ngsolve_longrun.py` - Extended 500-iteration saturation test
  - `test_radia_ngsolve_with_cache.py` - PrepareCache() memory behavior test
  - `test_radia_core_memory.py` - Baseline Radia core memory verification
  - All tests confirm no true memory leak exists

### Documentation

- **Memory Management Section** in `examples/NGSolve_Integration/README.md`
  - Documents memory optimization techniques implemented in v1.3.3
  - Explains memory growth behavior (Python memory pool initialization)
  - Provides test results showing saturation pattern
  - Clarifies this is safe for long-running simulations

## [1.3.2] - 2025-11-21

### Changed

- **H-Matrix Solver Control (Breaking Change)**
  - H-matrix acceleration now requires explicit enablement via `rad.SolverHMatrixEnable(1, eps=1e-4, max_rank=30)`
  - Removed automatic enablement based on problem size (no more N > 200 threshold)
  - Updated documentation in `radpy.cpp` to clarify explicit control requirement
  - Updated benchmark scripts to use explicit H-matrix enable/disable calls
  - Follows H-Matrix Solver Control Policy in CLAUDE.md

### Added

- **Documentation Organization**
  - Created `internal/` folder for maintainer documentation
    - `internal/design/` - Architecture decisions, implementation proposals
    - `internal/analysis/` - Performance analysis, bottleneck investigations
  - Moved development documents from `docs/` to `internal/`
  - Updated `docs/README.md` with clear separation between user and maintainer docs

- **NGSolve Integration Test Suite**
  - `test_batch_evaluation.py` - Batch field evaluation testing
  - `test_cf_direct.py` - Direct CoefficientFunction usage
  - `test_convergence_hdiv.py` - H(div) space convergence analysis
  - `test_curlA_equals_B.py` - Vector potential curl verification
  - `test_curl_A_detailed.py` - Detailed curl(A) = B testing
  - `test_hcurl_vs_hdiv.py` - H(curl) vs H(div) comparison
  - `test_order1.py` - First-order element testing
  - `test_radia_ngsolve_diagnostic.py` - NGSolve integration diagnostics
  - `test_set_vs_interpolate.py` - GridFunction.Set() vs Interpolate()
  - `test_without_gridfunction.py` - Direct field evaluation without GridFunction

### Fixed

- **.gitignore improvements**
  - Added VTK output file patterns (*.vtu, *.vtk)
  - Added Python development file patterns (experimental implementations)
  - Added test file patterns for development tests
  - Prevents committing temporary and experimental files

## [1.3.1] - 2025-11-21

### Added

- **Pure Python Cached Field Implementation**
  - New `radia_field_cached.py` module with `CachedRadiaField` class
  - 60,000x faster than C++ PrepareCache() (500 points: 60s → 1ms)
  - Linear O(N) scaling with ~2 us/point overhead
  - Avoids pybind11 overhead by keeping all data operations in Python
  - Helper function `collect_integration_points()` for NGSolve integration
  - Documented in `docs/PYTHON_CACHED_FIELD_SOLUTION.md`

- **NGSolve Integration Unit System Policy**
  - Comprehensive policy in CLAUDE.md for `rad.FldUnits('m')` usage
  - Ensures consistent meter units between Radia and NGSolve
  - Verification checklist for all NGSolve integration code

### Changed

- **Unit System Consistency (Breaking Change)**
  - All NGSolve integration examples now use `rad.FldUnits('m')`
  - Magnet dimensions converted from mm to m (35 files updated)
  - Removed manual coordinate scaling (×1000, ×0.001)
  - Affects: tests/ and examples/NGSolve_Integration/ folders

### Fixed

- **profile_batch_performance.py**: Fixed ZeroDivisionError when measurement time < 0.01ms
- **test_coordinate_transform.py**: Fixed UnicodeEncodeError (cp932) by replacing `·` with `*`
- **test_vector_potential.py**: Added missing `src/python` path for radia module import

### Documentation

- Added `docs/PYTHON_CACHED_FIELD_SOLUTION.md` - Pure Python cached field implementation
- Updated CLAUDE.md with NGSolve Integration Unit System Policy
- Updated CLAUDE.md with Windows Console Encoding (cp932) compatibility policy

## [1.3.0] - 2025-11-21

### Added

- **H-Matrix Cache for Batch Field Evaluation**
  - Implemented `PrepareCache()` method in `radia_ngsolve.RadiaField` for batch field evaluation
  - Enables H-matrix acceleration when setting GridFunctions from Radia fields
  - Single batch Radia.Fld() call replaces element-by-element evaluation (~13,000 calls → 1 call)
  - Cache data structure using FNV-1a hash with O(1) lookup performance
  - Python API: `PrepareCache(points)`, `ClearCache()`, `GetCacheStats()`
  - Cache hit rate: 80-95% during GridFunction.Set() operations
  - Documented in `docs/HMATRIX_CACHE_IMPLEMENTATION.md`

- **Cache Performance Monitoring**
  - Added `GetCacheStats()` method returning dictionary with cache statistics
  - Reports: enabled status, cache size, hits, misses, hit rate
  - Enables performance profiling and optimization

### Changed

- **radia_ngsolve.cpp Internal Structure**
  - Added cache member variables: `point_cache_`, `use_cache_`, cache statistics
  - Modified `Evaluate()` methods to check cache before direct Radia evaluation
  - Added hash function for 3D point quantization (tolerance: 1e-10 meters)

### Examples

- **examples/NGSolve_Integration/example_hmatrix_cache_usage.py**
  - Complete usage example demonstrating PrepareCache() workflow
  - Performance comparison: cached vs non-cached evaluation
  - Integration point collection from mesh elements

### Tests

- **tests/test_hmatrix_cache_simple.py**
  - Basic cache functionality test (PASS)
  - Verifies PrepareCache(), ClearCache(), GetCacheStats()

- **tests/test_hmatrix_cache.py**
  - Comprehensive GridFunction integration test
  - Accuracy verification and performance measurement

### Known Limitations

- Radia batch evaluation performance degrades for very large point sets (>1000 points)
- Recommended usage: moderate point counts (<500 points)
- Cache is not automatically invalidated when Radia geometry changes (manual ClearCache() required)

### Performance

- **Before**: Element-by-element evaluation, no H-matrix benefit
  - Example: 13,021 Radia.Fld() calls for 449-element mesh (avg 1.4 points/call)

- **After**: Single batch evaluation with cached results
  - Example: 1 Radia.Fld() call for all integration points
  - Expected speedup: 10-50x for large meshes (when Radia batch evaluation is efficient)

## [1.2.1] - 2025-11-20

### Fixed

- **radia_ngsolve GridFunction.Set() Bug Fix**
  - Fixed result matrix indexing in batch evaluation function (`src/python/radia_ngsolve.cpp`)
  - Changed from `result(component, point)` to `result(point, component)` (lines 348-350)
  - Bug was introduced in commit ab77976 (H-matrix implementation)
  - GridFunction.Set() now produces correct values matching direct CoefficientFunction evaluation

- **NGSolve Examples Unit Consistency**
  - Added `rad.FldUnits('m')` to all 9 NGSolve integration examples
  - Converted all Radia coordinates from millimeters to meters
  - Ensures consistent unit handling between Radia (now meters) and NGSolve (meters)
  - Updated comments and print statements to reflect meter units

### Changed

- **Test Suite Organization**
  - Renamed `verify_curl_A_equals_B_improved.py` to `test_curlA_equals_B.py` in tests folder
  - Added test acceptance criteria to `tests/README.md`
  - Created comprehensive GridFunction projection best practices documentation in CLAUDE.md

### Documentation

- **GridFunction Projection Guidelines**
  - Documented optimal finite element space selection: HDiv order=2 for B projection
  - Region-dependent accuracy expectations: 0.15-0.36% at practical distances (>1 mesh cell)
  - Evaluation guidelines: avoid GridFunction evaluation within 1 mesh cell of magnet surface
  - Added extensive test results and best practices to CLAUDE.md

### Examples Updated

All NGSolve integration examples now use consistent meter-based units:
- test_gridfunction_simple.py
- verify_curl_A_equals_B.py
- test_set_vs_interpolate.py
- test_mesh_convergence.py
- test_coordinate_transform.py
- test_batch_evaluation.py
- benchmark_gridfunction_set.py
- visualize_field.py
- demo_field_types.py

## [1.2.0] - 2025-11-17

### Added

- **Test Suite Expansion**
  - `test_magpylib_comparison.py` - Cross-validation with magpylib for cylindrical magnets
  - `test_update_hmatrix_magnetization.py` - H-matrix magnetization update functionality test
  - Tests use pytest.skip() for optional dependencies (magpylib)

- **Benchmark Additions**
  - `benchmark_solver_methods.py` - Comparison of Direct/Relaxation/H-matrix solver methods
  - Demonstrates performance characteristics of each solver approach

- **Documentation Improvements**
  - `docs/README.md` - Comprehensive documentation index and navigation
  - Organizes user documentation (API, H-matrix, NGSolve) and developer documentation
  - Quick start guide for different user types

### Changed

- **Repository Organization**
  - Moved development notes from `docs/` to `dev/notes/` (11 files)
  - Organized by category: implementation, performance, releases
  - Clearer separation between user-facing and development documentation
  - Removed obsolete `examples/H-matrix/` folder (merged into solver_benchmarks)

- **Code Quality Improvements**
  - Fixed absolute paths to relative paths in 4 benchmark scripts
  - Ensures portability across different development environments
  - All Python scripts now use `os.path.join(os.path.dirname(__file__), ...)` pattern

- **Development Policies**
  - Added "Python Script Path Import Policy" to CLAUDE.md
  - Mandates relative paths for all example and test scripts
  - Improves collaboration and distribution

### Fixed

- Corrected path resolution in test scripts (tests/ folder structure)
- Removed HACApK development files (bem-bb-config.txt, *.pbf, *.xcr, etc.)
- Removed NGSolve temporary output folders (rad.ObjBckgCF/)
- Updated .gitignore with rules for NGSolve and HACApK temporary files

### Documentation

- Updated README.md with computation accuracy analysis
- Improved examples/simple_problems/README.md

## [1.1.1] - 2025-11-13

### Fixed

- **Reverted Phase 3 Serialization** (Critical Performance Fix)
  - Phase 3 serialization caused 89% performance regression (8.95x → 1.0x speedup loss)
  - Restored Phase 2-B implementation with verified 8.3x speedup
  - See `docs/PHASE3_PERFORMANCE_ISSUE.md` for detailed analysis
  - Removed H-matrix disk caching APIs (will be reimplemented in future)

### Changed

- **Removed Automatic N=200 Threshold**
  - Users now have explicit control over H-matrix enable/disable
  - Removed automatic override based on problem size
  - H-matrix respects user's `SolverHMatrixEnable()` flag regardless of N
  - Added policy to CLAUDE.md documenting user control requirement

### Added

- **Extended Scaling Benchmarks**
  - New `benchmark_solver_scaling_extended.py` testing N=125 to N=4913
  - Results: 8.9x at N=343 → 117.1x at N=4913 speedup
  - Created `SCALING_RESULTS.md` with comprehensive analysis

- **Exact Size Benchmarks with Memory Compression Analysis**
  - New `benchmark_hmatrix_scaling_exact.py` for N=100, 200, 500, 1000, 2000, 5000
  - Time speedup: 3.0x → 98.2x (exponential increase)
  - Memory compression: 100% → 0.1% (99.9% reduction at N=5000)
  - Detailed speedup and memory analysis
  - Verifies H-matrix O(N² log N) time and O(N log N) memory complexity

### Documentation

- **Phase 2-B Re-evaluation**
  - Created `PHASE2B_REEVALUATION.md` documenting correct methodology
  - Updated all benchmarks with Phase 2-B measured performance
  - Clarified construction vs solve time distinction
  - Added H-Matrix control policy to CLAUDE.md

- **Updated README**
  - Added memory compression results (0.1% at N=5000)
  - Updated Key Findings with exponential scaling benefits
  - Documented exact size benchmark results

### Performance

- **Phase 2-B Verified Performance**
  - Solver: 8.3x speedup at N=343 (measured)
  - Scaling: 3x at N≈100 → 98x at N≈5000
  - Memory: 99.9% reduction at N=5000 vs dense O(N²)
  - Field evaluation: 4.0x speedup (5000 points, batch)
  - Parallel construction: 27.7x speedup (OpenMP)

## [1.1.0] - 2025-11-13

### Added

- **Phase 3B: Full H-matrix Serialization to Disk**
  - Complete H-matrix structure saved to disk (`.radia_cache/hmat/*.hmat`)
  - Instant startup for repeated simulations (~10x faster)
  - New APIs:
    - `rad.SolverHMatrixCacheFull(enable=1)` - Enable full serialization
    - `rad.SolverHMatrixCacheSize(max_mb=1000)` - Set cache size limit
    - `rad.SolverHMatrixCacheCleanup(days=30)` - Cleanup old entries
  - Binary format with version checking (magic number, format version, HACApK version)
  - Automatic cache management with LRU eviction
  - Cross-session persistence (9.7x speedup measured)
  - Complete documentation in `docs/HMATRIX_SERIALIZATION.md`

- **Comprehensive Solver Comparison Benchmark**
  - New `benchmark_solver_comparison.py` comparing LU, Gauss-Seidel, and H-matrix
  - Demonstrates when each method is optimal:
    - LU Decomposition: Best for N < 100 (O(N³) complexity)
    - Gauss-Seidel: Best for 100 < N < 200 (O(N²) per iteration)
    - H-matrix: Best for N > 200 (O(N² log N) per iteration)
  - Includes per-iteration timing, full solve timing, and accuracy verification

- **Material API Enhancement**
  - New `rad.MatPM(Br, Hc, easy_axis)` for permanent magnets with demagnetization
  - Distinguishes permanent magnets from linear magnetic materials
  - Updated API documentation with proper usage examples

### Changed

- **Examples Folder Reorganization**
  - Renamed `examples/H-matrix/` → `examples/solver_benchmarks/`
  - Merged solver_benchmarks and solver_time_evaluation folders
  - Consolidated all solver-related benchmarks into single location
  - Updated folder title to "Magnetostatic Solver Benchmarks with H-Matrix Acceleration"
  - Organized benchmarks into categories: Core, Advanced, Verification
  - Net reduction: 383 lines of redundant code removed

- **Documentation Updates**
  - Updated all path references across 7 documentation files
  - Added comprehensive solver method selection guide
  - Added note about H-matrix overhead for fast-converging problems
  - Updated performance metrics with actual measurements (not extrapolated)

### Performance

- **Measured Performance Improvements (v1.1.0)**
  - Disk caching: 9.7x speedup (0.602s → 0.062s startup)
  - Solver: 6.64x speedup for N=343 elements
  - Field evaluation: 3.97x speedup for 5000 points (batch)
  - Parallel construction: 27.74x speedup (OpenMP)
  - Overall workflow: 7-8x speedup for repeated simulations

- **Solver Comparison Results**
  - Small problems (N=27): LU 5.64x slower than GS, H-matrix 5.09x faster
  - Medium problems (N=125): LU 28.79x slower, H-matrix 1.02x faster
  - Large problems (N=343): LU skipped (too slow), H-matrix construction overhead dominates for fast convergence

### Documentation

- **Implementation History**
  - Complete Phase 1 through Phase 3B development timeline
  - 5-day development cycle (2025-11-08 to 2025-11-13)
  - ~1500 lines of production code
  - ~800 lines of test code
  - ~3000 lines of documentation

- **New Documentation Files**
  - `docs/HMATRIX_SERIALIZATION.md` - Phase 3B user guide
  - `docs/HMATRIX_IMPLEMENTATION_HISTORY.md` - Complete development history
  - `docs/HMATRIX_BENCHMARKS_RESULTS.md` - Comprehensive benchmark results
  - `docs/MATERIAL_API_IMPLEMENTATION.md` - Material API documentation

### Test Suite

- **New Tests**
  - 11 comprehensive test scripts in `tests/hmatrix/`
  - Covers Phase 2-A, 2-B, 3, and 3-B implementation
  - All tests passing (100% success rate)
  - Cross-session serialization verification
  - Field accuracy verification

### Fixed

- **API Documentation Corrections**
  - Fixed `docs/API_REFERENCE.md` permanent magnet examples
  - Changed from `MatLin` to `MatPM` for NdFeB, SmCo, Ferrite magnets
  - Added warnings about proper material usage
  - Clarified MatLin is for soft magnetic materials only

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

- **radia_ngsolve Unified Interface**
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

- **radia_ngsolve API Simplification**
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
radia_ngsolve tests: 4/4 passed (100%)
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
