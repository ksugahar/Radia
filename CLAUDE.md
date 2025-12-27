# Claude Code - Radia Project Development Guidelines

This document contains development guidelines and refactoring policies for the Radia project when working with Claude Code.

## Build Policy: MSVC + Intel MKL

### Compiler Requirement (2025-12-27)

**CRITICAL**: Use **MSVC (Visual Studio C++ compiler)** with **Intel MKL** for building Radia.

**Policy**:
- **Use MSVC** for building both `radia.pyd` and `radia_ngsolve.pyd`
- **Use Intel MKL** for BLAS/LAPACK operations (mkl_rt.lib)
- **radia_ngsolve REQUIRES MSVC** due to ABI compatibility with MSVC-compiled NGSolve libraries
- Intel oneAPI compiler (icx-cl) is NOT compatible with NGSolve linking

**Build Command**:
```powershell
powershell.exe -ExecutionPolicy Bypass -File "BuildMSVC.ps1"
# Or for clean rebuild:
powershell.exe -ExecutionPolicy Bypass -File "BuildMSVC.ps1" -Rebuild
```

**Why MSVC instead of Intel Compiler**:
1. **NGSolve Compatibility**: NGSolve is compiled with MSVC; Intel compiler produces LLVM bitcode objects that are incompatible with MSVC-compiled libraries
2. **radia_ngsolve Linking**: The `add_ngsolve_python_module` CMake function requires MSVC ABI compatibility
3. **Intel MKL Still Used**: mkl_rt.lib works with both MSVC and Intel compilers, providing fast BLAS/LAPACK

**Rationale**:
1. **NGSolve Integration**: Required for RadiaField CoefficientFunction support
2. **MKL Performance**: Intel MKL provides optimized BLAS/LAPACK regardless of compiler
3. **Build Simplicity**: Single build produces both radia.pyd and radia_ngsolve.pyd

**Required Software**:
- Visual Studio 2022 (MSVC compiler)
- Intel oneAPI Base Toolkit (for Intel MKL only, NOT the compiler)

---

## BLAS/LAPACK Policy: Intel MKL Only

### OpenBLAS Dropped (2025-12-28)

**Policy**: Radia uses **Intel MKL only** for BLAS/LAPACK operations. **OpenBLAS is NOT supported**.

**Rationale**:
1. **Performance**: Intel MKL is faster than OpenBLAS on Intel CPUs
2. **OpenMP Integration**: Intel MKL uses Intel OpenMP (libiomp5md.dll) which integrates well with parallel Radia code
3. **Simplicity**: Single BLAS library reduces complexity and potential conflicts

**Removed Files** (2025-12-28):
- `src/ext/openblas/` - OpenBLAS headers, libraries, and binaries
- `src/radia/libopenblas.dll` - OpenBLAS runtime DLL

**Required Intel MKL DLLs** (auto-copied by BuildMSVC.ps1):
- `mkl_rt.*.dll` - MKL SDL (Single Dynamic Library) runtime
- `mkl_core.*.dll` - MKL core
- `mkl_intel_thread.*.dll` - MKL threading
- `mkl_def.*.dll`, `mkl_avx2.*.dll` - CPU kernels
- `mkl_vml_def.*.dll`, `mkl_vml_avx2.*.dll` - Vector math library
- `libiomp5md.dll` - Intel OpenMP runtime
- `libmmd.dll`, `svml_dispmd.dll` - Intel compiler runtime

**Note**: The DLL patterns use wildcards (e.g., `mkl_rt.*.dll`) for version-agnostic compatibility with future Intel oneAPI releases.

---

## OpenMP Policy: Intel OpenMP Only

### Use Intel OpenMP Instead of MSVC OpenMP (2025-12-28)

**Policy**: Radia uses **Intel OpenMP (libiomp5md.dll)** only. **MSVC OpenMP (vcomp140.dll) is NOT used**.

**Rationale**:
1. **Intel MKL Compatibility**: Intel MKL uses Intel OpenMP internally; mixing with MSVC OpenMP causes conflicts
2. **Performance**: Intel OpenMP provides better threading performance on Intel CPUs
3. **Consistency**: Single OpenMP runtime avoids DLL conflicts and undefined behavior

**CMake Configuration**:
- Use standard `/openmp` compiler flag for MSVC
- Link against `libiomp5md.lib` directly instead of CMake's `OpenMP::OpenMP_CXX`
- CMake detects Intel OpenMP and links it instead of MSVC's vcomp140.dll

**Required DLLs** (auto-copied by BuildMSVC.ps1):
- `libiomp5md.dll` - Intel OpenMP runtime (from Intel oneAPI compiler directory)

**Verification**:
After building, verify only Intel OpenMP is loaded:
```python
import psutil, os
process = psutil.Process(os.getpid())
for dll in process.memory_maps():
    if 'omp' in dll.path.lower():
        print(dll.path)
# Should show: libiomp5md.dll
# Should NOT show: vcomp140.dll
```

**Note**: If both `libiomp5md.dll` and `vcomp140.dll` are loaded, there is a configuration error that will cause OpenMP parallelization to malfunction.

---

## Radia Solver Methods: MMM and MSC

Radia supports two solver methods:

### MMM (Magnetic Moment Method) - Tetrahedra
- Used for **tetrahedral elements** (4 faces)
- **3 DOF per element**: Magnetization vector (Mx, My, Mz)
- Represents magnetic objects as distributions of magnetic dipoles

### MSC (Magnetic Surface Charge) - Hexahedra
- Used for **hexahedral elements** (6 faces)
- **6 DOF per element**: Surface charge density (sigma) per face
- Computes field from surface charges using solid angle integration
- Use `ObjPolyhdr()` with `HEX_FACES` for hexahedral elements

**Note**: Radia does NOT use BEM (Boundary Element Method). The MSC method uses surface charges but differs from classical BEM.

### Mixed Hex/Tetra Element Support (2025-12-26)

Radia supports **mixed meshes** containing both hexahedral (6DOF) and tetrahedral (3DOF) elements.

**Solver Compatibility**:

| Solver | Mixed Elements | Notes |
|--------|----------------|-------|
| LU (Method 0) | **Supported** | Dense LU with variable DOF blocks |
| BiCGSTAB (Method 1) | **Supported** | Iterative solver with variable DOF |
| HACApK (Method 2) | NOT Supported | HACApK cluster tree requires uniform DOF |

**Note**: HACApK's cluster tree algorithm assumes fixed DOF per element. For mixed meshes, use LU (small problems) or BiCGSTAB (large problems).

**Interaction Matrix Blocks**:
- **3x3 block** (tetra-tetra): Standard demagnetization tensor
- **6x6 block** (hex-hex): Surface charge interaction (MSC)
- **3x6 block** (tetra from hex): H-field at tetra center from hex face charges
- **6x3 block** (hex from tetra): Normal dot N_matrix at hex eval points

**Usage**:
```python
import radia as rad

rad.FldUnits('m')

# Create mixed container with hex and tetra elements
HEX_FACES = [[1,4,3,2], [5,6,7,8], [1,2,6,5], [3,4,8,7], [1,5,8,4], [2,3,7,6]]
hex_vertices = [[0,0,0], [0.1,0,0], [0.1,0.1,0], [0,0.1,0],
                [0,0,0.1], [0.1,0,0.1], [0.1,0.1,0.1], [0,0.1,0.1]]
hex_obj = rad.ObjPolyhdr(hex_vertices, HEX_FACES, [0, 0, 0])   # 6DOF MSC
tetra_obj = rad.ObjPolyhdr(tetra_vertices, TETRA_FACES, [0, 0, 0])  # 3DOF MMM

container = rad.ObjCnt([hex_obj, tetra_obj])
mat = rad.MatSatIsoTab(BH_DATA)
rad.MatApl(container, mat)

# Solve with LU or BiCGSTAB (NOT HACApK)
rad.Solve(container, 0.001, 100, 0)  # Method 0 = LU
```

**Implementation Details**:
- `RADIA_MSC_SUPPORT` compile flag enables mixed element support
- Variable DOF offset arrays: `m_elemDOF`, `m_elemDOFOffset`, `m_totalDOF`
- Block interaction computation in `SetupInteractMatrix_VariableDOF()`

---

## Memory Management

### Exception Safety

All functions that allocate memory with `new` must follow this pattern:

```cpp
Type* ptr = nullptr;
try {
	ptr = new Type(...);
	Handle h(ptr);
	ptr = nullptr;  // Ownership transferred to handle
	...
}
catch(...) {
	if(ptr) delete ptr;  // Cleanup if exception before ownership transfer
	Initialize();
	return 0;
}
```

**Key Points**:
- Initialize raw pointers to `nullptr` before `try` block
- Set to `nullptr` immediately after ownership transfer
- Clean up in `catch(...)` block if pointer is still non-null

### RAII (Resource Acquisition Is Initialization)

Prefer RAII containers over manual memory management:

```cpp
// Good - RAII with std::vector
std::vector<radTPolygon> polygons;

// Avoid - Manual memory management
radTPolygon* polygons = new radTPolygon[n];  // Requires manual delete[]
```

---

## Unit System Policy

### Always Use Meters (SI Units)

**Policy**:
- **All examples** in `examples/` folder MUST use `rad.FldUnits('m')`
- **NGSolve integration** ALWAYS requires `rad.FldUnits('m')`

**Rationale**:
- Radia default: millimeters (mm)
- NGSolve default: meters (m)
- Without `rad.FldUnits('m')`, coordinates are off by 1000x

**Correct workflow**:
```python
import radia as rad
rad.FldUnits('m')  # REQUIRED for NGSolve integration

# Hexahedral magnet using ObjPolyhdr (8 vertices, magnetization in A/m)
HEX_FACES = [[1,4,3,2], [5,6,7,8], [1,2,6,5], [3,4,8,7], [1,5,8,4], [2,3,7,6]]
vertices = [[-0.02,-0.02,-0.03], [0.02,-0.02,-0.03], [0.02,0.02,-0.03], [-0.02,0.02,-0.03],
            [-0.02,-0.02,0.03], [0.02,-0.02,0.03], [0.02,0.02,0.03], [-0.02,0.02,0.03]]
magnet = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 954930])  # meters, A/m
```

---

## Radia Field Computation Limitations

### rad.Fld() Accuracy Inside Magnets

**Important Limitation**: `rad.Fld()` does **NOT** accurately compute field values **inside** permanent magnets.

**Rationale**:
- Radia MMM is designed for field calculation in **air regions** (outside magnetic materials)
- Inside magnets, `rad.Fld()` returns inaccurate values (known limitation, not a bug)

**Testing Strategy**:
- X **Avoid**: Direct comparison of `rad.Fld()` inside magnets
- OK **Use**: Large magnet with small mesh region (field approximately uniform)

### Vector Potential A Field Limitation (2025-12-27)

**Important Limitation**: `rad.Fld(obj, 'a', point)` returns `[0, 0, 0]` for `ObjPolyhdr` (MSC method).

**Current Status**:
- Vector potential A computation is **NOT implemented** for ObjPolyhdr/MSC elements
- B and H field computation works correctly for all element types
- A field returns zero for both hexahedral and tetrahedral ObjPolyhdr elements

**Affected Scripts**:
- `examples/ngsolve_integration/verify_curl_A_equals_B/` - Currently fails (curl(A) = 0)

**Future Implementation**:
The vector potential A can be computed from magnetization M using:
```
A(r) = (mu_0 / 4*pi) * integral( M(r') x (r - r') / |r - r'|^3 ) dV'
```
This requires volume integration over the magnetized element.

---

## Material Specification

### MatLin - Linear Materials

`rad.MatLin()` defines **linear magnetic materials** (soft magnetic materials, NOT permanent magnets).

**IMPORTANT**: MatLin is for **linear materials only**. For permanent magnets, use `ObjPolyhdr()` with magnetization vector.

### API Forms (Industry Standard)

```python
# Form 1: Isotropic linear material
mat = rad.MatLin(mu_r)  # Relative permeability (mu_r >= 1)

# Form 2: Anisotropic linear material with easy axis
mat = rad.MatLin([mu_r_par, mu_r_perp], [ex, ey, ez])
```

**Parameters**:
- **mu_r**: Relative permeability (industry standard, mu_r >= 1)
- **[mu_r_par, mu_r_perp]**: Parallel and perpendicular relative permeabilities
- **[ex, ey, ez]**: Easy axis direction vector (does NOT need normalization)

**Important Notes**:
1. **Linear materials ONLY**: MatLin is for soft magnetic materials (iron, steel, mu-metal, etc.)
2. **Permanent magnets**: Do NOT use MatLin - define magnetization directly in `ObjPolyhdr(vertices, HEX_FACES, [Mx,My,Mz])`
3. **Isotropic materials**: **ALWAYS prefer single-argument form `MatLin(mu_r)`** for isotropic materials.
4. **Easy axis**: For anisotropic materials, the easy axis vector must have significant magnitude (e.g., `[0, 0, 1]`)

**Example**:
```python
import radia as rad
rad.FldUnits('m')

# Define hexahedral element using ObjPolyhdr
HEX_FACES = [[1,4,3,2], [5,6,7,8], [1,2,6,5], [3,4,8,7], [1,5,8,4], [2,3,7,6]]

# Soft iron cube (isotropic, mu_r=4000)
iron_vertices = [[0,0,0], [0.1,0,0], [0.1,0.1,0], [0,0.1,0],
                 [0,0,0.1], [0.1,0,0.1], [0.1,0.1,0.1], [0,0.1,0.1]]
cube = rad.ObjPolyhdr(iron_vertices, HEX_FACES, [0, 0, 0])  # Zero magnetization
mat = rad.MatLin(4000)  # mu_r = 4000
rad.MatApl(cube, mat)

# Anisotropic material with easy axis in z-direction
iron_vertices2 = [[0.2,0,0], [0.3,0,0], [0.3,0.1,0], [0.2,0.1,0],
                  [0.2,0,0.1], [0.3,0,0.1], [0.3,0.1,0.1], [0.2,0.1,0.1]]
cube2 = rad.ObjPolyhdr(iron_vertices2, HEX_FACES, [0, 0, 0])
mat2 = rad.MatLin([5001, 101], [0, 0, 1])  # Easy axis along z
rad.MatApl(cube2, mat2)
```

### MatSatIsoTab - Nonlinear Materials (B-H Curve)

`rad.MatSatIsoTab()` defines **nonlinear isotropic magnetic materials** using a B-H curve.

**API (Industry Standard)**:

```python
# B-H curve: [[H1, B1], [H2, B2], ...] where H is in A/m and B is in Tesla
mat = rad.MatSatIsoTab(BH_data)
```

**Example**:
```python
# Steel B-H curve data
BH_DATA = [
    [0.0, 0.0],
    [100.0, 0.1],
    [500.0, 0.8],
    [1000.0, 1.2],
    [5000.0, 1.7],
    [50000.0, 2.0],
]

mat = rad.MatSatIsoTab(BH_DATA)
rad.MatApl(steel_obj, mat)
```

**Note**: Radia internally converts B-H to M-H using: M = B/mu_0 - H

### Permanent Magnet Materials

Radia provides several APIs for permanent magnet materials:

#### Method 1: ObjPolyhdr with Magnetization (Recommended for Fixed PM)

For permanent magnets where demagnetization is negligible, specify magnetization directly in `ObjPolyhdr`:

```python
import radia as rad
rad.FldUnits('m')

# Define hexahedral vertices (8 corners)
HEX_FACES = [[1,4,3,2], [5,6,7,8], [1,2,6,5], [3,4,8,7], [1,5,8,4], [2,3,7,6]]
vertices = [
    [-0.05, -0.05, -0.05],
    [0.05, -0.05, -0.05],
    [0.05, 0.05, -0.05],
    [-0.05, 0.05, -0.05],
    [-0.05, -0.05, 0.05],
    [0.05, -0.05, 0.05],
    [0.05, 0.05, 0.05],
    [-0.05, 0.05, 0.05],
]

# NdFeB magnet: Br = 1.2 T = 954930 A/m (= Br / mu_0)
Mr = 954930  # A/m
pm = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, Mr])

# Compute field (NO Solve needed for fixed PM)
B = rad.Fld(pm, 'b', [0, 0, 0.1])  # Field at z=0.1m
```

**Key Points**:
- Use `ObjPolyhdr` for arbitrary hexahedral shapes (8 vertices)
- Magnetization is specified in A/m: Mr = Br / mu_0
- **No `Solve()` required** for fixed magnetization permanent magnets
- Only call `Solve()` when soft iron is present in the model

#### Method 2: PM Material Classes (Future Demagnetization Support)

Three material classes are available for future demagnetization implementation:

```python
# MatMagFixed - Fixed magnetization (no demagnetization)
# Magnetization [Mx, My, Mz] in A/m
mat = rad.MatMagFixed([0, 0, 954930])

# MatMagLinear - Linear demagnetization (Br/Hc model)
# Br [T], Hc [A/m], easy axis [ex, ey, ez]
mat = rad.MatMagLinear(1.2, 955000, [0, 0, 1])

# MatMagCurve - User-defined B-H demagnetization curve
# [[H1,B1], [H2,B2], ...], easy axis [ex, ey, ez]
BH_curve = [[0.0, 1.2], [-500000, 0.6], [-955000, 0.0]]
mat = rad.MatMagCurve(BH_curve, [0, 0, 1])
```

**Note**: Currently all three material classes behave as fixed magnetization. Full demagnetization implementation is planned for future versions.

#### Permanent Magnet + Soft Iron Interaction

When combining permanent magnets with soft iron, use `Solve()` to calculate the induced magnetization in the iron:

```python
import radia as rad
rad.FldUnits('m')

HEX_FACES = [[1,4,3,2], [5,6,7,8], [1,2,6,5], [3,4,8,7], [1,5,8,4], [2,3,7,6]]

# PM magnet (with fixed magnetization)
pm_vertices = [...]  # 8 vertex coordinates
pm = rad.ObjPolyhdr(pm_vertices, HEX_FACES, [0, 0, 954930])

# Soft iron yoke (zero initial magnetization)
iron_vertices = [...]  # 8 vertex coordinates
iron = rad.ObjPolyhdr(iron_vertices, HEX_FACES, [0, 0, 0])
mat_iron = rad.MatLin(1000)  # mu_r = 1000
rad.MatApl(iron, mat_iron)

# Create assembly and solve
assembly = rad.ObjCnt([pm, iron])
result = rad.Solve(assembly, 0.0001, 1000, 0)  # LU solver

# Now compute fields
B = rad.Fld(assembly, 'b', [0, 0, 0.1])
```

---

## Windows Console Encoding (cp932) Compatibility

**Policy**: **NEVER use Unicode mathematical symbols** in print statements.

**Forbidden Unicode → ASCII Replacements**:

| Unicode | Symbol | ASCII | Example |
|---------|--------|-------|---------|
| `\u00b2` | ² | `^2` | `N²` → `N^2` |
| `\u00b3` | ³ | `^3` | `N³` → `N^3` |
| `\u2192` | → | `->` | `A → B` → `A -> B` |
| `\u2248` | ≈ | `~=` | `x ≈ 2` → `x ~= 2` |
| `\u2264` | ≤ | `<=` | `N ≤ 100` → `N <= 100` |
| `\u2265` | ≥ | `>=` | `N ≥ 250` → `N >= 250` |

**Rationale**: Windows console (cmd.exe) defaults to cp932 encoding in Japanese environments, causing `UnicodeEncodeError` for Unicode symbols.

---

## Python Script Path Import Policy

**Policy**: Use relative paths for module imports (not absolute paths).

```python
# ✓ CORRECT - Relative path
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
import radia as rad

# ✗ WRONG - Absolute path
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")
import radia as rad
```

**Path patterns**:
- Examples folder: `'../../build/Release'`
- Tests folder: `'../build/Release'`

---

## File Organization Policies

### Mesh File Preservation

**Policy**:
- **NEVER DELETE** mesh files (`.bdf`, `.nas`, `.msh`, `.vtk`)
- **NEVER DELETE** Cubit journal files (`.jou`, `.journal`)
- **NEVER DELETE** mesh generation scripts

**Rationale**: Mesh files are difficult to recreate without original CAD or mesh generation tools.

### Cubit Mesh Generation

**Policy**: Use `cubit_mesh_export` utilities for Cubit-based workflows.

**Requirements**:
- Journal files (`.jou`) MUST define blocks before export
- Use `cubit_mesh_export` for NASTRAN format conversion

**Correct workflow**:
```python
# In Cubit journal file (.jou):
# 1. Create geometry
# 2. Generate mesh
# 3. Define blocks (REQUIRED)
block 1 volume 1
block 1 element type hex8

# 4. Export using cubit_mesh_export
export nastran "geometry.bdf" dimension 3 overwrite

# In Python script:
from cubit_mesh_export import export_cubit_mesh
export_cubit_mesh('geometry.bdf', blocks={'1': {'material': 'NdFeB', 'Mr': [0, 0, 1.2]}})
```

**Common mistakes**:
```python
# WRONG - No block definition in .jou file
export nastran "geometry.bdf"  # Missing block assignment!

# WRONG - Block defined after export
export nastran "geometry.bdf"
block 1 volume 1  # Too late!
```

### VTK Export Policy

All example scripts should export VTK files with the same basename as the script.

---

## H-Matrix Acceleration Policy (HACApK ACA+)

### Policy: Do NOT Implement Custom H-Matrix - Use HACApK Only

**Date**: 2025-12-19
**Status**: Development Policy

**CRITICAL POLICY**:
1. **Do NOT implement custom H-matrix algorithms** (ACA, ACA+, or any low-rank approximation)
2. **If H-matrix acceleration is needed**, use the HACApK library at `src/ext/HACApK_LH-Cimplm/`
3. **Custom implementations have been removed** - rad_hmatrix*.cpp/h files were deleted (2025-12-18)

**Rationale**:
- Custom ACA implementations are prone to bugs and difficult to validate
- HACApK is a proven, MIT-licensed library with established correctness
- Benchmarks showed NO speedup for typical Radia use cases (single compact objects)

**Reference Implementation (for future use)**:
- **Library**: HACApK (Hierarchical Approximation with ACA+ for Krylov methods)
- **Source**: `src/ext/HACApK_LH-Cimplm/` (C implementation)
- **Original**: ppOpen-HPC project (MIT License)
- **Key Files**:
  - `cHACApK_base.h` - Core H-matrix structures and cluster tree
  - `cHACApK_base.c` - ACA+ implementation
  - `mpi_stub.h` - MPI stub for single-process execution (created for Radia)

**Why HACApK ACA+** (if needed in future):
1. **Proven Algorithm**: ACA+ is a well-established low-rank approximation method
2. **MIT License**: Compatible with Radia's license
3. **C Implementation**: Direct integration possible (no Fortran dependency)
4. **MPI Stub Available**: Can run without MPI dependency using provided stub

**Current Status (2025-12-23)**:

HACApK H-matrix solver is **implemented and available** as Method 2.

**When to Use HACApK**:
- Large-scale problems (N > 1,000 elements, >6,000 DOF)
- Memory-constrained environments (H-matrix uses O(N log N) memory vs O(N^2) for dense)
- NOT recommended for small problems where LU is faster

**Current Solver Methods** (v1.3.15):

| Method | Name | Description | Use Case |
|--------|------|-------------|----------|
| 0 | LU | Dense LU decomposition with LAPACK dgesv | Small problems (N < 500), guaranteed convergence |
| 1 | BiCGSTAB | Iterative BiCGSTAB with Jacobi preconditioner | General purpose, medium problems |
| 2 | HACApK | BiCGSTAB with H-matrix acceleration (ACA+) | Large problems (N > 1000), memory efficiency |

**Note**: Original Radia used Implicit SS (Gauss-Seidel) which had slow convergence for nonlinear materials. This was replaced with BiCGSTAB in v1.3.13. HACApK was added in v1.3.15.

### HACApK Parameters

Configure HACApK parameters using `rad.SetHACApKParams(eps, leaf_size, eta)` or `rad.SetHMatrixEpsilon(eps)`:

| Function | Description |
|----------|-------------|
| `SetHACApKParams(eps, leaf_size, eta)` | Set all H-matrix parameters |
| `SetHMatrixEpsilon(eps)` | Set only ACA tolerance (ELF-compatible: `magic.set_hmatrix_epsilon`) |
| `GetHACApKStats()` | Get H-matrix statistics after solve |

| Parameter | Description | Default | Range |
|-----------|-------------|---------|-------|
| **eps** | ACA tolerance (lower = more accurate, higher = faster) | 1e-4 | 1e-6 to 1e-2 |
| **leaf_size** | Minimum cluster size for leaf nodes | 10 | 5 to 50 |
| **eta** | Admissibility parameter (higher = more low-rank blocks) | 2.0 | 1.0 to 3.0 |

**Parameter Rationale**:

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `eps` | 1e-4 | Balance between accuracy and compression. ELF-compatible default. |
| `leaf_size` | 10 | Minimum cluster size. Smaller values allow deeper tree but increase overhead. 10 is ELF-compatible and provides good balance for typical problems. Values < 5 rarely improve performance; values > 50 reduce H-matrix benefit. |
| `eta` | 2.0 | Standard admissibility: clusters are "well-separated" when `dist(c1,c2) >= eta * max(diam(c1), diam(c2))`. eta=2.0 is conservative and ELF-compatible. Lower values (1.0-1.5) allow more aggressive compression but may reduce accuracy. |

**ACA Tolerance (eps) Guidelines**:

| eps | Accuracy | Speed | Use Case |
|-----|----------|-------|----------|
| 1e-6 | Highest | Slowest | Validation, reference solutions |
| 1e-4 | Good | Good | **Default, ELF-compatible** |
| 1e-3 | Moderate | Fast | Quick iterations, exploration |
| 1e-2 | Lower | Fastest | Initial guess, previews |

**Recommended Usage**:
```python
import radia as rad

# Method 0: Dense LU (direct solver, default)
rad.Solve(container, 0.0001, 1000, 0)

# Method 1: BiCGSTAB (iterative solver, better for large problems)
rad.Solve(container, 0.0001, 1000, 1)

# Method 2: HACApK (H-matrix accelerated BiCGSTAB)
rad.SetHACApKParams(1e-4, 10, 2.0)  # eps=1e-4 (ELF-compatible)
rad.Solve(container, 0.0001, 1000, 2)

# ELF-compatible: Set only epsilon (other params unchanged)
rad.SetHMatrixEpsilon(1e-3)  # magic.set_hmatrix_epsilon(1e-3)
rad.Solve(container, 0.0001, 1000, 2)

# Under-relaxation for difficult nonlinear problems
rad.SetRelaxParam(0.3)  # 30% damping (0.0 = full step, 0.0-1.0)
rad.Solve(container, 0.0001, 1000, 1)
rad.SetRelaxParam(0.0)  # Reset to full step
```

**Benchmark Script Usage**:
```bash
# Run HACApK benchmark with default eps=1e-4
python benchmark_hexahedron_msc.py --hacapk 5 10 15

# Run with custom ACA tolerance
python benchmark_hexahedron_msc.py --hacapk --eps 1e-3 5 10 15
```

**Documentation**:
- `docs/HMATRIX_EVALUATION.md` - Full evaluation report
- `examples/cube_uniform_field/nonlinear/` - Benchmark scripts and results

---

## NGSolve Integration Best Practices

### NGSolve Version Requirement

**CRITICAL**: Use NGSolve **6.2.2405** only.

**Version Constraint**: `ngsolve==6.2.2405`

**Reason**: NGSolve 6.2.2406+ has a regression bug in Periodic Boundary Conditions.
The `Identify()` information is lost during mesh generation with `Glue()`.

**Reference**: https://forum.ngsolve.org/t/ngsolve-periodic-boundary-condition-regression-bug-report/3805

**Installation**:
```bash
pip install radia[ngsolve]  # Installs with correct NGSolve version constraint
# OR
pip install ngsolve==6.2.2405
```

### Recommended Configuration

```python
fes = HDiv(mesh, order=2)  # Best accuracy
B_gf = GridFunction(fes)
B_gf.Set(radia_ngsolve.RadiaField(radia_obj, 'b'))
```

**Evaluation guidelines**:
- Evaluate GridFunction at distances > 1 mesh cell from magnet surface
- Use CoefficientFunction directly for maximum accuracy near boundaries
- Avoid GridFunction evaluation within 1 mesh cell of magnet surface

---

## PyPI Package Release Policy

### Version Management (Automated by Claude Code)

Claude Code is responsible for:
- Maintaining version numbers in `pyproject.toml` and `src/radia/__init__.py`
- Following semantic versioning (MAJOR.MINOR.PATCH)
- Updating `CHANGELOG.md` with release notes
- Committing and pushing version bump changes

### PyPI Upload (Manual by User)

**IMPORTANT**: Claude Code does NOT execute PyPI upload. The user must manually run the upload script.

**Workflow**:
1. Claude Code prepares the release (version bump, CHANGELOG, commit, push)
2. Claude Code runs `python -m build` to verify package builds correctly
3. **User manually executes** `Publish_to_PyPI.ps1` with their API token

```powershell
# Set PyPI API token (keep secure!)
$env:PYPI_TOKEN = "pypi-AgEIcGl..."

# Run upload script
powershell.exe -ExecutionPolicy Bypass -File Publish_to_PyPI.ps1
```

**Security**: NEVER commit PyPI tokens to repository. NEVER ask user for their token.

### MKL DLL Bundling Policy (2025-12-27)

**Policy**: PyPI packages MUST include Intel MKL runtime DLLs for user convenience.

**Required DLLs** (Windows):
- `mkl_core.2.dll` - MKL core library
- `mkl_intel_thread.2.dll` - MKL Intel OpenMP threading
- `mkl_def.2.dll` - MKL definitions
- `mkl_rt.2.dll` - MKL runtime (optional, for dynamic linking)
- `libiomp5md.dll` - Intel OpenMP runtime

**Source Location** (Intel oneAPI):
```
C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin\
C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\
```

**Package Structure**:
```
src/radia/
  __init__.py
  radia.pyd           # C++ extension
  mkl_core.2.dll      # MKL runtime (bundled)
  mkl_intel_thread.2.dll
  mkl_def.2.dll
  libiomp5md.dll      # Intel OpenMP
```

**Rationale**:
- Users should NOT need to install Intel oneAPI to use Radia
- Bundling DLLs ensures consistent behavior across environments
- Reduces support burden from "DLL not found" errors

**Build Integration**:
- `BuildWithIntel.ps1` should copy required DLLs to `src/radia/`
- `setup.py` / `pyproject.toml` should include DLLs in `package_data`

**License Consideration**:
Intel MKL runtime DLLs are redistributable under Intel oneAPI EULA.
Include appropriate license notices in the package.

---

## Nastran Mesh Import Unification (2025-11-23)

### Migration: nastran_reader.py → nastran_mesh_import.py

**Date**: 2025-11-23
**Status**: Complete

### Changes

**Removed**:
- `src/python/nastran_reader.py` - Legacy Nastran reader (deprecated)

**Enhanced**:
- `src/python/nastran_mesh_import.py` - Unified Nastran import module

### Supported Element Types

`nastran_mesh_import.py` now supports all major 3D element types:

| Element Type | Nastran Card | Nodes | Status |
|--------------|--------------|-------|--------|
| Hexahedron | CHEXA | 8 | ✓ Supported |
| Wedge/Prism | CPENTA | 6 | ✓ Supported |
| Pyramid | CPYRAM | 5 | ✓ Supported |
| Tetrahedron | CTETRA | 4 | ✓ Supported |
| Triangle (Surface) | CTRIA3 | 3 | ✓ Supported |

### CTRIA3 Surface Mesh Support

**Key Feature**: CTRIA3 elements are grouped by material ID (property ID).

- Each material ID creates **one polyhedron** from all its triangles
- Enables surface-based magnetic analysis
- Compatible with sphere.bdf (8 material groups, 7408 total faces)

**Usage**:
```python
from nastran_mesh_import import import_nastran_mesh, create_radia_from_nastran

# Read mesh
mesh_data = import_nastran_mesh('sphere.bdf', units='mm')

# Access triangle groups
tria_groups = mesh_data['tria_groups']
# Format: {material_id: {'faces': [[n1,n2,n3], ...], 'node_ids': set(...)}}

# Create Radia objects automatically
mag_obj = create_radia_from_nastran('sphere.bdf',
                                     material={'magnetization': [0, 0, 1.2]},
                                     units='mm')
```

### Migration Guide

**Before** (using nastran_reader.py):
```python
from nastran_reader import read_nastran_mesh, TETRA_FACES

mesh = read_nastran_mesh(nas_file)
nodes = mesh['nodes']  # numpy array
tetra_elements = mesh['tetra_elements']  # list
tria_groups = mesh['tria_groups']  # dict
```

**After** (using nastran_mesh_import.py):
```python
from nastran_mesh_import import import_nastran_mesh, create_radia_from_nastran
from netgen_mesh_import import TETRA_FACES, WEDGE_FACES, PYRAMID_FACES

# Option 1: Parse only
mesh = import_nastran_mesh(nas_file, units='mm')
vertices = mesh['vertices']  # list of [x,y,z]
tet_elements = mesh['tet_elements']  # list of vertex indices
tria_groups = mesh['tria_groups']  # dict (same format)

# Option 2: Create Radia objects directly (recommended)
mag_obj = create_radia_from_nastran(nas_file,
                                     material={'magnetization': [0, 0, 1.2]},
                                     units='mm')
```

### Affected Files

**Deprecated**:
- `examples/background_fields/sphere_nastran_analysis.py` - Marked as DEPRECATED, kept for reference

**Note**: If issues arise with Nastran import, refer to `nastran_mesh_import.py` as the single source of truth.

---

## Unit System Policy: No Hard-Coded Unit Conversions

### Requirement: Centralized Unit Control via rad.FldUnits() and radia_ngsolve

**Goal**: All unit conversions must be controlled through explicit API calls (`rad.FldUnits()` or `radia_ngsolve` constructor), never through hard-coded conversion factors in user code.

**Policy**:

**✓ ALLOWED - Explicit unit control**:
```python
# Method 1: Set Radia units globally
import radia as rad
rad.FldUnits('m')  # All Radia operations now use meters

# Create hexahedral magnet with ObjPolyhdr
HEX_FACES = [[1,4,3,2], [5,6,7,8], [1,2,6,5], [3,4,8,7], [1,5,8,4], [2,3,7,6]]
vertices = [[-0.05,-0.05,-0.05], [0.05,-0.05,-0.05], [0.05,0.05,-0.05], [-0.05,0.05,-0.05],
            [-0.05,-0.05,0.05], [0.05,-0.05,0.05], [0.05,0.05,0.05], [-0.05,0.05,0.05]]
magnet = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 954930])  # 0.1m cube, 1.2T equivalent

# Method 2: Specify units in radia_ngsolve constructor
from radia_ngsolve import RadiaField
B_cf = RadiaField(magnet, 'b', units='m')  # Explicitly use meters
```

**✗ FORBIDDEN - Hard-coded unit conversions**:
```python
# WRONG - Hard-coded mm to m conversion
for pt in obs_points:
    f.write(f'{pt[0]/1000.0} {pt[1]/1000.0} {pt[2]/1000.0}\n')  # ✗ DO NOT DO THIS

# WRONG - Hard-coded scaling factors
x_mm = x_m * 1000.0  # ✗ DO NOT DO THIS
field_m = field_mm / 1000.0  # ✗ DO NOT DO THIS
```

**Rationale**:
- **Single source of truth**: Units controlled by `rad.FldUnits()` only
- **Consistency**: All code uses same unit system set at initialization
- **Maintainability**: Changing units requires one line change, not searching for conversion factors
- **Error prevention**: Hard-coded conversions cause bugs when unit system changes

**Unit Detection**:

Use `rad.FldUnits()` without arguments to get current unit system:

```python
import radia as rad

# Get current units (returns multi-line string)
units_str = rad.FldUnits()
# Parse to detect length unit
if 'Length:  mm' in units_str:
    length_unit = 'mm'
    length_scale = 1.0  # No conversion needed
elif 'Length:  m' in units_str:
    length_unit = 'm'
    length_scale = 0.001  # mm to m
else:
    raise ValueError(f"Unknown length unit in: {units_str}")
```

**No Exceptions**:

All code, including `radia_vtk_export.py`, `radia_ngsolve.cpp`, `nastran_mesh_import.py`, must:
1. Query current unit system via `rad.FldUnits()`
2. Convert based on detected units, not hard-coded assumptions
3. Never assume Radia is using mm or m

This ensures code works regardless of user's `rad.FldUnits()` setting.

**Implementation Pattern**:

```python
# CORRECT - Use rad.FldUnits() to control units
import radia as rad

# Set unit system once at start
rad.FldUnits('mm')  # or 'm' for NGSolve integration

# Create hexahedral magnet with ObjPolyhdr (100mm cube)
HEX_FACES = [[1,4,3,2], [5,6,7,8], [1,2,6,5], [3,4,8,7], [1,5,8,4], [2,3,7,6]]
vertices = [[-50,-50,-50], [50,-50,-50], [50,50,-50], [-50,50,-50],
            [-50,-50,50], [50,-50,50], [50,50,50], [-50,50,50]]
magnet = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 954930])  # 100mm, A/m
field = rad.Fld(magnet, 'b', [50, 50, 50])  # 50mm point

# Export - automatically uses correct units
from radia_vtk_export import exportGeometryToVTK
exportGeometryToVTK(magnet, 'output.vtk')  # Handles units internally
```

**Migration Guidelines**:

When removing hard-coded unit conversions:

1. **Identify conversion factors**: Search for `/1000`, `*1000`, `0.001`, etc.
2. **Determine intended unit system**: mm or m?
3. **Add `rad.FldUnits()` at script start**: Set unit system explicitly
4. **Remove conversion factors**: Use values directly in chosen unit system
5. **Update comments**: Document which unit system is used

**Files to Check**:

When writing or modifying code:
- ✓ Check for hard-coded `*1000`, `/1000`, `*0.001`, `/0.001`
- ✓ Ensure `rad.FldUnits()` is called at script start
- ✓ Verify no manual coordinate scaling
- ✓ Use `radia_vtk_export.py` for VTK export (handles units)

**Example - Corrected Code**:

Before (hard-coded conversions):
```python
# BAD - Hard-coded unit conversion
x_range = np.linspace(-90, 90, 21)  # mm
for pt in obs_points:
    f.write(f'{pt[0]/1000.0} {pt[1]/1000.0} {pt[2]/1000.0}\n')  # Manual mm->m
```

After (unit-aware):
```python
# GOOD - Use rad.FldUnits() and exportGeometryToVTK
rad.FldUnits('m')  # Set to meters
x_range = np.linspace(-0.09, 0.09, 21)  # m (no conversion needed)

# Use radia_vtk_export for geometry
from radia_vtk_export import exportGeometryToVTK
exportGeometryToVTK(magnet, 'output.vtk')  # Automatic unit handling

# For field data, use same unit system
for pt in obs_points:
    f.write(f'{pt[0]} {pt[1]} {pt[2]}\n')  # Already in meters
```

---

## NGSolve Mesh Access Policy

### Centralized Mesh Access via netgen_mesh_import.py

**Policy**:
- **All NGSolve mesh access** MUST use functions from `src/radia/netgen_mesh_import.py`
- **NEVER** directly access `mesh.ngmesh.Points()`, `mesh.vertices[]`, or `el.vertices[].nr` in any script
- **ALWAYS** import mesh handling functions from `netgen_mesh_import.py`
- **NO EXCEPTIONS**: This applies to all scripts including examples, tests, and debugging code

**Enforcement**:
- Direct mesh access is a bug source due to index confusion
- All new code MUST use `extract_elements()` or `netgen_mesh_to_radia()`
- Existing code with direct access MUST be refactored

**Rationale**:

NGSolve has two different indexing schemes that cause off-by-one errors:

| Access Method | Indexing | Notes |
|--------------|----------|-------|
| `mesh.ngmesh.Points()[i]` | **1-indexed** | Index 0 raises error, valid: 1 to nv |
| `mesh.vertices[i]` | **0-indexed** | Valid: 0 to nv-1 |
| `el.vertices[i].nr` | Returns value for **0-indexed** `mesh.vertices[]` | Use with `mesh.vertices[]` only |

**Common Bug Pattern**:
```python
# WRONG - Using 0-indexed .nr with 1-indexed ngmesh.Points()
for v in el.vertices:
    pt = mesh.ngmesh.Points()[v.nr]  # Off-by-one error!

# CORRECT - Use 0-indexed consistently
for v in el.vertices:
    vertex = mesh.vertices[v.nr]
    pt = vertex.point
```

**Correct Usage**:

```python
# Import from centralized module
from netgen_mesh_import import netgen_mesh_to_radia, extract_elements, TETRA_FACES

# Option 1: Direct conversion to Radia (recommended)
radia_obj = netgen_mesh_to_radia(mesh,
                                  material={'magnetization': [0, 0, 0]},
                                  units='m',
                                  material_filter='magnetic')

# Option 2: Extract elements for custom processing
elements, _ = extract_elements(mesh, material_filter='magnetic')
for el in elements:
    vertices = el['vertices']  # Already extracted correctly
    # ...
```

**Module Location**: `src/radia/netgen_mesh_import.py`

**Available Functions**:
- `netgen_mesh_to_radia()`: Convert entire mesh to Radia geometry (recommended)
- `extract_elements()`: Extract element data for custom processing
- `compute_element_centroid()`: Compute centroid from vertex list
- `create_radia_tetrahedron()`: Create single Radia tetrahedron
- `create_radia_hexahedron()`: Create single Radia hexahedron

**Available Constants**:
- `TETRA_FACES`: 1-indexed face topology for tetrahedra
- `HEX_FACES`: 1-indexed face topology for hexahedra
- `WEDGE_FACES`: 1-indexed face topology for wedges
- `PYRAMID_FACES`: 1-indexed face topology for pyramids

---

## Example Script Naming Convention

### Requirement: Consistent snake_case Naming with Functional Prefixes

**Goal**: All example scripts in `examples/` folder must follow a consistent naming convention for easy identification and organization.

**Policy**:

**Naming pattern**: `<prefix>_<description>.py`

- Use **snake_case** (all lowercase with underscores)
- Use **functional prefix** to indicate script purpose
- Use **descriptive names** that explain what the script does

**Standard prefixes**:

| Prefix | Purpose | Example |
|--------|---------|---------|
| `demo_` | Educational demonstration of a feature | `demo_batch_evaluation.py` |
| `example_` | Complete working example | `example_hmatrix_cache_usage.py` |
| `benchmark_` | Performance measurement script | `benchmark_solver_scaling.py` |
| `test_` | Validation/verification script | `test_batch_evaluation.py` |
| `verify_` | Correctness verification | `verify_curl_A_equals_B.py` |
| `compare_` | Comparison between methods | `compare_radia_ngsolve_cube.py` |
| `visualize_` | Visualization script | `visualize_field.py` |
| `run_` | Runner/orchestrator script | `run_all_benchmarks.py` |
| (none) | Descriptive physical model name | `sphere_in_quadrupole.py`, `arc_current_with_magnet.py` |

**Naming rules**:

1. **✓ CORRECT - snake_case**:
   ```
   sphere_in_quadrupole.py
   benchmark_solver_scaling.py
   demo_batch_evaluation.py
   verify_curl_A_equals_B.py
   ```

2. **✗ INCORRECT - CamelCase or PascalCase**:
   ```
   Cubit2Nastran.py        # Should be: cubit_to_nastran.py
   York_cubit_mesh.py      # Should be: york_cubit_mesh.py (already correct case, but York should be lowercase)
   CompareResults.py       # Should be: compare_results.py
   ```

3. **✗ INCORRECT - No prefix for functional scripts**:
   ```
   accuracy.py             # Should be: verify_accuracy.py or benchmark_accuracy.py
   plot.py                 # Should be: visualize_results.py or plot_benchmark_results.py
   ```

**Directory-specific guidelines**:

| Directory | Typical Prefixes | Notes |
|-----------|------------------|-------|
| `examples/simple_problems/` | (none), `demo_` | Physical model names preferred |
| `examples/solver_benchmarks/` | `benchmark_`, `run_`, `plot_`, `verify_` | Performance focus |
| `examples/ngsolve_integration/` | `demo_`, `example_`, `test_`, `verify_` | Educational + validation |
| `examples/background_fields/` | (none), `compare_` | Physical model names |
| `examples/electromagnet/` | `main_`, `visualize_` | Workflow scripts |

**Migration checklist**:

When renaming files:
1. Use `git mv` to preserve history: `git mv OldName.py new_name.py`
2. Update imports in other files
3. Update README.md references
4. Update documentation
5. Commit with clear message: `"Rename OldName.py to new_name.py (naming convention)"`

**Examples of good names**:

```
# Physical models - descriptive, no prefix needed
sphere_in_quadrupole.py           # Clear physics description
arc_current_with_magnet.py        # Clear what it models
cubic_polyhedron_magnet.py        # Clear geometry + physics

# Functional scripts - prefix required
demo_batch_evaluation.py          # Demo of batch feature
benchmark_solver_scaling.py       # Benchmark solver performance
verify_curl_A_equals_B.py         # Verify Maxwell equation
compare_radia_ngsolve.py          # Compare two methods
visualize_field.py                # Visualize field data
run_all_benchmarks.py             # Orchestrator script
```

**Files to rename** (current violations):

1. `background_fields/Cubit2Nastran.py` → `background_fields/cubit_to_nastran.py`
2. `electromagnet/York_cubit_mesh.py` → `electromagnet/york_cubit_mesh.py`

**Rationale**:
- **Consistency**: Easy to scan and find scripts by purpose
- **Clarity**: Prefix immediately indicates script type
- **Python convention**: PEP 8 recommends snake_case for module names
- **Sorting**: Related scripts group together alphabetically

---

## Mesh Generation Tools

### Tool Selection by Element Type

| Element Type | Tool | Notes |
|--------------|------|-------|
| **Tetrahedral** | **Netgen** | Recommended. Uses `netgen.occ.Box` + `OCCGeometry.GenerateMesh()` |
| Tetrahedral | GMSH | Alternative. Export as .msh and import with nastran_mesh_import.py |
| Tetrahedral | Nastran | CTETRA elements from .bdf files |
| **Hexahedral** | **Cubit** | Recommended for complex hex mesh generation |
| Hexahedral | Nastran | CHEXA elements from .bdf files |

### Netgen Tetrahedral Mesh Generation

```python
from netgen.occ import Box, Pnt, OCCGeometry
from netgen.meshing import MeshingParameters

# Create box geometry
p1 = Pnt(-0.5, -0.5, -0.5)
p2 = Pnt(0.5, 0.5, 0.5)
box = Box(p1, p2)

# Generate mesh - MUST wrap in OCCGeometry first
geo = OCCGeometry(box)
mp = MeshingParameters(maxh=0.3)  # Maximum element size
mesh = geo.GenerateMesh(mp)

# Extract nodes and tetrahedra
nodes = [[p[0], p[1], p[2]] for p in mesh.Points()]
tetrahedra = [[v.nr - 1 for v in el.vertices] for el in mesh.Elements3D()]
```

**Important Notes**:
- Use `OCCGeometry(box).GenerateMesh(mp)`, NOT `box.GenerateMesh(mp)`
- Use `v.nr - 1` for vertex indices, NOT `int(v) - 1`

### Cubit Hexahedral Mesh Notes

- Use Cubit for complex hexahedral meshing
- Export as NASTRAN format (.bdf)
- Import with `nastran_mesh_import.import_nastran_mesh()`
- Ensure blocks are defined before export: `block 1 volume 1`

---

## MSC (Magnetic Surface Charge) Method

### Overview

Radia uses **MSC (Magnetic Surface Charge)** for all hexahedral elements:

| Element Type | Faces | DOF | API | Use Case |
|--------------|-------|-----|-----|----------|
| **Tetrahedron** | 4 triangular | 3 (Mx, My, Mz) | `ObjPolyhdr()` + `TETRA_FACES` | Complex curved geometry |
| **Hexahedron** | 6 quadrilateral | 6 (sigma per face) | `ObjPolyhdr()` + `HEX_FACES` | Permanent magnets, soft iron |

**Policy (2025-12-27)**:
- **ObjPolyhdr()** with HEX_FACES: Standard API for hexahedral elements (6 DOF MSC)
- **Tetrahedron**: 3 DOF (Mx, My, Mz) - MMM method with uniform magnetization
- **Hexahedron**: 6 DOF (sigma per face) - MSC method with surface charges
- 3 DOF hexahedron (MMM) is NOT supported - all hexahedra use 6 DOF MSC
- All meshes are expected to be generated externally (Netgen, GMSH, Cubit, etc.)

### Implementation

**Source files** (`src/core/`):
- `rad_polyhedron.cpp`: Element type detection and dispatch
- `B_comp_tetrahedron_MSC()`: 4-face tetrahedral field computation
- `B_comp_hexahedron_MSC()`: 6-face hexahedral field computation
- `FieldFromChargedTriangle()`: Analytic field from charged triangle (solid angle formula)
- `FieldFromQuadFace()`: Quad face split into 2 triangles
- `rad_poly_analytical.cpp`: `RadFieldFromTriangleFaceGlobal()` for triangle integration

**Key Features**:
- Uses **global coordinates** directly (no local coordinate transformations)
- Computes field using **solid angle integration** formula (van Oosterom & Strackee, 1983)
- Handles **outward normal orientation** automatically

### Quick Start

**Tetrahedral mesh (Netgen)**:

```python
import radia as rad
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

rad.FldUnits('m')

# Generate tetrahedral mesh
cube_solid = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
cube_solid.mat('magnetic')
geo = OCCGeometry(cube_solid)
mesh = Mesh(geo.GenerateMesh(maxh=0.3))

# Import to Radia (uses MSC automatically)
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='magnetic')

# Apply material and solve
mat = rad.MatLin(999)  # mu_r = 1000
rad.MatApl(mag_obj, mat)
rad.Solve(mag_obj, 0.0001, 1000, 1)  # Method 1 = BiCGSTAB
```

**Hexahedral element (ObjPolyhdr)**:

```python
import radia as rad
from netgen_mesh_import import HEX_FACES

rad.FldUnits('m')

# Create hexahedral element with explicit vertices
vertices = [[0,0,0], [1,0,0], [1,1,0], [0,1,0],
            [0,0,1], [1,0,1], [1,1,1], [0,1,1]]
hex_obj = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])

# Apply material and solve
mat = rad.MatLin(999)  # mu_r = 1000
rad.MatApl(hex_obj, mat)
rad.Solve(hex_obj, 0.0001, 1000, 1)
```

### Benchmark Results (2025-12-19)

**Nonlinear Material (H_ext = 50,000 A/m)**:

| Element Type | N_elem | Solver | Iterations | M_avg_z (A/m) |
|--------------|--------|--------|------------|---------------|
| Hexahedral MSC | 125 | LU | 3 | 173,400 |
| Tetrahedral MSC | ~200 | LU | 2 | ~190,000 |

**Notes**:
- Tetrahedron: 3 DOF (Mx, My, Mz) - uniform magnetization
- Hexahedron: 6 DOF (sigma per face) - surface charge density
- LU solver (Method 0) and BiCGSTAB (Method 1) both work

### Documentation

- [docs/MSC_QUICK_START.md](docs/MSC_QUICK_START.md): Quick start guide
- [examples/cube_uniform_field/nonlinear/](examples/cube_uniform_field/nonlinear/): Nonlinear benchmark

---

## PyPI Package Structure (v1.3.8+)

### Package Directory: src/radia (NOT src/python)

**Critical Requirement**: The package directory MUST be named `src/radia` for `pip install radia` to work correctly.

**Rationale**:
- `setuptools.find_packages(where="src")` discovers packages by directory name
- If the directory is named `python`, the wheel creates `python/` folder, not `radia/`
- Users would need `import python` instead of `import radia` - completely wrong!

**Correct Structure**:
```
src/
  radia/              # Package directory (import radia works)
    __init__.py       # Re-exports symbols from radia.pyd
    radia.pyd         # Core C++ extension module
    radia_ngsolve.pyd # Optional NGSolve integration
    *.py              # Python utility modules
```

**__init__.py Requirements**:
```python
# Radia Python package
__version__ = "1.3.8"

# Import all symbols from the C++ extension module
try:
    from radia.radia import *
except ImportError:
    from .radia import *
```

**Key Points**:
1. **Directory name = Package name**: `src/radia/` -> `import radia`
2. **C++ module import**: `__init__.py` must re-export symbols from `radia.pyd`
3. **setup.py package_data**: Use `"radia"` not `"python"` in package_data dict
4. **Version sync**: Keep version consistent across `__init__.py`, `setup.py`, `pyproject.toml`

### Build Policy: Always Use BuildRadiaInternal.ps1

**CRITICAL**: Always use the standard build script to ensure .pyd files are correctly copied.

**Policy**:
1. **ALWAYS use `BuildRadiaInternal.ps1`** for building - NEVER use manual cmake commands
2. The script automatically copies .pyd to `src/radia/radia.pyd` after build
3. The script verifies the copy was successful

**Standard Build Command**:
```powershell
# From project root directory:
powershell.exe -ExecutionPolicy Bypass -File BuildRadiaInternal.ps1

# With options:
powershell.exe -ExecutionPolicy Bypass -File BuildRadiaInternal.ps1 -Rebuild  # Clean rebuild
powershell.exe -ExecutionPolicy Bypass -File BuildRadiaInternal.ps1 -Verbose  # Verbose output
```

**Why This Matters**:
- CMake outputs `radia.cp312-win_amd64.pyd` (version-tagged name)
- Python package expects `src/radia/radia.pyd` (simple name)
- Manual cmake builds do NOT copy/rename the .pyd file
- Using old .pyd causes "AttributeError: module has no attribute" errors

**What the Script Does**:
1. Runs CMake configure and build
2. Finds `radia.cp312-win_amd64.pyd` in `build/Release/`
3. Copies it to `src/radia/radia.pyd` (with rename)
4. Verifies the copy with timestamp check

**NEVER Do This**:
```powershell
# WRONG - Does NOT copy .pyd to src/radia
cmake --build build --config Release --target radia

# WRONG - Manual copy forgets the rename
copy build\Release\radia.cp312-win_amd64.pyd src\radia\
```

**If You Must Use Manual Commands** (not recommended):
```powershell
# Build
cmake --build "s:\Radia\01_GitHub\build" --config Release --target radia

# MUST copy with rename
Copy-Item "build\Release\radia.cp312-win_amd64.pyd" "src\radia\radia.pyd" -Force

# Verify
(Get-Item "src\radia\radia.pyd").LastWriteTime
```

**Rationale**:
- `src/radia/` is the package source directory
- `python -m build` uses files from `src/radia/` to create wheels
- Ensures wheel contains the latest built binaries
- Avoids confusion about which .pyd is the "current" version

**Note**: The `.gitignore` excludes `.pyd` files, so these copied files are NOT committed to git.

**Testing Checklist** (before PyPI upload):
1. Build wheel: `python -m build`
2. Inspect wheel: `unzip -l dist/radia-x.y.z-py3-none-any.whl`
3. Verify structure: Contents should show `radia/radia.pyd`, NOT `python/radia.pyd`
4. Test install: `pip install dist/radia-x.y.z-py3-none-any.whl`
5. Test import: `python -c "import radia; print(radia.__version__)"`

---

## Benchmark Policy

### ベンチマーク実行ルール

1. **1ケース毎に実行**: 複数のソルバーやパラメータを比較する場合でも、1ケースずつ実行して結果を確認
2. **シーケンシャル実行のみ**: ベンチマークは並列実行せず、必ず1つずつ順番に実行する（メモリ測定の正確性のため）
3. **メモリ使用量を記録**: ピークメモリ使用量を測定し、結果JSONに含める

### メモリ測定方法

**psutil を使用**（C++拡張のメモリも含めて測定可能）:

```python
import psutil
import os

def get_peak_memory_mb():
    """Get peak memory usage in MB (Windows: peak_wset)"""
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    # Windows: peak_wset = peak working set size
    # Linux/Mac: rss = resident set size (ピークは別途追跡必要)
    if hasattr(mem_info, 'peak_wset'):
        return mem_info.peak_wset / (1024 * 1024)  # MB
    else:
        return mem_info.rss / (1024 * 1024)  # MB (fallback)
```

**注意**: `tracemalloc` はPython内部のメモリのみ追跡するため、C++拡張モジュール (radia.pyd) のメモリは測定されない。

### 結果JSONフォーマット

```json
{
  "n_elements": 1000,
  "solver_method": 0,
  "t_solve": 5.2,
  "peak_memory_mb": 450.0,
  "iterations": 6,
  "M_avg_z": 149846.0
}
```

---

## Naming Policy: External Project References

### Do NOT Use "ELF" or "ELF_MAGIC" in Radia Codebase

**Policy (2025-12-16)**:
- **Do NOT use "ELF" or "ELF_MAGIC"** as terminology in Radia source code, documentation, or comments
- The MSC (Magnetic Surface Charge) method in Radia is a **Radia-native implementation**
- References to external projects should be kept to academic citations only

**Rationale**:
- Radia is an independent project with its own implementation
- Avoid confusion with external codebases
- Maintain clear intellectual property boundaries

**Allowed**:
- Academic citations in documentation (e.g., "Reference: Yano & Sugahara, J. Magn. Soc. Jpn., 2023")
- General algorithm descriptions (e.g., "MSC method", "surface charge method")

**Not Allowed**:
- Variable/function names containing "ELF" (e.g., `ELF_MAGIC_compatible`)
- Comments like "following ELF_MAGIC convention"
- Documentation referring to "ELF_MAGIC format"

**Migration**: If existing code contains "ELF" references, replace with Radia-native terminology.

---

## Logging and Debugging Policy

### No Console Output from C++ Code

**Policy (2025-12-21)**:
- **Do NOT add printf/cout/cerr output in C++ code** for logging or debugging
- All user-facing output should be handled through **Python scripts**
- C++ code should be silent unless there is a critical error that requires immediate attention

**Rationale**:
- Python stdout/stderr can be captured and formatted properly
- C++ printf may not appear in Python console depending on buffering
- Keeps separation of concerns: C++ for computation, Python for user interaction

**Allowed**:
- Error messages via Radia's error handling system (`Send.ErrorMessage(...)`)
- Debug output controlled by compile-time flags (`#ifdef DEBUG_...`)

**Not Allowed**:
- Unconditional `printf()`, `fprintf(stderr, ...)`, `std::cout`, `std::cerr`
- Timing output from C++ code (use Python's `time.time()` instead)

**Debugging Approach**:
1. Use Python-side timing: `time.time()` before/after calls
2. Use conditional compile flags: `#ifdef RADIA_DEBUG_TIMING`
3. Write to file only when explicitly enabled by environment variable

---

**Last Updated**: 2025-12-21 (Added logging policy, DOF documentation: Tetra=3, Hexa=6)
**For**: Claude Code AI Assistant
**Project**: Radia Magnetic Field Computation