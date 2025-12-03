# Claude Code - Radia Project Development Guidelines

This document contains development guidelines and refactoring policies for the Radia project when working with Claude Code.

## Radia Solver Method: MMM (Magnetic Moment Method)

**IMPORTANT**: Radia uses the **MMM (Magnetic Moment Method)**, NOT BEM (Boundary Element Method).

**Terminology**:
- ✓ Correct: "Radia MMM solver", "MMM computation", "MMM field evaluation"
- ✗ Incorrect: "Radia BEM", "BEM solver", "boundary element method"

**Context**:
- MMM represents magnetic objects as distributions of magnetic moments
- This differs from BEM which uses surface integral equations
- All references to "BEM" in Radia-related documentation should be corrected to "MMM"

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
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.06], [0, 0, 1.2])  # meters
```

---

## Radia Field Computation Limitations

### rad.Fld() Accuracy Inside Magnets

**Important Limitation**: `rad.Fld()` does **NOT** accurately compute field values **inside** permanent magnets.

**Rationale**:
- Radia MMM is designed for field calculation in **air regions** (outside magnetic materials)
- Inside magnets, `rad.Fld()` returns inaccurate values (known limitation, not a bug)

**Testing Strategy**:
- ✗ **Avoid**: Direct comparison of `rad.Fld()` inside magnets
- ✓ **Use**: Large magnet with small mesh region (field approximately uniform)

---

## Material Specification (MatLin)

`rad.MatLin()` defines **linear magnetic materials** (soft magnetic materials, NOT permanent magnets).

**IMPORTANT**: MatLin is for **linear materials only** (materials with magnetic susceptibility). For permanent magnets, use `rad.ObjRecMag()` with magnetization vector.

### API Forms

```python
# Form 1: Isotropic linear material
mat = rad.MatLin(ksi)  # Single susceptibility value

# Form 2: Anisotropic linear material with easy axis
mat = rad.MatLin([ksi_par, ksi_perp], [ex, ey, ez])
```

**Parameters**:
- **ksi**: Isotropic magnetic susceptibility (χ = μr - 1)
- **[ksi_par, ksi_perp]**: Parallel and perpendicular susceptibilities
- **[ex, ey, ez]**: Easy axis direction vector (does NOT need normalization)

**Important Notes**:
1. **Linear materials ONLY**: MatLin is for soft magnetic materials (iron, steel, mu-metal, etc.)
2. **Permanent magnets**: Do NOT use MatLin - define magnetization directly in `rad.ObjRecMag([x,y,z], [dx,dy,dz], [Mx,My,Mz])`
3. **Isotropic materials**: Use single-argument form `MatLin(ksi)` or equal susceptibilities `MatLin([ksi, ksi], [ex, ey, ez])`
4. **Easy axis**: For anisotropic materials, easy axis must be specified as 3D vector

**Example**:
```python
import radia as rad
rad.FldUnits('m')

# Soft iron cube (isotropic, μr=4000)
cube = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 0])  # Zero magnetization
mat = rad.MatLin(3999.0)  # χ = μr - 1 = 3999
rad.MatApl(cube, mat)

# Anisotropic material with easy axis in z-direction
cube2 = rad.ObjRecMag([0.2, 0, 0], [0.1, 0.1, 0.1], [0, 0, 0])
mat2 = rad.MatLin([5000, 100], [0, 0, 1])  # Easy axis along z
rad.MatApl(cube2, mat2)
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

## H-Matrix Solver Control

**Policy**: Users have full explicit control - no automatic problem size threshold.

```cpp
// CORRECT - Respect user's choice
if(use_hmatrix) {
    return SetupInteractMatrix_HMatrix();
}

// WRONG - Don't override user's choice
if(use_hmatrix && AmOfMainElem >= 200) {
    // Use H-matrix
} else {
    use_hmatrix = false;  // Override user's choice!
}
```

**Rationale**: User decides when H-matrix is appropriate for their use case.

### H-Matrix Implementation Policy

**IMPORTANT**: The current H-matrix implementation (`rad_hmatrix_aca.cpp/h`) is based on **HACApK concepts** and is the **only authorized H-matrix implementation** for Radia.

**Policy**:
- Do NOT implement alternative H-matrix libraries (HODLR, H2Lib, STRUMPACK, etc.)
- Do NOT replace the current ACA-based implementation with other low-rank approximation methods
- All H-matrix improvements must be made within the existing `radTHMatrixACA` class

**Current Implementation Details**:
- Based on HACApK (ppOpenHPC-MATH-HACApK) concepts from ELF_MAGIC
- Uses Adaptive Cross Approximation (ACA+) for low-rank block compression
- Parameters tuned for MMM: eta=2.0, min_cluster_size=15, eps=1e-5
- Integrated with BiCGSTAB solver (Method 10)

**Rationale**:
- HACApK is proven to work with MMM (Magnetic Moment Method) in ELF_MAGIC
- Single implementation reduces maintenance burden and complexity
- Alternative H-matrix libraries would require significant testing and validation

---

## NGSolve Integration Best Practices

**Recommended configuration**:
```python
fes = HDiv(mesh, order=2)  # Best accuracy
B_gf = GridFunction(fes)
B_gf.Set(radia_ngsolve.RadiaField(radia_obj, 'b'))
```

**Evaluation guidelines**:
- ✓ Evaluate GridFunction at distances > 1 mesh cell from magnet surface
- ✓ Use CoefficientFunction directly for maximum accuracy near boundaries
- ✗ Avoid GridFunction evaluation within 1 mesh cell of magnet surface

---

## PyPI Package Release Policy

### Version Management (Automated by Claude Code)

Claude Code is responsible for:
- Maintaining version numbers in `pyproject.toml`
- Following semantic versioning (MAJOR.MINOR.PATCH)
- Updating `CHANGELOG.md` with release notes

### PyPI Upload (Manual by User)

```powershell
# Set PyPI API token (keep secure!)
$env:PYPI_TOKEN = "pypi-AgEIcGl..."

# Run upload script
powershell.exe -ExecutionPolicy Bypass -File Publish_to_PyPI.ps1
```

**Security**: NEVER commit PyPI tokens to repository.

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
magnet = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 1.2])  # 0.1m

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

# All subsequent operations use this unit system
magnet = rad.ObjRecMag([0, 0, 0], [100, 100, 100], [0, 0, 1.2])  # 100mm
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

## Tetrahedral Element Method Control

### Requirement: Explicit API for Tetrahedral Field Computation Method

**Goal**: Allow users to explicitly select tetrahedral element field computation method via API, not environment variables.

**Rationale**:
- Environment variables are not user-friendly and hard to discover
- Explicit API makes method selection clear in code
- Matches pattern of H-matrix solver control (`rad.SolverHMatrixEnable()`)

**API**:

```python
import radia as rad

# Set tetrahedral method
rad.SolverTetraMethod(method)
# method = 0: Original Radia method (default)
# method = 1: Analytical method (for high-permeability materials)
```

**Implementation**:

1. **Global variable** (`src/lib/radentry.cpp`):
   ```cpp
   static int g_TetrahedronMethod = 0;  // 0=original, 1=analytical
   ```

2. **C API** (`src/lib/radentry.cpp`, `src/lib/radentry.h`):
   ```cpp
   int CALL RadSolverTetraMethod(int method);
   int RadSolverGetTetraMethod();
   ```

3. **Python binding** (`src/python/radpy_pyapi.cpp`):
   ```python
   rad.SolverTetraMethod(method)
   ```

**Usage Example**:

```python
import radia as rad
from netgen_mesh_import import netgen_mesh_to_radia

# Set tetrahedral method to analytical
rad.SolverTetraMethod(1)

# Import tetrahedral mesh
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 12000]},
                                units='m')

# Solve
rad.Solve(mag_obj, 0.0001, 10000)
```

**Method Details**:

| Method | Name | Description | Use Case |
|--------|------|-------------|----------|
| 0 | Original | Radia's original tetrahedral method | Default, general purpose |
| 1 | Analytical | Analytical method with basis vector extraction | High permeability materials (μr > 100) |

**Notes**:
- Default is method=0 (original Radia method)
- Method setting is global and affects all tetrahedral elements
- Must be set before calling `rad.Solve()`
- Does not affect hexahedral or other element types

---

**Last Updated**: 2025-11-27
**For**: Claude Code AI Assistant
**Project**: Radia Magnetic Field Computation
