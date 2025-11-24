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
// ✓ CORRECT - Respect user's choice
if(use_hmatrix) {
    return SetupInteractMatrix_HMatrix();
}

// ✗ WRONG - Don't override user's choice
if(use_hmatrix && AmOfMainElem >= 200) {
    // Use H-matrix
} else {
    use_hmatrix = false;  // Override user's choice!
}
```

**Rationale**: User decides when H-matrix is appropriate for their use case.

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

**Last Updated**: 2025-11-23
**For**: Claude Code AI Assistant
**Project**: Radia Magnetic Field Computation
