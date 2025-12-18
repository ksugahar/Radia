# Claude Code - Radia Project Development Guidelines

This document contains development guidelines and refactoring policies for the Radia project when working with Claude Code.

## Radia Solver Method: MMM (Magnetic Moment Method)

**IMPORTANT**: Radia uses the **MMM (Magnetic Moment Method)**, NOT BEM (Boundary Element Method).

**Terminology**:
- Correct: "Radia MMM solver", "MMM computation", "MMM field evaluation"
- Incorrect: "Radia BEM", "BEM solver", "boundary element method"

**Context**:
- MMM represents magnetic objects as distributions of magnetic moments
- This differs from BEM which uses surface integral equations
- All references to "BEM" in Radia-related documentation should be corrected to "MMM"

---

## Polyhedral Element Field Computation (MSC)

### Magnetic Surface Charge (MSC) Method

Both tetrahedral and hexahedral polyhedral elements use the **Magnetic Surface Charge (MSC)** method for field computation. This method computes the magnetic field from each triangular face using closed-form formulas based on the surface magnetic charge density (sigma = M dot n).

**Supported Element Types**:
- **Tetrahedra (4 faces)**: 4 triangular faces
- **Hexahedra (6 faces)**: 6 quadrilateral faces, each split into 2 triangles

**Key Features**:
- Uses **global coordinates** directly (no local coordinate transformations)
- Computes field using **solid angle integration** formula
- Handles **outward normal orientation** per Radia convention

**Implementation** (`src/core/rad_polyhedron.cpp`):
- `B_comp_tetrahedron_MSC()`: 4-face tetrahedral field computation
- `B_comp_hexahedron_MSC()`: 6-face hexahedral field computation (quad -> 2 triangles)
- Uses `RadFieldFromTriangleFaceGlobal()` from `rad_poly_analytical.cpp`

**Usage**:

```python
import radia as rad
from netgen_mesh_import import netgen_mesh_to_radia

rad.FldUnits("m")

# Import tetrahedral mesh (automatically uses MSC method)
mag_obj = netgen_mesh_to_radia(mesh,
                                material={"magnetization": [0, 0, 0]},
                                units="m")

# Apply material and solve
mat = rad.MatLin(999.0)  # mu_r = 1000
rad.MatApl(mag_obj, mat)
rad.Solve(mag_obj, 0.0001, 1000)
```

**Notes**:
- Both tetrahedral and hexahedral meshes work with linear and nonlinear materials
- LU solver (Method 0) and BiCGSTAB (Method 1) both work with polyhedral elements

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

## H-Matrix Policy: Do NOT Implement Custom Algorithms

### Policy

**CRITICAL**: Do NOT implement custom H-matrix algorithms (ACA, ACA+, or any low-rank approximation).

**Rules**:
1. **No custom H-matrix code**: Never write new ACA/H-matrix implementations
2. **Use HACApK only**: If H-matrix is ever needed, use `src/ext/HACApK_LH-Cimplm/`
3. **Removed implementations**: rad_hmatrix*.cpp/h files were deleted (2025-12-18)

**Rationale**:
- Custom implementations are prone to bugs and difficult to validate
- HACApK is a proven, MIT-licensed library
- Benchmarks showed NO speedup for typical Radia use cases (single compact objects)

### radia_ngsolve.RadiaField

The `radia_ngsolve.RadiaField` CoefficientFunction uses **direct rad.Fld() calls**:
- Batch evaluation via `rad.Fld(obj, field_type, points_list)`
- Point caching with `PrepareCache()` for repeated evaluations
- **NO H-matrix acceleration** is implemented or planned

**See**: Main `CLAUDE.md` for full H-matrix policy details.

---

**Last Updated**: 2025-12-19
**For**: Claude Code AI Assistant
**Project**: Radia Magnetic Field Computation
