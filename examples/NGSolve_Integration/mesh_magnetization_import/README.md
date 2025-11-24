# Mesh Import Tests for Radia

This folder contains tests for importing external mesh files into Radia for magnetic field computation.

## Overview

Radia supports three methods for creating magnetic geometry:

1. **Built-in primitives** (`rad.ObjRecMag`, `rad.ObjCylMag`, etc.) - Fastest, most accurate
2. **Tetrahedral mesh import** (from NGSolve/Netgen) - NGSolve integration
3. **Hexahedral mesh import** (from Cubit, Nastran format) - Professional mesh tools

This test suite verifies that all three methods produce consistent results when combined with `rad.MatLin` (linear materials) and `rad.ObjBckgCF` (background fields).

## Files

### NGSolve Reference Solutions

- `ngsolve_cube_uniform_field.py` - **High-precision NGSolve reference solution**
  - H-formulation with perturbation potential method
  - Magnetic cube (0.1m) with μᵣ = 100, H_ext = 1000 A/m
  - Multi-region mesh (magnetic cube + inner air + outer air)
  - Graded mesh refinement for accuracy
  - Saves results to .npz file for comparison
  - Based on: S:/ngsolve/NGSolve/2024_01_31_H-formulation/

- `compare_radia_ngsolve_cube.py` - **Compare Radia with NGSolve reference**
  - Loads NGSolve high-precision solution
  - Evaluates Radia ANALYTICAL method at same points
  - Reports field errors and convergence

### Demonstration Scripts

- `demo_tetrahedral_methods_comparison.py` - **NGSolve mesh import demo**
  - Generate tetrahedral mesh with NGSolve/Netgen
  - Import into Radia with ANALYTICAL method
  - Compare with built-in hexahedral reference

### Mesh Generation

- `generate_hex_mesh_cubit.py` - Generate hexahedral mesh using Cubit Python API
  - Creates 0.1m cube with 5×5×5 = 125 hex8 elements
  - Exports to Nastran (.bdf), Gmsh (.msh), VTK (.vtk) formats
  - Requires: Coreform Cubit 2025.3

### Mesh Files (Generated)

- `cube_hex.bdf` - Nastran format hexahedral mesh (primary)
- `cube_hex.msh` - Gmsh format (alternative)
- `cube_hex.vtk` - VTK format (for ParaView visualization)
- `cube_hex.cub5` - Cubit session file (for reference)

## Usage

### 1. Generate NGSolve High-Precision Reference Solution

```bash
cd S:/Radia/01_GitHub/examples/ngsolve_integration/mesh_magnetization_import
python ngsolve_cube_uniform_field.py
```

**Requirements:**
- NGSolve with Netgen installed
- Python packages: numpy, ngsolve

**Output:**
- `ngsolve_cube_uniform_field_results.npz` - Saved field data at test points
- Console output with field values and solver statistics

**What it does:**
- Creates 3D geometry: magnetic cube (μᵣ=100) + inner/outer air domains
- Generates graded tetrahedral mesh (fine near cube, coarse in far field)
- Solves H-formulation: ∇·(μ∇φ) = ∇·(μH_s) with perturbation potential
- Evaluates total field H_total = H_s + H_pert at test points
- Provides high-precision reference (CG solver with tol=1e-8)

### 2. Compare Radia with NGSolve Reference

```bash
python compare_radia_ngsolve_cube.py
```

**Requirements:**
- NGSolve reference solution must be generated first (step 1)
- Radia built with ANALYTICAL method support
- Environment: `RADIA_TETRA_METHOD='ANALYTICAL'` (set in script)

**Output:**
- Field comparison at identical test points
- Error percentages for magnitude and components
- Pass/fail assessment (target: <20% error)

### 3. Run NGSolve Mesh Import Demo

```bash
python demo_tetrahedral_methods_comparison.py
```

**What it does:**
- Generates tetrahedral mesh with NGSolve/Netgen
- Imports into Radia using ANALYTICAL method
- Compares with Radia's built-in hexahedral method

### 4. Generate Hexahedral Mesh (Optional, for Cubit users)

```bash
python generate_hex_mesh_cubit.py
```

**Requirements:**
- Coreform Cubit 2025.3 installed
- `S:/CoreformCubit/04_GitHub/cubit_mesh_export.py` available

**Output:**
- `cube_hex.bdf` (35KB) - 216 vertices, 125 hex elements
- `cube_hex.msh`, `cube_hex.vtk`, `cube_hex.cub5`

**Expected Results:**

| Method | |H| (A/m) | Error vs Built-in | Status |
|--------|----------|-------------------|--------|
| Built-in primitive | ~0.990 | 0% (reference) | ✓ |
| Hexahedral mesh | ~0.972 | ~2% | ✓ |
| Tetrahedral mesh | ~1.056 | ~7% | ✓ |

All three methods should produce similar field values, confirming that mesh import works correctly.

## Test Physics

### Problem Setup

- **Geometry**: 0.1m × 0.1m × 0.1m cube centered at origin
- **Material**: Linear magnetic material (μᵣ = 10, χ = 9)
- **Background field**: H₀ = 1 A/m in z-direction (uniform)
- **Test point**: [0.2, 0, 0] m (far field)

### Expected Behavior

Inside the cube:
- H_inside ≈ H₀ / μᵣ ≈ 0.1 A/m (field screening by high permeability)

Outside the cube:
- H_outside ≈ H₀ + perturbation
- Far field (~0.2m away): H ≈ 0.99 A/m (slight perturbation)

### Mesh Discretization Effects

- **Built-in**: Exact analytical geometry → most accurate
- **Hexahedral mesh**: 5×5×5 structured grid → ~2% error (excellent)
- **Tetrahedral mesh**: Unstructured mesh → ~7% error (acceptable)

Errors are due to:
1. Discretization of continuous geometry into finite elements
2. Numerical integration over element faces
3. Mesh resolution (finer mesh → better accuracy)

## Technical Details

### Hexahedral Mesh Import

Uses `nastran_mesh_import.py` to parse Nastran .bdf files:

```python
from nastran_mesh_import import create_radia_from_nastran

cube = create_radia_from_nastran(
    'cube_hex.bdf',
    material={'magnetization': [0, 0, 0]},
    units='m',
    combine=True
)
```

**Supported formats:**
- GRID/GRID* cards (node definitions, fixed-width and long format)
- CHEXA cards (8-node hexahedral elements with continuation lines)
- CTETRA cards (4-node tetrahedral elements)

**Node ordering:**
- Nastran CHEXA8 standard node ordering
- Automatically mapped to Radia HEX_FACES topology

### Tetrahedral Mesh Import

Uses `netgen_mesh_import.py` with NGSolve mesh:

```python
from netgen.occ import Box, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import create_radia_from_mesh

geo = OCCGeometry(Box((-0.05, -0.05, -0.05), (0.05, 0.05, 0.05)))
mesh = Mesh(geo.GenerateMesh(maxh=0.03))

cube = create_radia_from_mesh(
    mesh,
    material={'magnetization': [0, 0, 0]},
    combine=True
)
```

**Features:**
- Direct integration with NGSolve mesh
- Automatic face winding correction (outward normals)
- Supports both tetrahedral and hexahedral NGSolve meshes

## Troubleshooting

### Error: "cube_hex.bdf not found"

**Solution:** Run `generate_hex_mesh_cubit.py` first to create the mesh file.

### Error: "Can not find at least three vertex points..."

**Cause:** Duplicate vertices or incorrect node ordering in mesh file.

**Solution:** Check that:
1. GRID* long format is parsed correctly (16-character fields)
2. Continuation lines are handled properly ('+' markers)
3. All vertices are unique

### Error: Cubit module not found

**Cause:** Cubit Python API not in system path.

**Solution:** Ensure `sys.path.append("C:/Program Files/Coreform Cubit 2025.3/bin")` is correct.

### Large errors (>20%)

**Possible causes:**
1. Face winding direction incorrect (check TETRA_FACES/HEX_FACES)
2. Mesh too coarse (reduce `maxh` parameter)
3. Units mismatch (ensure `rad.FldUnits('m')` for NGSolve integration)

## Related Documentation

- `src/python/nastran_mesh_import.py` - Nastran mesh parser
- `src/python/netgen_mesh_import.py` - NGSolve mesh importer
- `examples/ngsolve_integration/h_formulation/` - H-formulation comparison (planned)
- `examples/ngsolve_integration/magnetization_import/` - Magnetization import (planned)

## Development Notes

### Nastran CHEXA8 Node Ordering

Nastran standard (ISO/IEC 10303-238):
```
     7 ----------- 6
    /|            /|
   / |           / |
  4 ----------- 5  |
  |  |          |  |
  |  3 ---------|--2
  | /           | /
  |/            |/
  0 ----------- 1
```

Nodes 0-3: Bottom face (z-)
Nodes 4-7: Top face (z+)

This ordering is automatically mapped to Radia's HEX_FACES topology by `nastran_mesh_import.py`.

### Version History

- v0.1.0 (2025-11-22): Initial hexahedral mesh import support
  - GRID* long format parsing
  - CHEXA continuation line handling
  - Cubit mesh generation script

---

**Author**: Radia Development Team
**Last Updated**: 2025-11-22
