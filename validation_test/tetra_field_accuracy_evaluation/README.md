# Tetrahedral Field Accuracy Validation

This folder contains validation-class scripts to evaluate the accuracy of
magnetic field computation from tetrahedral meshes in Radia using the MSC
(Magnetic Surface Charge) method. The human-facing notebook is
`docs/tetra_field_accuracy_evaluation/tetra_field_accuracy_validation.ipynb`.

## Summary of Results

### Key Finding: Tetrahedral MSC Field Computation is ACCURATE

The evaluation confirms that `rad.Fld()` from tetrahedral (ObjTetrahedron) elements
produces **identical results** to hexahedral (ObjHexahedron) elements when given
the same magnetization distribution.

| Metric | Result |
|--------|--------|
| Average Error | 0.0000% |
| Maximum Error | 0.0000% |
| Test Points | 13 (all outside magnet) |
| Tetra Elements | 390 (maxh=0.25m) |
| Hexa Elements | 64 (ObjHexahedron) |

## Test Configuration

- **Geometry**: 1.0m x 1.0m x 1.0m cube centered at origin
- **Magnetization**: Uniform M = [0, 0, 1.0e6] A/m
- **Test Points**: On-axis (z=0.6 to 3.0m) and off-axis locations
- **Reference**: Hexahedral mesh using ObjHexahedron (MSC method)

## Scripts

### [analytical_reference.py](analytical_reference.py)

Compares tetrahedral MSC field computation with hexahedral (analytical) reference
using **uniform magnetization** (permanent magnet, no solver).

**Usage:**
```bash
python analytical_reference.py
```

**Output:**
- Field values at 13 test points
- Error comparison (tetra vs hexa)
- Results saved to `analytical_reference_results.json`

**Result:** 0% error - Tetrahedral and hexahedral produce identical B fields.

### [evaluate_radia_solvers.py](evaluate_radia_solvers.py)

Compares tetrahedral (Netgen) and hexahedral (ObjHexahedron) meshes with
**linear material and external field** - both solved using Radia's BiCGSTAB solver.

**Usage:**
```bash
python evaluate_radia_solvers.py
```

**Output:**
- Magnetization from Radia solver (M_avg_z)
- Field comparison at 13 test points
- Results saved to `solver_comparison_results.json`

**Result:** < 5% error - Both mesh types solve correctly and produce similar B fields.

| Mesh Type | Elements | M_avg_z (A/m) | Field Error vs Hexa |
|-----------|----------|---------------|---------------------|
| Hexahedral | 125 (n_div=5) | 173,357 | - (reference) |
| Tetrahedral | 200 (maxh=0.3) | 187,788 | 2.8% average, 5.0% max |

### [evaluate_tetra_field.py](evaluate_tetra_field.py)

Evaluates tetrahedral field accuracy using NGSolve H-formulation magnetization.

**Note:** This script has known issues with NGSolve magnetization extraction.
Use `evaluate_radia_solvers.py` for Radia solver comparison.

### [ngsolve_reference.py](ngsolve_reference.py)

Computes the magnetic field from a uniformly magnetized cube using NGSolve's A-formulation (vector potential) as an independent reference solution.

### [compare_all.py](compare_all.py)

Loads results from `evaluate_tetra_field.py` and `ngsolve_reference.py` and creates a comprehensive comparison table of Radia tetrahedral MSC, hexahedral, and NGSolve reference solutions.

## Conclusions

1. **rad.Fld from tetrahedral mesh is ACCURATE**
   - Field values match hexahedral reference exactly (0% error for uniform M)
   - MSC method implementation is correct

2. **Radia solver works with tetrahedral meshes**
   - Both LU and BiCGSTAB converge for linear materials
   - Field error < 5% compared to hexahedral reference
   - Mesh discretization causes small differences in magnetization

3. **Both tetrahedral (Netgen) and hexahedral (ObjHexahedron) are validated**
   - Tetrahedral: Use `rad.ObjTetrahedron()` or `netgen_mesh_import.netgen_mesh_to_radia()`
   - Hexahedral: Use `rad.ObjHexahedron()` for simple geometries

## Related Benchmarks

- `examples/cube_uniform_field/linear/`: Linear material benchmarks
- `examples/cube_uniform_field/nonlinear/`: Nonlinear material benchmarks

## Technical Details

### MSC Method Implementation

The tetrahedral field computation uses the Magnetic Surface Charge (MSC) method:

1. Each tetrahedral element has 4 triangular faces
2. Surface charge density: sigma = M dot n
3. Field computed using solid angle integration formula
4. Implementation in `src/core/rad_polyhedron.cpp`

### NGSolve Mesh Access Policy (MANDATORY)

**CRITICAL**: All NGSolve mesh access MUST use functions from `netgen_mesh_import.py`.

**Why?** NGSolve has two different indexing schemes that cause off-by-one errors:

| Access Method | Indexing | Notes |
|--------------|----------|-------|
| `mesh.ngmesh.Points()[i]` | **1-indexed** | Index 0 raises error |
| `mesh.vertices[i]` | **0-indexed** | Valid: 0 to nv-1 |
| `el.vertices[i].nr` | Returns **0-indexed** | Use with `mesh.vertices[]` only |

**Correct Usage:**
```python
from netgen_mesh_import import extract_elements, compute_element_centroid

# Extract elements with correct indexing
elements, _ = extract_elements(mesh, material_filter='magnetic')
for el in elements:
    vertices = el['vertices']  # Already correctly extracted
    centroid = compute_element_centroid(vertices)
```

**NEVER** directly access `mesh.ngmesh.Points()`, `mesh.vertices[]`, or `el.vertices[].nr`.

### Face Topology (TETRA_FACES)

```python
from netgen_mesh_import import TETRA_FACES

TETRA_FACES = [
    [1, 3, 2],  # Face 0: v0-v2-v1
    [1, 2, 4],  # Face 1: v0-v1-v3
    [2, 3, 4],  # Face 2: v1-v2-v3
    [3, 1, 4]   # Face 3: v2-v0-v3
]
```

Note: 1-indexed for Radia, with reversed winding for outward normals.

---

**Date**: 2025-12-13
**Version**: Radia v1.3.14+
