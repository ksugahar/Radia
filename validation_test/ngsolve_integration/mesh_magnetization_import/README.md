# NGSolve Magnetization Import

This folder contains scripts for importing magnetization distributions into Radia
tetrahedral meshes and evaluating magnetic fields.

## NGSolve Version Requirement

**CRITICAL**: Use NGSolve **6.2.2601** or later.

```bash
pip install ngsolve>=6.2.2601
```

**Note**: NGSolve 6.2.2406~6.2.2501 had a regression bug in Periodic Boundary Conditions
where `Identify()` information was lost during mesh generation with `Glue()`. This is fixed in 6.2.2601+.

**Reference**: https://forum.ngsolve.org/t/ngsolve-periodic-boundary-condition-regression-bug-report/3805

## Summary of Results

### Key Findings

1. **Radia fixed-magnetization tetrahedral field evaluation is accurate**
   (< 1% error vs analytical dipole)
2. **NGSolve H-formulation with Kelvin transform is ACCURATE** (< 0.15% error with correct version)
3. **Use `netgen_mesh_import` module** for correct Netgen -> Radia mesh transfer

### Verified Test: Analytical Magnetization -> Radia

| Metric | Result |
|--------|--------|
| Average Error | 0.98% |
| Maximum Error | 1.03% |
| Test Points | 12 (all outside sphere) |
| Mesh Elements | 4496 tetrahedra |
| Magnetization | Uniform M = [0, 0, 1000] A/m |

**Conclusion**: Radia's fixed-magnetization field evaluator correctly computes
B field from tetrahedral meshes
when given correct magnetization values.

### NGSolve H-formulation with Kelvin Transform

Using NGSolve 6.2.2405 with H-formulation and Kelvin transform for infinite domain:

| Location | NGSolve | Analytical | Error |
|----------|---------|------------|-------|
| Origin (0,0,0) | -0.970587 A/m | -0.970588 A/m | **0.000%** |
| (0.7, 0, 0) exterior | -0.353201 A/m | -0.353713 A/m | **0.145%** |

| Region | Max Error | RMS Error |
|--------|-----------|-----------|
| Interior (|x| < 0.5 m) | 3.32e-05 A/m | **0.001%** |
| Exterior (|x| >= 0.5 m) | 1.39e-03 A/m | 6.27e-04 A/m |

**Test Configuration**:
- Sphere radius: 0.5 m
- Kelvin radius: 1.0 m
- Relative permeability: mu_r = 100
- Mesh elements: 1,370,773
- DOFs: 6,310,708

### NGSolve Version Comparison

| Version | Origin Error | (0.7,0,0) Error | Status |
|---------|--------------|-----------------|--------|
| 6.2.2405 | 0.000% | 0.145% | **OK** |
| 6.2.2406 | 3.030% | 182.713% | **BUG** |

**Root Cause**: In 6.2.2406+, `Identify()` periodic BC information is lost during
mesh generation with `Glue()`, causing incorrect coupling between interior and
Kelvin-transformed exterior domains

## Files

### Working Scripts

- **`sphere_analytical_to_radia.py`** - **VERIFIED: < 1% error**
  - Creates sphere mesh with Netgen
  - Assigns analytical uniform magnetization
  - Compares Radia B field with dipole formula
  - Uses `netgen_mesh_to_radia()` for correct mesh transfer

### Reference Scripts

- `ngsolve_cube_uniform_field.py` - NGSolve H-formulation reference
- `verified_ngsolve_to_radia.py` - Full NGSolve -> Radia pipeline test

## Usage

### Verified Test (Analytical Magnetization)

```bash
cd validation_test/ngsolve_integration/mesh_magnetization_import
python sphere_analytical_to_radia.py
```

**Output**:
```
Field Comparison (Analytical dipole vs Radia fixed-magnet field):
  Average error: 0.9832%
  Maximum error: 1.0311%
  [PASS] Radia fixed-magnet field matches analytical dipole field (< 5% error)
```

## Key Implementation Notes

### 1. Use `netgen_mesh_import` Module (MANDATORY)

**CRITICAL POLICY**: All NGSolve mesh access MUST use functions from `netgen_mesh_import.py`.

**Why?** NGSolve has TWO different indexing schemes that cause bugs:

| Access Method | Indexing | Valid Range |
|--------------|----------|-------------|
| `mesh.ngmesh.Points()[i]` | **1-indexed** | 1 to nv |
| `mesh.vertices[i]` | **0-indexed** | 0 to nv-1 |
| `el.vertices[i].nr` | Returns **0-indexed** value | Use with `mesh.vertices[]` only |

**Common Bug Pattern (DO NOT DO THIS):**
```python
# WRONG - Using 0-indexed .nr with 1-indexed ngmesh.Points()
for v in el.vertices:
    pt = mesh.ngmesh.Points()[v.nr]  # Off-by-one error!
```

**Correct Usage:**
```python
from netgen_mesh_import import netgen_mesh_to_radia, extract_elements, compute_element_centroid

# Radia always uses meters (compatible with NGSolve SI units)

# Option 1: Direct conversion (recommended)
mag_obj = netgen_mesh_to_radia(
    mesh,
    material={'magnetization': [0, 0, M_z]},
    units='m',
    material_filter='magnetic'
)

# Option 2: Custom processing with extract_elements
elements, _ = extract_elements(mesh, material_filter='magnetic')
for el in elements:
    vertices = el['vertices']  # Correctly extracted
    centroid = compute_element_centroid(vertices)
```

**Available Functions in `netgen_mesh_import.py`:**
- `netgen_mesh_to_radia()`: Convert entire mesh to Radia geometry (recommended)
- `extract_elements()`: Extract element data for custom processing
- `compute_element_centroid()`: Compute centroid from vertex list
- `create_radia_tetrahedron()`: Create single Radia tetrahedron
- `create_radia_hexahedron()`: Create single Radia hexahedron
- `TETRA_FACES`, `HEX_FACES`, `WEDGE_FACES`, `PYRAMID_FACES`: Face topology constants

### 2. Radia Magnetization is in A/m (NOT Tesla)

```python
# Correct: A/m
magnetization = [0, 0, 1000.0]  # 1000 A/m

# Wrong: Tesla
# magnetization = [0, 0, 1.2]  # This would be interpreted as 1.2 A/m
```

### 3. Radia Always Uses Meters

Radia always uses meters, which is compatible with NGSolve/Netgen SI units.

## Physics

### Analytical Solution for Uniformly Magnetized Sphere

Outside sphere (dipole field):
```
B = (mu_0/4pi) * (3(m*r)r/r^5 - m/r^3)
where m = (4/3)*pi*a^3*M is the magnetic dipole moment
```

Inside sphere:
```
B = (2/3)*mu_0*M
```

### Radia Fixed-Magnet Surface-Charge Field

Radia computes B field from surface charge density:
```
sigma = M dot n  (on each face)
```
Using closed-form solid angle integration for each triangular face.

## Related Folders

- `../../tetra_field_accuracy_evaluation/` - More tetrahedral accuracy tests
- `../../cube_uniform_field/` - Radia solver benchmarks

---

**Date**: 2025-12-13
**Version**: Radia v1.3.14+
