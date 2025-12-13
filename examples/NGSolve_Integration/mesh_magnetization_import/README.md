# NGSolve Magnetization Import

This folder contains scripts for importing magnetization distributions into Radia
tetrahedral meshes and evaluating magnetic fields.

## NGSolve Version Requirement

**CRITICAL**: Use NGSolve **6.2.2405** only.

```bash
pip install ngsolve==6.2.2405
```

**Reason**: NGSolve 6.2.2406+ has a regression bug in Periodic Boundary Conditions.
The `Identify()` information is lost during mesh generation with `Glue()`.

**Reference**: https://forum.ngsolve.org/t/ngsolve-periodic-boundary-condition-regression-bug-report/3805

## Summary of Results

### Key Findings

1. **Radia tetrahedral MSC method is ACCURATE** (< 1% error vs analytical dipole)
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

**Conclusion**: Radia's MSC method correctly computes B field from tetrahedral meshes
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
cd examples/ngsolve_integration/mesh_magnetization_import
python sphere_analytical_to_radia.py
```

**Output**:
```
Field Comparison (Analytical dipole vs Radia MSC):
  Average error: 0.9832%
  Maximum error: 1.0311%
  [PASS] Radia MSC matches analytical dipole field (< 5% error)
```

## Key Implementation Notes

### 1. Use `netgen_mesh_import` Module

```python
from netgen_mesh_import import netgen_mesh_to_radia

rad.FldUnits('m')  # REQUIRED for NGSolve integration

mag_obj = netgen_mesh_to_radia(
    mesh,
    material={'magnetization': [0, 0, M_z]},
    units='m',
    material_filter='magnetic'
)
```

### 2. Radia Magnetization is in A/m (NOT Tesla)

```python
# Correct: A/m
magnetization = [0, 0, 1000.0]  # 1000 A/m

# Wrong: Tesla
# magnetization = [0, 0, 1.2]  # This would be interpreted as 1.2 A/m
```

### 3. Always Use `rad.FldUnits('m')` with NGSolve

Netgen uses meters, Radia defaults to millimeters. Always set:
```python
rad.FldUnits('m')
```

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

### Radia MSC (Magnetic Surface Charge) Method

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
