# Electromagnet Example (mu_r = 1000)

This example validates Radia's MSC (Magnetic Surface Charge) hexahedron implementation against ELF_MAGIC.

## Model Description

- **Geometry**: C-type electromagnet with racetrack coil
- **Elements**: 52 hexahedral elements (full model), 26 elements (x-mirror)
- **DOF**: 312 (52 elements x 6 faces)
- **Material**: Linear soft iron, mu_r = 1000
- **Coil**: Racetrack coil, -2000 A

### Coil Parameters

| Parameter | Value |
|-----------|-------|
| Center | (0, 131.25, 0) mm |
| Inner radius | 5 mm |
| Outer radius | 40 mm |
| Straight section | 50 mm, 62.5 mm |
| Height | 105 mm |
| Current | -2000 A |

## Validation Results

### Interaction Matrix Comparison

**Radia vs ELF_MAGIC EIEM2 (full model, 52 elements)**

| Metric | Value |
|--------|-------|
| Max \|Radia - ELF\| | 7.99e-07 |
| Relative max diff | 0.0000% |
| RMS difference | 2.13e-08 |

**Result: Matrices match exactly.**

The diagonal block [1,1] comparison:

```
ELF:
  [ -7.4083  -2.7877  -6.7522  -2.7877  -3.1643  -3.1643]
  [ -1.0277  -0.6309  -1.0277  -0.4066  -0.4842  -0.4842]
  [ -6.7522  -2.7877  -7.4083  -2.7877  -3.1643  -3.1643]
  [ -1.0277  -0.4066  -1.0277  -0.6309  -0.4842  -0.4842]
  [ -1.3484  -0.5598  -1.3484  -0.5598  -0.8688  -0.6016]
  [ -1.3484  -0.5598  -1.3484  -0.5598  -0.6016  -0.8688]

Radia:
  [ -7.4083  -2.7877  -6.7522  -2.7877  -3.1643  -3.1643]
  [ -1.0277  -0.6309  -1.0277  -0.4066  -0.4842  -0.4842]
  [ -6.7522  -2.7877  -7.4083  -2.7877  -3.1643  -3.1643]
  [ -1.0277  -0.4066  -1.0277  -0.6309  -0.4842  -0.4842]
  [ -1.3484  -0.5598  -1.3484  -0.5598  -0.8688  -0.6016]
  [ -1.3484  -0.5598  -1.3484  -0.5598  -0.6016  -0.8688]
```

### Field Comparison at Origin

| Solver | Bz (mT) |
|--------|---------|
| ELF EIEM2 | -228.12 |
| ELF EIEM1 | -226.36 |
| Radia | -226.24 |

**Difference from ELF EIEM2**: 1.88 mT (0.82%)

### Cause of Field Difference

The interaction matrices are identical, so the field difference comes from **coil modeling**:

| Aspect | ELF_MAGIC | Radia |
|--------|-----------|-------|
| Coil model | Segmented (discretized) | Analytical racetrack |
| Arc sections | Approximated by straight segments | Exact arc integration |

The 0.82% difference is due to the coil discretization in ELF vs the analytical coil in Radia.

## EIEM Evaluation Points

Radia implements **EIEM2** (Yano-Sugahara method):

- **Evaluation point**: Midpoint between face center and element center
- **Formula**: `EvalPt = 0.5 * (FaceCenter + ElementCenter)`

| Method | Evaluation Point | ELF Bz (mT) |
|--------|------------------|-------------|
| EIEM0 | Element centroid | -924.33 |
| EIEM1 | Face center | -226.36 |
| EIEM2 | Midpoint (face center + element center) | -228.12 |
| EIEM3 | Different scheme | -226.12 |

## Directory Structure

```
mu=1000/
  README.md           # This file
  full/               # Full model (52 elements, no symmetry)
    compare_matrix_api.py    # Matrix comparison using GetInteractMatrix API
    compare_coil_field.py    # Coil field analysis
  x-mirror/           # X-mirror model (26 elements + TrfMlt)
    trace_eval_points.py     # Evaluation point tracing
    compare_matrices_direct.py
```

## Usage

### Compare interaction matrices

```python
python full/compare_matrix_api.py
```

### Analyze coil field differences

```python
python full/compare_coil_field.py
```

## TrfMlt REMOVED (2026-01-31)

**IMPORTANT**: `TrfMlt`, `TrfPlSym`, `TrfZerPara`, and `TrfZerPerp` have been **REMOVED** from Radia.

The shared-DOF design in TrfMlt was fundamentally incompatible with MSC 6DOF hexahedra. Element-based management (IMA) is the correct approach.

**Note**: The test scripts in this directory (`x-mirror/`, `z-mirror/`, etc.) reference the old TrfMlt API and will not work with current Radia. Use IMA symmetry instead.

## IMA (Image) Symmetry (2026-01-30)

**Radia now supports IMA symmetry for MSC hexahedra** - an alternative to TrfMlt that works correctly with plane symmetry.

### Usage

```python
import radia as rad

rad.FldUnits('m')

# Build full model geometry
hex_objects = [rad.ObjHexahedron(verts, [0,0,0]) for verts in all_vertices]
for h in hex_objects:
    rad.MatApl(h, rad.MatLin(mu_r))
container = rad.ObjCnt([coil] + hex_objects)

# Enable IMA x-mirror (half model)
intrc = rad.PreRelax(container, container)
n_ima = rad.SetIMASymmetry(intrc, 'x')  # Returns 26 (half of 52)
rad.BuildIMAMatrix(intrc)

# Solve with reduced DOF (156 instead of 312)
rad.Solve(container, 0.0001, 100, 0)
B = rad.Fld(container, 'b', [0, 0, 0])  # Same result as full model
```

### Validation Results

| Model | Elements | DOF | Bz (mT) |
|-------|----------|-----|---------|
| Full model | 52 | 312 | -226.24 |
| **IMA x-mirror** | 26 | 156 | **-226.24** |
| TrfMlt (broken) | 26 | 156 | Incorrect |

IMA produces **identical results** to the full model with half the DOF.

## Conclusions

1. **Radia's MSC interaction matrix matches ELF_MAGIC EIEM2 exactly** (full model)
2. **Field difference (0.82%) is due to coil modeling**, not MSC implementation
3. Radia's hexahedral MSC implementation is validated against ELF_MAGIC
4. **TrfMlt has been REMOVED** - use IMA symmetry for plane symmetry
5. **IMA symmetry works correctly** with MSC hexahedra for plane symmetry

## References

- ELF_MAGIC: `S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\`
- Yano-Sugahara MSC method (EIEM2 evaluation points)
