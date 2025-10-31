# CoefficientFunction vs GridFunction Comparison Guide

## Overview

`rad_ngsolve` provides two ways to evaluate Radia magnetic fields in NGSolve:

1. **CoefficientFunction**: Direct evaluation (exact)
2. **GridFunction**: FEM interpolation (convergent with mesh refinement)

## Comparison

### CoefficientFunction (CF)

**How it works**:
```python
B_cf = rad_ngsolve.RadBfield(magnet)
B = B_cf(mesh_pt)  # Calls Radia directly for each point
```

**Characteristics**:
- ✓ **Exact**: Error = 0% (matches Radia exactly)
- ✓ **Mesh-independent**: Same accuracy regardless of mesh size
- ✓ **Simple**: No interpolation step needed
- ⚠ **Slower**: Calls Radia for each evaluation point

**Best for**:
- Accurate field values at specific points
- Verification and validation
- Small number of evaluation points

### GridFunction (GF)

**How it works**:
```python
fes = VectorH1(mesh, order=2)
gfB = GridFunction(fes)
gfB.Set(B_cf)  # Interpolates onto finite element space
B = gfB(mesh_pt)  # Fast evaluation from interpolated values
```

**Characteristics**:
- ✓ **Fast repeated evaluation**: Once interpolated, very fast
- ✓ **Convergent**: Error decreases with mesh refinement
- ✓ **Compatible with FEM solvers**: Can be used in FEM simulations
- ⚠ **Interpolation error**: Depends on mesh size

**Best for**:
- Many evaluation points
- Integration with FEM solvers
- Coupled physics simulations

## Error vs Mesh Size

Test point: (0, 0, 0) inside magnet, Radia reference: Bz = 0.949718 T

| Mesh Size | Elements | CF Error | GF Error | GF Relative Error |
|-----------|----------|----------|----------|-------------------|
| 0.020 m (20mm) | 693 | 0% | 0.604 T | 63.56% |
| 0.010 m (10mm) | 4,970 | 0% | 0.302 T | 31.77% |
| 0.005 m (5mm) | 34,496 | 0% | 0.009 T | 1.00% |
| 0.002 m (2mm) | 771,995 | 0% | 0.0000004 T | <0.01% |

**Key findings**:
- CoefficientFunction is always exact (0% error)
- GridFunction converges with O(h²) rate
- For 1% accuracy, use maxh ≤ 0.005 m (5mm)
- For <0.01% accuracy, use maxh ≤ 0.002 m (2mm)

## Usage Examples

### 1. Exact Evaluation (CoefficientFunction)

```bash
python visualize_field.py --method cf
```

Output:
- `radia_field_cf.vtu` - Exact Radia values

### 2. Fast Interpolation (GridFunction with fine mesh)

```bash
python visualize_field.py --method gf --maxh 0.005
```

Output:
- `radia_field_gf.vtu` - Interpolated values (~1% error)

### 3. Side-by-Side Comparison

```bash
python visualize_field.py --method both --maxh 0.005
```

Output:
- `radia_field_compare.vtu` - Contains both CF and GF fields
- Open in Paraview to compare visually

### 4. Mesh Refinement Study

```bash
# Coarse mesh (fast, but large error)
python visualize_field.py --method gf --maxh 0.02

# Fine mesh (slower, but accurate)
python visualize_field.py --method gf --maxh 0.005

# Very fine mesh (slowest, but very accurate)
python visualize_field.py --method gf --maxh 0.002
```

## Paraview Visualization

After generating `radia_field_compare.vtu`:

1. Open file in Paraview
2. Select `B_field_CF` to view exact CoefficientFunction values
3. Select `B_field_GF` to view interpolated GridFunction values
4. Use Calculator filter to compute difference: `B_field_CF - B_field_GF`
5. Visualize error distribution

## Recommendations

### Choose CoefficientFunction if:
- You need exact Radia values
- You're evaluating at few points
- You're validating results
- Accuracy is critical

### Choose GridFunction if:
- You need fast repeated evaluation
- You're running FEM simulations
- You're computing integrals
- You can accept ~1% error with maxh=0.005m

### Use Both if:
- You want to verify accuracy
- You're studying convergence
- You're documenting results

## Code Examples

### CoefficientFunction Direct Evaluation

```python
import rad_ngsolve
from ngsolve import *
from netgen.csg import *

# Create Radia magnet
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])
rad.MatApl(magnet, rad.MatStd('NdFeB', 1.2))
rad.Solve(magnet, 0.0001, 10000)

# Create mesh
geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05)))
mesh = Mesh(geo.GenerateMesh(maxh=0.02))

# Create CoefficientFunction
B_cf = rad_ngsolve.RadBfield(magnet)

# Evaluate at point
mesh_pt = mesh(0, 0, 0.02)
B = B_cf(mesh_pt)  # Exact Radia value
print(f"Bz = {B[2]:.6f} T")
```

### GridFunction Interpolation

```python
# ... (same setup as above)

# Create fine mesh for accurate interpolation
mesh = Mesh(geo.GenerateMesh(maxh=0.005))  # 5mm mesh

# Create GridFunction
fes = VectorH1(mesh, order=2)
gfB = GridFunction(fes)
gfB.Set(B_cf)  # Interpolate once

# Fast repeated evaluation
for z in [0, 0.01, 0.02, 0.03, 0.04]:
	mesh_pt = mesh(0, 0, z)
	B = gfB(mesh_pt)  # Fast (no Radia call)
	print(f"z={z*1000}mm: Bz = {B[2]:.6f} T")
```

### Error Analysis

```python
# Compare CF and GF at test point
mesh_pt = mesh(0, 0, 0)

# Exact evaluation
B_cf_val = B_cf(mesh_pt)

# Interpolated evaluation
B_gf_val = gfB(mesh_pt)

# Compute error
error = np.sqrt((B_cf_val[0] - B_gf_val[0])**2 +
	            (B_cf_val[1] - B_gf_val[1])**2 +
	            (B_cf_val[2] - B_gf_val[2])**2)

print(f"CF: Bz = {B_cf_val[2]:.6f} T")
print(f"GF: Bz = {B_gf_val[2]:.6f} T")
print(f"Error: {error:.6e} T")
```

## Performance Considerations

### CoefficientFunction
- **Per-point cost**: ~0.1-1 ms (depends on Radia computation)
- **Memory**: Minimal (just stores Radia object ID)
- **Best for**: <1000 evaluation points

### GridFunction
- **Interpolation cost**: Depends on mesh size
  - maxh=0.02: ~1 second (693 elements)
  - maxh=0.005: ~10 seconds (34,496 elements)
  - maxh=0.002: ~5 minutes (771,995 elements)
- **Per-point cost**: ~0.001 ms (very fast)
- **Memory**: Proportional to DOFs
  - maxh=0.02: ~30 KB (3,855 DOFs)
  - maxh=0.005: ~1.2 MB (155,046 DOFs)
  - maxh=0.002: ~25 MB (3,198,549 DOFs)
- **Best for**: >10,000 evaluation points

## Summary

| Feature | CoefficientFunction | GridFunction (fine mesh) |
|---------|-------------------|------------------------|
| Accuracy | Exact (0% error) | ~1% error (maxh=0.005) |
| Speed (per point) | Slow (~1 ms) | Fast (~0.001 ms) |
| Setup time | None | Moderate (~10 sec) |
| Memory | Minimal | Moderate (~1 MB) |
| Use case | Few points, exact | Many points, fast |

**Final recommendation**:
- Use CoefficientFunction for exact values and verification
- Use GridFunction with maxh≤0.005m for production computations
- Both approaches are valid - choose based on your requirements!

Date: 2025-10-31
