# Nonlinear Material Benchmark

This benchmark tests Radia's nonlinear solver with saturable BH curves.

**Reference**: ELF/MAGIC benchmark at `S:\ELF_MAGIC\01_GitHub\examples\cube_uniform_field\nonlinear`

## Problem Setup

- **Geometry**: 1.0 m x 1.0 m x 1.0 m cube centered at origin
- **Mesh**: 10x10x10 hexahedral elements (1000 elements)
- **Material**: Nonlinear BH curve (soft iron with saturation)
- **External Field**: Hz = 50,000 A/m (uniform field along z-axis)
- **Nonlinear Solver**: Newton-Raphson method

## BH Curve (Soft Iron)

```
H [A/m]     B [T]
0           0.0
100         0.1
200         0.3
500         0.8
1000        1.2
2000        1.5
5000        1.7
10000       1.8
50000       2.0
100000      2.1
```

- Initial permeability: mu_r ~= 800
- Saturation onset: B ~= 1.5 T

## BH Curve Definition in Radia

The nonlinear material is defined using `rad.MatSatIsoTab()` with [H, M] pairs:

```python
import numpy as np
import radia as rad

# BH curve data: [H (A/m), B (T)]
bh_data = [
    [0.0, 0.0],
    [100.0, 0.1],
    [200.0, 0.3],
    [500.0, 0.8],
    [1000.0, 1.2],
    [2000.0, 1.5],
    [5000.0, 1.7],
    [10000.0, 1.8],
    [50000.0, 2.0],
    [100000.0, 2.1],
]

# Convert to [H, M] format
# B = mu_0 * (H + M), so M = B/mu_0 - H
mu_0 = 4 * np.pi * 1e-7
hm_data = [[h, b/mu_0 - h] for h, b in bh_data]

mat = rad.MatSatIsoTab(hm_data)
```

## Reference Results

### Radia Results (this benchmark)

From Radia benchmark with 6x6x6 mesh, Hz = 50,000 A/m:

| Material | Time [s] | Iterations | Mz_avg [A/m] |
|----------|----------|------------|--------------|
| Linear (mu_r=800) | 0.03 | 118 | 157,168 |
| Nonlinear (BH curve) | 0.03 | 116 | 156,815 |

**Key Result**: Linear and nonlinear materials give similar Mz because the
BH curve initial slope corresponds to mu_r ~= 800.

### ELF/MAGIC Reference

From ELF/MAGIC benchmark with 10x10x10 mesh, Hz = 50,000 A/m:

| Solver | Time [s] | Newton Iter | Mz_center [A/m] | Residual | Status |
|--------|----------|-------------|-----------------|----------|--------|
| Dense LU | 11.30 | 8 | 341,651 | 8.2e-4 | **OK** |
| BiCGSTAB+H-matrix | 15.08 | 8 | 341,716 | 8.2e-4 | **OK** |

**Note**: ELF/MAGIC shows higher Mz due to differences in:
1. Newton-Raphson solver vs Radia's relaxation method
2. Mesh refinement (10x10x10 vs 6x6x6)
3. BH curve interpolation method

## H-Matrix Acceleration

For large-scale nonlinear problems (> 3000 elements), H-matrix acceleration significantly reduces computation time:

```python
rad.SolverHMatrixEnable(1)
```

### Crossover Point

- **< 3000 elements**: Dense LU is faster
- **> 3000 elements**: H-matrix is faster
  - N=15 (3,375 elements): H-matrix 1.29x faster
  - N=20 (8,000 elements): H-matrix 1.87x faster

## Running the Benchmark

```bash
cd examples/cube_uniform_field/nonlinear
python cube_nonlinear_hmatrix.py
```

## Expected Results

With 50,000 A/m external field:
- Interior H field: Hz ~= 50,000 A/m
- Interior magnetization: Mz ~= 157,000 A/m
- Interior B field: Bz = mu_0 * (H + M) ~= 0.2 T
- For the given BH curve, this is below saturation (B_sat ~= 2.0 T)

## Tetrahedral Mesh Support (2025-12-03)

Radia now supports tetrahedral meshes for nonlinear materials. Key findings:

### Solver Method Requirements

| Mesh Type | Method 0/4 (Relaxation) | Method 9 (LU) | Method 10 (BiCGSTAB) |
|-----------|-------------------------|---------------|----------------------|
| Hexahedral | OK | OK | OK |
| Tetrahedral | **NOT CONVERGED** | OK | OK |

**Important**: For tetrahedral meshes with nonlinear materials, use Method 9 or 10:

```python
# Method 10 (BiCGSTAB) - recommended for large problems
res = rad.Solve(grp, 0.001, 1000, 10)

# Method 9 (Direct LU) - recommended for small problems
res = rad.Solve(grp, 0.001, 1000, 9)
```

### Mesh Convergence Study

External field: B0 = 1.0 T, evaluation point: z = 60mm (outside cube)

| Mesh Type | Elements | Iterations | Bz_pert (mT) | Error vs Hex 8^3 |
|-----------|----------|------------|--------------|------------------|
| **Hexahedral** |
| Hex 4x4x4 | 64 | 7 | -391.7 | 7.7% |
| Hex 6x6x6 | 216 | 8 | -363.1 | 0.2% |
| Hex 8x8x8 | 512 | 10 | -363.7 | (reference) |
| **Tetrahedral (Netgen)** |
| maxh=60mm | 28 | 9 | -343.6 | 5.5% |
| maxh=40mm | 104 | 15 | -382.9 | 5.3% |
| maxh=30mm | 201 | 14 | -357.6 | 1.7% |
| maxh=25mm | 390 | 15 | -364.9 | 0.35% |
| maxh=20mm | 642 | 17 | -367.3 | 1.0% |

**Key Finding**: Tetrahedral mesh with ~400 elements achieves <1% accuracy compared to hexahedral reference.

### Tetrahedral Usage Example

```python
import radia as rad
from netgen.occ import Box, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

rad.FldUnits('m')
rad.SolverTetraMethod(0)  # Use polygon-based method

# Create tetrahedral mesh
cube_geom = Box((-0.05, -0.05, -0.05), (0.05, 0.05, 0.05))
geo = OCCGeometry(cube_geom)
ngmesh = Mesh(geo.GenerateMesh(maxh=0.025))  # ~400 elements

# Import to Radia
cube = netgen_mesh_to_radia(ngmesh, material={'magnetization': [0,0,0]}, units='m')

# Apply nonlinear material
mat = rad.MatSatIsoTab(hm_data)  # [H, M] format!
rad.MatApl(cube, mat)

# External field and solve
bg = rad.ObjBckg([0.0, 0.0, 1.0])  # 1 T
grp = rad.ObjCnt([cube, bg])
res = rad.Solve(grp, 0.001, 1000, 10)  # Method 10 (BiCGSTAB)
```

## Notes

1. **Strong external field required**: 50,000 A/m ensures saturation effects are visible
2. **Newton-Raphson**: Required for nonlinear problems (RLXM relaxation is insufficient)
3. **H-matrix recommendation**: Use for N > 3000 elements for best performance
4. **Tetrahedral solver**: Use Method 9 (LU) or Method 10 (BiCGSTAB) for tetrahedral nonlinear

## Files

- `cube_nonlinear_hmatrix.py` - Material comparison benchmark (linear vs nonlinear)
- `benchmark_nonlinear_tetra_vs_hex.py` - Tetrahedral vs hexahedral comparison (2025-12-03)
- `README.md` - This documentation

## Physics Notes

### Nonlinear Iteration

For nonlinear materials, Radia uses an iterative relaxation method:

1. Initialize magnetization M = 0
2. Compute H field from current M distribution
3. Update M using material BH curve: M = M(H)
4. Repeat until |delta M| < tolerance

### Saturation Effects

At high H fields, the material saturates:
- Low H: B ~ mu_r * mu_0 * H (linear regime)
- High H: B -> B_sat (saturation regime)
- Saturation magnetization for iron: ~2.1 T
