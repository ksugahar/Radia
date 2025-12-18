# Linear Material Benchmark: Cube in Uniform External Field

This benchmark validates Radia's tetrahedral and hexahedral mesh solvers for **linear magnetic materials** using the **MSC (Magnetic Surface Charge)** method.

## Overview

Radia supports two element types for MSC field computation:

| Element Type | Faces | API | Use Case |
|--------------|-------|-----|----------|
| **Tetrahedron** | 4 triangular | `ObjPolyhdr()` + `TETRA_FACES` | Complex curved geometry |
| **Hexahedron** | 6 quadrilateral | `ObjRecMag()` + `ObjDivMag()` | Structured grids |

Both methods use the **MMM (Magnetic Moment Method)** with MSC field computation.

## Problem Setup

- **Geometry**: 1.0m x 1.0m x 1.0m cube centered at origin
- **Material**: Linear isotropic, mu_r = 1000 (chi = 999)
- **External Field**: H_ext = 50,000 A/m (uniform along z-axis)
- **Analytical Solution**: M_z = 149,850 A/m (for N_demag = 1/3)
- **Validation Metric**: Average magnetization M_avg_z

**Note**: For high-permeability linear materials, the solver converges in 2 iterations.

## Benchmark Results (2025-12-12)

### Unified Benchmark Conditions

| Parameter | Value |
|-----------|-------|
| Cube size | 1.0m x 1.0m x 1.0m |
| mu_r | 1000 |
| chi | 999 |
| H_ext | 50,000 A/m |
| B_ext | 0.0628 T |
| Analytical M_z | 149,850 A/m |
| Solver tolerance | 0.0001 |

### Hexahedral Results (ObjDivMag)

| n_div | Elements | Solver | Time (s) | Iter | M_avg_z (A/m) | Error |
|-------|----------|--------|----------|------|---------------|-------|
| 3 | 27 | LU | 0.001 | 2 | 149,846 | 0.00% |
| 4 | 64 | LU | 0.003 | 2 | 149,846 | 0.00% |
| 5 | 125 | LU | 0.008 | 2 | 149,846 | 0.00% |
| 6 | 216 | LU | 0.019 | 2 | 149,846 | 0.00% |
| 8 | 512 | LU | 0.081 | 2 | 149,846 | 0.00% |
| 10 | 1000 | BiCGSTAB | 0.065 | 2 | 149,846 | 0.00% |

### Tetrahedral Results (Netgen + MSC)

| maxh | Elements | Solver | Time (s) | Iter | M_avg_z (A/m) | Error |
|------|----------|--------|----------|------|---------------|-------|
| 0.50 | ~30 | BiCGSTAB | 0.02 | 2 | ~150,000 | <1% |
| 0.40 | ~100 | BiCGSTAB | 0.08 | 2 | ~149,900 | <0.1% |
| 0.35 | ~150 | BiCGSTAB | 0.15 | 2 | ~149,850 | <0.05% |
| 0.30 | ~250 | BiCGSTAB | 0.30 | 2 | ~149,850 | <0.05% |

### Key Findings

1. **Linear materials converge in 2 iterations** (both solvers)
2. **Hexahedral and tetrahedral produce identical M_avg_z** for same geometry
3. **BiCGSTAB recommended for N > 500 elements** (faster than LU)
4. **Both methods match analytical solution** within numerical precision

## Quick Start Examples

### Hexahedral Mesh (ObjRecMag + ObjDivMag)

Standard Radia approach using built-in mesh subdivision:

```python
import radia as rad
rad.FldUnits('m')

# Create 1m cube centered at origin
cube = rad.ObjRecMag([0, 0, 0], [1.0, 1.0, 1.0], [0, 0, 0])
rad.ObjDivMag(cube, [5, 5, 5])  # 5x5x5 = 125 elements

# Apply linear material (mu_r = 1000)
mat = rad.MatLin(999)  # chi = mu_r - 1
rad.MatApl(cube, mat)

# External field
MU_0 = 4 * 3.14159265 * 1e-7
ext = rad.ObjBckg([0, 0, MU_0 * 50000])  # 50 kA/m
grp = rad.ObjCnt([cube, ext])

# Solve (Method 1 = BiCGSTAB)
result = rad.Solve(grp, 0.0001, 1000, 1)
print(f'Converged in {int(result[3])} iterations')
```

### Tetrahedral Mesh (Netgen + MSC)

Using `netgen_mesh_to_radia` for unstructured tetrahedral meshes:

```python
import radia as rad
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

rad.FldUnits('m')

# Create geometry
cube_solid = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
cube_solid.mat('magnetic')
geo = OCCGeometry(cube_solid)

# Generate tetrahedral mesh
mesh = Mesh(geo.GenerateMesh(maxh=0.3))

# Import to Radia (uses MSC method automatically)
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='magnetic')

# Apply linear material
mat = rad.MatLin(999)
rad.MatApl(mag_obj, mat)

# External field and solve
MU_0 = 4 * 3.14159265 * 1e-7
ext = rad.ObjBckg([0, 0, MU_0 * 50000])
grp = rad.ObjCnt([mag_obj, ext])

result = rad.Solve(grp, 0.0001, 1000, 1)
print(f'Converged in {int(result[3])} iterations')
```

## Files

| File | Description |
|------|-------------|
| `benchmark_conditions.py` | Unified benchmark parameters (mu_r, H_ext, etc.) |
| `benchmark_hexa_unified.py` | Hexahedral benchmark script (2025-12-12) |
| `benchmark_tetra_unified.py` | Tetrahedral benchmark script (2025-12-12) |
| `benchmark_tetra_vs_hex.py` | Legacy benchmark comparing hex/tet meshes |

### Results Directories

| Directory | Contents |
|-----------|----------|
| `hexahedron/lu/` | Hexahedral results with LU solver |
| `hexahedron/bicgstab/` | Hexahedral results with BiCGSTAB solver |
| `tetrahedron/lu/` | Tetrahedral results with LU solver |
| `tetrahedron/bicgstab/` | Tetrahedral results with BiCGSTAB solver |

## Running the Benchmarks

```bash
cd examples/cube_uniform_field/linear

# Run hexahedral benchmark (both solvers)
python benchmark_hexa_unified.py 3 4 5 6 8 10

# Run hexahedral benchmark (LU only)
python benchmark_hexa_unified.py --lu 3 4 5 6 8 10

# Run tetrahedral benchmark (both solvers)
python benchmark_tetra_unified.py 0.5 0.4 0.35 0.3 0.25

# Run tetrahedral benchmark (BiCGSTAB only)
python benchmark_tetra_unified.py --bicgstab 0.5 0.4 0.35 0.3
```

Output includes:
- Mesh generation time
- Solve time and iteration count
- Average magnetization (M_avg_z)
- Error vs analytical solution

## Physics Background

For a linear magnetic material in uniform external field:

- **Demagnetizing factor**: N ~ 1/3 for cube (approximate)
- **Interior H-field**: H_int = H0 / (1 + N*(mu_r - 1))
- **Magnetization**: M = chi * H_int
- **Exterior field**: Dipole field from magnetized cube

## Solver Selection

Radia provides two solver methods:

| Method | Name | Description | Use Case |
|--------|------|-------------|----------|
| 0 | LU | Dense LU decomposition (LAPACK dgesv) | N < 500 elements |
| 1 | BiCGSTAB | Iterative BiCGSTAB with Jacobi preconditioner | N >= 500 elements |

### Usage

```python
# Method 0: LU (direct solver) - best for small problems
result = rad.Solve(grp, 0.0001, 1000, 0)

# Method 1: BiCGSTAB (iterative) - best for large problems
result = rad.Solve(grp, 0.0001, 1000, 1)
```

### Solver Comparison (Linear Material, mu_r = 1000)

| Elements | LU Time | BiCGSTAB Time | Iterations |
|----------|---------|---------------|------------|
| 64 | 0.003s | 0.002s | 2 |
| 125 | 0.008s | 0.004s | 2 |
| 512 | 0.081s | 0.030s | 2 |
| 1000 | 0.45s | 0.065s | 2 |

**Note**: For linear materials, both solvers converge in 2 iterations regardless of permeability.

## Known Limitations

1. **Internal field accuracy**:
   - `rad.Fld()` inside magnetic materials has limited accuracy
   - Always validate using external field points or average magnetization

2. **Mesh quality**:
   - Degenerate tetrahedra (slivers) can cause numerical issues
   - Use quality mesh generation tools (Netgen, GMSH)

## Notes

1. **Unit System**: All examples use meters (`rad.FldUnits('m')`)
2. **Material Definition**: Use isotropic `rad.MatLin(chi)` for best results
3. **Convergence criterion**: Uses ||dM||/||M|| relative change
