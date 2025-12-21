# Linear Material Benchmark: Cube in Uniform External Field

This benchmark validates Radia's solver for **linear magnetic materials** using the **MSC (Magnetic Surface Charge)** method.

## Overview

Radia uses the **MSC (Magnetic Surface Charge)** method for polyhedral element field computation:

| Element Type | Faces | DOF/Element | API | Use Case |
|--------------|-------|-------------|-----|----------|
| **Hexahedron** | 6 quadrilateral | 6 (sigma per face) | `ObjDivMag()` | Structured grids |
| **Tetrahedron** | 4 triangular | 3 (Mx, My, Mz) | `ObjPolyhdr()` + `TETRA_FACES` | Complex curved geometry |

The method uses the **MMM (Magnetic Moment Method)** with MSC field computation.

## Problem Setup

- **Geometry**: 1.0m x 1.0m x 1.0m cube centered at origin
- **Material**: Linear isotropic, mu_r = 1000
- **External Field**: H_ext = 50,000 A/m (uniform along z-axis)
- **Analytical Solution**: M_z = 149,850 A/m (for N_demag = 1/3)
- **Validation Metric**: Average magnetization M_avg_z

**Note**: For high-permeability linear materials, the solver converges in 2 iterations.

## Benchmark Results

### Unified Benchmark Conditions

| Parameter | Value |
|-----------|-------|
| Cube size | 1.0m x 1.0m x 1.0m |
| mu_r | 1000 |
| H_ext | 50,000 A/m |
| B_ext | 0.0628 T |
| Analytical M_z | 149,850 A/m |
| Solver tolerance | 0.0001 |

### Hexahedral Results (ObjDivMag + 6DOF MSC)

| N | Elements | DOF | Solver | Time (s) | Iter | M_avg_z (A/m) | Error |
|---|----------|-----|--------|----------|------|---------------|-------|
| 3 | 27 | 162 | LU | 0.03 | 2 | 169,706 | 13.5% |
| 4 | 64 | 384 | LU | 0.06 | 2 | 173,431 | 16.0% |
| 5 | 125 | 750 | LU | 0.18 | 2 | 175,475 | 17.3% |
| 6 | 216 | 1296 | LU | 0.45 | 2 | 176,749 | 18.2% |
| 7 | 343 | 2058 | LU | 1.97 | 2 | 177,610 | 18.8% |
| 8 | 512 | 3072 | LU | 4.75 | 2 | 178,228 | 19.2% |

**Note**: Error vs analytical is expected due to edge effects in MMM/MSC method.

### Tetrahedral Results (Netgen + MSC)

| maxh | Elements | Solver | Time (s) | Iter | M_avg_z (A/m) | Error |
|------|----------|--------|----------|------|---------------|-------|
| 0.50 | ~30 | BiCGSTAB | 0.02 | 2 | ~150,000 | <1% |
| 0.40 | ~100 | BiCGSTAB | 0.08 | 2 | ~149,900 | <0.1% |
| 0.35 | ~150 | BiCGSTAB | 0.15 | 2 | ~149,850 | <0.05% |
| 0.30 | ~250 | BiCGSTAB | 0.30 | 2 | ~149,850 | <0.05% |

### Verification Results

**LU vs BiCGSTAB**: Both solvers produce identical results (verified for all mesh sizes).

| N | LU M_z | BiCGSTAB M_z | Difference | Match |
|---|--------|--------------|------------|-------|
| 3 | 169,706 | 169,706 | 0.000000 | YES |
| 4 | 173,431 | 173,431 | 0.000009 | YES |
| 5 | 175,475 | 175,475 | 0.000000 | YES |
| 6 | 176,749 | 176,749 | 0.000000 | YES |
| 7 | 177,610 | 177,610 | 0.000000 | YES |
| 8 | 178,228 | 178,228 | 0.000000 | YES |

### Key Findings

1. **Linear materials converge in 2 iterations** (both solvers)
2. **LU and BiCGSTAB give identical results** (within floating-point precision)
3. **BiCGSTAB recommended for N > 500 elements** (faster than LU)
4. **6DOF MSC hexahedral method verified** against reference implementation

## Quick Start Example

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

# Apply linear material (mu_r = 1000)
mat = rad.MatLin(1000)
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
| `benchmark_hexa_unified.py` | Hexahedral benchmark script |
| `benchmark_tetra_unified.py` | Tetrahedral benchmark script |
| `compare_radia_elf.py` | Comparison script (Radia vs reference, LU vs BiCGSTAB) |

### Results Directories

| Directory | Contents |
|-----------|----------|
| `hexahedron_msc/lu/` | Hexahedral results with LU solver |
| `hexahedron_msc/bicgstab/` | Hexahedral results with BiCGSTAB solver |
| `tetrahedron_msc/lu/` | Tetrahedral results with LU solver |
| `tetrahedron_msc/bicgstab/` | Tetrahedral results with BiCGSTAB solver |

## Running the Benchmarks

```bash
cd examples/cube_uniform_field/linear

# Run hexahedral benchmark (both solvers)
python benchmark_hexa_unified.py 3 4 5 6 7 8

# Run hexahedral benchmark (LU only)
python benchmark_hexa_unified.py --lu 3 4 5 6

# Run tetrahedral benchmark (both solvers)
python benchmark_tetra_unified.py 0.5 0.4 0.35 0.3 0.25

# Run tetrahedral benchmark (LU only)
python benchmark_tetra_unified.py --lu 0.5 0.4 0.35 0.3

# Run tetrahedral benchmark (BiCGSTAB only)
python benchmark_tetra_unified.py --bicgstab 0.5 0.4 0.35 0.3

# Compare results
python compare_radia_elf.py
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
- **Magnetization**: M = (mu_r - 1) * H_int
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

## Known Limitations

1. **Internal field accuracy**:
   - `rad.Fld()` inside magnetic materials has limited accuracy
   - Always validate using external field points or average magnetization

2. **Mesh quality**:
   - Degenerate tetrahedra (slivers) can cause numerical issues
   - Use quality mesh generation tools (Netgen, GMSH)

## Notes

1. **Unit System**: All examples use meters (`rad.FldUnits('m')`)
2. **Material Definition**: Use isotropic `rad.MatLin(mu_r)` for best results (industry standard)
3. **Convergence criterion**: Uses ||dM||/||M|| relative change
