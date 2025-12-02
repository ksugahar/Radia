# Linear Material Benchmark: Cube in Uniform External Field

This benchmark validates Radia's tetrahedral and hexahedral mesh solvers for **linear magnetic materials**.

## Problem Setup

- **Geometry**: 100mm x 100mm x 100mm cube centered at origin
- **Material**: Linear isotropic, mu_r = 100 (chi = 99)
- **External Field**: B0 = 1.0 T (uniform along z-axis)
- **Validation Metric**: External field at points OUTSIDE the cube

**Important**: Radia MMM is designed for external field calculation. The proper validation metric is the field in air regions, not inside the magnetic material.

## Benchmark Results (2025-12-02)

### External Field Comparison

Reference: Hexahedral 8x8x8 mesh (512 elements)

| Mesh Type | Elements | Avg Error | Max Error | Notes |
|-----------|----------|-----------|-----------|-------|
| **Hexahedral Meshes** |
| Hex 1x1x1 | 1 | 3.33% | 9.27% | Single element |
| Hex 2x2x2 | 8 | 2.38% | 6.80% | |
| Hex 3x3x3 | 27 | 1.17% | 3.53% | |
| Hex 4x4x4 | 64 | 0.65% | 1.98% | |
| Hex 5x5x5 | 125 | 0.37% | 1.12% | |
| Hex 6x6x6 | 216 | 0.19% | 0.59% | Good convergence |
| **Tetrahedral Meshes (Manual 5-tetra cells)** |
| Tet 1^3x5 | 5 | 3.03% | 8.18% | Single cell |
| Tet 2^3x5 | 40 | 0.43% | 0.87% | Good accuracy |
| Tet 3^3x5 | 135 | 0.23% | 0.62% | Excellent |
| **Tetrahedral Meshes (Netgen/NGSolve)** |
| maxh=0.10 | 12 | 1.93% | 5.34% | |
| maxh=0.06 | 28 | 0.20% | 0.49% | Best accuracy |
| maxh=0.04 | 80 | 0.26% | 0.67% | |

### Key Findings

1. **Hexahedral mesh convergence**: Error decreases consistently with mesh refinement
2. **Tetrahedral mesh accuracy**: Comparable to hexahedral at similar element counts
3. **External field validation**: Proper metric for MMM validation

## Mesh Configurations

### Hexahedral Mesh (ObjRecMag + ObjDivMag)

Standard Radia approach using built-in mesh subdivision:

```python
import radia as rad
rad.FldUnits('m')

cube = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 0])
rad.ObjDivMag(cube, [5, 5, 5])  # 5x5x5 = 125 elements
mat = rad.MatLin(99.0)  # chi = mu_r - 1
rad.MatApl(cube, mat)
```

### Tetrahedral Mesh (Manual 5-tetra decomposition)

Each cubic cell is decomposed into 5 tetrahedra:

```python
rad.SolverTetraMethod(0)  # Method 0: polygon-based
# Create 5 tetrahedra per cubic cell using rad.ObjPolyhdr()
```

### Tetrahedral Mesh (Netgen import)

Using netgen_mesh_import for unstructured meshes:

```python
from netgen_mesh_import import netgen_mesh_to_radia
from netgen.occ import Box, OCCGeometry
from ngsolve import Mesh

cube_solid = Box((-0.05, -0.05, -0.05), (0.05, 0.05, 0.05))
geo = OCCGeometry(cube_solid)
mesh = Mesh(geo.GenerateMesh(maxh=0.06))

mag_obj = netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m')
```

## Files

| File | Description |
|------|-------------|
| `benchmark_external_field.py` | Main benchmark comparing hex/tet meshes |
| `compare_external_field.py` | External field comparison with NGSolve |
| `precision_evaluation.py` | Comprehensive precision evaluation |
| `cube_benchmark_ngsolve.vtu` | NGSolve solution (VTU format) |
| `cube_benchmark_radia.vtk` | Radia solution (VTK format) |

## Running the Benchmark

```bash
cd examples/cube_uniform_field/linear
python benchmark_external_field.py
```

Output will show:
- Hexahedral mesh convergence study (1x1x1 to 6x6x6)
- Tetrahedral mesh (manual 5-tetra decomposition)
- Tetrahedral mesh (Netgen import)
- Summary table with error statistics

## Physics Background

For a linear magnetic material in uniform external field:

- **Demagnetizing factor**: N ~ 1/3 for cube (approximate)
- **Interior H-field**: H_int = H0 / (1 + N*(mu_r - 1))
- **Magnetization**: M = chi * H_int
- **Exterior field**: Dipole field from magnetized cube

## Known Limitations

1. **High permeability**: For mu_r > 100, solver convergence may be slow
2. **Adjacent tetrahedra**: Dense tetrahedral meshes may show numerical issues
3. **Internal field**: `rad.Fld()` inside magnetic materials has limited accuracy

## Notes

1. **Unit System**: All examples use meters (`rad.FldUnits('m')`)
2. **Material Definition**: Use isotropic `rad.MatLin(chi)` for best results
3. **External field**: Evaluate at points outside the magnetic material for validation
