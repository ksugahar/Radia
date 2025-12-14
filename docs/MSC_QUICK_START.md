# MSC (Magnetic Surface Charge) Quick Start Guide

**Version:** 1.3.14
**Date:** 2025-12-15

## Overview

Radia supports mesh-based magnetic field computation using the **MSC (Magnetic Surface Charge)** method.
Each element has **3 DOF (degrees of freedom)**: Mx, My, Mz (magnetization components).

### Supported Element Types

| Element Type | Faces | API | Use Case |
|--------------|-------|-----|----------|
| **Hexahedron** | 6 quad | `ObjThckPgn()` or `ObjPolyhdr()` + `HEX_FACES` | Structured grids |
| **Tetrahedron** | 4 tri | `ObjPolyhdr()` + `TETRA_FACES` | Complex curved geometry |
| **Wedge/Prism** | 5 (3 quad + 2 tri) | `ObjPolyhdr()` + `WEDGE_FACES` | Hybrid meshes |
| **Pyramid** | 5 (4 tri + 1 quad) | `ObjPolyhdr()` + `PYRAMID_FACES` | Mesh transitions |

---

## Complete Examples

### 1. Hexahedral Mesh with ObjThckPgn (Recommended for Hexahedra)

`ObjThckPgn` creates a thick polygon (extruded 2D polygon) - the simplest way to create hexahedral elements.

```python
import sys
import os
sys.path.insert(0, 'path/to/build/Release')
import radia as rad
import numpy as np

rad.FldUnits('m')
rad.UtiDelAll()

# Physical constants
MU_0 = 4 * np.pi * 1e-7

# Problem setup
CUBE_SIZE = 1.0       # 1 m cube
H_EXT = 50000.0       # 50,000 A/m external field
n_div = 5             # 5x5x5 = 125 elements

# Create hexahedral mesh using ObjThckPgn
elements = []
elem_size = CUBE_SIZE / n_div
half_size = elem_size / 2.0

for ix in range(n_div):
    for iy in range(n_div):
        for iz in range(n_div):
            # Center of this sub-cube
            cx = (ix + 0.5) * elem_size - CUBE_SIZE / 2
            cy = (iy + 0.5) * elem_size - CUBE_SIZE / 2
            cz = (iz + 0.5) * elem_size - CUBE_SIZE / 2

            # Create 2D polygon (XY plane)
            polygon_vertices = [
                [cx - half_size, cy - half_size],
                [cx + half_size, cy - half_size],
                [cx + half_size, cy + half_size],
                [cx - half_size, cy + half_size],
            ]

            # ObjThckPgn(z_base, thickness, vertices_2d, axis, magnetization)
            z_base = cz - half_size
            obj = rad.ObjThckPgn(z_base, elem_size, polygon_vertices, 'z', [0, 0, 0])
            elements.append(obj)

print(f'Created {len(elements)} hexahedral elements')

# Create container and apply material
container = rad.ObjCnt(elements)
mat = rad.MatLin(999)  # mu_r = 1000
rad.MatApl(container, mat)

# External field
ext = rad.ObjBckg([0, 0, MU_0 * H_EXT])
grp = rad.ObjCnt([container, ext])

# Solve (Method 1 = BiCGSTAB)
result = rad.Solve(grp, 0.001, 1000, 1)
print(f'Converged in {int(result[3])} iterations')

# Get average magnetization
all_M = rad.ObjM(container)
M_list = [m[1] for m in all_M]
M_avg_z = np.mean([m[2] for m in M_list])
print(f'M_avg_z = {M_avg_z:.0f} A/m')
```

### 2. Tetrahedral Mesh with Netgen

```python
import sys
import os
sys.path.insert(0, 'path/to/build/Release')
sys.path.append('path/to/src/radia')

import radia as rad
import numpy as np

# IMPORTANT: Import ngsolve BEFORE radia modules
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia

rad.FldUnits('m')
rad.UtiDelAll()

MU_0 = 4 * np.pi * 1e-7
H_EXT = 50000.0

# Create geometry with Netgen
cube_solid = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
cube_solid.mat('magnetic')
geo = OCCGeometry(cube_solid)

# Generate tetrahedral mesh
ngmesh = geo.GenerateMesh(maxh=0.3)
mesh = Mesh(ngmesh)
print(f'Generated {mesh.ne} tetrahedral elements')

# Import to Radia (uses ObjPolyhdr + TETRA_FACES internally)
cube = netgen_mesh_to_radia(mesh,
                             material={'magnetization': [0, 0, 0]},
                             units='m',
                             material_filter='magnetic')

# Apply material and solve
mat = rad.MatLin(999)
rad.MatApl(cube, mat)

ext = rad.ObjBckg([0, 0, MU_0 * H_EXT])
grp = rad.ObjCnt([cube, ext])

result = rad.Solve(grp, 0.001, 1000, 1)
print(f'Converged in {int(result[3])} iterations')
```

### 3. Nonlinear Material with B-H Curve

For saturable magnetic materials, use `MatSatIsoTab` with an H-M table.

```python
import sys
import os
sys.path.insert(0, 'path/to/build/Release')
import radia as rad
import numpy as np

rad.FldUnits('m')
rad.UtiDelAll()

MU_0 = 4 * np.pi * 1e-7
H_EXT = 50000.0

# B-H curve data: [H (A/m), B (T)]
BH_DATA = [
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

# IMPORTANT: Convert B-H to H-M format for Radia
# Formula: M = B/mu_0 - H
HM_DATA = [[h, b/MU_0 - h] for h, b in BH_DATA]

# Create hexahedral mesh
cube = rad.ObjRecMag([0, 0, 0], [1.0, 1.0, 1.0], [0, 0, 0])
rad.ObjDivMag(cube, [10, 10, 10])  # 1000 elements

# Apply nonlinear material
mat = rad.MatSatIsoTab(HM_DATA)
rad.MatApl(cube, mat)

# External field and solve
ext = rad.ObjBckg([0, 0, MU_0 * H_EXT])
grp = rad.ObjCnt([cube, ext])

# Nonlinear solve (typically 3-5 iterations)
result = rad.Solve(grp, 0.001, 100, 1)
print(f'Converged in {int(result[3])} nonlinear iterations')

# Get magnetization
all_M = rad.ObjM(cube)
M_list = [m[1] for m in all_M]
M_avg_z = np.mean([m[2] for m in M_list])
print(f'M_avg_z = {M_avg_z:.0f} A/m')
```

### 4. Manual ObjPolyhdr with HEX_FACES

For custom hexahedral meshes (rotated, deformed), use `ObjPolyhdr` with `HEX_FACES`.

```python
import sys
import os
sys.path.insert(0, 'path/to/build/Release')
sys.path.append('path/to/src/radia')

import radia as rad
from netgen_mesh_import import HEX_FACES

rad.FldUnits('m')
rad.UtiDelAll()

# Define hexahedron vertices (can be arbitrary convex hexahedron)
# Vertex ordering: bottom face CCW (V1-V4), top face CCW (V5-V8)
vertices = [
    [0.0, 0.0, 0.0],  # V1 (index 0)
    [1.0, 0.0, 0.0],  # V2 (index 1)
    [1.0, 1.0, 0.0],  # V3 (index 2)
    [0.0, 1.0, 0.0],  # V4 (index 3)
    [0.0, 0.0, 1.0],  # V5 (index 4)
    [1.0, 0.0, 1.0],  # V6 (index 5)
    [1.0, 1.0, 1.0],  # V7 (index 6)
    [0.0, 1.0, 1.0],  # V8 (index 7)
]

# HEX_FACES defines 6 quadrilateral faces (1-indexed)
# Each face: [v1, v2, v3, v4] with outward normal (CCW when viewed from outside)
print(f'HEX_FACES = {HEX_FACES}')

# Create hexahedral element
hex_elem = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 1.0e6])

# Compute field at a point
B = rad.Fld(hex_elem, 'b', [0.5, 0.5, 1.5])
print(f'B at (0.5, 0.5, 1.5) = {B} T')
```

---

## API Reference

### ObjThckPgn - Thick Polygon (Extruded 2D Polygon)

```python
obj = rad.ObjThckPgn(z_base, thickness, vertices_2d, axis, magnetization)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `z_base` | float | Base position along extrusion axis |
| `thickness` | float | Extrusion length |
| `vertices_2d` | list[[x,y], ...] | 2D polygon vertices (CCW order) |
| `axis` | str | Extrusion axis: 'x', 'y', or 'z' |
| `magnetization` | [Mx, My, Mz] | Initial magnetization vector |

### ObjPolyhdr - General Polyhedron

```python
obj = rad.ObjPolyhdr(vertices, faces, magnetization)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `vertices` | list[[x,y,z], ...] | 3D vertex coordinates |
| `faces` | list[[v1,v2,...], ...] | Face vertex indices (1-indexed!) |
| `magnetization` | [Mx, My, Mz] | Initial magnetization vector |

**IMPORTANT**: Face vertex indices are **1-indexed** (Radia convention).

### Face Topology Constants

Import from `netgen_mesh_import`:

```python
from netgen_mesh_import import TETRA_FACES, HEX_FACES, WEDGE_FACES, PYRAMID_FACES
```

| Constant | Faces | Description |
|----------|-------|-------------|
| `TETRA_FACES` | 4 | Tetrahedron: 4 triangular faces |
| `HEX_FACES` | 6 | Hexahedron: 6 quadrilateral faces |
| `WEDGE_FACES` | 5 | Wedge/Prism: 2 tri + 3 quad |
| `PYRAMID_FACES` | 5 | Pyramid: 1 quad + 4 tri |

### MatSatIsoTab - Nonlinear Isotropic Material

```python
mat = rad.MatSatIsoTab(HM_data)
rad.MatApl(obj, mat)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `HM_data` | list[[H, M], ...] | H-M curve data points |

**B-H to H-M Conversion**:
```python
# B = mu_0 * (H + M) -> M = B/mu_0 - H
HM_DATA = [[h, b/MU_0 - h] for h, b in BH_DATA]
```

### MatLin - Linear Isotropic Material

```python
mat = rad.MatLin(chi)  # chi = mu_r - 1
rad.MatApl(obj, mat)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `chi` | float | Magnetic susceptibility (mu_r - 1) |

---

## Solver Selection Guide

| Problem Size | Element Count | Recommended Solver | Code |
|--------------|---------------|-------------------|------|
| Small | < 1,000 | LU (direct) | `rad.Solve(grp, tol, max_iter, 0)` |
| Medium | 1,000 - 10,000 | BiCGSTAB | `rad.Solve(grp, tol, max_iter, 1)` |
| Large | > 10,000 | BiCGSTAB | `rad.Solve(grp, tol, max_iter, 1)` |

**Typical Parameters**:
- `tol`: 0.001 (0.1% relative tolerance)
- `max_iter`: 100-1000

**Iteration Counts**:
- Linear materials: 1-2 iterations
- Nonlinear materials: 3-6 iterations

---

## Performance Comparison (ELF_MAGIC vs Radia)

Benchmark: 1m cube, 50,000 A/m external field, nonlinear BH curve.

### LU Solver

| N | Elements | Radia Time | ELF Time | Speedup |
|---|----------|------------|----------|---------|
| 10 | 1,000 | 1.95s | 3.73s | **1.9x** |
| 15 | 3,375 | 24.35s | 61.65s | **2.5x** |
| 20 | 8,000 | ~280s | 595.23s | **2.1x** |

### BiCGSTAB Solver

| N | Elements | Radia Time | ELF Time | Speedup |
|---|----------|------------|----------|---------|
| 10 | 1,000 | 0.55s | 4.05s | **7.4x** |
| 15 | 3,375 | 7.30s | 54.49s | **7.5x** |
| 20 | 8,000 | 51.81s | 343.01s | **6.6x** |

**Key Finding**: Radia BiCGSTAB is **6-7x faster** than ELF_MAGIC for the same mesh and material.

---

## Common Issues and Solutions

### 1. "DLL load failed" Error

**Problem**: ImportError when importing radia_ngsolve

**Solution**: Import ngsolve BEFORE radia_ngsolve:
```python
import radia as rad
import ngsolve  # Import BEFORE radia_ngsolve
from radia import radia_ngsolve
```

### 2. Coordinates Off by 1000x

**Problem**: Results seem completely wrong

**Solution**: Set units to meters:
```python
rad.FldUnits('m')  # REQUIRED for NGSolve/Netgen integration
```

### 3. Solver Not Converging

**Problem**: Solver doesn't converge for nonlinear materials

**Solution**:
1. Use BiCGSTAB (Method 1) instead of LU (Method 0)
2. Increase max iterations
3. Check B-H curve data (must be monotonically increasing)
4. Ensure H-M conversion is correct: `M = B/mu_0 - H`

### 4. ObjPolyhdr Face Indexing Error

**Problem**: "Invalid face indices" error

**Solution**: Face indices are **1-indexed** in Radia:
```python
# CORRECT (1-indexed)
faces = [[1, 2, 3], [1, 3, 4], ...]

# WRONG (0-indexed)
faces = [[0, 1, 2], [0, 2, 3], ...]
```

---

## Face Topology Diagrams

### Tetrahedron (TETRA_FACES)

```
        V4 (apex)
       /|\
      / | \
     /  |  \
    /   |   \
   V1---+---V3
    \   |   /
     \  |  /
      \ | /
       \|/
        V2

TETRA_FACES = [
    [1, 3, 2],  # Bottom: V1-V3-V2
    [1, 2, 4],  # Front:  V1-V2-V4
    [2, 3, 4],  # Right:  V2-V3-V4
    [3, 1, 4],  # Left:   V3-V1-V4
]
```

### Hexahedron (HEX_FACES)

```
        V8--------V7
       /|        /|
      / |       / |
     V5--------V6 |
     |  V4-----|--V3
     | /       | /
     |/        |/
     V1--------V2

HEX_FACES = [
    [1, 4, 3, 2],  # Bottom (Z-): V1-V4-V3-V2
    [5, 6, 7, 8],  # Top (Z+):    V5-V6-V7-V8
    [1, 2, 6, 5],  # Front (Y-):  V1-V2-V6-V5
    [3, 4, 8, 7],  # Back (Y+):   V3-V4-V8-V7
    [1, 5, 8, 4],  # Left (X-):   V1-V5-V8-V4
    [2, 3, 7, 6],  # Right (X+):  V2-V3-V7-V6
]
```

---

## Related Documentation

- [examples/cube_uniform_field/nonlinear/](../examples/cube_uniform_field/nonlinear/) - Benchmark examples
- [examples/cube_uniform_field/linear/](../examples/cube_uniform_field/linear/) - Linear material examples
- [src/radia/netgen_mesh_import.py](../src/radia/netgen_mesh_import.py) - Mesh import utilities

---

**Last Updated:** 2025-12-15
**Project:** Radia Magnetic Field Computation
