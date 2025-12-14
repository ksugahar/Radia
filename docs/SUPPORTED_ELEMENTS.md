# Supported Magnetic Element Types in Radia

**Version:** 1.3.14
**Date:** 2025-12-15

## Overview

Radia supports the following magnetic element types for magnetostatic field computation:

| Element Type | API | Method | Status |
|--------------|-----|--------|--------|
| **Axis-aligned Rectangular Block** | `rad.ObjRecMag()` + `rad.ObjDivMag()` | Analytical (8-vertex atan) | **Supported** |
| **Hexahedron (6 faces)** | `rad.ObjThckPgn()` or `rad.ObjPolyhdr()` + `HEX_FACES` | MSC | **Supported** |
| **Tetrahedron (4 faces)** | `rad.ObjPolyhdr()` + `TETRA_FACES` | MSC | **Supported** |
| **Wedge/Prism (5 faces)** | `rad.ObjPolyhdr()` + `WEDGE_FACES` | MSC | **Supported** |
| **Pyramid (5 faces)** | `rad.ObjPolyhdr()` + `PYRAMID_FACES` | MSC | **Supported** |
| General polyhedra (>6 faces) | `rad.ObjPolyhdr()` | - | **Not Supported** |

## Supported Elements in Detail

### 1. Axis-Aligned Rectangular Block (`radTRecMag`)

**API:**
```python
import radia as rad
rad.FldUnits('m')

# Create rectangular block
cube = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 1.2e6])

# Subdivide for relaxation
rad.ObjDivMag(cube, [5, 5, 5])  # 125 elements
```

**Characteristics:**
- Uses closed-form analytical formula (8-vertex atan integration)
- **Fastest** computation method
- Faces must be perpendicular to coordinate axes
- Best for structured rectangular grids

### 2. Hexahedron (6-Face MSC)

**Recommended API: ObjThckPgn**
```python
import radia as rad
rad.FldUnits('m')

# Create hexahedron using thick polygon (extruded rectangle)
polygon_2d = [[-0.5, -0.5], [0.5, -0.5], [0.5, 0.5], [-0.5, 0.5]]
hex_elem = rad.ObjThckPgn(-0.5, 1.0, polygon_2d, 'z', [0, 0, 1.2e6])
```

**Alternative API: ObjPolyhdr**
```python
from netgen_mesh_import import HEX_FACES

vertices = [[0,0,0], [1,0,0], [1,1,0], [0,1,0],
            [0,0,1], [1,0,1], [1,1,1], [0,1,1]]
hex_elem = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 1.2e6])
```

**Face Topology (`HEX_FACES`):**
```python
HEX_FACES = [
    [1, 4, 3, 2],  # Bottom (Z-)
    [5, 6, 7, 8],  # Top (Z+)
    [1, 2, 6, 5],  # Front (Y-)
    [3, 4, 8, 7],  # Back (Y+)
    [1, 5, 8, 4],  # Left (X-)
    [2, 3, 7, 6],  # Right (X+)
]
```

**Characteristics:**
- Uses MSC method with 6 quadrilateral faces
- Each quad split into 2 triangles internally
- ELF_MAGIC compatible format
- Suitable for structured/deformed hexahedral meshes

### 3. Tetrahedron (4-Face MSC)

**API:**
```python
import radia as rad
from netgen_mesh_import import TETRA_FACES, netgen_mesh_to_radia

rad.FldUnits('m')

# Direct creation
vertices = [[0,0,0], [1,0,0], [0.5,1,0], [0.5,0.5,1]]
tet = rad.ObjPolyhdr(vertices, TETRA_FACES, [0, 0, 1.2e6])

# Or import from Netgen mesh
mag_obj = netgen_mesh_to_radia(mesh, material={'magnetization': [0, 0, 0]}, units='m')
```

**Face Topology (`TETRA_FACES`):**
```python
TETRA_FACES = [
    [1, 3, 2],  # Bottom (outward normal down)
    [1, 2, 4],  # Front
    [2, 3, 4],  # Right
    [3, 1, 4],  # Left
]
```

**Characteristics:**
- Uses MSC (Magnetic Surface Charge) method
- 4 triangular faces per element
- Suitable for complex curved geometries
- Compatible with Netgen/GMSH mesh generators

### 4. Wedge/Prism (5-Face MSC)

**API:**
```python
import radia as rad
from netgen_mesh_import import WEDGE_FACES

rad.FldUnits('m')

# Wedge vertices: bottom triangle (V1-V3) + top triangle (V4-V6)
vertices = [[0,0,0], [1,0,0], [0.5,1,0],   # Bottom triangle
            [0,0,1], [1,0,1], [0.5,1,1]]   # Top triangle
wedge = rad.ObjPolyhdr(vertices, WEDGE_FACES, [0, 0, 1.2e6])
```

**Face Topology (`WEDGE_FACES`):**
```python
WEDGE_FACES = [
    [1, 3, 2],     # Bottom triangle
    [4, 5, 6],     # Top triangle
    [1, 2, 5, 4],  # Quad face 1
    [2, 3, 6, 5],  # Quad face 2
    [3, 1, 4, 6],  # Quad face 3
]
```

### 5. Pyramid (5-Face MSC)

**API:**
```python
import radia as rad
from netgen_mesh_import import PYRAMID_FACES

rad.FldUnits('m')

# Pyramid vertices: base quadrilateral (V1-V4) + apex (V5)
vertices = [[0,0,0], [1,0,0], [1,1,0], [0,1,0],  # Base
            [0.5,0.5,1]]                          # Apex
pyramid = rad.ObjPolyhdr(vertices, PYRAMID_FACES, [0, 0, 1.2e6])
```

**Face Topology (`PYRAMID_FACES`):**
```python
PYRAMID_FACES = [
    [1, 4, 3, 2],  # Base quadrilateral
    [1, 2, 5],     # Triangle 1
    [2, 3, 5],     # Triangle 2
    [3, 4, 5],     # Triangle 3
    [4, 1, 5],     # Triangle 4
]
```

---

## Mesh Import Utilities

### Netgen Mesh Import

```python
from netgen_mesh_import import netgen_mesh_to_radia
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh

# Create geometry and mesh
geo = OCCGeometry(Box(Pnt(-0.5,-0.5,-0.5), Pnt(0.5,0.5,0.5)))
mesh = Mesh(geo.GenerateMesh(maxh=0.3))

# Import to Radia (tetrahedral elements)
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m')
```

### Nastran Mesh Import

```python
from nastran_mesh_import import create_radia_from_nastran

# Import Nastran .bdf file
mag_obj = create_radia_from_nastran('model.bdf',
                                     material={'magnetization': [0, 0, 1.2e6]},
                                     units='m')
```

**Supported Nastran Elements:**
- `CTETRA` (4-node tetrahedron)
- `CHEXA` (8-node hexahedron)
- `CPENTA` (6-node wedge/prism)
- `CPYRAM` (5-node pyramid)
- `CTRIA3` (surface triangles - grouped by material ID)

---

## Solver Compatibility

All supported element types work with both solvers:

| Solver | Method 0 (LU) | Method 1 (BiCGSTAB) |
|--------|---------------|---------------------|
| Rectangular block | Yes | Yes |
| Tetrahedron MSC | Yes | Yes |
| Hexahedron MSC | Yes | Yes |
| Wedge MSC | Yes | Yes |
| Pyramid MSC | Yes | Yes |

**Recommended solver selection:**
- N < 1000 elements: Either solver (LU is more robust)
- N >= 1000 elements: BiCGSTAB (faster)

---

## Performance Comparison

| Element Type | Field Computation | Memory per Element |
|--------------|-------------------|-------------------|
| Rectangular (analytical) | Fastest | 3 DOF (Mx, My, Mz) |
| Tetrahedron MSC | Medium | 3 DOF |
| Hexahedron MSC | Medium | 3 DOF |
| Wedge/Prism MSC | Medium | 3 DOF |
| Pyramid MSC | Medium | 3 DOF |

---

**Last Updated:** 2025-12-15
**Project:** Radia Magnetic Field Computation
