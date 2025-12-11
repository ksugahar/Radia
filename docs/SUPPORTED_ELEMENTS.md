# Supported Magnetic Element Types in Radia

**Version:** 1.3.13
**Date:** 2025-12-11

## Overview

Radia currently supports the following magnetic element types for magnetostatic field computation:

| Element Type | API | Method | Status |
|--------------|-----|--------|--------|
| **Axis-aligned Rectangular Block** | `rad.ObjRecMag()` + `rad.ObjDivMag()` | Analytical (8-vertex atan) | **Supported** |
| **Tetrahedron (4 faces)** | `rad.ObjPolyhdr()` with `TETRA_FACES` | MSC (Magnetic Surface Charge) | **Supported** |
| **Hexahedron (6 faces)** | `rad.ObjPolyhdr()` with `HEX_FACES` | MSC (6-quad -> 12-tri) | **Supported** |
| General polyhedra (>6 faces) | `rad.ObjPolyhdr()` | - | **Not Supported** |
| Wedge/Prism | - | - | **Not Supported** |
| Pyramid | - | - | **Not Supported** |

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

**Source:** [rad_rectangular_block.cpp](../src/core/rad_rectangular_block.cpp)

### 2. Tetrahedron (4-Face MSC)

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
    [1, 3, 2],  # Face 0: bottom (outward normal down)
    [1, 2, 4],  # Face 1: front
    [2, 3, 4],  # Face 2: right
    [3, 1, 4],  # Face 3: left
]
```

**Characteristics:**
- Uses MSC (Magnetic Surface Charge) method
- 4 triangular faces per element
- Suitable for complex curved geometries
- Compatible with Netgen/GMSH mesh generators

**Source:** [rad_polyhedron.cpp](../src/core/rad_polyhedron.cpp) - `B_comp_tetrahedron_MSC()`

### 3. Hexahedron (6-Face MSC)

**API:**
```python
import radia as rad
from netgen_mesh_import import HEX_FACES

rad.FldUnits('m')

# Direct creation
vertices = [[0,0,0], [1,0,0], [1,1,0], [0,1,0],
            [0,0,1], [1,0,1], [1,1,1], [0,1,1]]
hex_elem = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 1.2e6])
```

**Face Topology (`HEX_FACES`):**
```python
HEX_FACES = [
    [1, 4, 3, 2],  # Face 0: bottom (Z-)
    [5, 6, 7, 8],  # Face 1: top (Z+)
    [1, 2, 6, 5],  # Face 2: front (Y-)
    [3, 4, 8, 7],  # Face 3: back (Y+)
    [1, 5, 8, 4],  # Face 4: left (X-)
    [2, 3, 7, 6],  # Face 5: right (X+)
]
```

**Characteristics:**
- Uses MSC method with 6 quadrilateral faces
- Each quad split into 2 triangles: [V0,V1,V2] + [V0,V2,V3]
- ELF_MAGIC compatible format
- Suitable for structured/deformed hexahedral meshes

**Source:** [rad_polyhedron.cpp](../src/core/rad_polyhedron.cpp) - `B_comp_hexahedron_MSC()`

## Unsupported Elements

The following element types are **NOT currently supported**:

| Element Type | Reason |
|--------------|--------|
| **General polyhedra (>6 faces)** | MSC implementation limited to tet/hex |
| **Wedge/Prism (5 faces)** | Not implemented |
| **Pyramid (5 faces)** | Not implemented |
| **12-face triangular hexahedron** | Deprecated, use 6-face hex instead |

**If you need unsupported element types:**
1. Decompose into tetrahedra (always possible)
2. Use Netgen for automatic tetrahedral meshing
3. Contact developers for feature requests

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

# Import Nastran .bdf file (supports CTETRA, CHEXA)
mag_obj = create_radia_from_nastran('model.bdf',
                                     material={'magnetization': [0, 0, 1.2e6]},
                                     units='m')
```

**Supported Nastran Elements:**
- `CTETRA` (4-node tetrahedron) -> Radia tetrahedron MSC
- `CHEXA` (8-node hexahedron) -> Radia hexahedron MSC
- `CTRIA3` (surface triangles) -> Grouped by material ID

## Solver Compatibility

All supported element types work with both solvers:

| Solver | Method 0 (LU) | Method 1 (BiCGSTAB) |
|--------|---------------|---------------------|
| Rectangular block | Yes | Yes |
| Tetrahedron MSC | Yes | Yes |
| Hexahedron MSC | Yes | Yes |

**Recommended solver selection:**
- N < 500 elements: Either solver (LU is more robust)
- N >= 500 elements: BiCGSTAB (faster)

## Performance Comparison

| Element Type | Field Computation | Memory per Element |
|--------------|-------------------|-------------------|
| Rectangular (analytical) | Fastest | 3 DOF (Mx, My, Mz) |
| Tetrahedron MSC | Medium | 3 DOF |
| Hexahedron MSC | Medium | 3 DOF |

**Note:** Hexahedron MSC is slightly slower than rectangular analytical due to surface integral computation, but allows arbitrary orientation.

## Future Development

The following features are under consideration:

1. **HACApK H-matrix acceleration** - Currently under research
2. **Wedge/Prism elements** - Low priority
3. **Higher-order elements** - Not planned

---

**Last Updated:** 2025-12-11
**Project:** Radia Magnetic Field Computation
