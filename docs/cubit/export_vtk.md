# export_vtk

Export mesh to VTK Legacy format (.vtk) for visualization in ParaView etc.

## Syntax

```
radia_export vtk "filename.vtk" [order <1|2>] [dimension <2|3>] [overwrite]
```

No block assignment required. Sidesets are exported as surface cells.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `order 1` | yes | 1st order elements |
| `order 2` | | 2nd order (quadratic) elements via NetgenCurver |
| `dimension 3` | yes | 3D solid mesh |
| `dimension 2` | | 2D shell mesh (z forced to 0) |
| `overwrite` | off | Overwrite existing file |

## Supported Elements

### 1st Order

| Element | VTK Type | Code |
|---------|----------|------|
| Point | VTK_VERTEX | 1 |
| Line | VTK_LINE | 3 |
| Triangle | VTK_TRIANGLE | 5 |
| Quad | VTK_QUAD | 9 |
| Tetrahedron | VTK_TETRA | 10 |
| Hexahedron | VTK_HEXAHEDRON | 12 |
| Wedge | VTK_WEDGE | 13 |
| Pyramid | VTK_PYRAMID | 14 |

### 2nd Order

| Element | VTK Type | Code | Nodes |
|---------|----------|------|-------|
| Triangle6 | VTK_QUADRATIC_TRIANGLE | 22 | 6 |
| Quad8 | VTK_QUADRATIC_QUAD | 23 | 8 |
| Tetrahedron10 | VTK_QUADRATIC_TETRA | 24 | 10 |
| Hexahedron20 | VTK_QUADRATIC_HEXAHEDRON | 25 | 20 |
| Wedge15 | VTK_QUADRATIC_WEDGE | 26 | 15 |
| Pyramid13 | VTK_QUADRATIC_PYRAMID | 27 | 13 |

## File Format

VTK Legacy ASCII format (version 3.0), `UNSTRUCTURED_GRID` dataset:

```
# vtk DataFile Version 3.0
Radia Cubit Plugin (order 2)
ASCII
DATASET UNSTRUCTURED_GRID
POINTS <n> double
...
CELLS <n> <size>
...
CELL_TYPES <n>
...
CELL_DATA <n>
SCALARS GroupID int 1
LOOKUP_TABLE default
...
SCALARS GroupType int 1
LOOKUP_TABLE default
...
```

### Cell Data Fields

| Field | Description |
|-------|-------------|
| `GroupID` | Cubit block ID (volume elements) or sideset ID (surface elements) |
| `GroupType` | 0 = block element, 1 = sideset face element |

## Usage Examples

### Basic Export

```python
import cubit

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")

cubit.cmd('radia_export vtk "mesh.vtk" overwrite')
```

### 2nd Order with Curved Surfaces

```python
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.3")
cubit.cmd("mesh volume 1")

cubit.cmd('radia_export vtk "sphere_o2.vtk" order 2 overwrite')
```

## Limitations

- Order 3+ is not supported (VTK Legacy format has no standard cell types beyond quadratic)
- VTK_QUADRATIC_PYRAMID (type 27) cannot be imported by Gmsh

## See Also

- [export_Gmsh.md](export_Gmsh.md) — Gmsh v4.1 export (order 1-3)
- [Function_Reference.md](Function_Reference.md) — All plugin commands
