# export_Gmsh_ver4

Export mesh to Gmsh format version 4.1.

## Cubit Plugin (Recommended)

```
radia_export gmsh "mesh.msh" version 4.1
radia_export gmsh "mesh.msh" version 4.1 dim 2d
```

No block assignment or `#!python` required. See [export_Gmsh_ver2.md](export_Gmsh_ver2.md) for full plugin documentation.

---

## Plugin Command

```python
cubit.cmd('radia_export gmsh "mesh.msh" version 4 overwrite')
cubit.cmd('radia_export gmsh "mesh.msh" version 4 dimension 2 overwrite')
```

> **Note**: The old `cubit_mesh_export.export_Gmesh()` Python function has been replaced by the `radia_export gmsh` plugin command. The old Python module (`src/radia/cubit_mesh_export.py`) has been replaced by the C++ pybind11 module (`src/cubit_plugin/radia_cubit_pybind.cpp`).

### DIM Parameter Options

| Option | Description |
|--------|-------------|
| `dimension 3` (default) | 3D mode - no normal orientation applied |
| `dimension 2` | 2D mode - orient surface element normals to +z direction, z-coordinates set to 0 |

## Supported Elements

### 1st Order Elements

| Element Type | Gmsh Code | Nodes |
|--------------|-----------|-------|
| Point | 15 | 1 |
| Line | 1 | 2 |
| Triangle | 2 | 3 |
| Quadrilateral | 3 | 4 |
| Tetrahedron | 4 | 4 |
| Hexahedron | 5 | 8 |
| Wedge/Prism | 6 | 6 |
| Pyramid | 7 | 5 |

### 2nd Order Elements

| Element Type | Gmsh Code | Nodes |
|--------------|-----------|-------|
| Line3 | 8 | 3 |
| Triangle6 | 9 | 6 |
| Triangle7 | 42 | 7 |
| Quad8 | 16 | 8 |
| Quad9 | 10 | 9 |
| Tetrahedron10 | 11 | 10 |
| Tetrahedron11 | 35 | 11 |
| Hexahedron20 | 17 | 20 |
| Wedge15 | 18 | 15 |
| Pyramid13 | 19 | 13 |

## File Format

The Gmsh v4.1 format includes:

1. **$MeshFormat** - Version identifier (4.1 0 8)
2. **$PhysicalNames** - Block names with dimensions
3. **$Entities** - Geometry topology (vertices, curves, surfaces, volumes)
4. **$Nodes** - Node coordinates grouped by entity
5. **$Elements** - Elements grouped by entity block

## Usage Examples

### Basic 3D Export

```python
import cubit

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")

cubit.cmd('radia_export gmsh "mesh.msh" version 4 overwrite')
```

### 2D Export with Normal Orientation

```python
cubit.cmd("create surface rectangle width 1 height 1 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("mesh surface 1")

# Force 2D mode to ensure normals point in +z direction
cubit.cmd('radia_export gmsh "plate.msh" version 4 dimension 2 overwrite')
```

### 2nd Order Elements

```python
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 element type tetra10")  # Convert to 2nd order

cubit.cmd('radia_export gmsh "mesh_2nd_order.msh" version 4 overwrite')
```

## Differences from v2

| Feature | v2 (`export_Gmsh_ver2`) | v4 (`export_Gmsh_ver4`) |
|---------|------------------------|------------------------|
| Format version | 2.2 | 4.1 |
| $Entities section | No | Yes |
| Node grouping | Flat list | Entity blocks |
| Element grouping | Flat list | Entity blocks |
| DIM parameter | No | Yes |
| 2D normal control | No | Yes |

## $Entities Section

The v4 format includes geometry topology information:

- **Points**: Vertex coordinates
- **Curves**: Bounding box and vertex references
- **Surfaces**: Bounding box, physical tags, curve references
- **Volumes**: Bounding box, physical tags, surface references

This enables solvers to associate mesh elements with geometric entities.

## Compatibility

The output file is compatible with:
- Gmsh 4.x
- NGSolve/Netgen
- FEniCS
- Other solvers supporting Gmsh v4 format

## See Also

- [export_Gmsh_ver2](export_Gmsh_ver2.md) - Gmsh v2.2 format export
- [Cubit_Element_Order.md](Cubit_Element_Order.md) - How to control element order
