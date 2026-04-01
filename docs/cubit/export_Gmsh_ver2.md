# export_Gmsh_ver2

Export mesh to Gmsh format version 2.2.

## Cubit Plugin (Recommended)

The Radia Cubit plugin provides a native APREPRO command:

```
radia export gmsh "mesh.msh"
radia export gmsh "mesh.msh" version 4
radia export gmsh "mesh.msh" version 4 dimension 2
```

### Syntax

```
radia export gmsh <"filename"> [version <2|4>] [dimension <2|3>] [overwrite]
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `version 2` | yes | Gmsh v2.2 format |
| `version 4` | | Gmsh v4.1 format |
| `dimension 3` | yes | 3D mode |
| `dimension 2` | | 2D mode (orient surface normals to +z) |
| `overwrite` | off | Overwrite existing file |

**Advantages:** No block assignment, no `#!python`, pure APREPRO command.

**Installation:** Copy `radia_cubit.ccm` to `<Cubit install>/bin/plugins/`.

---

## Python API (via Plugin Command)

```python
cubit.cmd('radia export gmsh "mesh.msh" overwrite')
```

> **Note**: The old `cubit_mesh_export.export_Gmesh()` Python function has been replaced by the `radia export gmsh` plugin command. The old Python module (`src/radia/cubit_mesh_export.py`) has been replaced by the C++ pybind11 module (`src/cubit_plugin/radia_cubit_pybind.cpp`).

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

The Gmsh v2.2 format uses a simple flat structure:

1. **$MeshFormat** - Version identifier (2.2 0 8)
2. **$PhysicalNames** - Block names with dimensions
3. **$Nodes** - Node coordinates as a flat list
4. **$Elements** - Elements as a flat list

This format is simpler than v4 but lacks geometry topology information.

## Usage Examples

### Basic 3D Export

```python
import cubit

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")

cubit.cmd('radia export gmsh "mesh.msh" overwrite')
```

### 2nd Order Elements

```python
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 element type tetra10")  # Convert to 2nd order

cubit.cmd('radia export gmsh "mesh_2nd_order.msh" overwrite')
```

### Mixed Element Types

```python
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")

cubit.cmd('radia export gmsh "mixed.msh" overwrite')
```

## Comparison with v4

| Feature | v2 (`export_Gmsh_ver2`) | v4 (`export_Gmsh_ver4`) |
|---------|------------------------|------------------------|
| Format version | 2.2 | 4.1 |
| $Entities section | No | Yes |
| Node grouping | Flat list | Entity blocks |
| Element grouping | Flat list | Entity blocks |
| DIM parameter | No | Yes |
| 2D normal control | No | Yes |

**When to use v2:**
- Maximum compatibility with older software
- Simple mesh transfer without geometry information
- Smaller file size for large meshes

**When to use v4:**
- Need geometry topology ($Entities section)
- 2D meshes requiring normal orientation control
- Modern solvers (NGSolve, FEniCS)

## Compatibility

The output file is compatible with:
- Gmsh 2.x and later
- NGSolve/Netgen
- FEniCS
- Most FEM solvers supporting Gmsh format

## See Also

- [export_Gmsh_ver4](export_Gmsh_ver4.md) - Gmsh v4.1 format with $Entities section
- [Cubit_Element_Order.md](Cubit_Element_Order.md) - How to control element order
