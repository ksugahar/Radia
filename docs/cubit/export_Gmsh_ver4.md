# export_Gmsh_ver4

Export mesh to Gmsh format version 4.1.

## Synopsis

```python
cubit_mesh_export.export_Gmesh(cubit, FileName, version="4.1", DIM="auto")
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface object |
| `FileName` | str | required | Output file path for the .msh file |
| `DIM` | str | `"auto"` | Dimension mode (see below) |

### DIM Parameter Options

| Value | Description |
|-------|-------------|
| `"auto"` | Auto-detect dimension (3D if volume elements exist, else 2D) |
| `"2D"` | 2D mode - orient surface element normals to +z direction, z-coordinates set to 0 |
| `"3D"` | 3D mode - no normal orientation applied |

## Returns

Returns the `cubit` object for method chaining.

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
import cubit_mesh_export

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")

cubit_mesh_export.export_Gmesh(cubit, "mesh.msh")
```

### 2D Export with Normal Orientation

```python
cubit.cmd("create surface rectangle width 1 height 1 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("mesh surface 1")
cubit.cmd("block 1 add tri all")

# Force 2D mode to ensure normals point in +z direction
cubit_mesh_export.export_Gmesh(cubit, "plate.msh", DIM="2D")
```

### 2nd Order Elements

```python
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 element type tetra10")  # Convert to 2nd order

cubit_mesh_export.export_Gmesh(cubit, "mesh_2nd_order.msh")
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
