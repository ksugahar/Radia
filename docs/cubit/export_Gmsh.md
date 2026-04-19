# export_Gmsh

Export mesh to Gmsh v4.1 format with high-order element support (order 1-3).

## Syntax

```
radia_export gmsh "filename.msh" [order <1-3>] [dimension <2|3>] [overwrite]
```

No block assignment required — all meshed elements are exported automatically.
Sidesets are exported as surface elements. Nodesets as point elements.

> **IMPORTANT**: Use `radia_export gmsh`, NOT `export gmsh`.
> Cubit has no built-in `export gmsh` command — only `radia_export gmsh` is available.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `order 1` | yes | 1st order elements |
| `order 2` | | 2nd order (edge mid-nodes via NetgenCurver + ACIS geometry projection) |
| `order 3` | | 3rd order (edge + face nodes via NetgenCurver; wedge not supported) |
| `dimension 3` | yes | 3D mode |
| `dimension 2` | | 2D mode — orient surface normals to +z, z-coordinates set to 0 |
| `overwrite` | off | Overwrite existing file |

> **Order 4-5**: Not supported in Gmsh export (face/volume interior node extraction
> is unreliable, leading to negative Jacobians). Use `radia_export netgen` for order 4-5.

## Supported Elements

### 1st Order (order 1)

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

### 2nd Order (order 2) — Serendipity

| Element Type | Gmsh Code | Nodes |
|--------------|-----------|-------|
| Line3 | 8 | 3 |
| Triangle6 | 9 | 6 |
| Quad8 | 16 | 8 |
| Tetrahedron10 | 11 | 10 |
| Hexahedron20 | 17 | 20 |
| Wedge15 | 18 | 15 |
| Pyramid13 | 19 | 13 |

### 3rd Order (order 3) — Serendipity

| Element Type | Gmsh Code | Nodes |
|--------------|-----------|-------|
| Line4 | 26 | 4 |
| Triangle10 | 21 | 10 |
| Quad12 | 39 | 12 |
| Tetrahedron20 | 29 | 20 |
| Hexahedron32 | 99 | 32 |
| Pyramid21 | 125 | 21 |

> **Note**: Wedge/Prism is **not supported** at order 3 (Gmsh limitation:
> `FaceClosureFull` not implemented for prisms). Wedge elements fall back
> to linear when order 3 is requested.

## File Format

The output uses Gmsh v4.1 format:

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

cubit.cmd('radia_export gmsh "mesh.msh" overwrite')
```

### 2D Export with Normal Orientation

```python
cubit.cmd("create surface rectangle width 1 height 1 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("mesh surface 1")

cubit.cmd('radia_export gmsh "plate.msh" dimension 2 overwrite')
```

### High-Order Export (order 2)

```python
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.3")
cubit.cmd("mesh volume 1")

# No block assignment or element type change needed
cubit.cmd('radia_export gmsh "mesh_o2.msh" order 2 overwrite')
```

### 3rd Order Export

```python
cubit.cmd('radia_export gmsh "mesh_o3.msh" order 3 overwrite')
```

> Order 2+ generates a companion `.geo` file (with `Mesh.NumSubEdges=4`)
> for proper curved element display in the Gmsh GUI.

## $Entities Section

Geometry topology information:

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

- [export_NetgenMesh.md](export_NetgenMesh.md) - Netgen .vol export (order 1-5)
- [Function_Reference.md](Function_Reference.md) - All plugin commands
