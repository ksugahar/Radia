# Function Reference

Reference documentation for `cubit_mesh_export` module functions.

## Export Functions Overview

| Function | Format | 1st Order | 2nd Order | 3rd+ Order |
|----------|--------|-----------|-----------|------------|
| `export_exodus()` | Exodus II (.exo) | Yes | Yes | Yes |
| `export_netgen()` | Netgen mesh object | Yes | Yes (via Curve) | Yes (via Curve) |
| `export_netgen_with_names()` | Netgen mesh object | Yes | Yes (via Curve) | Yes (via Curve) |
| `export_vtk()` | VTK Legacy (.vtk) | Yes | Yes | No |
| `export_vtu()` | VTK XML (.vtu) | Yes | Yes | No |
| `export_gmesh()` | Gmsh v2.2/v4.1 | Yes | Yes | No |
| `export_nastran()` | Nastran BDF | Yes | No | No |
| `export_meg()` | MEG (ELF) | Yes | No | No |

---

## Exodus II Export

```python
export_exodus(cubit, FileName, large_model=False)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `FileName` | str | required | Output .exo file path |
| `large_model` | bool | False | Use 64-bit integers for large meshes |

**Features**: All element types, 1st/2nd order, nodesets, sidesets, block definitions.

[Full documentation](export_exodus.md) | [Examples](../../examples/cubit/exodus/)

---

## Gmsh Export

```python
export_gmesh(cubit, FileName, version="2.2", DIM="auto")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `FileName` | str | required | Output .msh file path |
| `version` | str | "2.2" | Format version: "2.2" or "4.1" |
| `DIM` | str | "auto" | Dimension mode (v4.1 only): "auto", "2D", or "3D" |

**DIM Options** (v4.1 only):
| Value | Description |
|-------|-------------|
| `"auto"` | Auto-detect (3D if volume elements exist) |
| `"2D"` | Orient normals to +z, z-coordinates set to 0 |
| `"3D"` | No normal orientation |

### v2.2 vs v4.1

| Feature | v2.2 | v4.1 |
|---------|------|------|
| $Entities section | No | Yes |
| DIM parameter | No | Yes |
| NGSolve/Netgen | **Supported** | Not recommended |
| Radia mesh import | **Supported** | Not supported |
| GMSH visualization | Supported | **Recommended** |

[v2.2 documentation](export_Gmsh_ver2.md) | [v4.1 documentation](export_Gmsh_ver4.md) | [Examples](../../examples/cubit/gmsh/)

---

## Nastran BDF Export

```python
export_nastran(cubit, FileName, DIM="3D", PYRAM=True)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `FileName` | str | required | Output .bdf file path |
| `DIM` | str | "3D" | "3D" or "2D" |
| `PYRAM` | bool | True | Pyramid handling |

**DIM Options**:
| Value | Elements |
|-------|----------|
| `"3D"` | CTETRA, CHEXA, CPENTA, CPYRAM |
| `"2D"` | CTRIA3, CQUAD4 (normals to +z) |

**PYRAM Options**:
| Value | Output | Use Case |
|-------|--------|----------|
| `True` | CPYRAM (5-node) | Standard Nastran |
| `False` | Degenerate CHEXA | JMAG compatibility |

**Limitation**: 1st order elements only.

[Full documentation](export_Nastran.md) | [Examples](../../examples/cubit/nastran/)

---

## MEG Export (ELF/MAGIC)

```python
export_meg(cubit, FileName, DIM='T', MGR2=None)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `FileName` | str | required | Output .meg file path |
| `DIM` | str | 'T' | 'T', 'K', or 'R' |
| `MGR2` | list | None | Spatial nodes [[x,y,z], ...] |

**DIM Options**:
| Value | Description | Coordinate System |
|-------|-------------|-------------------|
| `'T'` | 3D | X, Y, Z |
| `'K'` | 2D Planar | X, Y (Z=0) |
| `'R'` | Axisymmetric | R (X), Z |

**Block Names = ELF Element Names**:
| DIM | Tri | Quad | Tet | Wedge | Hex |
|-----|-----|------|-----|-------|-----|
| `'T'` | - | - | MMB4T | MMB6T | MMB8T |
| `'K'` | MMB3K | MMB4K | - | - | - |
| `'R'` | MMB3R | MMB4R | - | - | - |

**Limitation**: 1st order elements only.

[Full documentation](export_meg.md) | [Examples](../../examples/cubit/meg/)

---

## VTK Export

### Legacy Format (.vtk)

```python
export_vtk(cubit, FileName)
```

### XML Format (.vtu)

```python
export_vtu(cubit, FileName)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `FileName` | str | required | Output file path |

**VTK vs VTU**:
| Feature | VTK Legacy | VTK XML |
|---------|------------|---------|
| Format | ASCII text | XML |
| Metadata | Basic | BlockID, NodeID |
| ParaView | Compatible | **Recommended** |

**Auto-detection**: Element order detected from node count.

[Full documentation](export_VTK.md) | [Examples](../../examples/cubit/vtk/)

---

## Netgen Export

### Standard Export

```python
export_netgen(cubit, geometry_file=None, geometry=None)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `geometry_file` | str | None | Path to STEP/BREP/IGES file for Curve() support |
| `geometry` | OCCGeometry | None | OCC geometry object (takes precedence over geometry_file) |

**Returns**: `netgen.meshing.Mesh` object.

**Use case**: Simple geometries (cylinder, sphere, torus, cone) or any geometry loaded from STEP.

### Name-based Export (Complex Geometry)

```python
export_netgen_with_names(cubit, geometry)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface |
| `geometry` | OCCGeometry | required | OCC geometry with named faces |

**Returns**: `netgen.meshing.Mesh` object with exact OCC face mapping.

**Use case**: Complex geometries with Boolean operations (e.g., brick with cylindrical hole). Requires prior face naming with `name_occ_faces()`.

### Choosing the Right Workflow

| Geometry | Recommended Function |
|----------|---------------------|
| Simple shape (cylinder, sphere, etc.) | `export_netgen()` |
| Complex shape (Boolean ops) | `export_netgen_with_names()` |
| 2nd order only (no Curve(3+)) | `export_gmsh_v2()` + `ReadGmsh` |

[Full documentation](export_NetgenMesh.md) | [Examples](../../examples/cubit/netgen/)

---

## OCC Face Naming Utility

```python
name_occ_faces(shape, prefix="occ_face_")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `shape` | OCC shape | required | OCC geometry shape object |
| `prefix` | str | "occ_face_" | Prefix for face names |

**Returns**: Number of faces named.

**Use case**: Assign unique names to OCC faces before STEP export, so that `export_netgen_with_names()` can map Cubit surfaces to OCC faces after reimport.

---

## SetGeomInfo UV Utilities

These functions set UV parameters on Netgen mesh elements, enabling `mesh.Curve(order)` for high-order geometry approximation on curved surfaces.

### Surface Detection Functions

```python
set_cylinder_geominfo(ngmesh, radius, height, center=(0,0,0), axis='z', tol=0.01)
set_sphere_geominfo(ngmesh, radius, center=(0,0,0), tol=0.01)
set_torus_geominfo(ngmesh, major_radius, minor_radius, center=(0,0,0), axis='z', tol=0.01)
set_cone_geominfo(ngmesh, base_radius, height, center=(0,0,0), axis='z', tol=0.01)
```

**Returns**: Number of vertex geominfo entries modified.

| Parameter | Description |
|-----------|-------------|
| `ngmesh` | Netgen mesh object |
| `radius` / `base_radius` | Surface radius |
| `height` | Cylinder/cone height |
| `major_radius`, `minor_radius` | Torus radii |
| `center` | Center coordinates |
| `axis` | 'x', 'y', or 'z' |
| `tol` | Tolerance for detecting surface vertices |

### UV Computation Functions

```python
compute_cylinder_uv(x, y, z, radius, height, center=(0,0,0), axis='z')
compute_sphere_uv(x, y, z, radius, center=(0,0,0))
compute_torus_uv(x, y, z, major_radius, minor_radius, center=(0,0,0), axis='z')
compute_cone_uv(x, y, z, base_radius, height, center=(0,0,0), axis='z')
```

**Returns**: `(u, v)` tuple of UV parameters for the given point.

[Full documentation](export_NetgenMesh.md) | [Examples](../../examples/cubit/netgen/)

---

## Technical Guides

| Document | Description |
|----------|-------------|
| [Cubit_Element_Order.md](Cubit_Element_Order.md) | How to control element order (1st/2nd) in Cubit |

### Key Concepts

- **`get_connectivity()`** - Returns 1st order nodes only (corner nodes)
- **`get_expanded_connectivity()`** - Returns all nodes including mid-edge nodes
- **`block X element type tetra10`** - Cubit command to convert to 2nd order

### Internal Helper Functions

These functions are used internally by all export functions:

| Function | Description |
|----------|-------------|
| `_block_contains_geometry(cubit, block_id)` | Returns True if block contains geometry (volume/surface/curve/vertex) |
| `_get_block_elements(cubit, block_id, elem_type)` | Gets mesh elements from block, supporting both geometry and mesh element blocks |
| `_warn_mixed_element_types_in_blocks(cubit)` | Warns if blocks contain multiple 3D or 2D element types |

### Block Types

Blocks can contain either mesh elements or geometry:

| Block Contains | Elements Returned by `_get_block_elements()` |
|----------------|----------------------------------------------|
| Volume | tet, hex, wedge, pyramid (3D only) |
| Surface | tri, quad (2D only) |
| Curve | edge (1D only) |
| Vertex | node (0D only) |
| Mesh elements | Only the registered element types |

---

## Bug Fix Records

| Document | Description |
|----------|-------------|
| [BUGFIX_cubit_mesh_export_vtk_indexing.md](BUGFIX_cubit_mesh_export_vtk_indexing.md) | VTK node indexing fix for non-contiguous node IDs |
