# Function Reference

Reference documentation for the Radia Cubit plugin mesh export commands and the `radia_cubit_mesh` module.

## Cubit Plugin Commands (Recommended)

Native APREPRO commands — no Python, no block assignment required:

| Command | Format | Block Required |
|---------|--------|---------------|
| `export radia_nastran "f.bdf" [dimension <2\|3>] [nopyramid]` | Nastran BDF | **No** |
| `radia export gmsh "f.msh" [version <2\|4>] [dimension <2\|3>]` | Gmsh v2.2/v4.1 | **No** |
| `radia export meg "f.meg"` | MEG (ELF/MAGIC) | **No** |

**Installation:** Copy `radia_cubit.ccm` to `<Cubit install>/bin/plugins/` or set `CUBIT_PLUGIN_DIR`.

**Build:**

```bash
cmake -G "Visual Studio 17 2022" -A x64 \
  -DCubit_DIR="C:/Program Files/Coreform Cubit 2025.3/cmake" \
  ../src/cubit_plugin
cmake --build . --config Release
```

---

## Python API Functions (Legacy)

> **Note**: The `cubit_mesh_export` Python module has been replaced by the `radia_cubit_mesh` C++ pybind11 module (`src/cubit_plugin/radia_cubit_pybind.cpp`). For file exports, use the Cubit Plugin commands above (`cubit.cmd('radia export ...')`). For NGSolve mesh extraction, use `radia_cubit_mesh.extract_curved_mesh()`.

| Function | Format | 1st Order | 2nd Order | 3rd+ Order |
|----------|--------|-----------|-----------|------------|
| `cubit.cmd('export mesh ...')` | Exodus II (.exo) | Yes | Yes | Yes |
| `radia_cubit_mesh.extract_curved_mesh()` | Netgen mesh object | Yes | Yes (via Curve) | Yes (via Curve) |
| `cubit.cmd('radia export gmsh ...')` | Gmsh v2.2/v4.1 | Yes | Yes | No |
| `cubit.cmd('export radia_nastran ...')` | Nastran BDF | Yes | No | No |
| `cubit.cmd('radia export meg ...')` | MEG (ELF) | Yes | No | No |

> Note: The `radia export` plugin commands require no block assignment. The `radia_cubit_mesh.extract_curved_mesh()` function requires blocks for NGSolve mesh extraction.

---

## Exodus II Export

```python
cubit.cmd('export mesh "output.exo" overwrite')
```

Exodus II is Cubit's native format. Use Cubit's built-in export command directly.

**Features**: All element types, 1st/2nd order, nodesets, sidesets, block definitions.

[Full documentation](export_exodus.md) | [Examples](../../examples/cubit_mesh_export/exodus/)

---

## Gmsh Export

```python
cubit.cmd('radia export gmsh "mesh.msh" overwrite')              # v2.2 (default)
cubit.cmd('radia export gmsh "mesh.msh" version 4 overwrite')    # v4.1
```

### v2.2 vs v4.1

| Feature | v2.2 | v4.1 |
|---------|------|------|
| $Entities section | No | Yes |
| DIM parameter | No | Yes |
| NGSolve/Netgen | **Supported** | Not recommended |
| Radia mesh import | **Supported** | Not supported |
| GMSH visualization | Supported | **Recommended** |

[v2.2 documentation](export_Gmsh_ver2.md) | [v4.1 documentation](export_Gmsh_ver4.md) | [Examples](../../examples/cubit_mesh_export/gmsh/)

---

## Nastran BDF Export

```python
cubit.cmd('export radia_nastran "mesh.bdf" dimension 3 overwrite')
```

**DIM Options**:
| Value | Elements |
|-------|----------|
| `dimension 3` | CTETRA, CHEXA, CPENTA, CPYRAM |
| `dimension 2` | CTRIA3, CQUAD4 (normals to +z) |

**PYRAM Options**:
| Option | Output | Use Case |
|--------|--------|----------|
| (default) | CPYRAM (5-node) | Standard Nastran |
| `nopyramid` | Degenerate CHEXA | JMAG compatibility |

**Limitation**: 1st order elements only.

[Full documentation](export_Nastran.md) | [Examples](../../examples/cubit_mesh_export/nastran/)

---

## MEG Export (ELF/MAGIC)

```python
cubit.cmd('radia export meg "mesh.meg" overwrite')
```

**Block Names = ELF Element Names**:
| DIM | Tri | Quad | Tet | Wedge | Hex |
|-----|-----|------|-----|-------|-----|
| `'T'` | - | - | MMB4T | MMB6T | MMB8T |
| `'K'` | MMB3K | MMB4K | - | - | - |
| `'R'` | MMB3R | MMB4R | - | - | - |

**Limitation**: 1st order elements only.

[Full documentation](export_meg.md) | [Examples](../../examples/cubit_mesh_export/meg/)

---

## Netgen Export (with High-Order Curving)

```python
import radia_cubit_mesh
ngmesh = radia_cubit_mesh.extract_curved_mesh(order=2)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `order` | int | 2 | Curve order for high-order elements |

**Returns**: `netgen.meshing.Mesh` object with high-order curving applied.

**Features**: Automatically detects curved surfaces (cylinders, spheres, tori, cones), sets UV parameters, and applies `mesh.Curve(order)`.

| Geometry | Recommended Approach |
|----------|---------------------|
| Simple shape (cylinder, sphere, etc.) | `radia_cubit_mesh.extract_curved_mesh(order=N)` |
| 2nd order only (no Curve(3+)) | `cubit.cmd('radia export gmsh ...')` with 2nd-order blocks + `ReadGmsh` |

[Full documentation](export_NetgenMesh.md) | [Examples](../../examples/cubit_mesh_export/netgen/)

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

