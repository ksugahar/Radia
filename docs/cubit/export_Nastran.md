# export_Nastran

Export mesh to NX Nastran bulk data format.

## Cubit Plugin (Recommended)

The Radia Cubit plugin provides a native APREPRO command — no Python or block assignment required:

```
radia_export nastran "mesh.bdf"
radia_export nastran "mesh.bdf" dimension 2
radia_export nastran "mesh.bdf" nopyramid
radia_export nastran "mesh.bdf" dimension 2 nopyramid overwrite
```

### Syntax

```
radia_export nastran <"filename"> [order <1|2>] [dimension <2|3>] [nopyramid] [overwrite]
```

> **IMPORTANT**: Use `radia_export nastran`, NOT `export nastran`.
> Cubit has a built-in `export nastran` command with different format and no order 2 support.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `order 1` | yes | 1st order elements (CTETRA, CHEXA, CPENTA, CPYRAM) |
| `order 2` | | 2nd order elements (CTETRA10, CHEXA20, via NetgenCurver) |
| `dimension 3` | yes | 3D solid mesh |
| `dimension 2` | | 2D shell mesh (CTRIA3/CTRIA6, CQUAD4/CQUAD8), normals oriented to +z |
| `nopyramid` | off | Convert pyramid elements to degenerate CHEXA (JMAG compatible) |
| `overwrite` | off | Overwrite existing file without warning |

**Advantages over Python API:**
- No `block` assignment needed — exports all mesh elements automatically
- No `#!python` block — pure APREPRO/journal command
- Usable directly in `.jou` files

**Installation:** Copy `cubit_mesh_export.ccm` to `<Cubit install>/bin/plugins/` or set `CUBIT_PLUGIN_DIR`.

**Build:** See [src/cubit_plugin/CMakeLists.txt](../../src/cubit_plugin/CMakeLists.txt).

---

## Usage from Python

```python
cubit.cmd('radia_export nastran "mesh.bdf" order 2 overwrite')
cubit.cmd('radia_export nastran "mesh.bdf" dimension 2 nopyramid overwrite')
```

**Background**: In hybrid meshes combining tetrahedra and hexahedra, pyramid elements are required at the interface between element types. However, some solvers (e.g., JMAG) cannot import pyramid elements from Nastran files and interpret them as degenerate hexahedra. Use `nopyramid` for compatibility with such solvers.

## Supported Elements

### 3D Elements — Order 1

| Cubit Element | Nastran Card | Nodes |
|---------------|--------------|-------|
| Tetrahedron | CTETRA | 4 |
| Hexahedron | CHEXA | 8 |
| Wedge | CPENTA | 6 |
| Pyramid | CPYRAM | 5 |

### 3D Elements — Order 2 (via NetgenCurver)

| Cubit Element | Nastran Card | Nodes |
|---------------|--------------|-------|
| Tetrahedron | CTETRA | 10 |
| Hexahedron | CHEXA | 20 |
| Wedge | CPENTA | 15 |
| Pyramid | CPYRAM | 13 |

### 2D Elements

| Order | Triangle | Quadrilateral |
|-------|----------|---------------|
| 1 | CTRIA3 (3) | CQUAD4 (4) |
| 2 | CTRIA6 (6) | CQUAD8 (8) |

### 1D / 0D Elements

| Cubit Element | Nastran Card | Nodes |
|---------------|--------------|-------|
| Edge/Bar | CROD | 2 |
| Node | CMASS | 1 |

## File Format

The output file follows NX Nastran bulk data format:

1. **Header comments** - File info, timestamp, warnings
2. **Executive Control** - SOL card placeholder
3. **Case Control** - Analysis setup placeholder
4. **Bulk Data** - GRID and element cards

### Example Output Structure

```
$ CUBIT NX Nastran Translator
$ File: mesh.nas
$ Time Stamp: 29-Nov-24 at 12:00:00
$
BEGIN BULK
GRID,1,0,0.0,0.0,0.0,0
GRID,2,0,1.0,0.0,0.0,0
...
CTETRA,1,1,1,2,3,4
...
ENDDATA
```

## Usage Examples

### Basic 3D Export

```python
import cubit

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")

cubit.cmd('radia_export nastran "mesh.bdf" dimension 3 overwrite')
```

### 2D Plate Export

```python
cubit.cmd("create surface rectangle width 1 height 1 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("mesh surface 1")

# Export as 2D - normals oriented to +z
cubit.cmd('radia_export nastran "plate.bdf" dimension 2 overwrite')
```

### Handling Pyramids

```python
# Export pyramids as CPYRAM (default)
cubit.cmd('radia_export nastran "mesh.bdf" overwrite')

# Convert pyramids to degenerate hex (for solver compatibility)
cubit.cmd('radia_export nastran "mesh.bdf" nopyramid overwrite')
```

## 2D Mode Normal Orientation

When `DIM="2D"` is specified:
- Surface element normals are checked and reoriented to point in +z direction
- Z-coordinates are set to 0 for all nodes
- This ensures consistent normal direction for shell element analysis

## Compatibility

The output file is compatible with:
- NX Nastran
- MSC Nastran
- OptiStruct
- Other Nastran-compatible solvers

## See Also

- [export_Gmsh.md](export_Gmsh.md) — Gmsh v4.1 export (order 1-3)
- [export_NetgenMesh.md](export_NetgenMesh.md) — Netgen .vol export (order 1-5)
- [Function_Reference.md](Function_Reference.md) — All plugin commands
