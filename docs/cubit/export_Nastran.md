# export_Nastran

Export mesh to NX Nastran bulk data format.

## Synopsis

```python
cubit_mesh_export.export_Nastran(cubit, FileName, DIM="3D", PYRAM=True)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface object |
| `FileName` | str | required | Output file path for the .nas or .bdf file |
| `DIM` | str | `"3D"` | Dimension mode (see below) |
| `PYRAM` | bool | `True` | Pyramid element handling (see below) |

### DIM Parameter Options

| Value | Description |
|-------|-------------|
| `"3D"` | 3D mode - exports volume elements (CTETRA, CHEXA, CPENTA, CPYRAM) |
| `"2D"` | 2D mode - exports surface elements (CTRIA3, CQUAD4), orients normals to +z |

### PYRAM Parameter Options

| Value | Description |
|-------|-------------|
| `True` | Export pyramids as CPYRAM elements |
| `False` | Convert pyramids to degenerate CHEXA (8-node hex with duplicate nodes) |

**Background**: In hybrid meshes combining tetrahedra and hexahedra, pyramid elements are required at the interface between element types. However, some solvers (e.g., JMAG) cannot import pyramid elements from Nastran files and interpret them as degenerate hexahedra. Use `PYRAM=False` for compatibility with such solvers.

## Returns

Returns the `cubit` object for method chaining.

## Supported Elements

### 3D Elements (DIM="3D")

| Cubit Element | Nastran Card | Nodes | Description |
|---------------|--------------|-------|-------------|
| Tetrahedron | CTETRA | 4 | 4-node tetrahedral solid |
| Hexahedron | CHEXA | 8 | 8-node hexahedral solid |
| Wedge | CPENTA | 6 | 6-node pentahedral (prism) solid |
| Pyramid | CPYRAM | 5 | 5-node pyramid solid |

### 2D Elements (DIM="2D")

| Cubit Element | Nastran Card | Nodes | Description |
|---------------|--------------|-------|-------------|
| Triangle | CTRIA3 | 3 | 3-node triangular shell |
| Quadrilateral | CQUAD4 | 4 | 4-node quadrilateral shell |

### 1D Elements

| Cubit Element | Nastran Card | Nodes | Description |
|---------------|--------------|-------|-------------|
| Edge/Bar | CROD | 2 | 2-node rod/bar element |

### 0D Elements

| Cubit Element | Nastran Card | Nodes | Description |
|---------------|--------------|-------|-------------|
| Node | CMASS | 1 | Point mass element |

## Limitations

- **First-order elements only**: Second-order elements (TETRA10, HEX20, etc.) are not supported
- Elements are exported as 1st order regardless of Cubit element type setting

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
import cubit_mesh_export

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")

cubit_mesh_export.export_Nastran(cubit, "mesh.nas", DIM="3D")
```

### 2D Plate Export

```python
cubit.cmd("create surface rectangle width 1 height 1 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("mesh surface 1")
cubit.cmd("block 1 add tri all")

# Export as 2D - normals oriented to +z
cubit_mesh_export.export_Nastran(cubit, "plate.nas", DIM="2D")
```

### Handling Pyramids

```python
# Export pyramids as CPYRAM (default)
cubit_mesh_export.export_Nastran(cubit, "mesh.nas", PYRAM=True)

# Convert pyramids to degenerate hex (for solver compatibility)
cubit_mesh_export.export_Nastran(cubit, "mesh.nas", PYRAM=False)
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

- [Cubit_Element_Order.md](Cubit_Element_Order.md) - Element order control (note: Nastran export is 1st order only)
- [export_Gmsh_ver4](export_Gmsh_ver4.md) - Alternative format with 2nd order support
