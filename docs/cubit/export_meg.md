# export_meg

Export mesh to ELF/MAGIC MEG format.

## Synopsis

```python
cubit_mesh_export.export_meg(cubit, FileName, DIM='T', MGR2=None)
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cubit` | object | required | Cubit Python interface object |
| `FileName` | str | required | Output file path for the .meg file |
| `DIM` | str | `'T'` | Dimension mode (see below) |
| `MGR2` | list | `None` | Optional spatial nodes for boundary conditions |

### DIM Parameter Options

| Value | Description |
|-------|-------------|
| `'T'` | 3D mode (three-dimensional analysis) |
| `'R'` | Axisymmetric mode (2D axisymmetric analysis) |
| `'K'` | 2D mode (plane stress/strain analysis) |

### MGR2 Parameter

Optional list of spatial node coordinates for additional boundary condition nodes:
```python
MGR2 = [[x1, y1, z1], [x2, y2, z2], ...]
```

These nodes are written as MGR2 records in the MEG file.

## Returns

Returns the `cubit` object for method chaining.

## Supported Elements

### 3D Elements (DIM='T')

| Cubit Element | MEG Element | Nodes | Description |
|---------------|-------------|-------|-------------|
| Tetrahedron | TET | 4 | 4-node tetrahedral solid |
| Hexahedron | HEX | 8 | 8-node hexahedral solid |
| Wedge | WDG | 6 | 6-node wedge/prism solid |
| Pyramid | PYR | 5 | 5-node pyramid solid |

### 2D Elements (DIM='K' or 'R')

| Cubit Element | MEG Element | Nodes | Description |
|---------------|-------------|-------|-------------|
| Triangle | TRI | 3 | 3-node triangular element |
| Quadrilateral | QUA | 4 | 4-node quadrilateral element |

### 1D Elements

| Cubit Element | MEG Element | Nodes | Description |
|---------------|-------------|-------|-------------|
| Edge | BAR | 2 | 2-node bar/edge element |

## Limitations

- **First-order elements only**: Second-order elements are not supported
- Elements are exported with 1st order connectivity only

## File Format

The MEG format is used by ELF/MAGIC electromagnetic solver:

```
BOOK  MEP  3.50
* ELF/MESH VERSION 7.3.0
* SOLVER = ELF/MAGIC
MGSC 0.001
* NODE
MGR1 1 0 0.0 0.0 0.0
MGR1 2 0 1.0 0.0 0.0
...
* ELEMENT
MGE1 1 1 TET 1 2 3 4
...
MGEND
```

### Record Types

| Record | Description |
|--------|-------------|
| BOOK | File format identifier |
| MGSC | Scale factor |
| MGR1 | Node definition (id, flag, x, y, z) |
| MGR2 | Spatial node for boundary conditions |
| MGE1 | Element definition (id, material, type, nodes...) |
| MGEND | End of file marker |

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

cubit_mesh_export.export_meg(cubit, "mesh.meg", DIM='T')
```

### 2D Planar Export

```python
cubit.cmd("create surface rectangle width 1 height 1 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("mesh surface 1")
cubit.cmd("block 1 add tri all")

cubit_mesh_export.export_meg(cubit, "plate.meg", DIM='K')
```

### Axisymmetric Export

```python
cubit.cmd("create surface rectangle width 1 height 1 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("mesh surface 1")
cubit.cmd("block 1 add tri all")

cubit_mesh_export.export_meg(cubit, "axisym.meg", DIM='R')
```

### With Spatial Nodes

```python
# Define spatial nodes for boundary conditions
spatial_nodes = [
    [0.0, 0.0, 10.0],  # Far-field node 1
    [0.0, 0.0, -10.0]  # Far-field node 2
]

cubit_mesh_export.export_meg(cubit, "mesh.meg", DIM='T', MGR2=spatial_nodes)
```

## Compatibility

The output file is compatible with:
- ELF/MAGIC electromagnetic solver
- ELF/MESH mesh editor

## See Also

- [Cubit_Element_Order.md](Cubit_Element_Order.md) - Element order control (note: MEG export is 1st order only)
- [export_Nastran](export_Nastran.md) - Alternative format for structural analysis
