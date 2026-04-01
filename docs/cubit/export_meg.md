# export_meg

Export mesh to ELF/MAGIC MEG format.

## Cubit Plugin (Recommended)

```
radia export meg "mesh.meg"
```

No block assignment or `#!python` required.

**Installation:** Copy `radia_cubit.ccm` to `<Cubit install>/bin/plugins/`.

---

## Plugin Command

```python
cubit.cmd('radia export meg "mesh.meg" overwrite')
```

> **Note**: The old `cubit_mesh_export.export_meg()` Python function has been replaced by the `radia export meg` plugin command. The old Python module (`src/radia/cubit_mesh_export.py`) has been replaced by the C++ pybind11 module (`src/cubit_plugin/radia_cubit_pybind.cpp`).

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

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")

cubit.cmd('radia export meg "mesh.meg" overwrite')
```

### 2D Planar Export

```python
cubit.cmd("create surface rectangle width 1 height 1 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("mesh surface 1")
cubit.cmd("block 1 add tri all")

cubit.cmd('radia export meg "plate.meg" overwrite')
```

### Axisymmetric Export

```python
cubit.cmd("create surface rectangle width 1 height 1 zplane")
cubit.cmd("surface 1 scheme trimesh")
cubit.cmd("mesh surface 1")
cubit.cmd("block 1 add tri all")

cubit.cmd('radia export meg "axisym.meg" overwrite')
```

## Compatibility

The output file is compatible with:
- ELF/MAGIC electromagnetic solver
- ELF/MESH mesh editor

## See Also

- [Cubit_Element_Order.md](Cubit_Element_Order.md) - Element order control (note: MEG export is 1st order only)
- [export_Nastran](export_Nastran.md) - Alternative format for structural analysis
