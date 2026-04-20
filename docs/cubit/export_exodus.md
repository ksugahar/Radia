# export_exodus

Export mesh to Exodus II format (.exo, .e, .g).

## Overview

Exodus II is Cubit's native mesh format and the standard format for SIERRA tools developed at Sandia National Laboratories. This function provides a consistent API wrapper around Cubit's built-in export functionality.

## Command

Exodus II export uses Cubit's built-in export command directly:

```python
cubit.cmd('export mesh "output.exo" overwrite')
```

> **Note**: The old `cubit_mesh_export.export_exodus()` Python function has been removed. Use Cubit's native `export mesh` command instead.

## Why no Radia plugin for Exodus?

Exodus II is Cubit's native format and already supports arbitrary
element orders through Cubit's own block element-type settings
(`HEX20`, `HEX27`, `TET10`, etc).  The Radia plugin's main value on
other formats is the `NetgenCurver`-based ACIS projection for high
order curving — Exodus does not need that because Cubit already
produces curved-element Exodus files from its own high-order meshes.
`radia_export exodus` does not exist; use `export mesh` directly.

## Supported Elements

### 3D Elements
| Element | 1st Order | 2nd Order | 3rd Order |
|---------|-----------|-----------|-----------|
| Tetrahedron | TET4 | TET10 | TET11 |
| Hexahedron | HEX8 | HEX20 | HEX27 |
| Wedge | WEDGE6 | WEDGE15 | - |
| Pyramid | PYRAMID5 | PYRAMID13 | - |

### 2D Elements
| Element | 1st Order | 2nd Order | 3rd Order |
|---------|-----------|-----------|-----------|
| Triangle | TRI3 | TRI6 | TRI7 |
| Quadrilateral | QUAD4 | QUAD8 | QUAD9 |

### 1D/0D Elements
- BAR2, BAR3 (edge elements)
- NODE (point elements)

## Usage Examples

### Basic Usage

```python
import os, sys
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)
import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Create geometry and mesh
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.2")
cubit.cmd("mesh volume 1")

# Define block
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'solid'")

# Export to Exodus
cubit.cmd('export mesh "output.exo" overwrite')
```

### 2nd Order Elements

```python
# Create mesh
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")

# Add mesh elements to block (not geometry!)
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 name 'sphere'")

# Convert to 2nd order
cubit.cmd("block 1 element type tetra10")

# Export
cubit.cmd('export mesh "sphere_2nd_order.exo" overwrite')
```

### With Nodesets and Sidesets

```python
# Create geometry
cubit.cmd("create brick x 2 y 1 z 1")
cubit.cmd("volume 1 scheme map")
cubit.cmd("mesh volume 1")

# Define block
cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'solid'")

# Define nodesets (boundary conditions)
cubit.cmd("nodeset 1 add node in surface 1")
cubit.cmd("nodeset 1 name 'fixed'")
cubit.cmd("nodeset 2 add node in surface 2")
cubit.cmd("nodeset 2 name 'load'")

# Define sidesets (surface loads)
cubit.cmd("sideset 1 add surface 1")
cubit.cmd("sideset 1 name 'inlet'")
cubit.cmd("sideset 2 add surface 2")
cubit.cmd("sideset 2 name 'outlet'")

# Export (nodesets and sidesets included automatically)
cubit.cmd('export mesh "with_bc.exo" overwrite')
```

### Large Model Support

For meshes with more than 2^31 (~2 billion) elements or nodes:

```python
cubit.cmd('export mesh "large_model.exo" overwrite')
```

> For large model support (64-bit integers), consult Cubit's documentation for the `large` option.

## Output Summary

The function prints a summary of the exported mesh:

```
Exodus export: output.exo
--------------------------------------------------
Nodes: 1234
Tetrahedra: 5678

Blocks: 1
  Block 1: "solid"

Nodesets: 2
  Nodeset 1: "fixed"
  Nodeset 2: "load"

Sidesets: 2
  Sideset 1: "inlet"
  Sideset 2: "outlet"
--------------------------------------------------
```

## Comparison with Other Formats

| Feature | Exodus II | Gmsh | VTK | Nastran |
|---------|-----------|------|-----|---------|
| Native to Cubit | Yes | No | No | No |
| High-order elements | Yes (all) | Yes (2nd) | Yes (2nd) | No |
| Nodesets | Yes | No | No | No |
| Sidesets | Yes | No | No | No |
| Block names | Yes | Yes | Yes | Yes |
| Large model (64-bit) | Yes | No | No | No |

## Compatible Software

- **ParaView** - With Exodus reader plugin
- **VisIt** - Native support
- **SIERRA** - Sandia's simulation tools
- **Cubit** - Can re-import for further processing
- **Abaqus** - Via conversion tools

## Notes

1. **Nodesets and Sidesets**: Unlike other export functions in this module, Exodus natively supports nodesets and sidesets. They are automatically included in the export.

2. **Element Order**: For 2nd/3rd order elements, blocks must contain mesh elements (not geometry). Use `block X element type tetra10` to convert.

3. **File Extensions**: Common extensions are `.exo` (Exodus), `.e` (Exodus), and `.g` (Genesis input format).

## See Also

- [Cubit_Element_Order.md](Cubit_Element_Order.md) - How to control element order
- [Function_Reference.md](Function_Reference.md) - All export functions
