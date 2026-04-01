# Exodus II Export Examples

Export Cubit mesh to Exodus II format (.exo).

## Overview

Exodus II is Cubit's native mesh format, providing full feature support:

- All element types (tet, hex, wedge, pyramid, tri, quad)
- 1st and 2nd order elements
- Nodesets and Sidesets for boundary conditions
- Block/material definitions
- Large model support (64-bit integers)

## Usage

```python
import cubit
import radia_cubit_mesh

cubit.cmd("block 1 add tet all")
cubit.cmd('export mesh "output.exo" overwrite')
```

### Large Model Option

For meshes exceeding 2^31 elements/nodes:

```python
cubit.cmd('export mesh "large.exo" overwrite')
```

## Sample Files

| File | Description |
|------|-------------|
| `cube_tet.exo` | Simple tetrahedral mesh |
| `cube_hex.exo` | Simple hexahedral mesh |
| `sphere_2nd_order.exo` | 2nd order mesh (TET10, TRI6) |
| `mixed_with_bc.exo` | Mixed elements with nodesets/sidesets |
| `large_model.exo` | Large model with 64-bit integers |

## Regenerate Samples

```bash
"${CUBIT_PATH:-C:/Program Files/Coreform Cubit 2025.3/bin}/python3/python.exe" exodus_export_example.py
```

## Compatible Software

- ParaView (with Exodus plugin)
- Cubit (import mesh)
- SIERRA tools
- VisIt

## See Also

- [docs/cubit/export_exodus.md](../../../docs/cubit/export_exodus.md) - Full documentation
