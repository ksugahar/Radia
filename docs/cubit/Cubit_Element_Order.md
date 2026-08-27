# Cubit Element Order Control

## Overview

This document explains how to create and work with 1st and 2nd order elements in Coreform Cubit, and how they interact with the Radia Cubit plugin export commands (`export <fmt>`).

## Creating 2nd Order Elements

### Correct Method: Block Element Type

The correct way to create 2nd order elements in Cubit is:

1. **Create mesh** (generates 1st order elements by default)
2. **Add elements to a block**
3. **Set the element type** to convert to 2nd order

```python
# Step 1: Create mesh
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")

# Step 2: Add elements to block
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 name 'sphere'")

# Step 3: Convert to 2nd order
cubit.cmd("block 1 element type tetra10")
```

### Element Type Commands

| 1st Order | 2nd Order | Command |
|-----------|-----------|---------|
| TET4 | TET10 | `block X element type tetra10` |
| HEX8 | HEX20 | `block X element type hex20` |
| HEX8 | HEX27 | `block X element type hex27` |
| WEDGE6 | WEDGE15 | `block X element type wedge15` |
| PYRAMID5 | PYRAMID13 | `block X element type pyramid13` |
| TRI3 | TRI6 | `block X element type tri6` |
| QUAD4 | QUAD8 | `block X element type quad8` |
| QUAD4 | QUAD9 | `block X element type quad9` |
| EDGE2 | EDGE3 | `block X element type bar3` |

## Python API: get_connectivity vs get_expanded_connectivity

Cubit provides two functions to retrieve element nodes:

### get_connectivity()

Returns **only the corner (1st order) nodes** of an element.

```python
tet_id = cubit.get_block_tets(1)[0]
nodes = cubit.get_connectivity("tet", tet_id)
print(len(nodes))  # Always 4 for tetrahedra
```

### get_expanded_connectivity()

Returns **all nodes including mid-edge nodes** for 2nd order elements.

```python
tet_id = cubit.get_block_tets(1)[0]
nodes = cubit.get_expanded_connectivity("tet", tet_id)
print(len(nodes))  # 4 for TET4, 10 for TET10
```

### Comparison

```python
# For a 2nd order tetrahedron (TET10):
cubit.cmd("block 1 element type tetra10")

tet_id = cubit.get_block_tets(1)[0]

# get_connectivity returns 4 nodes (corners only)
nodes_1st = cubit.get_connectivity("tet", tet_id)
print(f"get_connectivity: {len(nodes_1st)} nodes")  # 4

# get_expanded_connectivity returns 10 nodes (corners + mid-edge)
nodes_2nd = cubit.get_expanded_connectivity("tet", tet_id)
print(f"get_expanded_connectivity: {len(nodes_2nd)} nodes")  # 10
```

## Impact on Export Functions

All plugin `export <fmt>` commands internally extract 1st order elements only, then use
**NetgenCurver** (CallbackGeometry + ACIS projection) to generate high-order nodes
at the requested order. Block element type settings (`tetra10`, etc.) are ignored.

```python
# The order parameter controls high-order curving — not the block element type
cubit.cmd('export netgen "mesh.vol" order 3 overwrite')
cubit.cmd('export gmsh "mesh.msh" order 2 overwrite')
```

## Design Philosophy

1. **Cubit's role**: Generate high-quality 1st order mesh (topology)
2. **NetgenCurver's role**: Add high-order nodes projected onto ACIS CAD surfaces

This separation provides:
- Arbitrary curve orders (1-5 for .vol, 1-3 for .msh, 1-2 for .bdf/.vtk)
- Nodes placed exactly on CAD geometry
- No Cubit-Netgen node ordering conversion needed for high-order nodes

## Complete Example

```python
import os, sys
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Create geometry and mesh
cubit.cmd("reset")
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.3")
cubit.cmd("mesh volume 1")

# Add to blocks
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 name 'sphere'")
cubit.cmd("block 2 add tri all in surface all")
cubit.cmd("block 2 name 'boundary'")

# Check 1st order
tet_id = cubit.get_block_tets(1)[0]
print(f"1st order - get_connectivity: {len(cubit.get_connectivity('tet', tet_id))}")
print(f"1st order - get_expanded_connectivity: {len(cubit.get_expanded_connectivity('tet', tet_id))}")

# Export 1st order Gmsh
cubit.cmd('export gmsh "sphere_1st.msh" overwrite')

# Convert to 2nd order
cubit.cmd("block 1 element type tetra10")
cubit.cmd("block 2 element type tri6")

# Check 2nd order
print(f"2nd order - get_connectivity: {len(cubit.get_connectivity('tet', tet_id))}")
print(f"2nd order - get_expanded_connectivity: {len(cubit.get_expanded_connectivity('tet', tet_id))}")

# Export 2nd order Gmsh
cubit.cmd('export gmsh "sphere_2nd.msh" overwrite')
```

Output:
```
1st order - get_connectivity: 4
1st order - get_expanded_connectivity: 4
2nd order - get_connectivity: 4
2nd order - get_expanded_connectivity: 10
```

## Summary Table

| Export Function | High-Order Method | Max Order |
|----------------|-------------------|-----------|
| `export netgen "f.vol" order N` | NetgenCurver (compact_netgen) | 1-5 |
| `export gmsh "f.msh" order N` | NetgenCurver | 1-3 |
| `export nastran_bdf "f.bdf" order N` | NetgenCurver | 1-2 |
| `export vtk "f.vtk" order N` | NetgenCurver | 1-2 |

> **Note**: use `export nastran_bdf` (not Cubit's built-in `export nastran`).
> `export jmag_nastran` remains a deprecated compatibility alias.

## See Also

- [export_NetgenMesh.md](export_NetgenMesh.md) — Netgen .vol export (order 1-5)
- [Function_Reference.md](Function_Reference.md) — All plugin commands
