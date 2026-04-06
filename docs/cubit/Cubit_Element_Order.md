# Cubit Element Order Control

## Overview

This document explains how to create and work with 1st and 2nd order elements in Coreform Cubit, and how they interact with the radia Cubit plugin export commands and the `radia_cubit_mesh` module.

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

### radia_cubit_mesh.extract_curved_mesh()

Uses `get_connectivity()` to export **only 1st order elements**:

- TET4 ↁE4 nodes (exported)
- TET10 ↁE4 nodes (only corners exported)

This is intentional! Netgen's `mesh.Curve(order)` generates high-order nodes from geometry.

```python
# Even with TET10, exports as 1st order for Netgen
cubit.cmd("block 1 element type tetra10")
import radia_cubit_mesh
ngmesh = radia_cubit_mesh.extract_curved_mesh(order=2)
# ngmesh contains 4-node tets, mesh.Curve() adds high-order nodes
```

### radia export gmsh

Uses `get_expanded_connectivity()` to export 2nd order elements:

- TET4 ↁEGmsh type 4
- TET10 ↁEGmsh type 11

## Design Philosophy

### Why extract_curved_mesh Uses 1st Order Base Mesh

The design philosophy for `radia_cubit_mesh.extract_curved_mesh()`:

1. **Cubit's role**: Generate high-quality 1st order mesh (topology)
2. **Netgen's role**: Add high-order nodes based on CAD geometry

This separation provides:
- Arbitrary curve orders (2, 3, 4, 5, ...) via `mesh.Curve(order)`
- Nodes placed exactly on CAD geometry
- No need for Cubit-Netgen node ordering conversion for high-order nodes

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
cubit.cmd('radia export gmsh "sphere_1st.msh" overwrite')

# Convert to 2nd order
cubit.cmd("block 1 element type tetra10")
cubit.cmd("block 2 element type tri6")

# Check 2nd order
print(f"2nd order - get_connectivity: {len(cubit.get_connectivity('tet', tet_id))}")
print(f"2nd order - get_expanded_connectivity: {len(cubit.get_expanded_connectivity('tet', tet_id))}")

# Export 2nd order Gmsh
cubit.cmd('radia export gmsh "sphere_2nd.msh" overwrite')
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
| `export gmsh "f.msh" order N` | NetgenCurver | 1-4 |
| `export radia_nastran "f.bdf" order N` | NetgenCurver | 1-2 |
| `export vtk "f.vtk" order N` | NetgenCurver | 1-2 |
| `extract_curved_mesh(cubit, order=N)` | CallbackGeometry + BuildCurvedElements | 1-5 |

> **Note**: `export radia_nastran` (NOT `export nastran`). Cubit has a built-in `export nastran` with different format.

## See Also

- [export_NetgenMesh.md](export_NetgenMesh.md) - Cubit to Netgen mesh export with high-order curving (`radia_cubit_mesh.extract_curved_mesh`)
- [Cubit Documentation](https://coreform.com/products/coreform-cubit/documentation/)
