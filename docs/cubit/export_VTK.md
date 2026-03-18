# VTK Export - Cubit to VTK Mesh Export

## Overview

The `cubit_mesh_export` module provides two VTK export functions:

| Function | Format | Description |
|----------|--------|-------------|
| `export_vtk()` | Legacy VTK (.vtk) | ASCII format, widely compatible |
| `export_vtu()` | VTK XML (.vtu) | Modern XML format, recommended for ParaView |

Both functions automatically detect element order (1st or 2nd) based on node count, supporting mixed-order meshes without user configuration.

---

# export_vtk - Legacy VTK Format

`export_vtk()` exports Cubit mesh data to Legacy VTK (Visualization Toolkit) format. This function automatically detects element order (1st or 2nd) based on node count, supporting mixed-order meshes without user configuration.

## Key Features

- **Automatic Order Detection**: No need to specify element order - detected from node count
- **Mixed-Order Support**: 1st and 2nd order elements can coexist in the same file
- **All Element Types**: Tet, Hex, Wedge, Pyramid, Tri, Quad, Edge, Point
- **VTK Legacy Format**: ASCII format compatible with ParaView and other visualization tools

## Usage

### Basic Usage

```python
import cubit
import cubit_mesh_export

# Create mesh in Cubit
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")

# Export to VTK (auto-detects 1st order)
cubit_mesh_export.export_vtk(cubit, "mesh.vtk")
```

### 2nd Order Elements

```python
# Create mesh
cubit.cmd("create sphere radius 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")

# Add to block and convert to 2nd order
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 name 'sphere'")
cubit.cmd("block 1 element type tetra10")  # Convert to TET10

cubit.cmd("block 2 add tri all in surface all")
cubit.cmd("block 2 name 'boundary'")
cubit.cmd("block 2 element type tri6")  # Convert to TRI6

# Export (auto-detects 2nd order)
cubit_mesh_export.export_vtk(cubit, "sphere_2nd_order.vtk")
```

## Function Signature

```python
def export_vtk(cubit, FileName: str):
    """Export mesh to Legacy VTK format.

    Args:
        cubit: Cubit Python interface object
        FileName: Output file path for the .vtk file

    Returns:
        cubit: The cubit object (for method chaining)
    """
```

## Automatic Order Detection

The function automatically determines element order by examining the node count of each element:

| Element Type | 1st Order Nodes | 2nd Order Nodes | VTK Type (1st/2nd) |
|--------------|-----------------|-----------------|-------------------|
| Tetrahedron | 4 | 10 | 10 / 24 |
| Hexahedron | 8 | 20 | 12 / 25 |
| Wedge/Prism | 6 | 15 | 13 / 26 |
| Pyramid | 5 | 13 | 14 / 27 |
| Triangle | 3 | 6 | 5 / 22 |
| Quadrilateral | 4 | 8 | 9 / 23 |
| Edge | 2 | 3 | 3 / 21 |
| Point | 1 | - | 1 |

### How It Works

1. **Pre-scan**: Before writing, all elements are scanned to determine their node counts
2. **CELLS size**: Calculated dynamically based on actual node counts
3. **CELL_TYPES**: Each element gets the appropriate VTK type based on its node count

This approach allows mixed-order meshes where some elements are 1st order and others are 2nd order.

## Creating 2nd Order Elements in Cubit

To create 2nd order elements in Cubit:

1. **Create mesh first** (1st order by default)
2. **Add elements to block**
3. **Set element type** to convert to 2nd order

```python
# Step 1: Create mesh
cubit.cmd("mesh volume 1")

# Step 2: Add to block
cubit.cmd("block 1 add tet all")

# Step 3: Convert to 2nd order
cubit.cmd("block 1 element type tetra10")
```

### Available 2nd Order Element Types

| 1st Order | 2nd Order | Cubit Command |
|-----------|-----------|---------------|
| TET4 | TET10 | `block X element type tetra10` |
| HEX8 | HEX20 | `block X element type hex20` |
| WEDGE6 | WEDGE15 | `block X element type wedge15` |
| PYRAMID5 | PYRAMID13 | `block X element type pyramid13` |
| TRI3 | TRI6 | `block X element type tri6` |
| QUAD4 | QUAD8 | `block X element type quad8` |
| EDGE2 | EDGE3 | `block X element type bar3` |

## VTK File Format

The exported file follows VTK Legacy ASCII format:

```
# vtk DataFile Version 3.0
Unstructured Grid filename.vtk
ASCII
DATASET UNSTRUCTURED_GRID
POINTS <n> float
<x1> <y1> <z1>
<x2> <y2> <z2>
...
CELLS <num_cells> <total_size>
<num_nodes> <n1> <n2> ...
...
CELL_TYPES <num_cells>
<type1>
<type2>
...
CELL_DATA <num_cells>
SCALARS scalars float
LOOKUP_TABLE default
<value1>
<value2>
...
```

### CELL_DATA Values

The scalar values in CELL_DATA identify element types:

| Value | Element Type |
|-------|--------------|
| 1 | Tetrahedron |
| 2 | Hexahedron |
| 3 | Wedge |
| 4 | Pyramid |
| 5 | Triangle |
| 6 | Quadrilateral |
| 0 | Edge |
| -1 | Point |

## Node Ordering

### 2nd Order Hex (HEX20) - Cubit to VTK

Cubit and VTK use different node ordering for HEX20 elements:

```
Cubit:  [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
VTK:    [0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15]
```

The function handles this conversion automatically.

### 2nd Order Wedge (WEDGE15) - Cubit to VTK

```
Cubit:  [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
VTK:    [0,1,2,3,4,5,6,7,8,12,13,14,9,10,11]
```

## Comparison with Previous Versions

### Before v1.5.1

```python
# User had to specify ORDER parameter
cubit_mesh_export.export_vtk(cubit, "mesh.vtk", ORDER="2nd")
cubit_mesh_export.export_vtk(cubit, "mesh.vtk", ORDER="1st")
```

### v1.5.1 and Later

```python
# Automatic detection - no ORDER parameter needed
cubit_mesh_export.export_vtk(cubit, "mesh.vtk")
```

Benefits of automatic detection:
- Simpler API - no need to specify element order
- Supports mixed-order meshes automatically
- Reduces user errors from incorrect ORDER specification

## Examples

### Complete 3D Example

```python
import os, sys
cubit_path = os.environ.get("CUBIT_PATH")
if cubit_path:
    sys.path.append(cubit_path)

import cubit
cubit.init(['cubit', '-nojournal', '-batch'])

import cubit_mesh_export

# Create geometry
cubit.cmd("reset")
cubit.cmd("create brick x 2 y 2 z 2")

# Mesh with tet elements
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("volume 1 size 0.5")
cubit.cmd("mesh volume 1")

# Define blocks
cubit.cmd("block 1 add tet all in volume 1")
cubit.cmd("block 1 name 'solid'")

cubit.cmd("block 2 add tri all in surface all")
cubit.cmd("block 2 name 'boundary'")

# Export 1st order mesh
cubit_mesh_export.export_vtk(cubit, "brick_1st_order.vtk")

# Convert to 2nd order and export
cubit.cmd("block 1 element type tetra10")
cubit.cmd("block 2 element type tri6")
cubit_mesh_export.export_vtk(cubit, "brick_2nd_order.vtk")

print("Done!")
```

### Visualization in ParaView

1. Open the `.vtk` file in ParaView
2. Apply the filter to visualize
3. Use "Cell Data" scalars to color by element type
4. For 2nd order elements, enable "Nonlinear Subdivision Level" for curved display

## Troubleshooting

### Empty VTK File

**Cause**: No elements added to blocks.

**Solution**: Ensure elements are added to blocks before export:
```python
cubit.cmd("block 1 add tet all in volume 1")
```

### Wrong Element Order

**Cause**: Element type not set after adding to block.

**Solution**: Set element type after adding to block:
```python
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 element type tetra10")  # Must be after adding
```

### Node Indexing Issues

If node indices appear wrong, ensure you're using the latest version. Version 1.5.0 fixed a critical node indexing bug for non-contiguous node IDs.

---

# export_vtu - VTK XML Format

`export_vtu()` exports Cubit mesh data to VTK XML format (.vtu). This is the modern VTK format recommended for ParaView and other visualization tools.

## Key Features

- **XML Format**: Modern, extensible format with better metadata support
- **Automatic Order Detection**: Same as `export_vtk()`
- **BlockID Cell Data**: Each cell includes its source block ID
- **NodeID Point Data**: Original Cubit node IDs are preserved
- **ParaView Recommended**: This is the preferred format for ParaView

## Usage

### Basic Usage

```python
import cubit
import cubit_mesh_export

# Create mesh in Cubit
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 scheme tetmesh")
cubit.cmd("mesh volume 1")
cubit.cmd("block 1 add tet all")

# Export to VTU
cubit_mesh_export.export_vtu(cubit, "mesh.vtu")
```

### With 2nd Order Elements

```python
# Create mesh and convert to 2nd order
cubit.cmd("block 1 add tet all")
cubit.cmd("block 1 element type tetra10")

# Export (auto-detects 2nd order)
cubit_mesh_export.export_vtu(cubit, "mesh_2nd_order.vtu")
```

## Function Signature

```python
def export_vtu(cubit, FileName: str, binary: bool = False):
    """Export mesh to VTK XML format (VTU).

    Args:
        cubit: Cubit Python interface object
        FileName: Output file path for the .vtu file
        binary: If True, write binary data. Default: False (ASCII)

    Returns:
        cubit: The cubit object (for method chaining)
    """
```

## VTU File Structure

```xml
<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="N" NumberOfCells="M">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          ...coordinates...
        </DataArray>
      </Points>
      <Cells>
        <DataArray Name="connectivity" ...>...</DataArray>
        <DataArray Name="offsets" ...>...</DataArray>
        <DataArray Name="types" ...>...</DataArray>
      </Cells>
      <CellData Scalars="BlockID">
        <DataArray Name="BlockID" ...>...</DataArray>
      </CellData>
      <PointData Scalars="NodeID">
        <DataArray Name="NodeID" ...>...</DataArray>
      </PointData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
```

## Additional Data Arrays

### BlockID (Cell Data)

Each cell includes its source Cubit block ID. This allows:
- Material region identification in ParaView
- Post-processing by block
- Color mapping by block

### NodeID (Point Data)

Original Cubit node IDs are preserved as point data:
- Correlation with Cubit model
- Debugging and verification
- Reference to original mesh

### Visualization in ParaView

1. **Open the .vtu file** in ParaView
2. Click **Apply** to load the mesh
3. **Color by BlockID**: Select "BlockID" from the coloring dropdown to visualize different material regions
4. **View NodeID**: Select "NodeID" from the point data dropdown to see original Cubit node numbers
5. **2nd order elements**: Set "Nonlinear Subdivision Level" > 0 in Properties panel to see curved edges

## VTK XML vs Legacy Format Comparison

| Feature | Legacy (.vtk) | XML (.vtu) |
|---------|---------------|------------|
| Format | ASCII text | XML structured |
| Extensibility | Limited | Highly extensible |
| Compression | No | Supported (binary) |
| Metadata | Basic | Rich |
| ParaView | Compatible | Recommended |
| File size | Larger | Smaller (binary) |
| Human readable | Yes | Yes (ASCII mode) |

## When to Use Which

**Use `export_vtk()` when:**
- Maximum compatibility needed
- Simple visualization tasks
- Legacy workflow integration

**Use `export_vtu()` when:**
- Working with ParaView
- Need BlockID/NodeID data
- Modern workflow
- Large meshes (future binary support)

## Example: Multi-Block Export

```python
import cubit
import cubit_mesh_export

cubit.cmd("reset")

# Create two volumes with different mesh types
cubit.cmd("create brick x 1 y 1 z 1")
cubit.cmd("volume 1 move 0 0 0")
cubit.cmd("create sphere radius 0.5")
cubit.cmd("volume 2 move 2 0 0")

cubit.cmd("volume 1 scheme map")
cubit.cmd("volume 2 scheme tetmesh")
cubit.cmd("volume all size 0.2")
cubit.cmd("mesh volume all")

# Assign to different blocks
cubit.cmd("block 1 add hex all")
cubit.cmd("block 1 name 'cube'")
cubit.cmd("block 2 add tet all")
cubit.cmd("block 2 name 'sphere'")

# Export to VTU - BlockID will identify each region
cubit_mesh_export.export_vtu(cubit, "multi_block.vtu")

# In ParaView: Color by "BlockID" to see regions
```

---

## See Also

- [VTK File Formats](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf)
- [ParaView](https://www.paraview.org/)
- [export_Netgen](export_NetgenMesh.md) - For high-order (3rd+) elements via NGSolve

## Version History

| Version | Changes |
|---------|---------|
| 1.6.0 | Added `export_vtu()` for VTK XML format |
| 1.5.1 | Removed ORDER parameter, added automatic order detection |
| 1.5.0 | Fixed node indexing bug for non-contiguous node IDs |
| 1.4.x | Added ORDER parameter for 1st/2nd order selection |
