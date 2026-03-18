# VTK Export Examples

Export Cubit mesh to VTK format (.vtk, .vtu).

## VTK vs VTU Format

| Feature | VTK Legacy (.vtk) | VTK XML (.vtu) |
|---------|-------------------|----------------|
| Format | ASCII text | XML structured |
| Function | `export_vtk()` | `export_vtu()` |
| Readability | Human-readable | Human-readable (ASCII mode) |
| Extensibility | Limited | Highly extensible |
| Metadata | Basic | Rich (BlockID, NodeID) |
| File size | Larger | Smaller (binary option) |
| ParaView | Compatible | **Recommended** |

### When to Use Which

**Use VTK Legacy (.vtk) when:**
- Maximum compatibility with older software
- Simple visualization tasks
- Legacy workflow integration

**Use VTK XML (.vtu) when:**
- Working with ParaView (recommended format)
- Need BlockID/NodeID data for post-processing
- Modern workflow

## Usage

```python
import cubit
import cubit_mesh_export

cubit.cmd("block 1 add tet all")
cubit.cmd("block 2 add tri all")

# VTK Legacy format
cubit_mesh_export.export_vtk(cubit, "mesh.vtk")

# VTK XML format (recommended for ParaView)
cubit_mesh_export.export_vtu(cubit, "mesh.vtu")
```

## Additional Data in VTU Format

The VTU format includes extra data arrays:

| Data | Type | Description |
|------|------|-------------|
| **BlockID** | Cell Data | Source Cubit block ID for each element |
| **NodeID** | Point Data | Original Cubit node ID |

This enables:
- Material region identification
- Color mapping by block in ParaView
- Correlation with original Cubit model

## Automatic Order Detection

Both functions auto-detect element order (1st/2nd) from node count:

```python
# 1st order (auto-detected)
cubit.cmd("block 1 add tet all")
cubit_mesh_export.export_vtk(cubit, "mesh_1st.vtk")

# 2nd order (auto-detected)
cubit.cmd("block 1 element type tetra10")
cubit_mesh_export.export_vtk(cubit, "mesh_2nd.vtk")
```

## Sample Files

| File | Description |
|------|-------------|
| `sphere_1st_order.vtk/.vtu` | 1st order tet mesh |
| `sphere_2nd_order.vtk/.vtu` | 2nd order tet mesh (TET10) |
| `mixed_elements.vtk/.vtu` | Mixed hex + tet mesh |
| `circle_2d.vtk/.vtu` | 2D triangular mesh |

## Visualization in ParaView

1. File > Open > select .vtk or .vtu file
2. Click "Apply" in Properties panel
3. For VTU: Color by "BlockID" to visualize material regions
4. For 2nd order: Set "Nonlinear Subdivision Level" > 0 for curved display

## Regenerate Samples

```bash
"${CUBIT_PATH:-C:/Program Files/Coreform Cubit 2025.3/bin}/python3/python.exe" vtk_vtu_export_example.py
```

## See Also

- [docs/cubit/export_VTK.md](../../../docs/cubit/export_VTK.md) - Full documentation
