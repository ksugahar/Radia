# VTK Export

Demonstrates how to export Radia magnetic field data to VTK Structured Grid (VTS) format for 3D visualization in ParaView.

## Scripts

| File | Description |
|------|-------------|
| `demo_vts_export.py` | Creates a permanent magnet and exports B and H fields on a 3D grid to VTS format using `rad.FldVTS()` |

## Usage

```bash
python demo_vts_export.py
```

This generates `magnet_field.vts`, which can be opened in ParaView. Use the Glyph filter for vector arrows and color by `B_magnitude` to visualize the field.
