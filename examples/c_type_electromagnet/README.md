# Electromagnet Simulation - Complete Workflow

Complete 3D magnetostatic simulation of beam steering electromagnet with racetrack coil and magnetic yoke.

## Overview

This directory contains a complete electromagnet simulation workflow:

1. **Mesh Generation**: Generate Netgen mesh directly from Cubit (no intermediate files)
2. **Magnetostatic Simulation**: Solve field distribution with Radia
3. **Visualization**: Export field distribution for ParaView

**Workflow**: Cubit geometry -> export_netgen() -> Netgen mesh -> Radia

## Files

### Main Simulation

- **`main_simulation_workflow.py`** - Main simulation script
  - Loads magnetic yoke via Cubit -> Netgen direct export
  - Creates racetrack coil geometry
  - Solves magnetostatic problem
  - Exports field distribution (VTS)

### Mesh Generation

- **`york_cubit_mesh.py`** - Generate Netgen mesh from Cubit journal file
  - Input: `York.jou` (Cubit journal file)
  - Output: Netgen mesh (in-memory), `York.vtk` (visualization)
  - Creates hexahedral and pentahedral mesh for magnetic yoke

- **`York.jou`** - Cubit journal file with yoke geometry definition

## Complete Workflow

### Step 1: Run Simulation (includes mesh generation)

```bash
cd examples/electromagnet

# Run complete simulation (Cubit -> Netgen -> Radia)
python main_simulation_workflow.py
```

**Output**:
```
======================================================================
ELECTROMAGNET SIMULATION WORKFLOW
======================================================================

[Step 1/5] Importing yoke mesh via Cubit -> Netgen...
  Reading: York.jou
  Exporting to Netgen mesh...
  [OK] Yoke imported: ID=297

[Step 2/5] Creating racetrack coil...
  [OK] Coil created: ID=600
  Total current: -2000 A

[Step 3/5] Combining coil + yoke...
  [OK] Combined model: ID=605

[Step 4/5] Solving magnetostatics...
  [OK] Solution converged

[Step 5/5] Exporting field distribution...
  [OK] Field distribution exported to field_distribution.vts

======================================================================
SIMULATION COMPLETE
======================================================================
```

**Output Files**:
- `field_distribution.vts` - 3D magnetic field distribution (VTS format)

### Step 2: Visualize in ParaView

```bash
# Open magnetic field distribution
paraview field_distribution.vts
```

**In ParaView**:
1. Click "Apply"
2. Add **Glyph** filter:
   - Filters -> Common -> Glyph
   - Glyph Type: Arrow
   - Scalars: None
   - Vectors: B_field
   - Scale Mode: vector
   - Scale Factor: 0.1 (adjust for visibility)
3. Click "Apply" to show field vectors

**Visualization Tips**:
- Use **Slice** filter to view field on cutting planes
- Use **Contour** filter for field magnitude iso-surfaces
- Use **Calculator** to compute |B| magnitude

## Requirements

- **Coreform Cubit 2025.3+** - For hex mesh generation
- **cubit_mesh_export** - From S:\CoreformCubit\01_GitHub
- **NGSolve** - For Netgen mesh wrapper
- **ParaView** - For visualization

## Geometry Specifications

### Racetrack Coil

```python
Center: [0, 131.25, 0] mm
X dimensions: inner=5 mm, outer=40 mm
Y dimensions: inner=50 mm, outer=62.5 mm
Height: 105 mm
Turns: 105
Current: -2000 A
Current density: -0.544218 A/mm^2
Arc approximation: 3 segments
```

### Magnetic Yoke

**Source**: `York.jou` (Cubit journal file)

**Mesh**:
- Format: Netgen (direct from Cubit via export_netgen)
- Elements: 240 hexahedra + 48 pentahedra = 288 total
- Nodes: ~495

**Material**:
- Type: Nonlinear isotropic (saturation model)
- Applied via `rad.MatSatIsoFrm()`

## Troubleshooting

### "Cubit not found"

**Solution**: Install Coreform Cubit or adjust paths:
```python
CUBIT_PATH = "C:/Program Files/Coreform Cubit 2025.3/bin"
CUBIT_EXPORT_PATH = "S:/CoreformCubit/01_GitHub"
```

### "cubit_mesh_export not found"

**Solution**: Clone the Coreform Cubit Mesh Export repository:
```bash
git clone https://github.com/ksugahar/Coreform_Cubit_Mesh_Export S:/CoreformCubit/01_GitHub
```

### Solver returns NaN

**Causes**:
1. Geometry scale mismatch
2. Invalid polyhedra (degenerate elements)
3. Material property errors

**Solution**:
- Verify mesh quality in ParaView: `paraview York.vtk`
- Check coil and yoke bounding boxes overlap correctly

## Coordinate System

- **X**: Horizontal (perpendicular to beam)
- **Y**: Beam direction
- **Z**: Vertical

All dimensions in **millimeters (mm)**.
All magnetic field values in **Tesla (T)**.

## References

- **Radia**: https://github.com/ochubar/Radia
- **Coreform Cubit**: https://coreform.com/products/coreform-cubit/
- **Coreform Cubit Mesh Export**: https://github.com/ksugahar/Coreform_Cubit_Mesh_Export
- **ParaView**: https://www.paraview.org/

---

**Last Updated**: 2026-01-16
**Workflow**: Cubit -> Netgen (direct) -> Radia -> VTS -> ParaView
**Status**: Updated to use Cubit -> Netgen direct export (no Nastran)
