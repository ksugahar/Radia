# Electromagnet Simulation - Complete Workflow

Complete 3D magnetostatic simulation of beam steering electromagnet with racetrack coil and magnetic yoke.

## Overview

This directory contains a complete electromagnet simulation workflow:

1. **Mesh Generation**: Generate Nastran mesh from Cubit journal file
2. **Magnetostatic Simulation**: Solve field distribution with Radia
3. **Visualization**: Export geometry and fields for ParaView

## Files

### Main Simulation

- **`main_simulation_workflow.py`** - Main simulation script
  - Creates racetrack coil geometry
  - Loads magnetic yoke from Nastran mesh
  - Solves magnetostatic problem
  - Exports geometry (VTK) and field distribution (VTK)

### Mesh Generation

- **`york_cubit_mesh.py`** - Generate York.bdf from Cubit journal file
  - Input: `york.jou` (Cubit journal file)
  - Output: `York.bdf` (Nastran format), `York.vtk` (visualization)
  - Creates hexahedral and pentahedral mesh for magnetic yoke

- **`York.jou`** - Cubit journal file with yoke geometry definition

### Utilities (in `src/python/`)

- **`nastran_mesh_import.py`** - Nastran mesh import and conversion to Radia
- **`radia_vtk_export.py`** - VTK Legacy format export

## Complete Workflow

### Step 1: Generate Mesh from Cubit Journal

```bash
cd examples/electromagnet

# Generate Nastran mesh from Cubit journal
python york_cubit_mesh.py
```

**Output**:
- `York.bdf` - Nastran mesh for Radia simulation (569 vertices, 288 elements)
- `York.vtk` - VTK mesh for ParaView visualization
- `York.msh` - Gmsh format (optional)

**Requirements**:
- Coreform Cubit 2025.3 (or compatible version)
- `york.jou` - Cubit journal file with geometry definition

### Step 2: Run Magnetostatic Simulation

```bash
# Run complete simulation
python main_simulation_workflow.py
```

**Output**:
```
======================================================================
ELECTROMAGNET SIMULATION WORKFLOW
======================================================================

[Step 1/5] Importing yoke mesh from York.bdf...
[Nastran Import] Reading file: York.bdf
                 Vertices: 569
                 Hexahedral elements: 288
                 Wedge elements: 8
  [OK] Yoke imported: ID=297

[Step 2/5] Creating racetrack coil...
  [OK] Coil created: ID=600
  Total current: -2000 A

[Step 3/5] Combining coil + yoke...
  [OK] Combined model: ID=605

[Step 4/5] Exporting Radia_model.vtk...
  Polygons: 1848
  Points: 7376

[Step 5/5] Calculating magnetic field distribution...
  Grid: 21x31x21 = 13671 points
  Solving magnetostatics...
  [OK] Solution converged
  Calculating magnetic field...
  [OK] Field calculated at 13671 points
  [OK] Exported: field_distribution.vtk

======================================================================
SIMULATION COMPLETE
======================================================================
```

**Output Files**:
- `Radia_model.vtk` - Combined geometry (coil + yoke) in VTK Legacy format
- `field_distribution.vtk` - 3D magnetic field distribution with vectors

### Step 3: Visualize in ParaView

#### Method 1: Open Combined Geometry

```bash
# Open combined geometry
paraview Radia_model.vtk
```

**In ParaView**:
1. Click "Apply"
2. Adjust colors using "Radia_colours" field
3. Rotate view to inspect coil and yoke geometry

#### Method 2: Open Field Distribution

```bash
# Open magnetic field distribution
paraview field_distribution.vtk
```

**In ParaView**:
1. Click "Apply"
2. Add **Glyph** filter:
   - Filters → Common → Glyph
   - Glyph Type: Arrow
   - Scalars: None
   - Vectors: B_field
   - Scale Mode: vector
   - Scale Factor: 0.1 (adjust for visibility)
3. Click "Apply" to show field vectors

**Visualization Tips**:
- Use **Slice** filter to view field on cutting planes
- Use **Contour** filter for field magnitude iso-surfaces
- Use **Calculator** to compute |B| magnitude: `sqrt(B_field_X^2 + B_field_Y^2 + B_field_Z^2)`

## Geometry Specifications

### Racetrack Coil

```python
# From racetrack_coil_model.py
Center: [0, 131.25, 0] mm
X dimensions: inner=5 mm, outer=40 mm
Y dimensions: inner=50 mm, outer=62.5 mm
Height: 105 mm
Turns: 105
Current: -2000 A
Current density: -0.544218 A/mm^2
Arc approximation: 3 segments
```

**Bounding box**: X[-65, 65], Y[60, 202.5], Z[-52.5, 52.5] mm

### Magnetic Yoke

**Source**: `york.jou` (Cubit journal file)

**Mesh**:
- Format: Nastran bulk data (.bdf)
- Elements: 240 hexahedra + 48 pentahedra = 288 total
- Nodes: ~495

**Material**:
- Type: Linear isotropic
- Relative permeability: μr = 1000
- No remanent magnetization

## Solver Configuration

```python
# Magnetostatic solver settings
Precision: 0.01
Max iterations: 1000
Method: 4 (relaxation)
```

**Typical convergence**: ~20-30 iterations

## Field Calculation

**Field points**: Three positions along Z-axis
- Origin: [0, 0, 0]
- Z=100mm: [0, 0, 100]
- Z=500mm: [0, 0, 500]

**Field distribution grid**:
- Range: Geometry bbox + 50mm margin
- Resolution: 21 × 31 × 21 = 13,671 points
- Format: VTK STRUCTURED_POINTS with vector data

## File Formats

### VTK Legacy (ASCII)

**Current format** - Used for all outputs

**Advantages**:
- Human-readable ASCII
- Compatible with all VTK tools and ParaView versions
- Single combined file for geometry
- Point data format for field distribution

**Files**:
- `Radia_model.vtk` - Combined coil + yoke geometry (POLYDATA)
- `field_distribution.vtk` - Structured point cloud with vector field data (POLYDATA with VECTORS)

## Troubleshooting

### "York.bdf not found"

**Solution**: Run mesh generation first:
```bash
python york_cubit_mesh.py
```

This will generate `York.bdf` from `york.jou` using Cubit.

### "Cubit not found" (york_cubit_mesh.py)

**Solution**: Install Coreform Cubit or adjust path in `york_cubit_mesh.py`:
```python
sys.path.append("C:/Program Files/Coreform Cubit 2025.3/bin")
```

**Alternative**: Use pre-generated `York.bdf` (already provided in repository)

### Solver returns NaN

**Causes**:
1. Geometry scale mismatch
2. Invalid polyhedra (degenerate elements)
3. Material property errors

**Solution**:
- Verify mesh quality in ParaView: `paraview York_mesh.vtk`
- Check coil and yoke bounding boxes overlap correctly
- Verify material properties in `yoke_model.py`

### VTK export fails

**Error**: `ModuleNotFoundError: No module named 'radia_vtk_export'`

**Solution**: Ensure `src/python/radia_vtk_export.py` exists and path is correct:
```python
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))
```

### Field distribution file is huge

**Issue**: 13,671 points × 3 components = ~40KB (normal size)

If file is >10MB, grid resolution may be too high. Adjust in `main_simulation_workflow.py`:
```python
x_range = np.linspace(-100, 100, 21)  # Reduce number of points
y_range = np.linspace(0, 250, 31)
z_range = np.linspace(-100, 100, 21)
```

## Coordinate System

- **X**: Horizontal (perpendicular to beam)
- **Y**: Beam direction
- **Z**: Vertical

All dimensions in **millimeters (mm)**.

All magnetic field values in **Tesla (T)**.

## Performance Notes

**VTK export timing** (typical):
- Radia_model.vtk: ~1-2 seconds (combined coil + yoke, 1848 polygons)

**Field calculation timing**:
- 13,671 points: ~5-10 seconds (depending on CPU)

## Further Reading

- [Radia Python API](../../README.md)
- [VTK Export Utilities](../../src/python/radia_vtk_export.py)
- [Nastran Mesh Import](../../src/python/nastran_mesh_import.py)

## References

- **Radia**: https://github.com/ochubar/Radia
- **Coreform Cubit**: https://coreform.com/products/coreform-cubit/
- **ParaView**: https://www.paraview.org/
- **Nastran**: MSC Nastran Bulk Data format

---

**Last Updated**: 2025-11-26
**Workflow**: Cubit → Nastran → Radia → VTK → ParaView
**Status**: ✓ Working - All bugs fixed (face duplication resolved)
