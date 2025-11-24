# Electromagnet Simulation - Complete Workflow

Complete 3D magnetostatic simulation of beam steering electromagnet with racetrack coil and magnetic yoke.

## Overview

This directory contains a complete electromagnet simulation workflow:

1. **Mesh Generation**: Generate Nastran mesh from Cubit journal file
2. **Magnetostatic Simulation**: Solve field distribution with Radia
3. **Visualization**: Export geometry and fields for ParaView

## Files

### Main Simulation

- **`main_electromagnet_simulation.py`** - Main simulation script
  - Creates racetrack coil geometry
  - Loads magnetic yoke from Nastran mesh
  - Solves magnetostatic problem
  - Exports geometry (VTU) and field distribution (VTK)

### Mesh Generation

- **`York_cubit_mesh.py`** - Generate York.bdf from Cubit journal file
  - Input: `york.jou` (Cubit journal file)
  - Output: `York.bdf` (Nastran format), `York.vtk` (visualization)
  - Creates hexahedral and pentahedral mesh for magnetic yoke

- **`york.jou`** - Cubit journal file with yoke geometry definition

### Component Models

- **`racetrack_coil_model.py`** - Racetrack coil geometry
  - `create_racetrack_coil()`: Create coil with specified current
  - `get_coil_info()`: Get bounding box information

- **`yoke_model.py`** - Magnetic yoke from Nastran mesh
  - `create_yoke_from_nastran()`: Load and convert Nastran mesh to Radia

### Utilities

- **`nastran_reader.py`** - Low-level Nastran file parser (in `src/python/`)
- **`radia_vtk_export.py`** - VTK Legacy format export (in `src/python/`)
- **`radia_vtu_export.py`** - VTU (VTK XML) format export (in `src/python/`)

## Complete Workflow

### Step 1: Generate Mesh from Cubit Journal

```bash
cd examples/electromagnet

# Generate Nastran mesh from Cubit journal
python York_cubit_mesh.py
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
python main_electromagnet_simulation.py
```

**Output**:
```
======================================================================
MAIN ELECTROMAGNET SIMULATION
======================================================================

Creating racetrack coil...
  Current: -2000 A
  Turns: 105
  Current density: -0.544218 A/mm^2
  [OK] Coil created

Creating magnetic yoke from Nastran mesh...
  [OK] Created 288 polyhedra (240 hex + 48 penta)

Solving magnetostatics...
  [OK] Solver completed (iterations: 21)

Calculating magnetic field...
Position (mm)        Bx (mT)         By (mT)         Bz (mT)         |B| (mT)
(0, 0, 0)            -8579.596      0.941           3351.036        9210.804

Exporting geometries to VTU...
  [OK] Created: coil_geometry.vtu (0.15s)
  [OK] Created: yoke_geometry.vtu (1.23s)
  [OK] Created: main_electromagnet_simulation.vtk (combined, 1.45s)

Calculating field distribution...
  Grid resolution: 21 × 31 × 21 = 13671 points
  [OK] Created: field_distribution.vtk
======================================================================
```

**Output Files**:
- `coil_geometry.vtu` - Racetrack coil geometry (VTU format)
- `yoke_geometry.vtu` - Magnetic yoke geometry (VTU format)
- `main_electromagnet_simulation.vtk` - Combined geometry (VTK Legacy)
- `field_distribution.vtk` - 3D magnetic field distribution

### Step 3: Visualize in ParaView

#### Method 1: Open Individual Geometry Files

```bash
# Open coil and yoke separately
paraview coil_geometry.vtu yoke_geometry.vtu
```

**In ParaView**:
1. Click "Apply" for both geometries
2. Adjust colors using "Radia_colours" field
3. Rotate view to inspect geometry

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

### VTU (VTK XML Unstructured Grid)

**Modern format** - Recommended for ParaView

**Advantages**:
- XML-based, robust parsing
- Better support in modern ParaView
- Separate files for coil and yoke

**Files**:
- `coil_geometry.vtu`
- `yoke_geometry.vtu`

### VTK Legacy (ASCII)

**Legacy format** - For compatibility

**Advantages**:
- Human-readable ASCII
- Compatible with older VTK tools
- Single combined file

**Files**:
- `main_electromagnet_simulation.vtk` (combined coil + yoke)
- `field_distribution.vtk` (structured grid with vectors)

## Troubleshooting

### "York.bdf not found"

**Solution**: Run mesh generation first:
```bash
python York_cubit_mesh.py
```

This will generate `York.bdf` from `york.jou` using Cubit.

### "Cubit not found" (York_cubit_mesh.py)

**Solution**: Install Coreform Cubit or adjust path in `York_cubit_mesh.py`:
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

### VTU export fails

**Error**: `ModuleNotFoundError: No module named 'radia_vtu_export'`

**Solution**: Ensure `src/python/radia_vtu_export.py` exists and path is correct:
```python
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))
```

### Field distribution file is huge

**Issue**: 13,671 points × 3 components = ~40KB (normal size)

If file is >10MB, grid resolution may be too high. Adjust in `main_electromagnet_simulation.py`:
```python
nx, ny, nz = 21, 31, 21  # Reduce these numbers
```

## Coordinate System

- **X**: Horizontal (perpendicular to beam)
- **Y**: Beam direction
- **Z**: Vertical

All dimensions in **millimeters (mm)**.

All magnetic field values in **Tesla (T)**.

## Performance Notes

**VTU export timing** (typical):
- Coil geometry: ~0.1-0.2 seconds
- Yoke geometry: ~1-2 seconds (288 elements)
- Combined VTK: ~1.5-3 seconds

**Field calculation timing**:
- 13,671 points: ~5-10 seconds (depending on CPU)

## Further Reading

- [Radia Python API](../../README.md)
- [VTK Export Utilities](../../src/python/README.md)
- [Nastran Format Details](./nastran_reader.py)

## References

- **Radia**: https://github.com/ochubar/Radia
- **Coreform Cubit**: https://coreform.com/products/coreform-cubit/
- **ParaView**: https://www.paraview.org/
- **Nastran**: MSC Nastran Bulk Data format

---

**Last Updated**: 2025-11-22
**Workflow**: Cubit → Nastran → Radia → VTU/VTK → ParaView
