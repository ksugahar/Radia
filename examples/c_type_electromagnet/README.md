# Electromagnet Examples

Racetrack electromagnet simulation with magnetic yoke from Nastran mesh.

## Overview

This directory contains a complete electromagnet simulation combining:
- Racetrack coil geometry (using `rad.ObjRaceTrk`)
- Magnetic yoke imported from Nastran .nas mesh file
- Magnetostatic solver
- VTK export for visualization

## Files

### magnet.py

**Complete electromagnet simulation**

Production-ready racetrack electromagnet simulation featuring:
- Racetrack coil matching York.nas yoke geometry
- Magnetic yoke imported from Nastran mesh file
- Automatic VTK export
- Magnetostatic solver for nonlinear magnetic analysis
- Field calculation at multiple points

**Usage:**
```bash
cd examples/electromagnet
python magnet.py
```

**Requirements:**
- radia (built from this project)
- numpy

**Features:**
- Works with or without magnetic yoke
- Automatic VTK export to `electromagnet.vtk`
- Solver convergence checking
- Direct Nastran mesh reading (no Cubit required)

**Output:**
```
======================================================================
RACETRACK ELECTROMAGNET SIMULATION
======================================================================

Creating racetrack coil...
  Current: -2000 A
  Turns: 105
  Current density: -0.544218 A/mm^2
  [OK] Racetrack coil created

Creating magnetic yoke from Nastran mesh...
  Nodes: 495
  Elements (CHEXA): 240
  Elements (CPENTA): 48
  Total elements: 288
  [OK] Created 288 polyhedra (240 hex + 48 penta)

Solving magnetostatics...
  Solver result: [0.00954..., 1.206..., 0.00219..., 23.0]
  [OK] Solver completed (iterations: 23.0)

Calculating magnetic field...
Position (mm)        Bx (mT)         By (mT)         Bz (mT)         |B| (mT)
(0, 0, 0)            -0.438314       -0.239402       205.171258      205.171866
(0, 0, 100)          -35.995225      9.777235        0.362883        37.301236
(0, 0, 500)          -0.000594       0.046554        -0.214708       0.219698

Exporting geometry to VTK...
  [OK] Created: electromagnet.vtk

Calculating field distribution...
  Grid resolution: 21 × 31 × 21 = 13671 points
  Field calculation range: bbox + 50mm
  [OK] Created: field_distribution.vtk
```

### racetrack_coil_model.py

**Racetrack coil geometry**

Creates a racetrack coil matching the York.nas yoke geometry.

**Function:**
```python
def create_racetrack_coil(current=-2000):
    """
    Create racetrack coil for electromagnet.

    Original Mathematica specification:
    coil = radObjRaceTrk[{0, 131.25, 0}, {5, 40}, {50, 62.5}, 105, 3, current/105/35]

    Args:
        current: Total current (A), default -2000 A

    Returns:
        tuple: (coil_object, coil_parameters)
    """
```

**Coil parameters:**
- Center: [0, 131.25, 0] mm
- X dimensions: inner=5 mm, outer=40 mm
- Y dimensions: inner=50 mm, outer=62.5 mm
- Height: 105 mm
- Turns: 105
- Current: -2000 A
- Arc segments: 3

**Test:**
```bash
python racetrack_coil_model.py
```

### yoke_model.py

**Magnetic yoke model from Nastran mesh**

Imports magnetic yoke geometry from Nastran .nas file.

**Function:**
```python
def create_yoke_from_nastran(nas_file, material_properties=None):
    """
    Create magnetic yoke from Nastran mesh file.

    Args:
        nas_file: Path to .nas file
        material_properties: Dictionary with material properties
            - permeability: [mu_x, mu_y] or scalar
            - type: 'linear' (default)

    Returns:
        tuple: (yoke_object, mesh_info)
    """
```

**Features:**
- Reads Nastran bulk data format (.nas/.bdf)
- Parses GRID (nodes), CHEXA (hexahedra), and CPENTA (pentahedra) entries
- Converts to Radia polyhedra
- Applies linear magnetic material (μr=1000)
- No symmetries applied (compatible with arbitrary coil geometry)

**Test:**
```bash
python yoke_model.py
```

### nastran_reader.py

**Nastran mesh file parser**

Low-level parser for Nastran bulk data format.

**Function:**
```python
def read_nastran_mesh(filename):
    """
    Read Nastran mesh file (.nas format).

    Returns:
        dict: Dictionary with 'nodes' and 'elements'
            - nodes: numpy array (N, 3) with node coordinates
            - elements: numpy array (M, 8) with element connectivity (1-indexed)
            - node_ids: dict mapping node_id to array index
    """
```

**Supported entries:**
- GRID: Node coordinates (fixed 8-character fields)
- CHEXA: 8-node hexahedron elements (with continuation line)
- CPENTA: 6-node pentahedron elements (single line)

**Test:**
```bash
python nastran_reader.py
```

### York.bdf

**Nastran mesh file for magnetic yoke**

- Format: Nastran bulk data
- Nodes: 495 (GRID entries)
- Elements: 288 total
  - 240 hexahedra (CHEXA entries)
  - 48 pentahedra (CPENTA entries)
- Bounding box: X[-25, 25], Y[-20, 162.5], Z[-105, 105] mm
- Designed for racetrack coil geometry

## Operating Modes

### Mode 1: Coils Only (No Yoke)

If York.nas is not found, the simulation runs with coils only:

```bash
python magnet.py
```

**Output:**
- Racetrack coil geometry created
- Field calculated (current source only)
- VTK file exported
- No magnetostatic solver needed

### Mode 2: Full Simulation (With Yoke)

If York.nas is available:

**Process:**
1. Creates racetrack coil
2. Imports magnetic yoke from Nastran mesh
3. Combines geometries
4. Runs magnetostatic solver
5. Exports complete geometry to VTK

**Solver:**
- Precision: 0.01
- Max iterations: 1000
- Method: Relaxation (method 4)
- Typical convergence: ~30 iterations

## Coil Specification

**Racetrack coil parameters:**
- Type: Racetrack (rectangular loop with rounded ends)
- Current: -2000 A
- Turns: 105
- Cross-sectional area: 35 mm^2
- Current density: -0.544218 A/mm^2

**Geometry:**
- Center: [0, 131.25, 0] mm
- Inner semi-width X: 5 mm
- Outer semi-width X: 40 mm
- Inner semi-length Y: 50 mm
- Outer semi-length Y: 62.5 mm
- Height Z: 105 mm
- Arc approximation: 3 segments

**Bounding box:**
- X: [-65, 65] mm
- Y: [60, 202.5] mm
- Z: [-52.5, 52.5] mm

## Visualization

### ParaView - Geometry

Open the geometry VTK file in ParaView:

```bash
# After running magnet.py
# Open electromagnet.vtk in ParaView
```

**File info:**
- Polygons: ~1760 (coil 80 + yoke 1680)
- Points: ~6944
- Format: VTK Legacy ASCII

### ParaView - Magnetic Field Distribution

Open the field distribution VTK file in ParaView:

```bash
# Open field_distribution.vtk in ParaView
```

**File info:**
- Grid: 21 × 31 × 21 = 13,671 points
- Range: Geometry bbox + 50mm margin
- Data: Magnetic field vectors (Bx, By, Bz) in Tesla
- Format: VTK STRUCTURED_POINTS

**Visualization steps in ParaView:**
1. Open `field_distribution.vtk`
2. Apply **Glyph** filter:
   - Glyph Type: Arrow
   - Scalars: B_field (magnitude)
   - Vectors: B_field
   - Scale Mode: vector
   - Scale Factor: adjust for visibility
3. Apply **Contour** filter for field magnitude iso-surfaces
4. Use **Slice** filter to view field on cutting planes

**Field calculation range:**
- X: [-115, 115] mm (geometry: [-65, 65] mm + 50mm margin)
- Y: [-70, 252.5] mm (geometry: [-20, 202.5] mm + 50mm margin)
- Z: [-155, 155] mm (geometry: [-105, 105] mm + 50mm margin)

### Field Point Values

The script also calculates field at specific points:
- Origin: [0, 0, 0] → ~205 mT (with yoke)
- Z=100mm: [0, 0, 100] → ~37 mT
- Z=500mm: [0, 0, 500] → ~0.22 mT

All field values in milli-Tesla (mT).

## Coordinate System

- X: Horizontal (perpendicular to beam)
- Y: Beam direction
- Z: Vertical

All dimensions in millimeters (mm).

## Nastran Format Details

### GRID Entry (Nodes)

Fixed-width format (8 characters per field):

```
GRID    ID      CP      X1      X2      X3
GRID    1               0.0     0.0     0.0
```

Fields:
- Positions 0-7: "GRID"
- Positions 8-15: Node ID
- Positions 24-31: X coordinate
- Positions 32-39: Y coordinate
- Positions 40-47: Z coordinate

### CHEXA Entry (Hexahedra)

Fixed-width format with continuation:

```
CHEXA   EID     PID     G1      G2      G3      G4      G5      G6      +
+       G7      G8
```

Fields:
- First line: Element ID, Property ID, Nodes 1-6
- Continuation line: Nodes 7-8

Node ordering matches standard CHEXA convention.

### CPENTA Entry (Pentahedra)

Fixed-width format (single line):

```
CPENTA  EID     PID     G1      G2      G3      G4      G5      G6
```

Fields:
- All on one line: Element ID, Property ID, Nodes 1-6
- Nodes 1-3: Bottom triangle
- Nodes 4-6: Top triangle

Node ordering matches standard CPENTA convention.

## Troubleshooting

### "File not found: York.bdf"

**Solution:** This is normal if the yoke mesh file is not available. The simulation will run in coils-only mode.

To enable full simulation:
1. Ensure `York.bdf` is in the `examples/electromagnet/` directory
2. Verify the file contains valid Nastran bulk data (GRID, CHEXA, and CPENTA entries)

### Solver returns NaN

**Causes:**
1. Geometry scale mismatch between coil and yoke
2. Symmetries incompatible with coil geometry
3. Invalid polyhedra (degenerate elements)

**Solution:**
- Verify coil and yoke bounding boxes are compatible
- Disable symmetries if needed (see yoke_model.py:135-137)
- Check mesh quality in Nastran file

### OpenGL viewer fails

**Solution:** This is expected on some systems. Use the VTK file with ParaView instead:

```bash
# Open electromagnet.vtk in ParaView
```

### Unicode encoding errors (Windows)

**Issue:** Console cannot display superscript characters (², ³, etc.)

**Solution:** Code uses ASCII alternatives (mm^2, mm^3) for console output.

## Design Notes

### Nastran vs Cubit

This implementation uses **direct Nastran mesh reading** instead of Cubit Python API:

**Advantages:**
- No Cubit installation required
- Faster startup (no Cubit initialization)
- Platform independent
- Easier to distribute

**Previous approach:**
- Required Cubit 2023.11 installed
- Used Cubit Python API to read .jou journal files
- Generated intermediate .wls Wolfram Language files

### Symmetry Considerations

The yoke model originally applied symmetries:
- Y-Z plane symmetry: `rad.TrfZerPerp(yoke, [0, 0, 0], [0, 1, 0])`
- X=0 plane symmetry: `rad.TrfZerPara(yoke, [0, 0, 0], [1, 0, 0])`

These are **disabled by default** (see yoke_model.py:135-137) to ensure compatibility with arbitrary coil geometries. Enable only if your geometry truly has these symmetries.

## Further Reading

- [src/python/README.md](../../src/python/README.md) - Radia Python utilities
- [examples/complex_coil_geometry/](../complex_coil_geometry/) - CoilBuilder examples
- [README_BUILD.md](../../README_BUILD.md) - Build instructions

## References

- **Nastran format**: MSC Nastran Bulk Data specification
- **Radia**: https://github.com/ochubar/Radia
- **ParaView**: https://www.paraview.org/

---

**Last Updated**: 2025-10-30
