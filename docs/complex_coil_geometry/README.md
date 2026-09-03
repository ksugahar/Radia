# Complex Coil Geometry Docs Showcase

Result-bearing docs showcase demonstrating complex multi-segment coil geometry construction using the modern CoilBuilder API.

## Overview

This directory contains examples of building complex coil paths with straight and arc segments. The examples use the `coil_builder` module which provides a fluent interface for constructing multi-segment coils with automatic state tracking.

## Files

### coil_model.py

**Coil model definition module**

Defines the 8-segment beam steering coil geometry as a reusable module:
- `create_beam_steering_coil()` - Creates the coil and returns (coil_object, parameters)
- Imported by the notebook and the two docs-local inspection helpers

**Usage:**
```python
from coil_model import create_beam_steering_coil

# Create coil
coil, params = create_beam_steering_coil()
```

**Test:**
```bash
cd docs/complex_coil_geometry
python coil_model.py
```

### visualize_coils.py

**Coil geometry visualization and shape verification**

Loads the coil model and visualizes it to verify the geometry is correct:
- Imports coil from coil_model.py (single source of truth)
- Displays coil parameters and bounding box
- Calculates magnetic field at test points for verification
- Reports representative field probes and an axis profile
- Purpose: **Verify coil shape before field calculation**

**Usage:**
```bash
cd docs/complex_coil_geometry
python visualize_coils.py
```

**Requirements:**
- radia (built from this project)
- numpy
- scipy (for rotation matrices)

**Coil specification:**
- Current: 1265 A
- Cross-section: 122×122 mm
- 8 segments total:
  1. Straight: 32.9 mm
  2. Arc: R=121 mm, φ=64.6°, tilt=90°
  3. Straight: 1018.5 mm, tilt=90°
  4. Arc: R=121 mm, φ=115.4°, tilt=-90°
  5. Straight: 906.9 mm, tilt=90°
  6. Arc: R=121 mm, φ=115.4°, tilt=-90°
  7. Straight: 1018.5 mm, tilt=90°
  8. Arc: R=121 mm, φ=64.6°, tilt=-90°

**Output:**
- Console: Magnetic field values at test points
- Console: Field along Z-axis profile

### field_map.py

**3D magnetic field distribution calculation**

Loads the coil model and calculates magnetic field distribution on a 3D grid:
- Imports coil from coil_model.py (same model as visualize_coils.py)
- Fixed evaluation grid sized to the 8-segment coil extent
- Structured grid field calculation
- Vector field export (Bx, By, Bz) in GMSH `.msh v4.1`
- Purpose: **Calculate field distribution for analysis**

**Usage:**
```bash
cd docs/complex_coil_geometry
python field_map.py
```

**Output:**
- File: `field_map.msh` - checked GMSH mesh with field data
- Grid size: Configurable (default: 31×51×31 = 49,011 points)
- Coverage: fixed grid in field_map.py covering the full coil extent
- Field data:
  - **Vector field**: B (Bx, By, Bz) in mT

**Visualization in GMSH:**
- Open `field_map.msh` and select the saved B-field view
- Use vector glyphs, clipping planes, and field ranges as appropriate

## CoilBuilder API

The `coil_builder` module provides a modern, elegant interface for building complex coil geometries.

### Key Features

- **Fluent method chaining** - Readable, declarative code
- **Automatic position/orientation tracking** - No manual state management
- **Type-safe with abstract base classes** - `isinstance()` instead of string comparison
- **~75% less code** than manual segment tracking
- **Automatic cross-section transformation with tilt**
- **Direct conversion to Radia objects**

**Location:** `src/radia/coil_builder.py`

### Basic Example

```python
from radia.coil_builder import CoilBuilder
import radia as rad

# Radia always uses meters
mm = 1e-3  # unit conversion factor

# Simple racetrack coil
coil_segments = (CoilBuilder(current=1000)
	.set_start([0, -50*mm, 0])
	.set_cross_section(width=20*mm, height=20*mm)
	.add_straight(length=100*mm)
	.add_arc(radius=50*mm, arc_angle=180, tilt=90)
	.add_straight(length=100*mm)
	.add_arc(radius=50*mm, arc_angle=180, tilt=90)
	.to_radia())

# Combine into Radia container
coils = rad.ObjCnt(coil_segments)
```

### Complex Example (8-Segment Beam Steering Magnet)

```python
# Radia always uses meters
mm = 1e-3

coil = (CoilBuilder(current=1265)
	.set_start([218*mm, -16.4*mm, -81*mm])
	.set_cross_section(122*mm, 122*mm)
	.add_straight(32.9*mm)
	.add_arc(121*mm, 64.6, tilt=90)
	.add_straight(1018.5*mm, tilt=90)
	.add_arc(121*mm, 115.4, tilt=-90)
	.add_straight(906.9*mm, tilt=90)
	.add_arc(121*mm, 115.4, tilt=-90)
	.add_straight(1018.5*mm, tilt=90)
	.add_arc(121*mm, 64.6, tilt=-90)
	.to_radia())
```

For the complete example with field calculation and visualization, see `visualize_coils.py`.

## Design Overview

### Architecture

The CoilBuilder uses modern design patterns for clean, maintainable code:

```
CoilSegment (ABC)
├── StraightSegment
│   ├── Geometry: rectangular conductor
│   ├── end_pos: linear progression
│   └── to_radia(): ObjRecCur + transforms
│
└── ArcSegment
	├── Geometry: toroidal sector
	├── end_pos: rotation around center
	└── to_radia(): ObjArcCur + transforms

CoilBuilder
├── State management: position, orientation, cross-section
├── Fluent methods: add_straight(), add_arc()
└── Conversion: to_radia()
```

### Design Patterns Used

1. **Builder Pattern** - CoilBuilder fluent interface
2. **Strategy Pattern** - CoilSegment hierarchy
3. **Template Method** - Abstract base class
4. **Property-Based Design** - `@property` decorators for computed values

### Code Comparison

**Old approach (manual tracking):**
```python
COIL = []
x0 = array([0, 0, 0])
V = eye(3).T

# Segment 1
COIL.append(cCOIL('GCE', I, x0, V, W, H, [], 0, 0, L))

# Must manually update for next segment
x0 = COIL[-1].x1
V = COIL[-1].V1
COIL.append(cCOIL('ARC', I, x0, V, W, H, R, phi, 90, 0))

# ... repeat for each segment ...
```

**New approach (automatic tracking):**
```python
mm = 1e-3
segments = (CoilBuilder(current=1000)
	.set_start([0, 0, 0])
	.set_cross_section(20*mm, 20*mm)
	.add_straight(100*mm)
	.add_arc(50*mm, 180, tilt=90)
	.to_radia())
```

**Benefits:**
- **75% less code** - From ~60 lines to ~15 lines for 8 segments
- **Zero manual tracking** - Automatic state management
- **Type-safe** - No string comparison
- **Readable** - Reads like a description
- **Hard to make mistakes** - No forgotten state updates

## Critical Bug Fixes

### Orientation Matrix Convention

The implementation uses **row vector** format for orientation matrices:

```python
# Each ROW is a basis vector
orientation[0, :]  # X-axis
orientation[1, :]  # Y-axis
orientation[2, :]  # Z-axis

# Arc center calculation
self.arc_center = self.start_pos - self.radius * self.orientation[0, :]

# Straight segment end position
self.end_pos = self.start_pos + self.length * self.orientation[1, :]
```

This matches the original COIL.py convention and ensures correct geometric transformations, especially for tilted segments.

### Tilt Transformation

Tilt is applied **only once** in the subclass constructors (`StraightSegment`, `ArcSegment`), not in the base class. This prevents double-application of tilt transformations.

## Quick Start

```bash
cd docs/complex_coil_geometry
python visualize_coils.py
```

This creates an 8-segment coil and calculates magnetic field at multiple points.

**What it demonstrates:**
- Complex multi-segment coil construction
- Fluent CoilBuilder API
- Automatic state tracking
- Magnetic field calculation

Run `python field_map.py` when a checked GMSH field artifact is also needed.

## Field Calculations

The examples demonstrate field calculation at multiple points:

```python
import radia as rad

# Radia always uses meters
mm = 1e-3

# Single point
B = rad.Fld(obj, 'b', [0, 0, 100*mm])  # Returns [Bx, By, Bz] in Tesla

# Multiple points
for z_mm in [0, 100, 500]:
	B = rad.Fld(obj, 'b', [0, 0, z_mm * mm])
	print(f"At z={z_mm} mm: B={B}")
```

**Example output:**
```
Position (mm)             Bx (mT)         By (mT)         Bz (mT)
----------------------------------------------------------------------
(0, 0, 0)                 1.580230        0.000000        0.000000
(0, 0, 100)               2.161479        0.000000        0.000000
(0, 0, 500)               2.443147        0.000000        0.000000
```

## Visualization Options

### GMSH Field Export

```bash
python field_map.py
```

## Coordinate System

- X: Horizontal (left-right)
- Y: Horizontal (front-back)
- Z: Vertical (up)

All dimensions in meters (Radia always uses meters). Design values use `mm = 1e-3` conversion.
All fields in Tesla (T).

## Troubleshooting

### "No module named 'radia'"

**Solution:** Install Radia in the active Python environment:
```bash
python -m pip install radia
```

### "No module named 'coil_builder'"

**Solution:** Use the public package import `from radia.coil_builder import CoilBuilder`.
For repository development, install the checkout in editable mode instead of
adding `src/radia` to `sys.path`.

### "No module named 'scipy'"

**Solution:** Install scipy:
```bash
pip install scipy
```

## Further Reading

- [src/radia/README.md](../../src/radia/README.md) - Visualization utilities and CoilBuilder documentation
- [tests/README.md](../../tests/README.md) - Test suite
- [Build from source](../../README.md#build-from-source) - Build instructions

---

**Last Updated**: 2026-09-03
