# Samarium-Cobalt (SmCo) Magnet Array

Hexagonal array of cylindrical SmCo permanent magnets.

**Ported from**: `2023_10_01_サマコバ/magnet.nb` (Mathematica notebook)

## Overview

This example simulates a hexagonal close-packed array of cylindrical Samarium-Cobalt (SmCo) permanent magnets mounted on a base plate. The configuration is commonly used for:
- Magnetic field generation
- Undulator/wiggler magnet arrays
- Multipole magnet designs
- Halbach array configurations

## Features

- **Hexagonal packing**: Efficient close-packed arrangement
- **Cylindrical magnets**: Approximated as K-sided polyhedra
- **Base plate**: Non-magnetic mounting structure
- **Field calculation**: Magnetic field at arbitrary points
- **VTK export**: 3D visualization in ParaView

## Usage

```bash
cd examples/smco_magnet_array
python smco_array.py
```

## Configuration

Default parameters (can be modified in `main()` function):

```python
mag_radius=0.005,        # 5 mm magnet radius
mag_height=0.01,         # 10 mm magnet height
mag_M=[0, 0, 1],         # 1 T vertical magnetization
spacing=0.01,            # 10 mm center-to-center spacing
array_radius=0.06,       # 60 mm array radius
base_plate_height=0.005  # 5 mm base plate thickness
```

Note: Each magnet is created using `rad.ObjCylMag()` with nseg=16 segments for cylindrical approximation.

## Output

```
======================================================================
SMCO MAGNET ARRAY SIMULATION
======================================================================

Creating SmCo magnet array...
  Magnet radius: 5.00 mm
  Magnet height: 10.00 mm
  Magnetization: [0, 0, 1] T
  Array radius: 60.00 mm
  Magnet spacing: 10.00 mm

  Creating base plate...
    Creating meshed disk: 6 radial × 24 angular = 144 elements
  [OK] Created 125 magnets in hexagonal array

Calculating magnetic field...
Position (m)              Bx (mT)         By (mT)         Bz (mT)         |B| (mT)
(0.000, 0.000, 0.020)     0.001707        0.000447        75.778921       75.778921
(0.000, 0.000, 0.050)     0.000020        -0.000057       42.502116       42.502116
(0.030, 0.000, 0.020)     15.668715       -0.000156       90.325789       91.674734

Exporting to VTK...
  [OK] Created: smco_array.vtk

Summary:
  Number of magnets: 125
  Array radius: 60.00 mm
  Magnet radius: 5.00 mm
  Magnet height: 10.00 mm
```

## Geometry Details

### Hexagonal Array Pattern

The magnets are arranged in a hexagonal close-packed pattern:
- **Grid spacing**: `spacing` parameter (default 10 mm)
- **Offset**: Every other row offset by `spacing/2`
- **Y-spacing**: `spacing * √3/2` for hexagonal geometry
- **Coverage**: Only magnets within `array_radius` are created

### Magnet Implementation

Each cylindrical magnet is created using `rad.ObjCylMag()`:
- **Function**: `rad.ObjCylMag([x,y,z], radius, height, nseg, axis, magnetization)`
- **Approximation**: nseg-sided prism with polygon base
- **Default nseg=16**: Good balance between accuracy and computation time
- **Axis**: 'z' (vertical orientation)

### Coordinate System

- **X-Y plane**: Magnet array plane
- **Z-axis**: Vertical (magnetization direction)
- **Origin**: Center of base plate

## Field Calculation Results

Field strength decreases with distance from array:
- **20 mm above center**: ~76 mT
- **50 mm above center**: ~43 mT
- **30 mm off-axis, 20 mm above**: ~92 mT

## Visualization

### ParaView

Open `smco_array.vtk` in ParaView:

```bash
# After running smco_array.py
# Open smco_array.vtk in ParaView
```

**File info**:
- Polygons: ~3090 (125 magnets × ~18 faces + base plate)
- Points: ~15312
- Base plate: Gray (144 meshed elements)
- Magnets: Blue (16-sided cylinders)

### View Settings

Recommended ParaView settings:
- **Color By**: Solid Color or object ID
- **Opacity**: 0.8 to see internal structure
- **Lighting**: Enable ambient light

## Code Structure

### Main Function

```python
def create_smco_magnet_array(
    mag_radius=0.005,      # Magnet radius (m)
    mag_height=0.01,       # Magnet height (m)
    mag_M=[0, 0, 1],       # Magnetization (T)
    spacing=0.01,          # Magnet spacing (m)
    array_radius=0.06,     # Array radius (m)
    base_plate_height=0.005  # Base plate height (m)
):
    """Create hexagonal array of SmCo magnets."""
```

### Helper Function

```python
def create_meshed_disk(R, H, n_radial, n_angular, x0=0, y0=0, z0=0):
    """
    Create meshed circular disk using hexahedral elements.

    Uses pentahedra for center ring and hexahedra for outer rings.
    """
```

### Magnet Creation

Magnets are created directly using `rad.ObjCylMag()`:

```python
# rad.ObjCylMag([x,y,z], radius, height, nseg, axis, magnetization)
magnet = rad.ObjCylMag([x, y, z], mag_radius, mag_height, 16, 'z', mag_M)
```

## Hexagonal Grid Algorithm

The hexagonal packing is achieved using:

```python
for nx in range(-20, 21):
    for ny in range(-20, 21):
        # Offset every other row by half spacing
        x = nx * spacing + (ny % 2) * spacing / 2

        # Vertical spacing for hexagonal pattern
        y = ny * spacing * np.sqrt(3) / 2

        # Only create magnets within array radius
        if x**2 + y**2 < array_radius**2:
            # Create magnet at (x, y)
```

This creates a close-packed hexagonal arrangement with optimal packing density.

## Material Properties

### SmCo (Samarium-Cobalt)

- **Type**: Rare-earth permanent magnet
- **Composition**: Sm₂Co₁₇ or SmCo₅
- **Remanence**: ~0.9-1.15 T (varies by grade)
- **Coercivity**: Very high (~700-900 kA/m)
- **Temperature stability**: Excellent (up to 250-350°C)
- **Advantages**: High field, good temperature stability
- **Disadvantages**: Brittle, expensive

### Simulation Parameters

In this example:
- **Magnetization**: 1 T (vertical direction)
- **No demagnetization**: Simplified linear model
- **No temperature effects**: Room temperature assumed

## Modifications

### Change Array Size

```python
# Create larger array
geometry, info = create_smco_magnet_array(
    array_radius=0.10,  # 100 mm radius
    spacing=0.015,      # 15 mm spacing
)
```

### Change Magnetization Direction

```python
# Horizontal magnetization
geometry, info = create_smco_magnet_array(
    mag_M=[1, 0, 0],  # X-direction
)
```

## Troubleshooting

### Too Many Magnets

If the simulation is slow:
1. Reduce `array_radius`
2. Increase `spacing`
3. Reduce base plate mesh resolution (n_radial, n_angular)

### Field Calculation Slow

Field calculation time scales with number of magnets:
- 125 magnets: ~1 second per point
- 500 magnets: ~4 seconds per point

### Memory Issues

Each magnet is represented by an nseg-sided prism:
- nseg=16: ~18 faces per magnet
- Consider reducing array_radius or increasing spacing

## Comparison with Original Mathematica

### Ported Features

✓ Hexagonal array generation
✓ Cylindrical magnet approximation
✓ Base plate
✓ Field calculation
✓ Geometry visualization

### Differences

- **Visualization**: VTK export instead of Mathematica Graphics3D
- **Units**: Explicit meters (m) instead of Mathematica's unit-free
- **Performance**: Python/NumPy typically faster for large arrays

## Further Reading

- [examples/complex_coil_geometry/](../complex_coil_geometry/) - CoilBuilder examples
- [examples/electromagnet/](../electromagnet/) - Yoke + coil simulation
- [src/python/README.md](../../src/python/README.md) - Radia utilities

## References

- **SmCo magnets**: https://en.wikipedia.org/wiki/Samarium-cobalt_magnet
- **Hexagonal packing**: https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres
- **Radia**: https://github.com/ochubar/Radia

---

**Created**: 2025-10-30
**Ported from**: Mathematica notebook (2023-10-01)
