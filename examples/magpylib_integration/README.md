# Magpylib + Radia Integration

This folder contains examples demonstrating how to use [magpylib](https://github.com/magpylib/magpylib) to define external magnetic fields for Radia MMM computations.

## Overview

magpylib is an open-source Python library for magnetic field computation using analytical expressions. By combining magpylib with Radia, you can:

1. Use magpylib's analytical permanent magnet and current source models as background fields
2. Solve for magnetization in soft magnetic materials subject to these fields
3. Compute the total field (magpylib source + induced field from magnetized material)

## Requirements

```bash
pip install magpylib>=5.0
pip install radia  # or build from source
```

## CRITICAL: ObjBckgCF Callback Unit Behavior

**ObjBckgCF passes positions to the callback in Radia's internal units (millimeters), regardless of the `rad.FldUnits()` setting.**

This is a critical detail that can cause ~1000x errors in field evaluation if not handled correctly.

### Correct Implementation

```python
# CORRECT - Always convert mm to m for magpylib
def magpylib_field(pos):
    """Callback for ObjBckgCF.

    IMPORTANT: pos is ALWAYS in mm (Radia internal units),
    regardless of FldUnits() setting. magpylib expects meters.
    """
    pos_m = [p * 0.001 for p in pos]  # mm -> m
    B = magpy_source.getB(pos_m)
    return [float(B[0]), float(B[1]), float(B[2])]
```

### Incorrect Implementation (DO NOT USE)

```python
# WRONG - This assumes FldUnits() affects callback coordinates
def magpylib_field(pos):
    B = magpy_source.getB(pos)  # pos is in mm, magpylib expects m!
    return [float(B[0]), float(B[1]), float(B[2])]
```

## Key Function: `magpylib_to_radia_callback()`

The adapter function converts a magpylib source to a Radia `ObjBckgCF` callback:

```python
def magpylib_to_radia_callback(magpy_source):
    """
    Create a Radia ObjBckgCF callback from a magpylib source.

    CRITICAL: ObjBckgCF passes positions in mm (Radia internal units),
    regardless of the FldUnits() setting. magpylib expects meters.
    This function handles the mm -> m conversion automatically.

    Parameters:
    - magpy_source: Any magpylib source (Cuboid, Circle, Collection, etc.)

    Returns:
    - Callback function: pos [mm] -> [Bx, By, Bz] in Tesla
    """
    def field_callback(pos):
        pos_m = [p * 0.001 for p in pos]  # mm -> m conversion
        B = magpy_source.getB(pos_m)
        return [float(B[0]), float(B[1]), float(B[2])]

    return field_callback
```

## Usage Example

```python
import radia as rad
import magpylib as magpy
import numpy as np

rad.FldUnits('m')
rad.UtiDelAll()

# Create magpylib source (permanent magnet)
magnet = magpy.magnet.Cuboid(
    magnetization=(0, 0, 1e6),  # A/m (approximately 1.25 T remanence)
    dimension=(0.05, 0.05, 0.02),  # 50x50x20 mm in meters
    position=(0, 0, 0.05)  # 50mm above origin in meters
)

# Create callback with mm -> m conversion
def field_callback(pos):
    """ObjBckgCF callback - pos is in mm, convert to m for magpylib."""
    pos_m = [p * 0.001 for p in pos]
    B = magnet.getB(pos_m)
    return [float(B[0]), float(B[1]), float(B[2])]

bg_field = rad.ObjBckgCF(field_callback)

# Create soft iron cube
cube = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.04], [0, 0, 0])
rad.ObjDivMag(cube, [4, 4, 4])

# Apply material
mat = rad.MatLin(999.0)  # mu_r = 1000
rad.MatApl(cube, mat)

# Solve
system = rad.ObjCnt([cube, bg_field])
rad.Solve(system, 0.0001, 1000)

# Get results
M = rad.ObjM(cube)
```

## Examples

### demo_magpylib_integration.py

Simple demo comparing ObjBckgCF (magpylib) with ObjBckg (uniform field).

### sphere_in_halbach_cylinder.py

Demonstrates three scenarios:
1. Single permanent magnet above iron cube
2. Helmholtz coil pair creating uniform field
3. Quadrupole magnet configuration

## Notes

1. **Unit Conversion**: ObjBckgCF always passes positions in mm. Always convert to meters for magpylib.

2. **Return Type**: The callback must return Python native floats, not numpy.float64. Use `float()` conversion.

3. **Field Evaluation**: Radia's `rad.Fld()` is most accurate OUTSIDE magnetic materials. Evaluate at external observation points.

## License

BSD 2-Clause License (same as magpylib)
