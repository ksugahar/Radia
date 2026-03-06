# Magpylib + Radia Integration

This folder contains examples demonstrating how to use [magpylib](https://github.com/magpylib/magpylib) to define external magnetic fields for Radia MMM computations.

## Overview

magpylib is an open-source Python library for magnetic field computation using analytical expressions. By combining magpylib with Radia, you can:

1. Use magpylib's analytical permanent magnet and current source models as background fields
2. Solve for magnetization in soft magnetic materials subject to these fields
3. Compute the total field (magpylib source + induced field from magnetized material)

## Requirements

\Requirement already satisfied: radia in c:\program files\python312\lib\site-packages (1.3.9)
Requirement already satisfied: numpy>=1.20 in c:\program files\python312\lib\site-packages (from radia) (2.3.5)
## CRITICAL: ObjBckgCF Callback Unit Behavior

**ObjBckgCF passes positions to the callback in Radia's internal units (millimeters), regardless of the \ setting.**

This is a critical detail that can cause ~1000x errors in field evaluation if not handled correctly.

### Correct Implementation

### Incorrect Implementation (DO NOT USE)

## Key Function: 
The adapter function converts a magpylib source to a Radia \ callback:

## Examples

### demo_magpylib_integration.py

Simple demo comparing ObjBckgCF (magpylib) with ObjBckg (uniform field).
Verifies that the unit conversion fix produces correct results.

### sphere_in_halbach_cylinder.py

Demonstrates three scenarios with cuboid magnets:
1. Single permanent magnet above iron cube
2. Helmholtz coil pair creating uniform field
3. Quadrupole magnet configuration

### cylinder_magnet_examples.py

Comprehensive cylindrical magnet examples:

1. **Axially Magnetized Cylinder** - Most common PM configuration
   - Magnetization along cylinder axis
   - Creates strong axial field on axis

2. **Diametrically Magnetized Cylinder** - Transverse field
   - Magnetization perpendicular to axis
   - Used for magnetic couplings

3. **Halbach Cylinder (k=2 Dipole)** - Strong uniform field in bore
   - Creates ~0.66 T uniform field in center
   - Minimal external field
   - 12 CylinderSegment elements

4. **Ring Magnet Pair** - Axial gradient field
   - Two opposing ring magnets
   - Zero field at center, strong gradient
   - Used for magnetic levitation

## Supported magpylib Sources

All magpylib source types work with ObjBckgCF:

| Source Type | magpylib Class | Use Case |
|-------------|----------------|----------|
| Cuboid | \ | Rectangular magnets |
| Cylinder | \ | Axial/diametric magnets |
| CylinderSegment | \ | Ring magnets, Halbach arrays |
| Sphere | \ | Spherical magnets |
| Circle (current) | \ | Coils, Helmholtz pairs |
| Line (current) | \ | Wire segments |
| Collection | \ | Combined sources |

## Usage Example

## Notes

1. **Unit Conversion**: ObjBckgCF always passes positions in mm. Always convert to meters for magpylib.

2. **Return Type**: The callback must return Python native floats, not numpy.float64. Use \ conversion.

3. **Field Evaluation**: Radia's \ is most accurate OUTSIDE magnetic materials. Evaluate at external observation points.

4. **Halbach Arrays**: Use \ with varying magnetization angles. More segments = more uniform field.

## License

BSD 2-Clause License (same as magpylib)
