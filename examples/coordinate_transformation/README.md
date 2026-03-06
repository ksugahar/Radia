# Coordinate Transformation in rad_ngsolve.RadiaField

## Overview

This feature allows you to define a local coordinate system for Radia field calculations in NGSolve. This is useful when:
- Your magnet is positioned/rotated relative to the global coordinate system
- You want to work in a local coordinate frame
- You need to transform field results between different coordinate systems

## New Parameters

### `RadiaField` Constructor

```python
import rad_ngsolve

# Basic usage (no transformation)
B_cf = rad_ngsolve.RadiaField(radia_obj, field_type='b')

# With coordinate transformation
B_cf = rad_ngsolve.RadiaField(
    radia_obj,
    field_type='b',
    origin=[x, y, z],      # Translation vector in meters
    u_axis=[ux, uy, uz],   # Local u-axis (will be normalized)
    v_axis=[vx, vy, vz],   # Local v-axis (will be normalized)
    w_axis=[wx, wy, wz]    # Local w-axis (will be normalized)
)
```

### Parameters

- **origin** (optional): Translation vector `[x, y, z]` in meters
  - Default: `[0, 0, 0]` (no translation)
  - Defines the origin of the local coordinate system in global coordinates

- **u_axis** (optional): Local u-axis `[ux, uy, uz]`
  - Default: `[1, 0, 0]` (global x-axis)
  - Will be automatically normalized to unit length

- **v_axis** (optional): Local v-axis `[vx, vy, vz]`
  - Default: `[0, 1, 0]` (global y-axis)
  - Will be automatically normalized to unit length

- **w_axis** (optional): Local w-axis `[wx, wy, wz]`
  - Default: `[0, 0, 1]` (global z-axis)
  - Will be automatically normalized to unit length

## Coordinate Transformation

### Forward Transformation (Global → Local)

When evaluating the field at a point `p_global` in NGSolve:

1. **Translation**: `p' = p_global - origin`
2. **Rotation to local frame**:
   ```
   p_local[0] = u_axis · p'
   p_local[1] = v_axis · p'
   p_local[2] = w_axis · p'
   ```
3. **Radia calculation**: Field calculated at `p_local` (converted to mm)

### Inverse Transformation (Local → Global)

The field result is transformed back to global coordinates:

```
F_global = u_axis * F_local[0] + v_axis * F_local[1] + w_axis * F_local[2]
```

Where:
- `F_local` is the field from Radia in local coordinates
- `F_global` is the field in NGSolve's global coordinates

## Examples

### Example 1: Translation Only

Move the coordinate origin to `[10mm, 0, 0]`:

```python
import radia as rad
import rad_ngsolve

# Create magnet at origin in Radia
dipole = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
rad.MatApl(dipole, rad.MatStd('NdFeB', 1.2))

# Create field with translated origin
B_cf = rad_ngsolve.RadiaField(
    dipole,
    field_type='b',
    origin=[0.01, 0, 0]  # 10mm in meters
)

# Now when NGSolve evaluates at global point [0.02, 0, 0],
# it will calculate field at [0.01, 0, 0] in local coords
# which is [10, 0, 0] mm for Radia
```

### Example 2: 90° Rotation Around Z-Axis

Rotate the coordinate system 90° counterclockwise around z:

```python
B_cf = rad_ngsolve.RadiaField(
    dipole,
    field_type='b',
    u_axis=[0, 1, 0],    # Global y becomes local u
    v_axis=[-1, 0, 0],   # Global -x becomes local v
    w_axis=[0, 0, 1]     # Global z stays as local w
)
```

### Example 3: Non-Normalized Axes (Auto-Normalization)

The axes don't need to be unit vectors - they will be automatically normalized:

```python
B_cf = rad_ngsolve.RadiaField(
    dipole,
    field_type='b',
    u_axis=[2, 0, 0],    # Length = 2 → normalized to [1, 0, 0]
    v_axis=[0, 3, 0],    # Length = 3 → normalized to [0, 1, 0]
    w_axis=[0, 0, 0.5]   # Length = 0.5 → normalized to [0, 0, 1]
)
```

### Example 4: Combined Translation and Rotation

Translate by `[5mm, 5mm, 0]` and rotate 45° in xy-plane:

```python
import math

B_cf = rad_ngsolve.RadiaField(
    dipole,
    field_type='b',
    origin=[0.005, 0.005, 0],           # Translation
    u_axis=[1, 1, 0],                    # 45° in xy-plane
    v_axis=[-1, 1, 0],                   # Perpendicular to u
    w_axis=[0, 0, 1]                     # z unchanged
)
```

### Example 5: Vector Potential with Transformation

The coordinate transformation works with all field types:

```python
# Magnetic flux density
B_cf = rad_ngsolve.RadiaField(dipole, 'b', origin=[0.01, 0, 0])

# Magnetic field
H_cf = rad_ngsolve.RadiaField(dipole, 'h', origin=[0.01, 0, 0])

# Vector potential
A_cf = rad_ngsolve.RadiaField(dipole, 'a', origin=[0.01, 0, 0])

# Magnetization
M_cf = rad_ngsolve.RadiaField(dipole, 'm', origin=[0.01, 0, 0])
```

## Technical Details

### Rotation Matrix

The transformation is defined by the rotation matrix:

```
R = [u_axis | v_axis | w_axis]
```

Where each axis is a column vector (after normalization).

### Coordinate Transformation Formulas

**Point transformation (global → local)**:
```
p_local = R^T * (p_global - origin)
```

**Field transformation (local → global)**:
```
F_global = R * F_local
```

### Unit Conversion

The coordinate transformation is applied **before** unit conversion:
1. NGSolve point in meters
2. Apply transformation (still in meters)
3. Convert to millimeters for Radia
4. Calculate field
5. Transform field back to global
6. Apply unit scaling (e.g., for vector potential)

## Notes

- All axes are automatically normalized to unit length
- The default is an identity transformation (no change)
- The `use_transform` attribute indicates whether transformation is active
- Works with all field types: 'b', 'h', 'a', 'm'
- Origin is specified in meters (NGSolve units)
- Axes can be specified as lists or tuples

## Version History

- **v0.07** (2025-11-07): Added coordinate transformation support
  - Added `origin`, `u_axis`, `v_axis`, `w_axis` parameters
  - Automatic normalization of axis vectors
  - Forward and inverse coordinate transformations

## See Also

- [Radia Documentation](http://radia.forge.jp/)
- [NGSolve Documentation](https://ngsolve.org/)
- [Python examples](../ngsolve_integration/)
