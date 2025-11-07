# Radia Vector Potential to NGSolve CoefficientFunction

This example demonstrates the integration between Radia's vector potential A and NGSolve's CoefficientFunction framework, verifying that **curl(A) = B**.

## Quick Start

```python
import radia as rad
import rad_ngsolve
from ngsolve import curl

# Create magnet
dipole = rad.ObjRecMag([0, 0, 0], [20, 10, 10], [1.0, 0, 0])

# Get vector potential and magnetic field
A_cf = rad_ngsolve.RadiaField(dipole, 'a')  # Vector potential
B_cf = rad_ngsolve.RadiaField(dipole, 'b')  # Magnetic flux density

# Verify curl(A) = B
curl_A_cf = curl(A_cf)
# curl_A_cf should match B_cf

# With coordinate transformation (NEW in v0.07)
A_cf_rotated = rad_ngsolve.RadiaField(
    dipole, 'a',
    origin=[0.01, 0, 0],           # Translation
    u_axis=[0.707, 0.707, 0],      # 45° rotation
    v_axis=[-0.707, 0.707, 0],
    w_axis=[0, 0, 1]
)
```

## Overview

When working with magnetic fields, the vector potential A is fundamental:
- **B = curl(A)**: Magnetic flux density is the curl of the vector potential
- **A** satisfies: div(A) = 0 (Coulomb gauge)

This example shows how to:
1. Extract both B and A from Radia magnetostatic calculations
2. Pass A to NGSolve as a CoefficientFunction
3. Compute curl(A) numerically in NGSolve
4. Verify that curl(A) matches B from Radia
5. **NEW**: Use coordinate transformations for arbitrary magnet placement

## Files

- `verify_curl_A_equals_B.py`: Main verification script (curl(A) = B)
- `test_coordinate_transform.py`: Coordinate transformation examples
- `test_mesh_convergence.py`: Mesh convergence analysis

## Usage

### Basic verification
```bash
python verify_curl_A_equals_B.py
```

### Coordinate transformation examples
```bash
python test_coordinate_transform.py
```

### Mesh convergence test
```bash
python test_mesh_convergence.py
```

## Requirements

- Python 3.12+
- Radia (with NGSolve integration)
- NGSolve
- NumPy

## How It Works

### 1. Radia Background Field with A

The Radia background field callback can now return both B and A:

```python
def radia_field_with_A(coords):
    """Return both B and A from Radia"""
    x, y, z = coords  # mm

    # Get B field (Tesla)
    B = rad.Fld(magnet, 'b', [x, y, z])

    # Get vector potential (T*m)
    A = rad.Fld(magnet, 'a', [x, y, z])

    return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)
```

### 2. Extract A as NGSolve CoefficientFunction

```python
import rad_ngsolve

# Get vector potential as CF
A_cf = rad_ngsolve.RadiaField(bg_field, 'a')

# Get B field for comparison
B_cf = rad_ngsolve.RadiaField(bg_field, 'b')
```

### 3. Compute curl(A) in NGSolve

```python
from ngsolve import curl

# Compute curl of vector potential
curl_A_cf = curl(A_cf)

# Evaluate at a point
mip = mesh(x, y, z)
curl_A_value = curl_A_cf(mip)
B_value = B_cf(mip)

# Verify: curl_A_value ≈ B_value
```

## Expected Output

```
================================================================================
Verification: curl(A) = B for Radia Background Field
================================================================================

[Step 1] Creating Radia rectangular magnet
--------------------------------------------------------------------------------
  Magnet ID: 1
  Center: [0, 0, 0] mm
  Dimensions: [20, 20, 30] mm
  Magnetization: [0, 0, 1000] A/m

[Step 2] Creating background field wrapper for B and A
--------------------------------------------------------------------------------
  Background field ID: 2
  Test point: [0.025, 0.015, 0.04] m
  B from Radia: [0.12345 0.23456 0.34567] T
  A from Radia: [0.00012 0.00023 0.00034] T*m

[Step 3] Creating NGSolve mesh and CoefficientFunction
--------------------------------------------------------------------------------
  Mesh created: 1234 vertices, 5678 elements
  A as CoefficientFunction: <class 'ngsolve.fem.CoefficientFunction'>

[Step 4] Computing curl(A) and comparing with B
--------------------------------------------------------------------------------

  Point (m)                curl(A) (T)              B (T)                 Error
  ----------------------------------------------------------------------------
  (0.025, 0.015, 0.03)  [ 0.12345,  0.23456,  0.34567]  [ 0.12345,  0.23456,  0.34567]  1.23e-06
  (0.025, 0.015, 0.04)  [ 0.12345,  0.23456,  0.34567]  [ 0.12345,  0.23456,  0.34567]  9.87e-07
  ...

[Step 5] Verification Result
--------------------------------------------------------------------------------
  [SUCCESS] curl(A) = B verified!
  Maximum error: 1.23e-06 T (tolerance: 1.00e-04 T)

[Step 6] Creating visualization
--------------------------------------------------------------------------------
  FE space created: 12345 DOFs
  GridFunctions created for A, curl(A), and B
  [OK] VTK output saved: curl_A_verification.vtu
       Open in ParaView to visualize A, curl(A), and B

================================================================================
Verification complete!
================================================================================
```

## Visualization in ParaView

The script generates `curl_A_verification.vtu` which contains:
- **A_vector_potential**: Vector potential field
- **curl_A**: Numerically computed curl(A)
- **B_field**: Direct B field from Radia

Load in ParaView and compare the fields visually.

## Technical Details

### Units
- **Radia**: mm (coordinates), Tesla (B), T·m (A), A/m (M)
- **NGSolve**: meters (coordinates), Tesla (B), T·m (A)
- Conversion handled automatically by `rad_ngsolve` module

### Backward Compatibility

The callback can still return just B for backward compatibility:

```python
def simple_field(coords):
    return [Bx, By, Bz]  # Old format still works
```

### Dictionary Format

To provide A, use dictionary format:

```python
def field_with_A(coords):
    return {
        'B': [Bx, By, Bz],  # Tesla
        'A': [Ax, Ay, Az]   # T·m
    }
```

## Mathematical Background

The vector potential A satisfies:

1. **B = ∇ × A** (definition of B)
2. **∇ · A = 0** (Coulomb gauge)
3. **∇² A = -μ₀ J** (Poisson equation for A)

For magnetostatics with no currents outside the magnet:
- Inside magnet: **∇² A = -μ₀ ∇ × M**
- Outside magnet: **∇² A = 0**

## Coordinate Transformation (New in v0.07)

### Overview

RadiaField now supports coordinate transformations, allowing you to define magnets in a local coordinate system and place them arbitrarily in the global NGSolve mesh.

### Basic Usage

```python
import rad_ngsolve

# With translation only
A_cf = rad_ngsolve.RadiaField(
    magnet,
    field_type='a',
    origin=[0.01, 0.005, 0]  # Translation in meters
)

# With rotation (45° around z-axis)
import math
cos45 = math.cos(math.radians(45))
sin45 = math.sin(math.radians(45))

A_cf = rad_ngsolve.RadiaField(
    magnet,
    field_type='a',
    u_axis=[cos45, sin45, 0],   # Rotated u-axis
    v_axis=[-sin45, cos45, 0],  # Rotated v-axis
    w_axis=[0, 0, 1]             # z-axis unchanged
)

# Combined translation + rotation
A_cf = rad_ngsolve.RadiaField(
    magnet,
    field_type='a',
    origin=[0.01, 0.005, 0],
    u_axis=[cos45, sin45, 0],
    v_axis=[-sin45, cos45, 0],
    w_axis=[0, 0, 1]
)
```

### How It Works

1. **Forward transformation** (Global → Local):
   - Point `p_global` in NGSolve mesh
   - Translate: `p' = p_global - origin`
   - Rotate: `p_local = [u·p', v·p', w·p']`
   - Evaluate Radia field at `p_local`

2. **Inverse transformation** (Local → Global):
   - Field `F_local` from Radia
   - Transform: `F_global = u*F_local[0] + v*F_local[1] + w*F_local[2]`
   - Return `F_global` to NGSolve

3. **Invariance of curl**:
   - `curl(A)` is invariant under orthogonal transformations
   - Therefore: `curl(A_global) = B_global` remains valid

### Automatic Normalization

Axis vectors are automatically normalized:

```python
A_cf = rad_ngsolve.RadiaField(
    magnet, 'a',
    u_axis=[2, 0, 0],    # Length = 2 → normalized to [1, 0, 0]
    v_axis=[0, 3, 0],    # Length = 3 → normalized to [0, 1, 0]
    w_axis=[0, 0, 0.5]   # Length = 0.5 → normalized to [0, 0, 1]
)
```

### Physical Interpretation

- Vector potential A transforms as a vector under rotations
- The curl operator commutes with orthogonal transformations
- Physical fields (B, H, M) transform consistently
- Gauge invariance is preserved

### Use Cases

- **Rotated magnets**: Define magnet in natural orientation, rotate to placement
- **Multi-magnet systems**: Each magnet in its own local frame
- **Rotating machinery**: Time-dependent coordinate transformations
- **Complex geometries**: Simplify magnet definition in local coordinates

### Example: Dipole Rotated 45°

```python
# Create dipole along local x-axis
dipole = rad.ObjRecMag(
    [0, 0, 0],      # Center
    [20, 10, 10],   # Elongated along x
    [1.0, 0, 0]     # Magnetization along x
)

# Rotate 45° around z-axis
cos45 = 0.7071
sin45 = 0.7071

A_cf = rad_ngsolve.RadiaField(
    dipole, 'a',
    u_axis=[cos45, sin45, 0],
    v_axis=[-sin45, cos45, 0],
    w_axis=[0, 0, 1]
)

# curl(A) will give correctly rotated B field
curl_A = curl(A_cf)
```

### See Examples

- `test_coordinate_transform.py`: Comprehensive coordinate transformation tests
- `../coordinate_transformation/`: More transformation examples

## See Also

- `examples/Radia_to_NGSolve_CoefficientFunction`: B field integration
- `examples/NGSolve_CoefficientFunction_to_Radia_BackgroundField`: NGSolve to Radia
- `examples/coordinate_transformation/`: Coordinate transformation documentation
