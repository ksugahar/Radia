# Radia Vector Potential to NGSolve CoefficientFunction

This example demonstrates the integration between Radia's vector potential A and NGSolve's CoefficientFunction framework, verifying that **curl(A) = B**.

## Overview

When working with magnetic fields, the vector potential A is fundamental:
- **B = curl(A)**: Magnetic flux density is the curl of the vector potential
- **A** satisfies: div(A) = 0 (Coulomb gauge)

This example shows how to:
1. Extract both B and A from Radia magnetostatic calculations
2. Pass A to NGSolve as a CoefficientFunction
3. Compute curl(A) numerically in NGSolve
4. Verify that curl(A) matches B from Radia

## Files

- `verify_curl_A_equals_B.py`: Main verification script

## Usage

```bash
python verify_curl_A_equals_B.py
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

## See Also

- `examples/Radia_to_NGSolve_CoefficientFunction`: B field integration
- `examples/NGSolve_CoefficientFunction_to_Radia_BackgroundField`: NGSolve to Radia
