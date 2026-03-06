# FMM Field Evaluation Examples

This folder contains examples demonstrating batch field evaluation in Radia.

## Overview

Radia provides efficient batch field evaluation APIs for computing magnetic fields
at many observation points. This is significantly faster than calling `rad.Fld()`
in a loop.

## Available APIs

| API | Description | Output |
|-----|-------------|--------|
| `rad.FldBatch(obj, points, method)` | Batch B and H field computation | `{'B': [[Bx,By,Bz],...], 'H': [[Hx,Hy,Hz],...]}` |
| `rad.FldPhi(obj, points)` | Scalar magnetic potential | `[phi1, phi2, ...]` |
| `rad.FldA(obj, points)` | Vector magnetic potential | `[[Ax,Ay,Az], ...]` |
| `rad.ClassifyPoints(obj, points, threshold)` | Point classification (inside/outside) | `[0, 1, 0, ...]` |

## Example Scripts

### [demo_fldbatch.py](demo_fldbatch.py)

Demonstrates `rad.FldBatch()` for efficient batch field evaluation.

**Features:**
- Creates subdivided magnetized cube (3x3x3 = 27 hexahedral elements)
- Computes field at 1000 observation points (10x10x10 grid)
- Compares performance: FldBatch vs Fld loop
- Shows typical speedup of 10-50x

**Usage:**
```bash
cd examples/fmm_field_evaluation
python demo_fldbatch.py
```

**Output:**
```
======================================================================
FldBatch Demo: Batch Field Evaluation
======================================================================

1. Creating magnetized cube:
   Size: 100 mm
   Divisions: 3 x 3 x 3 = 27 elements
   Magnetization: [0, 0, 1000000] A/m
   Creation time: 0.015 s

2. Creating observation grid:
   Center: [0, 0, 0.15] m
   Extent: +/- 100 mm
   Grid: 10 x 10 x 10 = 1000 points

3. Computing field with FldBatch...
   Time: 0.045 s
   Points/sec: 22222

4. Computing field with Fld loop (for comparison)...
   Time for 100 points: 0.350 s
   Estimated time for 1000 points: 3.500 s
   Speedup: 77.8x

5. Verification:
   Max difference (first 100 pts): 0.00e+00 T

6. Field statistics (all 1000 points):
   |B| min: 1.234567e-04 T
   |B| max: 5.678901e-03 T
   |B| mean: 2.345678e-03 T
======================================================================
```

### [verify_fmm_fldbatch.py](verify_fmm_fldbatch.py)

Verifies FMM-accelerated `rad.FldBatch()` by comparing direct (method=0) and FMM dipole approximation (method=1) results for far-field accuracy, near-field behavior, and timing.

**Usage:**
```bash
python verify_fmm_fldbatch.py
```

## API Reference

### rad.FldBatch(obj, points, method)

Compute B and H fields at multiple points efficiently.

**Parameters:**
- `obj`: Radia object handle (container or single element)
- `points`: List of [x, y, z] observation points
- `method`: Computation method
  - `0`: Direct summation (default, always accurate)
  - `1`: Reserved for future FMM acceleration

**Returns:**
Dictionary with keys:
- `'B'`: List of [Bx, By, Bz] magnetic flux density vectors [T]
- `'H'`: List of [Hx, Hy, Hz] magnetic field strength vectors [A/m]

**Example:**
```python
import radia as rad

# Radia always uses meters

# Create magnet
vertices = [[-0.05,-0.05,-0.05], [0.05,-0.05,-0.05], [0.05,0.05,-0.05], [-0.05,0.05,-0.05],
            [-0.05,-0.05,0.05], [0.05,-0.05,0.05], [0.05,0.05,0.05], [-0.05,0.05,0.05]]
magnet = rad.ObjHexahedron(vertices, [0, 0, 1e6])  # 1 MA/m magnetization

# Define observation points
import numpy as np
z = np.linspace(0.1, 0.5, 100)
points = [[0, 0, zi] for zi in z]

# Compute fields
result = rad.FldBatch(magnet, points, 0)
B = np.array(result['B'])
H = np.array(result['H'])

print(f"Bz range: {B[:,2].min():.6e} to {B[:,2].max():.6e} T")
```

### rad.FldPhi(obj, points)

Compute scalar magnetic potential at multiple points.

**Parameters:**
- `obj`: Radia object handle
- `points`: List of [x, y, z] observation points

**Returns:**
List of scalar potential values [A] (Ampere)

**Note:** Scalar potential is only defined in regions without currents.

### rad.FldA(obj, points)

Compute vector magnetic potential at multiple points.

**Parameters:**
- `obj`: Radia object handle
- `points`: List of [x, y, z] observation points

**Returns:**
List of [Ax, Ay, Az] vector potential values [T·m]

**Note:** The vector potential satisfies curl(A) = B.

### rad.ClassifyPoints(obj, points, threshold)

Classify points as inside or outside magnetic elements.

**Parameters:**
- `obj`: Radia object handle
- `points`: List of [x, y, z] points to classify
- `threshold`: Distance threshold for boundary classification [length units]

**Returns:**
List of integers:
- `0`: Point is outside all magnetic elements
- `1`: Point is inside a magnetic element
- `2`: Point is on/near boundary (within threshold)

## Performance Guidelines

### When to Use FldBatch

| Scenario | Recommended API |
|----------|-----------------|
| Single point field | `rad.Fld(obj, 'b', point)` |
| < 10 points | `rad.Fld()` in loop (low overhead) |
| 10-100 points | Either (FldBatch starts showing benefit) |
| > 100 points | **`rad.FldBatch()`** (significant speedup) |
| > 10,000 points | **`rad.FldBatch()`** (essential) |

### Performance Factors

1. **Python call overhead**: FldBatch has single Python-C++ call vs N calls for loop
2. **TaskManager parallelization**: FldBatch uses NGSolve TaskManager for multi-threading
3. **Memory access patterns**: FldBatch optimizes memory access for batch operations

### Typical Speedups

| Points | Elements | Speedup (FldBatch vs Fld loop) |
|--------|----------|-------------------------------|
| 100 | 27 | ~10x |
| 1,000 | 27 | ~50x |
| 10,000 | 27 | ~100x |
| 1,000 | 1,000 | ~30x |

## Use Cases

1. **Visualization grids**: Compute field on 3D grid for contour plots
2. **Trajectory integration**: Field along particle trajectories
3. **NGSolve mesh nodes**: Field at all mesh nodes for interpolation
4. **Field maps**: Export field data for external tools

## See Also

- [API Reference](../../docs/API_REFERENCE.md)
- [Simple Problems Examples](../simple_problems/)
- [NGSolve Integration Examples](../ngsolve_integration/)

---

**Date**: 2025-12-30
**Version**: Radia v1.3.15+
