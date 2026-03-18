# FMM Field Evaluation Examples

This folder contains examples demonstrating batch field evaluation in Radia.

## Overview

Radia provides efficient batch field evaluation via the unified `rad.Fld()` API.
When given points with shape `(N, 3)`, it automatically uses batch evaluation
which is significantly faster than calling `rad.Fld()` with single points in a loop.

## API

```python
import radia as rad
import numpy as np

# Single point: shape (3,)
B = rad.Fld(obj, 'b', np.array([0, 0, 0.1]))

# Batch: shape (N, 3) - auto-detected, parallelized
points = np.array([[0, 0, 0.1], [0, 0, 0.2], [0, 0, 0.3]])
B_batch = rad.Fld(obj, 'b', points)   # returns (N, 3) array
H_batch = rad.Fld(obj, 'h', points)   # returns (N, 3) array
A_batch = rad.Fld(obj, 'a', points)   # returns (N, 3) array
phi = rad.Fld(obj, 'phi', points)     # returns (N,) array
```

| Field type | Returns | Units |
|------------|---------|-------|
| `'b'` | B field | Tesla |
| `'h'` | H field | A/m |
| `'a'` | Vector potential A | T*m |
| `'phi'` | Scalar potential | A |

## Example Scripts

### [demo_fldbatch.py](demo_fldbatch.py)

Demonstrates batch field evaluation with the unified `rad.Fld()` API.

**Features:**
- Creates subdivided magnetized cube (3x3x3 = 27 hexahedral elements)
- Computes field at 1000 observation points (10x10x10 grid)
- Compares performance: batch vs single-point loop
- Shows typical speedup of 10-50x

**Usage:**
```bash
cd examples/fmm_field_evaluation
python demo_fldbatch.py
```

### [verify_fmm_fldbatch.py](verify_fmm_fldbatch.py)

Verifies batch `rad.Fld()` by comparing results with single-point evaluation.

## Performance Guidelines

### When to Use Batch Evaluation

| Scenario | Recommended API |
|----------|-----------------|
| Single point field | `rad.Fld(obj, 'b', point)` with shape (3,) |
| < 10 points | `rad.Fld()` in loop (low overhead) |
| 10-100 points | Batch starts showing benefit |
| > 100 points | **Batch** (significant speedup) |
| > 10,000 points | **Batch** (essential) |

### Performance Factors

1. **Python call overhead**: Batch has single Python-C++ call vs N calls for loop
2. **TaskManager parallelization**: Batch uses NGSolve TaskManager for multi-threading
3. **Memory access patterns**: Batch optimizes memory access for batch operations

### Typical Speedups

| Points | Elements | Speedup (batch vs single-point loop) |
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

- [API Reference](../../docs/api/API_REFERENCE.md)
- [Simple Problems Examples](../simple_problems/)
- [NGSolve Integration Examples](../ngsolve_integration/)
