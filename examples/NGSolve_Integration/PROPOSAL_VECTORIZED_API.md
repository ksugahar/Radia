# Proposal: Vectorized Field Evaluation API

**Date**: 2025-11-08
**Goal**: Accelerate GridFunction.Set() by reducing Python call overhead

---

## Problem Statement

Current performance bottleneck in `GridFunction.Set()`:

```python
# Current implementation in rad_ngsolve.cpp
for each vertex in mesh:  # M iterations
    coords = [x*1000, y*1000, z*1000]
    field = rad.Fld(obj, 'b', coords)  # Python call each time!
```

**Overhead per call**: ~10-50 µs (Python ↔ C++ transition)
**Total for M=5000 vertices**: 50-250 ms just for Python calls!

---

## Proposed Solution

### Add Vectorized Field Evaluation

**New Python API**:
```python
# Evaluate field at multiple points in one call
coords_array = [[x1, y1, z1], [x2, y2, z2], ...]  # M×3 array
fields = rad.FldVec(obj, 'b', coords_array)  # Single Python call
# Returns: M×3 array of field values
```

### Implementation Path

#### 1. Python API (already exists!)

`radentry.cpp:1011` already supports multiple points:
```cpp
int CALL RadFld(double* pB, int* pNb, int Obj, char* ID, double* pCoord, int Np)
{
    // ...
    FieldArbitraryPointsArray((long)Obj, ID, PointsArray, (long)Np);
    // ...
}
```

**The C API already supports batch evaluation!**

#### 2. Modify rad_ngsolve.cpp

Current implementation:
```cpp
virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
                     FlatVector<> result) const override
{
    // Single point evaluation
    py::object field_result = rad.attr("Fld")(radia_obj, field_type, coords);
}
```

Proposed addition:
```cpp
// Add batch evaluation method
virtual void Evaluate(const ngfem::SIMD_BaseMappedIntegrationRule& mir,
                     ngfem::BareSliceMatrix<SIMD<double>> values) const override
{
    // Batch evaluation for all points in mir
    size_t npts = mir.Size();

    // Collect all points
    py::array_t<double> coords_array({(int)npts, 3});
    auto coords_buf = coords_array.mutable_unchecked<2>();

    for (size_t i = 0; i < npts; i++) {
        auto pnt = mir[i].GetPoint();
        // ... coordinate transformation ...
        coords_buf(i, 0) = p_local[0] * 1000.0;
        coords_buf(i, 1) = p_local[1] * 1000.0;
        coords_buf(i, 2) = p_local[2] * 1000.0;
    }

    // Single Python call for all points!
    py::module_ rad = py::module_::import("radia");
    py::object result = rad.attr("Fld")(radia_obj, field_type, coords_array);

    // Extract results (single loop, no Python calls)
    // ...
}
```

---

## Expected Performance Improvement

### Current Performance

| M vertices | Python calls | Python overhead | Compute time | Total |
|-----------|-------------|-----------------|--------------|-------|
| 135 | 135 | ~5 ms | ~10 ms | 15 ms |
| 5034 | 5034 | ~250 ms | 7100 ms | 7350 ms |

### After Vectorization

| M vertices | Python calls | Python overhead | Compute time | Total |
|-----------|-------------|-----------------|--------------|-------|
| 135 | 1 | ~0.05 ms | ~10 ms | **10 ms** |
| 5034 | 1 | ~0.05 ms | 7100 ms | **7100 ms** |

**Speedup**:
- Small meshes: ~1.5x
- Large meshes: ~1.03x (overhead becomes negligible)

### With OpenMP Parallelization

The existing OpenMP loop (`#pragma omp parallel for if(Np > 100)` in radapl3.cpp:167) will now be effective:

| M vertices | Threads | Compute time (parallel) | Total |
|-----------|---------|------------------------|-------|
| 5034 | 8 cores | ~900 ms | **900 ms** |

**Combined speedup**: **~8x** for large meshes!

---

## Implementation Steps

### Step 1: Verify C API works

Test existing batch evaluation:
```python
import radia as rad
magnet = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1.2])

# Test batch evaluation (should work already)
points = [[10,0,0], [0,10,0], [0,0,10]]  # 3 points
fields = rad.Fld(magnet, 'b', points)  # Does this work?
```

### Step 2: Modify rad_ngsolve.cpp

Add SIMD evaluation method for batch processing.

### Step 3: Benchmark

Compare old vs new implementation.

---

## Alternative: Fast Multipole Method (FMM)

For true O(M) complexity, implement FMM:
- Current: O(M × N) per evaluation
- With FMM: O(M + N) per evaluation

**Effort**: Large (weeks of development)
**Benefit**: Fundamental scaling improvement

---

## Recommendation

1. **Immediate**: Implement vectorized API (small change, 2-8x speedup)
2. **Future**: Consider FMM for large-scale problems (N > 1000, M > 10000)

---

**Status**: Proposal - awaiting approval for implementation
