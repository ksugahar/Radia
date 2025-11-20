# PrepareCache() Optimization Proposal

## Problem Analysis

### Current Performance

**Measured with profile_batch_performance.py:**
- Radia.Fld() batch: 2000 points in **1ms** (0.5 us/point) ✓ FAST
- PrepareCache(): 500 points in **>60s** (>120 ms/point) ✗ VERY SLOW

**Bottleneck:** C++ PrepareCache() implementation

### Root Cause

PrepareCache() makes thousands of pybind11 API calls in loops:

```cpp
for (size_t i = 0; i < npts; i++) {
    py::list coords;          // Python object creation
    coords.append(val1);      // Python method call
    coords.append(val2);      // Python method call
    coords.append(val3);      // Python method call
    radia_points.append(coords);  // Python method call
}
// Total: npts * 4 = 2000 * 4 = 8000 Python API calls
```

For 3360 points: **13,440+ pybind11 API calls** → 60+ seconds

**pybind11 overhead dominates**, not Radia evaluation.

## Solution 1: Pure Python PrepareCache() [RECOMMENDED]

Implement PrepareCache() entirely in Python, store results in C++ cache.

### Architecture

```python
# Python side (fast for list operations)
def PrepareCache_Python(cf, points_meters):
    import radia as rad

    # 1. Build Radia points (Python list ops are fast)
    radia_points = [[x*1000, y*1000, z*1000] for x,y,z in points_meters]

    # 2. Single batch Radia call (fast: 2000 pts in 1ms)
    results = rad.Fld(cf.radia_obj, cf.field_type, radia_points)

    # 3. Store in C++ cache via new SetCache() method
    cf._SetCacheData(points_meters, results)  # Single C++ call
```

```cpp
// C++ side (minimal pybind11 calls)
void _SetCacheData(py::list points, py::list results) {
    // Extract all data in bulk (one Python API call per array)
    std::vector<std::array<double,3>> pts = extract_points(points);
    std::vector<std::array<double,3>> res = extract_results(results);

    // Process with pure C++ (no Python API calls)
    for (size_t i = 0; i < pts.size(); i++) {
        uint64_t hash = hash_point(pts[i][0], pts[i][1], pts[i][2]);
        point_cache_[hash] = res[i];
    }
}
```

### Benefits

- **10-100x faster**: Python list ops + single C++ call
- **Simpler**: Less C++ code, easier to maintain
- **Flexible**: Easy to add features in Python

### Estimated Performance

- 3360 points: ~10ms (vs current 60s)
- 10000 points: ~30ms (linear scaling)

## Solution 2: NumPy Array Interface [COMPLEX]

Use NumPy arrays for bulk data transfer.

```cpp
void PrepareCache(py::array_t<double> points_array) {
    auto pts = points_array.unchecked<2>();  // Fast NumPy access
    // Process without Python API calls
}
```

### Benefits

- Even faster for very large arrays
- Zero-copy data transfer

### Drawbacks

- Requires NumPy dependency
- More complex implementation
- NumPy C API learning curve

## Solution 3: Pre-allocated Buffer [MODERATE]

Allocate Radia points list once, reuse.

```cpp
// Pre-allocate and reuse
py::list radia_points = py::list(npts);
for (size_t i = 0; i < npts; i++) {
    radia_points[i] = py::make_tuple(x, y, z);  // Direct assignment, no append
}
```

### Benefits

- Reduces allocations
- Still uses pybind11

### Drawbacks

- Still has pybind11 overhead
- Complex index management

## Recommendation

**Implement Solution 1: Pure Python PrepareCache()**

**Reasons:**
1. Simplest implementation (mostly Python code)
2. Leverages Python's fast list operations
3. Minimizes pybind11 overhead (single C++ call)
4. Easy to test and maintain
5. Expected 100-1000x speedup

**Implementation Plan:**
1. Add `_SetCacheData(points, results)` C++ method
2. Implement `PrepareCache()` in Python wrapper
3. Keep C++ `PrepareCache()` for backward compatibility (mark deprecated)
4. Test with 3000-5000 points

**Migration Path:**
```python
# Old API (slow)
A_cf.PrepareCache(points)  # C++ implementation

# New API (fast)
import radia_ngsolve_fast
radia_ngsolve_fast.prepare_cache(A_cf, points)  # Python implementation
```

## Conclusion

The bottleneck is **pybind11 overhead in loops**, not Radia evaluation.

**Solution: Move list-heavy operations to Python side, keep only caching in C++.**

Expected result: **3360 points in <100ms** (600x improvement)
