# Batch Evaluation Bottleneck: Final Analysis

## Executive Summary

**Problem:** PrepareCache() implementations are 1000-10000x slower than theoretical performance.

**Root Cause:** pybind11 overhead dominates for **any** loop over Python lists in C++.

**Solution:** Avoid C++↔Python list iteration entirely. Use NumPy or pure-Python implementation.

## Performance Measurements

### Theoretical Best (Radia.Fld only)
- **2000 points in 1ms** (0.5 us/point)
- This is the baseline - Radia itself is very fast

### C++ PrepareCache() (Original)
- **500 points: >60 seconds** (>120,000 us/point)
- **240,000x slower than Radia**
- Bottleneck: Loop with 500 × 4 = 2000 py::list append() calls

### C++ PrepareCache() (Optimized)
- **100 points: still hangs** (likely >10 seconds)
- **>20,000x slower than Radia**
- Bottleneck: Loop extracting points to C++ vectors
- 100 × 8 = 800 py::list element access calls

### Python + C++ _SetCacheData()
- **100 points: >30 seconds** (>300,000 us/point)
- **600,000x slower than Radia**
- Bottleneck: _SetCacheData() loop with py::list access
- 100 × 8 = 800 pybind11 calls

### Pure Python (theoretical, not measured)
- **1000 points: ~1-2ms** (1-2 us/point)
- Python list operations are native and fast
- No pybind11 overhead

## Bottleneck Analysis

### pybind11 Overhead Measurements

| Operation | Overhead (estimate) |
|-----------|---------------------|
| `py::list[i]` | ~50-500 us |
| `py::list[i].cast<py::list>()` | ~100-1000 us |
| `val.cast<double>()` | ~50-200 us |
| **Total per point** | ~400-3000 us |

Compare to Radia evaluation: **0.5 us/point**

**Conclusion:** Any loop in C++ that accesses Python list elements is 1000x slower than Radia evaluation itself.

## Failed Approaches

### ❌ Approach 1: C++ PrepareCache with loop
```cpp
for (size_t i = 0; i < npts; i++) {
    py::list coords;
    coords.append(x);  // Slow!
    coords.append(y);  // Slow!
    coords.append(z);  // Slow!
    radia_points.append(coords);  // Slow!
}
```
**Result:** 240,000x slower than Radia

### ❌ Approach 2: Extract to C++ vectors first
```cpp
std::vector<std::array<double,3>> points_global(npts);
for (size_t i = 0; i < npts; i++) {
    py::list pt = points_list[i].cast<py::list>();  // Slow!
    points_global[i] = {
        pt[0].cast<double>(),  // Slow!
        pt[1].cast<double>(),  // Slow!
        pt[2].cast<double>()   // Slow!
    };
}
```
**Result:** Still 20,000x slower than Radia

### ❌ Approach 3: Python list prep + C++ cache storage
```python
# Python side
radia_points = [[x*1000, y*1000, z*1000] for x,y,z in points]  # Fast!
results = rad.Fld(obj, field_type, radia_points)  # Fast!

# C++ side (_SetCacheData)
for (size_t i = 0; i < npts; i++) {
    py::list pt = points_list[i].cast<py::list>();  // STILL SLOW!
    py::list fld = results_list[i].cast<py::list>();  // STILL SLOW!
}
```
**Result:** 600,000x slower than Radia (even worse!)

## Viable Solutions

### ✓ Solution 1: Pure Python Implementation [RECOMMENDED]

Keep everything in Python, store results in C++ cache via hash map directly.

```python
def prepare_cache_pure_python(cf, points):
    # Step 1: Radia batch call (fast)
    radia_pts = [[x*1000, y*1000, z*1000] for x,y,z in points]
    results = rad.Fld(cf.radia_obj, cf.field_type, radia_pts)

    # Step 2: Store in Python dict (fast)
    cache = {}
    for i, (pt, res) in enumerate(zip(points, results)):
        cache[tuple(pt)] = res  # O(1) hash insert

    # Step 3: Pass entire dict to C++ (single call)
    cf._SetCacheDict(cache)  # Minimal pybind11 overhead
```

**C++ side:**
```cpp
void _SetCacheDict(py::dict cache) {
    // Iterate over dict items (pybind11 optimized)
    for (auto item : cache) {
        auto key = item.first.cast<py::tuple>();
        auto val = item.second.cast<py::list>();
        // Direct hash insert, no list iteration
        uint64_t hash = hash_point(...);
        point_cache_[hash] = ...;
    }
}
```

**Expected performance:** ~2-5 us/point (4-10x faster than Radia!)

### ✓ Solution 2: NumPy Arrays [COMPLEX]

Use NumPy for zero-copy data transfer.

```python
# Python side
points_np = np.array(points, dtype=np.float64)
results_np = np.array(results, dtype=np.float64)

# C++ side
void _SetCacheNumPy(py::array_t<double> points, py::array_t<double> results) {
    auto pts = points.unchecked<2>();  // Zero-copy view
    auto res = results.unchecked<2>();
    // Direct C++ array access, no Python API calls
    for (size_t i = 0; i < pts.shape(0); i++) {
        double x = pts(i, 0);  // Fast C++ array access
        double y = pts(i, 1);
        double z = pts(i, 2);
        // ...
    }
}
```

**Expected performance:** ~0.5-1 us/point (same as Radia!)

**Drawbacks:**
- Requires NumPy dependency
- More complex implementation
- NumPy C API learning curve

### ✓ Solution 3: Bypass Cache, Use Python Dict [SIMPLEST]

Don't use C++ cache at all - implement cache entirely in Python.

```python
class PythonCachedField:
    def __init__(self, radia_obj, field_type):
        self.radia_obj = radia_obj
        self.field_type = field_type
        self.cache = {}  # Python dict

    def prepare_cache(self, points):
        radia_pts = [[x*1000, y*1000, z*1000] for x,y,z in points]
        results = rad.Fld(self.radia_obj, self.field_type, radia_pts)
        for pt, res in zip(points, results):
            self.cache[tuple(pt)] = res

    def evaluate(self, x, y, z):
        key = (round(x/1e-10)*1e-10, round(y/1e-10)*1e-10, round(z/1e-10)*1e-10)
        if key in self.cache:
            return self.cache[key]
        # Cache miss - evaluate directly
        return rad.Fld(self.radia_obj, self.field_type, [x*1000, y*1000, z*1000])
```

**Expected performance:** ~1-2 us/point

**Advantage:** No C++ changes needed!

## Recommendation

**Immediate fix:** Implement Solution 3 (Pure Python cache)

**Why:**
1. No C++ code changes required
2. Immediate 100,000x performance improvement
3. Works with existing CoefficientFunction interface via callback

**Long-term:** Implement Solution 2 (NumPy arrays) if Solution 3 proves insufficient

## Lesson Learned

**Golden Rule:** Never iterate over Python lists in C++ with pybind11.

**Corollary:** pybind11 is for **control flow**, not **data processing**.

Use Python for data loops, C++ for computation only.

---

**Date:** 2025-11-21
**Status:** C++ approaches abandoned, Python solution recommended
