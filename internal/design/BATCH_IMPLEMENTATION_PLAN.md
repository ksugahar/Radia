# Batch Evaluation Implementation Plan

## Overview

Implement `PrepareCache()` method in `RadiaFieldCF` class to enable H-matrix acceleration by batching all integration point evaluations.

**Target file**: `src/python/radia_ngsolve.cpp`

**Estimated time**: 4-5 hours total implementation + testing

---

## Phase 1: Add Cache Infrastructure

### 1.1 Add Cache Members to RadiaFieldCF Class

**Location**: `src/python/radia_ngsolve.cpp`, `RadiaFieldCF` class definition

**Add these member variables**:

```cpp
class RadiaFieldCF : public CoefficientFunction {
private:
    // Existing members
    PyObjectWrapper radia_obj_;
    std::string field_type_;

    // NEW: Batch evaluation cache
    std::unordered_map<uint64_t, std::array<double,3>> point_cache_;
    bool use_cache_;
    double cache_tolerance_;  // For hash quantization (e.g., 1e-10)
    size_t cache_hits_;
    size_t cache_misses_;

    // Hash function for 3D points
    uint64_t HashPoint(double x, double y, double z) const {
        // Quantize to tolerance grid
        int64_t ix = static_cast<int64_t>(x / cache_tolerance_ + 0.5);
        int64_t iy = static_cast<int64_t>(y / cache_tolerance_ + 0.5);
        int64_t iz = static_cast<int64_t>(z / cache_tolerance_ + 0.5);

        // Simple hash combination
        uint64_t hx = static_cast<uint64_t>(ix) * 73856093;
        uint64_t hy = static_cast<uint64_t>(iy) * 19349663;
        uint64_t hz = static_cast<uint64_t>(iz) * 83492791;
        return hx ^ hy ^ hz;
    }
};
```

### 1.2 Initialize Cache in Constructor

**Modify constructor**:

```cpp
RadiaFieldCF(PyObjectWrapper radia_obj, std::string field_type)
    : CoefficientFunction(3, false),  // 3D vector field
      radia_obj_(radia_obj),
      field_type_(field_type),
      use_cache_(false),              // NEW
      cache_tolerance_(1e-10),        // NEW
      cache_hits_(0),                 // NEW
      cache_misses_(0)                // NEW
{
    SetDimensions(Array<int>({3}));
}
```

---

## Phase 2: Implement PrepareCache() Method

### 2.1 PrepareCache() Implementation

**Add this method to RadiaFieldCF class**:

```cpp
void PrepareCache(shared_ptr<ngcomp::MeshAccess> mesh, int integration_order = -1) {
    py::gil_scoped_acquire acquire;

    // Clear any existing cache
    point_cache_.clear();
    cache_hits_ = 0;
    cache_misses_ = 0;

    // Step 1: Collect all integration points from all elements
    std::vector<std::array<double,3>> all_points;

    for (auto el : mesh->Elements(VOL)) {
        auto& trafo = mesh->GetTrafo(el, false);

        // Determine integration rule
        int order = (integration_order > 0) ? integration_order : 2 * mesh->GetFE(el, VOL).Order();
        const IntegrationRule& ir = SelectIntegrationRule(trafo.GetElementType(), order);

        // Get physical coordinates of integration points
        for (int i = 0; i < ir.Size(); i++) {
            MappedIntegrationPoint<3,3> mip(ir[i], trafo);
            auto pt = mip.GetPoint();
            all_points.push_back({pt(0), pt(1), pt(2)});
        }
    }

    if (all_points.empty()) {
        std::cout << "[PrepareCache] No integration points found" << std::endl;
        return;
    }

    std::cout << "[PrepareCache] Collected " << all_points.size() << " integration points" << std::endl;

    // Step 2: Single batch evaluation via Radia
    py::list points_list;
    for (const auto& pt : all_points) {
        py::list coord;
        coord.append(pt[0]);
        coord.append(pt[1]);
        coord.append(pt[2]);
        points_list.append(coord);
    }

    // Call rad.FldBatch()
    py::object radia_mod = py::module_::import("radia");
    py::object fld_batch = radia_mod.attr("FldBatch");

    // FldBatch(obj, field_type, points, use_hmatrix=1)
    py::object result = fld_batch(radia_obj_.get(), field_type_, points_list, 1);

    // Step 3: Store in cache
    py::list result_list = result.cast<py::list>();

    for (size_t i = 0; i < all_points.size(); i++) {
        const auto& pt = all_points[i];
        py::list field = result_list[i].cast<py::list>();

        std::array<double,3> field_val = {
            field[0].cast<double>(),
            field[1].cast<double>(),
            field[2].cast<double>()
        };

        uint64_t hash = HashPoint(pt[0], pt[1], pt[2]);
        point_cache_[hash] = field_val;
    }

    use_cache_ = true;
    std::cout << "[PrepareCache] Cached " << point_cache_.size()
              << " unique points for field type: " << field_type_ << std::endl;
}
```

---

## Phase 3: Modify Evaluate() to Use Cache

### 3.1 Add EvaluateFromCache() Helper

**Add private helper method**:

```cpp
void EvaluateFromCache(const BaseMappedIntegrationRule& mir,
                       BareSliceMatrix<> result) const {
    for (size_t i = 0; i < mir.Size(); i++) {
        auto pt = mir[i].GetPoint();
        uint64_t hash = HashPoint(pt(0), pt(1), pt(2));

        auto it = point_cache_.find(hash);
        if (it != point_cache_.end()) {
            // Cache hit
            result(0, i) = it->second[0];
            result(1, i) = it->second[1];
            result(2, i) = it->second[2];
            cache_hits_++;
        } else {
            // Cache miss - fallback to direct evaluation
            cache_misses_++;
            EvaluateSinglePoint(pt(0), pt(1), pt(2), result.Col(i));
        }
    }
}
```

### 3.2 Modify Main Evaluate() Method

**Update existing Evaluate() method**:

```cpp
virtual void Evaluate(const BaseMappedIntegrationRule& mir,
                      BareSliceMatrix<> result) const override {
    // Fast path: use cache if available
    if (use_cache_) {
        EvaluateFromCache(mir, result);
        return;
    }

    // Standard path: batch evaluation per element (existing code)
    EvaluateBatch(mir, result);
}
```

### 3.3 Add Cache Statistics Method

```cpp
void PrintCacheStats() const {
    if (!use_cache_) {
        std::cout << "[Cache] Not enabled" << std::endl;
        return;
    }

    std::cout << "[Cache] Statistics:" << std::endl;
    std::cout << "  Entries: " << point_cache_.size() << std::endl;
    std::cout << "  Hits: " << cache_hits_ << std::endl;
    std::cout << "  Misses: " << cache_misses_ << std::endl;

    if (cache_hits_ + cache_misses_ > 0) {
        double hit_rate = 100.0 * cache_hits_ / (cache_hits_ + cache_misses_);
        std::cout << "  Hit rate: " << hit_rate << "%" << std::endl;
    }
}
```

---

## Phase 4: Export to Python via pybind11

### 4.1 Add Python Bindings

**Location**: End of `src/python/radia_ngsolve.cpp`, in `PYBIND11_MODULE` section

**Add these bindings**:

```cpp
PYBIND11_MODULE(radia_ngsolve, m) {
    m.doc() = "Radia-NGSolve integration module";

    // Existing RadiaField binding
    py::class_<RadiaFieldCF, shared_ptr<RadiaFieldCF>, CoefficientFunction>(m, "RadiaField")
        .def(py::init<PyObjectWrapper, std::string>(),
             py::arg("radia_obj"),
             py::arg("field_type") = "b",
             "Create Radia field CoefficientFunction")

        // NEW: PrepareCache method
        .def("PrepareCache", &RadiaFieldCF::PrepareCache,
             py::arg("mesh"),
             py::arg("integration_order") = -1,
             "Pre-compute field at all integration points for H-matrix acceleration.\n"
             "\n"
             "Parameters:\n"
             "  mesh: NGSolve mesh object\n"
             "  integration_order: Integration rule order (default: 2*element_order)\n"
             "\n"
             "Example:\n"
             "  B_cf = radia_ngsolve.RadiaField(magnet, 'b')\n"
             "  B_cf.PrepareCache(mesh)\n"
             "  B_gf.Set(B_cf)  # Fast: uses cached values\n")

        // NEW: Cache statistics
        .def("PrintCacheStats", &RadiaFieldCF::PrintCacheStats,
             "Print cache hit/miss statistics")

        // NEW: Clear cache
        .def("ClearCache", [](RadiaFieldCF& self) {
            self.use_cache_ = false;
            self.point_cache_.clear();
            self.cache_hits_ = 0;
            self.cache_misses_ = 0;
        }, "Clear cached field values");
}
```

---

## Phase 5: Testing and Validation

### 5.1 Create Test Script

**File**: `tests/test_batch_evaluation.py`

```python
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

import radia as rad
import radia_ngsolve
from ngsolve import *
from netgen.occ import *
import time
import numpy as np

print("="*70)
print("Batch Evaluation Performance Test")
print("="*70)

# Create magnet array (N=125 elements)
rad.UtiDelAll()
n = 5
elements = []
for i in range(n):
    for j in range(n):
        for k in range(n):
            x = (i - n/2 + 0.5) * 0.02
            y = (j - n/2 + 0.5) * 0.02
            z = (k - n/2 + 0.5) * 0.02
            elem = rad.ObjRecMag([x, y, z], [0.02, 0.02, 0.02], [0, 0, 1.2])
            elements.append(elem)

magnet = rad.ObjCnt(elements)
print(f"\n[Setup] Created magnet array: {n**3} elements")

# Create mesh
box = Box((0.015, 0.015, 0.015), (0.063, 0.063, 0.063))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.010))
print(f"[Setup] Mesh: {mesh.ne} elements, {mesh.nv} vertices")

# Enable H-matrix
rad.SetHMatrixFieldEval(1, 1e-6)

# Test 1: Standard evaluation (element-by-element)
print("\n" + "="*70)
print("TEST 1: Standard GridFunction.Set() (element-by-element)")
print("="*70)

fes = HDiv(mesh, order=2)
B_gf_standard = GridFunction(fes)
B_cf_standard = radia_ngsolve.RadiaField(magnet, 'b')

t0 = time.time()
B_gf_standard.Set(B_cf_standard)
t_standard = time.time() - t0

print(f"  Time: {t_standard*1000:.1f} ms")

# Test 2: Batch evaluation with PrepareCache()
print("\n" + "="*70)
print("TEST 2: Optimized GridFunction.Set() with PrepareCache()")
print("="*70)

B_gf_batch = GridFunction(fes)
B_cf_batch = radia_ngsolve.RadiaField(magnet, 'b')

t0 = time.time()
B_cf_batch.PrepareCache(mesh)
t_cache_prep = time.time() - t0

t0 = time.time()
B_gf_batch.Set(B_cf_batch)
t_batch_set = time.time() - t0

t_batch_total = t_cache_prep + t_batch_set

print(f"  PrepareCache time: {t_cache_prep*1000:.1f} ms")
print(f"  Set() time: {t_batch_set*1000:.1f} ms")
print(f"  Total time: {t_batch_total*1000:.1f} ms")

B_cf_batch.PrintCacheStats()

# Verify accuracy
test_points = [
    (0.030, 0.020, 0.040),
    (0.040, 0.040, 0.050),
    (0.050, 0.030, 0.060),
]

print("\n" + "="*70)
print("ACCURACY VERIFICATION")
print("="*70)

errors = []
for pt in test_points:
    mip = mesh(*pt)
    B_std = np.array(B_gf_standard(mip))
    B_bat = np.array(B_gf_batch(mip))

    error = np.linalg.norm(B_bat - B_std)
    B_norm = np.linalg.norm(B_std)
    rel_error = error / B_norm * 100 if B_norm > 0 else 0
    errors.append(rel_error)

    print(f"  Point {pt}: Error = {rel_error:.6f}%")

# Performance summary
print("\n" + "="*70)
print("PERFORMANCE SUMMARY")
print("="*70)

speedup = t_standard / t_batch_total

print(f"  Standard method: {t_standard*1000:.1f} ms")
print(f"  Batch method: {t_batch_total*1000:.1f} ms")
print(f"  Speedup: {speedup:.1f}x")
print(f"  Mean accuracy error: {np.mean(errors):.6f}%")

if speedup > 2.0 and np.mean(errors) < 1.0:
    print("\n[SUCCESS] Batch evaluation provides significant speedup with good accuracy!")
else:
    print(f"\n[INFO] Results recorded")
```

### 5.2 Test with Different Problem Sizes

Create additional test with N=1000 elements to verify scalability.

---

## Implementation Timeline

| Phase | Task | Estimated Time | Files Modified |
|-------|------|----------------|----------------|
| 1 | Add cache infrastructure | 30 min | `radia_ngsolve.cpp` |
| 2 | Implement PrepareCache() | 1.5 hours | `radia_ngsolve.cpp` |
| 3 | Modify Evaluate() | 1 hour | `radia_ngsolve.cpp` |
| 4 | Add Python bindings | 30 min | `radia_ngsolve.cpp` |
| 5 | Testing and validation | 1.5 hours | `test_batch_evaluation.py` |
| **Total** | | **~5 hours** | |

---

## Expected Results

**Performance (N=125 elements, ~7500 integration points)**:
- Standard: ~1000 ms (no H-matrix benefit)
- Batch: ~50 ms (full H-matrix benefit)
- **Expected speedup: 20x**

**Performance (N=1000 elements, ~60000 integration points)**:
- Standard: ~10000 ms
- Batch: ~100 ms
- **Expected speedup: 100x**

**Accuracy**:
- Should match standard method exactly (same Radia evaluations)
- Expected difference: < 1e-12 (numerical precision)

---

## Completion Criteria

1. ✅ Code compiles without errors
2. ✅ PrepareCache() successfully collects all integration points
3. ✅ Single FldBatch() call completes with H-matrix
4. ✅ Cache hit rate > 99%
5. ✅ Speedup > 10x for N=125 problem
6. ✅ Accuracy error < 1e-6% vs standard method
7. ✅ Python API works as documented

---

**Last Updated**: 2025-11-20
