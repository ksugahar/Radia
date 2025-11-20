# Batch Evaluation for H-Matrix Acceleration in radia_ngsolve

## Problem Statement

Current implementation of `GridFunction.Set(coefficient_function)` in radia_ngsolve:

### How it works now:
1. NGSolve calls `RadiaFieldCF::Evaluate()` for **each mesh element**
2. Each call evaluates ~10-20 points (one element's integration points)
3. Even with H-matrix enabled, each call has overhead >> computation time
4. **H-matrix speedup is NOT realized**

### Performance impact:
- For mesh with 500 elements:
  - 500 calls to `Evaluate()`
  - Each call evaluates 15 points
  - Total: 7500 points evaluated in 500 batches
  - H-matrix overhead: ~500 × (setup cost)
  - **Result: No speedup, possibly slower than direct evaluation**

## Proposed Solution

### Batch Evaluation with PrepareCache():

```cpp
// User code:
rad.SetHMatrixFieldEval(1, 1e-6)  // Enable H-matrix
cf = radia_ngsolve.RadiaField(magnet, 'b')
cf.PrepareCache(mesh)  // NEW: Pre-compute all values
gf.Set(cf)  // Fast: returns cached values
```

### Implementation plan:

#### 1. Add cache to RadiaFieldCF class:
```cpp
class RadiaFieldCF : public CoefficientFunction {
    // ...existing members...
    
    // Batch evaluation cache
    std::map<std::tuple<double,double,double>, std::array<double,3>> point_cache_;
    bool use_cache_;
    
public:
    void PrepareCache(py::object py_mesh);
    void ClearCache();
};
```

#### 2. Implement PrepareCache():
```cpp
void RadiaFieldCF::PrepareCache(py::object py_mesh) {
    // Step 1: Collect ALL integration points from mesh
    std::vector<std::array<double,3>> all_points;
    
    // Iterate over all mesh elements
    // For each element:
    //   - Get integration rule
    //   - Map integration points to global coordinates
    //   - Add to all_points
    
    // Step 2: Batch evaluate using rad.FldBatch()
    py::module_ rad = py::module_::import("radia");
    py::list points_list;
    for (auto& pt : all_points) {
        py::list coords;
        coords.append(pt[0] * 1000.0);  // m -> mm
        coords.append(pt[1] * 1000.0);
        coords.append(pt[2] * 1000.0);
        points_list.append(coords);
    }
    
    // Single batch call - full H-matrix speedup!
    int use_hmatrix_flag = use_hmatrix.is_none() ? -1 : 
                           use_hmatrix.cast<int>();
    py::object results = rad.attr("FldBatch")(
        radia_obj, field_type, points_list, use_hmatrix_flag
    );
    
    // Step 3: Store in cache
    py::list results_list = results.cast<py::list>();
    for (size_t i = 0; i < all_points.size(); i++) {
        auto& pt = all_points[i];
        py::list field = results_list[i].cast<py::list>();
        
        std::array<double,3> value = {
            field[0].cast<double>(),
            field[1].cast<double>(),
            field[2].cast<double>()
        };
        
        auto key = std::make_tuple(pt[0], pt[1], pt[2]);
        point_cache_[key] = value;
    }
    
    use_cache_ = true;
}
```

#### 3. Modify Evaluate() to use cache:
```cpp
void RadiaFieldCF::Evaluate(const BaseMappedIntegrationRule &mir,
                            BareSliceMatrix<> result) const {
    if (use_cache_) {
        // Fast path: return cached values
        for (size_t i = 0; i < mir.Size(); i++) {
            auto pt = mir[i].GetPoint();
            auto key = std::make_tuple(pt[0], pt[1], pt[2]);
            
            auto it = point_cache_.find(key);
            if (it != point_cache_.end()) {
                result(i, 0) = it->second[0];
                result(i, 1) = it->second[1];
                result(i, 2) = it->second[2];
            } else {
                // Point not in cache - evaluate directly
                // (shouldn't happen if PrepareCache was called correctly)
                EvaluateDirect(mir[i], result, i);
            }
        }
    } else {
        // Standard path: batch evaluation as before
        EvaluateBatch(mir, result);
    }
}
```

## Expected Performance

### Current (element-by-element):
- Mesh: 500 elements × 15 points = 7500 total points
- Calls to Radia: 500 calls
- H-matrix setup overhead: 500× (wasted)
- **Time: ~1000 ms** (no H-matrix benefit)

### Optimized (batch with PrepareCache):
- Mesh: same 7500 points
- Calls to Radia: **1 call** (PrepareCache)
- H-matrix setup overhead: 1× (efficient!)
- **Time: ~50 ms** (full H-matrix benefit)
- **Speedup: 20x**

### Scalability:
For larger problems (N >> 1000):
- Element-by-element: O(N_elem) × overhead → No speedup
- Batch: O(1) × overhead → Full H-matrix speedup (O(N log N))
- **Expected speedup: 50-100x for N > 5000**

## Usage Example

```python
import radia as rad
from ngsolve import *
from netgen.occ import *
import radia_ngsolve

# Create Radia geometry (N=1000 elements)
rad.FldUnits('m')
magnet = create_large_magnet()  # 1000+ elements

# Create NGSolve mesh
mesh = Mesh(...)
fes = HCurl(mesh, order=2)
gf = GridFunction(fes)

# Enable H-matrix
rad.SetHMatrixFieldEval(1, 1e-6)

# Create CoefficientFunction
B_cf = radia_ngsolve.RadiaField(magnet, 'b')

# SLOW way (current):
# gf.Set(B_cf)  # 500 element calls, no H-matrix benefit, ~1000 ms

# FAST way (proposed):
B_cf.PrepareCache(mesh)  # Single batch evaluation, ~50 ms
gf.Set(B_cf)  # Returns cached values, ~1 ms
# Total: ~51 ms (20x faster)
```

## Implementation Status

- [x] Problem identified and analyzed
- [x] Solution proposed
- [ ] C++ implementation (PrepareCache)
- [ ] Testing and benchmarking
- [ ] Documentation
- [ ] Integration with existing code

## Alternative Approaches Considered

### 1. Python-only solution:
**Problem**: Can't intercept NGSolve's internal GridFunction.Set() calls
**Verdict**: Not feasible without C++ changes

### 2. Custom GridFunction setter:
**Problem**: Would require reimplementing L² projection
**Verdict**: Too complex, error-prone

### 3. Batch evaluation in Evaluate():
**Current**: Already does batch evaluation per element
**Problem**: Still called N_elem times by NGSolve
**Verdict**: Not sufficient

## Recommendation

Implement PrepareCache() in C++ as proposed. This is:
- ✅ Clean API (user explicitly enables optimization)
- ✅ No breaking changes (optional feature)
- ✅ Maximum performance gain
- ✅ Works with existing NGSolve infrastructure

---

**Status**: Proposal  
**Priority**: High (enables H-matrix speedup in coupled simulations)  
**Effort**: ~1-2 days implementation + testing

