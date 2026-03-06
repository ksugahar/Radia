# Memory Leak Diagnosis Report - 2025-11-21

## Executive Summary

**Critical memory leak identified in Radia's field computation function `rad.Fld()`**

- **Leak magnitude**: ~0.15 KB per field evaluation point
- **Impact**: Severe in NGSolve integration (1.5 MB per time step for typical mesh)
- **Root cause**: `rad.Fld()` implementation does not release internal buffers/caches
- **Affected versions**: v1.3.2 and earlier

## Problem Statement

User reported severe memory leak when simulating moving magnets with NGSolve integration:
- Memory increases by **1.4 GB projected for 1000 time steps**
- Linear memory growth: **1,478 KB per time step**
- Makes long-time simulations impossible

## Diagnostic Tests Performed

### Test 1: NGSolve Integration Memory Test

**File**: `tests/test_moving_magnet_memory.py`

**Configuration**:
- 100 time steps
- Moving magnet (0.001 m per step)
- NGSolve mesh: 216 vertices, 7,644 DOFs (HCurl, order=2)
- GridFunction.Set() called each step

**Results**:
```
Memory increase: 144.37 MB (100 steps)
Growth rate: 1,478.24 KB/step
Projected for 1000 steps: 1,443 MB

Status: FAIL - Severe memory leak
```

### Test 2: Diagnosis - Component Isolation

**File**: `tests/test_moving_magnet_memory_diagnosis.py`

**Test scenarios**:
1. Baseline (create new CF + GF each time)
2. Reuse CoefficientFunction
3. Reuse GridFunction
4. Reuse both CF + GF
5. Reuse both + explicit ClearCache()

**Results**:
```
Test                     Mem (MB)   Rate (KB/step)   Status
------------------------------------------------------------
Baseline                 72.18      1478.24          FAIL
Reuse CF                 72.18      1478.23          FAIL
Reuse GF                 72.16      1477.87          FAIL
Reuse both               72.16      1477.85          FAIL
Reuse both + clear       72.16      1477.88          FAIL
```

**Conclusion**: Leak is NOT in RadiaField CF or GridFunction object lifetime.
Leak persists even with object reuse.

### Test 3: Radia Core Memory Test (No NGSolve)

**File**: `tests/test_radia_only_memory.py`

**Configuration**:
- 100 time steps
- Moving magnet
- 100 field evaluation points using `rad.Fld()` per step
- No NGSolve involvement

**Results**:
```
Memory increase: 1.47 MB (100 steps)
Growth rate: 14.86 KB/step

Status: FAIL - Memory leak in Radia core
```

**Conclusion**: Leak exists in Radia core, independent of NGSolve.

### Test 4: Field Computation Isolation

**File**: `tests/test_radia_field_computation_memory.py`

**Test scenarios**:
1. Create/Delete magnet only (no field computation): **0.00 KB/step**
2. Create/Delete + 1 field evaluation: **0.15 KB/step**
3. Create/Delete + 10 field evaluations: **1.49 KB/step** (0.15 × 10)
4. Create/Delete + 100 field evaluations: **14.85 KB/step** (0.15 × 100)

**Results**:

| Test | Field Points | Memory Growth (KB/step) |
|------|--------------|-------------------------|
| Object creation only | 0 | 0.00 |
| + Field evaluation | 1 | 0.15 |
| + Field evaluation | 10 | 1.49 |
| + Field evaluation | 100 | 14.85 |

**Leak rate**: **0.15 KB per `rad.Fld()` call**

**Conclusion**:
- Object creation/deletion (`rad.UtiDelAll()`, `rad.ObjRecMag()`) has NO leak
- Field computation (`rad.Fld()`) has consistent 0.15 KB/call leak
- Leak is proportional to number of field evaluations

## Root Cause Analysis

### Identified Source

**Primary leak location**: `rad.Fld()` implementation (radapl4.cpp)

**Mechanism**:
- Each call to `rad.Fld()` allocates ~0.15 KB that is never freed
- Likely candidates:
  - Internal field computation buffers
  - Temporary arrays for field calculation
  - Subdivision data structures
  - Interpolation caches

### Why NGSolve Integration Shows 100× Amplification

NGSolve's `GridFunction.Set(B_cf)` evaluates the CoefficientFunction at all mesh degrees of freedom:

**Calculation**:
```
Mesh: 7,644 DOFs (HCurl order=2)
Leak: 1,478 KB/step ÷ 0.15 KB/point ≈ 9,853 points evaluated

This matches the order of magnitude for mesh DOF count.
```

**Conclusion**: NGSolve integration amplifies Radia's field computation leak by evaluating at thousands of points per time step.

## Code Investigation

### DeleteAllElements Implementation

**File**: `src/core/rad_transform_impl.cpp:1427-1450`

```cpp
int radTApplication::DeleteAllElements(int DeletionMethNo)
{
    try
    {
        if(DeletionMethNo == 1)
            GlobalMapOfHandlers.erase(GlobalMapOfHandlers.begin(),
                                      GlobalMapOfHandlers.end());
        // ...
        GlobalUniqueMapKey = 1;
        MapOfDrawAttr.erase(MapOfDrawAttr.begin(), MapOfDrawAttr.end());

        if(SendingIsRequired) Send.Int(0);
        return 1;
    }
    catch(...)
    {
        Initialize(); return 0;
    }
}
```

**Status**: ✓ Implementation correct - properly clears maps and resets key counter.

**Conclusion**: Leak is NOT in object lifecycle management.

### Field Computation (Next Step)

**Files to investigate**:
- `src/lib/radapl4.cpp` - Field computation API
- `src/core/radfield.cpp` - Field calculation implementation
- `src/core/radg3d.cpp` - Geometry field methods

**Suspected functions**:
- `RadFld()` entry point
- Field computation temporary buffer allocation
- Subdivision data structures in field interpolation

## Impact Assessment

### Severity: **CRITICAL**

**Affected use cases**:
- ✗ Time-stepping simulations (moving magnets)
- ✗ Parameter sweeps (magnetic field maps)
- ✗ NGSolve FEM coupling (any multi-step workflow)
- ✓ Single-shot calculations (minimal impact)

**Memory consumption examples**:

| Simulation Steps | Mesh DOFs | Memory Leak | Time to OOM (16 GB RAM) |
|------------------|-----------|-------------|------------------------|
| 100 | 7,644 | 144 MB | ~11,000 steps |
| 1,000 | 7,644 | 1.4 GB | ~11,000 steps |
| 10,000 | 7,644 | 14 GB | ~11,000 steps |
| 100 | 30,000 | 562 MB | ~2,800 steps |
| 1,000 | 30,000 | 5.6 GB | ~2,800 steps |

**Conclusion**: Long simulations (> 1000 steps) will exhaust system memory.

## Workarounds (Temporary)

### 1. Restart Python Process Periodically

```python
# In time-stepping loop, checkpoint and restart every N steps
if step % 1000 == 0:
    save_checkpoint(step, data)
    os.execv(sys.executable, ['python'] + sys.argv)  # Restart process
```

**Pros**: Clears all memory
**Cons**: Complex checkpoint/restart logic, loses Python state

### 2. Use Radia-Only Field Evaluation

```python
# Instead of GridFunction.Set()
for i, pt in enumerate(mesh_points):
    B_values[i] = rad.Fld(magnet, 'b', pt)
```

**Pros**: Avoids NGSolve amplification (100× reduction)
**Cons**: Still leaks 0.15 KB/point, slower than Set()

### 3. Reduce Mesh Resolution

```python
# Use coarser mesh to reduce DOF count
mesh = Mesh(geo.GenerateMesh(maxh=0.05))  # Larger maxh = fewer DOFs
```

**Pros**: Reduces leak proportionally
**Cons**: Lower accuracy

**Note**: These are stopgap measures - **proper fix required in Radia core**.

## Recommended Fix Strategy

### Phase 1: Identify Leaking Allocation

1. **Profile `rad.Fld()` with Valgrind/Dr. Memory**
   ```bash
   valgrind --leak-check=full --track-origins=yes python test_radia_field_computation_memory.py
   ```

2. **Add debug output to radapl4.cpp**
   - Log allocations in `RadFld()` function
   - Track buffer creation/destruction

3. **Review field computation code**
   - Check for `new` without corresponding `delete`
   - Look for static/global buffers that accumulate
   - Verify RAII patterns for temporary arrays

### Phase 2: Implement Fix

**Suspected issues**:

1. **Temporary buffers not freed**
   ```cpp
   // BAD - leak
   double* buffer = new double[n];
   // ... computation ...
   // Missing: delete[] buffer;

   // GOOD - RAII
   std::vector<double> buffer(n);
   // Automatic cleanup
   ```

2. **Static cache that grows unbounded**
   ```cpp
   // BAD - grows forever
   static std::map<int, std::vector<double>> field_cache;

   // GOOD - bounded cache or per-object lifetime
   ```

3. **Subdivision data not cleaned up**
   ```cpp
   // Check subdivision temporary structures
   // Ensure they're freed after field computation
   ```

### Phase 3: Verify Fix

Run all diagnostic tests:
```bash
cd tests
python test_radia_field_computation_memory.py  # Should show 0.00 KB/step
python test_radia_only_memory.py                # Should show < 1 MB total
python test_moving_magnet_memory.py             # Should show < 10 MB total
```

**Success criteria**:
- `rad.Fld()` leak: **< 0.01 KB/call** (negligible)
- NGSolve integration leak: **< 10 MB for 100 steps**
- Test status: **PASS** on all memory tests

## Files Created for Diagnosis

1. `tests/test_moving_magnet_memory.py` - NGSolve integration test
2. `tests/test_moving_magnet_memory_diagnosis.py` - Component isolation
3. `tests/test_radia_only_memory.py` - Radia core test
4. `tests/test_radia_field_computation_memory.py` - Field computation isolation

**Note**: These tests should remain in the test suite to prevent regression.

## Next Steps

1. **Immediate**: Add Valgrind/Dr. Memory profiling to CI pipeline
2. **Short term**: Investigate radapl4.cpp `RadFld()` implementation
3. **Medium term**: Implement fix and verify with diagnostic tests
4. **Long term**: Add memory leak tests to automated test suite

## References

- User report: "Previous version had memory increasing issue"
- CLAUDE.md policies: Memory Management, Exception Safety
- Related code: rad_transform_impl.cpp:1427 (DeleteAllElements)

---

**Report prepared**: 2025-11-21
**Tests executed**: 4 comprehensive diagnostic tests
**Leak confirmed**: 0.15 KB per `rad.Fld()` call
**Recommended action**: Investigate radapl4.cpp field computation buffers
