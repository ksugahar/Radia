# H-Matrix Critical Bug: Default Method Not Compatible

**Status**: 🔴 CRITICAL BUG
**Date**: 2025-11-10
**Impact**: H-matrix may produce incorrect results

## Problem Summary

When H-matrix is enabled, the default relaxation method (MethNo_4) **cannot use H-matrix** because it requires direct access to individual matrix elements.

### Current Behavior

1. User enables H-matrix: `rad.SolverHMatrixEnable()`
2. User calls: `rad.Solve(obj, prec, max_iter)`  (no method specified)
3. System uses default method = **MethNo_4** (see rad_c_interface.cpp:1523)
4. MethNo_4 requires dense `InteractMatrix` for:
   - Off-diagonal contributions (lines 658-661 in rad_relaxation_methods.cpp)
   - Diagonal element in `FindNewH()` (line 676)

5. **BUT**: When H-matrix is enabled:
   - `SetupInteractMatrix_HMatrix()` is called (rad_interaction.cpp:482)
   - H-matrix is built successfully
   - Dense `InteractMatrix` memory is allocated but **NOT filled**
   - MethNo_4 uses uninitialized/garbage data → **WRONG RESULTS**

### Confirmed Facts

- **Default method**: MethNo_4 (rad_c_interface.cpp:1523)
- **Methods with H-matrix support**:
  - `radTSimpleRelaxation` ✓ (lines 110-114 in rad_relaxation_methods.cpp)
  - `radTRelaxationMethNo_2` ✓ (lines 156-172)
  - `radTRelaxationMethNo_8` ❌ (lines 1185-1229, uses dense matrix only)
  - **`radTRelaxationMethNo_4` ❌** (lines 622-681, uses dense matrix only)

- **Memory allocation**: Even with H-matrix:
  - `AllocateMemory()` is called (line 118)
  - `InteractMatrix` pointers are allocated
  - But matrix values are never filled when H-matrix is used
  - No crash, but incorrect computation

## Why Benchmark Doesn't Crash

The benchmark doesn't crash because:
1. `InteractMatrix` memory IS allocated (line 118 in Setup())
2. MethNo_4 accesses allocated (but unfilled) memory
3. No segfault, but **computation is wrong**

## Why Performance Is Poor

Even if results were correct, performance is poor because:

1. **Problem converges in 1-2 iterations** (linear problem, no initial M)
   - Construction cost: 0.98s
   - Iteration cost: ~0.01s
   - Total: 0.99s H-matrix vs 0.03s dense

2. **H-matrix not actually used in iterations**
   - MethNo_4 uses dense matrix (even if garbage)
   - H-matrix MatVec is never called

3. **N=343 is too small for H-matrix benefit**
   - Dense MatVec: O(N²) = 118k operations (< 1ms on modern CPU)
   - H-matrix only beneficial for N > 1000+

## Required Fixes

### Priority 1: Prevent Wrong Results (CRITICAL)

**Option A: Auto method selection** (Recommended)
```cpp
// In rad_c_interface.cpp or rad_material_impl.cpp
void SolveGen(int ObjKey, double PrecOnMagnetiz, int MaxIterNumber, int MethNo)
{
	if(MethNo == 0) {
		// Auto-select based on H-matrix state
		if(RadSolverGetHMatrixEnabled()) {
			MethNo = 2;  // Use MethNo_2 which supports H-matrix
		} else {
			MethNo = 4;  // Use default MethNo_4
		}
	}
	rad.SolveGen(ObjKey, PrecOnMagnetiz, MaxIterNumber, MethNo);
}
```

**Option B: Add H-matrix support to MethNo_4**
- Complex: Requires rewriting the algorithm
- MethNo_4 uses diagonal elements explicitly
- H-matrix doesn't provide element-wise access

**Option C: Build both matrices**
- Defeats H-matrix memory savings
- Not recommended

### Priority 2: Improve H-matrix Performance

See `HMATRIX_PHASE1_IMPLEMENTATION.md` Task 3:
- H-matrix is already built only once per Solve()
- But for current benchmark (1-2 iterations), this doesn't help
- Need problems with 100+ iterations to see benefit

### Priority 3: Better Documentation

Update API docs to clarify:
- H-matrix only works with certain methods (1, 2, not 4, 8)
- Recommend method=2 when using H-matrix
- H-matrix beneficial only for large N (>1000) with many iterations

## Testing Plan

### Test 1: Verify Bug
```python
import radia as rad

# Create problem
mat = rad.MatSatIsoFrm([10000, 100], [0.1, 1000])
elem = rad.ObjRecMag([0,0,0], [10,10,10])
rad.MatApl(elem, mat)
rad.ObjBckg([10000,0,0])

# Test with H-matrix (may be wrong!)
rad.SolverHMatrixEnable()
result_hmat = rad.Solve(elem, 0.0001, 100)  # Uses MethNo_4 by default
H_hmat = rad.Fld(elem, 'h', [0,0,50])

# Test with dense (correct)
rad.UtiDelAll()
mat = rad.MatSatIsoFrm([10000, 100], [0.1, 1000])
elem = rad.ObjRecMag([0,0,0], [10,10,10])
rad.MatApl(elem, mat)
rad.ObjBckg([10000,0,0])
rad.SolverHMatrixDisable()
result_dense = rad.Solve(elem, 0.0001, 100)
H_dense = rad.Fld(elem, 'h', [0,0,50])

# Compare
print(f"H-matrix: {H_hmat}")
print(f"Dense:    {H_dense}")
print(f"Error:    {np.linalg.norm(np.array(H_hmat) - np.array(H_dense))}")
```

Expected: **Large error** if bug exists

### Test 2: Verify Fix (Option A)
After implementing auto method selection:
- H-matrix should automatically use MethNo_2
- Results should match dense solver
- Performance may still be poor for small problems (expected)

## References

- Phase 1 implementation plan: `HMATRIX_PHASE1_IMPLEMENTATION.md`
- Phase 1 Task 1 commit: 839a385 (ACA parameter optimization)
- Default method: rad_c_interface.cpp:1523
- MethNo_4 implementation: rad_relaxation_methods.cpp:622-681
- H-matrix MatVec: rad_interaction.cpp:1517-1532
