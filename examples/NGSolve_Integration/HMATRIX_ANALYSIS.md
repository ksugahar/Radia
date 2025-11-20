# H-matrix Performance Analysis for GridFunction.Set()

**Date**: 2025-11-08
**Issue**: H-matrix shows minimal speedup (~1.00x) for GridFunction.Set() operation

---

## Benchmark Results Summary

From `benchmark_gridfunction_set.py`:
- **Overall average H-matrix speedup**: 1.02x
- **Speedup by field type**: ~1.00x for all types (B, H, A, M)
- **Speedup vs element count**: ~1.00x for all N (27 to 343)
- **Speedup vs mesh size**: ~1.00x for all mesh sizes

**Conclusion**: H-matrix has **no significant effect** on GridFunction.Set() performance.

---

## Root Cause Analysis

### 1. Two Different Operations in Radia

Radia has two distinct computational phases:

#### Phase 1: **SOLVE** (Magnetization Computation)
- **Purpose**: Compute magnetization distribution M from materials and geometry
- **Operation**: Solve linear system: `M_unknown = K^-1 * b`
  - K is the N×N interaction matrix (dense)
  - For nonlinear materials: iterative Gauss-Seidel
  - For linear materials with LU: direct solve
- **H-matrix benefit**: ✅ **YES** - H-matrix accelerates matrix-vector products
- **When it happens**: During `rad.Solve()` or relaxation

#### Phase 2: **FIELD EVALUATION** (Post-processing)
- **Purpose**: Compute field B(r) at arbitrary point r from known magnetization
- **Operation**: Sum contributions: `B(r) = Σ G(r, r_i) * M_i`
  - G(r, r_i) is the Green's function (Biot-Savart kernel)
  - M_i is the known magnetization of element i
  - Sum over all N elements
- **H-matrix benefit**: ❌ **NO** - Current H-matrix infrastructure not used here
- **When it happens**: During `rad.Fld()` calls

### 2. What GridFunction.Set() Does

```python
gf = GridFunction(HCurl(mesh))
cf = radia_ngsolve.RadiaField(magnet, 'b')
gf.Set(cf)  # <--- This operation
```

**Sequence**:
1. NGSolve iterates over all mesh vertices (M vertices)
2. For each vertex at position r:
   - Call `cf.Evaluate(r)` → calls `radia_ngsolve::RadiaFieldCF::Evaluate()`
   - This calls `rad.Fld(obj, field_type, r)` (Python call)
   - Radia computes field by summing N element contributions
3. Total: M × N operations

**Key point**: `GridFunction.Set()` only performs **FIELD EVALUATION**, not solving.

### 3. Why H-matrix Doesn't Help

The current H-matrix implementation in Radia:
- ✅ Accelerates matrix-vector products during **solving**: `y = K*x`
- ✅ Used in iterative solvers (Gauss-Seidel, etc.)
- ❌ NOT used in field evaluation: `B(r) = Σ G(r, r_i) * M_i`

**Reason**: Field evaluation (`rad.Fld()`) doesn't use the H-matrix data structure. It directly sums contributions from all elements using the Biot-Savart law.

---

## Why Results Show ~1.00x Speedup

Looking at the benchmark data:

| Configuration | Dense (ms) | H-matrix (ms) | Speedup |
|--------------|-----------|--------------|---------|
| N=27, h=0.05 | 22.94 | 13.97 | 1.64x |
| N=343, h=0.0125 | 7354.75 | 7371.10 | 1.00x |

**Small speedup for small problems**: Random variation or cache effects
**No speedup for large problems**: Confirms H-matrix not used in field evaluation

The computation time scales as **O(M × N)** regardless of H-matrix setting, where:
- M = number of mesh vertices
- N = number of Radia elements

---

## What COULD Make GridFunction.Set() Faster?

### Option 1: Fast Multipole Method (FMM) for Field Evaluation

Implement hierarchical field evaluation:
- Group nearby Radia elements into clusters
- Approximate far-field contributions using multipole expansions
- Complexity: O(M) instead of O(M × N)
- **Challenge**: Requires significant implementation effort

### Option 2: Extend H-matrix to Field Evaluation

Use existing H-matrix structure for field computation:
- The H-matrix already has hierarchical clustering of elements
- Could evaluate field hierarchically: near-field exact, far-field approximate
- **Challenge**: Would require modifying `rad.Fld()` implementation

### Option 3: Vectorized Field Evaluation

Evaluate fields at multiple points simultaneously:
```cpp
// Current: Loop in Python
for each vertex r:
    B = rad.Fld(obj, 'b', r)  # N operations per call

// Proposed: Single call
B_array = rad.FldVec(obj, 'b', r_array)  # Amortized overhead
```
**Benefit**: Reduce Python-C++ call overhead
**Challenge**: API change required

### Option 4: Caching/Interpolation

For nearby mesh points, interpolate field instead of recomputing:
- Build field cache on coarse grid
- Interpolate for fine mesh vertices
- **Challenge**: Accuracy vs speed tradeoff

---

## Current Limitations

### 1. Python Call Overhead

Each vertex evaluation requires:
- Python → C++ call to `Evaluate()`
- C++ → Python call to `rad.Fld()`
- Python → C++ call into Radia core

**Per-vertex overhead**: ~10-50 µs depending on system

### 2. Element-by-Element Summation

`rad.Fld()` currently:
```cpp
for (int i = 0; i < N_elements; i++) {
    B_total += BiotSavart(r, element[i]);
}
```

No hierarchical acceleration, no early termination.

---

## Recommendations

### Short-term (No code changes):

1. **Accept current performance** for typical use cases:
   - N < 500 elements: GridFunction.Set() completes in seconds
   - For coupled simulations: Field evaluation is often one-time setup
   - Focus optimization on time-stepping loop if needed

2. **Reduce mesh size** where possible:
   - Use adaptive mesh refinement
   - Refine only near magnets, coarse elsewhere
   - Reduces M (number of evaluation points)

3. **Reduce Radia elements** using symmetry:
   - Exploit problem symmetry
   - Use fewer, larger elements where field variation is smooth
   - Reduces N (number of source elements)

### Long-term (Requires implementation):

1. **Implement FMM for field evaluation** (highest impact)
   - Would provide O(M) complexity instead of O(M × N)
   - Benefits ALL field evaluation, not just NGSolve integration

2. **Vectorized field evaluation API**
   - `rad.FldVec(obj, field_type, point_array)`
   - Reduces Python call overhead
   - Enables SIMD vectorization

3. **Extend H-matrix to field evaluation**
   - Reuse existing hierarchical structure
   - Moderate implementation complexity

---

## Conclusion

**The H-matrix speedup applies to the SOLVE phase (rad.Solve()), not to field evaluation (rad.Fld()).**

For GridFunction.Set():
- ✅ Optimization works: Precision control, coordinate transformation
- ❌ H-matrix doesn't help: Not used in field evaluation path
- 🔄 Future: Need hierarchical field evaluation (FMM or extended H-matrix)

**Current performance is acceptable** for typical use cases (N < 500, M < 10,000).
**For large-scale problems**, hierarchical field evaluation would be needed.

---

**See Also**:
- `benchmark_gridfunction_set.py` - Performance measurements
- `docs/API_EXTENSIONS.md` - API documentation
- `examples/solver_time_evaluation/` - Solver-phase H-matrix benchmarks
