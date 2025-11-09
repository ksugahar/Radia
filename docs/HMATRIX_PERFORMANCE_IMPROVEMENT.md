# H-Matrix Performance Improvement Plan

## Current Performance Issues

### Benchmark Results (N=343, Nonlinear Material + Background Field)

| Method | Time (s) | Memory | Speedup |
|--------|----------|--------|---------|
| Standard Relaxation | 0.0277 | 8 MB (dense) | 1.00x (baseline) |
| **H-matrix Relaxation** | **1.751** | **11 MB** | **0.02x (50x SLOWER!)** |

**Critical Problems:**
1. **Construction overhead: 99%** (1.74s construction / 1.75s total)
2. **No compression: 148%** (11 MB H-matrix vs 8 MB dense matrix)
3. **Result: 50x SLOWER than standard method**

---

## Root Cause Analysis

### 1. Construction Overhead (99% of runtime)

**Current behavior:**
```
Building 9 H-matrices (3x3 tensor components) in parallel...
Construction time: 1.74186 s
```

**Problem:**
- Each relaxation solve rebuilds ALL 9 H-matrices from scratch
- No reuse between iterations
- Parallel construction helps but overhead still dominates

**Impact:** For N < 500, construction time exceeds computation savings

### 2. No Compression (148% memory usage)

**Current behavior:**
```
Compression ratio: 148.286%
Dense would be: 8 MB
H-matrix memory: 11 MB
```

**Problem:**
- ACA (Adaptive Cross Approximation) creates too many full-rank blocks
- Low-rank approximation not effective
- Admissibility criterion too strict

**Likely causes:**
- `eta` parameter (admissibility) = 0.8 (too conservative)
- `eps` tolerance (ACA) = 1e-6 (too tight)
- Near-field interactions dominate for small problems

### 3. Cache Not Utilized

**Current behavior:**
```
Caching symmetry transformations for 343 elements...
Cached 343 transformation lists
```

**Problem:**
- Symmetry transformations cached but H-matrices rebuilt every solve
- No H-matrix reuse between relaxation iterations
- Magnetization updates don't use incremental updates

---

## Proposed Improvements

### Priority 1: Reduce Construction Overhead

#### 1.1 H-Matrix Reuse Between Iterations
```cpp
// Current: Rebuild every iteration
for (iter = 0; iter < max_iter; iter++) {
    BuildHMatrix();  // ← Expensive! 1.74s
    SolveIteration();
}

// Proposed: Build once, update magnetization
BuildHMatrix();  // Once only
for (iter = 0; iter < max_iter; iter++) {
    UpdateMagnetization();  // Fast update
    SolveIteration();
}
```

**Expected gain:** 50-100x speedup (from 1.75s → 0.02-0.03s)

#### 1.2 Lazy H-Matrix Construction
```cpp
// Only build H-matrix if problem size > threshold
if (n_elements < 100) {
    use_dense_solver();  // Faster for small N
} else {
    use_hmatrix_solver();
}
```

**Expected gain:** Automatic selection of fastest method

### Priority 2: Improve Compression

#### 2.1 Relax ACA Parameters
```cpp
// Current (too strict):
eps = 1e-6;      // ACA tolerance
eta = 0.8;       // Admissibility
max_rank = 50;   // Maximum rank

// Proposed (balanced accuracy/compression):
eps = 1e-4;      // Relax tolerance (still < 1% error)
eta = 1.5;       // More aggressive clustering
max_rank = 30;   // Lower rank for better compression
```

**Expected gain:**
- Compression ratio: 148% → 30-50%
- Memory: 11 MB → 2-4 MB
- Construction: 1.74s → 0.5-1.0s

#### 2.2 Adaptive Parameters by Problem Size
```cpp
if (n_elements < 200) {
    eps = 1e-3;  // Coarse approximation OK
    eta = 2.0;   // Very aggressive
} else {
    eps = 1e-4;  // Standard
    eta = 1.5;
}
```

### Priority 3: Cache Strategy

#### 3.1 Persistent H-Matrix Cache
```cpp
// Cache H-matrix across multiple Solve() calls
static HMatrixCache g_hmatrix_cache;

void Solve(...) {
    if (!g_hmatrix_cache.is_valid(geometry_hash)) {
        g_hmatrix_cache.build(geometry);
    }
    g_hmatrix_cache.solve_with_current_magnetization();
}
```

**Expected gain:** Near-zero overhead for repeated solves

#### 3.2 Incremental Magnetization Updates
Instead of rebuilding H-matrix, update only affected blocks when magnetization changes.

---

## Implementation Roadmap

### Phase 1: Quick Wins (High Impact, Low Risk)
**Target: 10-20x speedup**

1. **H-Matrix Reuse in Relaxation Loop** (Week 1)
   - Build once before iteration loop
   - Update magnetization in-place
   - Expected: 1.75s → 0.1s

2. **Relax ACA Parameters** (Week 1)
   - Change eps: 1e-6 → 1e-4
   - Change eta: 0.8 → 1.5
   - Expected: Better compression, faster construction

3. **Add Automatic Method Selection** (Week 1)
   - Use dense for N < 100
   - Use H-matrix for N ≥ 100
   - Expected: Optimal performance at all scales

### Phase 2: Medium-Term Improvements (Medium Impact)
**Target: Additional 2-5x speedup**

4. **Persistent H-Matrix Cache** (Week 2-3)
   - Cache between Solve() calls
   - Validate with geometry hash
   - Expected: Near-zero overhead for repeated solves

5. **Adaptive Parameters by N** (Week 2)
   - Different parameters for different problem sizes
   - Expected: Better scaling

### Phase 3: Advanced Optimizations (Long-Term)
**Target: Additional 2-3x speedup**

6. **Incremental Magnetization Updates** (Week 4-6)
   - Update affected blocks only
   - Avoid full H-matrix rebuild

7. **GPU Acceleration** (Future)
   - Offload ACA and matrix operations to GPU
   - Expected: 5-10x additional speedup for large N

---

## Expected Performance After Improvements

### Target Performance (N=343)

| Metric | Current | Target | Improvement |
|--------|---------|--------|-------------|
| Construction | 1.74 s | 0.05 s | **35x faster** |
| Total Time | 1.75 s | 0.01 s | **175x faster** |
| Memory | 11 MB | 3 MB | **3.7x less** |
| Compression | 148% | 40% | **3.7x better** |
| Speedup vs Dense | 0.02x | **3-5x** | **150-250x improvement** |

### Validation Criteria

✅ **Success:**
- H-matrix faster than dense for N > 100
- Compression ratio < 50% (better than dense)
- Construction overhead < 20% of total time
- Accuracy maintained (< 1% error)

⚠️ **Acceptable:**
- H-matrix competitive with dense for N > 200
- Compression ratio < 80%

❌ **Failure:**
- H-matrix still slower than dense for any N
- Compression ratio > 100%

---

## Testing Strategy

1. **Benchmark Suite** (`examples/solver_benchmarks/`)
   - Test N = 27, 125, 343, 1000
   - Compare: Dense vs H-matrix
   - Measure: Time, memory, accuracy

2. **Accuracy Verification**
   - Field error < 1% vs dense method
   - Magnetization convergence

3. **Performance Regression Tests**
   - Ensure improvements don't break correctness
   - CI/CD integration

---

## References

- CHANGELOG.md: v1.0.10 H-matrix implementation
- examples/solver_benchmarks/: Current benchmark results
- src/ext/HACApK_LH-Cimplm/: H-matrix library
- src/core/rad_intrc_hmat.cpp: Relaxation solver integration

---

**Status:** 🔴 Critical Performance Issue
**Priority:** P0 (Blocker for production use)
**Owner:** TBD
**Last Updated:** 2025-11-10
