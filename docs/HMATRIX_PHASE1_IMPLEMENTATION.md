# H-Matrix Phase 1 Implementation Plan

## Overview

Phase 1の目標：**10-20倍の高速化**を達成

## Task List

### Task 1: ACAパラメータ緩和（最優先・最も簡単）
**Expected gain**: 構築時間 1.74s → 0.5-1.0s

#### 1.1 デフォルトパラメータの変更
**File**: `src/core/rad_intrc_hmat.h` or where config is defined

**Current values:**
```cpp
eps = 1e-6;      // ACA tolerance (too strict)
eta = 0.8;       // Admissibility parameter (too conservative)
max_rank = 50;   // Maximum rank
```

**Proposed values:**
```cpp
eps = 1e-4;      // Relax tolerance (still < 1% error)
eta = 1.5;       // More aggressive clustering
max_rank = 30;   // Lower rank for better compression
```

**Implementation:**
1. Find `radTHMatrixSolverConfig` struct definition
2. Change default values
3. Test with benchmark_solver_methods.py
4. Verify accuracy < 1%

**Test command:**
```bash
python examples/solver_benchmarks/benchmark_solver_methods.py
```

**Success criteria:**
- Compression ratio: 148% → < 80%
- Construction time: 1.74s → < 1.0s
- Accuracy: < 1% error

---

### Task 2: 自動手法選択（簡単・大きな効果）
**Expected gain**: Small problems (N<100) use dense solver

#### 2.1 N < 100で密行列ソルバーを使用
**File**: `src/core/rad_interaction.cpp` (CompRelaxInt method)

**Current behavior:**
- Always builds H-matrix if enabled

**Proposed logic:**
```cpp
if (g_use_hmatrix_relaxation && n_elements >= 100) {
    // Use H-matrix
    if (hmat_intrct == nullptr) {
        hmat_intrct = new radTHMatrixInteraction(this, hmatrix_config);
    }
    hmat_intrct->CompRelaxInt(...);
} else {
    // Use dense solver
    CompRelaxInt_Dense(...);  // Standard method
}
```

**Implementation steps:**
1. Find where g_use_hmatrix_relaxation is checked
2. Add `n_elements >= 100` condition
3. Ensure dense solver fallback works
4. Test with N=27, 125, 343

**Success criteria:**
- N=27: Uses dense (fast)
- N=125, 343: Uses H-matrix (if beneficial)
- No crashes or errors

---

### Task 3: H-Matrix再利用（最重要・中難度）
**Expected gain**: 1.75s → 0.01-0.03s (50-100x speedup!)

#### 3.1 現状の問題
**Current code location**: Likely in `rad_interaction.cpp` or `rad_relaxation_methods.cpp`

**Current behavior** (pseudocode):
```cpp
AutoRelax(precision, max_iter) {
    for (iter = 0; iter < max_iter; iter++) {
        // Problem: H-matrix rebuilt EVERY iteration!
        if (use_hmatrix) {
            Build H-matrix;  // 1.74s ← EXPENSIVE!
        }

        Compute new magnetization;

        if (converged) break;
    }
}
```

**Root cause**: H-matrix depends on geometry (element positions), NOT magnetization.
Once built, it can be reused for all iterations.

#### 3.2 提案する実装

**Option A: Build once before loop** (Recommended - Simple)
```cpp
AutoRelax(precision, max_iter) {
    // Build H-matrix ONCE before iteration loop
    if (use_hmatrix) {
        if (hmat_intrct == nullptr || !hmat_intrct->IsBuilt()) {
            hmat_intrct->BuildHMatrix();  // Only once!
        }
    }

    // Iteration loop
    for (iter = 0; iter < max_iter; iter++) {
        // Use existing H-matrix (no rebuild!)
        Compute new magnetization using hmat_intrct;

        if (converged) break;
    }
}
```

**Option B: Lazy build with validation** (More robust)
```cpp
AutoRelax(precision, max_iter) {
    for (iter = 0; iter < max_iter; iter++) {
        if (use_hmatrix) {
            // Only build if not built OR geometry changed
            if (!hmat_intrct->IsBuilt() || geometry_changed) {
                hmat_intrct->BuildHMatrix();
            }
            // else: Reuse existing H-matrix
        }

        Compute new magnetization;

        if (converged) break;
    }
}
```

#### 3.3 実装ステップ

**Step 1**: Find relaxation loop location
```bash
grep -rn "MaxIterNumber" src/core/
grep -rn "DefineNewMagnetizations" src/core/
```

**Step 2**: Add `IsBuilt()` check to radTHMatrixInteraction
```cpp
class radTHMatrixInteraction {
    bool is_built;  // Already exists!

public:
    bool IsBuilt() const { return is_built; }
    void Invalidate() { is_built = false; }  // For geometry changes
};
```

**Step 3**: Move BuildHMatrix call outside loop
- Find current location (inside loop?)
- Move to before loop
- Add is_built check

**Step 4**: Test thoroughly
```bash
python examples/solver_benchmarks/benchmark_solver_methods.py
python tests/test_relaxation.py
```

**Success criteria:**
- BuildHMatrix called ONCE per Solve(), not once per iteration
- Total time: 1.75s → ~0.05s (construction) + ~0.01s (iterations)
- Accuracy maintained
- No memory leaks

---

## Implementation Order

### Week 1 (Quick Wins)
1. **Day 1-2**: Task 1 (ACA parameters) - 簡単
2. **Day 2-3**: Task 2 (Auto selection) - 簡単
3. **Day 4-5**: Task 3 (H-matrix reuse) - 中難度だが最重要

**Total expected:** 10-50x speedup

---

## Code Locations to Investigate

### Primary files:
1. `src/core/rad_relaxation_methods.cpp`
   - `radTRelaxationMethNo_8::AutoRelax()`
   - Iteration loop location

2. `src/core/rad_interaction.cpp`
   - `radTInteraction::CompRelaxInt()`
   - H-matrix usage

3. `src/core/rad_intrc_hmat.cpp`
   - `radTHMatrixInteraction::BuildHMatrix()`
   - Construction logic

4. `src/core/rad_intrc_hmat.h`
   - `radTHMatrixSolverConfig`
   - Default parameters

### Configuration:
```cpp
struct radTHMatrixSolverConfig {
    double eps;           // ACA tolerance
    double eta;           // Admissibility
    int max_rank;         // Maximum rank
    int min_cluster_size; // Minimum cluster size
};
```

---

## Testing Strategy

### Unit tests:
```bash
# Accuracy test
python tests/test_relaxation.py

# Performance benchmarks
python examples/solver_benchmarks/benchmark_solver_methods.py
python examples/solver_benchmarks/benchmark_hmatrix_field.py
```

### Validation:
1. **Accuracy**: Field error < 1% vs dense
2. **Performance**:
   - N=27: Dense faster (auto-selection)
   - N=125: H-matrix 2-5x faster
   - N=343: H-matrix 5-10x faster
3. **Memory**: Compression ratio < 80%

---

## Rollback Plan

If any task causes problems:

1. **Git branches**: Create feature branches
   ```bash
   git checkout -b feature/hmatrix-aca-params
   git checkout -b feature/hmatrix-auto-select
   git checkout -b feature/hmatrix-reuse
   ```

2. **Incremental commits**: One task per commit

3. **Quick rollback**: `git revert` if needed

---

## Next Steps

**Immediate actions:**
1. ✅ Create implementation plan (this document)
2. ⬜ Investigate code locations (grep commands above)
3. ⬜ Start with Task 1 (easiest, good warm-up)
4. ⬜ Progress to Task 2 and 3

**Status**: 🟡 Ready to start implementation
**Priority**: P0
**Owner**: TBD
**Created**: 2025-11-10
