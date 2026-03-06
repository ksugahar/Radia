# Radia Nonlinear Iteration Unification Plan

## Overview

This document describes a **safe, incremental migration plan** to unify the nonlinear iteration code across all 3 solver methods in Radia, following the ELF architecture pattern.

**Critical Constraint**: A previous migration attempt "destroyed" Radia. This plan is designed with multiple safety checkpoints and rollback capability.

## Current State Analysis

### Code Duplication (~1800 lines)

| Solver | Method | Lines | Size |
|--------|--------|-------|------|
| LU | 0 | 489-1136 | ~650 lines |
| BiCGSTAB | 1 | 1501-2100+ | ~600 lines |
| HACApK | 2 | 2317-2858 | ~540 lines |

### Identical Code Blocks (Duplicated in All 3)

1. **Initialization Phase** (~100 lines each):
   - Get totalDOF, FlatMagn, FlatField, FlatExtern arrays
   - Initialize chi using `GetInitialChi_ELF_Style()`
   - Detect linear materials (`all_materials_linear` flag)
   - Cache polyhedron pointers (`polyCache`)
   - Initialize H field and NewFieldArray

2. **Pre-Solve Setup** (~50 lines each):
   - Store old values (OldMagn, OldBnorm)
   - Build/update system matrix diagonal with 1/chi

3. **Post-Solve Update** (~150 lines each):
   - Update element magnetization from FlatMagn
   - Handle 3DOF (tetrahedra) and 6DOF (hexahedra MSC) elements
   - Chi update using `ComputeChiDualMethod()`
   - B-field convergence check (ELF mucal2 style)

4. **Convergence Check** (~30 lines each):
   - Linear material shortcut (converge in 1 iteration)
   - Relative change check: `if(rel_change <= PrecOnMagnetiz) break`

### Unique Code Per Solver

| Solver | Unique Code |
|--------|-------------|
| LU | Matrix copy, `dgesv()` call (~20 lines) |
| BiCGSTAB | `SolveBiCGSTAB_VariableDOF()` call (~10 lines) |
| HACApK | H-matrix build, `SolveBiCGSTAB_HMatrix_VariableDOF()` call (~30 lines) |

## Migration Strategy: "Wrapper First, Replace Later"

### Philosophy

1. **DO NOT modify existing solver code initially**
2. Create a NEW unified function that CALLS existing solvers
3. Verify the wrapper produces identical results
4. Only then refactor internals

### Why This Is Safe

- Existing code remains unchanged until verified
- Easy rollback: just use old functions directly
- Each step can be tested independently
- No "big bang" refactoring

## Implementation Phases

### Phase 0: Baseline Verification (BEFORE ANY CHANGES)

**Goal**: Establish reference results that MUST be matched after each change.

**Actions**:
1. Run full benchmark suite (nonlinear hex N=5,10,15, tetra maxh=0.35,0.25,0.20)
2. Record exact M_avg_z values for all 3 solvers
3. Save as `baseline_results_YYYYMMDD.json`

**Checkpoint**: Results file exists, all tests pass

---

### Phase 1: Add Unified Entry Point (NO FUNCTIONAL CHANGE)

**Goal**: Create a dispatcher function that routes to existing implementations.

**File**: `rad_relaxation_methods.cpp`

**New Function**:
```cpp
// Add to radTRelaxationMethNo_0 (or base class)
int AutoRelax_Unified(double PrecOnMagnetiz, int MaxIterNumber,
                      char MagnResetIsNotNeeded, int solverMethod)
{
    switch(solverMethod) {
        case 0: return AutoRelax_VariableDOF(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
        case 1: return static_cast<radTRelaxationMethNo_1*>(this)->AutoRelax_VariableDOF(...);
        case 2: return static_cast<radTRelaxationMethNo_2*>(this)->AutoRelax_VariableDOF(...);
        default: return AutoRelax_VariableDOF(PrecOnMagnetiz, MaxIterNumber, MagnResetIsNotNeeded);
    }
}
```

**Checkpoint**:
- Compile succeeds
- Calling `AutoRelax_Unified(tol, max, false, 0)` gives EXACT same results as `AutoRelax_VariableDOF()`
- Benchmark results match baseline

---

### Phase 2: Extract Common Initialization

**Goal**: Move initialization code to shared helper functions.

**New Helper Functions**:
```cpp
// In rad_relaxation_methods.cpp (static or member functions)

struct NonlinearContext {
    int totalDOF;
    std::vector<double> FlatMagn;
    std::vector<double> FlatField;
    std::vector<double> FlatExtern;
    std::vector<double> CurrentChiArray;
    std::vector<double> OldMagn;
    std::vector<double> NewFieldArray;
    std::vector<radTg3dRelax*> polyCache;
    bool all_materials_linear;
    // ... other shared state
};

bool InitializeNonlinearContext(NonlinearContext& ctx);
void StoreOldValues(NonlinearContext& ctx);
void UpdateMagnetization(NonlinearContext& ctx);
void UpdateChi(NonlinearContext& ctx, double relax_param);
bool CheckConvergence(NonlinearContext& ctx, double tol, double& rel_change);
```

**Migration Steps**:
1. Create `NonlinearContext` struct
2. Add `InitializeNonlinearContext()` - copy code from LU solver
3. Modify LU solver to use helper (test)
4. Modify BiCGSTAB solver to use same helper (test)
5. Modify HACApK solver to use same helper (test)

**Checkpoint**: After EACH step, run benchmark and verify EXACT match with baseline

---

### Phase 3: Extract Common Post-Solve Update

**Goal**: Unify the magnetization update and chi update logic.

**New Helper Functions**:
```cpp
void UpdateElementMagnetization(NonlinearContext& ctx);  // 3DOF/6DOF handling
void ComputeNewChi(NonlinearContext& ctx);               // Uses ComputeChiDualMethod
double ComputeConvergenceMetric(NonlinearContext& ctx);  // B-field based
```

**Migration Steps**:
1. Add `UpdateElementMagnetization()` - extract from LU solver
2. Test LU solver with extracted function
3. Apply to BiCGSTAB (test)
4. Apply to HACApK (test)
5. Repeat for `ComputeNewChi()` and `ComputeConvergenceMetric()`

**Checkpoint**: After EACH step, benchmark must match baseline

---

### Phase 4: Create Linear Solve Dispatcher

**Goal**: Abstract the linear solve step.

**New Interface**:
```cpp
enum class LinearSolverType { LU = 0, BiCGSTAB = 1, HACApK = 2 };

// Abstract linear solve - returns number of iterations (0 for direct solvers)
int SolveLinearStep(NonlinearContext& ctx, LinearSolverType solver,
                    double* SystemMatrix,  // For LU
                    double bicg_tol,       // For iterative
                    int bicg_max_iter);
```

**Implementation**:
```cpp
int SolveLinearStep(NonlinearContext& ctx, LinearSolverType solver, ...) {
    switch(solver) {
        case LinearSolverType::LU:
            // Call dgesv
            return 0;
        case LinearSolverType::BiCGSTAB:
            return SolveBiCGSTAB_VariableDOF(...);
        case LinearSolverType::HACApK:
            return SolveBiCGSTAB_HMatrix_VariableDOF(...);
    }
}
```

**Checkpoint**: All 3 solvers work through dispatcher, results match baseline

---

### Phase 5: Create Unified Nonlinear Loop

**Goal**: Single implementation of Newton-Raphson iteration.

**New Function**:
```cpp
int AutoRelax_Unified_Impl(double PrecOnMagnetiz, int MaxIterNumber,
                           char MagnResetIsNotNeeded, LinearSolverType solver)
{
    NonlinearContext ctx;
    if(!InitializeNonlinearContext(ctx)) return 0;

    // Pre-build H-matrix if HACApK
    if(solver == LinearSolverType::HACApK) {
        BuildHMatrix(ctx);
    }

    int iterCount = 0;
    for(; iterCount < MaxIterNumber; iterCount++) {
        StoreOldValues(ctx);

        // The ONLY solver-specific part
        int lin_iter = SolveLinearStep(ctx, solver, ...);
        m_solve_bicg_iter += lin_iter;

        UpdateElementMagnetization(ctx);
        ComputeNewChi(ctx);

        double rel_change = ComputeConvergenceMetric(ctx);
        if(rel_change <= PrecOnMagnetiz) break;

        // Linear material shortcut
        if(ctx.all_materials_linear) break;
    }

    return iterCount + 1;
}
```

**Checkpoint**: Unified function produces IDENTICAL results for all 3 solvers

---

### Phase 6: Cleanup (OPTIONAL)

**Goal**: Remove duplicated code from original solvers.

**Actions**:
1. Mark old `AutoRelax_VariableDOF()` methods as deprecated
2. Have them call `AutoRelax_Unified_Impl()` internally
3. Eventually remove deprecated methods

**Note**: This phase is optional and can be deferred indefinitely. The unified implementation is the primary goal.

---

## Test Checkpoints

### After Each Phase

1. **Compile Test**: Code compiles without errors/warnings
2. **Unit Test**: Individual solver functions work correctly
3. **Benchmark Test**: Results match baseline EXACTLY (M_avg_z within 1e-6)
4. **Regression Test**: All existing examples still work

### Baseline Test Cases

| Test Case | Solver | Expected M_avg_z |
|-----------|--------|------------------|
| Hex N=10 Nonlinear | LU | 716281 |
| Hex N=10 Nonlinear | BiCGSTAB | 716307 |
| Hex N=10 Nonlinear | HACApK | 716362 |
| Linear Hex N=10 | All | 716110 (1 iteration) |

### Rollback Procedure

If any checkpoint fails:
1. `git stash` current changes
2. Verify baseline still works
3. Analyze what broke
4. Start phase again with corrected approach

---

## Risk Mitigation

### Known Risks

1. **Pointer aliasing**: Different solvers may have subtle differences in how they access shared data
   - Mitigation: Use `NonlinearContext` struct to encapsulate all state

2. **Initialization order**: Some arrays may need specific initialization order
   - Mitigation: Preserve exact order from working LU solver

3. **Memory management**: Different solvers may have different allocation patterns
   - Mitigation: Keep memory allocation in solver-specific code initially

4. **Thread safety**: HACApK uses OpenMP internally
   - Mitigation: Keep solver-specific threading in `SolveLinearStep()`

### Safety Rules

1. **NEVER** delete working code - only add new code initially
2. **ALWAYS** run benchmarks after each change
3. **KEEP** old functions callable until new is verified
4. **COMMIT** after each successful checkpoint
5. **DOCUMENT** any deviations from plan

---

## Timeline Estimate

| Phase | Effort | Risk |
|-------|--------|------|
| Phase 0 (Baseline) | 30 min | Low |
| Phase 1 (Dispatcher) | 1 hour | Low |
| Phase 2 (Init Extract) | 2 hours | Medium |
| Phase 3 (Update Extract) | 2 hours | Medium |
| Phase 4 (Linear Solve) | 1 hour | Low |
| Phase 5 (Unified Loop) | 2 hours | Medium |
| Phase 6 (Cleanup) | Optional | Low |

**Total**: ~8-10 hours of careful, incremental work

---

## Success Criteria

1. Single `AutoRelax_Unified_Impl()` function handles all 3 solvers
2. Benchmark results match baseline for all test cases
3. Code reduction of ~1200+ lines (keeping ~600 for unified implementation)
4. No regression in any existing functionality
5. Easier to add new solver methods in future (just add case to `SolveLinearStep`)

---

## Appendix: File Changes Summary

### Modified Files

- `src/core/rad_relaxation_methods.cpp` - Main implementation
- `src/core/rad_relaxation_methods.h` - New function declarations

### New Code (~300-400 lines)

- `NonlinearContext` struct
- Helper functions (5-6 functions)
- `SolveLinearStep()` dispatcher
- `AutoRelax_Unified_Impl()` main function

### Removed Code (Phase 6, ~1200 lines)

- Duplicated initialization in BiCGSTAB/HACApK
- Duplicated update logic in BiCGSTAB/HACApK
- Duplicated convergence check in BiCGSTAB/HACApK

---

**Document Version**: 1.0
**Date**: 2025-12-28
**Author**: Claude Code (AI Assistant)
**Status**: DRAFT - Awaiting User Approval
