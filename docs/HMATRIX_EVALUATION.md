# H-Matrix Acceleration Evaluation for Radia

**Date**: 2025-12-05
**Status**: Evaluation Complete - NOT Recommended for Radia

## Executive Summary

H-matrix (Hierarchical Matrix) acceleration was evaluated for Radia's MMM (Magnetic Moment Method) solver. **Conclusion: H-matrix provides NO benefit for Radia's typical use cases** and is NOT recommended for integration.

## Benchmark Results

### Single Compact Object (Cube in Uniform Field)

| N_div | Elements | DOF | LU Time | BiCGSTAB Time | H-matrix BiCGSTAB |
|-------|----------|-----|---------|---------------|-------------------|
| 10 | 1,000 | 3,000 | 2.01s | 0.54s | 0.53s (**1.02x**) |
| 15 | 3,375 | 10,125 | 24.28s | 7.56s | 7.27s (**1.04x**) |
| 20 | 8,000 | 24,000 | 278.68s | 47.05s | 47.33s (**0.99x**) |

**Key Finding**: H-matrix shows **NO speedup** (0.97x-1.04x ratio) for compact objects.

## Why H-Matrix Does NOT Help Radia

### 1. Admissibility Criterion Failure

H-matrix partitions the interaction matrix into:
- **Far-field blocks**: Low-rank approximation (ACA)
- **Near-field blocks**: Dense (full) storage

The admissibility criterion is:
```
dist(cluster_i, cluster_j) >= eta * min(diam(cluster_i), diam(cluster_j))
```

For a **single compact object** (typical Radia use case):
- All elements are spatially close together
- `dist` is always small relative to `diam`
- **NO blocks satisfy the admissibility criterion**
- **ALL blocks remain dense** -> No compression benefit

### 2. Radia's Typical Use Cases

Radia is designed for:
- Single magnet assemblies (undulators, wigglers)
- Compact electromagnet designs
- Objects where elements are tightly packed

These geometries have **high element density** with **no spatial separation**, making H-matrix ineffective.

### 3. When H-Matrix WOULD Help

H-matrix is beneficial for:
- **Multiple separated objects** (e.g., array of magnets with large gaps)
- **Particle simulations** with distributed sources
- **BEM** for 3D geometries with large surface-to-volume ratio

These are **NOT typical Radia use cases**.

## HACApK Integration Assessment

### HACApK Library Overview
- **Source**: ppOpen-HPC project (MIT License)
- **Purpose**: H-matrix with ACA+ for large-scale problems
- **Design**: MPI-parallel, targeting supercomputers

### Integration Challenges

1. **MPI Dependency**
   - HACApK requires MPI (`#include <mpi.h>`)
   - Solution: MPI stub header created (`mpi_stub.h`)
   - Status: Solved

2. **Callback Interface**
   - HACApK requires `cHACApK_calc_entry_ij.h` for matrix element computation
   - Must provide Radia-specific implementation
   - Complexity: Medium

3. **Code Adaptation**
   - ~5000+ lines of C/Fortran code
   - Designed for distributed memory systems
   - Would require significant adaptation for single-process use

4. **Build System Integration**
   - CMake changes needed
   - Cross-platform testing required

### Effort vs. Benefit Analysis

| Aspect | Effort | Benefit |
|--------|--------|---------|
| MPI stub | Low | - |
| Callback implementation | Medium | - |
| Code adaptation | High | - |
| Build integration | Medium | - |
| Testing/debugging | High | - |
| **Total Effort** | **HIGH** | **ZERO** (for typical use cases) |

## Recommendation

**DO NOT integrate HACApK** for the following reasons:

1. **No benefit for typical use cases**: H-matrix provides zero speedup for compact objects
2. **High integration effort**: ~2-4 weeks of development time
3. **Code complexity**: Adds 5000+ lines of code to maintain
4. **Limited applicability**: Only helps for separated objects (rare in Radia)

### Alternative Optimizations

For Radia performance improvement, consider:

1. **OpenMP parallelization** (already implemented)
   - Matrix-vector product parallelization
   - LU decomposition parallelization

2. **BLAS/LAPACK optimization**
   - Intel MKL with optimized kernels (required)

3. **Preconditioner improvements**
   - Better Jacobi/block-Jacobi preconditioners
   - ILU preconditioner for BiCGSTAB

4. **Memory layout optimization**
   - Cache-friendly data structures
   - SIMD vectorization

## Files Created During Evaluation

- `src/ext/HACApK_LH-Cimplm/mpi_stub.h` - MPI stub for single-process

## References

1. HACApK: https://github.com/RIKENGITHUB/ppOpen-HPC
2. Original paper: Ida, A., et al. "HACApK: An H-matrix library for Krylov methods"
3. Radia benchmarks: `examples/cube_uniform_field/nonlinear/`

---

**Conclusion**: H-matrix acceleration is NOT suitable for Radia. Focus on other optimization strategies.
