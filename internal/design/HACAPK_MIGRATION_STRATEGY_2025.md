# HACApK Migration Strategy (2025-12-05)

## Executive Summary

After detailed analysis of both Radia's current H-matrix implementation and the HACApK library,
a direct migration to HACApK is **NOT recommended** due to:

1. **MPI Dependency** - HACApK requires MPI even for single-process use
2. **Complexity Overhead** - Significant wrapper code needed for C++ integration
3. **Limited Benefit for Typical Radia Use Cases** - Single compact objects show no speedup

Instead, we recommend **improving the current implementation** with lessons learned from HACApK.

## Analysis Results

### HACApK Library Characteristics

**Location**: `src/ext/HACApK/` and `src/ext/HACApK_LH-Cimplm/`

**Requirements**:
- MPI (required, line 43: `#include <mpi.h>`)
- OpenMP (optional but highly recommended)
- BLAS/LAPACK (for DGEMV, DGEMM operations)

**Key Functions**:
- `cHACApK_acaplus()` - ACA+ low-rank approximation
- `cHACApK_generate_cbitree()` - Cluster tree construction
- `cHACApK_generate_frame_blrleaf()` - H-matrix frame generation

**Matrix Element Computation**:
- Uses callback `cHACApK_calc_entry_ij()` for on-demand computation
- Enables true matrix-free operation (no O(N^2) pre-computation)

### ELF_MAGIC's Approach

**Key Finding**: ELF_MAGIC does NOT directly use HACApK for its DLL interface.

Instead, `m_dll_hmatrix.f90` implements:
- Custom H-matrix structure (`st_hmatrix`)
- Block-based storage with ACA compression
- `hmatrix_build_from_dense()` - builds from pre-computed matrix
- **Important**: ACA is temporarily disabled (line 164-168) due to issues

### Radia's Current Implementation

**Location**: `src/core/rad_hmatrix_aca.cpp`

**Characteristics**:
- Custom ACA+ implementation
- OpenMP parallelization
- Cluster tree with bounding boxes
- Standard admissibility: `dist >= eta * min(diam1, diam2)`

**Problem Identified**:
The `ComputeEntry()` function reads from pre-computed dense matrix:
```cpp
void radTHMatrixACA::ComputeEntry(int i, int j, double* entry_3x3) const
{
    double* mat = m_interaction->IntrcMat;  // Requires O(N^2) pre-computation!
    int base_idx = (i * n + j) * 9;
    for(int k = 0; k < 9; ++k) {
        entry_3x3[k] = mat[base_idx + k];
    }
}
```

This defeats the memory benefit of H-matrix (matrix-free approach).

## Why HACApK Migration is NOT Recommended

### 1. MPI Requirement

HACApK is designed for distributed computing. Even single-process use requires MPI initialization:
```c
MPI_Comm_rank(comm, &irank);
MPI_Comm_split(comm, iclr, ikey, &commn);
```

Radia is a single-process library. Adding MPI dependency would:
- Require MPI library installation for all users
- Add complexity for basic use cases
- Not provide benefit for typical workstations

### 2. No Benefit for Compact Objects

Benchmark results (2025-12-05) show H-matrix provides **no speedup** for single compact cubes:

| N | Elements | Dense Time | H-matrix Time | Speedup |
|---|----------|------------|---------------|---------|
| 10 | 1000 | 0.55s | 0.55s | 1.01x |
| 15 | 3375 | 7.55s | 7.60s | 0.99x |
| 20 | 8000 | 48.98s | 50.28s | 0.97x |

**Reason**: All elements are close together, no clusters satisfy admissibility criterion.

### 3. Integration Complexity

Wrapping HACApK for C++ use requires:
- MPI stub or initialization code
- Fortran-C interoperability layer (if using Fortran HACApK)
- Callback registration for `calc_entry_ij`
- Memory management bridging

## Recommended Approach

### Phase 1: Optimize Current Implementation (Short-term)

1. **Fix the core issue** - Implement matrix-free `ComputeEntry()`:
   ```cpp
   void radTHMatrixACA::ComputeEntry(int i, int j, double* entry_3x3) const
   {
       // Compute interaction matrix element on-demand
       // instead of reading from pre-computed matrix
       m_interaction->ComputeInteractionBlock(i, j, entry_3x3);
   }
   ```

2. **Improve admissibility criterion** - Use HACApK's proven parameters:
   - `eta = 2.0` (current)
   - `min_cluster_size = 15` (from HACApK param[21])

3. **Add BLAS Level 3 operations** - Use OpenBLAS for MatVec:
   - Dense blocks: `dgemv()`
   - Low-rank blocks: `dgemm()` for U*V^T*x

### Phase 2: Selective HACApK Algorithms (Medium-term)

Extract and adapt specific HACApK algorithms (without MPI dependency):

1. **ACA+ Algorithm** - Port `cHACApK_acaplus()` to C++
2. **Cluster Tree** - Use HACApK's proven cluster construction
3. **Parameter Tuning** - Copy HACApK's optimal parameters

### Phase 3: Full Matrix-Free Implementation (Long-term)

For problems where O(N^2) matrix storage is prohibitive:

1. Implement on-demand interaction computation
2. Use HACApK-style H-matrix structure
3. Consider optional MPI support for cluster computing

## Immediate Actions

1. **Document current H-matrix limitations**
   - Add warning in user documentation
   - Recommend Dense BiCGSTAB for N < 5000 compact objects

2. **Update CLAUDE.md**
   - H-matrix policy: not beneficial for compact single objects
   - Recommend solver selection based on geometry

3. **Add `SolverHMatrixStats()` API**
   - Return compression ratio, number of low-rank blocks
   - Help users understand H-matrix effectiveness

## Files to Track

| File | Purpose | Status |
|------|---------|--------|
| `rad_hmatrix_aca.cpp` | Current H-matrix implementation | Keep, optimize |
| `src/ext/HACApK/` | HACApK C library | Reference only |
| `internal/design/HMATRIX_IMPLEMENTATION_REVIEW_2025.md` | Implementation review | Complete |
| `internal/design/HMATRIX_BICGSTAB_FIX_2025.md` | BiCGSTAB fix documentation | Complete |

## Conclusion

**Do NOT migrate to HACApK library directly.**

Instead:
1. Keep current `rad_hmatrix_aca.cpp` implementation
2. Fix matrix-free computation
3. Extract useful algorithms from HACApK (ACA+, clustering)
4. Document that H-matrix is beneficial only for spatially distributed geometries

---

**Author**: Claude Code
**Date**: 2025-12-05
**Status**: Strategy Defined
