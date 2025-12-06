# H-matrix Implementation Review (2025-12-05)

## Summary

Based on ELF_MAGIC's BICGSTAB_COMPARISON.md analysis, this document reviews the differences between Radia's custom ACA implementation and HACApK library, with recommendations for future development.

## Current Implementation Comparison

### Radia: Custom ACA (`rad_hmatrix_aca.cpp`)

**Location**: `src/core/rad_hmatrix_aca.cpp/h`

**Features**:
- Custom cluster tree construction (bisection by longest axis)
- Standard admissibility criterion: `dist >= eta * min(diam1, diam2)`
- ACA+ low-rank approximation
- OpenMP parallelization for block computation
- Parameters from HACApK: eta=2.0, min_cluster_size=15

**Code structure** (lines 329-356):
```cpp
bool radTHMatrixACA::IsAdmissible(const radTCluster* c1, const radTCluster* c2) const
{
    // Parameters from ELF_MAGIC m_HACApK_base.f90
    const int min_cluster_size = 15;  // param(21)
    if(c1->size < min_cluster_size || c2->size < min_cluster_size) return false;

    double dist = c1->bbox.distance(c2->bbox);
    double min_diam = std::min(c1->bbox.diameter(), c2->bbox.diameter());

    return dist >= m_eta * min_diam;  // eta=2.0 from param(51)
}
```

### ELF_MAGIC: HACApK Library

**Location**: `S:\ELF_MAGIC\01_GitHub\src\ppOpenHPC-MATH-HACApK\src\HACApK_1.0.0\`

**Features**:
- Proven H-matrix library from ppOpen-HPC project
- Matrix-free ACA+ (computes matrix elements on-demand)
- Integrated BiCGSTAB solver
- MPI parallelization support (for supercomputer use)
- BLAS Level 3 optimizations for MatVec

**Key function**: `cHACApK_acaplus()` in `cHACApK_base.c`:700

## Key Differences

| Aspect | Radia Custom ACA | HACApK Library |
|--------|------------------|----------------|
| Matrix storage | Pre-computed full matrix | Matrix-free (on-demand) |
| BLAS usage | Custom loops + OpenMP | Intel MKL (dgemv, dgemm) |
| BiCGSTAB | Separate implementation | Integrated solver |
| Parallelization | OpenMP (shared memory) | MPI + OpenMP |
| Cluster tree | Simple bisection | Octree-based |
| Memory efficiency | Dense matrix backup | True matrix-free |

## Bug Fix Applied (2025-12-05)

**Problem**: `DenseMatVec` function never called H-matrix routines despite comment claiming H-matrix support.

**Fix**: Added `IsHMatrixEnabled()` check and `HMatrixMatVec()` call in `rad_relaxation_methods.cpp`:566-594.

**Result**: Mathematically correct (identical Bz values), but no speedup for compact cube geometry.

## Why No Speedup for Compact Cubes

For a single compact cube mesh:
1. All elements are adjacent/close
2. No clusters satisfy `dist >= eta * min(diam1, diam2)`
3. All blocks become dense (no low-rank compression)
4. H-matrix overhead (cluster tree, block bookkeeping) adds cost

**Benchmark results** (linear material, BiCGSTAB Method 1):

| N | Elements | Dense Time | H-matrix Time | Speedup |
|---|----------|------------|---------------|---------|
| 10 | 1000 | 0.55s | 0.55s | 1.01x |
| 15 | 3375 | 7.55s | 7.60s | 0.99x |
| 20 | 8000 | 48.98s | 50.28s | 0.97x |

## When H-matrix IS Beneficial

From BICGSTAB_COMPARISON.md crossover analysis:

- **N < 12**: Dense LU is fastest
- **12 <= N < 25**: Dense BiCGSTAB most efficient
- **25 <= N < 50**: Problem-dependent (geometry matters)
- **N >= 50**: HACApK BiCGSTAB should win

**Effective geometries**:
1. Spatially distributed magnets (undulators, magnet arrays)
2. Long aspect ratio objects (beamline magnets)
3. Assemblies with air gaps

## Recommendations

### Short-term (Current state)

1. Document that H-matrix is NOT beneficial for compact single objects
2. Add `SolverHMatrixStats()` API to report compression ratio
3. Recommend Dense BiCGSTAB for N < 5000 compact objects

### Medium-term (Migration to HACApK)

1. Replace custom ACA with HACApK C library calls
2. Use HACApK's `cHACApK_acaplus()` for ACA+ compression
3. Use HACApK's integrated BiCGSTAB solver
4. Keep Radia's interaction matrix as fallback

### Long-term (Full HACApK Integration)

1. Adopt matrix-free approach (compute N_ij on-demand)
2. Use HACApK's cluster tree and admissibility
3. Support both shared-memory and distributed computing

## Files to Modify

| File | Current State | Required Change |
|------|---------------|-----------------|
| `rad_hmatrix_aca.cpp` | Custom ACA | Wrap HACApK calls |
| `rad_interaction.cpp` | H-matrix build | Use HACApK build |
| `rad_relaxation_methods.cpp` | BiCGSTAB | Option for HACApK solver |
| `CMakeLists.txt` | No HACApK link | Link HACApK library |

## Reference

- ELF_MAGIC BICGSTAB_COMPARISON.md: `S:\ELF_MAGIC\01_GitHub\examples\cube_uniform_field\radia_comparison\BICGSTAB_COMPARISON.md`
- HACApK library: `src/ext/HACApK_LH-Cimplm/`
- Bug fix document: `internal/design/HMATRIX_BICGSTAB_FIX_2025.md`

---

**Author**: Claude Code
**Date**: 2025-12-05
**Status**: Review Complete
