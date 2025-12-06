# H-matrix BiCGSTAB Integration Fix (2025-12-05)

## Summary

Fixed a critical bug where H-matrix acceleration was not being used in the BiCGSTAB solver (Method 1).

## Bug Description

The `DenseMatVec` function in `rad_relaxation_methods.cpp` had a comment claiming H-matrix support but never actually called the H-matrix matvec functions.

**Original code (buggy)**:
```cpp
void radTRelaxationMethNo_1::DenseMatVec(...) {
    // Uses H-matrix if available, otherwise dense matrix  <- LIE!
    if(IntrcMat != nullptr) {
        // Always O(N^2) dense loop - H-matrix never used!
        for(int i = 0; i < n_elem; i++) {
            for(int j = 0; j < n_elem; j++) { ... }
        }
    }
}
```

## Fix Applied

Added proper H-matrix integration to `DenseMatVec` function (lines 566-594):

```cpp
// Check if H-matrix is available and enabled
if(IntrctPtr->IsHMatrixEnabled())
{
    // H-matrix accelerated matrix-vector product: O(N log N)
    IntrctPtr->HMatrixMatVec(x.data(), y.data());

    // Add diagonal term (1/chi) * x
    #pragma omp parallel for if(n_elem > 50)
    for(int i = 0; i < n_elem; i++) {
        // ... diagonal contribution
    }
}
else if(IntrcMat != nullptr)
{
    // Fallback: Dense matrix-vector product O(N^2)
    // ... original code
}
```

## Verification

Benchmark results confirm the fix is **mathematically correct**:

| N | Elements | Dense Bz | H-matrix Bz | Difference |
|---|----------|----------|-------------|------------|
| 8 | 512 | 2.69071912 | 2.69071912 | 0 |
| 10 | 1000 | 2.56547500 | 2.56547500 | 0 |
| 12 | 1728 | 2.50013625 | 2.50013625 | 0 |
| 15 | 3375 | 2.31384211 | 2.31384211 | 0 |
| 18 | 5832 | 2.43170848 | 2.43170848 | 0 |
| 20 | 8000 | 2.42536806 | 2.42536806 | 0 |

Iteration counts are also identical between Dense and H-matrix BiCGSTAB.

## Performance Observation

H-matrix does not show speedup for single compact cube mesh:

| N | Elements | Dense Time | H-matrix Time | Speedup |
|---|----------|------------|---------------|---------|
| 10 | 1000 | 0.55s | 0.55s | 1.01x |
| 15 | 3375 | 7.55s | 7.60s | 0.99x |
| 20 | 8000 | 48.98s | 50.28s | 0.97x |

**Reason**: For a single compact cube mesh, all elements are close together. The H-matrix admissibility criterion `dist(c1, c2) >= eta * min(diam(c1), diam(c2))` is never satisfied because no clusters are "well-separated". All blocks remain dense.

## When H-matrix is Beneficial

H-matrix acceleration is most effective for:

1. **Spatially distributed geometries** (multiple magnets, undulators)
2. **Long aspect ratio objects** (beamline magnets, coils)
3. **Assemblies with air gaps** between components

For single compact objects like a subdivided cube, dense BiCGSTAB may be more efficient due to H-matrix overhead.

## Files Modified

- `src/core/rad_relaxation_methods.cpp` (lines 566-594)

## Related Files

- `src/core/rad_interaction.h` - `IsHMatrixEnabled()`, `HMatrixMatVec()` declarations
- `src/core/rad_hmatrix_aca.cpp` - ACA implementation (unchanged)

## Discovered By

ELF_MAGIC benchmark comparison (2025-12-05)

Bug report: `S:\ELF_MAGIC\01_GitHub\examples\cube_uniform_field\radia_comparison\RADIA_HMATRIX_BUG_REPORT.md`
