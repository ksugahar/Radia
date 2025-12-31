# Matrix Storage Format: Row-Major vs Column-Major

This document analyzes the optimal matrix storage format for each solver method in Radia.

## Summary

| Solver | Optimal Format | Reason |
|--------|---------------|--------|
| **LU (dgesv)** | **Column-Major** | LAPACK is Fortran-based, expects column-major |
| **BiCGSTAB** | **Row-Major** | BLAS dgemv can use either, but row-major is cache-friendly for C++ |
| **HACApK** | **Row-Major** | C implementation, row-major access pattern |

## Current Implementation

Radia currently uses **column-major** format for all solvers:
- `m_flatInteractMatrix` is stored in column-major (Fortran/LAPACK style)
- Element `A(i,j)` is at index `[j * totalDOF + i]`

This is **optimal for LU** but **suboptimal for BiCGSTAB and HACApK**.

## Detailed Analysis

### 1. LU Decomposition (LAPACK dgesv)

**Optimal: Column-Major**

LAPACK's `dgesv` is written in Fortran and expects column-major storage:
```cpp
// Column-major: A(i,j) at [j * n + i]
dgesv_(&n, &nrhs, A_colmajor, &n, ipiv, b, &n, &info);
```

If row-major data is passed, the solver effectively solves A^T * x = b (transposed problem).

**Performance**: No conversion needed when storing column-major.

### 2. BiCGSTAB Iterative Solver

**Optimal: Row-Major** (for C/C++ implementations)

BiCGSTAB's main operation is matrix-vector product `y = A * x`:
```cpp
// Row-major: A(i,j) at [i * n + j]
for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
        sum += A[i * n + j] * x[j];  // Sequential memory access in inner loop
    }
    y[i] = sum;
}
```

**Why row-major is faster**:
- Inner loop accesses `A[i*n + 0], A[i*n + 1], ...` sequentially
- CPU cache prefetching works efficiently
- Each row fits in cache line (spatial locality)

**Column-major access pattern** (current implementation):
```cpp
// Column-major: A(i,j) at [j * n + i]
for (int i = 0; i < n; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
        sum += A[j * n + i] * x[j];  // Strided memory access (stride = n)
    }
    y[i] = sum;
}
```
- Inner loop accesses `A[0*n + i], A[1*n + i], A[2*n + i], ...`
- Strided access causes cache misses
- Performance degradation: 2x-5x slower for large matrices

**Note**: Using BLAS `dgemv` can mitigate this:
```cpp
// Column-major with BLAS
cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1.0, A, n, x, 1, 0.0, y, 1);
```
BLAS internally optimizes for the storage format, but there's still some overhead compared to native row-major.

### 3. HACApK (H-Matrix)

**Optimal: Row-Major**

HACApK is implemented in C (src/ext/HACApK_LH-Cimplm/) and uses row-major storage internally:
```c
// cHACApK_base.c: Row-major access
for (int i = 0; i < ndl; i++) {
    for (int j = 0; j < ndt; j++) {
        waa[i * ndt + j] = ...;  // Row-major storage
    }
}
```

**Current overhead**:
- Radia stores column-major in `m_flatInteractMatrix`
- HACApK requires row-major for its internal H-matrix blocks
- Conversion happens during H-matrix construction (O(n_block^2) per block)

## Recommendations

### Option A: Keep Column-Major (Current)
- Pro: Optimal for LU (most accurate solver)
- Pro: No code changes needed
- Con: Suboptimal for BiCGSTAB and HACApK
- Con: ~50% slower for large iterative problems

### Option B: Use Row-Major for Iterative Solvers
- Store `m_flatInteractMatrix` in row-major
- Transpose to column-major only when calling LU
- Pro: Faster BiCGSTAB and HACApK
- Con: O(n^2) transpose for LU

### Option C: Maintain Both Formats
- Store both row-major and column-major versions
- Pro: Optimal for all solvers
- Con: 2x memory usage for interaction matrix

### Recommendation

**Keep Column-Major** (Option A) for now:
1. LU is used for small problems where performance matters most
2. BiCGSTAB with BLAS dgemv is still efficient
3. HACApK's overhead is amortized over multiple iterations
4. Memory reduction is more important than marginal speed improvement

## Memory Layout Visualization

```
Column-Major (Fortran/LAPACK):     Row-Major (C/C++):
Memory address increases →          Memory address increases →

[A00] [A10] [A20] ...               [A00] [A01] [A02] ...
[A01] [A11] [A21] ...               [A10] [A11] [A12] ...
[A02] [A12] [A22] ...               [A20] [A21] [A22] ...

Element A(i,j) at [j*n + i]         Element A(i,j) at [i*n + j]
```

## References

1. LAPACK Users' Guide: https://www.netlib.org/lapack/lug/
2. Intel MKL Developer Reference: https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/
3. "Anatomy of High-Performance Matrix Multiplication" by Goto & van de Geijn (2008)
