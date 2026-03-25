# Tip: COO extraction for small-scale BEM dense matrices

## Background

I am using `ngsolve.bem` for self-inductance extraction of conductor
surfaces (LaplaceSL, saddle point EFIE). For our problem sizes
(N = 1000-5000 surface DOFs), I need the full dense SL matrix for
LU factorization and the energy inner product `L = mu_0 * J^T @ SL @ J`.

I understand that ngsolve.bem is designed for large-scale problems
where FMM + iterative solvers avoid forming the dense matrix entirely.
However, for small-to-medium BEM problems where the dense matrix IS
needed, I found that `COO()` extraction is much faster than `ToDense()`.

## Observation

With `use_fmm=False`, `LaplaceSL` stores the assembled matrix as
`SparseMatrixdouble` with 100% fill (as expected -- BEM matrices are
dense by nature). Extracting the dense NumPy array:

| Method | Time (N=5085) | Notes |
|--------|---------------|-------|
| Operator creation | 22 s | BEM integral assembly |
| `mat.COO()` + scipy `toarray()` | **0.18 s** | Direct data extraction |
| `mat.ToDense().NumPy()` | **144 s** | N column-by-column MatVecs |

Both methods produce bit-identical results.

## Why the difference

Looking at the NGSolve source (`basematrix.cpp`), `ToDense()` is
implemented in `BaseMatrix` as N column-by-column MatVecs:

```cpp
// basematrix.cpp
template <typename TSCAL>
Matrix<TSCAL> BaseMatrix :: ToDense() const
{
  for (int i = 0; i < fx.Size(); i++)
    {
      fx = 0.0;
      fx(i) = 1;
      Mult (vecx, vecy);
      dmat.Col(i) = fy;
    }
  return dmat;
}
```

This is a reasonable generic implementation -- it works for any
`BaseMatrix` subclass, including FMM operators that have no stored
matrix. For FEM sparse matrices (with << N nonzeros per row), the
overhead is small.

For BEM with `use_fmm=False`, however, the matrix is already fully
assembled in CSR format (100% fill). The N MatVecs each traverse the
full N entries per row through CSR indirect addressing, making the
total cost O(N^2) MatVec operations on an O(N^2) matrix.

`COO()` simply returns the stored triplets, which scipy converts
to dense in a single pass.

## Workaround

```python
from scipy.sparse import coo_matrix

# Fast: COO extraction (~0.18s at N=5085)
rows, cols, vals = mat.COO()
SL = coo_matrix((vals, (rows, cols)),
                shape=(mat.height, mat.width)).toarray()

# Slow: ToDense (~144s at N=5085)
# SL = mat.ToDense().NumPy()
```

## Suggestion

Would it be useful to add a `SparseMatrix::ToDense()` override that
directly copies from the stored CSR data instead of going through
N MatVecs? For small-scale BEM users who need the dense matrix
(circuit extraction, model order reduction, etc.), this would be a
nice quality-of-life improvement.

Alternatively, if the design intent is that dense extraction should
always go through `COO()`, perhaps a note in the documentation would
help users discover this path.

## Environment

- NGSolve 6.2.2602
- Windows Server 2022, Python 3.12, 8-core Xeon

Thank you for the excellent BEM implementation -- the LaplaceSL
operator itself is fast and accurate for our inductance work.
