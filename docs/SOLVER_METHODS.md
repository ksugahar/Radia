# Radia Solver Methods

**Version:** 1.3.13
**Date:** 2025-12-11

This document describes the available solver methods in Radia.

## Available Methods

| Method | Name | Complexity | Best For |
|--------|------|------------|----------|
| **LU Direct** | `'lu'` or `'direct'` or `0` | O(N^3) | Small problems (N < 500), guaranteed convergence |
| **BiCGSTAB** (Default) | `'bicgstab'` or `'iterative'` or `1` | O(N^2 * k) | General purpose, large problems |

**Note:** Both solvers support linear and nonlinear materials. They use pure Newton-Raphson iteration for nonlinear convergence.

## Solver Architecture (v1.3.13)

### Historical Background

- **Original Radia**: Used Implicit SS (Successive Substitution / Gauss-Seidel) method
  - Slow convergence for high-permeability nonlinear materials
  - Could require hundreds of iterations

- **Current Radia (v1.3.13+)**: Replaced with modern solvers
  - **BiCGSTAB** iterative solver for better convergence
  - **LU direct** solver for guaranteed convergence
  - Both use pure Newton-Raphson iteration (no Gauss-Seidel M(H) correction)

### Performance (v1.3.13)

After OpenBLAS and OpenMP optimization:

| N_elem | LU Time | BiCGSTAB Time | Notes |
|--------|---------|---------------|-------|
| 104 | 0.13s | 0.13s | Linear material |
| 200 | 0.44s | 0.43s | Linear material |
| 390 | 1.72s | 1.71s | Linear material |

**Comparison with ELF_MAGIC (nonlinear material):**
- BiCGSTAB: **0.4x-0.6x** ratio (Radia is faster)
- LU: 7-11x ratio (more iterations in Radia due to stricter tolerance)

## Usage

### Basic Usage

```python
import radia as rad

# Create geometry
cube = rad.ObjRecMag([0, 0, 0], [40, 40, 40], [0, 0, 0])
rad.ObjDivMag(cube, [5, 5, 5])
mat = rad.MatLin(999.0)  # mu_r = 1000
rad.MatApl(cube, mat)

# Apply external field
ext_field = rad.ObjBckg([0, 0, 1.0])  # 1 T in SI
grp = rad.ObjCnt([cube, ext_field])

# Solve using default method (BiCGSTAB)
res = rad.Solve(grp, 0.0001, 1000)

# Or specify method by name
res = rad.Solve(grp, 0.0001, 1000, 'bicgstab')    # BiCGSTAB (default)
res = rad.Solve(grp, 0.0001, 1000, 'lu')          # LU decomposition
```

### Method Selection Guide

```
Problem size?
  |-- N < 500    --> Either method (LU is more robust)
  |-- N >= 500   --> 'bicgstab' (faster)

Material type?
  |-- Linear (MatLin)      --> Any solver works
  |-- Nonlinear (MatSatIso, etc.) --> Any solver works
```

## LU Direct Solver (Method 0)

Direct solver using LU decomposition with LAPACK `dgesv`.

**Pros:**
- Exact solution per nonlinear iteration
- Always converges (no divergence risk)
- Stable for all materials

**Cons:**
- O(N^3) time complexity
- O(N^2) memory usage

**Best for:** Small problems (N < 500), validation/debugging

```python
res = rad.Solve(grp, 0.0001, 100, 'lu')     # By name
res = rad.Solve(grp, 0.0001, 100, 0)        # By number
```

## BiCGSTAB Solver (Method 1, Default)

BiCGSTAB (Biconjugate Gradient Stabilized) iterative solver with Jacobi preconditioning.

**Pros:**
- O(N^2 * k) time complexity (k = iterations)
- Fast for medium to large problems
- Good convergence with preconditioning

**Cons:**
- May not converge for ill-conditioned problems

**Best for:** General magnetostatic problems, tetrahedral/hexahedral meshes

```python
res = rad.Solve(grp, 0.0001, 1000)              # Default (BiCGSTAB)
res = rad.Solve(grp, 0.0001, 1000, 'bicgstab')  # By name
res = rad.Solve(grp, 0.0001, 1000, 1)           # By number
```

## Convergence Tolerance

The `PrecOnMagnetiz` parameter controls convergence:

```python
rad.Solve(grp, tol, max_iter, method)
#              ^^^
#              Relative tolerance: ||dM||/||M||
```

**Comparison with ELF_MAGIC:**
- ELF default: `0.01` (1%)
- Radia benchmark: `0.0001` (0.01%)

**Recommendation:** Use `0.01` for typical applications, `0.0001` for high precision.

## H-Matrix Acceleration

**Status:** Under research (not available in v1.3.13)

H-matrix (HACApK) acceleration was evaluated but found to provide **no benefit for typical Radia use cases** (single compact objects). This is because:

1. All elements are spatially close together
2. No blocks satisfy the admissibility criterion
3. All blocks remain dense (no compression benefit)

See [HMATRIX_EVALUATION.md](HMATRIX_EVALUATION.md) for details.

**Future:** H-matrix may be beneficial for:
- Multiple well-separated magnetic objects
- Large-scale problems with distributed geometry

## Technical Details

### System Equation

The MMM (Magnetic Moment Method) system equation:

```
(1/chi - N) * M = H_ext
```

Where:
- `chi`: Magnetic susceptibility tensor
- `N`: Interaction matrix (demagnetization coefficients)
- `M`: Magnetization vector
- `H_ext`: External field

### Nonlinear Iteration

For nonlinear materials, outer Newton-Raphson iteration:

1. Compute `chi(H)` from current field estimate
2. Solve linear system with current `chi`
3. Update magnetization: `M_new = solution`
4. Check convergence: `||M_new - M_old|| / ||M_new|| < tol`
5. Repeat until converged

### BLAS/LAPACK Optimization (v1.3.13)

- `cblas_ddot`, `cblas_dnrm2`, `cblas_daxpy`: Vector operations
- `cblas_dgemv`: Matrix-vector product
- `dgesv_`: LU decomposition with partial pivoting

### OpenMP Parallelization (v1.3.13)

- Interaction matrix O(N^2) construction is parallelized
- Speedup: Up to 240x for large problems

## Migration Notes

### From v1.3.6 or earlier

| Old Method | New Method | Notes |
|------------|------------|-------|
| `8` (Newton) | Removed | Use `0` (LU) or `1` (BiCGSTAB) |
| `9` (LU) | `0` or `'lu'` | Same functionality |
| `10` (BiCGSTAB) | `1` or `'bicgstab'` | Same functionality |

### From v1.3.12 (Implicit SS removal)

Method 2 (Implicit SS / Gauss-Seidel) was removed due to slow convergence for nonlinear materials. Use BiCGSTAB (Method 1) instead.

---

**Last Updated:** 2025-12-11
**Project:** Radia Magnetic Field Computation
