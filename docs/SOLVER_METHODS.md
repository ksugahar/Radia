# Radia Solver Methods

This document describes the available solver methods in Radia.

## Available Methods (v1.3.7+)

| Method | Name | Complexity | Linear | Nonlinear | Best For |
|--------|------|------------|--------|-----------|----------|
| **LU Direct** | `'lu'` or `'direct'` or `0` | O(N^3 * k) | Yes | Yes | Small problems (N < 500) |
| **BiCGSTAB** (Default) | `'bicgstab'` or `'iterative'` or `1` | O(N^2 * k) | Yes | Yes | General purpose |
| **BiCGSTAB + H-matrix** | `'bicgstab'` + `SolverHMatrixEnable()` | O(N log N * k) | Yes | Yes | Large problems (N > 1000) |

**Note:** All solvers support both linear and nonlinear materials. The Newton-Raphson method (former Method 8) has been removed in v1.3.7 - Newton-style M(H) updates are now integrated into both LU and BiCGSTAB solvers.

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
ext_field = rad.ObjBckg([0, 0, 1.0])  # 1 T
grp = rad.ObjCnt([cube, ext_field])

# Solve using default method (BiCGSTAB)
res = rad.Solve(grp, 0.0001, 1000)

# Or specify method by name
res = rad.Solve(grp, 0.0001, 1000, 'bicgstab')    # BiCGSTAB (default)
res = rad.Solve(grp, 0.0001, 1000, 'iterative')   # Same as 'bicgstab'
res = rad.Solve(grp, 0.0001, 1000, 'lu')          # LU decomposition
res = rad.Solve(grp, 0.0001, 1000, 'direct')      # Same as 'lu'

# Or by number
res = rad.Solve(grp, 0.0001, 1000, 0)  # LU Direct
res = rad.Solve(grp, 0.0001, 1000, 1)  # BiCGSTAB
```

### Method Selection Guide

```
Problem size?
  |-- N < 500    --> 'lu' (LU Direct) or 'bicgstab' (default)
  |-- N >= 500   --> 'bicgstab' (BiCGSTAB, default)
  |-- N > 1000   --> 'bicgstab' + SolverHMatrixEnable()

Material type?
  |-- Linear (MatLin)      --> Any solver works
  |-- Nonlinear (MatSatIso, MatSatIsoTab, etc.) --> Any solver works
```

**Note:** Both solvers have outer nonlinear iteration loops with Newton-style M(H) updates, so they handle nonlinear materials correctly. The choice depends mainly on problem size.

## LU Direct Solver (Method 0)

Direct solver using LU decomposition with partial pivoting. O(N^3) complexity per nonlinear iteration.

**Pros:**
- Exact solution per nonlinear iteration
- Always converges
- Stable for all materials
- Supports both linear and nonlinear materials

**Cons:**
- O(N^3) time complexity - slow for large N
- O(N^2) memory usage

**Best for:** Small problems (N < 500), validation/debugging

```python
res = rad.Solve(grp, 0.0001, 100, 'lu')     # By name
res = rad.Solve(grp, 0.0001, 100, 'direct') # Alias
res = rad.Solve(grp, 0.0001, 100, 0)        # By number
```

**Note:** For linear materials, LU converges in 1-2 outer iterations. For nonlinear materials, multiple outer iterations are needed for chi(H) to converge.

## BiCGSTAB (Method 1, Default)

BiCGSTAB (Biconjugate Gradient Stabilized) is an iterative solver with O(N^2 * k) complexity where k is the number of iterations.

**Pros:**
- Fast for medium to large problems
- Stable for high permeability materials
- Good convergence with Jacobi preconditioning
- Supports both linear and nonlinear materials

**Cons:**
- May not converge for ill-conditioned problems

**Best for:** General magnetostatic problems, tetrahedral meshes

```python
res = rad.Solve(grp, 0.0001, 1000)              # Default (BiCGSTAB)
res = rad.Solve(grp, 0.0001, 1000, 'bicgstab')  # By name
res = rad.Solve(grp, 0.0001, 1000, 'iterative') # Alias
res = rad.Solve(grp, 0.0001, 1000, 1)           # By number
```

## H-Matrix Acceleration

Enable H-matrix with HACApK ACA+ algorithm for BiCGSTAB:

```python
# Enable H-matrix
rad.SolverHMatrixEnable()
res = rad.Solve(grp, 0.0001, 1000, 'bicgstab')

# Disable H-matrix
rad.SolverHMatrixDisable()

# Check status
status = rad.SolverHMatrixStatus()  # 1 if enabled, 0 if disabled
```

**Pros:**
- O(N log N) per iteration instead of O(N^2)
- Reduced memory for large problems

**Cons:**
- Overhead for small problems

**Best for:** Large problems (N > 1000)

## Performance Benchmark

Results from 40mm soft iron cube (mu_r=1000) in 1T uniform field:

| N_elem | LU Time | BiCGSTAB Time | BiCGSTAB Iters |
|--------|---------|---------------|----------------|
| 27 | 0.005s | 0.0003s | 6 |
| 125 | 0.010s | 0.005s | 12 |
| 512 | 0.72s | 0.10s | 14 |
| 1000 | 6.58s | 0.39s | 16 |
| 1728 | 34.1s | 1.25s | 18 |

## Accuracy

Both LU and BiCGSTAB methods produce identical results:

| N_elem | LU Bz (T) | BiCGSTAB Bz (T) | Difference |
|--------|-----------|-----------------|------------|
| 27 | 14.067658 | 14.067658 | 0.0000% |
| 64 | 14.067658 | 14.067658 | 0.0000% |
| 125 | 14.067658 | 14.067658 | 0.0000% |

## Notes

1. **Simplified method numbering (v1.3.7):** Methods are now 0 (LU) and 1 (BiCGSTAB)
2. **Newton-Raphson removed (v1.3.7):** Newton-style M(H) updates are integrated into both solvers
3. **Default solver:** BiCGSTAB is the default
4. **Tetrahedral meshes:** All methods work correctly with tetrahedral elements
5. **Material types:**
   - Linear materials (MatLin): Any solver works; 'bicgstab' is fastest for large problems
   - Nonlinear materials (MatSatIso, MatSatIsoTab, MatLam): Any solver works; all produce identical results

## Why BiCGSTAB and H-matrix are Efficient for Nonlinear Materials

**Key insight: Nonlinearity only affects diagonal blocks (self-interaction terms N_ii).**

In the MMM (Magnetic Moment Method) system equation:

```
M = chi(H) * (N * M + H_ext)
```

The interaction matrix N has a special structure:
- **Off-diagonal blocks N_ij (i != j)**: Depend only on geometry - constant throughout solution
- **Diagonal blocks N_ii (self-interaction)**: Only these terms interact with chi(H)

This means:
1. **H-matrix compression remains valid**: Off-diagonal blocks can still be compressed with ACA+
2. **BiCGSTAB matvec is unchanged**: Same interaction matrix used for all iterations
3. **Nonlinear update is O(N)**: Only chi * N_ii diagonal update is nonlinear, not O(N^2)

**Practical implications:**
- H-matrix acceleration works for nonlinear materials (same compression as linear case)
- BiCGSTAB convergence is similar for linear and nonlinear materials
- HACApK-BiCGSTAB combination is efficient for large nonlinear problems

**Reference:** This property is well-known in computational electromagnetics and is exploited by many fast solvers including HACApK.

## Technical Details: Nonlinear Material Handling

Both solvers use Newton-style M(H) updates for nonlinear materials:

1. **Outer iteration loop**: After each linear system solve, apply Newton-style correction
2. **Gauss-Seidel update**: For each element i:
   - Compute quasi-external field: sum of contributions from all OTHER elements + external field
   - Solve local equation: H = (I - chi*Nii)^{-1} * (QuasiExtField + Mr)
   - Apply material's M(H) function directly: M = M(H)
3. **Convergence check**: Monitor change in magnetization between iterations

This hybrid approach combines the efficiency of LU/BiCGSTAB with the accuracy of Newton-Raphson:
- LU/BiCGSTAB provide a good initial guess for M
- Newton-style M(H) update ensures correct nonlinear behavior
- Both solvers produce identical results for both linear and nonlinear materials

**Solver Comparison** (v1.3.7+):

| Method | Number | Inner Method | Nonlinear Update | Best For |
|--------|--------|--------------|------------------|----------|
| LU | 0 | LU decomposition | M = M(H) | Small problems, validation |
| BiCGSTAB | 1 | BiCGSTAB iteration | M = M(H) | General purpose, large problems |

Both solvers produce identical results for both linear and nonlinear materials.

## Migration from v1.3.6

If you were using method numbers 8, 9, or 10:

| Old (v1.3.6) | New (v1.3.7+) | Notes |
|--------------|---------------|-------|
| `8` (Newton) | Removed | Use `0` (LU) or `1` (BiCGSTAB) - both have Newton-style M(H) updates |
| `9` (LU) | `0` or `'lu'` | Same functionality |
| `10` (BiCGSTAB) | `1` or `'bicgstab'` | Same functionality, now default |

**Recommended:** Use string names (`'lu'`, `'bicgstab'`, `'direct'`, `'iterative'`) for clarity.
