# Radia Solver Methods

This document describes the available solver methods in Radia.

## Available Methods

| Method | Name | Complexity | Linear | Nonlinear | Best For |
|--------|------|------------|--------|-----------|----------|
| **Newton-Raphson** | `'newton'` or `8` | O(N^2 * k) | Yes | Yes | Nonlinear materials |
| **LU Direct** | `'lu'` or `9` | O(N^3 * k) | Yes | Yes | Small problems (N < 500) |
| **BiCGSTAB** (Default) | `'bicgstab'` or `10` | O(N^2 * k) | Yes | Yes | General purpose |
| **BiCGSTAB + H-matrix** | `'bicgstab'` + `SolverHMatrixEnable()` | O(N log N * k) | Yes | Yes | Large problems (N > 1000) |

**Note:** All solvers now support both linear and nonlinear materials (v1.3.5+).

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
res = rad.Solve(grp, 0.0001, 1000, 'bicgstab')  # BiCGSTAB
res = rad.Solve(grp, 0.0001, 1000, 'lu')        # LU decomposition
res = rad.Solve(grp, 0.0001, 1000, 'newton')    # Newton-Raphson (nonlinear)

# Or by number (for backward compatibility)
res = rad.Solve(grp, 0.0001, 1000, 10)  # BiCGSTAB
res = rad.Solve(grp, 0.0001, 1, 9)      # LU (only 1 iteration needed)
res = rad.Solve(grp, 0.0001, 1000, 8)   # Newton-Raphson
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
                                                    'newton' is traditional choice
                                                    'bicgstab' also works well
```

**Note:** All solvers now have outer nonlinear iteration loops, so they all
handle nonlinear materials correctly. The choice depends mainly on problem size.

## Newton-Raphson (Method 8)

Newton-Raphson iterative solver for **nonlinear materials**. Uses local Jacobian for each element.

**Pros:**
- Handles nonlinear (saturable) materials correctly
- Uses instantaneous susceptibility at each iteration
- Good convergence for typical B-H curves

**Cons:**
- Slower than BiCGSTAB for linear materials
- Requires well-defined B-H curve

**Best for:** Nonlinear materials (MatSatIso, MatSatIsoTab, MatLam, etc.)

```python
# Nonlinear material with B-H curve
bh_curve = [
    [0, 0],
    [100, 0.5],
    [500, 1.2],
    [2000, 1.6],
    [10000, 1.9]
]
mat = rad.MatSatIsoTab(bh_curve)
rad.MatApl(cube, mat)

# Solve with Newton-Raphson
res = rad.Solve(grp, 0.0001, 1000, 'newton')
```

## BiCGSTAB (Default)

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
res = rad.Solve(grp, 0.0001, 1000)        # Uses default BiCGSTAB
res = rad.Solve(grp, 0.0001, 1000, 'bicgstab')  # Explicit
```

## LU Direct Solver

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
res = rad.Solve(grp, 0.0001, 100, 'lu')   # For nonlinear materials
res = rad.Solve(grp, 0.0001, 100, 9)      # Same as above
```

**Note:** For linear materials, LU converges in 1 outer iteration. For nonlinear
materials, multiple outer iterations are needed for chi(H) to converge.

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

Both LU and BiCGSTAB methods produce consistent results for linear materials:

| N_elem | LU Bz (T) | BiCGSTAB Bz (T) | Difference |
|--------|-----------|-----------------|------------|
| 27 | 0.1373162 | 0.1373157 | 0.0004% |
| 512 | 0.1406983 | 0.1406980 | 0.0002% |
| 1000 | 0.1411381 | 0.1411383 | 0.0001% |

## Notes

1. **Default change (v1.3.5):** BiCGSTAB is now the default solver
2. **Nonlinear support (v1.3.5):** All solvers now use Newton-style M(H) updates and produce identical results
3. **Tetrahedral meshes:** All methods work correctly with tetrahedral elements
4. **Material types:**
   - Linear materials (MatLin): Any solver works; 'bicgstab' is fastest for large problems
   - Nonlinear materials (MatSatIso, MatSatIsoTab, MatLam): Any solver works; all produce identical results
5. **Method numbers:** For backward compatibility, methods can also be specified by number (8=Newton, 9=LU, 10=BiCGSTAB)

## Technical Details: Nonlinear Material Handling

**All solvers (v1.3.5+)** now use Newton-style M(H) updates for nonlinear materials:

1. **Outer iteration loop**: After each linear system solve, apply Newton-style correction
2. **Gauss-Seidel update**: For each element i:
   - Compute quasi-external field: sum of contributions from all OTHER elements + external field
   - Solve local equation: H = (I - chi*Nii)^{-1} * (QuasiExtField + Mr)
   - Apply material's M(H) function directly: M = M(H)
3. **Convergence check**: Monitor change in magnetization between iterations

This hybrid approach combines the efficiency of LU/BiCGSTAB with the accuracy of Newton-Raphson:
- LU/BiCGSTAB provide a good initial guess for M
- Newton-style M(H) update ensures correct nonlinear behavior
- All three solvers now produce identical results for both linear and nonlinear materials

**Solver Comparison** (v1.3.5+):

| Solver | Inner Method | Nonlinear Update | Best For |
|--------|--------------|------------------|----------|
| Newton (8) | Gauss-Seidel | M = M(H) | General purpose |
| LU (9) | LU decomposition | M = M(H) | Small problems, validation |
| BiCGSTAB (10) | BiCGSTAB iteration | M = M(H) | Large problems |

All solvers produce identical results for both linear and nonlinear materials.
