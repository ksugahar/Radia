# Nonlinear Material Benchmark

Benchmark for Radia's nonlinear solver with saturable BH curves.

## Problem Setup

- **Geometry**: 1.0 m x 1.0 m x 1.0 m cube centered at origin
- **Mesh**: N x N x N hexahedral elements
- **Material**: Nonlinear BH curve (soft iron with saturation)
- **External Field**: Hz = 50,000 A/m (uniform field along z-axis)

## BH Curve (Soft Iron)

H [A/m]     B [T]
0           0.0
100         0.1
200         0.3
500         0.8
1000        1.2
2000        1.5
5000        1.7
10000       1.8
50000       2.0
100000      2.1

- Initial permeability: mu_r ~= 800
- Saturation onset: B ~= 1.5 T

## BH Curve Definition in Radia

The nonlinear material is defined using rad.MatSatIsoTab() with [H, M] pairs:

import numpy as np
import radia as rad

# BH curve data: [H (A/m), B (T)]
bh_data = [
    [0.0, 0.0],
    [100.0, 0.1],
    [200.0, 0.3],
    [500.0, 0.8],
    [1000.0, 1.2],
    [2000.0, 1.5],
    [5000.0, 1.7],
    [10000.0, 1.8],
    [50000.0, 2.0],
    [100000.0, 2.1],
]

# Convert to [H, M] format
# B = mu_0 * (H + M), so M = B/mu_0 - H
mu_0 = 4 * np.pi * 1e-7
hm_data = [[h, b/mu_0 - h] for h, b in bh_data]

mat = rad.MatSatIsoTab(hm_data)

## Solver Methods

Radia provides two solver methods:

| Method | Description | Recommendation |
|--------|-------------|----------------|
| 0 (LU) | Direct LU decomposition | Small problems (N < 500) |
| 1 (BiCGSTAB) | Iterative solver (default) | Large problems (N >= 500) |

# Method 0: Direct LU solver
res = rad.Solve(grp, 0.001, 1000, 0)

# Method 1: BiCGSTAB iterative solver (default)
res = rad.Solve(grp, 0.001, 1000, 1)

## Running the Benchmark

cd examples/cube_uniform_field/nonlinear
python benchmark_nonlinear_tetra_vs_hex.py

## Files

- benchmark_nonlinear_tetra_vs_hex.py - Tetrahedral vs hexahedral comparison
- README.md - This documentation

## Physics Notes

### Nonlinear Iteration

For nonlinear materials, Radia uses an iterative method:

1. Initialize magnetization M = 0
2. Compute H field from current M distribution
3. Update M using material BH curve: M = M(H)
4. Repeat until |delta M| < tolerance

### Saturation Effects

At high H fields, the material saturates:
- Low H: B ~ mu_r * mu_0 * H (linear regime)
- High H: B -> B_sat (saturation regime)
- Saturation magnetization for iron: ~2.1 T
