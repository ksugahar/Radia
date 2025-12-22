# Nonlinear Material Benchmark

Benchmark for Radia's nonlinear solver with saturable BH curves using MSC (Magnetic Surface Charge) method.

## Problem Setup

- **Geometry**: 1.0 m x 1.0 m x 1.0 m cube centered at origin
- **Material**: Nonlinear BH curve (soft iron with saturation)
- **External Field**: Hz = 200,000 A/m (uniform field along z-axis)

## Element Types

| Element | Mesh Tool | DOF | Benchmark Script |
|---------|-----------|-----|------------------|
| Hexahedron MSC | ObjPolyhdr | 6 (sigma per face) | `benchmark_hexahedron_msc.py` |
| Tetrahedron MSC | Netgen | 3 (Mx, My, Mz) | `benchmark_tetrahedron_msc_netgen.py` |

## BH Curve (Soft Iron)

```
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
```

- Initial permeability: mu_r ~= 800
- Saturation onset: B ~= 1.5 T

## Hexahedron MSC Benchmarks

Results stored in `hexahedron_msc/{lu,bicgstab}/` directories.

### LU Solver Results (H_ext = 200,000 A/m)

| N | Elements | DOF | Time (s) | Iter | M_avg_z (A/m) |
|---|----------|-----|----------|------|---------------|
| 3 | 27 | 162 | 0.049 | 5 | 679,134 |
| 5 | 125 | 750 | 0.299 | 6 | 702,129 |
| 8 | 512 | 3072 | 7.284 | 7 | 713,064 |
| 10 | 1000 | 6000 | 53.58 | 13 | 716,275 |

### BiCGSTAB Solver Results (H_ext = 200,000 A/m)

| N | Elements | DOF | Time (s) | Iter | M_avg_z (A/m) |
|---|----------|-----|----------|------|---------------|
| 3 | 27 | 162 | 0.036 | 4 | 678,577 |
| 5 | 125 | 750 | 1.135 | 8 | 701,259 |

**Note**: BiCGSTAB has convergence issues for N >= 8 with nonlinear materials.

## Tetrahedron MSC Benchmarks

Results stored in `tetrahedron_msc/{lu,bicgstab}/` directories.

### LU Solver Results (H_ext = 50,000 A/m)

| maxh | Elements | DOF | Time (s) | Iter | M_avg_z (A/m) |
|------|----------|-----|----------|------|---------------|
| 0.50 | 28 | 84 | 0.082 | 2 | 184,790 |
| 0.40 | 104 | 312 | 0.385 | 2 | 193,230 |
| 0.35 | 112 | 336 | 0.442 | 2 | 197,377 |
| 0.30 | 200 | 600 | 0.950 | 2 | 190,594 |

---

## Solver Methods

Radia provides two solver methods:

| Method | Description | Recommendation |
|--------|-------------|----------------|
| 0 (LU) | Direct LU decomposition | Small problems (N < 12) |
| 1 (BiCGSTAB) | Iterative solver (default) | Large problems (N >= 12) |

```python
# Method 0: Direct LU solver
res = rad.Solve(grp, 0.001, 1000, 0)

# Method 1: BiCGSTAB iterative solver (default)
res = rad.Solve(grp, 0.001, 1000, 1)
```

### When to Use Each Method

| Problem Size | Recommended Solver | Rationale |
|--------------|-------------------|-----------|
| N < 12 (< 1,728 elem) | LU (Method 0) | Fast, direct solution |
| N >= 12 | BiCGSTAB (Method 1) | Better O(N^2) vs O(N^3) scaling |
| Distributed geometry | BiCGSTAB + H-matrix | O(N log N) for well-separated clusters |

---

**Last Updated**: 2025-12-22
**Author**: Claude Code
