# Nonlinear Material Benchmark

Benchmark for Radia's nonlinear solver with saturable BH curves using MSC (Magnetic Surface Charge) method.

## Problem Setup

- **Geometry**: 1.0 m x 1.0 m x 1.0 m cube centered at origin
- **Material**: Nonlinear BH curve (soft iron with saturation)
- **External Field**: Hz = 50,000 A/m (uniform field along z-axis)

## Element Types

| Element | Mesh Tool | DOF | Benchmark Script |
|---------|-----------|-----|------------------|
| Hexahedron MSC | ObjPolyhdr | 3 (Mx, My, Mz) | `benchmark_hexahedron_msc.py` |
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

### LU Solver Results

| N | Elements | DOF | Time (s) | Iter | M_avg_z (A/m) |
|---|----------|-----|----------|------|---------------|
| 5 | 125 | 375 | 0.25 | 3 | 173,400 |

### BiCGSTAB Solver Results

| N | Elements | DOF | Time (s) | Iter | M_avg_z (A/m) |
|---|----------|-----|----------|------|---------------|
| 5 | 125 | 375 | 0.05 | 3 | 173,400 |

## Tetrahedron MSC Benchmarks

Results stored in `tetrahedron_msc/{lu,bicgstab}/` directories.

### LU Solver Results

| maxh | Elements | DOF | Time (s) | Iter | M_avg_z (A/m) |
|------|----------|-----|----------|------|---------------|
| 0.40 | ~60 | ~180 | 0.2 | 2 | ~175,000 |
| 0.30 | ~200 | ~600 | 1.5 | 2 | ~190,000 |
| 0.25 | ~350 | ~1050 | 4.5 | 2 | ~185,000 |

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

**Last Updated**: 2025-12-06
**Author**: Claude Code
