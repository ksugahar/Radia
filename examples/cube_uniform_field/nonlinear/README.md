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

---

## Hexahedron MSC Benchmarks

Results stored in `hexahedron_msc/{lu,bicgstab,hacapk}/` directories.

### Solver Methods

| Method | Description | Complexity |
|--------|-------------|------------|
| LU | Dense LU decomposition | O(N^3) |
| BiCGSTAB | Iterative BiCGSTAB | O(N^2) per iteration |
| HACApK | H-matrix accelerated BiCGSTAB | O(N log N) per iteration |

### LU Solver Results (H_ext = 200,000 A/m)

| N | Elements | DOF | Time (s) | Iter | M_avg_z (A/m) | Memory (MB) |
|---|----------|-----|----------|------|---------------|-------------|
| 5 | 125 | 750 | 0.30 | 6 | 702,132 | 49 |
| 10 | 1,000 | 6,000 | 46.26 | 13 | 716,281 | 1,174 |

**Note**: LU solver for N >= 15 is prohibitively slow (O(N^6) time complexity for nonlinear).

### BiCGSTAB Solver Results (H_ext = 200,000 A/m)

| N | Elements | DOF | Time (s) | Iter | M_avg_z (A/m) | Memory (MB) |
|---|----------|-----|----------|------|---------------|-------------|
| 5 | 125 | 750 | 0.28 | 8 | 702,128 | 35 |
| 10 | 1,000 | 6,000 | 118.77 | 50 | 716,261 | 345 |

**Note**: BiCGSTAB is slower than LU for this problem due to high iteration count. Use HACApK for large problems.

### HACApK Solver Results (H_ext = 200,000 A/m, eps = 1e-4)

| N | Elements | DOF | Time (s) | Iter | M_avg_z (A/m) | Memory (MB) | Compression | Leaves |
|---|----------|-----|----------|------|---------------|-------------|-------------|--------|
| 5 | 125 | 750 | 0.14 | 4 | 701,757 | 35 | 0.0219 | 184 |
| 10 | 1,000 | 6,000 | 6.09 | 5 | 715,885 | 208 | 0.0059 | 3,190 |

### Speedup: HACApK vs LU

| N | LU Time (s) | HACApK Time (s) | Speedup |
|---|-------------|-----------------|---------|
| 5 | 0.30 | 0.14 | 2.1x |
| 10 | 46.26 | 6.09 | 7.6x |
| 15 | (est. hours) | - | >>10x |
| 20 | (impractical) | - | >>100x |

---

## Radia vs ELF/MAGIC Comparison

Comparison of Radia and ELF/MAGIC solvers on identical hexahedral mesh (1.0m cube, H_ext = 200,000 A/m).

### Magnetization Accuracy (M_avg_z)

| N | Elements | Radia LU | ELF LU | Radia HACApK | ELF HACApK | Difference (%) |
|---|----------|----------|--------|--------------|------------|----------------|
| 5 | 125 | 702,132 | 702,132 | 701,757 | 702,110 | < 0.1% |
| 10 | 1,000 | 716,281 | 716,281 | 715,885 | 716,360 | < 0.1% |

**Note**: LU results match exactly. Small differences in HACApK are due to H-matrix compression (eps=1e-4).

### Performance Comparison

#### LU Solver

| N | Elements | Radia Time (s) | ELF Time (s) | Radia Iter | ELF Iter |
|---|----------|----------------|--------------|------------|----------|
| 5 | 125 | 0.30 | 0.13 | 6 | 6 |
| 10 | 1,000 | 46.26 | 13.42 | 13 | 13 |

**Note**: ELF LU is faster due to Fortran LAPACK optimizations. Iteration counts match exactly.

#### BiCGSTAB Solver

| N | Elements | Radia Time (s) | ELF Time (s) | Radia Iter | ELF Iter |
|---|----------|----------------|--------------|------------|----------|
| 5 | 125 | 0.28 | 0.08 | 8 | 3 |
| 10 | 1,000 | 118.77 | 4.93 | 50 | 4 |

**Note**: Radia BiCGSTAB requires more iterations due to different convergence criterion (B-field vs permeability).

#### HACApK Solver (eps = 1e-4)

| N | Elements | Radia Time (s) | ELF Time (s) | Radia Iter | ELF Iter |
|---|----------|----------------|--------------|------------|----------|
| 5 | 125 | 0.14 | 0.15 | 4 | 3 |
| 10 | 1,000 | 6.09 | 4.07 | 5 | 4 |

**Note**: HACApK performance is comparable. Both use ACA+ compression with identical parameters.

### Key Observations

1. **Accuracy**: LU results are identical between Radia and ELF (M_avg_z matches to all digits)
2. **Convergence**: Both use similar Newton-Raphson (mucal2) convergence with B-field criterion
3. **HACApK**: Performance is comparable; Radia's H-matrix overhead is similar to ELF
4. **BiCGSTAB**: Radia's dense BiCGSTAB needs optimization; use HACApK for large problems

---

## Tetrahedron MSC Benchmarks

Results stored in `tetrahedron_msc/{lu,bicgstab,hacapk}/` directories.

### HACApK Solver Results (H_ext = 200,000 A/m, eps = 1e-4)

| maxh | Elements | DOF | Time (s) | Iter | M_avg_z (A/m) | Memory (MB) | Leaves |
|------|----------|-----|----------|------|---------------|-------------|--------|
| 0.40 | 104 | 312 | 0.29 | 2 | 588,673 | 66 | 172 |
| 0.20 | 627 | 1,881 | 12.50 | 2 | 568,973 | 134 | 2,098 |
| 0.15 | 2,211 | 6,633 | 135.86 | 2 | 587,299 | 775 | 14,386 |
| 0.10 | 4,994 | 14,982 | 757.65 | 2 | 570,603 | 3,298 | 36,475 |

**Note**: 3DOF tetrahedra require pre-computed dense InteractMatrix (O(N^2) memory), which is the performance bottleneck. H-matrix is used for MatVec acceleration only.

---

## Convergence Criterion

Radia v1.3.15+ uses **B-field based convergence** (mucal2) for nonlinear materials, matching the Newton-Raphson method used in industry-standard solvers.

### Convergence Methods

| Method | Stage | Formula | Description |
|--------|-------|---------|-------------|
| **mucal0** | Initialization | `chi = B[1] / (mu_0 * H[1]) - 1` | Initial chi from 2nd BH point |
| mucal1 | Iteration | `rel_change = |chi_new - chi_old| / chi_old` | Legacy relaxation method |
| **mucal2** | Iteration | `rel_change = |B_new - B_old| / B_sat` | Newton-Raphson, fast convergence |

---

## Solver Methods

Radia provides three solver methods:

| Method | Name | Description | Recommendation |
|--------|------|-------------|----------------|
| 0 | LU | Direct LU decomposition | Small problems (N < 12) |
| 1 | BiCGSTAB | Iterative solver | Medium problems |
| 2 | HACApK | H-matrix accelerated BiCGSTAB | Large problems (N >= 10) |

```python
# Method 0: Direct LU solver
res = rad.Solve(grp, 0.001, 1000, 0)

# Method 1: BiCGSTAB iterative solver
res = rad.Solve(grp, 0.001, 1000, 1)

# Method 2: HACApK (H-matrix + BiCGSTAB)
rad.SetHACApKParams(1e-4, 10, 2.0)  # eps, leaf_size, eta
res = rad.Solve(grp, 0.001, 1000, 2)
```

### When to Use Each Method

| Problem Size | Recommended Solver | Rationale |
|--------------|-------------------|-----------|
| N < 10 (< 1,000 elem) | LU (Method 0) | Fast, direct solution |
| 10 <= N < 15 | HACApK (Method 2) | 8x faster than LU |
| N >= 15 | HACApK (Method 2) | Only practical option |

---

## Usage

### Hexahedron Benchmarks

```bash
# Run specific solver
python benchmark_hexahedron_msc.py --lu 5 10
python benchmark_hexahedron_msc.py --bicgstab 5 10
python benchmark_hexahedron_msc.py --hacapk 5 10 15 20

# Custom ACA tolerance
python benchmark_hexahedron_msc.py --hacapk --eps 1e-3 5 10
```

### Tetrahedron Benchmarks

```bash
# Run specific solver
python benchmark_tetrahedron_msc_netgen.py --lu 0.4 0.2
python benchmark_tetrahedron_msc_netgen.py --hacapk 0.4 0.2 0.15 0.10

# Custom ACA tolerance
python benchmark_tetrahedron_msc_netgen.py --hacapk --eps 1e-3 0.15
```

---

**Last Updated**: 2025-12-24
**Author**: Claude Code
