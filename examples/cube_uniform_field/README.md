# Cube in Uniform Field Benchmark

Radia magnetostatic solver benchmark for a soft iron cube in uniform external field.

## Test Environment

| Item | Specification |
|------|---------------|
| **CPU** | Intel Core i7-9700K (8 cores, 8 threads @ 3.60 GHz) |
| **Memory** | 128 GB DDR4 |
| **OS** | Windows 10/11 |
| **Parallelization** | NGSolve TaskManager (8 threads), MKL_NUM_THREADS=8 |
| **Compiler** | MSVC 2022 + Intel MKL |
| **Python** | 3.12 |

## Folder Structure

```
cube_uniform_field/
├── hexahedron/                 # Hexahedral benchmark
│   ├── benchmark_hex.py        # Hexahedral benchmark script
│   ├── linear/                 # Linear material results
│   │   ├── lu/
│   │   ├── bicgstab/
│   │   └── hacapk/
│   └── nonlinear/              # Nonlinear material results
│       ├── lu/
│       ├── bicgstab/
│       └── hacapk/
├── tetrahedron/                # Tetrahedral benchmark
│   ├── benchmark_tetra.py      # Tetrahedral benchmark script
│   ├── linear/                 # Linear material results
│   │   ├── lu/
│   │   ├── bicgstab/
│   │   └── hacapk/
│   └── nonlinear/              # Nonlinear material results
│       ├── lu/
│       ├── bicgstab/
│       └── hacapk/
├── benchmark_common.py         # Shared benchmark functions
├── run_benchmarks.py           # Run benchmarks from config in subprocesses for accurate memory measurement
├── compare_hmatrix_stats.py    # Compare H-matrix statistics between ELF and Radia leaf-by-leaf
├── profile_hacapk.py           # Profile HACApK H-matrix build to identify performance bottlenecks
├── test_objm_minimal.py        # Minimal test of ObjM with correct benchmark pattern (hex mesh)
├── benchmark_cube_block_jacobi.py  # Cube benchmark with Block Jacobi preconditioner and Newton-Raphson
└── README.md                   # This file
```

## Problem Description

- **Geometry**: 1.0 m x 1.0 m x 1.0 m soft iron cube (centered at origin)
- **External field**: H_z = 200,000 A/m
- **Material**:
  - Linear: Constant permeability mu_r = 1000 (chi = 999)
  - Nonlinear: Saturation BH curve (soft iron)

### Unified Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| H_ext | 200,000 A/m | External field |
| hmat_eps | 1e-4 | ACA+ compression tolerance |
| bicg_tol | 1e-4 | BiCGSTAB convergence |
| nonl_tol | 0.001 | Nonlinear convergence |

---

## Column Definitions

| Column | Description |
|--------|-------------|
| **M_avg_z** | Average magnetization in z-direction [A/m] |
| **Solver** | Solver method: LU, BiCGSTAB, or HACApK |
| **Memory** | Peak working set memory (Windows peak_wset) |
| **Compress** | H-matrix compression ratio (H-mat memory / Dense memory, HACApK only) |
| **Linear** | Number of linear (BiCGSTAB) iterations |
| **Nonl** | Number of nonlinear iterations |
| **MatBuild** | Matrix construction time: Dense matrix (LU/BiCGSTAB) or H-matrix with ACA+ (HACApK) |
| **LinSolve** | Linear system solve time |
| **Total** | Total solve time |

**Note**: For HACApK solver, "MatBuild" shows H-matrix construction time (including ACA+ fill). For LU/BiCGSTAB, it shows dense interaction matrix construction time.

---

## Benchmark Results (2025-12-31, 8-thread)

### 1. Hexahedral Linear Benchmark

#### N=5 (125 elements, 750 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 701,878 | LU | 51 MB | - | 0 | 1 | 0.07s | 0.40s | **0.48s** |
| 701,878 | BiCGSTAB | 46 MB | - | 0 | 1 | 0.03s | 3ms | **0.04s** |
| 701,878 | HACApK | 41 MB | 97% | 11 | 1 | 0.03s | 3ms | **0.04s** |

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 716,184 | LU | 933 MB | - | 0 | 1 | 1.0s | 1.0s | **2.1s** |
| 716,110 | BiCGSTAB | 632 MB | - | 0 | 1 | 1.0s | 0.2s | **1.3s** |
| 716,183 | HACApK | 193 MB | **50%** | 11 | 1 | 1.3s | 0.1s | **1.4s** |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 720,157 | LU | 9924 MB | - | 0 | 1 | 10.3s | 18.9s | **29.2s** |
| 720,159 | BiCGSTAB | 6723 MB | - | 0 | 1 | 10.2s | 4.7s | **15.0s** |
| 720,169 | HACApK | 911 MB | **26%** | 10 | 1 | 8.3s | 0.6s | **9.0s** |

#### N=20 (8,000 elements, 48,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| - | LU | - | - | - | - | - | - | **(~7min estimated)** |
| 721,938 | BiCGSTAB | 37481 MB | - | 0 | 1 | 57.7s | 24.9s | **82.6s** |
| 721,979 | HACApK | 2796 MB | **15%** | 16 | 1 | 31.7s | 3.4s | **35.2s** |

#### N=25 (15,625 elements, 93,750 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 723,047 | HACApK | 6300 MB | **9%** | 23 | 1 | 77.4s | 10.9s | **88.4s** |

---

### 2. Hexahedral Nonlinear Benchmark

**All solvers converge** for hexahedral nonlinear problems.

#### N=5 (125 elements, 750 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 702,132 | LU | 51 MB | - | 0 | 6 | 0.05s | 0.03s | **0.09s** |
| 702,148 | BiCGSTAB | 46 MB | - | 0 | 3 | 0.07s | 0.15s | **0.22s** |
| 702,143 | HACApK | 42 MB | 97% | 18 | 3 | 0.51s | 0.03s | **0.54s** |

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 716,281 | LU | 931 MB | - | 0 | 13 | 1.1s | 7.3s | **9.3s** |
| 716,334 | BiCGSTAB | 631 MB | - | 0 | 4 | 1.1s | 0.4s | **1.7s** |
| 716,335 | HACApK | 189 MB | **49%** | 19 | 4 | 1.4s | 0.3s | **1.7s** |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 719,832 | LU | 9925 MB | - | 0 | 35 | 12.5s | 726s | **767s** |
| 719,804 | BiCGSTAB | 6725 MB | - | 0 | 43 | 13.1s | 21.6s | **36s** |
| 719,809 | HACApK | 905 MB | **25%** | 72 | 43 | 9.4s | 6.4s | **16s** |

#### N=20 (8,000 elements, 48,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| - | LU | - | - | - | - | - | - | **(>1h estimated)** |
| 721,240 | BiCGSTAB | 37486 MB | - | 0 | 34 | 68.3s | 111s | **185s** |
| 721,285 | HACApK | 2794 MB | **15%** | 75 | 32 | 35.0s | 17.3s | **53s** |

#### N=25 (15,625 elements, 93,750 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 722,166 | HACApK | 6384 MB | **9%** | 148 | 26 | 79.1s | 67.9s | **148s** |

---

### 3. Tetrahedral Linear Benchmark

#### maxh=0.20m (627 elements, 1,881 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 816,133 | LU | 172 MB | - | 0 | 1 | 0.12s | 2.5s | 2.7s |
| 816,080 | BiCGSTAB | 140 MB | - | 0 | 1 | 0.12s | 0.08s | 0.4s |
| 816,089 | HACApK | 98 MB | 85% | 29 | 1 | 0.42s | 0.04s | **0.48s** |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 748,693 | LU | 1272 MB | - | 0 | 1 | 1.6s | 3.4s | 5.2s |
| 748,683 | BiCGSTAB | 922 MB | - | 0 | 1 | 1.6s | 1.0s | 3.0s |
| 748,705 | HACApK | 279 MB | **56%** | 25 | 1 | 1.71s | 0.33s | **2.06s** |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 754,617 | LU | 6105 MB | - | 0 | 1 | 8.9s | 34.6s | 45.6s |
| 754,394 | BiCGSTAB | 4378 MB | - | 0 | 1 | 9.0s | 4.6s | 15.4s |
| 754,549 | HACApK | 747 MB | **36%** | 29 | 1 | 6.41s | 1.29s | **7.77s** |

#### maxh=0.05m (33,974 elements, 101,922 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 737,397 | HACApK | 8758 MB | **10%** | 26 | 1 | 71.1s | 17.4s | **90.7s** |

---

### 4. Tetrahedral Nonlinear Benchmark

**All solvers converge** for tetrahedral nonlinear problems.

#### maxh=0.20m (627 elements, 1,881 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 748,822 | LU | 171 MB | - | 0 | 20 | 0.11s | 1.8s | 2.0s |
| 748,869 | BiCGSTAB | 140 MB | - | 0 | 14 | 0.12s | 0.24s | 0.4s |
| 748,941 | HACApK | 98 MB | 85% | 130 | 17 | 0.32s | 0.15s | **0.51s** |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 730,715 | LU | 1273 MB | - | 0 | 29 | 1.6s | 87.4s | 91.5s |
| 730,760 | BiCGSTAB | 922 MB | - | 0 | 21 | 1.7s | 4.0s | 6.1s |
| 730,596 | HACApK | 280 MB | **56%** | 134 | 36 | 1.70s | 2.0s | **3.72s** |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 730,996 | LU | 6108 MB | - | 0 | 39 | 9.0s | 1315s | 1336s |
| 731,003 | BiCGSTAB | 4380 MB | - | 0 | 41 | 9.3s | 36.3s | 47.7s |
| 731,063 | HACApK | 751 MB | **36%** | 206 | 41 | 6.42s | 9.9s | **16.7s** |

#### maxh=0.05m (33,974 elements, 101,922 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|------:|
| 727,285 | HACApK | 8772 MB | **10%** | 256 | 69 | 64.9s | 172.8s | **240s** |

---

## Key Findings

### Performance Summary

1. **HACApK vs LU**: HACApK is **48x faster** for hex N=15 nonlinear (16s vs 767s)
2. **HACApK vs BiCGSTAB**: HACApK is **3.5x faster** for hex N=20 nonlinear (53s vs 185s)
3. **Memory Efficiency**: HACApK uses **13x less memory** at hex N=20 (2.8 GB vs 37.5 GB)
4. **Compression Ratio**: Improves with problem size (49% at N=10, 25% at N=15, 15% at N=20)
5. **Tetra HACApK Optimization (v1.4.1)**: ELF-style face basis caching provides **12-20x speedup** in H-matrix build

### Solver Recommendations

| Element Type | Material | Problem Size | Recommended Solver |
|--------------|----------|--------------|-------------------|
| Hexahedral | Linear | DOF < 6,000 | Any solver |
| Hexahedral | Linear | DOF > 6,000 | **HACApK** (3x faster) |
| Hexahedral | Nonlinear | DOF < 6,000 | BiCGSTAB or HACApK |
| Hexahedral | Nonlinear | DOF > 6,000 | **HACApK** (33x faster than LU) |
| Tetrahedral | Linear | DOF < 7,000 | HACApK or BiCGSTAB |
| Tetrahedral | Linear | DOF > 7,000 | **HACApK** (2x faster than BiCGSTAB) |
| Tetrahedral | Nonlinear | DOF < 7,000 | HACApK or BiCGSTAB |
| Tetrahedral | Nonlinear | DOF > 7,000 | **HACApK** (3x faster than BiCGSTAB) |

**Note (v1.4.1)**: After ELF-style optimization, HACApK is now recommended for tetrahedral elements at all scales.

---

## H-matrix Statistics Summary

**Note**: Compression ratio = H-matrix memory / Dense memory. Lower is better.

### Hexahedral Elements (6 DOF per element)

| N | Elements | DOF | lowrank | dense | max_rank | H-mat [MB] | Dense [MB] | Compression |
|---|----------|-----|--------:|------:|---------:|-----------:|-----------:|------------:|
| 5 | 125 | 750 | 26 | 158 | 48 | 4 | 4 | 97% |
| 10 | 1,000 | 6,000 | 1,166 | 2,024 | 101 | 137 | 275 | **50%** |
| 15 | 3,375 | 20,250 | 8,202 | 9,442 | 85 | 805 | 3,129 | **26%** |
| 20 | 8,000 | 48,000 | 19,216 | 20,610 | 105 | 2,588 | 17,578 | **15%** |
| 25 | 15,625 | 93,750 | 31,450 | 50,400 | 96 | 5,940 | 67,055 | **9%** |

### Tetrahedral Elements (3 DOF per element)

| maxh | Elements | DOF | lowrank | dense | max_rank | H-mat [MB] | Dense [MB] | Compression |
|------|----------|-----|--------:|------:|---------:|-----------:|-----------:|------------:|
| 0.20 | 627 | 1,881 | 554 | 1,544 | 43 | 23 | 27 | 85% |
| 0.15 | 2,211 | 6,633 | 3,160 | 11,226 | 51 | 186 | 336 | **56%** |
| 0.10 | 4,994 | 14,982 | 12,074 | 24,401 | 52 | 623 | 1,712 | **36%** |
| 0.05 | 33,974 | 101,922 | 133,954 | 203,403 | 54 | 8,213 | 79,255 | **10%** |

---

## Computational Complexity

### Theoretical Complexity

| Solver | Time Complexity | Memory Complexity |
|--------|-----------------|-------------------|
| Dense LU | O(N^3) | O(N^2) |
| Dense BiCGSTAB | O(N^2) per iter | O(N^2) |
| BiCGSTAB+H-matrix | **O(N log N)** per iter | **O(N log N)** |

### Measured Scaling Exponents (Nonlinear, last 9 data points)

Power-law fit: `t = a * N^alpha` in log-log space.

#### Hexahedral Elements

| Solver | Time Exponent | Time R² | Memory Exponent | Memory R² |
|--------|--------------|---------|-----------------|-----------|
| LU | 3.48 | 0.996 | 1.91 | 1.000 |
| BiCGSTAB | 2.33 | 0.989 | 1.98 | 1.000 |
| **HACApK** | **1.58** | 0.985 | **1.28** | 0.999 |

#### Tetrahedral Elements

| Solver | Time Exponent | Time R² | Memory Exponent | Memory R² |
|--------|--------------|---------|-----------------|-----------|
| **HACApK** | **1.48** | 0.995 | **1.26** | 0.999 |

**Note**: HACApK memory exponent ~1.3 is consistent with O(N log N). The time exponent ~1.5 is higher than memory due to increasing BiCGSTAB iterations and nonlinear iteration count variation with problem size. See `hexahedron/plot_scaling.py` and `tetrahedron/plot_scaling.py` for details.

---

## H-Matrix and ACA+ Algorithm

### Hierarchical Matrix (H-Matrix)

The H-matrix approach partitions the interaction matrix into a hierarchy of blocks. Blocks representing well-separated clusters are approximated as low-rank matrices, while nearby interactions are stored as dense blocks.

**Cluster Tree Construction**:
1. Build binary space partitioning tree of element centroids
2. Leaf nodes contain at most `leaf_size` elements (default: 10)
3. Block admissibility: cluster pair (i,j) is low-rank if `dist(i,j) >= eta * max(diam(i), diam(j))` (eta=2.0)

### ACA+ (Adaptive Cross Approximation with Pivoting)

ACA+ is used to construct low-rank approximations of admissible blocks. It approximates a matrix A as:

```
A ~ U * V^T,  where U is m x k and V is n x k
```

**ACA+ Algorithm with Partial Pivoting**:

1. **Initialization**: Set residual R = A, rank k = 0
2. **Row Selection (Pivoting)**: Find row index i with maximum residual norm
3. **Column Selection**: In row i, find column j with maximum absolute value
4. **Update**:
   - u_k = R(:, j) / R(i, j)  (column vector)
   - v_k = R(i, :)^T         (row vector)
   - R = R - u_k * v_k^T     (rank-1 update)
5. **Convergence Check**: If ||u_k|| * ||v_k|| < eps * ||A||_F, stop
6. **Increment**: k = k + 1, goto step 2

**Key Features**:
- **Partial Pivoting**: Selects the row with largest residual norm, ensuring numerical stability
- **Adaptive Rank**: Automatically determines the rank needed to achieve tolerance `eps`
- **Efficiency**: O(k * (m + n)) operations per rank-1 update
- **No Full Matrix Required**: Only accesses matrix elements as needed

**Parameters**:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eps` | 1e-4 | ACA+ tolerance (convergence criterion) |
| `leaf_size` | 10 | Maximum elements per leaf cluster |
| `eta` | 2.0 | Admissibility parameter (higher = more low-rank blocks) |

**Typical Compression Results**:

| DOF | Compression | Memory Savings |
|-----|-------------|----------------|
| 6,000 | 50% | 2x |
| 20,000 | 26% | 4x |
| 48,000 | 15% | 7x |
| 94,000 | 9% | 11x |

---

## BH Curve (Nonlinear Material)

```
H [A/m]     B [T]     Notes
0           0.0
100         0.1       Initial mu_r ~ 800
200         0.3
500         0.8
1000        1.2
2000        1.5       Saturation begins
5000        1.7
10000       1.8
50000       2.0       Strong saturation
100000      2.1
```

---

## Usage

### Linear Benchmarks

```bash
cd hexahedron
python benchmark_hex.py --linear --lu --bicgstab --hacapk 10 15

cd tetrahedron
python benchmark_tetra.py --linear --lu --bicgstab --hacapk 0.20 0.15 0.10
```

### Nonlinear Benchmarks

```bash
cd hexahedron
python benchmark_hex.py --nonlinear --lu --bicgstab --hacapk 10 15
python benchmark_hex.py --nonlinear --hacapk 20

cd tetrahedron
python benchmark_tetra.py --nonlinear --lu --bicgstab --hacapk 0.20 0.15 0.10
python benchmark_tetra.py --nonlinear --hacapk 0.05
```

---

## Memory Measurement Notes

### Windows: peak_wset (Peak Working Set)

Memory measurements use Windows `peak_wset` via `psutil.Process().memory_info().peak_wset`.

**Dense matrix reference (for comparison):**
- Dense matrix memory = N^2 x 8 bytes (double precision)
- N=6,000 DOF: 275 MB dense matrix
- N=20,250 DOF: 3,129 MB dense matrix
- N=48,000 DOF: 17,578 MB dense matrix

### H-matrix Compression Ratio

The "Compression" column shows: **H-matrix memory / Dense matrix memory x 100%**

- **< 30%**: Excellent compression, H-matrix very beneficial
- **30-50%**: Good compression
- **50-80%**: Moderate compression
- **> 80%**: Poor compression (problem too small for H-matrix)

---

**Last Updated**: 2026-03-03 (Added measured scaling exponents for hexahedral and tetrahedral benchmarks)
