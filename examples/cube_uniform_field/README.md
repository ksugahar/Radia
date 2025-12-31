# Cube in Uniform Field Benchmark

Radia magnetostatic solver benchmark for a soft iron cube in uniform external field.

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

## Benchmark Results (2025-12-31, 8-thread)

**System**: 8-core CPU, OMP_NUM_THREADS=8, MKL_NUM_THREADS=8

### 1. Hexahedral Linear Benchmark

#### N=5 (125 elements, 750 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 701,878 | LU | 51 MB | - | 0 | 1 | 0.07s | - | 0.40s | **0.48s** |
| 701,878 | BiCGSTAB | 46 MB | - | 0 | 1 | 0.03s | - | 3ms | **0.04s** |
| 701,878 | HACApK | 41 MB | 97% | 11 | 1 | 0.1ms | 0.03s | 3ms | **0.04s** |

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 716,184 | LU | 933 MB | - | 0 | 1 | 1.0s | - | 1.0s | **2.1s** |
| 716,110 | BiCGSTAB | 632 MB | - | 0 | 1 | 1.0s | - | 0.2s | **1.3s** |
| 716,183 | HACApK | 193 MB | **50%** | 11 | 1 | 0.6ms | 1.3s | 0.1s | **1.4s** |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 720,157 | LU | 9924 MB | - | 0 | 1 | 10.3s | - | 18.9s | **29.2s** |
| 720,159 | BiCGSTAB | 6723 MB | - | 0 | 1 | 10.2s | - | 4.7s | **15.0s** |
| 720,169 | HACApK | 911 MB | **26%** | 10 | 1 | 1.9ms | 8.3s | 0.6s | **9.0s** |

#### N=20 (8,000 elements, 48,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| - | LU | - | - | - | - | - | - | - | **(~7min estimated)** |
| 721,938 | BiCGSTAB | 37481 MB | - | 0 | 1 | 57.7s | - | 24.9s | **82.6s** |
| 721,979 | HACApK | 2796 MB | **15%** | 16 | 1 | 6.3ms | 31.7s | 3.4s | **35.2s** |

#### N=25 (15,625 elements, 93,750 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 723,047 | HACApK | 6286 MB | **9%** | 23 | 1 | 14.2ms | 136.3s | 9.3s | **146.7s** |

---

### 2. Hexahedral Nonlinear Benchmark

**All solvers converge** for hexahedral nonlinear problems.

#### N=5 (125 elements, 750 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 702,114 | LU | 51 MB | - | 0 | 6 | 0.07s | - | 0.42s | **0.49s** |
| 702,114 | BiCGSTAB | 46 MB | - | 0 | 3 | 0.03s | - | 0.01s | **0.03s** |
| 702,114 | HACApK | 42 MB | 97% | 21 | 3 | 0.1ms | 0.03s | 3ms | **0.04s** |

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 716,281 | LU | 933 MB | - | 0 | 13 | 1.1s | - | 5.7s | **7.7s** |
| 716,316 | BiCGSTAB | 632 MB | - | 0 | 5 | 1.0s | - | 0.9s | **2.0s** |
| 716,353 | HACApK | 193 MB | **50%** | 31 | 4 | 0.6ms | 1.3s | 0.3s | **1.6s** |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 719,832 | LU | 9926 MB | - | 0 | 35 | 10.7s | - | 587s | **598s** |
| 719,902 | BiCGSTAB | 6724 MB | - | 0 | 27 | 10.5s | - | 24s | **34s** |
| 719,883 | HACApK | 912 MB | **26%** | 85 | 29 | 1.9ms | 8.7s | 5.5s | **14s** |

#### N=20 (8,000 elements, 48,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| - | LU | - | - | - | - | - | - | - | **(~4h estimated)** |
| 721,269 | BiCGSTAB | 37485 MB | - | 0 | 30 | 57.8s | - | 166s | **224s** |
| 721,299 | HACApK | 2799 MB | **15%** | 120 | 30 | 6.5ms | 32.1s | 24.6s | **57s** |

---

### 3. Tetrahedral Linear Benchmark

#### maxh=0.20m (627 elements, 1,881 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 816,133 | LU | 172 MB | - | 0 | 1 | 0.12s | - | 2.5s | 2.7s |
| 816,080 | BiCGSTAB | 140 MB | - | 0 | 1 | 0.12s | - | 0.08s | 0.4s |
| 816,089 | HACApK | 98 MB | 85% | 29 | 1 | 4.0ms | 0.42s | 0.04s | **0.48s** |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 748,693 | LU | 1272 MB | - | 0 | 1 | 1.6s | - | 3.4s | 5.2s |
| 748,683 | BiCGSTAB | 922 MB | - | 0 | 1 | 1.6s | - | 1.0s | 3.0s |
| 748,705 | HACApK | 279 MB | **56%** | 25 | 1 | 1.2ms | 1.71s | 0.33s | **2.06s** |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 754,617 | LU | 6105 MB | - | 0 | 1 | 8.9s | - | 34.6s | 45.6s |
| 754,394 | BiCGSTAB | 4378 MB | - | 0 | 1 | 9.0s | - | 4.6s | 15.4s |
| 754,549 | HACApK | 747 MB | **36%** | 29 | 1 | 2.4ms | 6.41s | 1.29s | **7.77s**

#### maxh=0.05m (33,974 elements, 101,922 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 737,397 | HACApK | 8715 MB | **10%** | 26 | 1 | 24ms | 2158s | 14.6s | 2175s |

---

### 4. Tetrahedral Nonlinear Benchmark

**All solvers converge** for tetrahedral nonlinear problems.

#### maxh=0.20m (627 elements, 1,881 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 748,822 | LU | 171 MB | - | 0 | 20 | 0.11s | - | 1.8s | 2.0s |
| 748,869 | BiCGSTAB | 140 MB | - | 0 | 14 | 0.12s | - | 0.24s | 0.4s |
| 748,941 | HACApK | 98 MB | 85% | 130 | 17 | 0.3ms | 0.32s | 0.15s | **0.51s** |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 730,715 | LU | 1273 MB | - | 0 | 29 | 1.6s | - | 87.4s | 91.5s |
| 730,760 | BiCGSTAB | 922 MB | - | 0 | 21 | 1.7s | - | 4.0s | 6.1s |
| 730,596 | HACApK | 280 MB | **56%** | 134 | 36 | 1.0ms | 1.70s | 2.0s | **3.72s** |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 730,996 | LU | 6108 MB | - | 0 | 39 | 9.0s | - | 1315s | 1336s |
| 731,003 | BiCGSTAB | 4380 MB | - | 0 | 41 | 9.3s | - | 36.3s | 47.7s |
| 731,063 | HACApK | 751 MB | **36%** | 206 | 41 | 2.6ms | 6.42s | 9.9s | **16.7s**

#### maxh=0.05m (33,974 elements, 101,922 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 727,285 | HACApK | 8734 MB | **10%** | 256 | 69 | 23ms | 2204s | 152s | 2357s |

---

## Key Findings

### Performance Summary

1. **HACApK vs LU**: HACApK is **43x faster** for hex N=15 nonlinear (14s vs 598s)
2. **HACApK vs BiCGSTAB**: HACApK is **3.9x faster** for hex N=20 nonlinear (57s vs 224s)
3. **Memory Efficiency**: HACApK uses **13x less memory** at hex N=20 (2.8 GB vs 37.5 GB)
4. **Compression Ratio**: Improves with problem size (50% at N=10, 26% at N=15, 15% at N=20)
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

| Solver | Time Complexity | Memory Complexity |
|--------|-----------------|-------------------|
| Dense LU | O(N^3) | O(N^2) |
| Dense BiCGSTAB | O(N^2) per iter | O(N^2) |
| BiCGSTAB+H-matrix | **O(N log N)** per iter | **O(N log N)** |

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

**Last Updated**: 2025-12-31 (8-thread verified benchmark: N=5, N=10, N=15, N=20)
