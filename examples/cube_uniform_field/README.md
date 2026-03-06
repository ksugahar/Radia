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

## Benchmark Results (2025-12-30)

### 1. Hexahedral Linear Benchmark

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 716,184 | LU | 925 MB | - | 0 | 1 | 1.9s | - | 1.1s | 3.1s |
| 716,110 | BiCGSTAB | 628 MB | - | 0 | 1 | 1.9s | - | 0.2s | 2.2s |
| 716,183 | HACApK | 187 MB | **50%** | 11 | 1 | 0.7ms | 2.1s | 0.1s | 2.3s |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 720,157 | LU | 9917 MB | - | 0 | 1 | 20.5s | - | 24.0s | 46.0s |
| 720,159 | BiCGSTAB | 6720 MB | - | 0 | 1 | 20.4s | - | 2.6s | 23.8s |
| 720,169 | HACApK | 906 MB | **26%** | 10 | 1 | 5.6ms | 14.9s | 0.6s | 15.6s |

#### N=20 (8,000 elements, 48,000 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 721,979 | HACApK | 2788 MB | **15%** | 16 | 1 | 8.5ms | 55.3s | 2.9s | 58.6s |

---

### 2. Hexahedral Nonlinear Benchmark

**All solvers converge** for hexahedral nonlinear problems.

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 716,281 | LU | 925 MB | - | 0 | 13 | 1.9s | - | 8.8s | 11.6s |
| 716,316 | BiCGSTAB | 629 MB | - | 0 | 5 | 1.9s | - | 1.2s | 3.1s |
| 716,353 | HACApK | 188 MB | **50%** | 31 | 4 | 0.7ms | 2.1s | 0.3s | 2.4s |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 719,832 | LU | 9920 MB | - | 0 | 35 | 21.0s | - | 876s | 927s |
| 719,838 | BiCGSTAB | 6721 MB | - | 0 | 34 | 25.5s | - | 36.8s | 63.3s |
| 719,873 | HACApK | 906 MB | **26%** | 91 | 31 | 3.5ms | 18.6s | 9.4s | 28.1s |

#### N=20 (8,000 elements, 48,000 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 721,305 | HACApK | 2792 MB | **15%** | 118 | 30 | 8.8ms | 55.0s | 23.2s | 78.6s |

---

### 3. Tetrahedral Linear Benchmark

#### maxh=0.20m (627 elements, 1,881 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 816,133 | LU | 171 MB | - | 0 | 1 | 0.15s | - | 0.26s | 0.45s |
| 816,080 | BiCGSTAB | 140 MB | - | 0 | 1 | 0.15s | - | 0.07s | 0.23s |
| 816,089 | HACApK | 98 MB | 85% | 29 | 1 | 0.3ms | 3.9s | 0.04s | 3.9s |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 748,693 | LU | 1272 MB | - | 0 | 1 | 1.9s | - | 3.1s | 5.3s |
| 748,683 | BiCGSTAB | 922 MB | - | 0 | 1 | 1.9s | - | 1.0s | 3.0s |
| 748,705 | HACApK | 279 MB | **56%** | 25 | 1 | 0.9ms | 17.4s | 0.3s | 17.8s |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 754,617 | LU | 6105 MB | - | 0 | 1 | 9.4s | - | 34.0s | 46.4s |
| 754,394 | BiCGSTAB | 4378 MB | - | 0 | 1 | 9.4s | - | 4.6s | 16.3s |
| 754,549 | HACApK | 748 MB | **36%** | 29 | 1 | 6.3ms | 66.9s | 1.3s | 68.3s |

#### maxh=0.05m (33,974 elements, 101,922 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 737,397 | HACApK | 8716 MB | **10%** | 26 | 1 | 24ms | 1355s | 14.6s | 1372s |

---

### 4. Tetrahedral Nonlinear Benchmark

**All solvers converge** for tetrahedral nonlinear problems.

#### maxh=0.20m (627 elements, 1,881 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 748,822 | LU | 172 MB | - | 0 | 20 | 0.15s | - | 1.8s | 2.1s |
| 748,869 | BiCGSTAB | 140 MB | - | 0 | 14 | 0.16s | - | 0.26s | 0.42s |
| 748,941 | HACApK | 98 MB | 85% | 130 | 17 | 0.3ms | 3.8s | 0.18s | 4.0s |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 730,715 | LU | 1273 MB | - | 0 | 29 | 1.9s | - | 87.9s | 92.3s |
| 730,760 | BiCGSTAB | 922 MB | - | 0 | 21 | 1.9s | - | 4.1s | 6.2s |
| 730,596 | HACApK | 280 MB | **56%** | 134 | 36 | 1.0ms | 23.6s | 2.0s | 25.6s |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 730,996 | LU | 6107 MB | - | 0 | 39 | 9.4s | - | 1341s | 1370s |
| 731,005 | BiCGSTAB | 4380 MB | - | 0 | 45 | 9.7s | - | 40.5s | 52.7s |
| 731,063 | HACApK | 750 MB | **36%** | 206 | 41 | 2.6ms | 68.9s | 11.2s | 80.3s |

#### maxh=0.05m (33,974 elements, 101,922 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 727,285 | HACApK | 8730 MB | **10%** | 256 | 69 | 23ms | 1452s | 152s | 1606s |

---

## Key Findings

### Performance Summary

1. **HACApK vs LU**: HACApK is **33x faster** for hex N=15 nonlinear (28s vs 927s)
2. **HACApK vs BiCGSTAB**: HACApK is **2.3x faster** for hex N=15 nonlinear (28s vs 63s)
3. **Memory Efficiency**: HACApK uses **11x less memory** at hex N=15 (906 MB vs 9920 MB)
4. **Compression Ratio**: Improves with problem size (50% at N=10, 26% at N=15, 15% at N=20)

### Solver Recommendations

| Element Type | Material | Problem Size | Recommended Solver |
|--------------|----------|--------------|-------------------|
| Hexahedral | Linear | DOF < 6,000 | Any solver |
| Hexahedral | Linear | DOF > 6,000 | **HACApK** (3x faster) |
| Hexahedral | Nonlinear | DOF < 6,000 | BiCGSTAB or HACApK |
| Hexahedral | Nonlinear | DOF > 6,000 | **HACApK** (33x faster than LU) |
| Tetrahedral | Linear | DOF < 7,000 | BiCGSTAB (fastest) |
| Tetrahedral | Linear | DOF > 7,000 | BiCGSTAB (H-matrix build overhead) |
| Tetrahedral | Nonlinear | DOF < 7,000 | BiCGSTAB (fastest) |
| Tetrahedral | Nonlinear | DOF > 7,000 | BiCGSTAB or HACApK |

**Note**: For tetrahedral elements, H-matrix build time dominates at smaller scales. BiCGSTAB is often faster for DOF < 15,000.

---

## H-matrix Statistics Summary

**Note**: Compression ratio = H-matrix memory / Dense memory. Lower is better.

### Hexahedral Elements (6 DOF per element)

| N | Elements | DOF | lowrank | dense | max_rank | H-mat [MB] | Dense [MB] | Compression |
|---|----------|-----|--------:|------:|---------:|-----------:|-----------:|------------:|
| 10 | 1,000 | 6,000 | 1,166 | 2,024 | 101 | 137 | 275 | **50%** |
| 15 | 3,375 | 20,250 | 8,202 | 9,442 | 85 | 805 | 3,129 | **26%** |
| 20 | 8,000 | 48,000 | 19,216 | 20,610 | 105 | 2,588 | 17,578 | **15%** |

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

**Last Updated**: 2025-12-30 (All benchmark results updated from JSON files)
