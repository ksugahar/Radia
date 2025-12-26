# Cube in Uniform Field Benchmark

Radia magnetostatic solver benchmark for a soft iron cube in uniform external field.

## Folder Structure

```
cube_uniform_field/
├── linear/                     # Linear material (constant mu_r)
│   ├── benchmark_hexahedron.py # Hexahedral benchmark script
│   ├── benchmark_tetrahedron.py# Tetrahedral benchmark script
│   ├── hexahedron/             # Hexahedral results
│   │   ├── lu/
│   │   ├── bicgstab/
│   │   └── hacapk/
│   └── tetrahedron/            # Tetrahedral results
│       ├── lu/
│       ├── bicgstab/
│       └── hacapk/
├── nonlinear/                  # Nonlinear material (saturation BH curve)
│   ├── benchmark_common.py     # Shared benchmark functions
│   ├── generate_readme.py      # README generation script
│   ├── hexahedron/             # Hexahedral benchmark
│   │   ├── benchmark_hex.py
│   │   ├── lu/
│   │   ├── bicgstab/
│   │   └── hacapk/
│   └── tetrahedron/            # Tetrahedral benchmark
│       ├── benchmark_tetra.py
│       ├── lu/
│       ├── bicgstab/
│       └── hacapk/
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

### Analytical Solution (Linear)

M_analytical = chi * H_ext / (1 + chi * N_demag) = 598,204 A/m

Where N_demag = 1/3 for a cube.

---

## Linear Benchmark Results (2025-12-26)

### Hexahedral Benchmark (Linear) - Detailed Analysis

#### N=5 (125 elements, 750 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Iter | H-mat Build [s] | Total [s] |
|--------|-------------:|------------:|-----:|----------------:|----------:|
| Dense LU | 701,900 | 34 | 2 | - | 1.38 |
| Dense BiCGSTAB | 701,907 | 33 | 4 | - | 0.15 |
| H-matrix | 701,897 | 33 | 3 | 0.07 | 0.08 |

**H-matrix stats (N=5)**: lowrank=26, dense=158, max_rank=46, leaves=184

#### N=10 (1,000 elements, 6,000 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Iter | H-mat Build [s] | Total [s] |
|--------|-------------:|------------:|-----:|----------------:|----------:|
| Dense LU | 716,184 | 60 | 2 | - | 9.6 |
| Dense BiCGSTAB | 716,184 | 37 | 3 | - | 9.7 |
| H-matrix | 716,192 | 41 | 2 | 3.97 | 4.1 |

**H-matrix stats (N=10)**: lowrank=1,166, dense=2,024, max_rank=90, leaves=3,190

#### N=15 (3,375 elements, 20,250 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Iter | H-mat Build [s] | Total [s] |
|--------|-------------:|------------:|-----:|----------------:|----------:|
| Dense LU | 720,157 | 123 | 2 | - | 147 |
| Dense BiCGSTAB | 720,156 | 52 | 3 | - | 133 |
| H-matrix | 720,186 | 58 | 2 | 26.5 | **28** |

**H-matrix stats (N=15)**: lowrank=8,262, dense=9,124, max_rank=79, leaves=17,386, compression=0.32%

#### N=20 (8,000 elements, 48,000 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Iter | H-mat Build [s] | Total [s] |
|--------|-------------:|------------:|-----:|----------------:|----------:|
| Dense LU | - | OOM | - | - | - |
| Dense BiCGSTAB | - | OOM | - | - | - |
| H-matrix | 722,006 | 92 | 3 | 87.2 | **95** |

**H-matrix stats (N=20)**: lowrank=17,976, dense=20,656, max_rank=99, leaves=38,632

---

### Tetrahedral Benchmark (Linear, Netgen mesh) - Detailed Analysis

#### maxh=0.35m (112 elements, 336 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Iter | H-mat Build [s] | Total [s] |
|--------|-------------:|------------:|-----:|----------------:|----------:|
| Dense LU | 908,535 | 64 | 2 | - | 0.78 |
| Dense BiCGSTAB | 908,291 | 65 | 2 | - | 0.19 |
| H-matrix | 908,290 | 65 | 2 | 0.01 | 0.17 |

**H-matrix stats**: lowrank=2, dense=254, max_rank=19, leaves=256

#### maxh=0.25m (390 elements, 1,170 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Iter | H-mat Build [s] | Total [s] |
|--------|-------------:|------------:|-----:|----------------:|----------:|
| Dense LU | 765,437 | 69 | 2 | - | 2.2 |
| Dense BiCGSTAB | 765,448 | 66 | 2 | - | 2.1 |
| H-matrix | 765,458 | 68 | 2 | 0.01 | 1.9 |

**H-matrix stats**: lowrank=290, dense=869, max_rank=36, leaves=1,159

#### maxh=0.20m (627 elements, 1,881 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Iter | H-mat Build [s] | Total [s] |
|--------|-------------:|------------:|-----:|----------------:|----------:|
| Dense LU | 816,133 | 71 | 2 | - | 5.8 |
| Dense BiCGSTAB | 816,056 | 67 | 3 | - | 5.6 |
| H-matrix | 816,165 | 70 | 3 | 0.03 | 5.0 |

**H-matrix stats**: lowrank=554, dense=1,544, max_rank=41, leaves=2,098

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Iter | H-mat Build [s] | Total [s] |
|--------|-------------:|------------:|-----:|----------------:|----------:|
| Dense LU | 748,693 | 90 | 2 | - | 79 |
| Dense BiCGSTAB | 748,683 | 75 | 2 | - | 76 |
| H-matrix | 748,718 | 79 | 2 | 0.41 | 79 |

**H-matrix stats**: lowrank=3,160, dense=11,226, max_rank=46, leaves=14,386

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Iter | H-mat Build [s] | Total [s] |
|--------|-------------:|------------:|-----:|----------------:|----------:|
| Dense LU | 754,617 | 104 | 2 | - | 561 |
| Dense BiCGSTAB | 754,394 | 88 | 2 | - | 530 |
| H-matrix | 754,533 | 95 | 2 | 1.17 | 499 |

**H-matrix stats**: lowrank=12,074, dense=24,401, max_rank=52, leaves=36,475

---

## Nonlinear Benchmark Results (2025-12-26)

### Hexahedral Benchmark (Nonlinear) - Detailed Analysis

#### N=5 (125 elements, 750 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Nonl Iter | H-mat Build [s] | Total [s] | Converged |
|--------|-------------:|------------:|----------:|----------------:|----------:|-----------|
| Dense LU | 702,132 | 34 | 6 | - | 0.69 | Yes |
| Dense BiCGSTAB | 702,114 | 32 | 3 | - | 0.23 | Yes |
| H-matrix | 701,751 | 33 | 4 | 0.07 | 0.26 | Yes |

**H-matrix stats (N=5)**: lowrank=26, dense=158, max_rank=46, leaves=184

#### N=10 (1,000 elements, 6,000 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Nonl Iter | H-mat Build [s] | Total [s] | Converged |
|--------|-------------:|------------:|----------:|----------------:|----------:|-----------|
| Dense LU | 716,281 | 61 | 13 | - | 18.5 | Yes |
| Dense BiCGSTAB | 716,318 | 37 | 5 | - | 13.2 | Yes |
| H-matrix | 715,883 | 39 | 5 | 3.93 | 5.0 | Yes |

**H-matrix stats (N=10)**: lowrank=1,166, dense=2,024, max_rank=90, leaves=3,190

#### N=15 (3,375 elements, 20,250 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Nonl Iter | H-mat Build [s] | Total [s] | Converged |
|--------|-------------:|------------:|----------:|----------------:|----------:|-----------|
| Dense LU | 719,832 | 124 | 35 | - | 1021 | Yes |
| Dense BiCGSTAB | 719,875 | 53 | 29 | - | 348 | Yes |
| H-matrix | 719,399 | 59 | 26 | 26.3 | **36** | Yes |

**H-matrix stats (N=15)**: lowrank=8,262, dense=9,124, max_rank=79, leaves=17,386, compression=0.32%

#### N=20 (8,000 elements, 48,000 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Nonl Iter | H-mat Build [s] | Total [s] | Converged |
|--------|-------------:|------------:|----------:|----------------:|----------:|-----------|
| Dense LU | - | OOM | - | - | - | - |
| Dense BiCGSTAB | - | OOM | - | - | - | - |
| H-matrix | 720,829 | 95 | 25 | 87.0 | **121** | Yes |

**H-matrix stats (N=20)**: lowrank=17,976, dense=20,656, max_rank=99, leaves=38,632

---

### Tetrahedral Benchmark (Nonlinear, Netgen mesh) - Detailed Analysis

**Note**: Only LU converges for nonlinear tetrahedral problems. BiCGSTAB and HACApK do not converge within 100 iterations.

#### maxh=0.35m (176 elements, 528 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Nonl Iter | H-mat Build [s] | Total [s] | Converged |
|--------|-------------:|------------:|----------:|----------------:|----------:|-----------|
| Dense LU | 724,820 | 65 | 12 | - | 0.70 | Yes |
| Dense BiCGSTAB | 649,818 | 65 | 100 | - | 1.00 | No |
| H-matrix | 717,152 | 67 | 100 | 0.01 | 1.21 | No |

#### maxh=0.30m (200 elements, 600 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Nonl Iter | H-mat Build [s] | Total [s] | Converged |
|--------|-------------:|------------:|----------:|----------------:|----------:|-----------|
| Dense LU | 726,577 | 66 | 14 | - | 0.87 | Yes |
| Dense BiCGSTAB | 639,258 | 66 | 100 | - | 1.18 | No |
| H-matrix | 715,943 | 67 | 100 | 0.01 | 1.30 | No |

#### maxh=0.25m (389 elements, 1,167 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Nonl Iter | H-mat Build [s] | Total [s] | Converged |
|--------|-------------:|------------:|----------:|----------------:|----------:|-----------|
| Dense LU | 726,958 | 68 | 17 | - | 2.25 | Yes |
| Dense BiCGSTAB | 655,929 | 66 | 100 | - | 5.95 | No |
| H-matrix | 714,299 | 69 | 100 | 0.01 | 3.24 | No |

#### maxh=0.20m (627 elements, 1,881 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Nonl Iter | H-mat Build [s] | Total [s] | Converged |
|--------|-------------:|------------:|----------:|----------------:|----------:|-----------|
| Dense LU | 729,107 | 71 | 21 | - | 6.60 | Yes |
| Dense BiCGSTAB | 674,866 | 68 | 100 | - | 25.9 | No |
| H-matrix | 714,631 | 71 | 100 | 0.03 | 9.80 | No |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Nonl Iter | H-mat Build [s] | Total [s] | Converged |
|--------|-------------:|------------:|----------:|----------------:|----------:|-----------|
| Dense LU | 730,715 | 92 | 29 | - | 181 | Yes |
| Dense BiCGSTAB | 710,935 | 76 | 100 | - | 415 | No |
| H-matrix | 706,867 | 82 | 100 | 0.43 | 119 | No |

**H-matrix stats**: lowrank=3,160, dense=11,226, max_rank=46, leaves=14,386

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| Solver | M_avg_z [A/m] | Memory [MB] | Nonl Iter | H-mat Build [s] | Total [s] | Converged |
|--------|-------------:|------------:|----------:|----------------:|----------:|-----------|
| Dense LU | 730,996 | 104 | 39 | - | 1813 | Yes |
| Dense BiCGSTAB | - | - | 100 | - | - | No |
| H-matrix | 691,313 | 95 | 100 | 1.50 | 421 | No |

**H-matrix stats**: lowrank=12,074, dense=24,401, max_rank=52, leaves=36,475

---

## Key Findings

### Linear Material

1. **All solvers converge**: LU, BiCGSTAB, HACApK all work for both element types
2. **H-matrix speedup at N=15**: 5x faster than LU (28s vs 147s)
3. **N=20 requires H-matrix**: Dense solvers run out of memory (OOM)
4. **H-matrix build time dominates**: At N=15, H-matrix build is 26.5s of 28s total

### Hexahedral Elements (Nonlinear - Good Convergence)

1. **All solvers converge**: LU, BiCGSTAB, HACApK all work
2. **H-matrix is fastest at N=15**: 28x faster than LU (36s vs 1021s)
3. **H-matrix is 10x faster than BiCGSTAB at N=15**: 36s vs 348s
4. **N=20 requires H-matrix**: Dense solvers run out of memory

### Tetrahedral Elements (Nonlinear - Convergence Issues)

1. **Only LU converges**: BiCGSTAB and HACApK fail to converge within 100 iterations
2. **Different behavior from ELF**: ELF converges with all solvers for tetrahedra
3. **Under investigation**: Radia tetrahedral nonlinear solver needs debugging
4. **Residual stays high**: BiCGSTAB/HACApK residual ~0.5-0.6 at 100 iterations

### Solver Recommendations

| Element Type | Material | Problem Size | Recommended Solver |
|--------------|----------|--------------|-------------------|
| Hexahedral | Linear/Nonlinear | DOF < 6,000 | BiCGSTAB |
| Hexahedral | Linear/Nonlinear | DOF 6,000-20,000 | HACApK |
| Hexahedral | Linear/Nonlinear | DOF > 20,000 | **HACApK required** |
| Tetrahedral | Linear | All sizes | Any solver |
| Tetrahedral | Nonlinear | All sizes | **LU only** (temporary) |

---

## H-matrix Statistics Summary

### Hexahedral Elements (6 DOF per element)

| N | Elements | DOF | lowrank | dense | max_rank | leaves | build [s] | compression |
|---|----------|-----|--------:|------:|---------:|-------:|----------:|------------:|
| 5 | 125 | 750 | 26 | 158 | 46 | 184 | 0.07 | - |
| 10 | 1,000 | 6,000 | 1,166 | 2,024 | 90 | 3,190 | 3.97 | - |
| 15 | 3,375 | 20,250 | 8,262 | 9,124 | 79 | 17,386 | 26.5 | 0.32% |
| 20 | 8,000 | 48,000 | 17,976 | 20,656 | 99 | 38,632 | 87.2 | - |

### Tetrahedral Elements (3 DOF per element)

| maxh | Elements | DOF | lowrank | dense | max_rank | leaves | build [s] |
|------|----------|-----|--------:|------:|---------:|-------:|----------:|
| 0.35 | 112 | 336 | 2 | 254 | 19 | 256 | 0.01 |
| 0.25 | 390 | 1,170 | 290 | 869 | 36 | 1,159 | 0.01 |
| 0.20 | 627 | 1,881 | 554 | 1,544 | 41 | 2,098 | 0.03 |
| 0.15 | 2,211 | 6,633 | 3,160 | 11,226 | 46 | 14,386 | 0.41 |
| 0.10 | 4,994 | 14,982 | 12,074 | 24,401 | 52 | 36,475 | 1.17 |

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
cd linear
python benchmark_hexahedron.py --lu --bicgstab --hacapk 5 10 15
python benchmark_hexahedron.py --hacapk 20
python benchmark_tetrahedron.py --lu --bicgstab --hacapk 0.35 0.25 0.20 0.15 0.10
```

### Nonlinear Benchmarks

```bash
cd nonlinear/hexahedron
python benchmark_hex.py --lu 5 10 15
python benchmark_hex.py --bicgstab 5 10 15
python benchmark_hex.py --hacapk 5 10 15 20

cd nonlinear/tetrahedron
python benchmark_tetra.py --lu 0.30 0.25 0.20 0.15 0.10
python benchmark_tetra.py --bicgstab 0.30 0.25 0.20 0.15 0.10
python benchmark_tetra.py --hacapk 0.30 0.25 0.20 0.15 0.10
```

---

**Last Updated**: 2025-12-26
