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

### Hexahedral Benchmark (Linear)

#### N=5 (125 elements, 750 DOF)

| M_avg_z | Solver | Memory | Iter | Time |
|--------:|--------|-------:|-----:|-----:|
| 701,900 | Dense LU | 34 MB | 2 | 1.38s |
| 701,907 | Dense BiCGSTAB | 33 MB | 4 | 0.15s |
| 701,897 | H-matrix | 33 MB | 3 | 0.08s |

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Iter | Time |
|--------:|--------|-------:|-----:|-----:|
| 716,184 | Dense LU | 60 MB | 2 | 9.6s |
| 716,184 | Dense BiCGSTAB | 37 MB | 3 | 9.7s |
| 716,192 | H-matrix | 41 MB | 2 | 4.1s |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Iter | Time |
|--------:|--------|-------:|-----:|-----:|
| 720,157 | Dense LU | 123 MB | 2 | 147s |
| 720,156 | Dense BiCGSTAB | 52 MB | 3 | 133s |
| 720,186 | H-matrix | 58 MB | 2 | **28s** |

#### N=20 (8,000 elements, 48,000 DOF)

| M_avg_z | Solver | Memory | Iter | Time |
|--------:|--------|-------:|-----:|-----:|
| - | Dense LU | OOM | - | - |
| - | Dense BiCGSTAB | OOM | - | - |
| 722,006 | H-matrix | 92 MB | 3 | **95s** |

### Tetrahedral Benchmark (Linear, Netgen mesh)

#### maxh=0.35m (112 elements, 336 DOF)

| M_avg_z | Solver | Memory | Iter | Time |
|--------:|--------|-------:|-----:|-----:|
| 908,535 | Dense LU | 64 MB | 2 | 0.78s |
| 908,291 | Dense BiCGSTAB | 65 MB | 2 | 0.19s |
| 908,290 | H-matrix | 65 MB | 2 | 0.17s |

#### maxh=0.25m (390 elements, 1,170 DOF)

| M_avg_z | Solver | Memory | Iter | Time |
|--------:|--------|-------:|-----:|-----:|
| 765,437 | Dense LU | 69 MB | 2 | 2.2s |
| 765,448 | Dense BiCGSTAB | 66 MB | 2 | 2.1s |
| 765,458 | H-matrix | 68 MB | 2 | 1.9s |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Iter | Time |
|--------:|--------|-------:|-----:|-----:|
| 748,693 | Dense LU | 90 MB | 2 | 79s |
| 748,683 | Dense BiCGSTAB | 75 MB | 2 | 76s |
| 748,718 | H-matrix | 79 MB | 2 | 79s |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Iter | Time |
|--------:|--------|-------:|-----:|-----:|
| 754,617 | Dense LU | 104 MB | 2 | 561s |
| 754,394 | Dense BiCGSTAB | 88 MB | 2 | 530s |
| 754,533 | H-matrix | 95 MB | 2 | 499s |

---

## Nonlinear Benchmark Results (2025-12-26)

### Hexahedral Benchmark (Nonlinear)

#### N=5 (125 elements, 750 DOF)

| M_avg_z | Solver | Memory | Nonl | Time |
|--------:|--------|-------:|-----:|-----:|
| 702,132 | Dense LU | 34 MB | 6 | 0.69s |
| 702,114 | Dense BiCGSTAB | 32 MB | 3 | 0.23s |
| 701,751 | H-matrix | 33 MB | 4 | 0.26s |

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Nonl | Time |
|--------:|--------|-------:|-----:|-----:|
| 716,281 | Dense LU | 61 MB | 13 | 18.5s |
| 716,318 | Dense BiCGSTAB | 37 MB | 5 | 13.2s |
| 715,883 | H-matrix | 39 MB | 5 | 5.0s |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Nonl | Time |
|--------:|--------|-------:|-----:|-----:|
| 719,832 | Dense LU | 124 MB | 35 | 1021s |
| 719,875 | Dense BiCGSTAB | 53 MB | 29 | 348s |
| 719,399 | H-matrix | 59 MB | 26 | **36s** |

#### N=20 (8,000 elements, 48,000 DOF)

| M_avg_z | Solver | Memory | Nonl | Time |
|--------:|--------|-------:|-----:|-----:|
| - | Dense LU | OOM | - | - |
| - | Dense BiCGSTAB | OOM | - | - |
| 720,829 | H-matrix | 95 MB | 25 | **121s** |

### Tetrahedral Benchmark (Nonlinear, Netgen mesh)

**Note**: Only LU converges. BiCGSTAB and HACApK do not converge within 100 iterations.

#### maxh=0.30m (200 elements, 600 DOF)

| M_avg_z | Solver | Memory | Nonl | Time | Converged |
|--------:|--------|-------:|-----:|-----:|-----------|
| 726,577 | Dense LU | 66 MB | 14 | 0.87s | Yes |
| 639,258 | Dense BiCGSTAB | 66 MB | 100 | 1.18s | No |
| 715,943 | H-matrix | 67 MB | 100 | 1.30s | No |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Nonl | Time | Converged |
|--------:|--------|-------:|-----:|-----:|-----------|
| 730,715 | Dense LU | 92 MB | 29 | 181s | Yes |
| 710,935 | Dense BiCGSTAB | 76 MB | 100 | 415s | No |
| 706,867 | H-matrix | 82 MB | 100 | 119s | No |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Nonl | Time | Converged |
|--------:|--------|-------:|-----:|-----:|-----------|
| 730,996 | Dense LU | 104 MB | 39 | 1813s | Yes |
| - | Dense BiCGSTAB | - | 100 | - | No |
| 691,313 | H-matrix | 95 MB | 100 | 421s | No |

---

## Key Findings

### Linear Material

1. **All solvers converge**: LU, BiCGSTAB, HACApK all work for both element types
2. **H-matrix speedup**: 5x faster than LU at N=15 for hexahedra
3. **N=20 requires H-matrix**: Dense solvers run out of memory

### Hexahedral Elements (Nonlinear - Good Convergence)

1. **All solvers converge**: LU, BiCGSTAB, HACApK all work
2. **H-matrix is fastest**: 28x faster than LU at N=15, 10x faster than BiCGSTAB
3. **N=20 requires H-matrix**: Dense solvers run out of memory

### Tetrahedral Elements (Nonlinear - Convergence Issues)

1. **Only LU converges**: BiCGSTAB and HACApK fail to converge
2. **Different behavior from ELF**: ELF converges with all solvers
3. **Under investigation**: Radia tetrahedral nonlinear solver needs debugging

### Solver Recommendations

| Element Type | Material | Problem Size | Recommended Solver |
|--------------|----------|--------------|-------------------|
| Hexahedral | Linear/Nonlinear | DOF < 6,000 | BiCGSTAB |
| Hexahedral | Linear/Nonlinear | DOF 6,000-20,000 | HACApK |
| Hexahedral | Linear/Nonlinear | DOF > 20,000 | **HACApK required** |
| Tetrahedral | Linear | All sizes | Any solver |
| Tetrahedral | Nonlinear | All sizes | **LU only** (temporary) |

---

## Computational Complexity

| Solver | Time | Memory |
|--------|------|--------|
| Dense LU | O(N^3) | O(N^2) |
| Dense BiCGSTAB | O(N^2)/iter | O(N^2) |
| BiCGSTAB+H-matrix | **O(N log N)** | **O(N log N)** |

---

## H-matrix Statistics (Hexahedral)

| N | DOF | lowrank | dense | max_rank | leaves |
|---|-----|---------|-------|----------|--------|
| 5 | 750 | 26 | 158 | 46 | 184 |
| 10 | 6,000 | 1,166 | 2,024 | 90 | 3,190 |
| 15 | 20,250 | 8,262 | 9,124 | 79 | 17,386 |
| 20 | 48,000 | 17,976 | 20,656 | 99 | 38,632 |

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
