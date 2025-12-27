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

## Benchmark Results

### 1. Hexahedral Linear Benchmark (2025-12-27)

**M_analytical = 598,204 A/m** (mu_r = 1000, N_demag = 1/3)

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 716,184 | LU | 61 MB | - | 0 | 2 | - | - | - | - | 10.2s |
| 716,200 | BiCGSTAB | 38 MB | - | 16 | 3 | - | - | - | - | 10.2s |
| 716,192 | HACApK | 41 MB | **47%** | 12 | 2 | 4.01s | - | 4.01s | - | 4.28s |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 720,157 | LU | 123 MB | - | 0 | 2 | - | - | - | - | 155.8s |
| 720,159 | BiCGSTAB | 52 MB | - | 11 | 2 | - | - | - | - | 132.6s |
| 720,186 | HACApK | 57 MB | **24%** | 11 | 2 | 26.0s | - | 26.0s | - | 27.9s |

---

### 2. Hexahedral Nonlinear Benchmark (2025-12-27)

**All solvers converge** for hexahedral nonlinear problems.

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 716,281 | LU | 60 MB | - | 0 | 13 | 0 | 9.19s | - | 9.19s | 18.1s |
| 716,318 | BiCGSTAB | 38 MB | - | 41 | 5 | 0 | - | - | 4.83s | 13.2s |
| 716,362 | HACApK | 211 MB | **48%** | 31 | 4 | 4.02s | - | 4.02s | 0.57s | 4.59s |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 719,832 | LU | 124 MB | - | 0 | 35 | 0 | 991s | - | 991s | 1098s |
| 719,875 | BiCGSTAB | 53 MB | - | 87 | 29 | 0 | - | - | 279s | 372s |
| 719,916 | HACApK | 1246 MB | **25%** | 84 | 28 | 26.1s | - | 26.1s | 12.1s | 38.3s |

#### N=20 (8,000 elements, 48,000 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 721,323 | HACApK | 4798 MB | **14%** | 110 | 30 | 87.8s | - | 87.8s | 40.7s | 128s |

---

### 3. Tetrahedral Linear Benchmark (2025-12-27)

#### maxh=0.30m (200 elements, 600 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 751,153 | LU | 66 MB | - | 0 | 2 | 0 | 0.219s | - | 0.219s | 1.88s |
| 751,162 | BiCGSTAB | 66 MB | - | 45 | 3 | 0 | - | - | 0.172s | 0.81s |
| 751,189 | HACApK | 66 MB | 0.85% | 30 | 3 | 0.007s | - | 0.007s | 0.002s | 0.95s |

#### maxh=0.25m (390 elements, 1,170 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 765,437 | LU | 69 MB | - | 0 | 2 | 0 | 0.257s | - | 0.257s | 2.63s |
| 765,444 | BiCGSTAB | 66 MB | - | 29 | 4 | 0 | - | - | 0.223s | 2.64s |
| 765,454 | HACApK | 68 MB | 1.59% | 29 | 4 | 0.158s | - | 0.158s | 0.014s | 3.20s |

#### maxh=0.20m (627 elements, 1,881 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 816,133 | LU | 72 MB | - | 0 | 2 | 0 | 0.388s | - | 0.388s | 6.54s |
| 816,056 | BiCGSTAB | 67 MB | - | 31 | 3 | 0 | - | - | 0.392s | 6.61s |
| 816,165 | HACApK | 70 MB | 0% | 31 | 3 | 0.175s | - | 0.175s | 0.047s | 4.65s |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 748,693 | LU | 90 MB | - | 0 | 2 | 0 | 6.78s | - | 6.78s | 87.0s |
| 748,677 | BiCGSTAB | 76 MB | - | 31 | 6 | 0 | - | - | 6.28s | 86.7s |
| 748,712 | HACApK | 82 MB | 0% | 31 | 6 | 0.397s | - | 0.397s | 0.476s | 57.5s |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 754,617 | LU | 103 MB | - | 0 | 2 | 0 | 73.2s | - | 73.2s | 674s |
| 754,466 | BiCGSTAB | 88 MB | - | 37 | 6 | 0 | - | - | 154.5s | 562s |
| 754,618 | HACApK | 93 MB | 0% | 36 | 5 | 1.19s | - | 1.19s | 1.73s | 557s |

---

### 4. Tetrahedral Nonlinear Benchmark (2025-12-27)

**All solvers converge** for tetrahedral nonlinear problems.

#### maxh=0.30m (200 elements, 600 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|--------:|---------:|------:|
| 726,577 | LU | 66 MB | - | 0 | 14 | 0 | 0.052s | - | 0.052s | 1.89s |
| 726,721 | BiCGSTAB | 66 MB | - | 171 | 12 | 0 | - | - | 0.063s | 0.62s |
| 726,810 | HACApK | 70 MB | 104% | 105 | 9 | 0.20s | - | 0.20s | 0.04s | 0.24s |

#### maxh=0.25m (390 elements, 1,170 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 727,050 | LU | 69 MB | - | 0 | 15 | 0 | 0.081s | - | 0.081s | 2.58s |
| 727,454 | BiCGSTAB | 66 MB | - | 93 | 7 | 0 | - | - | 0.131s | 2.43s |
| 727,463 | HACApK | 83 MB | 98% | 95 | 7 | 0.07s | - | 0.07s | 0.07s | 0.15s |

#### maxh=0.20m (627 elements, 1,881 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 748,822 | LU | 72 MB | - | 0 | 20 | 0 | 0.169s | - | 0.169s | 7.81s |
| 748,869 | BiCGSTAB | 68 MB | - | 118 | 14 | 0 | - | - | 0.346s | 6.70s |
| 748,858 | HACApK | 106 MB | 83% | 128 | 17 | 0.19s | - | 0.19s | 0.33s | 0.53s |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 730,715 | LU | 91 MB | - | 0 | 29 | 0 | 5.66s | - | 5.66s | 174.5s |
| 730,591 | BiCGSTAB | 77 MB | - | 134 | 36 | 0 | - | - | 2.10s | 101.5s |
| 730,616 | HACApK | 438 MB | **54%** | 134 | 36 | 1.97s | - | 1.97s | 3.08s | 5.17s |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | LU | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|-------:|---------:|---------:|------:|
| 730,996 | LU | 104 MB | - | 0 | 39 | 0 | 1435s | - | 1435s | 1823s |
| 731,003 | BiCGSTAB | 90 MB | - | 217 | 41 | 0 | - | - | 507s | 875s |
| 731,054 | HACApK | 1575 MB | **35%** | 205 | 38 | 7.08s | - | 7.08s | 12.5s | 20.2s |

---

## Key Findings

### Hexahedral Elements (Linear)

1. **All solvers converge**: LU, BiCGSTAB, HACApK all work
2. **HACApK is 5.6x faster than LU at N=15**: 27.9s vs 155.8s
3. **HACApK is 4.7x faster than BiCGSTAB at N=15**: 27.9s vs 132.6s

### Hexahedral Elements (Nonlinear)

1. **All solvers converge**: LU, BiCGSTAB, HACApK all work
2. **HACApK is 29x faster than LU at N=15**: 38.3s vs 1098s
3. **HACApK is 10x faster than BiCGSTAB at N=15**: 38.3s vs 372s
4. **N=20 solved only by HACApK**: 128s for 48,000 DOF
5. **M_avg_z consistent**: ~720,000 A/m across all solvers

### Tetrahedral Elements (Nonlinear)

1. **All solvers converge**: LU, BiCGSTAB, HACApK all work
2. **HACApK is fastest at large scale**: 90x faster than LU at maxh=0.10 (20.2s vs 1823s)
3. **HACApK is 43x faster than BiCGSTAB**: 20.2s vs 875s at maxh=0.10
4. **M_avg_z consistent**: ~731,000 A/m across all solvers

### Solver Recommendations

| Element Type | Material | Problem Size | Recommended Solver |
|--------------|----------|--------------|-------------------|
| Hexahedral | Linear | DOF < 6,000 | Any solver |
| Hexahedral | Linear | DOF > 6,000 | **HACApK** (5x faster) |
| Hexahedral | Nonlinear | DOF < 1,000 | Any solver |
| Hexahedral | Nonlinear | DOF 1,000-10,000 | BiCGSTAB or HACApK |
| Hexahedral | Nonlinear | DOF > 10,000 | **HACApK** (29x faster) |
| Tetrahedral | Nonlinear | DOF < 2,000 | Any solver |
| Tetrahedral | Nonlinear | DOF 2,000-10,000 | HACApK |
| Tetrahedral | Nonlinear | DOF > 10,000 | **HACApK** (90x faster) |

---

## H-matrix Statistics Summary (2025-12-27 Updated)

**Note**: Compression ratio = H-matrix memory / Dense memory. Lower is better.

### Hexahedral Elements (6 DOF per element)

| N | Elements | DOF | lowrank | dense | max_rank | H-mat [MB] | Dense [MB] | Compression |
|---|----------|-----|--------:|------:|---------:|-----------:|-----------:|------------:|
| 10 | 1,000 | 6,000 | 1,166 | 2,024 | 90 | 130 | 275 | **48%** |
| 15 | 3,375 | 20,250 | 8,262 | 9,124 | 79 | 766 | 3,129 | **25%** |
| 20 | 8,000 | 48,000 | 17,976 | 20,656 | 99 | 2,429 | 17,578 | **14%** |

### Tetrahedral Elements (3 DOF per element)

| maxh | Elements | DOF | lowrank | dense | max_rank | H-mat [MB] | Dense [MB] | Compression |
|------|----------|-----|--------:|------:|---------:|-----------:|-----------:|------------:|
| 0.30 | 200 | 600 | 40 | 414 | 33 | 2.9 | 2.8 | 104% |
| 0.25 | 390 | 1,170 | 290 | 869 | 36 | 10.3 | 10.4 | 98% |
| 0.20 | 627 | 1,881 | 554 | 1,544 | 41 | 22 | 27 | 83% |
| 0.15 | 2,211 | 6,633 | 3,160 | 11,226 | 46 | 180 | 336 | **54%** |
| 0.10 | 4,994 | 14,982 | 12,074 | 24,401 | 52 | 597 | 1,712 | **35%** |

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
python benchmark_tetrahedron.py --lu --bicgstab --hacapk 0.30 0.25 0.20 0.15 0.10
```

### Nonlinear Benchmarks

```bash
cd nonlinear/hexahedron
python benchmark_hex.py --lu --bicgstab --hacapk 5 10 15
python benchmark_hex.py --hacapk 20

cd nonlinear/tetrahedron
python benchmark_tetra.py --lu --bicgstab --hacapk 0.30 0.25 0.20 0.15 0.10
```

---

---

## Unit System Notes

### Radia Internal Unit System

**IMPORTANT**: Radia ALWAYS uses millimeters (mm) internally, regardless of `FldUnits()` setting.

| Setting | Coordinate Input | Field Output | Internal Storage |
|---------|------------------|--------------|------------------|
| `FldUnits('mm')` | mm | T (for B), A/m (for H) | mm |
| `FldUnits('m')` | m (scaled to mm internally) | T, A/m | mm |

### Vector Potential A Unit Conversion

When using `FldUnits('m')` for NGSolve integration:

- **B field**: Returned correctly in Tesla (no scaling needed)
- **H field**: Returned correctly in A/m (no scaling needed)
- **A field**: Returned in T*mm (needs conversion for curl(A) = B in meters)

**Why A needs scaling:**

1. Radia computes A using mm-based geometry internally
2. A is dimensionally [T*length], so A_internal = T*mm
3. NGSolve differentiates in meters: `curl(A) = dA/dx_m`
4. For correct B = curl(A): `A_m = A_mm / 1000`

**In radia_ngsolve:**

```cpp
// Vector potential A unit scaling (in radia_ngsolve.cpp)
// Radia ALWAYS uses mm internally, so A is always in T*mm
// NGSolve differentiates in meters: curl(A) = dA/dx_m
// To get correct B = curl(A), we scale A by 0.001:
double scale = (field_type == "a") ? 0.001 : 1.0;
```

### Maxwell Relation Verification

The test `examples/ngsolve_integration/verify_curl_A_equals_B/` verifies:

```
B = curl(A)
```

With proper A field scaling, the `|curl(A)|/|B|` ratio should be ~1.0.

---

**Last Updated**: 2025-12-27
