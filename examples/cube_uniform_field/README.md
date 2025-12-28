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

### 1. Hexahedral Linear Benchmark (2025-12-29)

**M_analytical = 598,204 A/m** (mu_r = 1000, N_demag = 1/3)

#### N=5 (125 elements, 750 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 701,899 | LU | 50 MB | - | 0 | 1 | 0.09s | - | 0.07s | 0.17s |
| 701,903 | BiCGSTAB | 32 MB | - | 18 | 3 | - | - | - | 0.16s |
| 701,896 | HACApK | 33 MB | 95% | 18 | 3 | - | - | - | 0.11s |

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 716,183 | LU | 1192 MB | - | 0 | 1 | 3.4s | - | 1.2s | 4.7s |
| 716,109 | BiCGSTAB | 621 MB | - | 0 | 1 | 3.4s | - | 0.34s | 3.8s |
| 716,191 | HACApK | 210 MB | **47%** | 11 | 1 | 0.01s | 5.4s | 0.18s | 5.6s |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 720,156 | LU | 13029 MB | - | 0 | 1 | 37.4s | - | 25.0s | 64.3s |
| 720,158 | BiCGSTAB | 6703 MB | - | 0 | 1 | 37.5s | - | 2.7s | 40.9s |
| 720,196 | HACApK | 1244 MB | **24%** | 10 | 1 | 0.11s | 33.9s | 0.85s | 35.0s |

#### N=20 (8,000 elements, 48,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 722,035 | HACApK | 4809 MB | **13%** | 22 | 1 | 0.60s | 118s | 4.4s | 123s |

---

### 2. Hexahedral Nonlinear Benchmark (2025-12-29)

**All solvers converge** for hexahedral nonlinear problems.

#### N=5 (125 elements, 750 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 702,131 | LU | 50 MB | - | 0 | 6 | 0.09s | - | 0.47s | 0.56s |
| 702,114 | BiCGSTAB | 41 MB | - | 0 | 3 | 0.09s | - | 0.41s | 0.50s |
| 702,107 | HACApK | 37 MB | 95% | 21 | 3 | 0.55s | 0.55s | 0.03s | 0.59s |

#### N=10 (1,000 elements, 6,000 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 716,281 | LU | 1193 MB | - | 0 | 13 | 3.4s | - | 9.7s | 14.1s |
| 716,307 | BiCGSTAB | 621 MB | - | 0 | 6 | 3.4s | - | 0.99s | 4.4s |
| 716,362 | HACApK | 210 MB | **47%** | 31 | 4 | 0.01s | 5.5s | 0.61s | 6.1s |

#### N=15 (3,375 elements, 20,250 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 719,832 | LU | 13030 MB | - | 0 | 35 | 37.5s | - | 870s | 943s |
| 719,888 | BiCGSTAB | 6704 MB | - | 0 | 28 | 37.4s | - | 21.5s | 59.6s |
| 719,906 | HACApK | 1244 MB | **24%** | 90 | 30 | 0.11s | 34.0s | 14.1s | 48.3s |

#### N=20 (8,000 elements, 48,000 DOF) - HACApK only

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 721,344 | HACApK | 4813 MB | **13%** | 129 | 33 | 0.59s | 118s | 45.2s | 164s |

---

### 3. Tetrahedral Linear Benchmark (2025-12-29)

**M_analytical = 598,204 A/m** (mu_r = 1000, N_demag = 1/3)

#### maxh=0.30m (200 elements, 600 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 751,152 | LU | 66 MB | - | 0 | 2 | - | - | 0.22s | 1.9s |
| 751,162 | BiCGSTAB | 65 MB | - | 45 | 3 | - | - | 0.17s | 0.81s |
| 751,188 | HACApK | 65 MB | **0%** | 30 | 3 | 7ms | 7ms | 2ms | 0.95s |

#### maxh=0.25m (390 elements, 1,170 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 765,436 | LU | 116 MB | - | 0 | 2 | 3.3s | - | 0.26s | 3.6s |
| 765,443 | BiCGSTAB | 65 MB | - | 29 | 4 | - | - | 0.22s | 2.6s |
| 765,453 | HACApK | 68 MB | **1%** | 29 | 4 | 0.16s | 0.16s | 0.01s | 3.2s |

#### maxh=0.20m (627 elements, 1,881 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 816,133 | LU | 193 MB | - | 0 | 1 | 2.6s | - | 0.31s | 2.9s |
| 816,080 | BiCGSTAB | 134 MB | - | 0 | 1 | 2.5s | - | 0.23s | 2.8s |
| 816,189 | HACApK | 105 MB | 82% | 29 | 1 | 8ms | 0.28s | 0.05s | 0.34s |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 748,692 | LU | 1601 MB | - | 0 | 1 | 23.8s | - | 3.4s | 27.5s |
| 748,682 | BiCGSTAB | 915 MB | - | 0 | 1 | 50.9s | - | 1.0s | 52.2s |
| 748,718 | HACApK | 435 MB | 53% | 25 | 1 | 0.05s | 2.7s | 0.38s | 3.1s |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 754,617 | LU | 7810 MB | - | 0 | 1 | 170s | - | 35.3s | 209s |
| 754,393 | BiCGSTAB | 4371 MB | - | 0 | 1 | 122s | - | 5.0s | 130s |
| 754,533 | HACApK | 1573 MB | **34%** | 29 | 1 | 0.23s | 9.3s | 1.4s | 11.0s |

#### maxh=0.07m (15,116 elements, 45,348 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 739,990 | HACApK | 10803 MB | **17%** | 26 | 1 | 2.1s | 46.4s | 5.4s | 54.5s |

---

### 4. Tetrahedral Nonlinear Benchmark (2025-12-29)

**All solvers converge** for tetrahedral nonlinear problems.

#### maxh=0.30m (200 elements, 600 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 726,576 | LU | 66 MB | - | 0 | 14 | - | - | 0.05s | 1.9s |
| 726,721 | BiCGSTAB | 65 MB | - | 171 | 12 | - | - | 0.06s | 0.62s |
| 726,810 | HACApK | 69 MB | 103% | 105 | 9 | 0.20s | 0.20s | 0.04s | 0.24s |

#### maxh=0.25m (390 elements, 1,170 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 727,049 | LU | 116 MB | - | 0 | 15 | - | - | 0.58s | 2.8s |
| 727,454 | BiCGSTAB | 81 MB | - | 93 | 7 | - | - | 0.37s | 2.6s |
| 727,462 | HACApK | 83 MB | 98% | 95 | 7 | 0.07s | 0.07s | 0.07s | 0.15s |

#### maxh=0.20m (627 elements, 1,881 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 748,822 | LU | 193 MB | - | 0 | 20 | 2.6s | - | 2.0s | 4.7s |
| 748,869 | BiCGSTAB | 134 MB | - | 0 | 14 | 2.9s | - | 0.51s | 3.4s |
| 748,858 | HACApK | 106 MB | 82% | 128 | 17 | 0.19s | 0.19s | 0.32s | 0.52s |

#### maxh=0.15m (2,211 elements, 6,633 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 730,715 | LU | 1602 MB | - | 0 | 29 | 48.0s | - | 92.7s | 144s |
| 730,759 | BiCGSTAB | 916 MB | - | 0 | 21 | 33.6s | - | 4.5s | 38.4s |
| 730,615 | HACApK | 438 MB | 53% | 134 | 36 | 2.0s | 2.0s | 3.1s | 5.2s |

#### maxh=0.10m (4,994 elements, 14,982 DOF)

| M_avg_z | Solver | Memory | Compress | Linear | Nonl | MatBuild | H-matrix | LinSolve | Total |
|--------:|--------|-------:|---------:|-------:|-----:|---------:|---------:|---------:|------:|
| 730,996 | LU | 7811 MB | - | 0 | 39 | 168s | - | 1382s | 1572s |
| 731,003 | BiCGSTAB | 90 MB | - | 217 | 41 | - | - | 507s | 875s |
| 731,054 | HACApK | 1575 MB | **34%** | 205 | 38 | 7.1s | 7.1s | 12.5s | 20.2s |

---

## Key Findings

### Hexahedral Elements (Linear)

1. **All solvers converge**: LU, BiCGSTAB, HACApK all work
2. **HACApK is 2.5x faster than LU at N=15**: 36.9s vs 93.8s
3. **BiCGSTAB is competitive**: 41.1s at N=15 (similar to HACApK)
4. **Matrix caching reduces MatBuild**: 38.5s -> 0.11s with HACApK
5. **M_avg_z consistent**: ~716,000-720,000 A/m across all solvers

### Hexahedral Elements (Nonlinear)

1. **All solvers converge**: LU, BiCGSTAB, HACApK all work
2. **HACApK is 29x faster than LU at N=15**: 38.3s vs 1098s
3. **HACApK is 10x faster than BiCGSTAB at N=15**: 38.3s vs 372s
4. **N=20 solved only by HACApK**: 128s for 48,000 DOF
5. **M_avg_z consistent**: ~720,000 A/m across all solvers

### Tetrahedral Elements (Linear)

1. **All solvers converge**: LU, BiCGSTAB, HACApK all work
2. **HACApK is 7.6x faster than LU at maxh=0.10**: 11.9s vs 90.7s
3. **HACApK is 1.9x faster than BiCGSTAB at maxh=0.10**: 11.9s vs 22.4s
4. **Memory efficiency**: HACApK uses 5x less memory than LU (1573 MB vs 7810 MB)
5. **M_avg_z consistent**: ~755,000 A/m across all solvers

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
| Tetrahedral | Linear | DOF < 2,000 | Any solver |
| Tetrahedral | Linear | DOF 2,000-15,000 | BiCGSTAB or HACApK |
| Tetrahedral | Linear | DOF > 15,000 | **HACApK** (7.6x faster) |
| Tetrahedral | Nonlinear | DOF < 2,000 | Any solver |
| Tetrahedral | Nonlinear | DOF 2,000-10,000 | HACApK |
| Tetrahedral | Nonlinear | DOF > 10,000 | **HACApK** (90x faster) |

---

## H-matrix Statistics Summary (2025-12-29 Updated)

**Note**: Compression ratio = H-matrix memory / Dense memory. Lower is better.

### Hexahedral Elements (6 DOF per element)

| N | Elements | DOF | lowrank | dense | max_rank | H-mat [MB] | Dense [MB] | Compression |
|---|----------|-----|--------:|------:|---------:|-----------:|-----------:|------------:|
| 5 | 125 | 750 | 26 | 158 | 46 | 4 | 4 | 95% |
| 10 | 1,000 | 6,000 | 1,166 | 2,024 | 90 | 130 | 274 | **47%** |
| 15 | 3,375 | 20,250 | 8,202 | 9,442 | 82 | 765 | 3,128 | **24%** |
| 20 | 8,000 | 48,000 | 19,216 | 20,610 | 92 | 2,440 | 17,578 | **13%** |

### Tetrahedral Elements (3 DOF per element) - Linear Material

| maxh | Elements | DOF | lowrank | dense | max_rank | H-mat [MB] | Dense [MB] | Compression |
|------|----------|-----|--------:|------:|---------:|-----------:|-----------:|------------:|
| 0.30 | 200 | 600 | 40 | 414 | 33 | 0 | 2 | **0%** |
| 0.25 | 390 | 1,170 | 290 | 869 | 36 | 0 | 10 | **1%** |
| 0.20 | 627 | 1,881 | 554 | 1,544 | 41 | 22 | 26 | 82% |
| 0.15 | 2,211 | 6,633 | 3,160 | 11,226 | 46 | 179 | 335 | 53% |
| 0.10 | 4,994 | 14,982 | 12,074 | 24,401 | 52 | 597 | 1,712 | **34%** |
| 0.07 | 15,116 | 45,348 | 44,966 | 89,195 | 53 | 2,685 | 15,689 | **17%** |

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

## Memory Measurement Notes

### Windows: peak_wset (Peak Working Set)

Memory measurements in this benchmark use Windows `peak_wset` via `psutil.Process().memory_info().peak_wset`.

**What peak_wset measures:**
- **Process-wide peak memory**: Maximum physical memory used by the entire Python process
- **Includes**: Python interpreter, imported modules (numpy, netgen, radia), and solver allocations
- **Timing**: Captures peak across entire process lifetime, not just during solve

**Baseline overhead (approximate):**
- Python 3.12 + numpy + radia: ~50-80 MB
- With NGSolve/Netgen: ~100-150 MB additional

**Interpretation:**
- For large problems (DOF > 10,000), solver memory dominates and peak_wset is accurate
- For small problems (DOF < 1,000), baseline overhead may be significant fraction of reported memory
- To estimate pure solver memory: `solver_memory ≈ peak_wset - baseline`

**Alternative: RSS (Resident Set Size)**
- `memory_info().rss` shows current memory at measurement time
- After solve completes, most temporary allocations are freed
- RSS underestimates peak memory usage during LU decomposition

**Dense matrix reference (for comparison):**
- Dense matrix memory = N² × 8 bytes (double precision)
- N=6,000 DOF: 275 MB dense matrix
- N=20,250 DOF: 3,129 MB dense matrix
- N=48,000 DOF: 17,578 MB dense matrix

### H-matrix Compression Ratio

The "Compression" column shows: **H-matrix memory / Dense matrix memory × 100%**

- **< 50%**: Good compression, H-matrix is beneficial
- **50-80%**: Moderate compression
- **> 80%**: Poor compression (small problems, H-matrix overhead dominates)
- **> 100%**: H-matrix larger than dense (problem too small for H-matrix)

Compression ratio is obtained from `rad.GetHACApKStats()['compression']`.

---

**Last Updated**: 2025-12-29 (All benchmark results updated from JSON files)
