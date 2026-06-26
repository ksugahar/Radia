# Radia Solver Architecture: Design Philosophy

**Date**: 2026-02-16
**Status**: Active Design Document

---

## Overview

Radia provides coupled electromagnetic analysis by combining specialized solvers,
each optimized for its domain. The key principle: **use the right tool for each physics**.

```
+------------------------------------------------------------------+
|                    Radia Solver Architecture                      |
+------------------------------------------------------------------+
|                                                                   |
|  Conductor (Coil)          Magnetic Material (Core)               |
|  ==================        ===========================            |
|  PEEC + SIBC (C++)         Radia MMM/MSC                         |
|  - L, R, C extraction      - Linear: MatLin(mu_r)                |
|  - Skin effect (SIBC)      - Nonlinear: MatSatIsoTab(BH)         |
|  - Node-segment topology   - Hex/Tet/Wedge elements              |
|  - MNA multi-port solver   - Volume mesh                         |
|  - FastHenry .inp parser                                         |
|                                                                   |
|        |                           |                              |
|        +--- Coupling (Delta_L) ----+                              |
|        |   Biot-Savart + A-field   |                              |
|        v                           v                              |
|  +--------------------------------------------------+            |
|  |  Coupled System:                                  |            |
|  |  Z_eff = diag(R + Zs) + jw * (L_air + Delta_L)  |            |
|  |                                                   |            |
|  |  Delta_L: from Radia Solve() of magnetic core     |            |
|  +--------------------------------------------------+            |
|                                                                   |
|  Optional: ngbem (NGSolve BEM)                                    |
|  =============================                                    |
|  - Linear materials ONLY                                          |
|  - Galerkin BEM (LaplaceSL)                                       |
|  - Stabilized low-frequency formulation                           |
|  - High-order elements                                            |
|                                                                   |
+------------------------------------------------------------------+
```

---

## Component Roles

### 1. Conductor (Coil): PEEC + SIBC (C++ MKL)

**Role**: Model coil windings, PCB traces, and conducting structures.

**Architecture**: Filament + Panel decomposition (FastImp-style, no Loop-Star needed).

**C++ Implementation** (`src/core/rad_peec_matrices.cpp`):
- LAPACK `zgesv_`/`zgetrf_`/`zgetrs_` for LU factorization (complex)
- Templated BiCGSTAB (`rad_bicgstab.h`) shared with MSC solver
- MKL `cblas_zgemm` for matrix multiplication
- MNA (Modified Nodal Analysis) multi-port solver

**Formulation**:
```
Filament-Panel block system (no Loop-Star transformation):
| R + jwL + Zs    jwM_LS  |   | I_filament |   | V |
| jwM_LS^T        P/(jw)  | * | Q_panel    | = | 0 |

MNA multi-port solver:
  Z_branch = diag(R_dc + Zs) + jw*L
  Y_branch = Z_branch^{-1}        (via LU or BiCGSTAB)
  Y_node = A_full * Y_branch * A_full^T
  Z_port from Y_reduced LU factorization

where:
  L = filament inductance (Neumann/Rosa-Grover analytical)
  R = DC resistance + surface impedance (SIBC)
  P = panel potential coefficients (Wilton formula)
  M_LS = filament-panel magnetic coupling
  A_full = node incidence matrix (from topology)
```

**SIBC (Surface Impedance Boundary Condition)**:
- Rectangular conductors: Dowell formula `F_R = xi * [sinh(2xi) + sin(2xi)] / [cosh(2xi) - cos(2xi)]`
- Circular conductors: Bessel function `Z = (k*l)/(2*pi*r*sigma) * I0(kr)/I1(kr)` (scipy.special.iv)
- Skin depth: `delta = sqrt(2/(omega*mu*sigma))`

**Node-Segment Topology API**:
```python
from peec_matrices import PyPEECBuilder
from peec_topology import PEECCircuitSolver

builder = PyPEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.1, 0, 0)
builder.add_connected_segment(n1, n2, w=1e-3, h=1e-3, sigma=5.8e7, nwinc=3, nhinc=3)
builder.add_port(n1, n2)
topo = builder.build_topology()

solver = PEECCircuitSolver(topo)
Z = solver.compute_port_impedance(freq=1e6)
```

**Multi-filament subdivision**: `nwinc`/`nhinc` parameters create parallel sub-filaments for skin/proximity effect.

**When to use**:
- Power electronics (DC - 1 MHz)
- WPT coils (6.78 MHz, 13.56 MHz)
- Transformer/inductor windings
- Any application needing SPICE models or circuit parameters

### 2. Magnetic Material (Core): Radia MMM/MSC

**Role**: Model magnetic cores (ferrite, iron, steel, permanent magnets).

**Why Radia MMM/MSC**:
- **Nonlinear materials supported** (B-H curve via MatSatIsoTab)
- Unbounded domain (no PML/ABC needed)
- Multiple solver methods (LU, BiCGSTAB, HACApK)
- IMA (Image Method) for symmetry exploitation (all solvers, `image='+x-z'`)
- NGSolve TaskManager parallelization (no OpenMP dependency)

**Element types**:

| Element | DOF | Method | Use Case |
|---------|-----|--------|----------|
| Hexahedron | 6 (sigma/face) | MSC | Structured meshes, PM blocks |
| Wedge | 5 (sigma/face) | MSC | Transition elements |
| Tetrahedron | 3 (Mx,My,Mz) | MMM | Complex geometry (Netgen) |

**Material models**:
- `MatLin(mu_r)`: Linear isotropic
- `MatLin([mu_par, mu_perp], axis)`: Linear anisotropic
- `MatSatIsoTab(BH_data)`: **Nonlinear** (B-H curve)
- `MatPlayHysteresis(K, eta, f_k_tables)`: **Vector hysteresis** (B-input Play, recommended)
- `MatEnergyHysteresis(K, eta, f_k_tables, eps)`: **Vector hysteresis** (energy-based Play)
- `MatPM(Br, Hc, axis)` or direct M via `ObjHexahedron(verts, M)`: Permanent magnet

**Key advantage**: **Only Radia can handle nonlinear materials** in the integral equation framework.
This includes **vector hysteresis** with Play operators -- a novel capability
for BEM.
FEM-BEM (ngbem) is limited to linear materials because BEM requires a known Green's function.

### 3. ngbem (NGSolve BEM): Linear Materials Only

**Role**: High-order Galerkin BEM for conductors and linear magnetic problems.

**Capabilities**:
- LaplaceSL on HDivSurface -> inductance L (PEEC)
- SingleLayerPotentialOperator on SurfaceL2 -> potential P
- Stabilized low-frequency formulation (Weggler 2026)
- High-order elements (order 0, 1, 2, ...)

**CRITICAL LIMITATION: Linear materials only**.

ngbem uses the Laplace (or Helmholtz) Green's function `G(r) = 1/(4*pi*r)`.
This requires the material properties to be spatially uniform within each domain.
Nonlinear materials (mu depends on H) cannot be handled by BEM because:
1. The Green's function assumes constant material properties
2. Nonlinear problems require iterative updates of material parameters
3. The boundary integral equation is only valid for linear, piecewise-homogeneous media

**Low-frequency stabilized formulation** (ngbem 2026):
- Block system `[A_k, Q_k; Q_k^T, k^2*V_k]` with O(1) condition number for all k
- Uses product space: HDivSurface x SurfaceL2
- Eliminates classical O(k^{-2}) blow-up at low frequencies

**When to use ngbem**:
- High-accuracy PEEC (Galerkin > Collocation)
- Linear eddy current problems (FEM interior + BEM exterior)
- Problems needing high-order convergence

**When NOT to use ngbem**:
- Nonlinear magnetic materials -> use Radia MMM/MSC
- Saturable cores -> use Radia MMM with MatSatIsoTab
- Problems with H-dependent permeability

### 4. Coupling: PEEC + MMM (CoupledPEECSolver)

**Role**: Connect conductor impedance with magnetic core response.

**Implementation**: `src/radia/peec_coupled.py` - CoupledPEECSolver class.

**Coupling mechanism**:
```
Z_eff(f) = diag(R + Zs(f)) + jw * (L_air + Delta_L)

Delta_L[i,j]:
  1. Unit current in segment j -> H-field via Biot-Savart (finite filament)
  2. H-field magnetizes material via rad.ObjBckg() + rad.Solve()
  3. Vector potential A from magnetized material: rad.Fld(mag_obj, 'a', point)
  4. Delta_L[i][j] = dot(A(center_i), dir_i) * length_i
```

**Key Property**: For linear materials, Delta_L is frequency-independent (computed once).

**Python API**:
```python
from peec_coupled import CoupledPEECSolver

solver = CoupledPEECSolver(topology_dict, magnetic_objects=[core_id])
solver.compute_coupling_matrix()  # N_seg Radia Solve calls
Z = solver.compute_port_impedance(freq)
Z_sweep = solver.frequency_sweep(freqs)
```

**FastHenry .magnetic block** support:
```
.magnetic
  type=box
  center=0.05,0.01,0.0
  size=0.06,0.01,0.01
  divisions=2,1,1
  mu_r=1000
.endmagnetic
```

---

## PEEC Solver Details

### C++ MNA Solver (rad_peec_matrices.cpp)

The PEEC MNA solver is implemented entirely in C++ using MKL LAPACK/BLAS:

| Operation | LU (Method 0) | BiCGSTAB (Method 1) |
|-----------|---------------|---------------------|
| Z_branch inversion | `zgesv_` | `bicgstab::DenseInvert<complex>` |
| Y_node assembly | `cblas_zgemm` | `cblas_zgemm` |
| Y_reduced factorization | `zgetrf_` | `zgetrf_` |
| Multi-RHS solve | `zgetrs_` | `zgetrs_` |

**Templated BiCGSTAB** (`rad_bicgstab.h`):
- Shared between MSC (real, `double`) and PEEC (complex, `std::complex<double>`)
- BLAS dispatch via C++ overloading in `radia::blas` namespace
- `cblas_d*` functions for real, `cblas_z*` for complex
- Non-LAPACK fallback for builds without MKL

**Solver method selection**:
```python
solver = PEECCircuitSolver(topo)
solver.set_solver_method(0)  # LU (default, LAPACK zgesv_)
solver.set_solver_method(1)  # BiCGSTAB (templated, shared with MSC)
solver.set_bicgstab_params(tol=1e-10, max_iter=1000)
```

### FastHenry .inp Parser

`src/radia/fasthenry_parser.py` parses FastHenry input files:

| Directive | Description |
|-----------|-------------|
| `.Units` | Length unit (m, cm, mm, um, in, mils) |
| `N<name>` | Node definition |
| `E<name>` | Segment (w, h, sigma, nwinc, nhinc) |
| `.external` | Port definition |
| `.freq` | Frequency sweep |
| `.default` | Default parameters |
| `.equiv` | Node merge |
| `.magnetic` | Magnetic material block |

**One-step solve**:
```python
from fasthenry_parser import FastHenryParser

parser = FastHenryParser()
parser.parse_file('inductor.inp')
result = parser.solve()  # Returns dict: freqs, Z_port, R, L, topology
```

---

## Design Decisions

### Why Not FEM for Everything?

NGSolve (FEM) excels at bounded domain problems with complex geometry.
But for electromagnetic component design, FEM has limitations:

| Aspect | FEM (NGSolve) | Radia (BEM/Integral) |
|--------|---------------|---------------------|
| Open boundary | Requires PML/ABC | Natural (Green's function) |
| Circuit parameters | Post-processing | Direct output (L, R, C, M) |
| Thin conductors | Mesh aspect ratio issues | PEEC (surface only) |
| SPICE export | Not supported | Built-in |
| Permanent magnets | Needs volume mesh | Analytical (ObjRecMag) |
| Nonlinear materials | **Supported** | **Supported (MMM/MSC)** |

### Why Not BEM for Everything?

BEM (ngbem) is powerful but has fundamental limitations:

| Aspect | BEM (ngbem) | Radia MMM/MSC |
|--------|-------------|---------------|
| Linear materials | **Supported** | **Supported** |
| Nonlinear materials | **NOT supported** | **Supported** |
| Green's function | Required (limits material models) | Uses interaction matrix |
| Iterative nonlinear solve | Not possible | BiCGSTAB with relaxation |
| High-order elements | Supported | Not needed (piecewise constant) |

### Why PEEC for Coils?

| Approach | Advantages | Limitations |
|----------|-----------|-------------|
| **PEEC + SIBC** | Direct L,R; surface mesh; SPICE output; skin effect | MQS only |
| FEM (NGSolve) | Full physics; complex geometry | Volume mesh; no direct circuit params |
| ngbem Galerkin | High accuracy; high-order | Linear only; no SIBC |

For power electronics and WPT, **PEEC + SIBC is the optimal choice**:
1. Engineers need L, R values directly (not fields)
2. Skin effect is well-modeled by SIBC (no volume mesh needed)
3. SPICE integration is essential for system-level simulation

---

## Solver Selection Guide

### By Application

| Application | Conductor | Core | Coupling |
|-------------|-----------|------|----------|
| **Transformer** | PEEC + SIBC | MMM (MatSatIsoTab) | CoupledPEECSolver |
| **WPT coil** | PEEC + SIBC | MMM (MatLin) | CoupledPEECSolver |
| **Induction heating** | PEEC + SIBC | MMM + ESIM (mu_eff) | CoupledPEECSolver |
| **PM motor** | - | MMM (direct M / MatPM + MatLin) | Radia Solve() |
| **EMC shielding** | PEEC | - (or MMM for mu-metal) | CoupledPEECSolver |
| **PCB trace** | PEEC + SIBC | - | - |

### By Material

| Material | Linear? | Solver | Material API |
|----------|---------|--------|-------------|
| Copper coil | N/A | PEEC + SIBC | sigma parameter |
| Ferrite core | Yes (usually) | MMM/MSC or ngbem | MatLin(mu_r) |
| Silicon steel | Nonlinear | **MMM/MSC only** | MatSatIsoTab(BH) |
| Soft iron (hysteresis) | Nonlinear | **MMM only** | MatPlayHysteresis / MatEnergyHysteresis + b_input_hantila |
| NdFeB PM | Fixed M | MMM/MSC | ObjHexahedron(v, M) |
| Aluminum shield | Linear | PEEC or ngbem | sigma parameter |
| Mu-metal | Nonlinear | **MMM/MSC only** | MatSatIsoTab(BH) |

### By Frequency

| Range | Conductor Model | Core Model | Notes |
|-------|----------------|------------|-------|
| DC | PEEC (R only) | MMM (static) | No skin effect |
| 50 Hz - 10 kHz | PEEC + SIBC | MMM | Lamination matters |
| 10 kHz - 1 MHz | PEEC + SIBC | MMM | Power electronics |
| 1 MHz - 100 MHz | PEEC + SIBC | MMM (ferrite, sigma~0) | WPT, RF |
| > 100 MHz | ngbem (full-wave) | - | EMC, antenna |

---

## Implementation Status

| Component | Status | File |
|-----------|--------|------|
| PEEC C++ MNA solver (LU) | **Implemented** | `src/core/rad_peec_matrices.cpp` |
| PEEC C++ BiCGSTAB | **Implemented** | `src/core/rad_bicgstab.h` |
| PEEC Node-segment topology | **Implemented** | `src/core/rad_peec_matrices.cpp` |
| PEEC Multi-filament (nwinc/nhinc) | **Implemented** | `src/core/rad_peec_matrices.cpp` |
| FastHenry .inp parser | **Implemented** | `src/radia/fasthenry_parser.py` |
| Coupled PEEC+MMM | **Implemented** | `src/radia/peec_coupled.py` |
| SIBC (Dowell) | **Implemented** | `src/core/rad_peec_surface_impedance.cpp` |
| SIBC (Bessel) | **Implemented** | Python: `scipy.special.iv` |
| Panel (capacitance) | **Implemented** | `src/core/rad_peec_matrices.cpp` |
| Radia MMM (nonlinear) | **Implemented** | `src/core/rad_relaxation_methods.cpp` |
| B-input Newton (hysteresis) | **Implemented** | `src/core/rad_relaxation_methods.cpp` |
| B-input Hantila hybrid (hysteresis) | **Implemented** | `src/core/rad_relaxation_methods.cpp` |
| Radia MSC (hex/tet/wedge) | **Implemented** | `src/core/rad_polyhedron.cpp` |
| Templated BiCGSTAB (real+complex) | **Implemented** | `src/core/rad_bicgstab.h` |
| ngbem Galerkin PEEC | **Implemented** | `src/radia/ngbem_peec.py` |
| ngbem FEM-BEM eddy current | **Implemented** | `src/radia/ngbem_eddy.py` |
| ngbem PEEC+MMM coupling | **Implemented** | `src/radia/ngbem_coupled.py` |

### Validation Results

| Test Suite | Tests | Status |
|------------|-------|--------|
| Topology (series/parallel/DC) | 4/4 | PASS |
| Coupling (2-port transformer) | 31/31 | PASS |
| FastHenry parser | 9/9 | PASS |
| Multi-filament (nwinc/nhinc) | 6/6 | PASS |
| Panel/Resonance | 23/23 | PASS |
| **Total** | **73/73** | **ALL PASS** |

BiCGSTAB vs LU numerical equivalence: < 1e-6% relative difference on all test cases.

---

## Key Principle: Complementary Solvers

```
Nonlinear core       -----> Radia MMM/MSC (ONLY option)
Linear core          -----> Radia MMM/MSC  OR  ngbem
Conductor coil       -----> PEEC + SIBC (BEST for circuits)
Eddy current (linear)-----> ngbem FEM eddy
High-order accuracy  -----> ngbem Galerkin BEM
SPICE extraction     -----> PEEC + Lanczos MOR (ONLY option)
Circuit parameters   -----> PEEC MNA solver (C++ LAPACK)
```

The fundamental constraint is:
- **BEM requires linear materials** (known Green's function)
- **Nonlinear materials require volume-based methods** (MMM/MSC or FEM)
- **Radia's unique value**: integral equation + nonlinear iteration in unbounded domains
- **PEEC's unique value**: direct circuit parameter extraction (L, R, C, M) from geometry

---

## Nonlinear Solver Options

### Standard Picard / Newton (B-H Curves)

For `MatLin` and `MatSatIsoTab` materials, `rad.Solve()` uses chi-based iteration:

```python
rad.SolverConfig(newton_method=False)   # Picard (default)
rad.SolverConfig(newton_method=True)    # Newton (faster, needs dM/dH)
rad.Solve(container, 1e-4, 100, 0)
```

### B-input Newton / Hantila (Hysteresis)

For `MatPlayHysteresis` or `MatEnergyHysteresis` materials, specialized B-input solvers are available.

```python
# B-input Newton (quadratic convergence, O(N^3) per iteration)
rad.SolverConfig(b_input_newton=True)

# B-input Hantila hybrid (recommended for hysteresis)
# Newton warmup + Hantila refinement
rad.SolverConfig(b_input_hantila=True, hantila_alpha=0)  # auto-alpha
rad.Solve(container, 1e-4, 5000, 0)
```

---

## Parallelization

Radia uses **NGSolve TaskManager** (work-stealing thread pool) for all parallelism.
No OpenMP dependency. See `src/core/rad_parallel.h` for the abstraction layer.

| Solver | Parallelization Strategy |
|--------|------------------------|
| LU (method=0) | `SuspendTaskManager` + `MKLThreadGuard` for multi-threaded `dgesv_` |
| BiCGSTAB (method=1) | `ParallelFor` for matrix-vector products |
| HACApK (method=2) | `ParallelFor` for H-matrix build, ACA+ compression, and BiCGSTAB |
| Field computation | `ParallelFor` for `Fld`, `FldLst`, analytical integrals |

Thread count: controlled by `ngsolve.SetNumThreads(n)` or TaskManager default (all cores).
Query via `rad.GetSolveStats()` → `num_threads`, `taskmanager_enabled`.
