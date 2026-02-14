# Radia Solver Architecture: Design Philosophy

**Date**: 2026-02-14
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
|  FastHenry PEEC + SIBC     Radia MMM/MSC                         |
|  - L, R extraction         - Linear: MatLin(mu_r)                |
|  - Skin effect (Dowell)    - Nonlinear: MatSatIsoTab(BH)         |
|  - Loop-Star decomp.       - Hex/Tet/Wedge elements              |
|  - Surface mesh only       - Volume mesh                         |
|                                                                   |
|        |                           |                              |
|        +--- Coupling (Delta_L) ----+                              |
|        |   Biot-Savart + A-field   |                              |
|        v                           v                              |
|  +--------------------------------------------------+            |
|  |  Coupled System:                                  |            |
|  |  Z = R + jw(L_air + Delta_L * (mu_eff(w) - 1))  |            |
|  |                                                   |            |
|  |  mu_eff(w): from ngbem FEM eddy current solution  |            |
|  |  <Hz> = volume average of Hz from FEM diffusion   |            |
|  |  mu_eff = mu_r * <Hz> / Hz_inc                    |            |
|  +--------------------------------------------------+            |
|                                                                   |
|  Optional: ngbem (NGSolve BEM)                                    |
|  =============================                                    |
|  - Linear materials ONLY                                          |
|  - Galerkin BEM (LaplaceSL)                                       |
|  - FEM-BEM coupling (Calderon)                                    |
|  - High-order elements                                            |
|                                                                   |
+------------------------------------------------------------------+
```

---

## Component Roles

### 1. Conductor (Coil): FastHenry PEEC + SIBC

**Role**: Model coil windings, PCB traces, and conducting structures.

**Why FastHenry approach**:
- PEEC naturally gives circuit parameters (L, R, C, M)
- Surface mesh only (no volume mesh needed for conductors)
- SIBC captures skin effect without meshing skin depth
- Direct SPICE netlist output

**Formulation**:
```
Loop-Star block system:
| R + jwL       M_LS^T  |   | I_loop |   | V_port |
| M_LS      P/(jw)      | * | Q_star | = | 0      |

where:
  L = inductance matrix (Laplace kernel, BEM)
  R = resistance (DC + skin effect via SIBC/Dowell)
  P = potential coefficients (capacitive)
  M_LS = divergence coupling
```

**SIBC (Surface Impedance Boundary Condition)**:
- Rectangular conductors: Dowell formula `F_R = xi * [sinh(2xi) + sin(2xi)] / [cosh(2xi) - cos(2xi)]`
- Circular conductors: Bessel function `Z = (k*l)/(2*pi*r*sigma) * J0(kr)/J1(kr)`
- Skin depth: `delta = sqrt(2/(omega*mu*sigma))`

**When to use**:
- Power electronics (DC - 1 MHz)
- WPT coils (6.78 MHz, 13.56 MHz)
- Transformer/inductor windings
- Any application needing SPICE models

### 2. Magnetic Material (Core): Radia MMM/MSC

**Role**: Model magnetic cores (ferrite, iron, steel, permanent magnets).

**Why Radia MMM/MSC**:
- **Nonlinear materials supported** (B-H curve via MatSatIsoTab)
- Unbounded domain (no PML/ABC needed)
- Multiple solver methods (LU, BiCGSTAB, HACApK)
- IMA (Image Method) for symmetry exploitation

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
- `MatMagFixed(M)`: Permanent magnet (fixed M)

**Key advantage**: **Only Radia can handle nonlinear materials** in the integral equation framework.
FEM-BEM (ngbem) is limited to linear materials because BEM requires a known Green's function.

### 3. ngbem (NGSolve BEM): Linear Materials Only

**Role**: High-order Galerkin BEM for conductors and linear magnetic problems.

**Capabilities**:
- LaplaceSL on HDivSurface -> inductance L (PEEC)
- SingleLayerPotentialOperator on SurfaceL2 -> potential P
- FEM-BEM coupling via Calderon projector (eddy current)
- High-order elements (order 0, 1, 2, ...)
- FMM acceleration for large problems

**CRITICAL LIMITATION: Linear materials only**.

ngbem uses the Laplace (or Helmholtz) Green's function `G(r) = 1/(4*pi*r)`.
This requires the material properties to be spatially uniform within each domain.
Nonlinear materials (mu depends on H) cannot be handled by BEM because:
1. The Green's function assumes constant material properties
2. Nonlinear problems require iterative updates of material parameters
3. The boundary integral equation is only valid for linear, piecewise-homogeneous media

**When to use ngbem**:
- High-accuracy PEEC (Galerkin > Collocation)
- Linear eddy current problems (FEM interior + BEM exterior)
- Problems needing high-order convergence

**When NOT to use ngbem**:
- Nonlinear magnetic materials -> use Radia MMM/MSC
- Saturable cores -> use Radia MMM with MatSatIsoTab
- Problems with H-dependent permeability

### 4. Coupling: PEEC + MMM (Delta_L)

**Role**: Connect conductor impedance with magnetic core response.

**Coupling mechanism**:
```
L_total = L_air + Delta_L * (mu_eff(omega) - 1)

Delta_L[i,j] = mu_0 * integral_core H_i(r) . H_j(r) dV
```

where `H_i` is the field from unit current in loop i.

**mu_eff(omega)** is computed from ngbem FEM eddy current solution:
```
1. Solve FEM diffusion: laplacian(Hz) = j*omega*mu*sigma*Hz
2. Volume average: <Hz> = Integrate(Hz, mesh) / Volume
3. mu_eff = mu_r * <Hz> / Hz_inc
```

**For non-conducting core** (ferrite, sigma ~ 0):
```
mu_eff = mu_r (constant, frequency-independent)
L_total = L_air + Delta_L * (mu_r - 1)
```
Core increases inductance at all frequencies.

**For conducting core** (steel, iron, sigma >> 0):
```
mu_eff(omega) computed from FEM -> decreases with frequency
```
Eddy currents cause:
- Low freq: mu_eff ~ mu_r, L increases (full permeability)
- Mid freq: mu_eff decreases, loss peak (d/delta ~ 1)
- High freq: mu_eff ~ 0, L ~ L_air (core fully shielded)

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

### Why FastHenry-Style PEEC for Coils?

| Approach | Advantages | Limitations |
|----------|-----------|-------------|
| **FastHenry PEEC + SIBC** | Direct L,R; surface mesh; SPICE output; skin effect | MQS only |
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
| **Transformer** | PEEC + SIBC | MMM (MatSatIsoTab) | Delta_L |
| **WPT coil** | PEEC + SIBC | MMM (MatLin) | Delta_L |
| **Induction heating** | PEEC + SIBC | MMM + FEM eddy (mu_eff) | Delta_L + eddy |
| **PM motor** | - | MMM (MatMagFixed + MatLin) | Radia Solve() |
| **EMC shielding** | PEEC | - (or MMM for mu-metal) | Delta_L |
| **PCB trace** | PEEC + SIBC | - | - |

### By Material

| Material | Linear? | Solver | Material API |
|----------|---------|--------|-------------|
| Copper coil | N/A | PEEC + SIBC | sigma parameter |
| Ferrite core | Yes (usually) | MMM/MSC or ngbem | MatLin(mu_r) |
| Silicon steel | Nonlinear | **MMM/MSC only** | MatSatIsoTab(BH) |
| NdFeB PM | Fixed M | MMM/MSC | ObjHexahedron(v, M) |
| Aluminum shield | Linear | PEEC or ngbem | sigma parameter |
| Mu-metal | Nonlinear | **MMM/MSC only** | MatSatIsoTab(BH) |

### By Frequency

| Range | Conductor Model | Core Model | Notes |
|-------|----------------|------------|-------|
| DC | PEEC (R only) | MMM (static) | No skin effect |
| 50 Hz - 10 kHz | PEEC + SIBC | MMM + FEM eddy | Lamination matters |
| 10 kHz - 1 MHz | PEEC + SIBC | MMM + FEM eddy | Power electronics |
| 1 MHz - 100 MHz | PEEC + SIBC | MMM (ferrite, sigma~0) | WPT, RF |
| > 100 MHz | ngbem (full-wave) | - | EMC, antenna |

---

## Implementation Status

| Component | Status | File |
|-----------|--------|------|
| PEEC Loop-Star (ngbem) | **Implemented** | `src/radia/ngbem_peec.py` |
| FEM-BEM eddy current | **Implemented** (FEM mode) | `src/radia/ngbem_eddy.py` |
| Coupled PEEC+MMM | **Implemented** | `src/radia/ngbem_coupled.py` |
| SIBC (Dowell) | **Implemented** | `src/core/rad_peec_surface_impedance.cpp` |
| SIBC (Bessel) | **Implemented** | `examples/.../validate_circular_coil_sibc.py` |
| FEM eddy current mu_eff | **Implemented** | `src/radia/ngbem_coupled.py` |
| Radia MMM (nonlinear) | **Implemented** | `src/core/rad_relaxation_methods.cpp` |
| Radia MSC (hex/tet/wedge) | **Implemented** | `src/core/rad_polyhedron.cpp` |
| FastHenry C++ PEEC | **Implemented** | `src/core/rad_peec_matrices.cpp` |
| ngbem Galerkin PEEC | **Implemented** | `src/radia/ngbem_peec.py` |

---

## Key Principle: Complementary Solvers

```
Nonlinear core       -----> Radia MMM/MSC (ONLY option)
Linear core          -----> Radia MMM/MSC  OR  ngbem
Conductor coil       -----> FastHenry PEEC + SIBC (BEST for circuits)
Eddy current (linear)-----> ngbem FEM eddy (mu_eff from FEM)
High-order accuracy  -----> ngbem Galerkin BEM
SPICE extraction     -----> PEEC + Lanczos MOR (ONLY option)
```

The fundamental constraint is:
- **BEM requires linear materials** (known Green's function)
- **Nonlinear materials require volume-based methods** (MMM/MSC or FEM)
- **Radia's unique value**: integral equation + nonlinear iteration in unbounded domains
