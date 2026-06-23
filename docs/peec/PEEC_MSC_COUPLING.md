# PEEC-MSC Coupled Solver: Theory and Implementation

**Status**: Design Document (2026-02-12)
**Purpose**: PEEC conductor + MSC magnetic material coupled analysis

---

## Table of Contents

1. [Overview](#overview)
2. [Governing Equations](#governing-equations)
3. [Coupling Mechanism](#coupling-mechanism)
4. [Matrix Structure](#matrix-structure)
5. [Relation to CplMag](#relation-to-cplmag)
6. [Symmetry and Reciprocity](#symmetry-and-reciprocity)
7. [Implementation Roadmap](#implementation-roadmap)

---

## Overview

### What is PEEC-MSC Coupling?

**PEEC-MSC** combines two formulations for unified electromagnetic analysis:

| Method | Physics | Unknowns | Elements |
|--------|---------|----------|----------|
| **PEEC** | Conductor (coil, winding) | Current density **J** | Surface mesh |
| **MSC** | Magnetic material (core, ferrite) | Surface charge σ or Magnetization **M** | Volume mesh (hex/tet) |

**Applications**:
- Inductive components (transformers, inductors with ferrite cores)
- Wireless power transfer (WPT) coils + ferrite shields
- Induction heating (coil + workpiece)
- Motors/generators (windings + iron cores)

### Why Coupling is Needed

**Magnetic coupling**:
1. **PEEC current → MSC magnetization**: Coil current generates external field, magnetizing the core
2. **MSC magnetization → PEEC voltage**: Core magnetization induces back-EMF in coil

Without coupling, we lose:
- Inductance change due to core saturation
- Core losses (eddy current, hysteresis)
- Frequency-dependent impedance from core material

---

## Governing Equations

### PEEC Equation (Conductor)

Loop-Star decomposition + SIBC (Surface Impedance Boundary Condition):

```
Z_peec(ω) · I = V_source + V_induced
```

Where:
- **Z_peec**: PEEC impedance matrix (Loop-Star basis, includes R, L, C)
- **I**: Loop current vector (Loop-Star basis)
- **V_source**: Applied voltage (port excitation)
- **V_induced**: Induced voltage from MSC magnetic field

**SIBC**: Skin effect via frequency-dependent surface impedance:
```
Z_s(ω) = R_s(ω) + jω L_internal(ω)
```

### MSC Equation (Magnetic Material)

For hexahedral elements (6-DOF MSC):
```
K_msc · σ = H_ext
```

For tetrahedral elements (3-DOF MMM):
```
K_mmm · M = H_ext
```

Where:
- **K_msc**: MSC interaction matrix (6x6 blocks for hex-hex)
- **σ**: Surface charge density vector (1 per face)
- **K_mmm**: MMM demagnetization tensor (3x3 blocks for tet-tet)
- **M**: Magnetization vector (Mx, My, Mz)
- **H_ext**: External magnetic field from PEEC currents

**Nonlinear materials**:
```
B = μ(H) · H  →  M = χ(H) · H = (μ_r(H) - 1) · H
```

Solved iteratively (Newton-Raphson with line search damping).

---

## Coupling Mechanism

### 1. PEEC Current → MSC External Field

**Biot-Savart Law**: PEEC surface current generates magnetic field at MSC element evaluation points.

For PEEC loop current **I_k** flowing in filament **C_k**:
```
H_ext(r) = (1/4π) ∫_{C_k} I_k · (dl × R) / |R|^3
```

Where **R = r - r'** is the vector from source point **r'** to observation point **r**.

**Matrix form**:
```
H_ext = B_pm · I
```

Where **B_pm** (PEEC-to-MSC) is the Biot-Savart coupling matrix:
- **B_pm[i, k]**: Field at MSC element i's evaluation point from PEEC loop k

**Evaluation points** (MSC):
- **Hexahedron**: Midpoint between face center and element center (EIEM2 convention)
- **Tetrahedron**: Element centroid

### 2. MSC Magnetization → PEEC Induced Voltage

**Faraday's Law**: Changing magnetic flux induces EMF in PEEC loop.

For MSC element magnetization **M_i**, induced voltage in PEEC loop **k**:
```
V_induced,k = -jω · Φ_k = -jω ∫_{S_k} B_i · dS
```

Where:
- **Φ_k**: Magnetic flux through loop k from MSC element i
- **B_i = μ_0 · M_i** (or **B_i = μ_0 · (H_ext + M_i)** for linear materials)

**Neumann Formula** (mutual inductance):
```
M_ki = (μ_0/4π) ∮_{C_k} ∮_{C_i} (dl_k · dl_i) / |R|
```

For MSC surface charges (hex/wedge), use **equivalent current loops** on element faces:
```
J_surface = ∇ × M  →  K_surface = σ × n
```

**Matrix form**:
```
V_induced = -jω · M_mp · σ  (for MSC)
V_induced = -jω · M_mp · M  (for MMM)
```

Where **M_mp** (MSC-to-PEEC) is the mutual inductance matrix:
- **M_mp[k, i]**: Flux linkage in PEEC loop k from MSC element i

### 3. Reciprocity

**Theoretical reciprocity**:
```
M_mp = B_pm^T  (exact for infinitesimal elements)
```

In practice:
- **Evaluation points differ**: PEEC uses loop filaments, MSC uses element centers/faces
- **Numerical asymmetry**: ~1-5% difference is typical
- **Symmetrization** (optional): `K_coupling = 0.5 * (B_pm + M_mp^T)`

---

## Matrix Structure

### Full Coupled System

```
┌─────────────┬─────────────┐ ┌───┐   ┌─────────┐
│  Z_peec     │  -jω·M_mp   │ │ I │   │ V_source│
│             │             │ │   │ = │         │
│  B_pm       │  K_msc      │ │ σ │   │    0    │
└─────────────┴─────────────┘ └───┘   └─────────┘
```

**Dimensions**:
- **Z_peec**: (n_loops × n_loops) - PEEC impedance
- **M_mp**: (n_loops × n_msc_dof) - MSC→PEEC coupling
- **B_pm**: (n_msc_dof × n_loops) - PEEC→MSC coupling
- **K_msc**: (n_msc_dof × n_msc_dof) - MSC interaction matrix

**n_msc_dof**:
- Hexahedron: 6 DOF/element (6 faces)
- Wedge / pyramid: 5 DOF/element (5 faces)
- Tetrahedron / RecMag: 3 DOF/element (MMM Mx, My, Mz)
- Mixed surface-charge meshes (hex+wedge+pyramid) use variable DOF per element.
- A single soft-iron `rad.Solve` that mixes MMM tet/RecMag elements with MSC hex/wedge/pyramid elements is rejected with `Radia::Error204`; split the solve or use the mesh-backed HDiv-VIM path.

### Linear vs Nonlinear

**Linear materials**:
```
K_msc = constant  →  Direct solve (LU, BiCGSTAB, HACApK)
```

**Nonlinear materials** (μ_r(H)):
```
Newton-Raphson iteration:
1. Compute H_ext from PEEC currents: H_ext = B_pm · I
2. Update K_msc from μ_r(H) via Jacobian
3. Solve coupled system
4. Check convergence: |ΔM| / |M| < tol
```

**Damping** (for difficult convergence):
```
M_new = M_old + α · ΔM  (0 < α ≤ 1)
```

Where α is line search damping factor.

---

## Relation to Legacy CplMag

The old C++ CplMag APIs (`CndLoop`, `CplMagCreate`, `CplMagSolve`, etc.) are **removed**.
Use the Python-based `CoupledPEECSolver` and `PEECBuilder` instead.

### PEEC-MSC Architecture

| Aspect | Description |
|--------|------------|
| Conductor | Arbitrary PEEC mesh via `PEECBuilder` |
| Magnetic | MSC (hex 6-DOF) + MMM (tet 3-DOF) |
| Coupling | General B_pm, M_mp matrices |
| Mesh | Netgen/Cubit mixed mesh |
| Solver | LU, BiCGSTAB, HACApK |

```python
from radia.peec_topology import PEECCircuitSolver
from radia.peec_coupled import CoupledPEECSolver

solver = CoupledPEECSolver(topology_dict, magnetic_objects=[core_id])
Z = solver.compute_port_impedance(freq)
```

---

## Symmetry and Reciprocity

### Theoretical Background

**Maxwell reciprocity theorem**: For linear media,
```
∫ E_1 · J_2 dV = ∫ E_2 · J_1 dV
```

Implies:
```
Z_12 = Z_21  (impedance matrix is symmetric)
```

For PEEC-MSC:
```
M_mp[k, i] = B_pm[i, k]  (flux linkage reciprocity)
```

### Numerical Asymmetry

**Causes**:
1. **Different evaluation points**:
   - PEEC: Loop filament path (line integral)
   - MSC: Element center or face midpoints (point evaluation)

2. **Discretization error**:
   - Neumann formula uses loop contours
   - Biot-Savart uses element positions
   - Non-matching quadrature rules

3. **MSC EIEM2 evaluation points**:
   - Face evaluation point: `0.5 * (face_center + element_center)`
   - NOT at face center → breaks perfect symmetry

**Typical asymmetry**: 1-5% difference in `|M_mp - B_pm^T|`

### Symmetrization Strategies

**Option 1: Average symmetrization**
```python
K_coupling = 0.5 * (B_pm + M_mp.T)
```
- Pros: Enforces exact symmetry
- Cons: Loses physical meaning if difference is large

**Option 2: Use B_pm only (Galerkin)**
```python
K_mp = B_pm.T
```
- Pros: Mathematically consistent (Galerkin method)
- Cons: Ignores actual flux computation

**Option 3: No symmetrization (current approach)**
```python
# Use computed matrices as-is
# Solver: LU or BiCGSTAB (handle non-symmetric)
```
- Pros: Physically accurate
- Cons: Cannot use symmetric solvers (Cholesky, CG)

**Recommendation**:
- **Validation phase**: Use Option 3 (no symmetrization) to verify physics
- **Production**: Use Option 1 or 2 if symmetric solver is needed (future optimization)

### Energy Conservation Test

**Check**: Total energy should be conserved in lossless case.

```python
# Input power from source
P_source = Re(V_source^H · I)

# Stored energy in PEEC
W_peec = 0.5 * Re(I^H · Z_peec · I)

# Stored energy in MSC
W_msc = 0.5 * Re(σ^H · K_msc · σ)

# Energy balance
P_source ≈ jω (W_peec + W_msc)  (for lossless case)
```

Asymmetry in coupling matrix affects energy balance accuracy.

---

## Implementation Roadmap

### Phase 1: PEEC Single Conductor (No MSC)

**Goal**: Verify PEEC solver independently.

**Test cases**:
1. Circular loop self-inductance (analytical: Neumann formula)
2. Rectangular conductor DC resistance (analytical: R = ρL/A)
3. Skin effect validation (compare to Dowell's formula)

**Deliverables**:
- `rad.PEECLoop()` API for simple loop
- `rad.PEECFromMesh()` API for general mesh
- Unit tests with analytical solutions

**Estimated effort**: 1-2 weeks

---

### Phase 2: MSC Single Element (No PEEC)

**Goal**: Verify MSC solver with external field.

**Test cases**:
1. Single hex in uniform field (compare to analytical M = χ·H)
2. Two-element mutual interaction (compare to analytical Neumann)
3. Nonlinear material (compare to CplMag or ELF/MAGIC)

**Deliverables**:
- External field application: `rad.ObjBckg(callback)`
- MSC solver with Newton-Raphson (already implemented)
- Validation examples

**Status**: ✅ Already implemented (MSC + HACApK complete)

---

### Phase 3: PEEC-MSC Coupling (Linear)

**Goal**: Couple PEEC and MSC for linear materials (constant μ_r).

**Test cases**:
1. Loop + single hex core (compare to CplMag)
2. Loop + multi-element core (compare to CplMag container)
3. Frequency sweep (verify impedance vs frequency)

**Implementation**:
```python
# Build coupling matrices
B_pm = build_peec_to_msc_coupling(peec_mesh, msc_mesh)
M_mp = build_msc_to_peec_coupling(msc_mesh, peec_mesh)

# Assemble coupled system
Z_coupled = [[Z_peec,  -1j*omega*M_mp],
             [B_pm,     K_msc          ]]

# Solve
solution = solve(Z_coupled, rhs)
I = solution[:n_peec]
sigma = solution[n_peec:]
```

**Deliverables**:
- Coupling matrix builders
- Coupled solver
- Comparison with CplMag

**Estimated effort**: 2-3 weeks

---

### Phase 4: PEEC-MSC Coupling (Nonlinear)

**Goal**: Extend to nonlinear magnetic materials (μ_r(H)).

**Test cases**:
1. Loop + nonlinear core (B-H curve)
2. Core saturation (verify M vs H_ext)
3. Convergence with Newton damping

**Implementation**:
```python
# Newton-Raphson iteration
for iter in range(max_iter):
    # Compute external field
    H_ext = B_pm @ I

    # Update MSC Jacobian from μ_r(H)
    K_msc, dK_dH = compute_msc_jacobian(sigma, H_ext, BH_curve)

    # Solve coupled system
    [[Z_peec,  -1j*omega*M_mp],  [[ΔI  ],   [[R_peec],
     [B_pm,     K_msc          ]]  [Δσ  ]] = -[[R_msc ]]

    # Line search damping
    alpha = line_search(...)
    I += alpha * ΔI
    sigma += alpha * Δσ

    # Check convergence
    if norm(Δσ) / norm(σ) < tol:
        break
```

**Deliverables**:
- Nonlinear coupled solver
- Newton-Raphson with damping
- Convergence diagnostics

**Estimated effort**: 2-3 weeks

---

### Phase 5: Symmetrization and Optimization

**Goal**: Improve performance and energy conservation.

**Tasks**:
1. Implement symmetrization options (average, Galerkin)
2. Energy conservation verification
3. Compare symmetric vs non-symmetric solvers

**Optional**:
- Cholesky solver for symmetric case
- Preconditioner tuning for BiCGSTAB

**Estimated effort**: 1 week

---

## Summary

**PEEC-MSC coupling** extends Radia's capabilities to:
- Arbitrary conductor geometries (not just loops)
- MSC hexahedral elements (6-DOF, higher accuracy than MMM)
- Mixed element meshes (hex + tet)
- Unified PEEC+MSC framework

**Development priority**:
1. ✅ MSC solver (complete)
2. 🔄 PEEC solver (next)
3. ⏳ Coupling (after PEEC validation)
4. ⏳ Symmetrization (optimization phase)

**Key design decisions**:
- Use non-symmetric matrices initially (verify physics first)
- Support both linear and nonlinear materials
- Generalize CplMag to arbitrary meshes
- Maintain backward compatibility

---

**References**:
1. Ruehli, A.E. (1974). "Equivalent Circuit Models for Three-Dimensional Multiconductor Systems." IEEE Trans. MTT.
2. Hollaus, K. et al. (2024). "A Nonlinear Effective Surface Impedance in a Magnetic Scalar Potential Formulation." IEEE Trans. Magnetics.
3. Radia Documentation: `PEEC_SURFACE_IMPEDANCE.md`, `NPORT_BLOCK_LANCZOS_SPICE.md`
