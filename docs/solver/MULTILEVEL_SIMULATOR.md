# Multi-Level Simulator Architecture

## Philosophy

Radia provides three levels of electromagnetic analysis within a single
repository. Each level trades speed for accuracy, and all share the same
geometry, coordinate system, and coupling interface — enabling seamless
cross-validation.

```
Level 1: PEEC          →  seconds   →  design exploration
Level 2: NGSBEM (BEM)  →  minutes   →  detailed analysis
Level 3: FEM           →  hours     →  final verification
```

## The Three Levels

### Level 1: PEEC (Partial Element Equivalent Circuit)

**Speed**: Sub-second to seconds
**Accuracy**: ±5% (filament approximation)
**Best for**: Design exploration, parametric sweeps, optimization loops

- Conductor modeled as 1D filaments with GMD self-inductance
- Neumann integral for mutual inductance
- Loop-Star decomposition for frequency response
- Magnetic core via Delta_L coupling (Radia MSC or analytical)
- Complex permeability for magnetic loss

**Modules**: `peec_matrices.pyd`, `peec_topology.py`, `peec_coupled.py`

**Typical workflow**:
```python
from fasthenry_parser import FastHenryParser
parser = FastHenryParser()
parser.parse_string(inp)           # FastHenry .inp format
result = parser.solve()            # L, R, Z(f) in one call
```

**Strengths**:
- 100+ frequency points in seconds
- Natural circuit representation (Z-parameters, S-parameters)
- Easy parametric studies (sweep core position, mu_r, geometry)

**Limitations**:
- No current redistribution (skin/proximity effect in conductor)
- Filament approximation limits accuracy for closely spaced conductors
- Static Delta_L (no eddy currents in core without FEM coupling)

### Level 2: NGSBEM (Boundary Element Method)

**Speed**: Seconds to minutes
**Accuracy**: ±1% (surface current distribution)
**Best for**: Detailed frequency response, validation of PEEC results

- Surface current on conductor boundary (HDivSurface, RWG basis)
- Laplace single-layer operator for inductance
- Helmholtz kernel for full-wave (future)
- Unbounded domain (no air mesh, no PML)
- Natural open boundary conditions

**Modules**: `ngsbem_peec.py`, `bem_inductance.py`, `ngsbem_coupled.py`

**Eddy current options**:
- Scalar FEM-BEM (`EddyCurrentFEMBEM`): mu_r=1 conductors (Al, Cu shields)
- Vector FEM-BEM (`VectorEddyCurrentFEMBEM`): any mu_r (steel cores)

**Typical workflow**:
```python
from ngsbem_peec import NGBEMPEECSolver
solver = NGBEMPEECSolver(mesh, order=0, sigma=5.8e7)
solver.assemble(intorder=6)
Z = solver.solve_frequency(1e6)
```

**Strengths**:
- Captures skin effect and proximity effect naturally
- Accurate mutual inductance (surface current, not filament)
- High-order elements available (order=1,2)
- Unbounded domain (exact for open structures)

**Limitations**:
- Dense BEM matrices (O(N^2) storage, O(N^3) solve)
- Assembly time dominated by singular integration
- Memory limit: ~10,000 surface DOFs for direct solve

### Level 3: FEM (Finite Element Method)

**Speed**: Minutes to hours
**Accuracy**: Reference solution
**Best for**: Final verification, nonlinear materials, complex geometry

- Volume discretization with air mesh
- A-formulation (H(curl)) or T-Omega formulation
- Nonlinear B-H curves, hysteresis (Play model)
- Kelvin transformation for open boundaries
- Adaptive mesh refinement

**Modules**: NGSolve (external), `esim_coupled_solver.py`

**Typical workflow**:
```python
from ngsolve import *
# ... standard NGSolve A-formulation or T-Omega
```

**Strengths**:
- Handles any geometry, any material
- Nonlinear iteration (Newton-Raphson)
- Adaptive refinement for error control
- Well-established theory and software ecosystem

**Limitations**:
- Requires air mesh (volume ratio problem for thin structures)
- Open boundary needs special treatment (Kelvin, PML, ABC)
- Mesh generation can be the bottleneck
- Slowest of the three levels

## Cross-Validation Matrix

The key advantage of the multi-level architecture: **any two levels
can validate each other** on the same geometry.

| Comparison | What it validates | Expected agreement |
|------------|------------------|--------------------|
| PEEC vs NGSBEM | Conductor self-inductance | ±5% |
| PEEC vs Analytical | Filament formula accuracy | ±0.1% |
| NGSBEM vs FEM | BEM vs volume method | ±1% |
| PEEC+MSC vs NGSBEM+MSC | Full coupled system | ±5% |
| All vs FastHenry | External tool validation | ±5% |

## When to Use Each Level

### Design Phase (Level 1: PEEC)

```
Question: "How does core position affect inductance?"
→ PEEC: sweep 100 positions in 10 seconds
→ Select top 5 candidates for Level 2
```

### Verification Phase (Level 2: NGSBEM)

```
Question: "What is the accurate frequency response?"
→ NGSBEM: compute Z(f) at 20 frequencies in 5 minutes
→ Compare with PEEC to validate trends
→ Flag any PEEC outliers for Level 3
```

### Final Sign-off (Level 3: FEM)

```
Question: "Does the design meet spec under saturation?"
→ FEM: nonlinear solve at operating point
→ Compare with Level 2 for linear validation
→ Report with confidence bounds
```

## Magnetic Core Coupling Across Levels

| Level | Core Method | Eddy Currents | Nonlinear |
|-------|------------|---------------|-----------|
| 1. PEEC | Radia MSC (Delta_L) | No | Yes |
| 2. NGSBEM | Scalar/Vector FEM-BEM | Yes | No |
| 3. FEM | Volume FEM | Yes | Yes |

The Radia MSC method bridges Level 1 and Level 2:
- At Level 1: Direct Delta_L computation (fast, column-by-column)
- At Level 2: Schur complement coupling preserving H-matrix (peec_msc_schur.py)

### Core Solver Selection

```
Core conducting (sigma > 0)?
 No  → Radia MSC ('radia')         [ferrite, laminated steel, nonlinear]
 Yes → mu_r > 1?
        No  → Scalar FEM-BEM ('fembem')     [Al/Cu shield]
        Yes → Vector FEM-BEM ('vector_fembem') [solid steel]
```

| Core Type | `core_model` | Eddy | Nonlinear | Module |
|-----------|-------------|------|-----------|--------|
| Ferrite (sigma~0) | `'radia'` | No | Yes | `peec_coupled.py` |
| Laminated steel | `'radia'` + complex mu | No | Yes | `peec_coupled.py` |
| Solid steel (mu_r>1, sigma>0) | `'vector_fembem'` | Yes | No | `ngsbem_eddy.py` |
| Al/Cu shield (mu_r=1) | `'fembem'` | Yes | No | `ngsbem_eddy.py` |
| Any, bounded domain | `'fem'` | Yes | No | NGSolve direct |

### Delta_L Computation APIs

```python
# Level 1a: FastHenry format (simplest)
from fasthenry_parser import FastHenryParser
result = FastHenryParser().parse_string(inp).solve()

# Level 1b: Python API (peec_coupled.py)
solver = CoupledPEECSolver(topo, [core_handle])
solver.compute_coupling_matrix()

# Level 1c: Standalone MSC (peec_msc_schur.py, H-matrix preserving)
schur = SchurComplementSolver()
schur.set_msc_system(N_matrix, dof_offset, inv_chi)
schur.solve(freq, V_source)

# Level 2: NGSBEM + Radia core (ngsbem_coupled.py)
coupled = CoupledPEECMMM(peec_solver, core_model='radia', radia_core=core)
coupled.compute_coupling_radia()
```

## Quick Start: Wire + Ferrite Core (Level 1)

```python
import radia as rad
import numpy as np
from peec_matrices import PEECBuilder
from peec_coupled import CoupledPEECSolver

rad.UtiDelAll()

# Conductor
builder = PEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.1, 0, 0)          # 100mm wire
builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
builder.add_port(n1, n2)
topo = builder.build_topology()

# Magnetic core
core = rad.ObjHexahedron([
    [0.02,0.005,-0.005],[0.08,0.005,-0.005],
    [0.08,0.015,-0.005],[0.02,0.015,-0.005],
    [0.02,0.005,0.005],[0.08,0.005,0.005],
    [0.08,0.015,0.005],[0.02,0.015,0.005]], [0,0,0])
rad.MatApl(core, rad.MatLin(999))

# Solve
solver = CoupledPEECSolver(topo, [core])
solver.compute_coupling_matrix(mu_r_real=1000)
Z = solver.compute_port_impedance(1e6)
L = np.imag(Z) / (2*np.pi*1e6)
```

## Common Pitfalls

1. **Coordinates in meters**: `60mm` = `0.06`, not `60`
2. **`rad.UtiDelAll()` first**: Radia keeps global state
3. **Hex vertex order**: bottom CCW (v0-v3), top CCW (v4-v7)
4. **NGSBEM surface mesh**: Use `Glue(wire.faces)`, not volume Box
5. **NGSBEM maxh**: `maxh <= min_cross_section / 2` for equilateral elements
6. **MSC sign**: the current multipole-moment MMM system is assembled by `BuildMomentSystemCore`; do not reconstruct it from the retired EIEM2 collocation sign convention.
7. **No Yano eval point in production**: the old midpoint eval point belonged to EIEM2 and was deleted. Current multipole-moment MMM samples the applied field at the element centroid and uses centroid field/gradient moment rows.
8. **Center-charge correction is internal**: mutual face-center cancellation lives inside the moment assembly / `CentroidFieldGradFromFace`, not in user-level PEEC coupling code.
9. **Loop port**: Split loop with two nodes at same position, not `add_port(n,n)`
10. **Nonlinear: Newton->Picard order**: Start Newton (fast), finish Picard (stable).
    Newton can excite zero-eigenvalue modes -> wrong solution. Use `keep_magnetization`:
    ```python
    rad.SolverConfig(newton_method=True)
    rad.Solve(obj, 1e-3, 10, 2)   # Newton phase
    rad.SolverConfig(newton_method=False, keep_magnetization=True)
    rad.Solve(obj, 1e-3, 100, 2)  # Picard continues
    ```
11. **BEM EFIE-SIBC wrong for finite Z_s**: The EFIE `Z_s*J + jw*mu0*SL(J) = -jw*A_inc`
    has SL eigenvalue R/3 (not R) for l=1 on sphere. Gives factor-of-3 error in Z_s term.
    Only correct for PEC (Z_s -> 0). Use **FEM-SIBC** (`fem_esim_3d.py`) instead.
    Fix requires MFIE (not available in ngsolve.bem).
13. **BEM formulation selection for MQS-SIBC**: Not all BEM formulations work in MQS.
    - **EFIE**: Only valid when `Z_s/(jw*mu0*R) < 0.1` (copper at high freq). Fails for steel.
    - **MFIE tangential** (n x H): Gives PEC solution only (no Z_s dependence).
    - **MFIE normal** (Sugahara): Supports finite Z_s via Biot-Savart + surface divergence,
      but requires custom matrix assembly (not in ngsolve.bem).
    - **PMCHWT-SIBC**: Impossible in MQS (`jw*eps*SL(M)` is O((kR)^2) ~ 10^{-14}).
    - **FEM-SIBC**: Recommended for all finite Z_s problems (total-field formulation).
12. **Scattered-field SIBC RHS**: Must include BOTH `-(jw/Z_s)*<A_inc, v>` AND
    `-<n x H_inc, v>` on the SIBC boundary. Missing the second term causes factor-of-3
    error. Total-field formulation (`fem_esim_3d.py`) does not have this issue.

## Implementation Status (2026-03-28)

| Component | Status | Validated |
|-----------|--------|-----------|
| PEEC (filament) | Working | L within 0.1% of Grover |
| PEEC + MSC coupling | Working | Delta_L within 0.03% of Radia |
| NGSBEM (Laplace SL) | Working | Loop L within 5% of PEEC |
| NGSBEM + eddy (scalar) | Working | Tested on Al cube |
| NGSBEM + eddy (vector) | Working | Tested on steel cube |
| Schur complement solver | Working | Framework validated |
| FEM (NGSolve) | Working | Via esim_coupled_solver |
| H-matrix (HACApK) | Working | Integrated with Radia core |
| Cross-validation scripts | Available | 4 verification scripts |

## File Map

```
src/radia/
├── peec_matrices.pyd        # Level 1: PEEC matrix assembly (C++)
├── peec_topology.py         # Level 1: Circuit topology solver
├── peec_coupled.py          # Level 1: PEEC + Radia core coupling
├── peec_msc_schur.py        # Level 1-2: Schur complement bridge
├── fasthenry_parser.py      # Level 1: FastHenry format input
├── ngsbem_peec.py           # Level 2: NGSBEM PEEC solver
├── ngsbem_coupled.py        # Level 2: NGSBEM + core coupling
├── ngsbem_eddy.py           # Level 2: Eddy current FEM-BEM
├── bem_inductance.py        # Level 2: Direct BEM inductance
└── esim_coupled_solver.py   # Level 3: FEM coupled solver

examples/peec_integration/
├── coupled/
│   ├── demo_schur_complement.py       # Schur complement demo
│   └── verify_multielement_msc.py     # Multi-element MSC validation
├── ngsbem_peec_demo/
│   ├── verify_peec_vs_ngsbem.py       # PEEC vs NGSBEM cross-validation
│   ├── verify_loop_peec_vs_ngsbem.py  # Loop inductance comparison
│   └── compute_L_final.py            # NGSBEM reference computation
└── applications/
    └── demo_fasthenry_magnetic_core.py # FastHenry + core demo
```
