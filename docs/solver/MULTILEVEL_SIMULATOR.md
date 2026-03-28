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

**Modules**: NGSolve (external), `esim_coupled_solver.py`, `mmm_ngsolve.py`

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
| MMMBuilder vs Radia | MSC kernel correctness | ±0.03% |
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

## Implementation Status (2026-03-28)

| Component | Status | Validated |
|-----------|--------|-----------|
| PEEC (filament) | Working | L within 0.1% of Grover |
| PEEC + MSC coupling | Working | Delta_L within 0.03% of Radia |
| NGSBEM (Laplace SL) | Working | Loop L within 5% of PEEC |
| NGSBEM + eddy (scalar) | Working | Tested on Al cube |
| NGSBEM + eddy (vector) | Working | Tested on steel cube |
| MMMBuilder (standalone MSC) | Working | 27-element validated |
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
├── mmm_core.pyd             # Level 1-2: Standalone MSC (C++)
├── ngsbem_peec.py           # Level 2: NGSBEM PEEC solver
├── ngsbem_coupled.py        # Level 2: NGSBEM + core coupling
├── ngsbem_eddy.py           # Level 2: Eddy current FEM-BEM
├── bem_inductance.py        # Level 2: Direct BEM inductance
├── esim_coupled_solver.py   # Level 3: FEM coupled solver
└── mmm_ngsolve.py           # Level 1-3: NGSolve mesh integration

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
