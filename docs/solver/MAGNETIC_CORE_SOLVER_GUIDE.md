# Magnetic Core Solver Selection Guide

## Overview

When coupling PEEC conductors with magnetic cores, several solver options
are available for the core's magnetic and eddy current response.
The choice depends on material properties, frequency range, and accuracy needs.

## Decision Matrix

| Core Property | Recommended Solver | `core_model` | Module |
|---------------|-------------------|-------------|--------|
| Linear, non-conducting (ferrite) | **Radia MSC** | `'radia'` | `ngsbem_coupled.py` |
| Nonlinear, non-conducting (iron, saturation) | **Radia MSC** | `'radia'` | `ngsbem_coupled.py` |
| Linear, conducting, mu_r=1 (Al, Cu shield) | **Scalar FEM-BEM** | `'fembem'` | `ngsbem_coupled.py` |
| Linear, conducting, mu_r>1 (steel) | **Vector FEM-BEM** | `'vector_fembem'` | `ngsbem_coupled.py` |
| Any conducting, bounded domain | **FEM Dirichlet** | `'fem'` | `ngsbem_coupled.py` |
| Static, no eddy currents | **None (static mu_r)** | `None` | `ngsbem_coupled.py` |

## Solver Details

### Radia MSC (`core_model='radia'`)

**Best for**: Magnetic materials without eddy currents (ferrite, laminated iron)

- **Method**: Magnetic Surface Charge (MSC) integral equation
- **DOF**: 6 per hexahedron (surface charge per face), 3 per tetrahedron (Mx,My,Mz)
- **Domain**: Unbounded (no air mesh)
- **Materials**: Linear or nonlinear (B-H curve, hysteresis via Play model)
- **Frequency**: Static Delta_L, frequency-independent for non-conducting cores
- **Acceleration**: HACApK (H-matrix with ACA+) for large problems
- **Standalone**: Also available via `MMMBuilder` / `MMMSolver` (mmm_core.pyd) for
  Radia-independent operation with Schur complement coupling (`peec_msc_schur.py`)

**Delta_L computation**:
```python
# Via ngsbem_coupled.py
coupled = CoupledPEECMMM(peec_solver, core_model='radia', radia_core=core_handle)
coupled.compute_coupling_radia(solver_method=2)  # 2 = HACApK

# Via peec_coupled.py (simpler, column-by-column Radia Solve)
solver = CoupledPEECSolver(topo, [core_handle], mu_r_imag=0)
solver.compute_coupling_matrix()

# Via peec_msc_schur.py (standalone, H-matrix preserving)
schur = SchurComplementSolver()
schur.set_msc_system(N_matrix, dof_offset, inv_chi)
schur.solve(freq, V_source)
```

**Validation** (2026-03-28):
- 1 hex (6 DOF): MMMBuilder/Radia ratio = 1.000002
- 27 hex (162 DOF): ratio = 0.999685 (0.03% error)

### Scalar FEM-BEM (`core_model='fembem'`)

**Best for**: Non-magnetic conductors (Al shields, Cu enclosures)

- **Method**: Calderon projector (Hz scalar formulation)
- **Restriction**: mu_r = 1 ONLY (scalar Hz cannot handle magnetic contrast)
- **Domain**: Unbounded (BEM exterior)
- **Materials**: Linear, conducting (sigma > 0)
- **Frequency**: Full frequency-dependent eddy currents
- **Module**: `ngbem_eddy.py` → `EddyCurrentFEMBEM`

### Vector FEM-BEM (`core_model='vector_fembem'`)

**Best for**: Magnetic conducting cores (steel, silicon steel)

- **Method**: Johnson-Nedelec A-formulation + Weggler stabilized BEM
- **DOF**: H(curl) volume + HDivSurface boundary
- **Domain**: Unbounded
- **Materials**: Any linear mu_r, conducting
- **Frequency**: Full eddy currents with arbitrary permeability
- **Module**: `ngbem_eddy.py` → `VectorEddyCurrentFEMBEM`

### FEM Dirichlet (`core_model='fem'`)

**Best for**: Quick validation, bounded domain problems

- **Method**: FEM with H_inc Dirichlet BC on boundary
- **Domain**: Bounded (must truncate domain)
- **Materials**: Any linear mu_r, conducting
- **Frequency**: Eddy currents, bounded accuracy
- **Module**: Direct NGSolve H1/H(curl)

## Comparison: Radia MSC vs NGSBEM for Magnetic Cores

| Aspect | Radia MSC | NGSBEM (vector FEM-BEM) |
|--------|-----------|------------------------|
| **Eddy currents** | No | Yes |
| **Nonlinear** | Yes (B-H, hysteresis) | No (linear only) |
| **Domain** | Unbounded | Unbounded |
| **Mesh type** | Volume (hex/tet/wedge) | Volume (tet) + surface BEM |
| **Air mesh** | Not needed | Not needed (BEM exterior) |
| **Acceleration** | H-matrix (HACApK) | FMM or dense |
| **Typical DOF** | 6/hex, 3/tet | ~100-1000 per element |
| **Best regime** | DC / low freq / nonlinear | AC / eddy current dominated |

## When to Use Which

1. **Ferrite core (high mu_r, sigma ≈ 0)**: Use `'radia'`. Eddy currents negligible,
   MSC gives exact unbounded solution.

2. **Laminated steel (high mu_r, low effective sigma)**: Use `'radia'` with
   complex mu (mu_r_imag for loss). Lamination suppresses eddy currents.

3. **Solid steel (high mu_r, high sigma)**: Use `'vector_fembem'`.
   Eddy currents and magnetic permeability both matter.

4. **Aluminum shield (mu_r=1, high sigma)**: Use `'fembem'`.
   Fastest, scalar formulation sufficient.

5. **Nonlinear core with saturation**: Use `'radia'`. Only solver supporting
   B-H curves and hysteresis (Play model).

6. **Frequency sweep with conducting core**: Use `'vector_fembem'` or `'fem'`.
   Radia MSC gives static Delta_L only (no frequency-dependent eddy currents).
