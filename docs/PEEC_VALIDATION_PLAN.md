# PEEC Solver Validation Plan

**Status**: Phases 1-3 Complete (2026-02-16)
**Purpose**: Systematic validation of PEEC solver and PEEC-MSC coupling

---

## Table of Contents

1. [Validation Philosophy](#validation-philosophy)
2. [Phase 1: PEEC Topology & MNA Solver](#phase-1-peec-topology--mna-solver)
3. [Phase 2: Multi-filament & Skin Effect](#phase-2-multi-filament--skin-effect)
4. [Phase 3: FastHenry Parser & Coupled PEEC+MMM](#phase-3-fasthenry-parser--coupled-peecmmm)
5. [Phase 4: Panel/Capacitance & Resonance](#phase-4-panelcapacitance--resonance)
6. [Phase 5: Complex Geometries](#phase-5-complex-geometries)
7. [Acceptance Criteria](#acceptance-criteria)

---

## Validation Philosophy

### Pyramid of Validation

```
         +---------------------+
         |  Complex Geometries |  Phase 5: Real applications
         |  (Transformers, WPT)|
         +---------------------+
               ^
               |
         +---------------------+
         | Panel + Resonance   |  Phase 4: Capacitance, LC circuits
         | FastHenry + Coupled |
         +---------------------+
               ^
               |
      +--------+--------+
      |                 |
+-------------+   +-------------+
| Multi-fil.  |   |  FastHenry  |  Phase 2-3: nwinc/nhinc, parser
| Skin effect |   |  Coupled    |
+-------------+   +-------------+
               ^
               |
         +---------------------+
         | Topology + MNA      |  Phase 1: Node-segment, multi-port
         | Series/Parallel/DC  |
         +---------------------+
```

### Validation Principles

1. **Bottom-up**: Validate simple cases before complex
2. **Analytical first**: Compare to closed-form solutions when available
3. **Cross-verification**: Compare C++ (LAPACK) vs Python (scipy) paths
4. **Convergence study**: Mesh refinement, solver tolerance
5. **Energy conservation**: Verify power balance

---

## Phase 1: PEEC Topology & MNA Solver

**Status**: COMPLETE (4/4 tests pass)
**Goal**: Verify node-segment topology and MNA multi-port solver

### Test 1.1: Series Wire (2 segments)

```python
from peec_matrices import PyPEECBuilder
from peec_topology import PEECCircuitSolver

builder = PyPEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.05, 0, 0)
n3 = builder.add_node_at(0.1, 0, 0)
builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
builder.add_connected_segment(n2, n3, 1e-3, 1e-3, sigma=5.8e7)
builder.add_port(n1, n3)
topo = builder.build_topology()

solver = PEECCircuitSolver(topo)
Z = solver.compute_port_impedance(freq=1e6)
# Verify: Z_series = Z1 + Z2 + 2*jw*M12
```

**Result**: 0.00% error vs legacy `create_wire` method.

### Test 1.2: Parallel Wires

```python
# Two segments sharing same nodes -> parallel circuit
builder = PyPEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.1, 0, 0)
builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
n3 = builder.add_node_at(0, 0.005, 0)
n4 = builder.add_node_at(0.1, 0.005, 0)
builder.add_connected_segment(n3, n4, 1e-3, 1e-3, sigma=5.8e7)
# Use .equiv to merge nodes...
```

**Result**: Z_parallel = Z_single/2 (0.00% error).

### Test 1.3: Series Analytical

Z_series = sum(R) + jw*(L11 + L22 + 2*M12)

**Result**: 0.00% error vs analytical formula.

### Test 1.4: DC Resistance

Series: R_total = R1 + R2
Parallel: R_total = R1*R2/(R1+R2)

**Result**: 0.00% error.

**Validation script**: `examples/peec_integration/validation/validate_topology.py`

---

## Phase 2: Multi-filament & Skin Effect

**Status**: COMPLETE (6/6 tests pass)
**Goal**: Verify nwinc/nhinc subdivision and skin effect modeling

### Test 2.1: Filament Count

Verify `nwinc=3, nhinc=3` creates 9 sub-filaments.

**Result**: PASSED.

### Test 2.2: DC Resistance Invariance

R_3x3 = R_1x1 (parallel reduction is exact).

**Result**: 0.00% error.

### Test 2.3: Inductance Reduction

L_3x3 < L_1x1 due to mutual coupling between sub-filaments.

**Result**: L_3x3 = 80.6 nH < L_1x1 = 82.1 nH (1.76% reduction). PASSED.

### Test 2.4: AC Resistance Ratio

R_ac/R_dc at 1 MHz for 3mm x 3mm Cu with 5x5 filaments.

**Result**: R_ac/R_dc = 2.85 (physically reasonable). PASSED.

### Test 2.5: Convergence

L_eff converges as nwinc/nhinc increases.

**Result**: Monotonic convergence confirmed. PASSED.

### Test 2.6: Series + Multi-filament

Combined topology with series chain + multi-filament segments.

**Result**: PASSED.

**Validation script**: `examples/peec_integration/validation/validate_multifilament.py`

---

## Phase 3: FastHenry Parser & Coupled PEEC+MMM

**Status**: COMPLETE (9/9 parser tests, coupling tests pass)

### FastHenry Parser Tests (9/9)

| Test | Description | Result |
|------|-------------|--------|
| Parser directives | .Units, N, E, .external, .freq | PASSED |
| .equiv | Node merge | PASSED |
| Single wire | Parsed vs manual builder (0.00% error) | PASSED |
| Parallel wires | R_parallel = R_single/2 | PASSED |
| Multi-filament | nwinc/nhinc from .inp matches manual | PASSED |
| Series chain | 4 segments in series | PASSED |
| Frequency sweep | .freq -> Z(f) sweep | PASSED |
| .default params | Inheritance of default parameters | PASSED |
| Continuation lines | Line continuation with '+' | PASSED |

**Validation script**: `examples/peec_integration/validation/validate_fasthenry.py`

### Coupled PEEC+MMM Tests

| Test | Description | Result |
|------|-------------|--------|
| Biot-Savart | Finite filament formula accuracy | PASSED |
| mu_r=1 | No coupling with air material | PASSED |
| High-mu | L increases with magnetic material | PASSED |
| Symmetry | Delta_L matrix symmetry | PASSED |
| Freq sweep | Z(f) physically reasonable | PASSED |
| .magnetic box | FastHenry box block parsing | PASSED |
| .magnetic hex | FastHenry hexahedron block parsing | PASSED |
| Coupled solve | Full parser -> coupled solve | PASSED |

**Validation script**: `examples/peec_integration/validation/validate_coupled.py`

---

## Phase 4: Panel/Capacitance & Resonance

**Status**: COMPLETE (23/23 tests pass)

### Coupling Tests (31/31)

2-port transformer coupling coefficient validation.

**Validation script**: `examples/peec_integration/validation/validate_coupling.py`

### Panel/Resonance Tests (23/23)

Panel potential coefficients and LC resonance validation.

**Validation script**: `examples/peec_integration/validation/validate_panel_resonance.py`

---

## Phase 5: Complex Geometries

**Status**: Future
**Goal**: Real-world applications

### Test 5.1: Planar Transformer

- Primary winding: PCB spiral (PEEC mesh)
- Secondary winding: PCB spiral (PEEC mesh)
- Core: Ferrite E-core (MSC hex mesh from Cubit)

### Test 5.2: Wireless Power Transfer (WPT)

- Transmitter coil: Litz wire (PEEC mesh)
- Receiver coil: Litz wire (PEEC mesh)
- Ferrite shields: MSC hex mesh
- Frequency: 6.78 MHz or 13.56 MHz

---

## C++ Solver Verification

### BiCGSTAB vs LU Numerical Equivalence

The templated BiCGSTAB solver (`rad_bicgstab.h`) produces results identical to LU:

| Test Case | Max Relative Difference |
|-----------|------------------------|
| Simple 2-segment wire | 3.40e-14% |
| 3x3 multi-filament (9 filaments) | 6.66e-13% |
| 2-port transformer (6 segments) | 6.07e-10% |
| Frequency sweep with panels | 7.42e-18% |

**All differences are at machine precision level** (< 1e-6%).

### Python scipy Removed

The PEEC MNA solver is **100% C++ LAPACK**. No Python scipy fallback exists.
All 73 validation tests pass with the C++ solver.

---

## Acceptance Criteria

### Error Tolerances

| Test Category | Acceptance Criterion |
|---------------|---------------------|
| **Topology analytical** | Error = 0.00% |
| **Multi-filament DC** | Error = 0.00% |
| **BiCGSTAB vs LU** | Error < 1e-6% |
| **FastHenry parser** | Exact match vs manual |
| **Coupled PEEC+MMM** | Physically reasonable |

### Convergence Requirements

| Solver | Max Iterations | Residual |
|--------|---------------|----------|
| **PEEC LU** | Direct (no iteration) | Exact |
| **PEEC BiCGSTAB** | < 1000 | < 1e-10 |
| **MSC Linear** | < 100 | < 1e-4 |
| **MSC Nonlinear** | < 50 | < 1e-3 |

### Test Summary

| Phase | Tests | Status |
|-------|-------|--------|
| 1: Topology & MNA | 4/4 | PASS |
| 2: Multi-filament | 6/6 | PASS |
| 3: FastHenry + Coupled | 9/9 + coupling | PASS |
| 4: Panel/Resonance + Coupling | 23/23 + 31/31 | PASS |
| 5: Complex geometries | - | Future |
| **Total** | **73/73** | **ALL PASS** |

---

**References**:
1. Grover, F.W. (2004). *Inductance Calculations*. Dover Publications.
2. Dowell, P.L. (1966). "Effects of eddy currents in transformer windings." Proc. IEE.
3. Zhu, Z., Song, B., White, J. (2005). "Algorithms in FastImp." IEEE Trans. TCAD.
4. Rosa, E.B. (1908). "The self and mutual inductances of linear conductors." NBS Bulletin.
