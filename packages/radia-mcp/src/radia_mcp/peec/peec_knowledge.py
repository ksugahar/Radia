"""

PUBLIC NOTEBOOK: docs/peec_integration/public_demo.ipynb -- human-facing PEEC
entry point.  It separates maintained representative demos from validation,
src-API, benchmark, experiment, smoke, and cleanup-review scripts in the
130-file docs/peec_integration/demos catalog.
SHOWCASE NOTEBOOK: docs/peec_integration/peec_showcase.ipynb -- Dowell continued-fraction + ngbem EFIE loop L + 4 paper figures (verified).
CATALOG NOTEBOOK: docs/peec_integration/examples_catalog.ipynb -- full source
text/hash/result-media ledger for docs/peec_integration/demos with migration lanes.
CLEANUP ROUTING: docs/peec_integration/cleanup_routing.ipynb -- routes the
63 cleanup-review scripts into validation_test (50), src API (5), and
distill-delete (8); deletion candidates are summarized in
memory/peec_integration_cleanup_routing.md.
VERIFICATION MIGRATION: docs/peec_integration/verification_migration.ipynb --
closed the PEEC verification corpus under
validation_test/peec_integration/verification, promoted the GMSH centerline
reader to radia.peec_mesh_import, and distilled/deleted check_funcs.py.
SHOWCASE NOTEBOOK: docs/peec/dowell_surface_impedance_demo.ipynb -- executable companion to docs/peec/PEEC_SURFACE_IMPEDANCE.md: F(xi)/G(xi) skin+proximity ratios vs normalized thickness, copper-foil Z(f) sweep, Dowell coth (H=0 BC) vs ESIM tanh (dH/dz=0 BC) boundary comparison, round-wire Bessel vs rect Dowell, and CF truncation accuracy (verified; F/G match deep_bar_resistance_factor/reactance_factor to machine precision).
PEEC knowledge base for the mcp-server-peec MCP server.

Covers: Loop-Star PEEC architecture, node-segment topology, circuit extraction,
multi-filament subdivision, SIBC surface impedance,
PRIMA model order reduction, SPICE netlist export, and ngsolve.bem integration.

Sources:
  - src/radia/peec_topology.py (PEECCircuitSolver, MNA)
  - src/radia/fasthenry_parser.py (FastHenry .inp parser)
  - src/radia/lanczos_reduction.py (PRIMA, Lanczos, SPICE export)
  - src/radia/esim_cell_problem.py (ESIM nonlinear Z_s)
  - src/core/rad_peec_matrices.h/cpp (C++ L, R, P, M_LS matrices)
  - A. E. Ruehli, IEEE MTT, 1974 (original PEEC formulation)
"""

PEEC_OVERVIEW = """
# PEEC (Partial Element Equivalent Circuit) Overview

## Architecture

PEEC extracts circuit parameters (L, R, C, M) from conductor geometry using
integral equations with the Laplace kernel G(r) = 1/(4*pi*r).

**Valid frequency range**: DC to ~100 MHz (Darwin/MQS approximation).

### Filament-Panel Decomposition (FastImp Style)

```
Surface Mesh -> Face -> Panel (Star: charge Q)  -> P matrix (capacitance)
             -> Edge -> Filament (Loop: current I) -> L, R matrices
```

Loop-Star basis transformation is NOT needed -- filaments and panels are
inherently separate in PEEC.

### System Equation

```
[Z_LL   Z_LS] [I_L]   [V_L]
[Z_SL   Z_SS] [I_S] = [V_S]

Z_LL = diag(R_dc + Z_s(f)) + jw*L    (loop impedance)
Z_SS = P / (jw)                        (panel admittance)
Z_LS = jw * M_LS                       (loop-star coupling)
Z_SL = Z_LS^T                          (reciprocity)
```

Where:
- L: mutual inductance matrix (Neumann formula, C++ `PEECMatrices`)
- R_dc: DC resistance (diagonal)
- Z_s(f): frequency-dependent surface impedance (Bessel/Dowell/ESIM)
- P: potential coefficient matrix (capacitance = P^{-1})
- M_LS: loop-star coupling matrix

**Loop-only mode** (no panels, backward compatible):
```
Z_branch = diag(R_dc + Z_s) + jw*L
Y_node = A * Z_branch^{-1} * A^T
Y_node * V = I_ext
```

### Matrix Assembly (C++)

Frequency-independent matrices (L, R_dc, P, M_LS) are computed in C++
(`rad_peec_matrices.cpp`). Frequency-dependent Z_s is computed in Python.

```python
from peec_matrices import PyPEECBuilder
builder = PyPEECBuilder()
# ... add segments ...
topo = builder.build_topology()
# topo contains: L, R, P, M_LS, segment_nodes, ports, etc.
```

### Benchmarks

SHOWCASE NOTEBOOK `docs/solver_benchmarks/peec_solver_benchmarks.ipynb`:
dense Ruehli L vs HACApK H-matrix (ACA+ compression + scaling) and the MNA
dense-vs-HACApK crossover, read live from the committed benchmark JSON; the
extracted |Z11| matches the dense reference. Corpus + JSON +
`findings_peec_mna_crossover.md` / `comparison_peec_dense_vs_hacapk.md` at
`docs/solver_benchmarks/`; executable benchmark drivers live under
`validation_test/solver_benchmarks/`.

### References

- A. E. Ruehli, "Equivalent Circuit Models...", IEEE MTT, 1974
- G. Vecchi, "Loop-Star decomposition...", IEEE TAP, 1999
- A. Odabasioglu et al., "PRIMA", IEEE TCAD, 1998
"""

PEEC_BUILDER_API = """
# PyPEECBuilder API (Node-Segment Topology)

## Basic Usage

```python
from peec_matrices import PyPEECBuilder

builder = PyPEECBuilder()

# Add nodes (returns node ID)
n1 = builder.add_node_at(0, 0, 0)           # x, y, z in meters
n2 = builder.add_node_at(0.05, 0, 0)
n3 = builder.add_node_at(0.1, 0, 0)

# Add segments between nodes
builder.add_connected_segment(n1, n2, 1e-3, 1e-3)  # width, height [m]
builder.add_connected_segment(n2, n3, 1e-3, 1e-3)

# Define port (positive node, negative node)
builder.add_port(n1, n3)

# Build topology -> returns dict with L, R, segment_nodes, ports, etc.
topo = builder.build_topology()
```

## add_connected_segment Parameters

```python
builder.add_connected_segment(
    node_from,       # Source node ID (from add_node_at)
    node_to,         # Destination node ID
    width,           # Cross-section width [m]
    height,          # Cross-section height [m]
    sigma=5.8e7,     # Conductivity [S/m] (default: copper)
    nwinc=1,         # Width subdivisions for multi-filament
    nhinc=1,         # Height subdivisions for multi-filament
)
```

## Multi-Port Extraction

```python
# 2-port transformer
builder.add_port(n1, n2)   # Port 1: primary
builder.add_port(n3, n4)   # Port 2: secondary

topo = builder.build_topology()
solver = PEECCircuitSolver(topo)

# Z-parameter matrix (2x2)
Z_mat = solver.compute_Z_matrix(freq=1e6)
# Z11 = V1/I1 (port 2 open), Z12 = V1/I2 (port 1 open), etc.

# Coupling coefficient
result = solver.compute_coupling_coefficient(freq=1e6)
k = result['k']       # k = M12 / sqrt(L1 * L2)
L1 = result['L1']     # Self-inductance port 1 [H]
L2 = result['L2']     # Self-inductance port 2 [H]
M = result['M']       # Mutual inductance [H]
```

## Topology Dict Contents

```python
topo = builder.build_topology()
# topo keys:
#   'L'              - (n_loop, n_loop) inductance matrix [H]
#   'R'              - (n_loop,) DC resistance array [Ohm]
#   'segment_nodes'  - (n_loop, 2) node connectivity [from, to]
#   'n_nodes'        - Number of nodes
#   'n_loop'         - Number of filament segments
#   'ports'          - List of (node_pos, node_neg, port_id)
#   'segment_centers'    - (n_loop, 3) segment centers [m]
#   'segment_directions' - (n_loop, 3) unit direction vectors
#   'segment_lengths'    - (n_loop,) segment lengths [m]
# Optional (when panels present):
#   'P'              - (n_star, n_star) potential coefficient matrix
#   'M_LS'           - (n_loop, n_star) loop-star coupling
#   'n_star'         - Number of panel (star) elements
#   'C_gathered'     - Gathered capacitance matrix
```
"""

PEEC_CIRCUIT_SOLVER = """
# PEECCircuitSolver (MNA Port Impedance)

## Basic Impedance Extraction

```python
from peec_topology import PEECCircuitSolver

solver = PEECCircuitSolver(topo)

# Single frequency
Z = solver.compute_port_impedance(freq=1e6)
L = Z.imag / (2 * np.pi * 1e6)   # Inductance [H]
R = Z.real                         # Resistance [Ohm]

# Frequency sweep
freqs = np.logspace(3, 9, 100)     # 1 kHz to 1 GHz
Z_sweep = solver.frequency_sweep(freqs)
```

## With Surface Impedance (Skin Effect)

```python
import numpy as np
from scipy.special import iv as besseli  # Modified Bessel I

mu0 = 4e-7 * np.pi

def circular_sibc(freq, radius, sigma, mu_r=1.0):
    \"\"\"Bessel SIBC for circular conductor.\"\"\"
    omega = 2 * np.pi * freq
    mu = mu0 * mu_r
    k = np.sqrt(1j * omega * mu * sigma)
    kr = k * radius
    Zs = (k / (2 * np.pi * radius * sigma)) * besseli(0, kr) / besseli(1, kr)
    return Zs

def Zs_func(freq):
    return np.array([circular_sibc(freq, 0.5e-3, 5.8e7)] * n_seg)

Z_sweep = solver.frequency_sweep(freqs, Zs_func=Zs_func)
```

## Solver Methods

```python
solver.set_solver_method(0)   # LU (LAPACK zgesv_, default)
solver.set_solver_method(1)   # BiCGSTAB (large systems)
solver.set_bicgstab_params(tol=1e-10, max_iter=1000)
```

## Resonance Detection (with panels)

```python
f_res = solver.get_resonant_frequencies(freq_min=1e3, freq_max=1e9)
# Returns list of frequencies where Im(Z) crosses zero
```
"""

PEEC_MULTI_FILAMENT = """
# Multi-Filament Subdivision (nwinc/nhinc)

## Purpose

Skin and proximity effects require current redistribution within
conductor cross-sections. Multi-filament splits each segment into
(nwinc * nhinc) sub-filaments with full mutual inductance coupling.

## Guidelines

| Regime | d/delta | nwinc x nhinc | Description |
|--------|---------|---------------|-------------|
| DC     | < 0.5   | 1 x 1         | Uniform current |
| Moderate | 2-5   | 3 x 3         | Onset of skin effect |
| Strong | > 5     | 5 x 5+        | Concentrated surface current |

Where d = conductor dimension, delta = skin depth.

## Usage

```python
builder = PyPEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.1, 0, 0)

# 3x3 = 9 sub-filaments per segment
builder.add_connected_segment(n1, n2, 3e-3, 3e-3,
                               sigma=5.8e7, nwinc=3, nhinc=3)
builder.add_port(n1, n2)
topo = builder.build_topology()
# topo['n_loop'] = 9 (9 filaments from 1 segment)
```

## Skin Depth Reference

```python
import numpy as np
mu0 = 4e-7 * np.pi
delta = np.sqrt(2 / (2 * np.pi * freq * mu0 * mu_r * sigma))
# Copper at 1 MHz: delta ~ 66 um
# Copper at 100 kHz: delta ~ 209 um
# Steel at 8 kHz: delta ~ 0.18 mm (mu_r=1000, sigma=10e6)
```
"""

PEEC_FASTHENRY = """
# FastHenry .inp Parser

## Parsing FastHenry Input Files

```python
from fasthenry_parser import FastHenryParser

parser = FastHenryParser()
parser.parse_file('inductor.inp')

# Get Radia PEEC builder
builder = parser.to_peec_builder()
topo = builder.build_topology()

# Get frequency sweep points
freqs = parser.get_frequencies()

# Solve conductor / shielded PEEC. .magnetic blocks are rejected by solve().
result = parser.solve(freqs)
```

## Supported Directives

| Directive | Format | Example |
|-----------|--------|---------|
| `.Units` | `.Units mm` | Length unit (m, mm, um, cm, in, mils) |
| `N<name>` | `N1 x=0 y=0 z=0` | Node definition |
| `E<name>` | `E1 N1 N2 w=1 h=1 nwinc=3 nhinc=3 sigma=5.8e7` | Segment |
| `.external` | `.external N1 N2` | Port definition |
| `.freq` | `.freq fmin=1e3 fmax=1e9 ndec=10` | Frequency sweep |
| `.default` | `.default sigma=5.8e7 w=0.5 h=0.1` | Default parameters |
| `.equiv` | `.equiv N1 N2 N3` | Node merge (equivalence) |
| `.magnetic` | Parsed, rejected by solve() | Magnetic material block |
| `.panel` | See below | Panel (capacitive) block |

## Magnetic Material Block (.magnetic)

The old PEEC magnetic-material coupling path is retired.  FastHenry
`.magnetic` blocks are still parsed so inputs can be diagnosed, but
`parser.solve()` raises and points you to HDiv-VIM / reduced FEM for
magnetic cores.

```
.magnetic
  type=box
  center=0.05,0.01,0.0
  size=0.06,0.01,0.01
  mu_r=1000
.endmagnetic
```

Supported types: `box` (hexahedron), `hexahedron` (8 vertices), `mesh` (file).

## FastHenry .inp File Format Example

```
* Spiral inductor on ferrite substrate
.Units mm
.default sigma=5.8e7 nwinc=3 nhinc=1

N1 x=0    y=0   z=0.5
N2 x=10   y=0   z=0.5
N3 x=10   y=10  z=0.5
N4 x=0    y=10  z=0.5

E1 N1 N2 w=0.5 h=0.1
E2 N2 N3 w=0.5 h=0.1
E3 N3 N4 w=0.5 h=0.1

.external N1 N4
.freq fmin=1e3 fmax=1e9 ndec=10

.magnetic
  type=box
  center=5,5,-0.5
  size=12,12,1
  mu_r=2000
.endmagnetic

.end
```
"""

PEEC_MAGNETIC_POLICY = """
# Magnetic Material Coupling Policy

## Overview

PEEC in Radia is conductor/shield circuit extraction.  Magnetic material
cores are no longer solved through the old PEEC magnetic-moment coupling
path.  Keep the magnetic state mesh-backed and use HDiv-VIM / reduced FEM,
then couple fields at the application layer.

## Usage

```python
from peec_matrices import PyPEECBuilder

# Build conductor topology
builder = PyPEECBuilder()
n1 = builder.add_node_at(0, 0, 0)
n2 = builder.add_node_at(0.1, 0, 0)
builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
builder.add_port(n1, n2)
topo = builder.build_topology()
```

For ferrite/iron cores:
1. solve the conductor PEEC problem or extract the current distribution;
2. solve the magnetic core with HDiv-VIM / reduced FEM;
3. exchange fields in the application workflow and validate against a
   result-bearing notebook or validation_test case.
"""

PEEC_SIBC = """
# SIBC (Surface Impedance Boundary Condition)

## Three SIBC Implementations

| Conductor Type | Method | Library | Use Case |
|----------------|--------|---------|----------|
| Circular wire | Bessel I0/I1 | `scipy.special.iv` | Round conductors |
| Rectangular (d << w) | Dowell formula | C++ `rad_peec_surface_impedance.cpp` | PCB traces, busbars |
| Nonlinear magnetic | ESIM cell problem | `esim_cell_problem.py` | Steel/ferrite workpieces |

## Bessel SIBC (Circular Conductor)

```python
import numpy as np
from scipy.special import iv as besseli  # Modified Bessel I_n

mu0 = 4e-7 * np.pi

def circular_sibc(freq, radius, sigma, mu_r=1.0):
    omega = 2 * np.pi * freq
    mu = mu0 * mu_r
    k = np.sqrt(1j * omega * mu * sigma)
    kr = k * radius
    Zs = (k / (2 * np.pi * radius * sigma)) * besseli(0, kr) / besseli(1, kr)
    return Zs

# CRITICAL: use scipy.special.iv (modified Bessel), NOT jv (regular Bessel)
```

## Dowell SIBC (Rectangular Conductor)

For thin rectangular conductors where thickness d << width w. The closed-form
Dowell factors live in Python:
`radia.em_material.EMMaterial.dowell_Zs(freq, R)` (tanh surface impedance) and
`radia_mcp.radia_ngsolve.solve.deep_bar_resistance_factor` /
`deep_bar_reactance_factor` (the single-conductor `F(xi)` / `G(xi)` ratios; the
C++ PEEC assembly in `rad_peec_matrices.cpp` consumes the per-segment `Z_s`).

SHOWCASE NOTEBOOK `docs/peec/dowell_surface_impedance_demo.ipynb` (companion to
`docs/peec/PEEC_SURFACE_IMPEDANCE.md`) runs these live: `F(xi)`/`G(xi)` vs
normalized thickness `xi=a/delta`, a copper-foil `Z(f)` sweep, the Dowell `coth`
(`H(a)=0`) vs ESIM `tanh` (`dH/dz(a)=0`) boundary comparison, a round-wire Bessel
cross-section comparison, and continued-fraction (1st/2nd-order) truncation
accuracy. (The full continued-fraction + PRIMA-ladder derivation is in
`docs/peec_integration/peec_showcase.ipynb`.)

## ESIM (Effective Surface Impedance Method)

For H-dependent (nonlinear) surface impedance on magnetic conductors:

```python
from radia.esim_cell_problem import ESIMCellProblem

# BH curve data: [[H (A/m), B (T)], ...]
BH_DATA = [[0, 0], [100, 0.1], [1000, 1.2], [50000, 2.0]]

esim = ESIMCellProblem(
    bh_curve=BH_DATA,
    sigma=10e6,        # Conductivity [S/m]
    freq=8000,         # Frequency [Hz]
)

# Solve cell problem for surface field H0
Z_s = esim.solve(H0=1000.0)
# Z_s = complex surface impedance [Ohm/m^2]
```

ESIM solves the 1D cell problem:
```
rho * d^2H/ds^2 + j*omega*mu(|H|)*H = 0   for s in [0, inf)
```

Supports complex permeability: mu = mu' - j*mu" (grain eddy current losses).

**Use for**: induction heating workpieces, nonlinear iron cores, lossy ferrite.
"""

PEEC_PRIMA = """
# PRIMA Model Order Reduction

## Overview

PRIMA (Passive Reduced-order Interconnect Macromodeling Algorithm) reduces
the full PEEC system to a small equivalent circuit suitable for SPICE.

Uses Lanczos tridiagonalization (NOT Arnoldi) -> RL ladder network directly.

### Reduction Pipeline

```
Full system: [Z_LL, Z_LS, Z_LM; ...] (N DOFs)
  -> Lanczos on L: L = Q @ T @ Q^T  (T tridiagonal, k << N)
  -> Transform: Z' = Q^T @ Z @ Q
  -> Schur complement: eliminate internal nodes
  -> SPICE netlist: RL ladder + mutual inductors
```

## Usage

```python
from radia.lanczos_reduction import (
    LanczosReducer, PRIMASchurExtractor, SPICEExtractionConfig
)

# Configure extraction
config = SPICEExtractionConfig(
    n_moments=20,       # Number of Lanczos vectors
    tol=1e-6,           # Truncation tolerance
)

# Reduce inductance matrix
reducer = LanczosReducer(tol=config.tol)
result = reducer.lanczos_symmetric(L, k=config.n_moments)
# result.T  = tridiagonal matrix (alpha, beta)
# result.Q  = orthonormal basis

# Schur complement for port extraction
extractor = PRIMASchurExtractor(topo, config)
circuit = extractor.extract()
# circuit: SparseCircuit with resistors, inductors, capacitors, ports
```

## SparseCircuit Structure

```python
circuit.resistors      # [(node1, node2, R), ...]
circuit.inductors      # [(node1, node2, L), ...]
circuit.capacitors     # [(node1, node2, C), ...]
circuit.mutual_inductors  # [(L1+, L1-, L2+, L2-, M), ...]
circuit.port_nodes     # [(pos, neg), ...] for each port
```

## SPICE Netlist Output

The SparseCircuit can be exported to SPICE format (LTspice, ngspice compatible):

```
* PEEC Reduced Model
R1 n1 n2 1.23e-3
L1 n2 n3 4.56e-9
K1 L1 L2 0.85
C1 n3 0 7.89e-12
.ends
```

## Key Classes

| Class | Purpose |
|-------|---------|
| `LanczosReducer` | Lanczos tridiagonalization with re-orthogonalization |
| `LanczosResult` | Result: Q (basis), alpha/beta (tridiagonal), T, rank |
| `PRIMASchurExtractor` | Port impedance extraction via Schur complement |
| `SPICEExtractionConfig` | Configuration (n_moments, tol) |
| `SparseCircuit` | Circuit netlist (R, L, C, M, ports) |
| Magnetic-material coupling | Retired from PEEC; use HDiv-VIM / reduced FEM |

## References

- Odabasioglu, Celik, Pileggi, "PRIMA", IEEE TCAD, 1998
- Gustavsen, Semlyen, "Vector Fitting", IEEE TPWRD, 1999
"""

PEEC_NGSOLVE_BEM = """
# ngsolve.bem Integration

## Frequency Range Overlap

| Range | Solver | Notes |
|-------|--------|-------|
| DC - 1 MHz | Radia PEEC + SIBC | Circuit extraction, SPICE netlist |
| DC - 1 MHz | ngsolve.bem (Weggler EFIE) | Low-freq stable EFIE |
| 1 MHz - GHz | ngsolve.bem (Helmholtz) | Full-wave BEM |

Both PEEC and ngsolve.bem cover the DC-1 MHz range. Choose based on output:
- **PEEC**: L, R, C, M extraction, SPICE netlist, model order reduction
- **ngsolve.bem**: field distribution, scattering, coupling to FEM

## Radia PEEC Unique Features (vs ngsolve.bem)

1. **Direct circuit extraction**: L, R, C, M matrices -> SPICE netlist
2. **Lanczos MOR**: PRIMA reduction to compact equivalent circuit
3. **FastHenry compatibility**: Parse existing .inp models
4. **Multi-filament**: nwinc/nhinc for skin/proximity within PEEC

## When to Use Which

| Task | Use |
|------|-----|
| Inductance/resistance extraction | PEEC |
| SPICE netlist for circuit simulator | PEEC + PRIMA |
| Inductor on ferrite substrate | PEEC for conductors + HDiv-VIM / reduced FEM for core |
| Field distribution visualization | ngsolve.bem |
| High-frequency scattering (> 1 MHz) | ngsolve.bem (Helmholtz) |
| Induction heating (P_total) | ngsolve.bem (Scalar BIE + SIBC) |
| Impedance of complex 3D geometry | Either (cross-validate) |
"""

PEEC_PITFALLS = """
# Common PEEC Pitfalls

## 1. Units: Always Meters

All coordinates in meters. Conductivity in S/m. Current density in A/m^2.
Do NOT convert units -- Radia always uses meters.

```python
# CORRECT
builder.add_node_at(0.001, 0, 0)  # 1 mm in meters
builder.add_connected_segment(n1, n2, 0.5e-3, 0.1e-3)  # 0.5mm x 0.1mm

# WRONG
builder.add_node_at(1, 0, 0)  # This is 1 meter, not 1 mm!
```

## 2. Bessel Function: iv (Modified), NOT jv

```python
from scipy.special import iv as besseli  # CORRECT: modified Bessel I
# from scipy.special import jv  # WRONG: regular Bessel J
```

MKL does not provide Bessel functions; use scipy.

## 3. Port Definition Required

Without `add_port()`, `compute_port_impedance()` returns 0.

## 4. Cleanup with UtiDelAll

PEEC no longer owns magnetic-material Radia objects.

## 5. Import Order (with Cubit)

If combining with Cubit panels, import order matters:
```python
import netgen       # First
import ngsolve      # Second
# ... setup cubit ...
import cubit        # Last (DLL conflict avoidance)
```

## 6. Surface Impedance Array Size

Zs array must have exactly n_loop elements (= number of filament segments).
With multi-filament (nwinc=3, nhinc=3), n_loop = 9 * n_segments.

## 7. Laplace Kernel Only

PEEC uses G(r) = 1/(4*pi*r). Do NOT add Helmholtz kernel (e^{-jkr}/r).
Wave propagation effects are NOT modeled. Use ngsolve.bem for high frequency.
"""


def get_peec_documentation(topic: str = "all") -> str:
    """Return PEEC documentation by topic."""
    topics = {
        "overview": PEEC_OVERVIEW,
        "builder": PEEC_BUILDER_API,
        "solver": PEEC_CIRCUIT_SOLVER,
        "multi_filament": PEEC_MULTI_FILAMENT,
        "fasthenry": PEEC_FASTHENRY,
        "magnetic": PEEC_MAGNETIC_POLICY,
        "sibc": PEEC_SIBC,
        "prima": PEEC_PRIMA,
        "ngsolve_bem": PEEC_NGSOLVE_BEM,
        "pitfalls": PEEC_PITFALLS,
    }

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    elif topic in topics:
        return topics[topic]
    else:
        return (
            f"Unknown topic: '{topic}'. "
            f"Available: all, {', '.join(topics.keys())}"
        )
