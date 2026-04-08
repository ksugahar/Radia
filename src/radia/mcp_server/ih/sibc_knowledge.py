"""
SIBC / ESIM knowledge for induction heating MCP server.

Covers: Surface Impedance Boundary Condition approaches for eddy current,
BEM (Scalar BIE) and FEM (scattered-field) formulations, Karl iteration,
ESIM nonlinear impedance, and validation results.

Updated: 2026-04-09 with FEM-SIBC BND trace investigation results.
"""

IH_SIBC_OVERVIEW = """
# SIBC for Induction Heating: Method Selection

## Two Validated Approaches

| Method | Formulation | Validated | Error | DOFs |
|--------|-------------|-----------|-------|------|
| **Scalar BIE + SIBC (BEM)** | H1 surface phi | pytest 4/4 | **<0.1% sphere** | 162-320 |
| **FEM scattered-field + hole** | HCurl A_scat | verify_sphere | **-2.7% sphere** | ~150k |

## Scalar BIE + SIBC (Recommended for P_total)

**File**: `src/radia/bem_sibc_solver.py` (ScalarBIESIBCSolver)

System: `(1/2*M - DL + gamma * SL * M^{-1} * K) phi = phi_inc`
- gamma = Z_s / (jw * mu0)
- M = H1 surface mass, K = H1 stiffness (Laplace-Beltrami)
- DL, SL = Laplace double/single layer (ngsolve.bem)
- Gauge: Lagrange multiplier for int(phi) dS = 0
- Unknown: phi = exterior scalar magnetic potential (H = -grad phi)
- H_t = |grad_s phi| = tangential H = surface current
- P = 0.5 * Re(Z_s) * H_t^2 * area
- **No C++ needed**: all existing ngsolve.bem operators

**Validated**:
- Sphere: <0.1% for ALL Z_s (PEC to transparent), pytest 4/4
- Cylinder workpiece + coil: ~7% mesh-dependent

**Usage**:
```python
from radia.bem_sibc_solver import ScalarBIESIBCSolver, compute_phi_inc_from_loop
solver = ScalarBIESIBCSolver(mesh_wp, order=1)
result = solver.solve(phi_inc, Z_s=Z_s, omega=omega)
P = result['P_density'] * result['area']
H_t = result['H_t_rms']
```

**Panel script**: `src/radia/panels/calc_heating_bem.py`
- Biot-Savart filament coil (vectorized, 117x faster than naive)
- --impedance-model esim (nonlinear BH) or linear (fixed mu_r)
- Karl iteration for ESIM convergence

## FEM Scattered-Field + Hole + Kelvin (For L, B, P)

**File**: `examples/cubit_panels/inductance/verify_sphere_sibc.py`

Decomposition: A = A_inc + A_scat
- A_inc: known analytically (uniform field, or Biot-Savart filament)
- Solve for A_scat only (HCurl, smaller perturbation)
- Hole approach: workpiece volume NOT meshed, wp_surface is air-side boundary
- **Kelvin transform REQUIRED** for correct screening (truncation gives 435% error)
- H_t from A_total = A_inc + A_scat on BND (A_inc dominates, correctly evaluated)

**Validated**: sphere R=10mm, steel 7kHz (ratio=3.0):
- FEM + Kelvin: H_t = 2.53 A/m (-1.7% vs analytical 2.58)
- FEM + Dirichlet (R_air=80mm): H_t = 13.8 A/m (+435%) ← WRONG

**RHS must include BOTH terms**:
```python
f(v) = -(jw/Z_s) * <A_inc, v>_sibc   # SIBC term
     + (-1)      * <n x H_inc, v>_sibc # incident field boundary term
```
Missing the second term causes factor-of-3 error.

## CRITICAL: What Does NOT Work (2026-04-09)

### FEM Total-Field + Internal Interface: BND Integral Returns ~0

`calc_fem_kelvin.py` and `fem_esim_3d.py` used total-field formulation with
workpiece as a separate volume (internal interface). This approach has a
**fundamental problem**:

```python
# This returns ~0 on internal interface (workpiece side trace):
At_sq = sum(gfu[i].real**2 + gfu[i].imag**2 for i in range(3))
int_At2 = Integrate(At_sq, mesh, BND, definedon=wp_region)  # ~1e-17
```

**Root cause** (verified 2026-04-09):
1. NGSolve evaluates BND integrals from the **workpiece side** of the interface
2. Workpiece volume has no source; gauge regularization suppresses A to ~0
3. Air-side |A_t| ~ 3.0, workpiece-side |A_t| ~ 2e-6 (ratio 1.5 million)
4. `gfu.Other()` is DG-only, not available for conforming HCurl
5. Robin BC IS assembled correctly (B field changes 2.5% with/without Robin)
6. But P_total extraction from `int|A_t|^2` is impossible

**Previous "validation"** (commit 3ef7555, 2026-03-28) compared FEM-SIBC vs
EFIE-SIBC (BEM). But EFIE-SIBC has -65% error on sphere (SL eigenvalue mismatch).
"10% agreement with a 65%-wrong method" is NOT validation.

### EFIE-SIBC (HDivSurface): Wrong for Finite Z_s

RT0 (order=0) has zero surface curl. SIBC requires curl_s(J) which vanishes
for RT0 basis. Results are wrong by 65% on sphere at Z_s/(jw*mu0*R) = 10.

### Hole + Total-Field: PEC

Without scattered-field decomposition, hole boundary has natural Neumann BC
(n x H = 0 = PEC). Robin penalty only strengthens PEC tendency.

## Summary: When to Use What

| Need | Method | Script |
|------|--------|--------|
| **P_total, H_t** | Scalar BIE (BEM) | calc_heating_bem.py |
| **L (inductance)** | FEM (Kelvin) | calc_fem_kelvin.py |
| **B distribution** | FEM (volume integral) | calc_fem_kelvin.py |
| **Coil optimization** | BEM for P, FEM for L | both |

BEM is fast for P_total (162 DOFs vs 110k). FEM is needed for field distribution
and inductance. Both are needed for coil optimization.
"""

IH_ESIM = """
# ESIM: Effective Surface Impedance Method

ESIM extends linear SIBC to nonlinear magnetic materials by solving a 1D
cell problem through the conductor depth.

## Linear SIBC (Baseline)
```python
Z_s = (1+j) * rho / delta
delta = sqrt(2*rho / (omega * mu0 * mu_r))
```
Fixed Z_s, no iteration. Fast but inaccurate for steel (mu depends on H).

## ESIM Cell Problem
Solves 1D BVP: `rho * d^2H/dz^2 + jw * mu(|H|) * H = 0`
- Boundary: H(0) = H_t (surface field), dH/dz(d) = 0 (center)
- Returns: Z_s(H_t) = E_t(0) / H_t(0) (surface impedance)
- Handles: nonlinear BH curve, complex mu, finite thickness, cylinder geometry

## Karl Iteration
```
1. Initial Z_s from ESIM at estimated H_t
2. Solve BEM/FEM with Z_s -> get H_t on workpiece surface
3. Update Z_s from ESIM cell problem at new H_t
4. Relaxation: Z_s = 0.5*Z_s_new + 0.5*Z_s_old
5. Converge when dZ/Z < 1e-3 (typically 4-6 iterations)
```

## Module
```python
from radia.esim_cell_problem import ESIMFiniteSlabSolver
esim = ESIMFiniteSlabSolver(half_thickness=R_wp, bh_curve=BH_DATA,
                            sigma=sigma, frequency=freq, geometry='cylinder')
sol = esim.solve(H_t_rms)
Z_s = sol['Z']        # Complex surface impedance
P_prime = sol['P_prime']  # Power density [W/m^2]
```
"""

IH_BIOT_SAVART = """
# Biot-Savart Coil Field for BEM-SIBC

## phi_inc Computation (Scalar Magnetic Potential)

BEM-SIBC needs phi_inc (incident scalar potential) on workpiece surface.
For filamentary coil:

```python
from radia.bem_sibc_solver import compute_phi_inc_from_loop
phi_inc = compute_phi_inc_from_loop(
    node_coords, loop_center=[0,0,0], loop_radius=R_coil,
    current=I, n_quad=30, gap_deg=5)
```

Algorithm: z-axis analytical + horizontal path integration with Biot-Savart.
Vectorized: h_segments_batch (117x faster than naive Python loop).

## A_inc Computation (Vector Potential, for FEM Scattered-Field)

```python
from radia.biot_savart import a_segments_batch
A_inc_nodes = a_segments_batch(coil_segs, node_coords, current=1.0)
```

Vectorized over observation points, loops over segments. 106x faster.
For Joachim: requesting BiotSavartCF with A-field output would eliminate
the need for pre-computed nodal values.

## H_inc Computation (for FEM RHS nxH_inc term)

```python
from radia.biot_savart import h_segments_batch
H_inc_nodes = h_segments_batch(coil_segs, node_coords, current=1.0)
```
"""

IH_SCREENING = """
# Screening Physics

## Dimensionless Parameter

The key dimensionless parameter is `Z_s / (jw * mu0 * a)` where `a` is the
workpiece characteristic size (radius for cylinder).

| Z_s / (jw*mu0*a) | Behavior | One-way accuracy | Example |
|-------------------|----------|-----------------|---------|
| < 0.3 | Weak screening | One-way OK (-11%) | Copper 1kHz |
| 0.3 - 3 | Transition | One-way unreliable | Steel 7kHz (ratio=3.0) |
| > 3 | Strong screening | **One-way fails (100x+ error)** | Steel at high freq |

One-way models use H_t = H_inc (PEC approximation).
For steel at 7kHz: H_t = 0.77 A/m, not 18 A/m. One-way overestimates P by 300x.
Two-way (Karl iteration with ESIM) is essential for magnetic materials.

## FEM Scattered-Field: Kelvin is REQUIRED for Screening

**CRITICAL (2026-04-09)**: FEM scattered-field with Dirichlet truncation
(R_air = 8*R_wp) gives H_t ≈ H_inc (no screening). Kelvin transform is
REQUIRED for correct screening at Z_s/(jw*mu0*a) >= 1.

Verified:
- **With Kelvin**: verify_sphere_sibc.py, ratio=10, -2.7% error ← CORRECT
- **Without Kelvin**: IH hole mesh, ratio=3.0, H_t=13.8 vs BEM 6.3 ← WRONG (2.2x)

The Dirichlet truncation reflects the scattered field back, contaminating
the screening effect. The scattered field from SIBC is small but its
incorrect reflection dominates at high screening ratios.

For BEM (Scalar BIE): no truncation issue. BEM naturally handles open boundary.

## Typical IH Parameters

| Material | sigma [S/m] | mu_r | f [Hz] | delta [mm] | ratio |
|----------|-------------|------|--------|------------|-------|
| Steel | 2e6 | 100 | 7000 | 0.43 | 3.0 |
| Copper | 5.8e7 | 1 | 1000 | 6.6 | 0.01 |
| Copper | 5.8e7 | 1 | 100000 | 0.66 | 0.14 |
| Aluminum | 3.5e7 | 1 | 7000 | 1.01 | 0.04 |

Steel always needs two-way. Copper/aluminum OK with one-way at low frequencies.
"""


def get_ih_sibc_documentation(topic="all"):
    """Return IH SIBC documentation by topic."""
    topics = {
        "overview": IH_SIBC_OVERVIEW,
        "esim": IH_ESIM,
        "biot_savart": IH_BIOT_SAVART,
        "screening": IH_SCREENING,
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
