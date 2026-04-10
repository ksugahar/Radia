"""
SIBC / ESIM knowledge for induction heating MCP server.

Covers: Surface Impedance Boundary Condition approaches for eddy current,
BEM (Scalar BIE) and FEM (scattered-field) formulations, Karl iteration,
ESIM nonlinear impedance, and validation results.

Updated: 2026-04-09 with specialcf.normal sign fix, scattered-field cancellation,
total-field hole+BND approach.
"""

IH_SIBC_OVERVIEW = """
# SIBC for Induction Heating: Method Selection

## Three Validated Approaches

| Method | Formulation | Validated | Error | DOFs |
|--------|-------------|-----------|-------|------|
| **Scalar BIE + SIBC (BEM)** | H1 surface phi | pytest 4/4 | **<0.1% sphere** | 162-320 |
| **FEM scattered-field + hole** | HCurl A_scat | verify_sphere | **-2.7% sphere** | ~150k |
| **FEM total-field + hole + BND** | HCurl A_total | coil+cylinder | **+8.2% vs BEM** | ~150k |

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

## FEM Total-Field + BND Integral (Recommended for Coil+Workpiece)

**File**: `src/radia/panels/calc_fem_kelvin.py`

Total-field formulation (no A_inc decomposition):
- Coil as volume source J in coil cross-section
- Two mesh approaches: interface (workpiece as volume) or hole (workpiece removed)
- **Interface is more accurate** (-2.9% vs BEM) than hole (+5.9%)

**H_t extraction**: BND integral with tangential projection (specialcf.normal):
```python
n = specialcf.normal(3)
A_sq = sum(gfu[i].real**2 + gfu[i].imag**2 for i in range(3))
A_dot_n_re = sum(gfu[i].real * n[i] for i in range(3))
A_dot_n_im = sum(gfu[i].imag * n[i] for i in range(3))
At_sq = A_sq - A_dot_n_re**2 - A_dot_n_im**2
H_t_rms = abs(1j*omega/Z_s) * sqrt(Integrate(At_sq, mesh, BND, definedon=sibc) / A_wp)
```

**CRITICAL: Pointwise evaluation FAILS** near SIBC boundary. Robin penalty
|jw/Z_s| ~ 1e9 creates huge normal A component (~1e6 T*m). Tangential
projection at interior points requires 10^14 cancellation. Only BND integral
works (NGSolve evaluates trace functions directly).

**Interface is the correct formulation. Hole (PEC) is wrong.**

SIBC models a thin conducting shell: the interior is transparent (air), and
surface current J_s = (jw/Z_s)*A_t flows on the interface. This maps directly
to Robin BC on an **internal interface** (workpiece volume meshed as air).

Hole approach = PEC baseline (natural BC: n x H = 0). Adding Robin penalty to
a PEC boundary is not SIBC — it perturbs the wrong physics. The hole error
(+5.9%) is a formulation error, not a mesh resolution issue.

| Approach | Physics | H_t vs BEM | P vs BEM |
|----------|---------|-----------|---------|
| **Interface (correct)** | Transparent + Robin | **-2.9%** | **-2.3%** |
| Hole (wrong) | PEC + Robin perturbation | +5.9% | +16.2% |

Interface approach: workpiece volume meshed (as air, nu=nu0), SIBC on internal
boundary "wp_surface". The Robin BC `+jw/Z_s * u.Trace() * v.Trace() * ds`
acts as a thin conducting shell on the shared face.

## FEM Scattered-Field + Hole (For Sphere / Analytical A_inc)

**File**: `examples/cubit_panels/inductance/verify_sphere_sibc.py`

Decomposition: A = A_inc + A_scat
- A_inc: known analytically (uniform field, or Biot-Savart filament)
- Solve for A_scat only (HCurl, smaller perturbation)
- Hole approach: workpiece volume NOT meshed, wp_surface is air-side boundary

**RHS (correct signs)**:
```python
n_cond = -specialcf.normal(3)  # outward from conductor
nxH = Cross(H_inc_cf, specialcf.normal(3))  # = n_cond x H_inc
f(v) = -(jw/Z_s) * <A_inc, v>_sibc - <n_cond x H_inc, v>_sibc
```

**Sign convention (CRITICAL)**:
- `specialcf.normal(3)` points outward from the mesh domain (air)
- For a hole: this points INTO the conductor = WRONG for SIBC
- **n_cond = -specialcf.normal(3)** = outward from conductor (SIBC convention)
- Alternatively: `Cross(H_inc, specialcf.normal(3))` gives n_cond x H_inc directly

**Subtraction cancellation limitation**: For coil+workpiece, A_scat approx -A_inc
on SIBC surface (85% cancellation for copper). H1 interpolation error (~1%) in
A_inc amplifies to ~12% in the residual A_total = A_inc + A_scat. Use total-field
approach instead for coil+workpiece cases.

**Validated**: sphere (analytical A_inc, H_inc): -2.7% vs analytical.

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

### Hole Approach: Wrong Physics for SIBC

Hole = workpiece removed from mesh. Natural BC is Neumann (n x H = 0 = PEC).
Adding Robin penalty perturbs PEC toward finite impedance, but the baseline
is wrong: PEC = total screening, while SIBC = partial screening (transparent
interior with surface current). The hole approach has a systematic +6% error
that does not decrease with mesh refinement.

### Scattered-Field + H1 Interpolation of A_inc: Cancellation Error

H1 interpolation of Biot-Savart A_inc at mesh vertices introduces ~1% error.
With 85% cancellation (A_scat approx -A_inc), the residual error is amplified
to ~12% in H_t. Fine mesh refinement barely helps (+26% -> +23%).

### Volume Form for RHS: Curl-Free Violation

`-int H_inc_interp . curl(v) dx` fails (-99%) because H1 interpolant of H_inc
is NOT curl-free. Integration by parts identity only holds for the exact field.
The H1 interpolant has inter-element discontinuities -> spurious curl contributions.

### Pointwise Evaluation on SIBC Boundary: Normal Component Explosion

Robin penalty |jw/Z_s| ~ 1e9 creates |A_n| ~ 1e6 on SIBC boundary while
|A_t| ~ 1e-8. Tangential projection at points near (not on) boundary requires
cancellation of 10^14 -> fails completely. Only BND integral works.

## calc_fem_kelvin.py Panel: Fixed (2026-04-09)

**Production panel** for Cubit workflow. Now supports hole approach with correct SIBC.

**Fixes applied**:
1. Robin sign: `-jw/Z_s` -> `+jw/Z_s` (hole = SIBC on external boundary of air)
2. Boundary name: `sibc` (hole) + `wp_surface` (legacy internal interface)
3. Tangential projection: `|A_t|^2 = |A|^2 - (A.n)^2` via `specialcf.normal(3)`
4. Energy-balance P extraction for internal interface fallback
5. Auto-detect: `is_hole = "workpiece" not in materials`

**Validated (copper, 7 kHz)**:
- **Interface (recommended)**: H_t = 15.81 (-2.9% vs BEM), P = 5.11e-6 (-2.3%)
- Hole: H_t = 17.25 (+5.9% vs BEM), P = 6.08e-6 (+16.2%)
**Before fix**: H_t ~ 0, P ~ 0 (BND integral lacked tangential projection)

**Cubit .jou for interface approach (recommended)**:
```python
# Keep workpiece as separate volume, label the shared face "wp_surface"
block 1 add volume <air_id>
block 1 name "air"
block 2 add volume <coil_id>
block 2 name "coil"
block 3 add volume <wp_id>
block 3 name "workpiece"
sideset 1 add surface <shared_face_ids>
sideset 1 name "wp_surface"
radia_export netgen "model.vol" order 2 overwrite
```

**Cubit .jou for hole approach** (if interface not possible):
```python
subtract volume <wp_id> from volume <air_id>
sideset 1 add surface <hole_face_ids>
sideset 1 name "sibc"
block 1 add volume <air_id>
block 1 name "air"
block 2 add volume <coil_id>
block 2 name "coil"
radia_export netgen "model.vol" order 2 overwrite
```

## Summary: When to Use What

| Need | Method | Script |
|------|--------|--------|
| **P_total, H_t (fast)** | Scalar BIE (BEM) | calc_heating_bem.py |
| **P_total, H_t, L, B (FEM)** | Total-field + interface + BND | **calc_fem_kelvin.py** |
| **Coil optimization** | BEM for P, FEM for L | both |
| **Validation (sphere)** | Scattered-field | verify_sphere_sibc.py |

BEM is fast for P_total (162 DOFs vs 150k). FEM is needed for field distribution
and inductance. Panel FEM interface approach gives -2.9% H_t vs BEM.
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
the need for pre-computed nodal values. (Discussed with Joachim previously.)

## specialcf.normal WORKS — Sign Convention is Critical (2026-04-09 CORRECTED)

`specialcf.normal(3)` DOES work in HCurl BND LinearForm. The previous claim
that it "returns near-zero" was due to a **sign error**.

`specialcf.normal(3)` points outward from the mesh domain (air). For a hole
(conductor removed), this is INWARD into the conductor — opposite to SIBC
convention (outward from conductor).

```python
# CORRECT: negate specialcf.normal for SIBC conductor-outward convention
n_mesh = specialcf.normal(3)         # outward from air = into conductor
n_cond = -n_mesh                      # outward from conductor (SIBC)

# n_cond x H: use Cross(H, n_mesh) which equals -(n_mesh x H) = n_cond x H
nxH = Cross(H_inc_cf, n_mesh)        # works for ANY geometry

# OR explicit cross product:
nxH = CF((n_cond[1]*H[2] - n_cond[2]*H[1],
          n_cond[2]*H[0] - n_cond[0]*H[2],
          n_cond[0]*H[1] - n_cond[1]*H[0]))
```

Verified on sphere: ||f_specialcf - f_explicit|| / ||f_explicit|| = 1.3e-4 (roundoff).
Both give H_t_rms = 88.34 A/m (-2.7% vs analytical 90.83 A/m).

## Biot-Savart Path Integration: Thick Coil Error

phi_inc from path integration (compute_phi_inc_from_loop) assumes a filamentary
coil. For thick coils (wire radius a_coil comparable to distance to workpiece),
the filament approximation introduces error. Use BiotSavartCF (multipole) or
T0 source/sink for accurate coil fields.

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

## FEM Open Boundary: Kelvin vs Dirichlet

Kelvin transform provides exact open boundary. Dirichlet truncation at finite
distance introduces small error but is much simpler.

For scattered-field copper cylinder (7 kHz):
- **Kelvin** (R_air=60mm, R_kelvin=120mm): P = 6.62e-6 W
- **Dirichlet** (R_air=100mm): P = 6.59e-6 W (negligible difference)

Kelvin IS important for the verify_sphere_sibc.py uniform-field case
(ratio=10, steel), where Dirichlet gives +435% error. But for typical
coil+workpiece (localized source, ratio < 1): Dirichlet at R_air ~ 3*R_coil is adequate.

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
