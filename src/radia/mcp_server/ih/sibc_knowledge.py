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

## Coupled BEM Coil Terminal Inductance Change (2026-04-12)

`bem_coupled_solver.CoupledBEMSolver` computes the **physically correct
coil terminal inductance change** caused by an SIBC workpiece, by
iterating coil EFIE + workpiece scalar BIE+SIBC with a per-DOF
back-reaction RHS.

### What it returns

```python
from radia.bem_coupled_solver import CoupledBEMSolver
solver = CoupledBEMSolver(mesh_coil, mesh_wp)
result = solver.solve(Z_s, omega, max_iter=10, tol=1e-3, relax=0.5)

L_air   = result['L_air']      # uncoupled coil-only inductance
L_total = result['L_total']    # coupled coil terminal inductance
Delta_L = result['Delta_L']    # = L_total - L_air, sign-correct
P_total = result['P_total']    # workpiece dissipation [W]
```

### Sign behavior (validated 2026-04-12)

| Workpiece | mu_r | Delta_L sign | Physics |
|---|---|---|---|
| Non-magnetic conductor (Cu, Al) | 1 | **negative** | Lenz screening |
| Weakly ferromagnetic | 2-10 | negative (smaller) | Lenz still dominates |
| Ferromagnetic (steel) | 100-1000 | **positive** | Skin-layer flux concentration |
| Asymptotic high freq | any | saturates | Cu @ 1 MHz: -1.02 nH (PEC limit) |

Frequency sweep (copper, mu_r=1, R_coil=30mm, R_wp=10mm, H_wp=20mm):

| freq    | delta    | L_air     | L_total   | Delta_L      |
|---------|----------|-----------|-----------|--------------|
| 100 Hz  | 6.61 mm  | 86.671 nH | 86.356 nH | -0.316 nH    |
| 1 kHz   | 2.09 mm  | 86.671 nH | 85.899 nH | -0.772 nH    |
| 10 kHz  | 0.66 mm  | 86.671 nH | 85.729 nH | -0.942 nH    |
| 100 kHz | 0.21 mm  | 86.671 nH | 85.673 nH | -0.998 nH    |
| 1 MHz   | 0.066 mm | 86.671 nH | 85.655 nH | -1.016 nH    |

mu_r sweep (steel sigma=2e6, half=5mm, f=50kHz) shows the sign change:

| mu_r | Delta_L      |
|------|--------------|
| 1    | -0.831 nH    |
| 10   | -0.462 nH    |
| 100  | **+0.342 nH** |
| 1000 | +1.304 nH    |

### Implementation key (do NOT use scalar rescale)

The back-reaction RHS is a **per-HDivSurface-DOF vector**:

    f_back[i] = int v_i.Trace() . A_wp dS_coil   (i = 0..n_J-1)

Built via a NGSolve LinearForm with A_wp as a CoefficientFunction sum
of M analytic 1/r kernels in (x, y, z). The previous (v1) implementation
used a scalar rescale `alpha * SL @ J_coil` and produced wrong-signed
Delta_L. Saved as `bem_coupled_solver_v1_buggy.py.bak` for reference.

### Limitations

- Linear SIBC (Dowell formula) only — nonlinear ESIM not yet supported
  in the coupled solver
- Scalar BIE limit captures only exterior scattered field; full mu_r
  flux concentration in the workpiece volume requires MFIE/PMCHWT
  (not yet implemented)
- Quantitative validation against FEM-Kelvin / analytical sphere is
  still pending; sign and trends are verified

### Wired into

- `calc_inductance.py::_run_coupled_bem`: called when `--workpiece` is
  set with `--impedance-model dowell`
- IH panel display shows `L (air)`, `delta L`, `L (eff)`, `R (added)`,
  iterations, skin depth

## Per-panel local curvature SIBC (2026-04-12)

NGSolve mesh-driven per-panel local radius for the SIBC cell problem.
Eliminates the need for the user to manually set ``half_thickness`` to
match the workpiece geometry.

### What it does

For each workpiece panel, compute a local radius from the **discrete
normal-angle** between adjacent panels:

```
edge-adjacent panel pair (i, j):
    angle_ij = arccos(n_i . n_j)
    dist_ij  = |c_j - c_i|
    R_ij     = dist_ij / angle_ij      (= 1 / local curvature)
R_local[i] = percentile_10({R_ij over all edge neighbors of i})
```

The percentile-10 aggregation picks up the maximum local curvature
direction (= principal direction with the smallest R) and is robust
against:

- cylinder axial neighbors (parallel normals -> R -> infinity, filtered)
- sliver triangles producing spurious tiny R (clipped)
- mesh discretization noise (median-like robustness)

### CLI / GUI

```bash
calc_inductance.py --workpiece sibc --impedance-model esim \
    --use-local-curvature
```

GUI (IH panel) -> "Per-panel curvature: on / off" combo, visible
whenever Workpiece SIBC mode is on. Off by default for backward
compatibility.

### Validation (2026-04-12)

| Geometry        | mesh maxh | R_local mean | expected | error |
|-----------------|-----------|--------------|----------|-------|
| Sphere R=25 mm  | R/4       | 22.79 mm     | 25 mm    | -8.8% |
| Sphere R=25 mm  | R/8       | 23.32 mm     | 25 mm    | -6.7% |
| Cylinder R=10 mm side | R/3 | 10.72 mm     | 10 mm    | +7.2% |
| Cylinder caps         | R/3 | 654.9 mm     | flat     | OK    |
| Flat plate            | -   | 507 mm       | flat     | OK    |

The chord-vs-arc discretization error decreases with mesh refinement.

### End-to-end demo (sphere workpiece, R=15 mm, steel mu_r=100, 50 kHz)

| Case                                | half_thickness | P_total [W] | error |
|-------------------------------------|---------------|-------------|-------|
| A. Wrong global R (user typo: 5 mm) | 5 mm          | 0.180       | **+189%** |
| B. Correct global R = 15 mm         | 15 mm         | 0.062       | reference |
| C. **Per-panel auto (--use-local-curvature)** | (any)  | **0.067**   | **+8%** |

The headline is **case A vs C**: the user no longer needs to know or
type the workpiece radius — the mesh tells the solver. A typo in
``half_thickness`` becomes harmless, and a non-trivial geometry
(ellipsoid, free-form workpiece) gets per-panel-correct Z_s without
any user input.

### Files

- `src/radia/panels/calc_inductance.py::_compute_panel_local_radii` —
  the discrete-normal-angle radius extractor
- `src/radia/panels/calc_inductance.py::_compute_wp_impedance_from_panels` —
  ``mesh_wp`` and ``use_local_curvature`` parameters; both ESIM and
  Dowell loops use per-panel R
- `src/radia/esim_cell_problem.py::ESIMFiniteSlabSolver` —
  ``set_radius(R)`` and ``solve(H0, R_local=...)`` for fast per-panel R
- `src/radia/radia_ih.py` — GUI combo "Per-panel curvature: on/off"

### Result keys (added)

- `wp_use_local_curvature`: bool
- `wp_R_local_min`: smallest per-panel R found [m]
- `wp_R_local_max`: largest per-panel R found [m]

### Limitations

1. The cell problem is still 1D radial — captures the maximum
   principal curvature only. For ellipsoid or doubly-curved
   surfaces (kappa_1 != kappa_2) it approximates as a cylinder of
   ``R = 1/kappa_max``. For a sphere this is exact.
2. The **coupled BEM solver** (`bem_coupled_solver.py`) still uses
   a scalar Z_s. Per-panel Z_s in the coupled solve is Phase 5.
3. The clamp ``R >= 0.5*half_thickness`` protects against sliver
   triangles producing spurious tiny R. Lower the bound for
   sub-millimeter workpieces.

### What used to be the case (do not re-search)

Before 2026-04-12, ``--esim-geometry local_curvature`` only meant
"use 1D radial Bessel with R = global half_thickness". The name was
misleading. Now ``--use-local-curvature`` is the per-panel mesh-driven
flag, distinct from the geometry mode.

See ``memory/sibc_per_panel_curvature.md`` for full implementation
notes; the older ``sibc_curvature_status.md`` is **superseded**.

## Phase 5 (2026-04-12): per-node Z_s in the COUPLED BEM solver

The same per-panel curvature is now wired into the **iterative coupled
BEM solver** ``bem_coupled_solver.CoupledBEMSolver`` by extending the
scalar BIE+SIBC core ``bem_sibc_solver.ScalarBIESIBCSolver`` to
accept a **per-H1-node Z_s ndarray**.

### Implementation

``ScalarBIESIBCSolver.solve(phi_inc, Z_s, omega)``:
- ``Z_s`` may be a complex scalar (legacy) OR a complex ndarray of
  length ``self.ndof``.
- For an ndarray, the system is::

      (1/2 M - DL + diag(gamma) @ SL @ M^-1 @ K) phi = M phi_inc
      gamma_i = Z_s_i / (jw mu_0)

  Each H1 row gets its own Robin coefficient. Verified to reproduce
  the scalar result to 7.6e-16 when the array is uniform.

``calc_inductance.py::_run_coupled_bem``:
- Builds per-panel R via ``_compute_panel_local_radii``
- Computes per-panel Z_s via the Dowell tanh formula
- Projects to H1 nodes by **vertex averaging** (P1 nodal projection)
- Passes the array to ``solver.solve(Z_s=ndarray, omega)``

### Validation: analytical sphere SIBC (Smythe)

Uniform external field H = H0 z_hat on a conducting sphere of radius
R. Smythe's closed-form result::

    J_s(theta) = (3/2) * j*omega*B0*R / (j*omega*mu_0*R + Z_s) * sin(theta)
    H_t_rms    = max|J_s| * sqrt(2/3)

Three SIBC paths on the same sphere mesh:

  1. **scalar(wrong R=5mm)**: legacy, user supplied wrong half_thickness
  2. **per-node(mesh extractor)**: production path
  3. **per-node(true R)**: machine-precision regression

| Sphere R=10mm        | xi    | scalar(wrong R) | **per-node(mesh)** | per-node(true) |
|----------------------|-------|-----------------|---------------------|----------------|
| Cu 1 MHz (asymptot.) | 72    | -0.65%          | -0.65%              | -0.65%         |
| Cu 1 kHz             | 2.3   | -0.53%          | -0.68%              | -0.68%         |
| **Cu 10 Hz**         | 0.48  | **+31.0%** ❌  | **+2.9%** ✓        | -0.79%         |
| Steel mu=100 1 kHz   | 8.9   | -0.92%          | -0.90%              | -0.90%         |

The Cu/10Hz case shows the headline: scalar with wrong R gives +31%
error, per-node recovers to +2.9% (**11x improvement**). In the
asymptotic regime Z_s is R-independent so all paths agree.

### Validation: prolate spheroid (spatially varying R)

Sphere alone has uniform R so it cannot test that per-panel
**variation** is captured correctly. The prolate spheroid
``x²/b² + y²/b² + z²/a² = 1`` (a > b) has principal radii::

    pole (cos t = 1):    R = b² / a   (tip-of-cigar radius)
    equator (cos t = 0): R = b        (equatorial radius)

For a/b = 4 the pole-vs-equator R ratio is 16x. The per-panel
extractor's percentile-10 picks the smaller principal radius.

| Cu cigar (a × b)         | xi    | scalar(mean R) | **per-node(mesh)** | per-node(analytic) |
|---------------------------|-------|----------------|---------------------|---------------------|
| 20 × 10 mm, **10 Hz**     | 0.48  | -11.2%         | **+3.1%**           | ground truth        |
| 20 × 10 mm, 1 MHz (asym.) | 318   | +0.00%         | +0.00%              | identical           |
| 40 × 10 mm, **100 Hz**    | 1.5   | -5.1%          | **+1.7%**           | ground truth        |

**This is the definitive validation**: spatially varying principal
curvature is correctly extracted from the mesh AND correctly fed
into the per-node SIBC solver. Per-node beats scalar by 4-7x in the
xi ~ 1 regime; identical to scalar in the asymptotic regime
(regression OK).

### Files

- `src/radia/bem_sibc_solver.py::ScalarBIESIBCSolver.solve` —
  per-node Z_s acceptance + Robin block diag(gamma) construction
- `src/radia/bem_coupled_solver.py::CoupledBEMSolver.solve` —
  passes through per-node Z_s to the SIBC core
- `src/radia/panels/calc_inductance.py::_dowell_Zs`,
  `_build_per_node_Zs`, `_run_coupled_bem(use_local_curvature=True)`
- `examples/cubit_panels/inductance/verify_per_node_sibc_sphere.py`
- `examples/cubit_panels/inductance/verify_per_node_sibc_spheroid.py`

### Cross-check vs FEM-Kelvin SIBC (validated 2026-04-12)

The coupled BEM was independently validated against the FEM-Kelvin
SIBC pipeline (`calc_fem_kelvin.py --impedance sibc`) on the same
`radia_model.vol` from `ih_sample.jou`:

| Material | mu_r | f      | L_BEM     | L_FEM     | diff      |
|----------|------|--------|-----------|-----------|-----------|
| copper   | 1    | 50 kHz | 84.31 nH  | 84.56 nH  | **+0.29%** |
| steel    | 100  | 50 kHz | 87.92 nH  | 89.43 nH  | +1.72%    |

L_air (coil only) = 87.81 nH. The 0.3% agreement on copper is the
strongest validation we have for the coupled BEM solver. Both methods
report the correct sign in both regimes (Lenz screening for copper,
flux concentration for ferromagnetic steel).

**Canonical regression script**:
`examples/cubit_panels/inductance/compare_bem_coupled_vs_fem_kelvin.py`

Run after any change to `bem_coupled_solver.py` to confirm the cross-
check is still tight.

## Summary: When to Use What

| Need | Method | Script |
|------|--------|--------|
| **P_total, H_t (fast)** | Scalar BIE (BEM) | calc_heating_bem.py |
| **L change vs workpiece** | **Coupled BEM (2026-04-12)** | **calc_inductance.py --workpiece** |
| **P_total, H_t, L, B (FEM)** | Total-field + interface + BND | **calc_fem_kelvin.py** |
| **Coil optimization** | BEM for L+P, FEM for field distribution | both |
| **Validation (sphere)** | Scattered-field | verify_sphere_sibc.py |

BEM is fast for P_total (162 DOFs vs 150k). The coupled BEM is the
right tool for "how does the coil terminal L change when I add a
workpiece" — it gets the sign right across all mu_r/frequency regimes.
FEM is needed for field distribution and full mu_r flux concentration.
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
