"""
SIBC / ESIM knowledge for induction heating MCP server.

=============================================================================
2026-06-24 PANEL PIPELINE NOTICE
=============================================================================
The IH panel is now a unified production panel with six EM methods plus
three integrated thermal follow-up methods:

  (a) PEEC inductance              -> calc_inductance.py --coil-solver peec
  (b) BEM-A inductance             -> calc_inductance.py --coil-solver bem-a
  (c) PEEC+BEM weak coupling       -> calc_inductance.py --coil-solver peec --vol <wp>
  (d) BEM-A+BEM weak coupling      -> calc_inductance.py --coil-solver bem-a --vol <wp>
  (e) PEEC + FEM wp + Kelvin       -> calc_fem_kelvin.py --formulation total
  (f) Full FEM A-V + wp SIBC       -> calc_fem_coilmesh.py
  (g-i) Thermal 3D/static, 3D/rotating, and 2D-axisym -> calc_heat*.py

The current GUI owns per-panel browse rows under a top-level Working
folder.  The workpiece mesh is `IHPanel.wp_vol`, not the hidden legacy
`AnalysisWindow._vol_edit` compatibility shim.

References below that mention retired helpers such as
`calc_heating_bem.py` are historical physics notes.  Current panel
entry points are the calc scripts listed above.

Historical context:
  2026-04-17: BEM modules moved into a research reference lane.
  2026-06-28: BEM modules promoted to radia.* APIs; runnable references moved
              to validation_test/induction_heating/bem_reference/.
  2026-04-18: T0 + A-V compound retired from panel (gap-corner 1/r cusps)
  2026-04-19: calc_peec_bem + calc_fem_coilmesh become the panel
              (A-V is BACK, this time gapped-torus-only + proper source/sink)
  2026-04-24: calc_fem_kelvin scattered formulation retired from the panel
  2026-06-24: docs refreshed for unified calc_inductance + working-folder UI
=============================================================================
"""

from ..common import load_prompt

IH_PEEC_FEM = """
# IH Architecture: PEEC+FEM (v4.6.0)

## Two Solver Paths

| Path | Coil | Workpiece | Panel Method |
|------|------|-----------|-------------|
| **PEEC+FEM** (default) | PEEC filaments (no mesh) | FEM-SIBC + Kelvin | "PEEC+FEM" |
| **ALL FEM** (reference) | FEM volume mesh | FEM-SIBC + Kelvin | "FEM" |

### Path 1: PEEC+FEM (Production, default)

```
STEP coil file -> auto centerline -> PEEC filaments (ms, no mesh)
    -> Biot-Savart source field
    -> FEM-SIBC + Kelvin on workpiece (sparse, seconds)
    -> Picard back-reaction (1 step)
```

Advantages:
- **No coil mesh** (STEP -> filaments, auto cross-section extraction)
- **Coil rotation/translation** = coordinate transform only (no remesh)
- Shifted Compact AMS preconditioner (HCurl p=1, TaskManager parallel)
- HACApK saddle-point for large coils (N > 3000)
- PRIMA MOR for broadband sweep (if needed)

Solver options in panel:
- Dense LU: fastest for N < 3000 (default)
- HACApK saddle-point: for N >= 3000
- PRIMA: for frequency sweep (rarely needed in IH)

Files:
- `calc_peec.py` (Layer 4 computation)
- `coil_from_cad.py` (STEP -> filaments)
- `peec_topology.py` (PEECCircuitSolver)
- `prima_hacapk.py` (PRIMA MOR)

### Path 2: ALL FEM (Reference, validation)

```
Cubit .vol (coil + workpiece + air + Kelvin) -> full volume FEM
    -> Omega or A-formulation + SIBC + Kelvin
    -> Direct solve (pardiso/bddc)
```

Advantages:
- Full physics (volume currents, nonlinear materials)
- Well-validated (2D axisym reference: 1.15%)
- No approximations beyond mesh discretization

Limitations:
- Requires full volume mesh (coil + air + Kelvin)
- Coil remesh required for geometry changes
- Slower for parametric studies

Files:
- `calc_fem_kelvin.py` (Layer 4 computation)
- `scalar_potential_solver.py` (Omega formulation)
- `kelvin_solver.py` (Kelvin helpers)

### When to Use Which

| Scenario | Recommended |
|----------|------------|
| IH coil design (fast iteration) | **PEEC+FEM** |
| Coil placement optimization | **PEEC+FEM** (rotation = coordinate transform) |
| Nonlinear workpiece (BH curve) | **PEEC+FEM** (ESIM in FEM workpiece) |
| Validation / reference | **ALL FEM** (full physics) |
| Complex coil geometry (not helix) | **ALL FEM** (until STEP extraction matures) |

### Workpiece Sub-options (both paths)

| Workpiece | Description |
|-----------|-------------|
| off | Coil self-impedance only |
| SIBC | Linear surface impedance (Cu, Al, fixed mu_r) |
| ESIM | Nonlinear 1D cell problem (steel, BH curve, freq-dependent mu) |
"""

IH_SIBC_OVERVIEW = """
# SIBC for Induction Heating: Method Selection

## Solver Accuracy (Smythe sphere analytical validation, 2026-04-13)

| Method | H_t error vs analytical | Conditions | Status |
|--------|------------------------|------------|--------|
| **FEM scattered-field** | **+/-3%** | All f, Cu+Steel | **Reference (A_inc=exact CF)** |
| **BEM Scalar BIE SIBC** | **-1.6%** | R/d>1.5, Cu+Steel | **Production solver** |
| **FEM total-field (hole)** | **L: +2.5%** | Kelvin, Cu 50kHz | **L production, H_t TBD** |
| FEM total-field (interface) | **-34% systematic** | All f, Cu+Steel | **WRONG: do not use** |
| BEM EFIE (HDivSurface) | -7% to -66% | Broken | Factor-3 eigenvalue error |

## Scalar BIE + SIBC (Production solver for H_t and P)

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
- **mu_r is handled correctly through Z_s** (no PMCHWT needed)
- Uses BOTH DL and SL -- captures full A+H physics (unlike EFIE)

**Validated (2026-04-13) vs exact Bessel solution** (no SIBC approx):
- Copper (R/d 1.5-33.8): **-1.0% to -1.6%**
- Steel mu_r=100 (R/d 2.8-62.8): **-1.5% to -1.7%**
- SIBC approximation itself fails at R/delta < 1 (low frequency limit)
- **EFIE (HDivSurface+LaplaceSL) is broken** (-7% to -66%): factor-3
  eigenvalue error from using A only (SL), dropping H contribution (curl SL)

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

## FEM Total-Field + Hole Approach (Production for L; H_t under validation)

**File**: `src/radia/panels/calc_fem_kelvin.py --formulation total`

Total-field formulation: solve for full A with Robin BC on the hole boundary.
- **Hole approach is REQUIRED**: workpiece subtracted from mesh, NOT meshed.
- **L (inductance) is correct**: +2.5% vs analytical torus (3D Kelvin, Cu 50kHz)
- **H_t extraction**: `H_t = |jw/Z_s| * |A_t|` from single-side BND trace.
  Hole approach has unambiguous trace (air side only). Under validation.
- 2D axisym validation: L < 1%, P < 2% vs full-resolution (Cu/Steel/Al, R/d=3-160)

**Interface approach (-34% error) is WRONG**: Robin BC on internal interface returns
values from the workpiece interior (A ~ 0). Do NOT use `--formulation total` with
workpiece as a meshed volume.

**For H_t/P with analytical A_inc (validation only)**: use `--formulation scattered`.

## FEM Scattered-Field (Reference solver for H_t and P)

**Files**:
- `src/radia/panels/calc_fem_kelvin.py --formulation scattered`
- `validation_test/eddy_current_analytical_validation/sphere_uniform_field.py`
- SHOWCASE NOTEBOOK `docs/eddy_current_analytical_validation/validation_suite.ipynb`
  (live Dodd-Deeds rod impedance + frequency sweep, plus the documented Smythe
  sphere FEM-scattered<3% / FEM-total~34% / BEM cross-validation table).

Decomposition: A = A_inc + A_scat.  A_inc is known (analytical CF or
free-space FEM solve). Solve for A_scat with Robin BC.

**Two modes**:

1. **Sphere / uniform field**: A_inc = CF((-B0/2*y, B0/2*x, 0)) -- exact,
   no discretization error. RHS includes both surface terms:
   ```python
   f(v) = -(jw/Z_s) * <A_inc, v>_sibc - <n_cond x H_inc, v>_sibc
   ```

2. **Coil problem (two-step)**: Step 1 solves free-space (no Robin) -> A_inc
   as GridFunction. Step 2 solves with Robin, RHS = -robin * A_inc on SIBC.
   Volume terms cancel because A_inc satisfies the source PDE.

**H_t extraction**: A_total = A_inc (CF) + A_scat (GF).  Since A_inc dominates
the tangential component, numerical errors in A_scat do not corrupt H_t.

**Validated against Smythe sphere (2026-04-13, 14 data points)**:
- Copper (R/delta 1.5 to 33.8): +0.1% to +2.0%
- Steel (R/delta 2.8 to 62.8): -2.9% to -1.0%
- **All conditions: error < 3%**

See `docs/eddy_current_analytical_validation/validation_suite.ipynb (SHOWCASE NOTEBOOK)` for full derivation and validation table.

**Limitation for coil problem**: When A_inc comes from a FEM solve (GridFunction,
not exact CF), subtraction cancellation (A_scat ~ -A_inc on SIBC, 85% for Cu)
may degrade H_t accuracy. Future: use Biot-Savart vector potential via
RadiaField('a') as exact CF for A_inc.

## CRITICAL: What Does NOT Work

### FEM Total-Field H_t: -34% Systematic Error (confirmed 2026-04-13)

`calc_fem_kelvin.py --formulation total` gives **-34% error on H_t** (and
therefore -56% on P) compared to Smythe analytical sphere solution. This
error is **constant across all frequencies and materials** (tested R/delta
1.5 to 62.8, copper and steel).

Root cause: BND trace evaluation of conforming HCurl GridFunction on the
SIBC boundary returns values biased toward the workpiece interior (where
A ~ 0 due to gauge regularization or simply because the conducting body
has no internal source). The bias is independent of Robin magnitude.

**L (inductance) is NOT affected** -- it uses a volume integral.

**Previous "validation"** (commit 3ef7555, 2026-03-28) compared FEM-SIBC vs
EFIE-SIBC (BEM). But EFIE-SIBC has -65% error on sphere (SL eigenvalue mismatch).
"10% agreement with a 65%-wrong method" is NOT validation.

### EFIE-SIBC (HDivSurface): Wrong for Finite Z_s

RT0 (order=0) has zero surface curl. SIBC requires curl_s(J) which vanishes
for RT0 basis. Results are wrong by 65% on sphere at Z_s/(jw*mu0*R) = 10.

### SIBC Implementation: Hole Approach is Correct (2026-04-14)

SIBC = Robin BC on the conductor surface.  The conductor interior is
NOT solved.  Use hole approach: subtract workpiece from mesh, apply
Robin BC on the hole boundary.  Validated in 2D axisymmetric Kelvin FEM:

- Full-resolution (eddy currents resolved at delta/5) vs SIBC (hole + Robin)
- L: < 1% for Cu, Steel (mu_r=100), Al at R/delta = 3 to 160
- P: < 2% for all cases
- Z_s: cylindrical Bessel `rho*gamma*I1(ga)/I0(ga)` (exact for solid cylinder)
- Interface approach (wp meshed as air + Robin) is WRONG: flux bypasses
  Robin through transparent interior (steel: 58% of correct L)

For H1/phi formulation: Robin term = `(jw/Z_s)/r * u * v * ds("wp_bnd")`
For 3D HCurl formulation: Robin term = `(jw/Z_s) * A_t . v_t * ds("sibc")`

Script: `validation_test/eddy_current_analytical_validation/reference_2d_axisym.py`

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

## calc_fem_kelvin.py Panel (2026-04-13)

**Production panel** for Cubit workflow. Two formulations:

- `--formulation total` (default): L is accurate, H_t/P has -34% systematic error
- `--formulation scattered`: L AND H_t/P are accurate (for A_inc = exact CF)

**EMMaterial integration (2026-04-13)**: `--material copper` now auto-sets
sigma=5.8e7 and mu_r=1 (via `EMMaterial.from_args()`). No need to pass
`--sigma` separately.

**For H_t/P accuracy in coil problems**: the two-step scattered solve uses
A_inc from free-space FEM (GridFunction). Accuracy depends on A_inc quality.
Future: Biot-Savart A via RadiaField('a') as exact CF.

**Cubit .jou for hole approach (REQUIRED)**:
```python
# Workpiece is subtracted from air (NOT meshed), SIBC on hole boundary.
# Validated in 2D axisym: L < 1%, P < 2% vs full-resolution eddy current.
# Interface approach (wp meshed as volume) is WRONG: flux bypasses Robin
# through transparent interior (steel: 58% of correct L).
subtract volume <wp_id> from volume <air_id> keep_tool
imprint volume <air_id> <coil_id>
merge volume <air_id> <coil_id>
# Do NOT mesh the workpiece. Do NOT add it to a block.
sideset 1 add surface in volume <wp_id>   # shared surfaces become air boundary
sideset 1 name "sibc"
block 1 add volume <air_id>
block 1 name "air"
block 2 add volume <coil_id>
block 2 name "coil"
export netgen "model.vol" order 1 overwrite
```

**Do NOT use interface approach** (workpiece meshed as separate volume).
The Robin BC on an internal interface has wrong-side trace evaluation:
NGSolve returns values from the workpiece interior (A ~ 0), not from
the air side. This causes -34% systematic H_t error (total-field) and
58% L error (steel). The hole approach eliminates the ambiguity because
the SIBC boundary is an external boundary of the air domain with a
single-side trace.

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

- Linear SIBC (Dowell formula) only -- nonlinear ESIM not yet supported
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
type the workpiece radius -- the mesh tells the solver. A typo in
``half_thickness`` becomes harmless, and a non-trivial geometry
(ellipsoid, free-form workpiece) gets per-panel-correct Z_s without
any user input.

### Files

- `src/radia/panels/calc_inductance.py::_compute_panel_local_radii` --
  the discrete-normal-angle radius extractor
- `src/radia/panels/calc_inductance.py::_compute_wp_impedance_from_panels` --
  ``mesh_wp`` and ``use_local_curvature`` parameters; both ESIM and
  Dowell loops use per-panel R
- `src/radia/esim_cell_problem.py::ESIMFiniteSlabSolver` --
  ``set_radius(R)`` and ``solve(H0, R_local=...)`` for fast per-panel R
- `src/radia/radia_ih.py` -- GUI combo "Per-panel curvature: on/off"

### Result keys (added)

- `wp_use_local_curvature`: bool
- `wp_R_local_min`: smallest per-panel R found [m]
- `wp_R_local_max`: largest per-panel R found [m]

### Limitations

1. The cell problem is still 1D radial -- captures the maximum
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

- `src/radia/bem_sibc_solver.py::ScalarBIESIBCSolver.solve` --
  per-node Z_s acceptance + Robin block diag(gamma) construction
- `src/radia/bem_coupled_solver.py::CoupledBEMSolver.solve` --
  passes through per-node Z_s to the SIBC core
- `src/radia/panels/calc_inductance.py::_dowell_Zs`,
  `_build_per_node_Zs`, `_run_coupled_bem(use_local_curvature=True)`
- `examples/cubit_panels/inductance/verify_per_node_sibc_sphere.py`
- `examples/cubit_panels/inductance/verify_per_node_sibc_spheroid.py`
- `examples/cubit_panels/inductance/verify_per_node_sibc_torus.py`

### Residual cleanup (commit ca5cb46, 2026-04-12)

After Phase 5 four small follow-ups landed:

**Torus validation script** added (the adversarial case for the
percentile-10 extractor). On a thin torus (R/a = 6) the mean per-
panel R error is +245%, but the integrated H_t_rms stays within 5%
of the analytical reference because vertex averaging in
``_build_per_node_Zs`` smooths the per-panel noise. Documented as
a known limitation; the right fix is a 1-ring vertex expansion or
analytic Frenet-Serret / shape-operator approach.

**ih_bem_sample.jou** workpiece changed from cylinder to sphere
(R=10mm). The cylinder is locally cylindrical everywhere so per-
panel had no visible effect. With a sphere, the per-node coupled
BEM gives a coil-terminal dL that is **10x smaller** than the
wrong-global-R scalar path -- the correction is now visible in the
IH panel output.

**Adaptive sliver clamp**: replaced the global
``R >= 0.5 * half_thickness`` clamp with a per-panel adaptive lower
bound ``R_floor[i] = 0.5 * panel_diameter[i]`` inside
``_compute_panel_local_radii``. Below this floor the chord-vs-arc
approximation breaks down so the discrete-normal-angle estimate is
geometrically meaningless. Removed the duplicate global clamps in
``_run_coupled_bem`` and ``_compute_wp_impedance_from_panels``;
the floor lives in the extractor only.

  Validation:
    Sphere / spheroid -- unchanged (regression OK)
    Torus thin (R/a=6)  -- H_t err 4.80% → 3.32%
    Torus fat (R/a=1.5) -- H_t err 3.34% → 0.48%

**``_build_per_node_Zs`` order > 1 projection**: added an optional
``fes`` parameter. When supplied AND ``fes.globalorder > 1``, the
per-panel Z_s is wrapped in a SurfaceL2(order=0) GridFunction and
L2-projected onto fes via ``GridFunction.Set``. Edge / face DOFs
get the right H1 mass-matrix solve instead of the panel-mean
fallback. Default ``fes=None`` keeps the legacy P1 vertex
averaging on the production code path.

### Known limitations as of 2026-04-12

1. **Torus extractor**: percentile-10 of 3 edge-neighbors fails on
   panels at z = ±a (top/bottom of poloidal circle). Fix is
   1-ring expansion or shape-operator; deferred to next session.

2. **MFIE / PMCHWT not implemented**: the SIBC scalar BIE captures
   the surface skin layer correctly but cannot reach the high-mu_r
   ferromagnetic flux concentration limit. For that an MFIE term
   with the gradient of the Green's function is needed.

3. **Order > 1 L2 projection** has not been exercised in production.
   The legacy vertex-averaging path is the only one tested for the
   coupled BEM end-to-end.

### Cross-check vs FEM-Kelvin SIBC (validated 2026-04-12)

The coupled BEM was independently validated against the FEM-Kelvin
SIBC pipeline (`calc_fem_kelvin.py --impedance sibc`) on the same
`radia_model.vol` from `ih_bem_sample.jou`:

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
| **L, B distribution (FEM)** | **Total-field + hole + Kelvin** | **calc_fem_kelvin.py** |
| **P_total, H_t (FEM)** | Total-field + hole (under validation) | calc_fem_kelvin.py |
| **Coil optimization** | BEM for L+P, FEM for field distribution | both |
| **Validation (sphere)** | Scattered-field | verify_sphere_sibc.py |

BEM is fast for P_total (162 DOFs vs 150k). The coupled BEM is the
right tool for "how does the coil terminal L change when I add a
workpiece" -- it gets the sign right across all mu_r/frequency regimes.
FEM is needed for field distribution and full mu_r flux concentration.
"""

IH_ESIM = """
# ESIM for Induction Heating Workpieces

The ESIM (Effective Surface Impedance Method) general technique — 1D cell
problem mathematics, Karl iteration, `radia.esim_cell_problem` module API
— now lives in `radia_mcp.radia_ngsolve.esim_knowledge` (callable as
`mcp-server-radia-ngsolve.esim(topic)`).  This topic only documents the
**IH-specific** application of ESIM (workpiece SIBC for steel/ferrite under
an induction coil).

For general ESIM theory, call:
- `esim(topic="overview")`        — when ESIM vs linear SIBC; nonlinear conductors
- `esim(topic="cell_problem")`    — 1D BVP and supported geometries
- `esim(topic="karl_iteration")`  — Picard relaxation + convergence pitfalls
- `esim(topic="module_api")`      — `ESIMFiniteSlabSolver` + BEM/FEM coupling examples

## IH-specific use of ESIM (this topic)

### Typical IH workpiece materials needing ESIM
- **Carbon steel** (S45C, S50C): saturates ~1.7 T, mu_r drops from ~500 (low H)
  to ~10 (saturation). Linear SIBC underpredicts P_wp by 30–60% near the surface.
- **Stainless ferritic** (SUS430): similar saturation but lower base mu_r.
- **Electrical steel** (35JN230, 50A1300): laminated → ESIM with anisotropic mu
  along the lamination plane (out-of-plane mu effectively zero).
- **Ferrite cores** (Mn-Zn, Ni-Zn): use complex mu = mu' - j mu''.  Karl
  iteration handles the loss term automatically.

### IH workflow integration

1. Get an initial Z_s estimate via `esim.solve(H0)` at a small seed |H_t|
   (~5 A/m) — the outer Karl loop will refresh it on iter 0.
2. Solve outer BEM/FEM with current Z_s -> read mesh-RMS H_t on the
   workpiece surface.
3. Refresh Z_s = esim.solve(H_t)['Z']; under-relax: Z_s = relax * Z_s_new
   + (1 - relax) * Z_s_old (relax = 0.5 default).
4. Repeat until |dZ_s| / |Z_s| < tol (default 1e-3, ~5-10 iter typical).
5. Final P_wp + L from the converged Z_s.

### Production scripts (v4.46+; Karl iteration wired in all four)

All four scripts accept `--impedance-model esim --bh-file <BH.txt>` and
report converged Z_s, iteration count, and convergence history in JSON:

| Script | Coil model | Workpiece model | Notes |
|---|---|---|---|
| `calc_inductance.py --coil-solver peec` | PEEC filament from STEP | scalar BEM-SIBC | Fast (~few s); 1-D Karl |
| `calc_inductance.py --coil-solver bem-a` | BEM-A (Weggler EFIE, RWG on .vol) | scalar BEM-SIBC | Same Karl loop |
| `calc_fem_kelvin.py --impedance esim` | PEEC filament source | FEM-HCurl A with Robin | PEEC+FEM+Kelvin |
| `calc_fem_coilmesh.py --impedance-model esim` | FEM A-V volumetric coil | FEM with Robin SIBC | Full FEM A-V; re-assembles per Karl iter |

#### Common CLI flags

```
--impedance-model {sibc|esim}    sibc = linear Dowell; esim = Karl
--bh-file BH_FILE                required for esim; 2-col [H[A/m] B[T]]
--esim-max-iter N                outer Karl cap (default 15)
--esim-tol T                     |dZ|/|Z| convergence (default 1e-3)
--esim-relax R                   under-relaxation (default 0.5)
--half-thickness D               ESIM cell radius / half-slab [m]
```

`calc_fem_kelvin` has the same physics but exposes `--max-iter` (not
`--esim-max-iter`) and currently does NOT take an `--esim-tol` knob.

#### JSON output schema (Karl diagnostics, all four scripts)

```json
{
  "impedance_model": "esim",
  "esim_iterations": 6,
  "esim_converged": true,
  "esim_history": [
    {"iteration": 0, "Z_s_abs": ..., "H_t_rms": ..., "dZ": ..., "t_solve": ...},
    ...
  ],
  "Z_s_wp_real": 0.0172,
  "Z_s_wp_imag": 0.0314,
  ...
}
```

Use `esim_history[-1].dZ` for a final convergence indicator; `esim_history`
length will equal `esim_iterations` (one entry per outer Karl iter).

### Geometry choice for IH workpieces

| Workpiece shape  | ESIM `geometry=` | Notes                                  |
|------------------|------------------|----------------------------------------|
| Bar / cylinder   | `'cylinder'`     | radius = workpiece radius              |
| Plate / coupon   | `'finite_slab'`  | half-thickness = sample / 2            |
| Pipe (thin wall) | `'slab'`         | wall is essentially infinite plate     |
| Coil-formed wire | `'cylinder'`     | wire radius (proximity effect ignored) |

### Pitfalls specific to IH

- **Don't use linear SIBC at low frequency** for steel (< 10 kHz) — the
  field penetrates several mm and BH curve dominates. Use ESIM from the
  start; the converged solution differs by an order of magnitude in P_wp.
- **BH file required for ESIM**: panel + calc scripts raise if
  `--bh-file` is empty when `--impedance-model esim` is selected
  (pre-validation in `radia_ih.py:_build_*_command`).
- **Karl iteration diverges with very high relax**: lower the relax flag
  for SUS430 / S50C above 30 kHz: `--esim-relax 0.3`.  (CLI default is
  0.5; matches `calc_fem_kelvin`'s historical setting.)
- **Single Z_s per workpiece is mesh-RMS averaged**.  For workpieces with
  large H_t variation (sharp corners, coil shadow effects), the current
  implementation under-resolves the spatial saturation pattern.
  Per-element ESIM is on the roadmap; track via the radia-mcp tool.
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

## specialcf.normal WORKS -- Sign Convention is Critical (2026-04-09 CORRECTED)

`specialcf.normal(3)` DOES work in HCurl BND LinearForm. The previous claim
that it "returns near-zero" was due to a **sign error**.

`specialcf.normal(3)` points outward from the mesh domain (air). For a hole
(conductor removed), this is INWARD into the conductor -- opposite to SIBC
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

## FEM Open Boundary: Kelvin Only (Policy 2026-04-14)

**POLICY**: Use Kelvin transformation for all FEM open boundary problems.
Do NOT use large Dirichlet truncation spheres.

Kelvin transform provides exact open boundary with no truncation error:
- 2-sphere Periodic Kelvin (Cubit or OCC): physical + mapped domain
- Periodic BC couples interior and exterior sphere DOFs
- GND vertex at exterior sphere center = Dirichlet 0 = physical infinity
- nu_kelvin = nu0 * (r'/R)^2, r' from exterior center

Validated: 3D FEM-Kelvin L = 90.71 nH vs analytical 88.5 nH (+2.5%, coarse mesh).
2D axisym Kelvin: L < 0.6% for Cu/Steel/Al at R/delta = 3-160.

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


IH_MATHEMATICA_VERIFICATION = load_prompt("ih", "mathematica_verification")

_IH_MATHEMATICA_VERIFICATION_INLINE_DEPRECATED = r"""
# Mathematica Verification Recipes for IH SIBC

Closed-form / symbolic derivations for SIBC quantities, suitable for
side-by-side validation of NGSolve / Radia numerical results. Use these
as the **first sanity check** before debugging FEM/BEM output.

Cross-reference: `radia_mcp.radia_ngsolve.analytical_formulas(topic=
"validation_use_cases")` lists the canonical analytical companion for
each IH analysis.

## Recipe 1: Leontovich Z_s symbolic derivation

Derive `Z_s = (1+j)/(sigma * delta)` from the 1D diffusion equation
inside a semi-infinite conductor. Verify against the limit `omega ->
inf`, `mur -> 1`.

```mathematica
(* 1D Helmholtz in conductor: H''(z) = j*omega*mu*sigma * H(z) *)
(* Solution: H(z) = H0 * Exp[-kc*z], decaying into conductor *)
kc[mu_, sigma_, omega_] := (1 - I)/delta /.
    delta -> Sqrt[2 / (omega * mu * sigma)];

(* Surface impedance: Z_s = E_t / H_t at z=0 *)
(* E from Maxwell: curl E = -dB/dt -> E = (1/sigma) * dH/dz *)
Zs[mu_, sigma_, omega_] := (kc[mu, sigma, omega] / sigma) /.
    Simplify[#, Assumptions -> {sigma > 0, mu > 0, omega > 0}] &;

(* Verify against textbook: Z_s = (1+j)*sqrt(omega*mu/(2*sigma)) *)
ZsTextbook[mu_, sigma_, omega_] := (1 + I)*Sqrt[omega*mu/(2*sigma)];

Simplify[Zs[mu, sigma, omega] - ZsTextbook[mu, sigma, omega],
         Assumptions -> {sigma > 0, mu > 0, omega > 0}]
(* Should return 0 *)
```

## Recipe 2: Smythe sphere induced dipole (BEM validation reference)

Smythe (1950) analytical solution for a conducting sphere in uniform
AC field. **The lab's primary SIBC validation case** (2026-04-13: 14
data points, BEM Scalar BIE: -1.6%, FEM scattered: +/- 3%).

```mathematica
(* Magnetic Reynolds number xi = R / delta *)
(* m-factor (induced dipole moment / applied field volume) *)
mFactor[a_, sigma_, mur_, omega_] := Module[
    {delta, k, ka, num, den},
    delta = Sqrt[2/(omega * mu0 * mur * sigma)];
    k = (1 + I)/delta;
    ka = k * a;
    (* Smythe Eq. (1) -- spherical Bessel based *)
    num = (mur + 2) * (Sinh[ka] - ka*Cosh[ka]) +
          ka^2 * (mur - 1) * Sinh[ka]/3;
    den = (mur - 1) * (Sinh[ka] - ka*Cosh[ka]) -
          ka^2 * Sinh[ka]/3;
    -(a^3) * (1 + 3 * num / den)
];

(* Induced surface current (a*sin(theta)*phi component) *)
JsSphere[a_, sigma_, mur_, omega_, H0_, theta_] := Module[
    {m, B0},
    B0 = mu0 * H0;
    m = mFactor[a, sigma, mur, omega];
    (3/2) * I*omega*B0*a / (I*omega*mu0*a + ZsApprox) *
        Sin[theta] /. ZsApprox -> (1 + I)*Sqrt[omega*mu0*mur/(2*sigma)]
];

(* Tangential H_t_rms (max at equator) for cross-checking BEM SIBC *)
HtRMS[a_, sigma_, mur_, omega_, H0_] :=
    Max[Abs[JsSphere[a, sigma, mur, omega, H0, theta]]
        /. theta -> Pi/2] * Sqrt[2/3];
```

Use this to validate `ScalarBIESIBCSolver.solve(...)['H_t_rms']` to
~1% accuracy across `R/delta = 1.5 .. 60`. Compare via:

```mathematica
mu0 = 4*Pi*1e-7;
HtRMS[0.01, 5.8e7, 1, 2*Pi*1e6, 1.0]   (* Cu sphere 1 MHz *)
HtRMS[0.01, 2e6, 100, 2*Pi*1e3, 1.0]   (* steel 1 kHz *)
```

## Recipe 3: Skin depth frequency scaling (sanity plot)

```mathematica
delta[mu_, sigma_, omega_] := Sqrt[2/(omega * mu * sigma)];

LogLogPlot[
    {delta[mu0, 5.8e7, 2*Pi*f] * 1000,        (* Cu, mm *)
     delta[100*mu0, 2e6, 2*Pi*f] * 1000,      (* steel mu_r=100, mm *)
     delta[mu0, 3.5e7, 2*Pi*f] * 1000},       (* Al, mm *)
    {f, 50, 1e6},
    PlotLegends -> {"Cu", "Steel mu_r=100", "Al"},
    AxesLabel -> {"Frequency [Hz]", "Skin depth [mm]"}
]
```

Use this when picking `--mesh-size` for FEM: the surface mesh should
resolve `delta/5` for full-resolution validation, or `delta/2` for
SIBC + hole approach.

## Recipe 4: Per-panel curvature SIBC correction (Mitzner check)

Symbolic version of the per-panel local curvature SIBC (lab 2026-04-12).
Use this to verify the numerical `R_local` extractor in
`_compute_panel_local_radii`.

```mathematica
(* Mitzner curvature-corrected Z_s for a sphere of radius R *)
ZsMitznerSphere[Zs_, delta_, R_] :=
    Zs * (1 + (1 + I)/2 * delta / R);

(* For a sphere R=10mm at 100 Hz, Cu *)
delta100Hz = delta[mu0, 5.8e7, 2*Pi*100];   (* ~6.6 mm *)
ZsLeontovich = (1 + I) * Sqrt[2*Pi*100 * mu0 / (2 * 5.8e7)];
ZsMitznerSphere[ZsLeontovich, delta100Hz, 0.01]
(* Result is ~33% different from ZsLeontovich -- the regime where the *)
(* lab's per-panel extractor gives +2.9% vs +31% for the scalar global *)
```

## Cross-MCP cross-check

After computing a result with `calc_inductance.py --workpiece sibc`,
recompute the same with Mathematica and compare:

| Lab compute                          | Mathematica reference                |
|--------------------------------------|--------------------------------------|
| `ScalarBIESIBCSolver.solve(...)['H_t_rms']` | `HtRMS[a, sigma, mur, omega, H0]`  |
| `--use-local-curvature` per-node Z_s | `ZsMitznerSphere[ZsLeo, delta, R]`   |
| `--impedance-model esim` after Karl  | (Mathematica `NDSolve` 1D BVP)       |

Run any of these as a sanity check before reporting numbers in a paper
or showing to a grant reviewer.

Reference MCP: call `mcp-server-mathematica` for the Wolfram Language
runtime to actually execute these recipes from a Cubit panel session.
"""


IH_NGSOLVE_RECIPES = load_prompt("ih", "ngsolve_recipes")

_IH_NGSOLVE_RECIPES_INLINE_DEPRECATED = r"""
# NGSolve Implementation Recipes for IH SIBC

Production-grade NGSolve code patterns extracted from the lab's
`calc_inductance.py`, `calc_fem_kelvin.py`, `bem_sibc_solver.py`,
and the supporting `radia.sparsesolv_ngsolve` (Compact AMS / COCR).

Use these as starting points when extending the IH pipeline. All
code follows the lab's policies:
- Verify-First Policy (FES check before physics solve)
- Compact HX preconditioner (HYPRE-free, TaskManager-native)
- Shifted preconditioner for air + conductor problems
- No fallback chains (fail-fast)

## Recipe 1: Scalar BIE + SIBC Robin term with COCR + ComplexCompactAMS

Production reference: `src/radia/bem_sibc_solver.py::ScalarBIESIBCSolver`.

```python
from ngsolve import (Mesh, H1, BilinearForm, LinearForm, GridFunction,
                     grad, ds, dx, InnerProduct, specialcf)
from ngsolve.bem import LaplaceDL, LaplaceSL, SingleLayerPotential
import radia.sparsesolv_ngsolve as ssn

# Surface H1 space (BND only)
fes = H1(mesh, order=1, definedon=mesh.Boundaries("wp_bnd"))
u, v = fes.TnT()

# System: (1/2 M - DL + gamma * SL * M^-1 * K) phi = M phi_inc
# gamma = Z_s / (jw * mu0)  -- complex
gamma = Zs / (1j * omega * mu0)

M = BilinearForm(fes); M += u*v*ds; M.Assemble()
K = BilinearForm(fes); K += grad(u).Trace() * grad(v).Trace() * ds; K.Assemble()
DL = LaplaceDL(fes, fes).Assemble()
SL = LaplaceSL(fes, fes).Assemble()

# Build A = 0.5 * M - DL + gamma * SL * M^-1 * K (Lagrange mult for gauge)
# ... assembly details in bem_sibc_solver.py

# COCR solver (Sogabe-Zhang 2007) for complex symmetric system
solver = ssn.COCRSolver(
    matrix=A, preconditioner=prec,
    tol=1e-7, maxiter=400
)
phi.vec.data = solver * f.vec
```

**When to use**: BEM Scalar BIE (no FEM volume), nonlinear ESIM
workpieces via Karl outer iteration. Reference: lab 2026-04-13
Smythe sphere validation (-1.6% accuracy).

## Recipe 2: Per-node Z_s array via per-panel curvature extraction

Production reference: `src/radia/panels/calc_inductance.py::
_compute_panel_local_radii` and `_build_per_node_Zs`.

```python
import numpy as np
from ngsolve import H1, GridFunction, SurfaceL2

def per_node_curvature_h1(mesh, surface_label, percentile=10):
    # Discrete normal-angle radius extractor on a workpiece surface.
    # Returns one R_local per H1 vertex on `surface_label`. Used to
    # compute per-node Z_s for the SIBC Robin block.
    fes_p1 = H1(mesh, order=1, definedon=mesh.Boundaries(surface_label))
    ndof = fes_p1.ndof
    R_local = np.zeros(ndof)

    # Build panel adjacency from BND elements
    panels = list(mesh.Boundaries(surface_label).Elements())
    # ... per-panel normal, centroid, neighbor pairs (see calc_inductance.py)

    # For each panel: angle_ij = arccos(n_i . n_j), dist_ij = |c_j - c_i|
    # R_ij = dist_ij / angle_ij ; R_panel[i] = percentile(R_ij, 10)
    R_panel = compute_per_panel_R(panels, percentile=percentile)

    # Adaptive sliver clamp: R_floor[i] = 0.5 * panel_diameter[i]
    R_panel = np.maximum(R_panel, 0.5 * panel_diameters)

    # Vertex averaging to project panel -> H1 nodes
    for vid in range(ndof):
        neighbor_panels = vertex_to_panels[vid]
        R_local[vid] = np.mean([R_panel[p] for p in neighbor_panels])

    return R_local

# Use per-node Z_s in ScalarBIESIBCSolver
R_local = per_node_curvature_h1(mesh, "wp_bnd")
Zs_per_node = (1 + 1j) / (sigma * skin_depth(omega, mu0*mur, sigma)) * \
              (1 + (1 + 1j)/2 * skin_depth(omega, mu0*mur, sigma) / R_local)
solver.solve(phi_inc, Z_s=Zs_per_node, omega=omega)   # ndarray accepted
```

**Validation** (lab 2026-04-12, prolate spheroid `a/b=4`):
- Cu 10 Hz, scalar(mean R): -11.2% error
- Cu 10 Hz, per-node(mesh): +3.1% error (4x improvement)

## Recipe 3: FEM + Kelvin + Robin BC (hole approach)

Production reference: `src/radia/panels/calc_fem_kelvin.py`.

**KEY**: workpiece is SUBTRACTED from mesh (hole approach), NOT meshed.
SIBC = Robin BC on the hole boundary. Avoids the -34% systematic error
of the interface approach.

```python
from ngsolve import (Mesh, HCurl, Periodic, BilinearForm, LinearForm,
                     GridFunction, curl, dx, ds, InnerProduct,
                     specialcf, x, y, z, Preconditioner, BVP)
import radia.sparsesolv_ngsolve as ssn

mesh = Mesh("model.vol")  # hole approach + Kelvin pair

# Verify-First Policy: check FES and Kelvin pair BEFORE solving
print("Materials:", mesh.GetMaterials())
print("Boundaries:", mesh.GetBoundaries())

fes_base = HCurl(mesh, order=1, dirichlet="gnd",
                  complex=True, gradientdomains={})
fes = Periodic(fes_base)  # Kelvin pair = master/slave
slaved = sum(fes_base.FreeDofs()) - sum(fes.FreeDofs())
assert slaved > 0, "Kelvin Periodic identification failed"

# Functional test: Set 1 on kelvin_int -> ratio on kelvin_ext should be 1.0
gfu_test = GridFunction(fes); gfu_test.vec[:] = 0
gfu_test.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
ratio = Integrate(gfu_test*gfu_test, mesh,
                   definedon=mesh.Boundaries("kelvin_ext")) / \
        Integrate(gfu_test*gfu_test, mesh,
                   definedon=mesh.Boundaries("kelvin_int"))
assert abs(ratio - 1.0) < 1e-3, f"Kelvin ratio = {ratio}, expected 1.0"

# Build system
u, v = fes.TnT()
nu = 1.0 / (mu0 * mur_air)
nu_kelvin = nu_kelvin_cf(mesh)  # = nu0 * (r'/R)^2

a = BilinearForm(fes, symmetric=False)
a += nu * curl(u) * curl(v) * dx("air")
a += nu_kelvin * curl(u) * curl(v) * dx("kelvin")
a += (1j * omega / Zs) * InnerProduct(u.Trace(), v.Trace()) * \
     ds(definedon=mesh.Boundaries("sibc"))  # Robin SIBC

# Shifted preconditioner for air + conductor (eps = 1e-6 * nu)
a_shifted = BilinearForm(fes, symmetric=False)
a_shifted += nu * curl(u) * curl(v) * dx
a_shifted += 1e-6 * nu * u * v * dx   # eps * mass shift on prec only

c = Preconditioner(a_shifted, "compactams")  # Compact HX (HYPRE-free)

# RHS: Biot-Savart from PEEC filaments
f = LinearForm(fes)
f += InnerProduct(H_inc_cf, v.Trace()) * \
     ds(definedon=mesh.Boundaries("sibc"))

# Solve via BiCGStab + Compact HX
A.Assemble(); c.Update(); f.Assemble()
inv = ssn.BiCGStab(matrix=A, c=c, tol=1e-7, maxiter=200)
gfu.vec.data = inv * f.vec
```

**Validated**: 3D FEM-Kelvin L = 90.71 nH vs analytical torus 88.5 nH
(+2.5% on coarse mesh, lab 2026-04-12). 2D axisym: L < 1%, P < 2%.

## Recipe 4: Verify-First Policy on a typical IH FES setup

**MANDATORY** before any 5-minute physics solve. Each check is <1 s.

```python
from ngsolve import (H1, HCurl, Periodic, GridFunction, Integrate,
                     CoefficientFunction as CF)

mesh = Mesh("model.vol")

# Check 1: materials and boundaries spelled as expected
print("Materials:", mesh.GetMaterials())   # expect: air, kelvin (NOT wp!)
print("Boundaries:", mesh.GetBoundaries())  # expect: sibc, gnd,
                                            #         kelvin_int, kelvin_ext

# Check 2: Periodic actually constrains slave DOFs
fes_base = HCurl(mesh, order=1, dirichlet="gnd", complex=True)
fes = Periodic(fes_base)
slaved = sum(fes_base.FreeDofs()) - sum(fes.FreeDofs())
print(f"Kelvin slaved DOFs: {slaved} (should be > 0)")
assert slaved > 0

# Check 3: Functional Kelvin pair test
gfu = GridFunction(fes); gfu.vec[:] = 0
gfu.Set(CF((1.0, 0, 0)), definedon=mesh.Boundaries("kelvin_int"))
norm_int = Integrate(InnerProduct(gfu, gfu), mesh,
                      definedon=mesh.Boundaries("kelvin_int"))
norm_ext = Integrate(InnerProduct(gfu, gfu), mesh,
                      definedon=mesh.Boundaries("kelvin_ext"))
ratio = norm_ext / max(norm_int, 1e-30)
print(f"Kelvin int->ext ratio: {ratio:.6f} (should be 1.0)")
assert abs(ratio - 1.0) < 1e-3

# Check 4: SIBC boundary has nonzero area
A_sibc = Integrate(CF(1), mesh, BND,
                   definedon=mesh.Boundaries("sibc"))
print(f"SIBC area: {A_sibc:.6e} m^2 (should match coil-facing wp surf)")
assert A_sibc > 1e-9

# Check 5: GND vertex is identified
gnd_dofs = [i for i in range(fes_base.ndof) if not fes_base.FreeDofs()[i]]
print(f"GND-constrained DOFs: {len(gnd_dofs)} (should be > 0)")
assert len(gnd_dofs) > 0
```

If any check fails, FIX THE GEOMETRY/LABELS first. Do NOT iterate
solver runs to discover an FES-level bug.

Real-world cost saving (lab 2026-04-25):
- Wrong-debug path: 10 minutes physics solve x N iterations
- Right path: 5 FES checks at <1 s each + 1 physics solve

## Cross-MCP recipe references

- `radia_mcp.matrix_solvers` -- Krylov solver selection (COCR vs BiCGStab)
- `radia_mcp.bem` -- ngsolve.bem operator catalog (LaplaceDL/SL etc.)
- `radia_mcp.fem` -- FEM formulation theory (A-Omega, T-Omega, gauging)
- `radia_mcp.radia_ngsolve.analytical_formulas('validation_use_cases')`
  -- closed-form reference for each FEM/BEM analysis
"""


_BEM_DEPRECATED_MSG = """
# BEM knowledge has moved

BEM (Scalar BIE, EFIE, coupled BEM-SIBC) is no longer part of the IH
production path.  Use `mcp-server-radia-ngsolve` with topic `ngsbem_inductance`
for BEM documentation.

BEM solver modules are importable as `radia.bem_inductance`,
`radia.bem_coupled_solver`, and `radia.ngsbem_*`.  Runnable reference
scripts are in `validation_test/induction_heating/bem_reference/`.
"""


def get_ih_sibc_documentation(topic="all"):
    """Return IH SIBC documentation by topic.

    Production topics: peec_fem, esim, screening.
    Deprecated (BEM): overview, biot_savart -> redirects to radia-ngsolve.
    """
    production_topics = {
        "peec_fem": IH_PEEC_FEM,
        "esim": IH_ESIM,
        "screening": IH_SCREENING,
        "mathematica_verification": IH_MATHEMATICA_VERIFICATION,
        "ngsolve_recipes": IH_NGSOLVE_RECIPES,
    }
    # Legacy topics kept for backward compat but redirect to ngsolve MCP
    deprecated_topics = {"overview", "biot_savart"}

    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(production_topics.values())
    elif topic in production_topics:
        return production_topics[topic]
    elif topic in deprecated_topics:
        return _BEM_DEPRECATED_MSG
    else:
        return (
            f"Unknown topic: '{topic}'. "
            f"Available: all, {', '.join(production_topics.keys())}. "
            f"BEM topics moved to mcp-server-radia-ngsolve (ngsbem_inductance)."
        )
