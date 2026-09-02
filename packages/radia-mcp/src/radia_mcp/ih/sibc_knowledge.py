"""
SIBC / ESIM knowledge for the induction-heating MCP server.

The human production interface is the Induction Heating block in the single
Radia Simulink library. It exposes the following validated headless solver
families through explicit initialization/configuration:

  (a) PEEC inductance              -> calc_inductance.py --coil-solver peec
  (b) BEM-A inductance             -> calc_inductance.py --coil-solver bem-a
  (c) PEEC+BEM weak coupling       -> calc_inductance.py --coil-solver peec --vol <wp>
  (d) BEM-A+BEM weak coupling      -> calc_inductance.py --coil-solver bem-a --vol <wp>
  (e) PEEC + FEM wp + Kelvin       -> calc_fem_kelvin.py --formulation total
  (f) Full FEM A-V + wp SIBC       -> calc_fem_coilmesh.py
  (g-i) Thermal 3D/static, 3D/rotating, and 2D-axisym -> calc_heat*.py

The native time-domain path uses separate readable Level-2 MATLAB Eddy and
Thermal S-Functions over checked standalone radia_mex object handles. Cubit's
embedded PySide toolbar owns mesh export only; it is not an analysis panel.
"""

from ..common import load_prompt

IH_PEEC_FEM = """
# IH Architecture: PEEC+FEM (v4.6.0)

## Two Solver Paths

| Path | Coil | Workpiece | Backend selection |
|------|------|-----------|-------------------|
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

Solver options in the application configuration:
- Dense LU: fastest for N < 3000 (default)
- HACApK saddle-point: for N >= 3000
- PRIMA: for frequency sweep (rarely needed in IH)

Files:
- `calc_peec.py` (headless computation)
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
- `calc_fem_kelvin.py` (headless computation)
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
- **BH file required for ESIM**: Simulink block and calc scripts raise if
  `--bh-file` is empty when `--impedance-model esim` is selected
  (pre-validation in `IHDesignSpec.build_command`).
- **Karl iteration diverges with very high relax**: lower the relax flag
  for SUS430 / S50C above 30 kHz: `--esim-relax 0.3`.  (CLI default is
  0.5; matches `calc_fem_kelvin`'s historical setting.)
- **Single Z_s per workpiece is mesh-RMS averaged**.  For workpieces with
  large H_t variation (sharp corners, coil shadow effects), the current
  implementation under-resolves the spatial saturation pattern.
  Per-element ESIM is on the roadmap; track via the radia-mcp tool.
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

IH_NGSOLVE_RECIPES = load_prompt("ih", "ngsolve_recipes")

_BEM_DEPRECATED_MSG = """
# BEM knowledge has moved

BEM (Scalar BIE, EFIE, coupled BEM-SIBC) documentation is owned by
`mcp-server-radia-ngsolve` topic `ngsbem_inductance` rather than duplicated
here. BEM remains a selectable IH solver family where its validated contract
applies.

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
