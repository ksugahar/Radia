"""ESIM (Effective Surface Impedance Method) — general knowledge.

ESIM is a 1D cell problem solved through conductor depth that returns a
nonlinear surface impedance Z_s(H_t) for use in BEM/FEM with surface
impedance boundary conditions (SIBC).

This module documents the GENERAL technique (cell problem mathematics,
Karl iteration, module API) without coupling to any specific application.
For application-specific use of ESIM (induction heating workpieces with
steel BH curves), see `radia_mcp.ih.ih_sibc(topic="esim")`.

Promoted from `radia_mcp.ih.sibc_knowledge.IH_ESIM` on 2026-04-24 — the
underlying technique is generally applicable to any nonlinear-magnetic
SIBC problem (induction heating, eddy-current brakes, magnetic shielding,
nonlinear core losses).
"""

OVERVIEW = """
# ESIM: Effective Surface Impedance Method (Overview)

ESIM extends linear SIBC to nonlinear magnetic materials by solving a 1D
cell problem through the conductor depth and returning a field-dependent
surface impedance Z_s(H_t).

## Linear SIBC (baseline; the limit ESIM extends)
```python
Z_s = (1 + 1j) * rho / delta
delta = sqrt(2 * rho / (omega * mu0 * mu_r))
```
Fixed Z_s, no iteration. Fast, accurate for non-magnetic conductors
(Cu, Al), inaccurate for steel/ferrite where mu depends on |H|.

## When to use ESIM
- Nonlinear soft magnetic conductor (steel, electrical steel, ferrite)
- Surface field amplitude varies enough to traverse the BH knee
- Skin depth « conductor thickness (otherwise solve full 2D/3D FEM)

## When linear SIBC is enough
- Cu / Al / Au workpieces (mu_r = 1)
- Steel at very high frequency where d/delta > 50 (deep saturation
  doesn't matter because field doesn't reach there anyway)
- Linear validation runs / sanity checks
"""

CELL_PROBLEM = """
# ESIM 1D Cell Problem

Solves a 1D BVP through the conductor depth (z = 0 surface, z = d center
for symmetric slab; z = R for cylinder centerline):

```
rho * d^2 H / dz^2 + j * omega * mu(|H|) * H = 0
```

## Boundary conditions
- z = 0:    H(0) = H_t        (surface tangential field, prescribed)
- z = d:    dH/dz(d) = 0      (zero derivative at center, symmetry)

## Material law
- `mu(|H|)` from BH curve table or analytical fit
- Complex mu allowed (mu = mu' - j mu''); covers grain eddy current
  losses + magnetic hysteresis losses inside the cell

## Output
- `Z_s(H_t) = E_t(0) / H_t(0)` — complex surface impedance
- `P_prime(H_t)` = 0.5 * Re(Z_s) * |H_t|^2 — surface power density [W/m^2]

## Geometries supported by `radia.esim_cell_problem`
- `'slab'`     — symmetric infinite plate, half-thickness d
- `'cylinder'` — solid cylinder, radius R
- `'finite_slab'` — slab with anti-symmetric BC at center (1-sided heating)

## Discretization
1D piecewise-linear FEM through the cell. Coarse mesh near surface (where
field is large), refined toward Newton iteration of the nonlinear system.
"""

KARL_ITERATION = """
# Karl Iteration (Picard relaxation for SIBC + ESIM)

ESIM gives Z_s as a function of H_t, but H_t is the unknown the BEM/FEM
solve produces.  Karl iteration breaks the chicken-and-egg with under-
relaxed fixed-point iteration:

```
1. Seed Z_s = esim.solve(H0)['Z']    # H0 = small estimate, e.g. 5 A/m
2. Solve outer BEM/FEM with current Z_s -> read mesh-RMS H_t
3. Z_s_new = esim.solve(H_t)['Z']
4. Z_s = relax * Z_s_new + (1 - relax) * Z_s_old   # relax = 0.5 default
5. dZ = |Z_s - Z_s_old| / |Z_s_old|
6. break if dZ < tol (tol = 1e-3 default; 5-10 iter typical for steel)
```

## CLI (v4.46+ production scripts)
- `--esim-max-iter N`   outer Karl cap (default 15)
- `--esim-tol T`        relative tolerance |dZ|/|Z| (default 1e-3)
- `--esim-relax R`      under-relaxation (default 0.5; lower if stiff)
- `--bh-file FILE`      required; 2-column [H[A/m] B[T]] table

`calc_fem_kelvin.py` uses `--max-iter` (legacy name) and currently
does NOT expose `--esim-tol`.  The other two Karl scripts
(`calc_inductance`, `calc_fem_coilmesh`) use the
new flag names.

## Convergence pitfalls
- **Diverging at low frequency**: linear SIBC initial guess is far
  from the saturated solution.  The seed `esim.solve(5.0)` at small
  H is preferred over copying the linear Dowell Z_s.
- **Oscillating**: relax too high.  Drop to `--esim-relax 0.3` for
  SUS430 / S45C above 30 kHz.  Anderson acceleration is on the
  roadmap.
- **max_iter=1 false-not-converged** (pre v4.46.1): the convergence
  check required `iteration > 0` which made any 1-iter run report
  `esim_converged = false`.  Fixed in v4.46.1 to accept iter 0 when
  the user explicitly asked for max_iter <= 1.
- **Stuck at coarse error**: outer BEM/FEM mesh too coarse - local
  H_t spikes drive ESIM into very nonlinear regime. Refine the
  surface mesh or upgrade to per-element Z_s.

## Local Z_s (current production) vs scalar legacy
- **Local production**: `LocalESIMSurfaceModel` evaluates
  `Z_s(f, |H_t(x)|)` at surface quadrature samples and
  `AssembleSurfaceImpedanceGram` projects those values against the
  surface-Omega basis.  `SolveLocalESIMSurfaceVIM` drives HCurl-only systems;
  `CoupledHDivHybridVIMSystem.solve_frequency_local_esim` places the same
  update around the full fixed-bulk-HDiv/HCurl solve.  Do not average these
  sample values into one modal diagonal.
- **LUT production path**: `BuildLocalESIMSurfaceLUT` stores a pickle-free
  frequency/field table.  Online updates interpolate it with zero cell solves,
  reject extrapolation, and verify a material/cell SHA-256 signature.
- **Scalar legacy**: one mesh-RMS field and one Z_s for the whole surface.
  It is useful only as a quick ablation because it under-resolves local
  saturation and edge-field variation.
- **Still open**: simultaneous ordinary bulk nonlinear B-H plus local ESIM,
  hysteretic/rotational surface state, and multidimensional corner cells.
"""

MODULE_API = """
# `radia.esim_cell_problem` — Module API

```python
from radia.esim_cell_problem import ESIMFiniteSlabSolver

esim = ESIMFiniteSlabSolver(
    half_thickness=R_wp,    # [m] cell depth
    bh_curve=BH_DATA,       # [(H, B), ...] table or callable mu(|H|)
    sigma=sigma,            # [S/m] electrical conductivity
    frequency=freq,         # [Hz] working frequency
    geometry='cylinder',    # 'slab' | 'cylinder' | 'finite_slab'
)
sol = esim.solve(H_t_rms)
Z_s = sol['Z']              # complex surface impedance [Ohm]
P_prime = sol['P_prime']    # surface power density [W/m^2]
H_profile = sol['H_z']      # H(z) profile inside the cell (debug)
```

## Coupling pattern (v4.46+ production scripts use this)

For reduced HCurl/HDiv-VIM use the current local API:

```python
cell = LocalESIMSurfaceModel(bh_curve=bh, sigma=sigma)
solution = mixed.solve_frequency_local_esim(
    cell,
    frequency_hz,
    mixed_galerkin_keep_blocks=("volume1", "surface"),
    mixed_galerkin_eliminate_blocks="volume",
)
assert solution.converged
assert solution.surface_impedance.diagnostics()["passive"]
```

This updates nonlinear skin impedance around a fixed bulk HDiv magnetic
operator and accepts exactly one physical excitation.  It is not yet a
simultaneous bulk nonlinear B-H solve.

```python
# Scalar Karl iteration loop (legacy/application compatibility)
esim = ESIMFiniteSlabSolver(half_thickness=R_wp, bh_curve=bh,
                            sigma=sigma, frequency=freq,
                            geometry='cylinder')
Z_s = complex(esim.solve(5.0, max_iter=5)['Z'])   # seed
history = []
for k in range(max_iter):
    res = solve_outer_BEM_or_FEM(Z_s=Z_s)         # bem.solve / a_bf re-assemble
    H_t_rms = float(res['H_t_rms'])               # mesh-RMS amplitude
    Z_s_old = Z_s
    Z_s_new = complex(esim.solve(max(H_t_rms, 1e-3))['Z'])
    Z_s = relax * Z_s_new + (1 - relax) * Z_s_old
    dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)
    history.append({"iteration": k, "Z_s_abs": abs(Z_s),
                    "H_t_rms": H_t_rms, "dZ": dZ})
    if dZ < tol and (k > 0 or max_iter <= 1):
        break
# One final outer solve at the converged Z_s so post-proc sees the
# matching residual.
final_res = solve_outer_BEM_or_FEM(Z_s=Z_s)
```

## Production sites (v4.46+)

| Outer model | Script (`src/radia/panels/`) | Karl loop location |
|---|---|---|
| Scalar BEM-SIBC (PEEC coil) | `calc_inductance.py` | `_solve_workpiece_weak_coupled` |
| Scalar BEM-SIBC (BEM-A coil) | `calc_inductance.py` (same) | (same) |
| FEM-HCurl with Robin (PEEC coil) | `calc_fem_kelvin.py` | `solve_fem_kelvin` |
| FEM A-V with Robin (FEM coil) | `calc_fem_coilmesh.py` | `solve_fem_coilmesh` |

All four return the same JSON Karl diagnostic schema:
`esim_iterations`, `esim_converged`, `esim_history` (list of per-iter
dicts), plus the final converged `Z_s_wp_real` / `Z_s_wp_imag`.

For application-specific guidance on choosing the right script (PEEC+BEM
vs Full FEM, half_thickness for non-cylindrical workpieces, etc.) call
`radia_mcp.ih.ih_sibc(topic='esim')`.
"""


EVRS_COUPLING = """
# Production HCurl EVRS + local ESIM-SIBC coupling

This is the current Radia production interpretation of the SIBC bridge.  It is
separate from the deprecated scalar Warburg closure found in older SIBC
documentation.

## Topology-aware space

Start from a high-order NGSolve `HCurl(p)` parent space and classify every
conductor face by its neighboring materials:

- conductor-air/exterior: eligible for surface-Omega/SIBC modes;
- conductor-conductor: retain conductive graph-cycle bridge modes;
- conductor-insulator: retain a non-SIBC blocking trace;
- ordinary conductor interior: compress by the Eddy-Visible Response Space
  (EVRS), then quotient directions carrying no independent `curl(T)` current.

Thus `p` controls parent-space admissibility near corners and edges, while the
retained EVRS rank controls reduced-model size.  An SIBC face is never
inferred merely from being a boundary face.

## Local ESIM surface operator

For a nonlinear magnetic conductor, solve the one-dimensional normal ESIM cell
at the local tangential-field magnitude and tabulate
`Z_s^ESIM(|H_t|,s)`.  Do not replace those values by their surface average.
For surface-current modes `K_m`, assemble

    [G_Gamma(Z_s)]_mn = int_Gamma Z_s(x,s) K_m(x) . K_n(x) dGamma.

Radia exposes this contract as `SurfaceImpedanceGram` and
`AssembleSurfaceImpedanceGram`.  The latter projects per-quadrature-sample
impedances and checks the surface basis and weights.  If
`Re Z_s(x,s) >= 0` pointwise, the Hermitian dissipative part of the Gram is
positive semidefinite.  Passing an unnamed array as if it were a modal diagonal
is rejected because it would hide the required surface projection.

`LocalESIMSurfaceModel` owns the validated BH curve and SciPy one-dimensional
cell controls. `SolveLocalESIMSurfaceVIM`, or
`TopologyAwareHybridVIM.solve_local_esim`, runs the fail-loud HCurl outer loop
for exactly one physical excitation.  For a fixed bulk HDiv magnetic operator,
`CoupledHDivHybridVIMSystem.solve_frequency_local_esim` runs the same update
around the complete magnetic/current mixed solve and returns
`CoupledHDivHCurlLocalESIMSolution`. The returned fields and
`SurfaceImpedanceGram` always belong to the same linearization point.

For design sweeps, `BuildLocalESIMSurfaceLUT` solves the cell grid offline and
`LocalESIMSurfaceLUT` persists `Z_s(f,|H_t|)` in a pickle-free NPZ archive.
Attach a loaded table with `LocalESIMSurfaceModel.with_lut`. The online Karl
iteration then reports `cell_solve_count = 0`. Interpolation is bilinear in
log-frequency and log-field for `log(|Z_s|)` and continuous impedance phase;
extrapolation is forbidden. A SHA-256 material/cell signature rejects reuse with
a different BH curve, conductivity, or numerical cell setup. Establish LUT
adequacy by field/frequency node refinement. Add state
dimensions, rather than forcing a two-dimensional LUT, for hysteresis,
temperature-dependent data not already represented by BH/sigma, rotational
magnetization, or multidimensional corner cells.
`ValidateLocalESIMSurfaceLUT` automates this check with direct cell solves at
field midpoints and at both stored and midpoint frequencies; choose the table
density from its reported maximum relative interpolation error.

The nonlinear solve uses an outer Karl iteration:

1. solve the reduced HCurl system for the current local surface Gram;
2. recover the solution-shaped `|H_t|` profile on conductor-air faces;
3. interpolate the converged ESIM cell table;
4. rebuild `G_Gamma(Z_s)` and repeat to the impedance-change tolerance.

The converged Gram is solved directly inside the fixed-operator
HDiv-MMM/HCurl coupling. Do not overclaim the scope:
`NgsolveBDMEddyBubbleVIM` currently builds the BDM-MMM response reduction from
scalar `mu_r`. A simultaneous constitutive iteration that updates an ordinary
nonlinear HDiv magnetic operator together with the local ESIM Gram is still an integration task.

## Exact bulk/surface Schur elimination

Partition the reduced operator into eliminated ordinary bulk coordinates `b`
and retained HDiv/bridge/SIBC coordinates `k`:

    A = [[A_bb, A_bk], [A_kb, A_kk]].

The frequency-specific trial and test maps are

    P_trial = [-A_bb^-1 A_bk; I],
    P_test  = [-A_bb^-T A_kb^T; I],
    P_test^T A P_trial = A_kk - A_kb A_bb^-1 A_bk.

For a nonzero eliminated-block right-hand side, reconstruction is affine:

    x = P_trial c + [A_bb^-1 f_b; 0],
    rhs_reduced = P_test^T f.

Omitting the particular lift loses the directly excited bulk field even though
the Schur matrix itself is correct.  `MixedGalerkinOrthogonalization.solve`
(a legacy identifier) implements the complete affine reconstruction.

## Recorded p=6 broadband gate

`validation_test/cln/evrs_esim_sibc_mixed_notched_p6.json` records an mdx run
on a re-entrant notched conductor from 1 kHz to 1 MHz:

- 3557 active parent-HCurl DoFs;
- Krylov depth 8: 10 parent-T responses -> 8 current-Gram bulk modes;
- 14 conductor-cycle modes + 3 exterior-SIBC modes;
- 25 retained eddy coordinates, 0.703% of the parent DoFs;
- local-ESIM port difference from depth 22: at most 4.95e-8;
- coupled bulk/surface port difference from the direct reduced solve: at most 6.44e-16;
- uniform ESIM vs local ESIM: 0.229%--0.655%;
- linear SIBC vs local ESIM: 6.98%--38.31%;
- the 4-by-96 LUT has maximum direct-midpoint relative error 8.17e-4;
- all 16 outer solves converged in 10--17 updates with passive local Grams;
- every online outer update used the LUT and performed zero cell solves.

The linear-to-ESIM difference is a constitutive-model comparison, not an error
estimate.  Use local ESIM-SIBC when the high-frequency layer is thin and locally
one-dimensional.  Retain volume EVRS/VIM for low-frequency penetration,
genuinely three-dimensional edge/corner current, or rotational/hysteretic state
outside the ESIM cell model.
"""


def get_esim_documentation(topic: str = "all") -> str:
    """Return ESIM documentation for the requested topic."""
    topics = {
        "overview": OVERVIEW,
        "cell_problem": CELL_PROBLEM,
        "karl_iteration": KARL_ITERATION,
        "module_api": MODULE_API,
        "evrs_coupling": EVRS_COUPLING,
    }
    if topic == "all":
        return "\n\n".join(topics[k] for k in ("overview", "cell_problem",
                                                "karl_iteration", "module_api",
                                                "evrs_coupling"))
    if topic not in topics:
        return (f"Unknown topic '{topic}'. Available: {', '.join(topics)}, "
                f"or 'all' for everything.")
    return topics[topic]
