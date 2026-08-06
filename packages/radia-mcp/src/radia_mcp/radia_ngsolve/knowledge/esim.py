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


def get_esim_documentation(topic: str = "all") -> str:
    """Return ESIM documentation for the requested topic."""
    topics = {
        "overview": OVERVIEW,
        "cell_problem": CELL_PROBLEM,
        "karl_iteration": KARL_ITERATION,
        "module_api": MODULE_API,
    }
    if topic == "all":
        return "\n\n".join(topics[k] for k in ("overview", "cell_problem",
                                                "karl_iteration", "module_api"))
    if topic not in topics:
        return (f"Unknown topic '{topic}'. Available: {', '.join(topics)}, "
                f"or 'all' for everything.")
    return topics[topic]
