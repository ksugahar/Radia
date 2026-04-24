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
1. Initial Z_s from ESIM at estimated H_t (e.g. 1 A/m, or from linear SIBC)
2. Solve outer BEM/FEM with current Z_s -> get H_t on conductor surface
3. Update Z_s from ESIM cell problem at the new H_t (per surface element)
4. Relaxation:  Z_s_new = (1 - alpha) * Z_s_old + alpha * Z_s_from_ESIM
                          (typical alpha = 0.5, lower for stiff problems)
5. Converge when ||dZ_s|| / ||Z_s|| < 1e-3 (typically 4-6 iterations for
   linear-ish steel, 10-20 for deep saturation)
```

## Convergence pitfalls
- **Diverging at low frequency**: linear SIBC initial guess is far from
  the saturated solution. Start with a single ESIM solve at average |H_t|.
- **Oscillating**: alpha too high. Drop to 0.3 or use Anderson acceleration.
- **Stuck at coarse error**: outer BEM/FEM mesh too coarse — local H_t
  spikes drive ESIM into very nonlinear regime. Refine the surface mesh
  or smooth Z_s spatially (per-node Z_s with vertex averaging).

## Per-element vs per-node Z_s
- Per-element (cheaper, FEM panel-quadrature compatible)
- Per-node (smoother, needed for high-curvature BEM, see
  `radia_mcp.radia_ngsolve.ngsbem_inductance` per-panel curvature SIBC)
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

## Coupling to BEM-SIBC

```python
# Per-surface-element Karl iteration loop
for it in range(max_iter):
    Z_s_per_elem = np.array([esim.solve(np.abs(H_t[i]))['Z']
                             for i in range(n_elem)])
    H_t_new = solve_bem_sibc(Z_s_per_elem)   # outer BEM solve
    if np.linalg.norm(H_t_new - H_t) / np.linalg.norm(H_t) < 1e-3:
        break
    H_t = (1 - alpha) * H_t + alpha * H_t_new
```

## Coupling to FEM-SIBC

The same Karl loop wraps the FEM solve.  See
`radia_mcp.ih.ih_sibc(topic='peec_fem')` for the production
PEEC + FEM induction-heating workflow that uses this exact pattern.
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
