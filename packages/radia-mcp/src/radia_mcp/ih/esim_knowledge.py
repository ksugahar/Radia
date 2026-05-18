"""
ESIM (Effective Surface Impedance Method) practical usage knowledge.

This module covers HOW TO USE ESIM in production -- CLI flags,
BH-curve file format, scalar-vs-per-element decision, Karl
iteration tuning, troubleshooting.  It is NOT paper material.

Source: docs/esim/USAGE.md (canonical) + relevant CLI / src code.
Authoring guide: every example in this module must be a copy-pastable
working command line; every flag listed must exist in the
src/radia/panels/calc_inductance.py argparse.
"""


ESIM_USAGE_OVERVIEW = """
# ESIM — when to use it

The 1-D nonlinear cell-problem solver returns Z_s(|H_t|) for a
ferromagnetic conductor with a BH curve.  ESIM **only** beats linear
SIBC in one regime: a ferromagnetic workpiece (typically steel /
electrical steel) driven hard enough that the surface H_t traverses
the BH knee.

## Decision table

| Workpiece                                      | Drive level         | Use |
|------------------------------------------------|---------------------|-----|
| Cu / Al / brass (mu_r ~ 1)                     | any                 | --impedance-model sibc (linear Dowell) |
| Steel, |H_t| stays well below knee (~ 100 A/m) | low                 | --impedance-model sibc --mu-r 100 (constant mu_r) |
| Steel, |H_t| straddles knee (~ 1 kA/m)         | mid                 | **--impedance-model esim --bh-file <bh.txt>** |
| Steel, deep saturation                         | high                | **--impedance-model esim --esim-per-panel** (BEM path only) |
| Lossy ferrite (Mn-Zn, Ni-Zn)                   | any                 | currently Python-API only (complex_mu kwarg); CLI does not yet expose |

## Cost vs benefit

| Path        | Linear SIBC | ESIM scalar | ESIM per-element |
|-------------|-------------|-------------|------------------|
| Outer solve | once        | 5-15 iter   | 5-60 iter (Karl) |
| Cost factor | 1           | 5-15x       | 5-60x + N_DOF cell calls |
| Accuracy    | depends on mu_r choice | good if |H_t| range narrow | needed if |H_t| varies > 3x across surface |

Rule of thumb: P_wp can change by up to 50 % moving from scalar to
per-element when surface H_t spans the BH knee.
"""


ESIM_USAGE_BH_FILE = """
# BH-curve file format

A BH curve is REQUIRED for --impedance-model esim.  Two-column
whitespace-separated text, ascending in H, including the origin:

  # comment lines start with #
  0.0        0.0
  10.0       0.00126
  50.0       0.00628
  100.0      0.0125
  500.0      0.061
  1000.0     0.12
  5000.0     0.55
  10000.0    1.20
  50000.0    1.85
  100000.0   1.95

Column 1: H [A/m]
Column 2: B [T]

Extra columns are ignored.

Sample BH curves ship under src/radia/panels/samples/, e.g.
  em_sample_bh.txt  (typical structural steel)

The cell solver interpolates with cubic spline (>= 4 points) or
linear (< 4 points); see esim_cell_problem.py::BHCurveInterpolator.
Out-of-range H is extrapolated.

For best results, include data points at H = 0, near the knee
(typical 500-2000 A/m for steel), and well past the knee
(>= 50000 A/m).
"""


ESIM_USAGE_INDUCTANCE_CLI = """
# calc_inductance.py with ESIM (BEM-SIBC weak coupling)

The PEEC-coil + BEM-workpiece path.  Workpiece BIE solved on a scalar
potential; SIBC enters as a complex Robin coefficient.

## Minimum invocation

  python src/radia/panels/calc_inductance.py \\
      --coil-solver peec --coil-step coil.step \\
      --vol workpiece.vol --wp-label sibc \\
      --sigma 2e6 --mu-r 100 --half-thickness 0.005 \\
      --frequency 100e3 --current 1.0 \\
      --impedance-model esim --bh-file em_sample_bh.txt

## ESIM-specific flags

| Flag                | Default | Meaning |
|---------------------|---------|---------|
| --impedance-model esim | sibc | Switch from linear Dowell to ESIM cell-problem |
| --bh-file <path>    | (none)  | Required when --impedance-model esim |
| --esim-max-iter     | 15      | Outer Karl iteration cap |
| --esim-tol          | 1e-3    | Convergence on max\\|dZ_s\\|/\\|Z_s\\| |
| --esim-relax        | 0.5     | Karl damping (under-relaxation); lower if oscillation |
| --esim-per-panel    | False   | Per-DOF Z_s mode (BEM-A path only; requires --wp-bem-backend intree-dense) |

## Workpiece flags

| Flag             | Required? | Meaning |
|------------------|-----------|---------|
| --vol            | yes       | Workpiece .vol from Cubit `radia_export netgen` |
| --wp-label sibc  | yes       | Workpiece boundary label (sideset) |
| --sigma          | yes       | Workpiece conductivity [S/m] |
| --mu-r           | yes       | Workpiece relative permeability (linear regime value) |
| --half-thickness | yes       | Cylinder radius for the 1-D cell problem [m] |

## PEEC coil flags

| Flag           | Default | Meaning |
|----------------|---------|---------|
| --coil-step    | (none)  | Coil STEP file (CAD) |
| --coil-solver peec | -- | Use PEEC perimeter filaments |
| --coil-sigma   | 5.8e7   | Coil conductivity [S/m] (Cu) |
| --peec-n-peri  | 16      | Perimeter filament count per cross-section |

## BEM workpiece backend

| Flag                | Default | Meaning |
|---------------------|---------|---------|
| --wp-bem-backend intree-dense | hacapk | Dense LU (required for --esim-per-panel) |
| --wp-bem-backend hacapk       | -- | ACA + GMRES (large meshes; scalar Z_s only) |
| --h1-order 1 or 2   | 1       | Basis order (P2 = curved Tri6 if .vol has curve_order 2) |
"""


ESIM_USAGE_FEM_KELVIN_CLI = """
# calc_fem_kelvin.py with ESIM (HCurl FEM workpiece + Kelvin)

Coil as PEEC filament line-integral RHS; workpiece as volumetric
HCurl A with Robin SIBC + Kelvin transformation for open boundary.

## Minimum invocation

  python src/radia/panels/calc_fem_kelvin.py \\
      --vol workpiece.vol --fes-order 1 \\
      --material custom --sigma 2e6 --mu-r 100 \\
      --half-thickness 0.005 \\
      --frequency 100e3 --current 1.0 \\
      --formulation total \\
      --impedance esim --bh-file em_sample_bh.txt \\
      --max-iter 15 \\
      --solver pardiso \\
      --peec-step coil.step --peec-sigma 5.8e7 \\
      --peec-n-peri 16 --peec-nwinc 3 --peec-nhinc 3 \\
      --require-kelvin

## ESIM-specific flags (note inconsistent names vs calc_inductance.py)

| Flag (calc_fem_kelvin.py)    | calc_inductance.py equivalent | Default |
|------------------------------|-------------------------------|---------|
| --impedance esim             | --impedance-model esim        | sibc    |
| --bh-file                    | (same)                        | none    |
| --max-iter                   | --esim-max-iter               | 15      |
| (no --esim-tol equivalent)   | --esim-tol                    | hard-coded 1e-3 |
| (per-DOF not in CLI yet)     | --esim-per-panel              | -- |

The flag-name inconsistency is a known wart -- both CLIs accept
their respective spellings unchanged.

## Solver choice

| --solver       | When |
|----------------|------|
| pardiso        | Default; fast sparse direct |
| ams            | HCurl p=1, low memory |
| bddc           | Preconditioned CG, p >= 2 recommended |
| iccg           | Generic fallback |
"""


ESIM_USAGE_FEM_COILMESH_CLI = """
# calc_fem_coilmesh.py with ESIM (full FEM A-V volumetric coil)

The highest-fidelity path: coil as volumetric mesh with H1 source
potential and HCurl A; workpiece SIBC as Robin BC.

## Minimum invocation

  python src/radia/panels/calc_fem_coilmesh.py \\
      --vol workpiece.vol \\
      --frequency 100e3 --current 1.0 \\
      --coil-sigma 5.8e7 --sigma 2e6 --mu-r 100 \\
      --half-thickness 0.005 \\
      --fes-order 1 \\
      --solver pardiso \\
      --sibc-bnd sibc --source-bnd source --sink-bnd sink \\
      --coil-mat coil \\
      --impedance-model esim --bh-file em_sample_bh.txt \\
      --esim-max-iter 15 --esim-tol 1e-3 \\
      --require-kelvin

## ESIM-specific flags

Same as calc_inductance.py (`--impedance-model esim`,
`--esim-max-iter`, `--esim-tol`, `--esim-relax`).  Per-element is
supported (--esim-per-panel) and operates on per-BND-DOF Z_s.

## Required workpiece + coil labels in the .vol

| Boundary  | Purpose |
|-----------|---------|
| sibc      | Workpiece surface (Robin BC) |
| source    | Coil source-cap face (Dirichlet phi=1) |
| sink      | Coil sink-cap face (Dirichlet phi=0) |
| kelvin_*  | Kelvin transformation pair (auto-detected) |

| Material  | Purpose |
|-----------|---------|
| coil      | Coil volume (sigma_coil applied here) |
| kelvin    | Kelvin sphere interior |
| air       | Air outside Kelvin |
"""


ESIM_USAGE_PER_ELEMENT = """
# Per-element (per-DOF) Z_s vs scalar Z_s

## What's the difference

| Aspect              | Scalar Karl (default) | Per-element Karl (--esim-per-panel) |
|---------------------|-----------------------|--------------------------------------|
| Z_s representation  | single complex value  | ndarray[ndof] complex per surface DOF |
| H_t source          | mesh-RMS of |H_t|     | per-DOF |H_t| from variational form |
| Cell solver calls per iter | 1               | N_DOF (typically 100-5000) |
| Captures saturation pattern? | NO          | YES |
| Convergence at default --esim-relax 0.5 | reliable | may fail in deep saturation |

## When per-element matters

Per-element pays off when |H_t| varies significantly across the
workpiece surface (typical IH: edge of workpiece sees 5-10x stronger
H_t than the far face).  Order-of-magnitude estimate:

- |H_t|_max / |H_t|_min < 2: scalar is fine; per-element is overkill
- 2 <= ratio < 5: per-element improves P_wp by 5-20 %
- ratio >= 5 AND ratio spans BH knee: scalar can under-estimate by 30-50 %

Headline benchmark (steel cylinder, 50 kHz, I = 100 A):
- |H_t| range: 250 - 2951 A/m (12x ratio, straddles BH knee at ~1 kA/m)
- Scalar P_wp = 30.6 W (converged in 15 iter)
- Per-element P_wp = 45.8 W (+48 %, NOT converged at 60 iter / alpha=0.3)

## When per-element doesn't converge

Per-element Karl with default --esim-relax 0.5 can fail to converge
to dZ < 1e-3 in 15 iter when surface DOFs straddle the BH knee.
Workarounds:

1. Lower --esim-relax (try 0.3 or 0.2) + raise --esim-max-iter (30-60)
2. Even at lower relax, may oscillate around 0.07-0.2 dZ_max.  P_wp
   typically stable to ~1 % despite dZ not reaching tol.
3. Anderson acceleration (planned, not yet implemented) is the
   intended fix for tight convergence.

## Backend restriction

--esim-per-panel requires --wp-bem-backend intree-dense (dense LU
on the BIE system matrix).  The HACApK ACA + GMRES backend does
NOT yet support per-DOF Z_s -- scalar only.  Adding HACApK per-DOF
is roadmap (v4.56+).
"""


ESIM_USAGE_CONVERGENCE = """
# Karl iteration convergence tuning

The outer Karl loop iterates on Z_s (scalar or per-DOF vector):

  Z_s^(k+1) = alpha * E(H_t(Z_s^(k))) + (1 - alpha) * Z_s^(k)

where alpha = --esim-relax and E is the cell-solver.  Picard
fixed-point with damping; converges when alpha * L < 1 where L is
the local Lipschitz of (E . H_t).

## Default

alpha = 0.5, max_iter = 15, tol = 1e-3.  Works for most engineering
IH cases.

## Symptom-driven tuning

| Symptom                                | Action |
|----------------------------------------|--------|
| dZ decreases monotonically in 5-8 iter | leave default |
| dZ oscillates after iter 1-2           | --esim-relax 0.3 |
| dZ slowly decreasing, > 15 iter needed | --esim-max-iter 30 |
| Deep saturation, dZ stalls > 0.1       | --esim-relax 0.2 + --esim-max-iter 60 |
| Even at low relax, oscillates          | Anderson accel (planned, not yet) |

## Reading the esim_history in JSON output

  "esim_history": [
    {"iteration": 0, "Z_s_abs": 0.0358, "H_t_rms": 247.3, "dZ": 1.0},
    {"iteration": 1, "Z_s_abs": 0.0352, "H_t_rms": 261.0, "dZ": 0.017},
    {"iteration": 2, "Z_s_abs": 0.0349, "H_t_rms": 268.4, "dZ": 0.008},
    ...
  ]

- dZ ratio between successive iter = empirical contraction factor;
  ~0.4 at alpha=0.5 (close to 0.5 = alpha) means L ~ 0.8 (good).
- If dZ stops decreasing or oscillates, L > 1/alpha (need lower alpha).

For per-element, the history has dZ_max instead of dZ:

  "esim_history": [
    {"iteration": 0, "Z_s_abs_mean": 0.027, "Z_s_abs_max": 0.031,
     "H_t_per_dof_mean": 514, "H_t_per_dof_max": 875, "dZ_max": 0.11},
    ...
  ]

dZ_max is max_i |Z_s_new[i] - Z_s_old[i]| / |Z_s_old[i]|.  If
dZ_max stops decreasing while Z_s_abs_mean stabilises, the global
average is converging but individual hotspot DOFs still oscillate.
For P_wp purposes this is usually OK (P_wp tracks the mean), but
local saturation pattern is unreliable.
"""


ESIM_USAGE_JSON_OUTPUT = """
# ESIM output JSON schema

calc_inductance.py / calc_fem_kelvin.py / calc_fem_coilmesh.py all
emit JSON to stdout.  ESIM-specific fields:

## Common fields

  "impedance_model": "esim",
  "esim_iterations": int,        # iter count to convergence (or max_iter)
  "esim_converged": bool,        # True iff dZ < tol reached
  "esim_history": [...],         # per-iter snapshots (see below)
  "Z_s_wp_real": float,          # final Re(Z_s), scalar or mean (per-element)
  "Z_s_wp_imag": float,          # final Im(Z_s)
  "skin_depth_wp_mm": float,     # delta for the linear-mu limit (diagnostic)

## Scalar Karl history entry

  {
    "iteration": 0,
    "Z_s_abs": float,            # |Z_s| at this iter
    "H_t_rms": float,            # mesh-RMS H_t [A/m]
    "dZ": float,                 # relative change vs prev iter
    "t_solve": float             # wall time for the inner BIE/FEM solve [s]
  }

## Per-element Karl history entry (--esim-per-panel)

  {
    "iteration": 0,
    "Z_s_abs_mean": float,       # mean |Z_s[i]| across all DOFs
    "Z_s_abs_max": float,        # max
    "Z_s_abs_min": float,        # min
    "H_t_per_dof_mean": float,   # mean |H_t| across all DOFs
    "H_t_per_dof_max": float,    # max
    "dZ_max": float,             # max relative change across DOFs
    "t_solve": float
  }

## Per-element ndarray output

When --esim-per-panel is used, the final Z_s array is exposed:

  "esim_per_panel": true,
  "esim_per_panel_Z_s_real": [float, ...],  # length = N_DOF
  "esim_per_panel_Z_s_imag": [float, ...],

Use this for spatial-pattern visualisation (e.g. plot Z_s magnitude
on the workpiece surface to see the saturation map).
"""


ESIM_USAGE_TROUBLESHOOTING = """
# Common ESIM errors and fixes

| Error                                              | Cause                         | Fix |
|----------------------------------------------------|-------------------------------|-----|
| `--impedance-model esim requires --bh-file`        | ESIM requested without BH     | Pass --bh-file <path> |
| `BH curve empty / not monotone in H`               | Malformed BH file             | Verify two-column ASCII, ascending H, includes (0, 0) |
| `ESIM:NOT-CONVERGED after N iter`                  | Karl hit max_iter             | Raise --esim-max-iter; lower --esim-relax |
| `--esim-per-panel ... wp-bem-backend hacapk`       | per-DOF not on HACApK backend | Use --wp-bem-backend intree-dense |
| `cell solver SCIPY_AVAILABLE False`                | scipy missing                 | pip install scipy |
| Karl converges but P_wp wildly off ref             | wrong --half-thickness        | Use min(R_wp, H_wp/2) for solid bulk |
| `BIE iv overflow` / NaN seed                       | very high xi (R/delta > 100)  | Cell solver uses thin-skin fallback automatically (v4.46.1+); upgrade radia |
| Per-element runs but stagnates around dZ_max=0.3   | BH-knee-straddling DOFs       | Try --esim-relax 0.2; if still stuck, P_wp is usually stable to ~1% anyway |

## Sanity checks

Before suspecting ESIM:
1. Run with `--impedance-model sibc` (linear) first.  Same workflow,
   no BH file.  If linear works and ESIM doesn't, the BH file or
   --half-thickness is wrong.
2. Check `skin_depth_wp_mm` in the output.  Should be << workpiece
   thickness for SIBC to apply.
3. Check `H_t_rms_A_per_m`.  If much smaller than the BH knee
   (e.g. 100 A/m for steel), linear SIBC is enough -- ESIM is
   over-engineered.
4. Inspect `esim_history` dZ trajectory.  Monotone decrease = good.
   Oscillation = lower --esim-relax.
"""


TOPICS = {
    "all":                None,
    "overview":           ESIM_USAGE_OVERVIEW,
    "bh_file":            ESIM_USAGE_BH_FILE,
    "inductance_cli":     ESIM_USAGE_INDUCTANCE_CLI,
    "fem_kelvin_cli":     ESIM_USAGE_FEM_KELVIN_CLI,
    "fem_coilmesh_cli":   ESIM_USAGE_FEM_COILMESH_CLI,
    "per_element":        ESIM_USAGE_PER_ELEMENT,
    "convergence":        ESIM_USAGE_CONVERGENCE,
    "json_output":        ESIM_USAGE_JSON_OUTPUT,
    "troubleshooting":    ESIM_USAGE_TROUBLESHOOTING,
}


def get_ih_esim_documentation(topic: str = "all") -> str:
    """Return ESIM usage knowledge for the requested topic."""
    topic = (topic or "all").strip().lower()
    if topic == "all":
        return "\n\n".join(
            f"## {key}\n{value}" if value else ""
            for key, value in TOPICS.items() if key != "all"
        )
    if topic in TOPICS:
        return TOPICS[topic]
    return (
        f"Unknown topic: {topic!r}.\n\n"
        f"Available topics: {', '.join(k for k in TOPICS if k != 'all')}.\n"
        f"Pass topic='all' for the full ESIM usage guide."
    )
