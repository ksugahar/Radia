"""

SHOWCASE NOTEBOOK: docs/ih_esim_benchmark/esim_showcase.ipynb -- Bessel cross-check (live) + committed digest/envelope/per-DOF/Karl figures.
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

## Primary reference (method origin)

ESIM is from K. Hollaus, V. Hanser, and M. Schoebinger, "A Nonlinear
Effective Surface Impedance in a Magnetic Scalar Potential Formulation,"
IEEE Trans. Magn., 2025, doi:10.1109/TMAG.2025.3613932 (bib key Hollaus2025;
official bibtex author field: "Hollaus, Karl and Hanser, Valentin and
Sch\\"{o}binger, Markus").  The Sugahara-lab IGTE 2026 paper re-casts this
method into a scalar BIE-SIBC; Karl Hollaus (TU Wien) is a co-author.
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
| --vol            | yes       | Workpiece .vol from Cubit `export netgen` |
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
| (no per-DOF equivalent)      | --esim-per-panel              | off     |
| (no Anderson equivalent)     | --esim-anderson-m <depth>     | 0 (off) |

The flag-name inconsistency is a known wart -- both CLIs accept
their respective spellings unchanged.  `--esim-per-panel` and
`--esim-anderson-m` are calc_inductance.py-only because they were
added during the BEM+SIBC convergence work (v4.66.0+) and have no
calc_fem_kelvin.py analogue.

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

Headline benchmark (steel cylinder, 50 kHz, I = 100 A; sweep_v2,
v4.67.0+ with triangle-P1-gradient extractor + Anderson m=5):
- |H_t| range: ~250 - ~2900 A/m (10x+ ratio, straddles BH knee at ~1 kA/m)
- Scalar P_wp     = 30.51 W (converged in 6 iter)
- Per-element P_wp = 18.75 W (converged in 7 iter)
- per/scalar = 0.615 ==> per-element predicts 38.5 % LESS than scalar
- Equivalently: scalar OVER-estimates by a factor 1.63

Earlier v4.47.2 - v4.66.x values (P_per = 45.4 W, +48 % "scalar
under-estimate") came from a Galerkin localization gradient
extractor that mis-placed the saturation hot-spot.  The sign of the
disagreement FLIPPED with the v4.67.0 fix.  Always cite the
v4.67.0+ sweep_v2 numbers, never the older Galerkin-era headline.

## When per-element doesn't converge

Per-element Karl with default --esim-relax 0.5 can fail to converge
to dZ < 1e-3 in 15 iter when surface DOFs straddle the BH knee.
Workarounds (in order of escalating cost):

1. Lower --esim-relax (try 0.3 or 0.2) + raise --esim-max-iter (30-60)
2. Even at lower relax, may oscillate around 0.07-0.2 dZ_max.  P_wp
   typically stable to ~1 % despite dZ not reaching tol.
3. **Anderson-type-II acceleration (v4.66.0+)**: enable with
   ``--esim-anderson-m <depth>`` where depth is the history length
   (typical 3-5).  Safeguarded against step-clip violations via
   ``esim_anderson.py::AndersonAccelerator`` (records ``n_clips``
   and ``n_restarts`` to the JSON output).  Hits dZ < 1e-3 on
   problems that Picard couldn't close even at relax=0.2.

### Per-DOF |H_t| extraction (v4.67.0+)

Per-element Karl needs the tangential H field SAMPLED PER DOF, not
just an L²-projected scalar.  The current extractor uses a
triangle-gradient of phi at each DOF:

  H_t_per_dof[k] = -grad_S(phi)|_{vertex k}

implemented in ``bem_sibc_solver.py::extract_H_t_per_dof_grad`` and
wired into the calc_inductance.py scalar + per-panel Karl paths
(lines ~805, ~882).  This replaced the earlier "Laplacian sample"
(which suffered from corner-vertex bias on the closed surface; see
commit 630527d4 for the symptom).

## Backend restriction

--esim-per-panel requires --wp-bem-backend intree-dense (dense LU
on the BIE system matrix).  The HACApK ACA + GMRES backend does
NOT yet support per-DOF Z_s -- scalar only.  Adding HACApK per-DOF
is roadmap.
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
| Even at low relax, oscillates          | --esim-anderson-m 3-5 (safeguarded Anderson-II, v4.66.0+) |

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

## Operating-regime stall map (dense (I, f) sweep, IGTE 2026 ESIM digest)

Empirical convergence behaviour on the IGTE workpiece (Ø50 mm x
H25 mm steel, BH src/radia/panels/samples/em_sample_bh.txt) with
default Karl settings (--esim-max-iter 30 --esim-relax 0.5
--esim-anderson-m 5), 108-case sweep (9 I_port x 6 f, scalar +
per-element each):

  Region                              | Behaviour
  ------------------------------------|----------------------------
  Bulk of low / medium I (47 of 49)   | converged 6-13 iter (median 8)
  I = 200 A, f = 500 kHz              | converged 26 iter (outlier)
  I = 500 A, f = 10 kHz               | converged 25 iter (outlier, max-gap case)
  I = 500 A, f = 20 kHz               | stall (iter=30 cap, conv=False)
  I = 500 A, f = 50 kHz               | stall
  I = 500 A, f = 100 kHz              | stall
  I = 500 A, f = 200 kHz              | stall
  I = 500 A, f = 500 kHz              | stall (spurious +20 % gap)

Interpretation: the I = 500 A high-frequency band is the only
operating region where per-element Karl + Anderson-II currently
fails to converge in 30 iterations.  Raising --esim-max-iter
beyond 30 does NOT help in this band (the per-DOF Z_s oscillates
around a non-fixed-point limit cycle, not slowly drifting).  The
underlying cause is suspected to be the locally non-Lipschitz
piecewise BH at deep saturation crossings; needs further study.

For converged cases the per-element vs scalar gap reaches:
  - max -49.8 % at (500 A, 10 kHz)    (converged, this is the
                                       digest panel-(b) case)
  - max -47.4 % at (200 A, 100 kHz)   (converged, second-worst)
  - max +28   % at (1 A, 500 kHz)     (sign flips below BH knee)

Reproducing the full sweep: examples/ih_esim_benchmark/sweep_f_I.py
(restart-safe after 1795e078 -- transient NAS-import flickers
auto-recover; see ih_esim "troubleshooting" topic).

Source-of-truth data: examples/ih_esim_benchmark/sweep_data_dense/*.json
(108 per-case JSONs + sweep_results.json, committed `847259d2`).
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
| Per-element stalls at iter=30 in I=500 A f>=20 kHz band | high-current high-freq limit cycle | See `convergence` topic, "Operating-regime stall map" -- not fixable by raising max_iter |
| `{"error": "No module named 'radia.<X>'"}` in JSON, returncode 0 | transient NAS-share import flicker under heavy 100bangoki OCR load | radia >= 4.92.0: calc_main now exits non-zero on any error (commit `5b88b67f`).  Sweep wrappers see returncode != 0 and can retry.  examples/ih_esim_benchmark/sweep_f_I.py has auto-recover for cached error JSONs (commit `1795e078`). |

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


ESIM_USAGE_SCALAR_VS_VECTOR_BEM = """
# Why scalar BIE-SIBC + ESIM is the right combination

Production guidance: for an MQS induction-heating workpiece with
Leontovich surface impedance, the minimal *and* most accurate
variational triple is

  scalar magnetic potential phi  (1 DOF / surface node)
+ Lagrange P2 basis on Tri6      (O(h^3) approximation)
+ isoparametric curved Tri6      (O(h^3) geometric error)

and the three orders match.  This is what `calc_inductance.py
--wp-bem-backend intree-dense --h1-order 2` configures by default
when the .vol has curve_order >= 2.

## Decision: which workpiece formulation?

| Path | When to use | When NOT to use |
|------|-------------|-----------------|
| scalar BIE (path A) | **closed workpiece + closed air** + MQS + SIBC valid | open strips, transformer windings (Omega_air not simply connected) |
| vector BEM-A (path B) | when coil topology requires n_peri-free Biot-Savart | NEVER for the workpiece itself if scalar BIE applies |
| FEM HCurl + Kelvin (path C) | magnetic core or thin gap inside air (BEM cannot resolve) | pure IH (coil + workpiece + vacuum) -- wasteful 100x DOF |
| FEM A-V coilmesh (path D) | when P_coil needs to be resolved volumetrically | publication headline -- use as cross-check only |

## Why scalar BIE is uniquely well-matched to Tri6 P2

The 6 Tri6 nodes carry **geometry AND basis** with the same Lagrange-P2
shape functions:

- isoparametric mapping x(u,v) = sum phi_i^P2 * x_i
- scalar field            phi(u,v) = sum phi_i^P2 * phi_i

No Piola transformation needed (unlike RWG / RT_0 in vector BEM-A,
where contravariant Piola couples det J non-trivially to the
basis).  No `n × A` evaluation needed (unlike HCurl FEM Robin BC).

End-to-end error rate `||phi - phi_h||_2 + ||x - x_h||_inf = O(h^3)`
on smooth Gamma + smooth Z_s.  The rates **match**, so refining h
gives full O(h^3) gain.  Other discretisations:

| Choice | basis-order | geom-order | match? |
|--------|-------------|------------|--------|
| scalar BIE Tri6 P2     | O(h^3) | O(h^3) | YES |
| scalar BIE Tri3 P1     | O(h^2) | O(h^2) | yes (lower) |
| scalar BIE Tri6 P1     | O(h^2) | O(h^3) | mismatch -- wasted curving |
| vector BEM-A RT_0 Tri6 | O(h)   | O(h^3) | mismatch -- DOF count grows but order stuck at 1 |
| FEM HCurl, Curve(2)    | O(h^3) | O(h^3) | yes but **100x more DOFs** + Kelvin truncation |

## Per-element Z_s is natural under scalar BIE

The BIE system matrix is `0.5 M - DL + gamma * SL * M^{-1} * K`
where `gamma = j*omega/Z_s`.  Switching from scalar to per-DOF Z_s
amounts to **row-scaling**: `A[i,:] = 0.5 M[i,:] - DL[i,:] +
gamma[i] * (SL M^{-1} K)[i,:]`.  3-line code change in the
assembler, diagonal-only re-application per Karl iteration.

In contrast, in FEM HCurl + Robin (path C) the Robin coefficient
`gamma * u.Trace() * v.Trace()` lives on a 2-D BND surface FES and
must be re-assembled as a full bilinear-form term per Karl
iteration -- 100x more expensive at the same accuracy.

## Per-iteration Karl cost comparison

Production gapped-torus + steel-cylinder benchmark (50 kHz, 1 A):

| Path | Outer DOFs | Per-iter cost | Karl scaling |
|------|-----------|---------------|--------------|
| A (scalar BIE, dense LU)     | 166-580      | 0.2 s | diagonal re-scaling only |
| A (scalar BIE, HACApK)       | 5000+        | 3-5 s | diagonal re-scaling only (scalar only -- HACApK per-DOF roadmap) |
| C (FEM-Kelvin, pardiso)      | 76000        | 12 s  | full re-assembly of Robin term |
| D (FEM-coilmesh, pardiso)    | 87500        | 25 s  | full re-assembly of Robin + coil-sigma term |

The scalar BIE path is the recommended PUBLICATION path; FEM paths
serve as cross-checks ([three-path consistency at <2 %](file:///S:/Radia/01_GitHub/docs/esim/CROSS_VALIDATION.md)).

## Reviewer-question pre-cooked answers

Q: "Why not vector EFIE/MFIE for completeness?"
A: MQS air is curl-free; vector BEM solves 6x the DOFs of the
   underlying scalar physics.  Used only as coil-source variant
   when topology demands.

Q: "Why scalar BIE in 2026 -- isn't this classical?"
A: The contribution is not the BIE.  It is the demonstration that
   scalar BIE + curved Tri6 + per-element ESIM is the **uniquely
   well-matched** discretisation for nonlinear SIBC IH, with
   measurable downstream cost (50-100x cheaper per Karl iter than
   FEM paths at the same convergence).

Q: "Does this generalize beyond IH?"
A: To any MQS thin-skin eddy-current problem on closed
   workpieces -- transformer cores in 3D, electric-motor rotors,
   magnetic shielding.  The scalar BIE applies whenever
   Omega_air is simply connected.

## Full discussion

See [`docs/esim/SCALAR_BIE_VS_VECTOR_BEM.md`](file:///S:/Radia/01_GitHub/docs/esim/SCALAR_BIE_VS_VECTOR_BEM.md)
for the full path-A vs path-B/C/D analysis, error-rate proofs, and
the historical Sauter-Schwab / Calderon-calculus reference list.
"""


ESIM_USAGE_HEADLINE_NUMBERS = """
# Headline numerical results -- for paper-style citation

The IGTE 2026 paper's primary numerical claim (sweep_v2, v4.67.0+,
2026-05-22 LAB benchmark):

## Per-element vs scalar Z_s (steel cylinder, BH knee)

Inputs:
  - Workpiece: cylindrical steel, R=5 mm, H=10 mm, sigma = 2e6 S/m,
                mu_r(linear) = 100, BH curve
                src/radia/panels/samples/em_sample_bh.txt (CEFC 2020),
                half_thickness = 5 mm
  - Coil: PEEC filament (16 perimeter), src/radia/panels/samples/
          ih_fem_kelvin_demo_coil.step
  - Frequency: 50 kHz
  - I_port: 100 A (chosen to push surface H_t past the BH knee)
  - Mesh: src/radia/panels/samples/ih_bem_sample_p1.vol
          (2150 BND tris, 1077 vertices)
  - Karl: --esim-max-iter 30 --esim-anderson-m 5 --esim-relax 0.5

Results (sweep_v2 = IGTE Fig. 1 source data):
  Quantity                  | Scalar Karl     | Per-element Karl       | Delta
  --------------------------|-----------------|------------------------|----------------
  Convergence (Anderson m=5)| yes @ iter 6    | yes @ iter 7           | both OK
  P_wp [W]                  | 30.51           | 18.75                  | per/scalar=0.615
                            |                 |                        | (= -38.5 %)
  L_total [nH]              | 90.70           | 92.88                  | +2.4 %
  delta_L [nH]              | -9.97           | -7.79                  | -22 % (less neg.)
  H_t_rms [A/m]             | 680 (scalar)    | 518 mean (per-DOF int) | --
  mean|Z_s|                 | 2.25e-2         | 2.97e-2                | per-elem > scalar

Engineering interpretation: at this drive level, hot-spot DOFs sit
past the BH knee, where the curve gives LOWER local Z_s (mu_r drops
sharply on the falling side).  Per-element resolves this; scalar
averages it away and over-estimates dissipation.  This is opposite
to what the older Galerkin-localization extractor (v4.47.2 - v4.66.x)
reported as +48 %; that result was the artifact of
|H_t|^2 ~ phi_i (K phi)_i sampling the surface Laplacian instead of
the gradient norm.  Always cite v4.67.0+ sweep_v2 numbers above.

Reproduction (from any machine with radia >= 4.67.0 installed):

  python src/radia/panels/calc_inductance.py \\
    --coil-step src/radia/panels/samples/ih_fem_kelvin_demo_coil.step \\
    --coil-solver peec \\
    --vol src/radia/panels/samples/ih_bem_sample_p1.vol --wp-label sibc \\
    --sigma 2e6 --mu-r 100 --half-thickness 0.005 \\
    --frequency 50000 --current 100.0 --coil-sigma 5.8e7 \\
    --impedance-model esim --bh-file src/radia/panels/samples/em_sample_bh.txt \\
    --esim-max-iter 30 --esim-tol 1e-3 --esim-relax 0.5 --esim-anderson-m 5 \\
    [--esim-per-panel] \\
    --h1-order 1 --wp-bem-backend intree-dense \\
    --output {scalar,per_panel}.json

For the full 108-case (I, f) sweep map (Fig. 1a of the IGTE 2026
digest) use examples/ih_esim_benchmark/sweep_f_I.py.  Grid:
I_port in {1, 2, 5, 10, 20, 50, 100, 200, 500} A x
f in {10, 20, 50, 100, 200, 500} kHz x {scalar, per_panel}
= 9 x 6 x 2 = 108 cases.

Frozen artifacts: examples/ih_esim_benchmark/sweep_data_dense/*.json
(108 per-case JSONs + sweep_results.json + side-wall |Z_s| at the
max-gap I=500 A / f=10 kHz case, committed `847259d2` on 2026-05-30).
Replaces the older sparse 32-case `sweep_v2` at C:/temp/igte_bench/.

## Dense-sweep max-gap point (IGTE 2026 digest panel (b))

Maximum per-element-vs-uniform P_wp gap over the dense grid:

  Case            Scalar P_wp    Per-element P_wp   Gap       Notes
  --------------  -------------  -----------------  --------  -----------------
  500 A, 10 kHz   162.85 W       81.83 W            -49.75 %  CONVERGED iter=6
                                                              (digest panel (b))
  200 A, 100 kHz  310.48 W       163.36 W           -47.39 %  CONVERGED iter=8

The 500 A / 10 kHz case has both the maximum P_wp gap AND the
maximum |Z_s| spatial spread (4.5--16.5 mOhm = 4.36x ratio) among
all CONVERGED cases.  Non-converged 500 A high-freq cases formally
show larger spreads (up to 9.3x at 500 A / 200 kHz) but the per-
element Karl limit cycle makes those Z_s values unreliable.

Source: examples/ih_esim_benchmark/sweep_data_dense/
        I500_f10k_per_panel.json (per-DOF Z_s) +
        examples/ih_esim_benchmark/sweep_data/
        I500_f10k_Zs_side_field.json (side-wall (theta, z, |Z_s|, R)).

## Three-path consistency (linear-mu screening)

Same geometry but I_port = 1 A (linear-mu regime):

  f [kHz]  PEEC-BEM  FEM-Kelvin  FEM-coilmesh  |Z_s| agreement
  ----     --------  ----------  ------------  ---------------
  10       1.103e-2  1.101e-2    1.097e-2      < 1 %
  50       2.572e-2  2.543e-2    2.321e-2      < 2 % (paths AB), -10 % (D)
  100      3.770e-2  3.722e-2    3.200e-2      < 2 % (paths AB), -18 % (D)
  500      1.314e-1  1.293e-1    8.55e-2       < 2 % (paths AB), -53 % (D)

FEM-coilmesh degradation is coil-mesh under-resolution
(delta_coil / h_coil), not a Karl-loop bug.  PEEC-BEM <-> FEM-Kelvin
remains <2 % across the full 10-500 kHz sweep.

## Bessel reference (cell-problem validation, linear-mu)

  xi = R/delta   Frequency   Re(Z_s)_anal   Re(Z_s)_num   rel.err
  4              1 kHz       7.50e-3        7.50e-3       < 1e-5
  14             10 kHz      2.49e-2        2.49e-2       < 1e-4
  45             100 kHz     7.95e-2        7.95e-2       4e-4
  140            1 MHz       2.50e-1        2.50e-1       1.3e-3

Solver: ESIMFiniteSlabSolver(geometry='cylinder', mu_r=100, n_nodes=2000).
Bessel reference: scipy.special.iv (modified Bessel I_0, I_1).
Reproducer: examples/ih_esim_benchmark/analytical_bessel_baseline.py.

Cite these tables verbatim in any paper / talk; values are
locked-in for radia >= 4.55.3.
"""


ESIM_USAGE_EM_TABLE_COUPLING = """
sigma(T) / mu(T) thermal-EM coupling via a precomputed Z_s(|H_t|, T)
table -- the engineering alternative to per-timestep Karl-FEM when the
EM problem is quasi-static on the heat timescale (EM 1/omega ~ 100 us
at 10 kHz vs heat tau ~ s -> ratio ~1e4).  Two CLI scripts (NOT yet
wired into the radia-ih panel).

## Step 1: build the table -- calc_em_table.py

Solves the 1-D ESIM cell problem at each node of a log-spaced |H_t|
grid x linear T grid, writes a .npz:

    python calc_em_table.py \\
        --material steel --frequency 50e3 \\
        --H-grid-min 1e2 --H-grid-max 1e5 --n-H 20 \\
        --T-grid-min 20  --T-grid-max 800 --n-T 20 \\
        --sigma-T-curve sigma_T_steel.csv \\
        --curie-temp 770 --curie-exp 0.5 \\
        --em-table-output em_table_steel_50kHz.npz

T-dependence flags (both OPT-IN; honest-uncertainty stance):
  --sigma-T-curve CSV   2-col (T_celsius, sigma_S_per_m).  Clamped
                        outside its T range (NO extrapolation).  Omit
                        -> material preset constant sigma at all T.
  --curie-temp Tc       Curie mu(T) collapse: scales the BH curve's
                        MAGNETIC part by f(T) = (1 - T/Tc)^beta,
                        clamped [0,1], so mu_r -> 1 at/above Tc.
                        Default 0 = BH T-invariant.
  --curie-exp beta      collapse exponent (default 0.5).

These are USER-CHOSEN functional forms (CSV curve, (Tc, beta)) --
T-dependent data is inherently uncertain, so it is supplied + clamped +
recorded in the table meta, never hardcoded false certainty (CLAUDE.md
"T-dependence: functionalize ...").  T-dependent BH SHAPE beyond the
Curie mu-collapse is a known gap.

.npz schema: H_grid, T_grid, Zs_re, Zs_im, q_surf (= 0.5 Re(Z_s) |H_t|^2),
meta (material, frequency, sigma_T_curve, curie_temp, curie_exp, ...).

## Step 2: runtime heat solve -- calc_heat_with_em_table.py

Backward-Euler heat; per timestep scales |H_t_ref(r)| by I(t)/I_ref
and bilinearly interpolates q_surf from the table at every surface DOF.
Dominant simplification: the SPATIAL shape of |H_t(r)| is frozen at the
reference solve; only the amplitude follows I(t).

Reference field --ht-source:
  kelvin (recommended)  --ht-sol <stem>_Jsurf.sol + --em-vol
                        <stem>_fem.vol from a calc_fem_kelvin run.
                        |J_s| = |H_t| on the SIBC face -- captures
                        workpiece backreaction.  One upstream FEM solve.
  biot                  --coil-step <coil.step>.  Walks the coil
                        centerline (same extractor as PEEC), evaluates
                        the incident Biot-Savart field at each surface
                        DOF, projects onto the local tangent plane.  NO
                        backreaction -> ~2x low for a good flat
                        conductor (the eddy-current image that doubles
                        H_t is absent).  --biot-image-factor 2.0 opts
                        into that doubling EXPLICITLY (never silent;
                        No-Fallbacks).  Needs no EM solve.

    python calc_heat_with_em_table.py \\
        --wp-vol workpiece.vol --em-table em_table_steel_50kHz.npz \\
        --ht-source kelvin --ht-sol run_Jsurf.sol --em-vol run_fem.vol \\
        --I-ref 1000 --coil-current 1200 --dt 0.5 --t-end 60

--coil-current-csv <t_s,I_A> gives a time-varying current (clamped).

## Tests
  tests/panels/test_em_table_curie.py       Curie mu(T) helpers
  tests/panels/test_heat_em_table_biot.py   biot |H_t| tangential
                                            projection + image-factor
                                            vs analytical circular loop

## When to prefer this over --impedance esim (Karl)
  - You need the WHOLE heat-up transient (many timesteps), not one EM
    operating point.
  - sigma(T) / mu(T) move enough during heating that a single Z_s is
    wrong, but re-solving the EM per step is too expensive.
  - The |H_t| spatial shape is roughly current-independent (true in the
    weak-redistribution regime; re-run the kelvin reference at a
    representative T if a Curie band sweeps across the part).
"""


TOPICS = {
    "all":                  None,
    "overview":             ESIM_USAGE_OVERVIEW,
    "em_table_coupling":    ESIM_USAGE_EM_TABLE_COUPLING,
    "bh_file":              ESIM_USAGE_BH_FILE,
    "inductance_cli":       ESIM_USAGE_INDUCTANCE_CLI,
    "fem_kelvin_cli":       ESIM_USAGE_FEM_KELVIN_CLI,
    "fem_coilmesh_cli":     ESIM_USAGE_FEM_COILMESH_CLI,
    "per_element":          ESIM_USAGE_PER_ELEMENT,
    "convergence":          ESIM_USAGE_CONVERGENCE,
    "json_output":          ESIM_USAGE_JSON_OUTPUT,
    "troubleshooting":      ESIM_USAGE_TROUBLESHOOTING,
    "scalar_vs_vector_bem": ESIM_USAGE_SCALAR_VS_VECTOR_BEM,
    "headline_numbers":     ESIM_USAGE_HEADLINE_NUMBERS,
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
