# Stream-Function coil-design panel (`calc_streamfunction.py`)

The FE-direct stream-function designer is exposed as a Layer-3 PySide6 panel
(`src/radia/radia_streamfunction.py`) over a headless Layer-4 calc
(`src/radia/panels/calc_streamfunction.py`).  ONE `argparse` surface drives
both; the panel is a composite `ModePanel` (a coil-`.vol` Browse + a Method
combo + a `QStackedWidget` over three sub-panels).

This file documents the panel/calc workflow and the boundary-condition /
contour / flux features.  The math core (the `(ACA+)+TSVD` least-norm solver,
`RegularizedTSVD`, the folded-Tikhonov Pareto front) is in
[`regularization.md`](regularization.md) and [`api.md`](api.md); the
single-stroke chaining theory is in [`single_stroke.md`](single_stroke.md).

## Three modes (`--method`)

| mode | does | key output |
|------|------|-----------|
| `design` | `target -> A psi = B` (folded-Tikhonov `RegularizedTSVD`) | `homogeneity_rms`, `peak_J` |
| `pareto` | `(homogeneity, peak current density)` front via `--pareto-lever {alpha, linf, geometry}` | `front[]` |
| `manufacture` | iso-contours -> orientation-consistent turns -> field-aware single-stroke wire -> STEP (`--step-output`) / PEEC (`--peec`) | wire + `loops_homogeneity_rms` + `peec` |

## I/O

- `--coil-vol` — a **standalone 2D surface** `.vol`.  `psi` is an `H1`
  GridFunction on it (Setup B: `definedon=coil.Boundaries('.*')`,
  `grad(v).Trace()`, `* ds`).
- `--eval-vol` — the evaluation region (surface **or** volume).
- `--target-cf` — the target field as a CoefficientFunction expression of
  `x,y,z`.  Scalar `-> Bz` (e.g. `"x"` = Gx, `"1"` = uniform), 3-vector
  `"(Bx,By,Bz)" -> ` full `B`.

## Current-confinement boundary condition (`--confine {off, on, abe}`)

On a **finite** former the contours run off the edges; closing them with a rim
chord injects a spurious edge current (LAB short cylinder Gx: single-current
rms 0.54, 42/42 contours open).  Confine the current to the patch:

| value | BC | use |
|-------|----|-----|
| `off` | none | contours that close on their own (full-ring solenoid, long former) |
| `on`  | `psi = 0` on every boundary edge (`dirichlet_bbnd`) | single-feed gradient/shim; **breaks** solenoids (both ends forced equal) |
| `abe` | **Abe edge-equipotential** — one free constant per physical boundary component + one ground | closes contours on **any** former, works for gradient **and** solenoid |

`abe` is the canonical current-potential BC (M. Abe, IEEE Trans. Magn.,
*DUCAS*; Appendix eq. 6 `T = R·T_IN`, constraints A-1/A-3).  It is implemented
as a DOF-reduction matrix `R` (`off`/`on` are the column-select special case);
a physical boundary edge (borders ONE surface element) is told apart from a CAD
**seam** (borders two) by element adjacency.  LAB short cylinder Gx: `abe`
closes the contours (`n_open 0`), gives 6.7× better separate-loops
single-current than `on` (0.022 vs 0.149), and does **not** break uniform.

> **Caveat.** `abe` is best for DESIGN + SEPARATE-TURNS + generality; it is
> **not** automatically best for the one-wire SINGLE-STROKE -- its edge
> equipotential makes a contour hug the boundary, so the field-aware chain can
> connect it worse than `on`.  The best single-stroke choice is target-specific.

## Contour drawing = flux-line drawing

The `psi` iso-contours are drawn by the magnetic-flux-line rule: between two
adjacent lines flows a fixed amount (current for `psi`, flux for the vector
potential `A_z`) -- Abe's "between nodes `i,j` flows `T_i − T_j`".  So
equal-`psi`-interval contouring already gives wire density `~ |grad psi| = |K|`,
the same density rule the flux-line **bubble system**
(Hirahatake/Noguchi/Igarashi/Yamashita, bubble `r ~ 1/sqrt|B|`) enforces.

- **`--contour-sub N` (order-p contour).**  The default marches the contour on
  the vertex (order-1) `psi`, dropping an order-2/3 design's edge DOFs.
  `--contour-sub 3` subdivides each surface triangle 3×3 and evaluates the
  full-order `psi` via `mesh.GetTrafo(el)` + `gfu(trafo(ip))` (the element-trafo
  MeshPoint dodges the boundary-point-eval-returns-0 quirk) -- the FE analogue
  of the analytical flux-line trace inside an element.  LAB Gx order 2:
  separate-turn `loops_homogeneity` 1.3e-4 → 1.1e-4, smoother wires.
- **`--flux-plot out.png` / `--flux-plane {x,y,z}` (bubble-system flux view).**
  Renders the designed coil's actual `B` field as flux lines on a cut-plane,
  seeded by the bubble system (line density `~ |B|`) -- a physical check that
  the coil produces the intended field (the four-lobe Gx gradient saddle).
- **`--steps-plot out.png` (per-step manufacturing view).** A 2×2 3D figure
  showing the coil at each stage: (1) the equal-current iso-contours
  (`N = --nlevels` turns -- this is how you set the number of lines),
  (2) the single-stroke (一筆書き) wire, (3) the sheet-metal (板金) `--distort`
  distorted wire, (4) the wire WITH thickness (太さ, `--wire-diam`, swept with a
  twist-free parallel-transport frame) + distortion.

## Single-stroke chain (`--chain {field_aware, nn}`)

The "extra lines" in the single-stroke view are the inter-loop **connectors**
(rungs): chaining `N` contour loops into one conductor needs `N−1` bridges, and
they carry the series current so their stray field is real, not cosmetic.
`field_aware` (default) keeps that field small two ways:

- **Cut placement.** Each loop's entry/exit cut is chosen by coordinate descent
  to minimise the **full** one-current wire error
  `min_I ||I·(loops + connectors) − B||` (not `||connectors||` alone -- that is
  worse on open contours, because the connectors are routed to *cancel* the
  rim-chord residual, not merely to be short).
- **Visit order.** The same wire-error objective is also minimised over a small
  set of candidate orders -- a nearest-neighbour seed and a 2-opt-shortened
  variant (which untangles the long "jump to a far lobe and back" crossings) --
  keeping whichever the cut-opt drives lowest.  The 2-opt **shortens** the
  rungs, which helps most cases (LAB Gx: connectors max 372→289 mm, delivered
  single-stroke +19…+70 %), but a length-optimal reorder can break the rungs'
  symmetric stray-field cancellation and *hurt* some cases (the documented
  "shorter rungs ≠ better field" trap, see [`single_stroke.md`](single_stroke.md)).
  Selecting the lower-wire-error order makes it **guaranteed never worse than
  nearest-neighbour** while capturing the real gains.

`field_aware` is never worse than `nn`; with closed contours (use `--confine
abe`) it reaches the separate-turns floor with **no** `--distort`.  The optional
sheet-metal `--distort` (single-current 3D control-grid Gauss-Newton) remains an
extra lever on top.

## End-to-end validation vs an independent codebase

`examples/stream_function/verify_coil_field_independent.py` designs an
MRI-gradient-scale coil (cylinder r = 0.15 m, L = 0.5 m, DSV r = 0.05 m) and
checks the field **two ways** -- the numpy straight-segment Biot-Savart used in
the designer **and** Radia's C++ `rad.ObjFlmCur` + `rad.Fld` (a separate
codebase).  They agree to 8–11 digits (uniform 3.5e-11, Gx 1.1e-8); the
`abe`-confined Gx coil reaches **1.0 %** nonlinearity on the short former,
cross-validated.  See the [examples README](../../examples/stream_function/README.md).

## Tests

- `tests/panels/test_streamfunction_golden.py` — the calc golden band (design /
  pareto levers / manufacture / field-aware chain / confine / order-p contour +
  bubble flux / cross-codebase), run via subprocess.
- `tests/panels/test_streamfunction_panel_qt.py` — headless PySide6 behaviour
  (mode combo, per-mode widget isolation, choice combos incl. `--confine abe`,
  `build_command` roundtrip, save/restore).
