# Stream-Function coil-design application (`calc_streamfunction.py`)

The FE-direct stream-function designer is exposed through the masked
`radia_simulink_library/Applications/Stream Function` block over the headless
`src/radia/panels/calc_streamfunction.py` calculation. One
`argparse`/`StreamFunctionDesignSpec` surface drives Python, MCP, and Simulink;
the block launches that contract on an explicit rising trigger and writes the
standard `run.log`, `solver_result.json`, and `result.json` artifacts.

Cubit exports the `.vol` inputs as a separate process. The block configuration
holds the coil/conductor `.vol` path and mode-specific inputs such as
`--eval-vol`. There is no Stream Function notebook workbench.

This file documents the application/calc workflow and the boundary-condition /
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

## MATLAB Optuna outer loop

Use `radia.stream.OptunaRunner` when geometry, manufacturing, regularization,
or other design settings must be searched outside the linear inverse. Each
Optuna trial writes an isolated application configuration and launches one
complete Stream Function analysis. The inner solve remains the production C++
ACA+ factorization with QR/TSVD recompression; Optuna does not replace or
reimplement it.

The masked `Stream Function Optuna` block runs the requested trial batch on one
rising trigger and reports completion, failure, Pareto-front, best-trial, and
elapsed-time signals. Python may be launched once per complete trial, never once
per Simulink time step. Keep `aca_eps` fixed while comparing physical designs.
Tune numerical tolerances, rank limits, and thread settings in a separate
calibration study constrained by accuracy and memory.

## I/O

- `--coil-vol` — a **standalone 2D surface** `.vol`.  `psi` is an `H1`
  GridFunction on it (Setup B: `definedon=coil.Boundaries('.*')`,
  `grad(v).Trace()`, `* ds`).  **Any surface, including multi-surface.**
  Because the FE-direct `psi` lives on `coil.Boundaries('.*')`, the same path
  designs a cylinder, a sphere, a plane, or a **multi-component** former with
  no special-casing — e.g. a **biplanar** coil (two parallel plates in ONE
  mesh).  The `abe` BC groups *each* disconnected component's boundary edges
  into its own free constant (so `ndof_free < ndof`), and the contours close
  on every component.  Locked by
  `validation_test/panels/test_streamfunction_golden.py::test_streamfunction_biplanar`
  (two 0.3 m plates at z = ±0.1 m, Gx target: design 1.8e-3, contours close on
  both plates).
- `--eval-vol` — the evaluation region (surface **or** volume).
- `--target-cf` — the target field as a CoefficientFunction expression of
  `x,y,z`.  Scalar `-> Bz` (e.g. `"x"` = Gx, `"1"` = uniform), 3-vector
  `"(Bx,By,Bz)" -> ` full `B`.
- `--target-harmonic` — the target as a **spherical harmonic** (see below).
  Give `--target-cf` **or** `--target-harmonic` (loud error on both).

## Spherical-harmonic target + field decomposition (`--target-harmonic`)

In a current-free region `Bz` is harmonic, so it expands in **real regular
solid harmonics** `R_l^m(x,y,z)` — homogeneous harmonic polynomials
(`∇²R = 0`).  These are the named MRI gradients/shims, and the designer speaks
them natively.

| target | `--target-harmonic` | polynomial |
|---|---|---|
| Gx, Gy, Gz | `X`, `Y`, `Z` | `x`, `y`, `z` |
| 2nd-order shims | `Z2`, `ZX`, `ZY`, `C2`, `S2` | `z²−(x²+y²)/2`, `xz`, `yz`, `x²−y²`, `xy` |
| 3rd-order shims | `Z3`, `Z2X`, `Z2Y`, `ZC2`, `ZS2`, `C3`, `S3` | … |
| 4th-order shims | `Z4`, `Z3X`, `Z3Y`, `Z2C2`, `Z2S2`, `ZC3`, `ZS3`, `C4`, `S4` | … |

The table spans **`l ≤ 4`** (`1+3+5+7+9 = 25` harmonics).  A name, an
`l=L,m=M` pair, a `(L,M)` pair, or a weighted sum is accepted — e.g. `Z2`
(pure 2nd-order shim), `X` (Gx), `Z2:1.0,Z:0.1` (Z2 with a Z offset),
`l=4,m=-4` (≡ `S4`).  The comma inside `l=L,m=M` / `(L,M)` is handled, so
those forms work alongside the comma-separated weighted-sum syntax.  It
generates the solid-harmonic `--target-cf` polynomial, so the whole pipeline
(design / pareto / manufacture / single-stroke) is unchanged.

**Achieved-field decomposition** (design mode, `result["harmonics"]`, depth
`--harmonic-lmax`, default 3, max 4): the designed `Bz` over the DSV is
least-squares decomposed into solid harmonics — the per-`(l,m)` field-RMS
**spectrum**, the LSQ **residual**, and (with a `--target-harmonic`) the
**purity** (target-harmonic field fraction — the standard gradient-coil
quality metric) and the largest **contaminant**.  One poly-string table is the
single source of truth for both target and analysis, so the round-trip is
exact.  Verified (gallery record `demo_shim_coil_purity.py`, order 2,
confine abe): Gx → pure `X` (purity 1.000, `Z2X` contaminant 7e-6); `Z2` →
pure `Z2`; the 4th-order `Z4` shim → purity 0.99983 with a named ~1 % `Z3`
contaminant (high-`l` shims are harder because an `l`-th harmonic's field
scales as `rˡ` over a fixed DSV).

## Active shielding — primary + shield coil (`--shield-vol`)

An actively-shielded gradient coil (Mansfield & Chapman 1986; Turner 1986)
adds a second, OUTER cylindrical surface (the **shield**) whose current cancels
the **stray field** outside the assembly — essential in MRI so the switching
gradients do not induce eddy currents in the cryostat.

- `--shield-vol SHIELD.vol` — the outer shield coil surface mesh.
- `--shield-eval-vol EXT.vol` — the **external region where the stray field
  must vanish** (a large annular shell OUTSIDE the shield).
- `--shield-weight w` *(default 1.0)* — weight on the stray-null constraint
  rows relative to the DSV target rows.

The primary and shield are designed JOINTLY: one stacked least-squares system
fits the target inside the DSV **and** nulls the field at the external sample
points, with a block-diagonal seminorm `diag(Sₚ, Sₛ)` — the `RegularizedTSVD`
machinery is unchanged, just a bigger system.  `design` mode only.

Reported (shielded run): `homogeneity_rms` (DSV), `stray_rms` (the HONEST
combined field at **independent** external measure points, relative to the DSV
target), `stray_fit_rms` (the circular constraint-point residual, info only),
and per-coil `peak_J_primary` / `peak_J_shield`.

**Coverage is everything (measured)**: the external null region must COVER the
exterior, not a thin mid-plane slice.  With a large full-length shell the
shielding is genuine — `demo_active_shield.py` measures **~86× (39 dB)** stray
reduction (100–1700× in the far field), DSV homogeneity preserved.  A too-small
region overfits locally (~4×) and makes the field WORSE at larger `z`; the
un-sampled gap between shield and null region stays unshielded (~2–3×).  Locked
by `test_streamfunction_active_shielding`.

## Design objective / regularizer (`--regularize {l2, h1, inductance}`)

The design solve is `min ‖Aψ−B‖² + α·ψᵀSψ` — the target field `Aψ=B` picks the
coil, the seminorm `S` picks *which* coil among the many that fit (the SF system
is hugely under-determined, `N ≫ M`).  `S` is the design objective:

| `--regularize` | `S` | objective | cost |
|----------------|-----|-----------|------|
| `l2` | `I` | min `‖ψ‖²` | sparse |
| `h1` *(default)* | surface-H1 stiffness | min `‖∇ψ‖²` — a **smoothness proxy** for current density | sparse |
| `inductance` | `μ₀ Cᵀ·SL·C` | **min ½ψᵀLψ — the physical magnetic stored energy** | dense (BEM) |

`inductance` is the canonical **minimum-stored-energy** gradient-coil objective
(Turner / Forbes target-field method): low stored energy = fast slew rate, the
defining gradient-coil figure of merit.  `L` is the true self-inductance of the
stream-function surface current `K = n̂×∇ψ`, assembled from the
`ngsolve.bem` Laplace single-layer operator `SL` (`C` is the discrete surface
rot mapping `H1 → HDivSurface`).  **Validated**: the BEM `L` matches the analytic
thin-torus inductance to −0.6 %, and `ψ=z` on the cylinder reproduces a solenoid
whose energy matches the Nagaoka coefficient (0.78 at `2R/ℓ=0.6`).  On the
canonical Gx case it gives the lowest stored energy **and** the lowest peak
current density of the three regularizers; the design reports `stored_energy_J`.
All three fold onto the **same** ACA+TSVD machinery (`RegularizedTSVD.
from_stiffness(base, S)`).

`inductance` is **dense** by nature — the inductance is a fully-coupled `1/r`
integral operator, so `L` is dense `N×N` (unlike the sparse `l2`/`h1`).  Use it
at moderate `N`; for very large meshes use `h1` (the sparse proxy) or a future
H-matrix `SL`.  Locked by
`validation_test/panels/test_streamfunction_golden.py::test_streamfunction_min_inductance`.

## IH-resonance design (`--target-inductance` / `--resonance-cap`)

A gradient coil wants **minimum** inductance (fast slew); an **induction-heating
work coil wants a SPECIFIC inductance** so the coil + tank capacitor resonate at
the inverter frequency, `L_target = 1/((2π f)² C)`.  The SAME inductance physics
serves both — minimise it (`--regularize inductance`) or TARGET it here.

The field design (`ψ`) is independent of the turn count; `L_coil` is set by the
turns (`~N²`).  So `--target-inductance <L_H>` (or `--resonance-cap <C_F>` with
`--peec-freq <f>`, which computes `L_target = 1/((2πf)²C)`) **searches `nlevels`
by bisection** — each candidate is a single-stroke → PEEC `L_coil` — for the
turn count whose coil resonates, then designs at that `nlevels` and reports
`resonance` (`nlevels`, achieved `L_coil`, `resonance_freq_Hz`, the achievable
`L_range_H`, and an `in_range` flag — an out-of-range `L_target` is **reported**,
not silently clamped).  LAB cylinder, `C = 22 nF`, `f = 200 kHz` → `L_target =
28.8 µH` → `nlevels = 13` → `L_coil = 30.3 µH`, coil resonates at **195 kHz**
(the integer-turn quantisation leaves a few %; the tank `C` trims the rest), and
the uniform field is still hit (4.5e-4).  Locked by
`validation_test/panels/test_streamfunction_golden.py::test_streamfunction_ih_resonance`;
connects the SF designer to the `radia-ih` work-coil pipeline.

### With `--greedy-turns` (the low-turn coil): report the required capacitor

When you ALSO ask for a **few-turn** coil with `--greedy-turns`, greedy owns the
turn count, so a resonance spec can no longer search `nlevels` — the two levers
both set the turns and would fight (`L_coil ~ N²` pins `N`).  Rather than
silently dropping the resonance target, the result reports
`required_cap_F = 1/((2π f)² L_coil)` — the tank capacitor that resonates the
**delivered** few-turn coil at `--peec-freq` — plus, for `--resonance-cap`, the
`resonance_freq_Hz` that the coil + the given capacitor actually reach (with its
relative error vs the operating frequency).  This is the **"few-turn uniform IH
coil → required capacitor"** answer: design the low-turn uniform coil, then size
the capacitor to `required_cap_F`; or drop `--greedy-turns` to let `nlevels` be
searched for `L_target`.  The `resonance` dict carries `mode = "from_greedy_coil"`
(vs `"search_nlevels"`).  Verified end-to-end (LAB cylinder + DSV, uniform `Bz`,
8 greedy turns, `47 nF` @ `200 kHz`): `L_coil = 25.3 µH`, `required_cap_F =
25.0 nF`, the `47 nF` cap resonates that coil at `146 kHz`.  Locked by
`validation_test/panels/test_streamfunction_golden.py::test_streamfunction_resonance_with_greedy_reports_cap`.

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
  (2) the single-stroke wire, (3) the sheet-metal `--distort` distorted wire,
  (4) the wire WITH thickness (`--wire-diam`, swept with a twist-free
  parallel-transport frame) + distortion.

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

`validation_test/stream_function/verify_coil_field_independent.py` designs an
MRI-gradient-scale coil (cylinder r = 0.15 m, L = 0.5 m, DSV r = 0.05 m) and
checks the field **two ways** -- the numpy straight-segment Biot-Savart used in
the designer **and** Radia's C++ `rad.ObjFlmCur` + `rad.Fld` (a separate
codebase).  They agree to 8–11 digits (uniform 3.5e-11, Gx 1.1e-8); the
`abe`-confined Gx coil reaches **1.0 %** nonlinearity on the short former,
cross-validated.  See the [demo ledger](examples.md).

## Scaling (large surface meshes)

The design solve is `min ‖Aψ−B‖² + α·ψᵀSψ`.  The regularisation seminorm `S`
(the surface-H1 stiffness) and the DOF-reduction matrix `R` (`Sᵢₙd = RᵀSR`) are
**naturally sparse** (the FE stiffness is ~2 % dense) and are kept sparse
end-to-end; `RegularizedTSVD.from_stiffness` factors `S` once (`splu`) and
back-solves the `k` ACA modes.  Densifying them (the historical `ToDense()`)
was an O(N²) time + memory wall — **13.3 s and 3.75 GB at N = 15 260 DOF**.
After the sparse fix the ACA fold is ~0.3 s and peak memory ~0.2 GB at the same
N; the design now scales **linearly** past 27 000 surface DOF (the remaining
cost is the Biot-Savart design-matrix assembly, itself linear in N).  Locked by
[`bench_sf_scaling.py`](../../validation_test/stream_function/bench_sf_scaling.py)
(+ committed `bench_sf_scaling.json` / `.png`).

## Validation

- `tests/test_simulink_application.py` and
  `tests/matlab/test_simulink_workflow.m` — the block runner, mask, and artifact
  contract.
- The stream-function golden validation lane locks the calc behavior (design /
  pareto levers / manufacture / field-aware chain / confine / order-p contour +
  bubble flux / cross-codebase), run via subprocess.
