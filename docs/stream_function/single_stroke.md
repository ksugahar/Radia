# Single-stroke chain construction

A real wound coil must be **one continuous conductor** driven by a
single current source.  The iso-contours of the stream function are
*N closed loops* — to manufacture, they must be connected into one
chain.  The connection segments carry the full current and produce a
**parasitic field** that was not part of the original SF design.

This page documents the three single-stroke methods shipped + the
empirical complexity tier framework that bounds their performance.

## The chain methods

All are in
[`demo_sf_to_peec_gx.py`](demo_sf_to_peec_gx.py),
selectable via
`--chain-method {field_aware, kuijpers, lobe, greedy, nn_blend}`.

> **Use the [`single-stroke-chain`](../../.claude/skills/single-stroke-chain/SKILL.md)
> skill** when joining a new coil's contours.  The single-stroke
> connection has no clean closed-form optimum — it is a reason-and-verify
> task (pick a method that respects the current-sign order → build the
> chain → measure DSV RMS → keep only if it beats the baseline → escalate
> to Path-A).  Five deterministic "improvements" have been tried and most
> made the field WORSE; the skill encodes the traps so you do not repeat
> them.

### `field_aware`: sign-order + azimuthal-min cuts *(default, recommended)*

Uses the `kuijpers` lobe/current-sign visiting ORDER (the dominant factor),
but instead of snapping every cut to a fixed φ, chooses each contour's cut
by coordinate descent to **minimise the azimuthal arc to its chain
neighbours** (axial `dz` is free — it carries no stray field).  This lets
the rungs' stray fields cancel more symmetrically over the DSV.

**Pros**: best DSV RMS found; beats `kuijpers` in every tested config.
**Cons**: slightly more, slightly scattered rungs than `kuijpers`.

Measured (default nphi=24 nz=40 nlevels=12): RMS **9.29 %**, x-axis
nonlin **7.20 %** — vs `kuijpers` 16.24 % / 9.73 %.  Robust: 30–54 % lower
RMS than `kuijpers` across nlevels 10/12/16 and nphi32/nz48.

**Key subtlety**: the azimuthal rung *total* is NOT the predictor —
`field_aware` has MORE azimuthal arc than `kuijpers` (12522 vs 7444 mm) yet
LOWER RMS, because field impact is the symmetric *cancellation* of the rung
stray fields, not their summed length.  And the geometry-only
`nn_blend` (12847 mm azimuthal, near-equal) gives RMS 0.65 because its
order interleaves the lobe signs.  The lever is the sign-respecting ORDER,
then symmetric cut placement — neither reducible to a single scalar metric
(hence the skill).

### `greedy`: global nearest-neighbour

For each polyline, find the closest unvisited polyline (in periodic
`(φ, z)` distance) and connect via a cylinder-surface geodesic.  Each
polyline is opened at the point closest to the previous chain end.

**Pros**: shortest total wire length.
**Cons**: chain criss-crosses the cylinder, visually obvious "wasted"
arcs (the user's flag 2026-05-30: 「等高線同士を無駄に接続している」).

Measured (planar Gx fingerprint baseline): RMS 21.78 %, x-axis
nonlin 11.40 %, length 22.65 m, 1350 segments.

### `lobe`: 4-quadrant Maxwell-pair topology

Classify each closed contour into one of four saddle quadrants
`(sign(x_centroid), sign(z_centroid))`.  Within each quadrant, sort by
polyline arc length (outer first) and chain via greedy NN.  Inter-
quadrant transitions follow a fixed order `(+x,+z) → (+x,-z) → (-x,-z)
→ (-x,+z)` with alternating direction so the exit of one quadrant is
near the entry of the next.

**Pros**: only 3 inter-quadrant arcs (vs N−1 for greedy).
**Cons**: slightly longer than greedy because forced quadrant ordering.

Measured: RMS 23.98 %, x-axis nonlin 9.67 %, length 24.67 m, 1350
segments.

### `kuijpers`: per-lobe cut line with straight rungs *(prior best, robust fallback)*

After Kuijpers, Jansen, Lomonova (Compumag 2023 [525]) Method-1.  Each
lobe is given a **fixed cut φ** (`+x` lobes at φ=0, `-x` at φ=π).
Every contour in a lobe is opened at the point closest to its lobe's
cut φ.  Contours within a lobe are sorted by the z coordinate of their
cut point; adjacent contours are connected by **one straight (φ, z)
"rung" at the cut**.  Inter-lobe transitions use the same straight
blend.

This is the Fig.4 rung pattern of Kuijpers et al.

**Pros**: cleanest visual structure (no criss-cross), best DSV RMS, 23%
fewer segments (rungs are 1 segment vs 8-segment geodesic helix).
**Cons**: slightly longer wire than greedy.

Measured: RMS **16.24 %** (-25 % vs greedy, -32 % vs lobe),
nonlin 9.73 %, length 23.99 m, **1040 segments**.

Confirms the paper's observation "deviation region = field error
region" — minimising the deviation minimises the field error.

## Avoiding the long bridges entirely: multi-wire

For a multi-lobe coil (Gx fingerprint) the long inter-lobe "rungs" the
single stroke needs are **structurally irreducible** — one continuous
wire MUST traverse all four lobes, so it must bridge between them.  Two
attempts to remove them both FAILED (2026-05-31):

  - **Reorder (boustrophedon / snake)** so consecutive lobes meet at the
    saddle: did NOT shorten the worst rung (still 291 mm) and made the
    field WORSE (9.3 % → 19.2 %) — the rung arrangement's symmetric
    stray-cancellation breaks.
  - **Cut the chain at the long rungs** into sub-wires: produces OPEN
    current paths (current can't start/stop mid-air; `div J ≠ 0`), so the
    Biot-Savart field is unphysical (9.3 % → 21 %).

> **Reordering CAN help — but only if you keep the better order.** The FE-direct
> calc (`calc_streamfunction.py`, `_field_aware_chain`) DOES reorder: it
> 2-opt-shortens the visit order to untangle the long crossings, then runs the
> field cut-opt and keeps whichever of {nearest-neighbour, 2-opt} order gives the
> lower **full-wire-error** (`min_I ||I·(loops+connectors) − B||`).  Measured on
> the FE-direct cylinder Gx (`--confine abe`/`off`, nlevels 10/12/16): 2-opt
> helps 5/6 cases by +19…+70 %, but HURTS one (abe nl=16, −78 %) — *exactly* the
> stray-cancellation break above.  Selecting the lower-error order discards that
> one regression while keeping the five gains, so it is **guaranteed never worse
> than nearest-neighbour**.  The lesson is unchanged: shortening rungs is not the
> objective (field cancellation is) — but a reorder *evaluated against the field
> objective and kept only when it wins* is a safe upside, not a trap.

The only physical way to have NO bridges is for every independent wire to
be a **closed loop** — i.e. drive each closed contour as its own
conductor (multiple current feeds).  That is exactly
[`demo_coil_design_gx.py`](demo_coil_design_gx.py):
independent saddle-shaped closed loops, **no bridges, DSV RMS 0.81 %** —
an order of magnitude better than any single stroke.

| Design | wires / feeds | bridges | DSV RMS |
|--------|---------------|---------|---------|
| single-stroke (`demo_sf_to_peec_gx.py`, field_aware) | 1 | yes (irreducible) | 9.3 % |
| multi-wire (`demo_coil_design_gx.py`, independent loops) | N (one per loop) | **none** | **0.81 %** |

**RULE**: if the wasteful long connections are unacceptable, switch to
the independent-closed-loop (multi-wire) design — do NOT try to cut or
reorder a single stroke to remove them.  The single stroke buys ONE feed
at the cost of the bridges + ~10× worse field; the multi-wire design buys
the best field + no bridges at the cost of N feeds.

## Compensating the single-stroke degradation

If you keep the single stroke but want to claw back the field quality the
bridges cost you, there are two levers:

### Path-A (re-contour, pure single stroke) — capped ~8 %

`--compensated-iter` folds the chain's parasitic field back into the SF
target and re-solves.  It is the only PURE single-stroke compensation (no
extra feeds), but it is capped: 9.3 % → 8.1 % (step-sensitive, oscillates,
best-psi-tracked).  A **pure** single stroke has uniform current and no
spare degrees of freedom, so it fundamentally cannot fully cancel the
fixed bridge field — and the Path-A trick of re-contouring is limited by
the chain being non-smooth in ψ.  (`--freeze-levels` to "smooth" it is a
NEGATIVE result — it removes the very level-drift that lets the best-psi
search explore, and finds nothing.)

### Shim loops via LS-OMP — the iterative compensation methodology

`--shim-loops K` (or `--shim-tol RMS`) keeps the single-stroke main wire and
adds **independent shim correction loops** that cancel the residual
`r = B_target − I_w·Bz_chain`.  The exact linear compensation is
`δI = A⁺ r` — a full independent VARYING-current distribution (= multi-wire,
0 % residual).  Shim loops realise its dominant part.  **The way to iterate
is Order-Recursive / Least-Squares OMP (`--shim-method ls_omp`, default)**:
at each step orthogonalise every candidate basis loop against the current
support and add the one whose ORTHOGONAL component (its actual least-squares
residual reduction, normalised by its norm) is largest; re-solve the LS over
the whole support; repeat.  This corrects the column-norm bias of plain OMP
(the SF basis loops have very different field magnitudes) and is the optimal
forward-greedy step.  The residual decreases **MONOTONICALLY** — guaranteed
convergence, no oscillation, no step tuning:

| feeds (1 + K) | LS-OMP RMS | plain OMP | top-K |
|---------------|------------|-----------|-------|
| 1 (no shim)   | 9.3 %      | 9.3 %     | 9.3 % |
| +3            | **5.5 %**  | 6.1 %     | 7.8 % |
| +5            | **3.9 %**  | 4.6 %     | 6.7 % |
| +10           | **2.3 %**  | 3.6 %     | 5.2 % |
| +20           | **1.1 %**  | 1.9 %     | 4.2 % |

LS-OMP reaches ~40 % lower RMS than plain OMP and ~4× lower than top-K for
the same feed count.  `--shim-tol` turns it into a spec-driven design: e.g.
`--shim-tol 0.015` adds loops until DSV RMS ≤ 1.5 % and reports the feed
count needed (16 here).  The shim currents are tiny (< 1 % of I_w at these K).

**Method survey (2026-05-31, web-search-informed)**: LS-OMP/ORMP was
benchmarked against the standard sparse-recovery upgrades of OMP —
Subspace Pursuit, CoSaMP, and convex L1/LASSO (the usual MRI-shim
sparse-coil tools).  For THIS problem — best K-term APPROXIMATION of a
dense residual with very few constraints (M = 25 ≪ N = 960) — the
forward-greedy LS-OMP wins; SP/CoSaMP/LASSO are tuned for exact sparse
recovery and their prune-back / shrinkage steps discard useful atoms here.
The decisive ingredient is the column-norm normalisation (orthogonalised
greedy), not a fancier combinatorial search.

The honest trade-off stands: **compensating the single-stroke degradation
costs independent feeds** — there is no free lunch within one
uniform-current wire — but OMP makes the iteration MONOTONE and gives the
minimum feeds for any target accuracy.  (A second equal-current SF coil —
"patch" the residual with another stream function and single-stroke it —
does NOT work: the residual is a high-frequency field whose SF correction
is a non-smooth stream function that does not contour into a clean
equal-current coil; it stalls at ~6–8 %.  The correction MUST be realised
as independent varying currents, which is exactly what OMP does.)

### Single-current "sheet-metal" coil distortion (bankin-ho) — no extra feeds

The shim trade-off above has one hidden assumption: the wire stays a **flat
planar pattern**.  Within that constraint cancelling the residual costs
independent feeds.  But a manufactured coil need not be flat.  `--distort`
keeps **one series current** and the stream-function contour LEVELS fixed,
and instead BENDS the single-stroke wire in 3D — a smooth low-dimensional
deformation field

```
d(x, y) -> (delta_x, delta_y, delta_z)     # an n_grid x n_grid bilinear
                                           # control field per active axis
```

the discrete realisation of an NGSolve VectorH1 mesh deformation.  This is
"bankin-ho" / sheet-metal forming — a geometric shape optimisation that frees
the wire to leave the plane.  **Trading geometric DOFs for current DOFs, a
single current cancels most of the single-stroke residual** (planar uniform
Bz benchmark, 33 contours, dense 21×21 eval grid; integrated `--distort` run,
flat single-stroke baseline 12290 ppm):

| `--distort-penalty` | eval MAE (ONE current) | max bend | feeds |
|---------------------|------------------------|----------|-------|
| (flat single stroke) | 12290 ppm | 0 mm | 1 |
| 1.0 (conservative)  | 2015 ppm  | 17 mm | 1 |
| 0.3 (balanced, default) | 1025 ppm | 24 mm | 1 |
| 0.1                 | 605 ppm   | 31 mm | 1 |
| 0.03 (aggressive)   | 340 ppm   | 34 mm | 1 |

Convergence is fast + monotone — penalty 0.1: 12290 → 3656 → 1060 → 687 →
613 ppm in 6 iters.  Gauss-Newton over the control-grid displacement field;
the step solves `(JᵀJ + (λ + λ_disp) I) δ = Jᵀr − λ_disp·d` with the
displacement penalty `λ_disp` RELATIVE to `mean(diag(JᵀJ))`.  Larger penalty
= smaller (more 3D-printable) bends; smaller penalty = lower MAE.  The series
current is re-fit (optimal single value) every step, and the honest
dense-grid MAE is tracked for the returned "best".

`--distort-comps {xyz, xy, z}` selects the deformed axes: full 3D
sheet-metal (`xyz`, best), pure out-of-plane bend (`z`), or in-plane
reroute (`xy`).  **In-plane (`xy`) alone is the weakest** — the connector
current it reroutes still sits in the DSV plane, so its parasitic Bz is only
moved around, not cancelled.  The out-of-plane lift is what lets a single
current actually cancel the parasitic field (a connector lifted toward/away
from the DSV changes its Bz weight).  An earlier contour-flow PoC that moved
each wire normal to itself IN-PLANE only managed a ~40 % reduction in 40
iters; full-3D Gauss-Newton gives an >80 % reduction in 5 (12290 → 2015).

**This is the single-current answer to the shim trade-off.**  The "no free
lunch within one uniform-current wire" caveat assumed a fixed planar
geometry; give the wire 3D shape freedom and one current reaches the
~340–2000 ppm class — manufacturable as a single 3D-printed bent conductor,
no extra power supplies.  Separate-feed shims (`--shim-loops`) still reach a
lower absolute floor (183 ppm at 10 feeds) but require independent supplies;
distortion trades a little accuracy for one feed + one printable part.  The
two compose: distort first (one part), then add a couple of feeds only if
the spec demands sub-500 ppm.

```bash
# single-current sheet-metal coil, ~600 ppm at ~31mm bends, 1 feed
python demo_planar_uniform_fem_psi.py --order 3 --nlevels 16 \
    --n-sample 121 --eval-n 21 --distort --distort-comps xyz \
    --distort-penalty 0.1
```

#### Cylinder port — `demo_sf_to_peec_gx.py --distort`

The same mechanism ports to a **cylinder** (the Gx fingerprint gradient coil,
HARD tier).  The cylinder analog of the planar out-of-plane z-lift is a
**radial bend** `δr` (out of the cylinder surface); azimuthal `s = a·δφ` and
axial `δz` reroute are also available.  `coil_distort_cyl` runs the same
Gauss-Newton with a displacement-Tikhonov penalty on a φ-periodic `(φ, z)`
control grid.

| stage | DSV RMS | note |
|-------|---------|------|
| continuous-SF ideal | ~0.4 % | design floor (24×40 grid) |
| → single-stroke (`field_aware`) | 8.5–9.3 % | hard-tier degradation |
| → **`--distort` (full r+s+z)** | **1.4 %** | ONE current, ~30 mm bend |

```bash
python demo_sf_to_peec_gx.py --nphi 24 --nz 40 --nlevels 12 \
    --distort --distort-comps rsz --distort-penalty 0.1
```

**Geometry-dependent lever (interesting reversal).**  On the *plane*, the
out-of-plane lift dominates and in-plane reroute is weak.  On the *cylinder*
it is the **opposite**: the in-surface reroute (`--distort-comps sz`,
azimuthal+axial) is the dominant lever and the radial bend alone
(`--distort-comps r`) is the weakest (and needs the largest displacement,
~50 mm).  Physical reason: cylinder wires already wrap in 3D, so repositioning
*along* the surface reshapes the current pattern more effectively than lifting
*off* it; for an internal DSV a radial move only weakly changes the wire-DSV
distance.  Full `r+s+z` combines both and wins.  The Gauss-Newton framework
picks the effective direction automatically — **the optimal sheet-metal
direction is geometry-dependent**.

> **The same reversal appears in surface-forming sheet-metal.**  This section
> bends the manufactured *wire* (ψ + contour levels FIXED).  A *distinct*
> sheet-metal lever FORMS the conductor *surface* and RE-SOLVES ψ to push the
> **(homogeneity, peak current density)** Pareto front — and it shows the same
> plane↔cylinder reversal (planar = out-of-surface bending; cylinder =
> in-surface), plus a finer rule: within the cylinder, AXIAL vs AZIMUTHAL
> in-surface forming is selected by the target's azimuthal order `m`
> (`Gx` m=1 → axial; `C2`=x²−y² ellipse m=2 → azimuthal helps too).  See
> [regularization.md § Pushing the front](regularization.md#pushing-the-homogeneity-peak-j-pareto-front)
> and `demo_pareto_deform.py` / `demo_pareto_cylinder_deform.py`.

**Compose with electric shims (`--distort --shim-loops K`).**  The geometric
bend and the separate-feed shims cancel *different* parts of the residual, so
they compose.  Critically, **the bend is far more feed-efficient**: one bent
wire (1 feed) reaches 1.4 % where ten electric shims alone (10 feeds) only
reach 2.3 %.  Adding shims on the *distorted* residual then refines it:

| coil | feeds | DSV RMS |
|------|-------|---------|
| single-stroke | 1 | 8.5 % |
| + 10 electric shims (no bend) | 11 | 2.3 % |
| **+ sheet-metal distort** | **1** | **1.4 %** |
| + distort **and** 10 shims | 11 | **1.0 %** |

So the practical recipe on the hard tier: **distort first** (one printable
part does the heavy lifting), then add a few electric shims only to clean the
remaining high-spatial-frequency residual.  Under `--distort`, the `[8]`
shim-only block is skipped and the shims are applied to the distorted residual
in `[9]` instead.

**ψ regularisation — `--regularize {tsvd, tikhonov, h1}`.**  The cylinder ψ
solve also takes a Tikhonov family: `tsvd` (ACA+TSVD mode-truncation,
default), `tikhonov` (ridge `(AᵀA + αI)ψ = AᵀB`, dense, `--alpha`), `h1`
(min-seminorm smoothest ψ via a graph Laplacian).  An α L-curve sweep shows
the relationship is **non-monotonic** (the contour topology jumps with α):
for the Gx fingerprint single-stroke, **TSVD mode-truncation is the best
regulariser** (8.45 %); the best ridge α≈1e-2·mean only ties it (8.97 %) and
H1 is worse (it over-spreads the current → longer connectors).  This is the
*opposite* ranking from the planar uniform case (where H1 was cleanest) —
another geometry-dependence.  Practical recipe: pick the regularisation that
minimises the **single-stroke RMS** (not the smoothest ψ), then distort.

### Min-inductance design ⊕ distort: same field, HALF the bend (2026-06-03)

The physical **min-inductance** objective (`calc_streamfunction.py --regularize
inductance`, `min ½ψᵀLψ` via the BEM self-inductance — see
[panel.md](panel.md#design-objective--regularizer---regularize-l2-h1-inductance))
gives a *smoother* surface current than `h1` / `l2`, so it single-strokes
cleaner AND needs a smaller sheet-metal deformation for the same delivered
field.  Measured on the cylinder Gx (`--confine abe`, nlevels 12, `--distort`):

| design | single-stroke (pre-distort) | post-distort | max bend |
|--------|-----------------------------|--------------|----------|
| `h1` | 4.85 % | 2.27 % | 9.8 mm |
| `inductance` | **2.61 %** | 2.26 % | **4.8 mm** |

Same final homogeneity, **half the deformation** — the min-inductance coil is
the more manufacturable 3D-printed / sheet-metal part.  So the maturity recipe
on the HARD tier is: **design min-inductance → abe (close the contours) →
field-aware single stroke → distort**.

**EASY tier is production-grade (sub-500 ppm, no shims).**  A uniform-Bz
(solenoid) target on the same cylinder has nested CLOSED rings (no irreducible
bridges), so the field-aware single stroke is a clean spiral: **254 ppm with no
distort, → 84 ppm after a 1.5 mm distort** (`--target-cf 1 --confine off
--nlevels 16 --distort`).  Locked by
`validation_test/panels/test_streamfunction_golden.py::
test_streamfunction_easy_tier_single_stroke_sub500ppm`.  The repository-backed
boundary is therefore explicit: **easy / nested targets reach the ~100 ppm
class single-stroke; the hard multi-lobe Gx is intrinsically ~2 %**, where
min-inductance + distort buys manufacturability (half the bend) rather than a
lower floor.

### Arbitrary curved formers (sphere) — FE-direct ψ ([`demo_sphere_fe_direct.py`](demo_sphere_fe_direct.py))

The basis-loop representation needs a structured `(φ, z)` grid, so it is
stuck on planes and cylinders.  **FE-direct ψ meshes ANY surface** — the
real payoff of the high-order (cubic) FE stream function.  On a *sphere*
former (NMR shim coil around the magnet), surface FEM done robustly (mesh
the solid, `H1(mesh, order=3, definedon=mesh.Boundaries(".*"))`, `ds` +
`grad(.).Trace()`, `mesh.Curve(3)` isoparametric) with the general kernel
`K = n̂ × ∇_s ψ` (the curved-surface generalisation of the planar `ẑ × ∇ψ`):

| target | continuous ψ | single-stroke | + sheet-metal distort |
|--------|--------------|---------------|-----------------------|
| uniform Bz (l=1) | cres 3e-15 (0 ppm) | **0.24 %** | — |
| Z2 shim (l=2,m=0) | cres 5e-14 | 4.3 % | **0.36 %** (1 current, ~2 mm bend) |

The full manufacturable pipeline runs on the curved former: FE-direct ψ →
single-stroke (latitude-ring spiral, real inter-turn spacing **10.5 mm** ≥
conductor width) → sphere sheet-metal distortion (radial `δr` + tangential)
→ 0.36 %.  Why it matters: the *design-surface shape is not a
manufacturability constraint* — any fabricable former (sphere / conformal /
3D-printed) holds a coil; the real constraints are the **wire pattern**
(single-stroke, min-spacing ≥ width, no self-cross), all handled here.  The
easy/hard split is set by **target complexity** (uniform l=1 → 0.24 % clean;
fingerprint → hard), NOT the surface shape.  Caveat: the supplied chainer is
for *axisymmetric* (m=0: Z1/Z2/Z3…) shims; m≠0 tesserals need a general
field-aware sphere chainer.

## Can single-stroke reach a 100 ppm-class spec? (No, for Gx)

For a high-precision target (~100 ppm = 0.01 % degradation), single-stroking
the **Gx multi-lobe** coil is fundamentally infeasible.  There is a
**resolution paradox** (measured on a dense 925-point eval grid):

| SF resolution | continuous-SF floor (ideal) | single-stroke | degradation |
|---------------|-----------------------------|---------------|-------------|
| 24×40, 12 lvl | 3968 ppm | 86 619 ppm  | 82 651 ppm |
| 36×60, 18 lvl | 996 ppm  | 816 406 ppm | 815 410 ppm |
| 48×80, 24 lvl | **528 ppm** | 730 873 ppm | 730 344 ppm |

- The **continuous-SF / multi-wire ideal improves with resolution**
  (3968 → 528 ppm; finer → approaches 100 ppm) — it is the only path to a
  100 ppm-class coil.
- The **single stroke gets dramatically WORSE with resolution**
  (8.7 % → 73-82 %), because a finer design has many more contours
  (68 → 893) and therefore many more bridges, which come to dominate the
  field.  The resolution you need for a 100 ppm IDEAL is exactly the
  resolution at which the single stroke is catastrophic.

Shim compensation (LS-OMP) converges toward the ideal but also overfits the
fit grid (eval RMS 2-3× the fit RMS), so closing an 80 % single-stroke
degradation to 100 ppm on a dense grid would need ~as many feeds as a full
multi-wire coil — i.e. no longer a single stroke.

**Conclusion**: a 100 ppm-class Gx coil MUST be the independent-closed-loop
(multi-wire) design; single-stroke is intrinsically a ~1-80 % degradation
technique for multi-lobe topologies.  100 ppm single-stroke IS feasible for
**EASY-tier** topologies — cylindrical Gz (a smooth helix is naturally one
wire, no bridges) and planar nested families — where there are no
irreducible inter-lobe bridges.  Match the realisation to the topology:
nested/helix → single stroke; multi-lobe + high precision → multi-wire.

### Demonstrated: EASY-tier (planar uniform Bz) to the ~200 ppm class

The established iterative pipeline (in
[`demo_planar_uniform_fem_psi.py`](demo_planar_uniform_fem_psi.py))
reaches the few-hundred-ppm class with only a handful of feeds:

```bash
python demo_planar_uniform_fem_psi.py --regularize h1 --order 3 \
    --nlevels 30 --n-sample 141 \
    --compensated-iter 60 --compensated-step 0.05 \
    --shim-grid 9 --shim-tol-ppm 200
```

| stage | dense-grid MAE |
|-------|----------------|
| continuous-SF ideal (H¹ order 3) | ~95 ppm |
| single stroke + Path-A | ~1375 ppm |
| + 9 LS-OMP shim loops (10 feeds) | **183 ppm** |

Why it works here (and not for Gx): uniform Bz has NESTED single-sign
contours → the single stroke is a clean spiral with NO irreducible
bridges, only short radial connectors; a finer spiral pitch
(`--nlevels 30`) shrinks them, Path-A (monotone on FE-direct ψ) folds the
rest into ψ, and a few overlapping shim loops (LS-OMP, monotone) close the
residual to the ~200 ppm class.  Note the honest evaluation: the
`--eval-n` dense grid is used for the ppm figure — the n_target FIT grid
over-reports because the massively-underdetermined SF solve over-fits it.
The shim currents are a few × I_w (the small plane loops couple weakly to
the offset target plane); a Helmholtz-pair-style shim geometry would lower
them.  The same pipeline is the template for the cylinder via
ngsolve.bem high-order (future work).

## Complexity tier framework

A coil's REACHABLE design quality is bounded by its TOPOLOGY:

| Tier   | Topology                  | Baseline    | + Path-A     | Path-A effectiveness         |
| ------ | ------------------------- | ----------- | ------------ | ---------------------------- |
| EASY   | axisym. / planar uniform  | RMS 2–3 %   | < 1 %        | useful (basis-loop)          |
|        |   + FE-direct H¹ ψ        | 2 %         | **0.47 %**   | **MONOTONE convergence**     |
| MEDIUM | cylindrical Gz            | already OK  | redundant    | smooth helix natural         |
| HARD   | cylindrical Gx fingerprint| 16 % (kuijpers) → **9.3 % (field_aware)** | **8.1 % (field_aware + Path-A)** | chain method matters; Path-A composes |
| HARDER | shielded / biplanar / 3D  | —           | —            | needs FE-direct or D-path    |

Past EASY needs FE-direct continuous ψ.  **The HARD tier "16 % ceiling"
was a `kuijpers`-method artifact, not a fundamental bound**: the
`field_aware` chain (sign-order + azimuthal-min cuts, 2026-05-31) reaches
9.3 % on the same SF design with no Path-A — so the connection METHOD, not
just the SF design, gates HARD-tier quality.  The two best ideas COMPOSE:
`field_aware` + Path-A (`--compensated-iter 40 --compensated-step 0.3`)
reaches **8.1 %**, roughly half the old `kuijpers` baseline.  Going lower
still likely needs B-spline SFD (Kuijpers Methods 2/3) or multivalued-
potential reformulation (Path D in `streamfunction(topic=single_stroke)`).

## Dead-end variants (do NOT re-try)

Six attempts on 2026-05-30 that all *reduced* field accuracy:

**Chain-orientation traps** (memory entry
`feedback_single_stroke_chain_orientation_traps`):

  - Sort by sign(centroid_x): RMS 16 → **43 %** (lobe reversal flips
    current direction for half the loops).
  - Force top/bottom cut on closed contours: RMS 16 → **30 %** (same
    flip).
  - Pair-aware reverse (`pos_y` outer-to-inner then `neg_y` reversed
    inner-to-outer): RMS 16 → **40 %**.
  - 8-sub-lobe (sx, sy, sz): RMS 16 → 35 %.

Common cause: any per-contour traversal-direction reversal flips
matplotlib's natural CCW orientation → current direction inverts →
field broken.

  **SAFE knob** = global within-lobe SORT KEY.
  **UNSAFE** = per-contour TRAVERSAL DIRECTION change.

**Path-A on basis-loop ψ** (memory entry
`feedback_path_a_naive_picard_negative`):

  - α=1.0 full step: 16 → 55 % (divergence).
  - α=0.5 / 0.1: oscillates 40–55 % (no contraction).
  - α=0.05 + N=60 iters: oscillates, best 15 % by accident.

Common cause: matplotlib contour topology jumps under small ψ
perturbations → `B_c(ψ)` not smooth → Picard does not contract.

## Why FE-direct ψ unsticks Path-A

`H¹` GridFunction is a CONTINUOUS function of its DOF vector → contour
family deforms SMOOTHLY → matplotlib emits a consistent contour count
and ordering → `B_c(ψ)` smooth → naive Picard CONTRACTS.

Empirical proof (planar uniform Bz):

  iter 40: residual 0.62 %
  iter 41: 0.60 %
  iter 42: 0.58 %
  ...
  iter 47: **0.49 %**       ← monotone, no backtracking

Same iteration on basis-loop ψ oscillates 1.7 % → 15 % around the best.

## Visualising the chain

The visualization modes (`contours`, `chain`, `step`) are tracked by the
`view_sf_coil_gx_gmsh.py` record in
[`theory.ipynb`](theory.ipynb). Keep new public
visualization entry points in docs notebooks or validation runners rather than
adding a fresh examples command path here.

`--mode contours`: shows the SF design's raw closed contours (no
connection arcs).  Each contour is a separate Physical Group in GMSH
so you can toggle them individually.

`--mode chain`: shows the single-stroke chain (= what PEEC sees).  Use
to inspect the rung pattern and identify any criss-cross artefacts.

`--mode step`: opens the multi-piece loft chain STEP from
`demo_sf_to_peec_gx.py --with-peec`.

(Pitfall: GMSH may restore the window to an off-screen monitor; the
viewer scripts force `General.GraphicsPositionX/Y` to (100, 100).  If
you still see no window, see `gmsh_usage(topic=pitfalls)` #9.)

## Cross-reference

  - Math: [theory.md](theory.md) Path-A section
  - Open extension paths: the paper outline (W:\02_学会資料\2025年度\2026_01_JIAM\streamfunction\) Path A/B/D
  - MCP topic: `streamfunction(topic=single_stroke)`,
    `streamfunction(topic=session_2026_05_30)` sections 2 (chain methods) +
    3 (Path-A representation) + 4 (dead-end variants).
