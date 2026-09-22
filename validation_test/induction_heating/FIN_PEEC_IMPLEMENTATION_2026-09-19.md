# Fin-capable PEEC: first conservative topology slice (2026-09-19)

This is an internal implementation/validation ledger, not a user workflow.
The production IH PEEC path remains the existing series-filament bundle.
No CLN code or CLN-derived snapshots are used by this work.

## Implemented

- `radia.peec_fin_topology` accepts CAD-derived perimeter station rings.
- Longitudinal branches exist between every pair of stations. Transverse
  branches exist only at selected interior stations. Every lane retains its
  own transition node; no artificial equipotential junction is introduced.
- The STEP adapter accepts only direct perimeter extraction
  (`step_uv`, `step_per_station`, `step_section_planes`). Equivalent-circle
  fallbacks are rejected. For an unknown section, the first-pass auto rule
  meshes every interior station; selective mesh removal is not yet claimed.
- A small dense MNA reference solve and a C++ PEEC assembly experiment both
  preserve node KCL on the hybrid graph. The terminal is a logical common
  node while each branch retains its physical CAD endpoints for L assembly.

## Not yet accepted as physical IH fin solver

The C++ experiment approximates a surface branch as a thin rectangular
segment. Its partial self/mutual inductances, internal-reactance split, and
surface resistance at a sharp fin have **not** been validated. There is no
SIBC coupling, no frequency/temperature-dependent skin-depth selection, no
thin-fin two-face coupling, and no branch-current back-reaction to the
workpiece BEM. It must not be enabled in `calc_inductance.py` or the Simulink
IH operator assembler on the strength of topology/KCL tests alone.

## First independent BEM-A comparison (not accepted)

`compare_beak_fin_bema.py --experimental-peec --n-peri 64` uses the same
6 mm synthetic STEP for both routes, copper conductivity 58 MS/m, 150 kHz,
and a 0.171 mm SIBC skin depth. BEM-A uses an independently generated
1156-triangle NGSolve surface mesh (1734 HDivSurface current DoFs); its
current-continuity residual is about 3.2e-16. The experimental PEEC graph
has 448 branches. The computed terminal/energy quantities were:

| Route | R (micro-ohm) | External L (nH) |
| --- | ---: | ---: |
| BEM-A impedance-EFIE | 34.372 | 1.12291 |
| Experimental rectangular-branch PEEC | 29.186 | 1.05442 |

PEEC is 15.1% low in R and 6.1% low in external L. These are **not**
equal and not a validation pass. The PEEC internal sheet reactance was
included in its solve but excluded from the reported external L, matching
the BEM-A result convention. Perimeter refinement from 16 to 64 lanes did
not close the R discrepancy. BEM-A mesh/PEEC axial convergence, spatial
current and tip/root loss maps, and the physical fin SIBC matrix remain
required. The original 60 mm fixture yielded ~9662 BEM-A surface triangles
even at 6 mm target maxh due to the 0.25 mm nose; a dense direct solve was
stopped after reaching ~6.8 GB working set. The 6 mm short fixture is a
bounded exploratory comparison, not a substitute for the full-length gate.

The current-distribution gate is stricter than terminal R/L. On the tracked
20 mm fixture at 150 kHz (BEM-A: 3520 triangles / 5280 current DoFs; PEEC:
64 perimeter lanes / 9 stations), terminal R and external L differ by only
+1.23% and -1.86%, respectively. Nevertheless, after mapping the central
30% of the BEM-A surface to the same perimeter lanes and removing only the
arbitrary global phasor, PEEC has 13.1% complex relative L2 error. It predicts
the mean beak |K| 23.6% high, the rounded-tip |K| 8.7% high, beak loss 33.6%
high, and tip loss 20.7% high. Therefore terminal agreement does not clear the
model. `beak_fin_bema_comparison_150kHz.json` records `accepted=false`; the
per-lane evidence is `beak_fin_current_profile_150kHz.csv`. Acceptance requires
profile L2 <= 5%, beak/tip mean |K| within 5%, and mean loss within 10%, after
both discretizations have independently converged.

The attached 2026-09-19 internal review also identifies the production-main
bulk impedance allocation and proximity iteration as non-SIBC-consistent.
Those paths must not be reused for this fin backend. The experimental branch
model uses a thin surface-band resistance and transverse MNA links, but its
rectangular Ruehli self term is not yet a validated surface partial element;
the failed local-current gate above confirms that further formulation work is
needed even where terminal quantities happen to agree.

### Review correction: cap convention is the dominant 6 mm discrepancy

The preceding rejection remains valid as a conservative gate, but its causal
interpretation is superseded here. The experimental `peec_fin_topology` path
does **not** use main's bulk `Zs_fil` or proximity iteration. It already uses
`R=Rs*length/dual_width`, `Zs=jR`, and transverse MNA branches at interior
stations. Review findings A--D apply to the existing production path, not this
experimental assembly.

Direct partition of BEM-A's sampled surface-current loss establishes the
terminal convention mismatch. At 6 mm, BEM-A total R is 34.372 micro-ohm:
6.531 micro-ohm is on the source/sink caps and 27.841 micro-ohm is lateral.
The experimental PEEC has ideal lossless equipotential terminal nodes and
gives 29.186 micro-ohm. Thus its difference from BEM-A total is -15.1%, but
from BEM-A lateral-only loss is +4.83%. The cap term alone is larger than the
original 5.19 micro-ohm discrepancy.

The 6/12/24 mm `cocr` sweep in `beak_fin_length_sweep_150kHz.json` reinforces
this: BEM-A cap R stays 6.53/6.76/6.90 micro-ohm, while total relative PEEC
error changes -15.1% -> -4.11% -> +2.24%. The absolute cap contribution is
approximately length-independent and its relative importance decays with
length. PEEC remains 4.8--8.3% above BEM-A lateral-only R, so discretization
and local-current convergence are not yet certified.

The local |K| comparison is also not yet a truth test at the nose: R_tip/delta
is only 1.47, violating a robust thin-skin asymptotic regime; 64 equal-arc
lanes place only about two lanes on the rounded tip; and the short-fixture
sampling window lies within end-effect penetration. BEM-A and thin-sheet PEEC
therefore require geometry/frequency and mesh convergence before applying the
numerical acceptance thresholds above. Neither current result validates the
physical nose current.

### Independent lateral-discretization convergence

`sibc2d_beak_reference.py` supplies the third, infinite-length reference
instead of treating either 3-D route as truth. Its resistance converges from
5.425 to 5.248 mOhm/m over 64--1024 perimeter panels, consistent with the
5.242 mOhm/m continuum value; rounded-tip mean |K|/mean(K) converges near 2.0.
A 256-panel regression requires R within 1% and the tip ratio within 5%.

Using the 12-to-24 mm difference quotient to cancel terminal effects:

| Route/refinement | lateral R (mOhm/m) | error vs 5.242 |
| --- | ---: | ---: |
| BEM-A maxh 3.0 mm | 4.977 | -5.05% |
| BEM-A maxh 1.5 mm | 5.027 | -4.10% |
| BEM-A maxh 0.75 mm | 5.105 | -2.61% |
| PEEC n=64, axial <=1.5 mm | 5.433 | +3.64% |
| PEEC n=128, axial <=1.5 mm | 5.335 | +1.78% |
| PEEC n=256, axial <=1.5 mm | 5.272 | +0.58% |

Thus BEM-A rises and PEEC falls toward the independent 2-D value. The earlier
4.8--8.3% lateral mismatch was shared discretization error, not evidence that
the experimental PEEC SIBC topology was physically wrong. This clears the
straight-prism **integrated lateral loss** formulation at the finest PEEC
point. It does not clear pointwise nose current or a real tapered/curved fin:
the R_tip/delta limitation, only about two 64-lane nose samples, and absence of
an independently converged 3-D tip profile remain explicit blockers. Full data
are in `beak_fin_discretization_convergence_150kHz.json`.

### Fin-delivery metrics replace pointwise K acceptance

The early 20 mm pointwise-profile rejection above is retained as provenance
but is no longer the acceptance definition. A rounded nose with R_tip/delta
near one and only a few lanes must not be judged by max(K), and terminal R/L
is insensitive to whether useful current reaches the fin. The primary metrics
are now beak current fraction, beak/tip integrated loss fraction, beak-current
centroid, and the H profile on x=5 mm (1 mm beyond the nose). Pointwise max(K)
is diagnostic only.

The independent 2-D n=1024 reference gives I_beak/I=0.30815,
P_beak/P=0.34167, P_tip/P=0.20231, beak-current centroid x=2.8531 mm,
tip mean |K|/mean(K)=2.014, and centre-probe |H|=41.784 A/m. For the finest
completed 24 mm comparison (BEM-A maxh=0.75 mm, PEEC n=256 and axial <=1.5
mm), BEM-A/PEEC respectively give current fraction 0.3040/0.3066, loss
fraction 0.3259/0.3365, centroid 2.8168/2.8309 mm, and probe H 41.06/41.45
A/m. Probe-line relative L2 error is 0.72%.

Acceptance requires current and loss fractions within 3%, centroid within one
perimeter-lane spacing, and probe-line H within 2%. Current fraction, centroid,
and H pass; beak loss fraction differs by 3.24%, narrowly failing the declared
limit. Therefore the strict combined result remains `accepted=false`, without
mischaracterising the well-converged integral and delivered-field quantities.
See `beak_fin_delivery_metrics_150kHz.json`.

### The pairwise 3% criterion is not decidable on this fixture (2026-09-21)

The 3.24% beak-loss-fraction difference was read as BEM-A's own discretisation
error, since against the independent 2-D n=1024 reference PEEC is 1.5% low and
BEM-A 4.6% low. Refining BEM-A does not close it. Over maxh 0.75, 0.55 and
0.45 mm its beak loss fraction is 0.325491, 0.323364 and 0.324678 -- -4.74%,
-5.36% and -4.97% against the 2-D value, non-monotone and flat near -5% -- and
the pairwise difference moves 3.00%, 3.68%, 3.26%, crossing the threshold in
both directions while neither route changes physically.

The reason is visible in the mesh: the surface face count grows only 5990 ->
7828 -> 9820 while the linear size falls by 1.67x, because the 0.25 mm tip
radius already drives the local refinement. Reducing maxh adds faces on the
flat body, not on the fin where this metric is measured.

Acceptance above already requires that both discretisations have independently
converged. That precondition is now measured to be unmet on the BEM-A side, so
the gate cannot be closed with the controls the comparison exposes. Closing it
needs either per-route scoring against the independent 2-D reference -- the
rule the mixed-Omega work adopted for the same reason -- or a fin-local
refinement control for the BEM-A surface mesh. Data:
`results/beak_fin_bema_refinement_20260921.json`.

### Geometry order is not the cause; basis order is (2026-09-21)

The BEM-A side of this comparison is geometrically first order.
`compare_beak_fin_bema.py` builds its surface through
`_extract_surface_mesh_filtered`, which returns a flat P1 surface mesh and
never calls `mesh.Curve`, so the 0.25 mm rounded tip is a polygon at every
`maxh`. `extract_surface_curved` (geom_order >= 2) exists but returns raw P2
node arrays for the SIBC-HACApK kernel and does not feed
`compute_inductance_source_sink`, which this comparison uses.

The measurement separates a geometric error from a discretisation one. Across
maxh 0.75, 0.55 and 0.45 mm the PEEC-over-BEM-A tip ratios are constant --
1.026, 1.026, 1.024 in mean |K| and 1.061, 1.061, 1.057 in mean loss -- while
the face count changes by 1.6x. A discretisation error shrinks under
refinement; this does not, because the tip's element size is already set by its
own curvature and `maxh` cannot reach it. BEM-A's tip loss fraction stays at
0.1846 +/- 0.0007 against the 2-D reference 0.2023, i.e. -8.7% and flat, where
PEEC is -4.0%.

This is not evidence against the SIBC model. The 2-D SIBC reference converges
to the continuum value and PEEC, which uses the same SIBC sheet, lands within
4% of it; the disagreement follows the tip geometry, not the impedance model.
The BEM-A solve is algebraically sound on the mesh it is given: R partition
closure -4.3e-15, current continuity about 3.2e-16.

The consequence for the gate is that BEM-A cannot serve as the converged
reference this comparison assumes until `compute_inductance_source_sink` is
given a curved surface, or this comparison is routed through the curved SIBC
path.

### Measured: basis order, not geometry order (2026-09-21)

The section above named geometry order. That code fact holds -- the coil
BEM-A is geometrically first order while the workpiece BEM-A reads
`vol_curve_order` and runs the curved SIBC kernel -- but measurement shows it
is not the cause, and this entry replaces that attribution.

`HDivSurface` also lives on the boundary of a curved volume mesh with the same
DOF count, so the same solve was run with the tip as an arc rather than a
polygon. Curving raises the tip boundary area from 3.14395e-05 to 3.15880e-05
m2 and converges by order 2, but R moves only from 125.81210 to 125.89625
micro-ohm: **+0.067%**. That cannot explain an 8.7% tip-loss deficit.

Basis order can. `compute_inductance_source_sink` defaults to `fes_order=0`,
RT0/RWG, the lowest order, and the comparison never varied it. One order step
on maxh 1.8 mm moves R from 123.68237 to 125.77627 micro-ohm, **+1.69%**,
about 25x the geometry-order effect, and `fes_order=1` at maxh 1.8 mm nearly
reproduces `fes_order=0` at maxh 0.75 mm. One order is therefore worth roughly
a 2.4x mesh refinement here, and R is still rising.

So the BEM-A side was never p-converged. The comparison only ever varied
`maxh`, which saturates because the tip element size is set by its own
curvature. Converging it needs `fes_order`; `hacapk_cocr` was RT0-only until
later the same day, when `_dof_cluster_coords` gave its cluster tree one point
per DOF for any order (and `compute_inductance_source_sink` gained `Compress`,
so the coil `.vol` is consumed whole and its curving order reaches the solve --
the same boundary the workpiece BEM already used).  The next entry records
what that order-1, curved run showed.

None of this indicts the SIBC model. The 2-D SIBC reference converges to the
continuum value and PEEC, on the same SIBC sheet, lands within 4% of it.

### Under every control the gate metric stands still (2026-09-21, later)

With `hacapk_cocr` now taking any order and the BEM-A boundary curved, the
gate was run at fes_order 1 and Curve(2) on maxh 0.75 (n_J 17970): R rose to
126.449 micro-ohm (+0.51%) but the beak loss fraction is 0.324837 and the tip
loss fraction 0.184909 -- unchanged from RT0 flat (0.325491 / 0.185637), so the
pairwise gap is 3.21% and the gate still reads `accepted=false`.  The tip is
not under-resolved in the observable either: 977 BEM-A triangles sit on the
tip inside the axial window, about 2.3 per PEEC lane per station.

So h (0.75 -> 0.45 mm), p (0 -> 1), and geometry (1 -> 2) have each been moved
and the delivery fractions did not follow any of them.  BEM-A is converged on
this fixture at about 3% below PEEC and about 5% (beak) / 8.6% (tip) below the
infinite-length 2-D reference, with PEEC itself 1.9% / 4.0% below it.  The 3%
criterion is therefore straddling a converged model difference, not a
discretisation error, and the two entries above that attributed it first to
geometry order and then to basis order are both superseded: geometry order was
measured too small, basis order moves R but not the gate metric.

Both routes sit below an infinite-length reference in the same order, so end
effects of the 24 mm fixture are the first candidate; the impedance-EFIE versus
thin-sheet PEEC current models are the second.  A length sweep of the
fractions discriminates the first; required gate 3 (an independent 3-D
reference) is what settles the second.  Data:
`results/beak_fin_bema_refinement_20260921.json`.

### The gate closes at 48 mm: it was the fixture length (2026-09-21)

Doubling the fixture converges both routes onto the independent 2-D reference
and the delivery gate passes.  The beak loss fraction goes from -4.93% (BEM-A)
and -1.88% (PEEC) of the reference at 24 mm to **-0.24%** and **+0.34%** at
48 mm; tip loss from -8.60% / -4.03% to -0.68% / +0.53%; current fraction to
+0.16% / -0.13%.  The two routes bracket the reference on all three fractions,
which is what two independently converged discretisations do around a true
value.  The pairwise difference falls from 3.21% to **0.58%**, and every
acceptance metric passes with at least five times margin: current fraction
0.28%, loss fraction 0.58% (limit 3%), centroid 0.0100 mm against a 0.0894 mm
lane spacing, probe-line H 0.17% (limit 2%).  Terminal agreement follows:
PEEC is +0.14% in R and -0.58% in L, against -3.01% and -1.06% at 24 mm.

The 48 mm run uses `fes_order=0`, so the basis order was not what the fractions
needed either -- only R responded to it.

Why length: the fractions are ratios measured mid-span, and the caps carry an
approximately length-independent loss (6.53, 6.76, 6.90 micro-ohm at 6, 12,
24 mm in the earlier sweep).  At 24 mm that end contribution still depresses
the mid-span beak and tip shares of both routes, by different amounts because
the two discretisations distribute the end current differently.  It is not a
discretisation error in either, which is why no refinement control reached it.

`beak_fin_delivery_metrics_150kHz.json` now records `accepted=true` on the
48 mm basis, with the 24 mm block retained as provenance.  Scope: a straight
constant-section fin in copper at 150 kHz with both discretisations
independently converged -- not tapered or curved fins, brazed fin solids,
filleted roots, nonlinear material, or the gate-4 back-reaction.  Data:
`results/beak_fin_length_resolution_20260921.json`.

### Gate 1 closes: the metrics are perimeter-limited and axially converged (2026-09-22)

The delivery gate had passed at exactly one discretisation, which says nothing
about convergence.  Refining each direction independently against a fixed
BEM-A reference settles it.  Perimeter, at 33 stations:

| `n_peri` | branches | R (uOhm) | beak current | beak loss | probe H | accepted |
|---|---|---|---|---|---|---|
| 32 | 2016 | 265.930 | 3.96% | 10.92% | 2.63% | no |
| 64 | 4032 | 257.828 | 0.55% | 3.29% | 0.77% | no |
| 128 | 8064 | 253.305 | 0.20% | 0.71% | 0.36% | **yes** |
| 256 | 16128 | 250.422 | 0.28% | 0.58% | 0.17% | **yes** |

R descends monotonically onto the BEM-A reference `250.074 uOhm` and is within
`0.14%` of it at `n_peri=256`.  Axial, at 256 lanes:

| `n_stations` | branches | R (uOhm) | beak current | beak loss | accepted |
|---|---|---|---|---|---|
| 5 | 1792 | 250.458 | 0.32% | 0.49% | **yes** |
| 9 | 3840 | 250.498 | 0.30% | 0.58% | **yes** |
| 17 | 7936 | 250.478 | 0.29% | 0.60% | **yes** |
| 33 | 16128 | 250.422 | 0.28% | 0.58% | **yes** |

The axial column is flat.  Five stations already give the converged answer,
and the geometry rejections -- station alignment, zero-area cells, a present
tip -- hold at every level.  Data: `results/beak_fin_refinement_gate_20260922.json`.

The accepted 48 mm basis is now a tracked fixture,
`tests/coil_from_cad/fixtures/beak_fin_48mm.step`, checked against the
generator by solid volume, area and bounding box rather than by STEP bytes.
Before this it existed only in a scratch directory, so the closed gate could
not be reproduced from the repository.

### The fin PEEC is 319x faster at the converged discretisation (2026-09-22)

The axial column above is not only free of information, it is where all the
cost is.  Measured on the 48 mm fixture at 256 lanes:

| `n_stations` | branches | seconds | resident | R (uOhm) |
|---|---|---|---|---|
| 33 | 16128 | 290.03 | 6.81 GB | 250.4218 |
| 17 | 7936 | 37.01 | 1.56 GB | 250.4783 |
| 9 | 3840 | 5.66 | 0.36 GB | 250.4977 |
| 5 | 1792 | **0.91** | **0.08 GB** | 250.4580 |

Thirty-three stations cost `290 s` and `6.8 GB` to move R by `0.014%` against
five.  The converged configuration is `n_peri=256, n_stations=5`, and it runs
in under a second.

A compressed PEEC was tried first and rejected.  `PEECCircuitSolver` already
carries a HACApK path, and routing the fin assembly through it avoids the
dense L fill; but its nodal saddle solve does not converge on this graph --
`8.4e-2` residual after 2000 outer matvecs -- and tightening the H-matrix
accuracy to match the inner tolerance did not rescue it.  The option is not
kept: a code path that does not converge is not an option, and the
discretisation result above removes the need for one.

Compression does not help the BEM-A reference either, which is now the whole
remaining cost.  On the same 48 mm system, `n_J = 16962`:

| solver | seconds | R (uOhm) | residual |
|---|---|---|---|
| `cocr` | 437.3 | 250.074338 | 1.945e-10 |
| `hacapk_cocr` | 459.8 | 250.074524 | 1.951e-10 |

The H-matrix path is `5%` slower for the same answer -- R agrees to `7.5e-7`,
L to `7.5e-9`.  At this size the dense COCR is already the right tool; the
compression overhead is not repaid.  That also corrects the cost attribution:
the original `1164 s` was `437 s` of BEM-A plus `290 s` of dense PEEC plus
meshing, so with the PEEC at `0.91 s` the comparison is about eight minutes
and is entirely its reference.

The reference mesh is not a lever either:

| `maxh` | faces | `n_J` | BEM-A (s) | R (uOhm) | vs finest |
|---|---|---|---|---|---|
| 2.00 mm | 8280 | 12420 | 217.5 | 244.660 | -2.17% |
| 1.50 mm | 8528 | 12792 | 231.9 | 246.629 | -1.38% |
| 1.00 mm | 9506 | 14259 | 291.1 | 248.355 | -0.69% |
| 0.75 mm | 11308 | 16962 | 474.1 | 250.074 | 0 |

A `2.7x` coarser `maxh` removes only `27%` of the faces, because the floor is
set by the `0.25 mm` tip radius rather than by `maxh`, and it buys `2.2x` at
the price of a reference that is itself `2.17%` from converged.  The gate
closes at every level, but the beak-loss error wanders -- `0.40%`, `0.22%`,
`1.59%`, `0.49%` -- which is the reference moving under it, not the PEEC.
A gate decided against an unconverged reference is not worth the minutes.
`0.75 mm` stays.  Data:
`results/beak_fin_reference_mesh_cost_20260922.json`.

So the speed answer is asymmetric and worth stating plainly: the production
fin path is now under a second, and the validation comparison is about eight
minutes of irreducible reference.  Both available levers on the reference --
compression and mesh -- were measured and rejected.

### Gate 2 measured: the branch mutual is filamentary (2026-09-22)

`MutualInductanceRectBar` averages the Neumann kernel over both cross-sections
and exists, in the kernel's own comment, because filamentary Neumann gives a
"spurious circulating current artifact ... for close parallel bars".  It is
reached only when both segments are sub-filaments of one parent, which
`add_connected_segment` never produces.  Every fin branch pair therefore takes
the filamentary formula, confirmed behaviourally: the built L matrix
reproduces Grover's equal, aligned, parallel filament closed form to machine
precision at four separations.

That is exactly the case the comment warns about.  At the accepted
discretisation the perimeter spacing is `0.0894 mm`, the band width is the
same `0.0894 mm`, and the sheet depth is the `0.1706 mm` skin depth: adjacent
bands touch, and their separation is smaller than their own depth.  Against
the cross-section-averaged value over the same two bars, **the filamentary
mutual is 5.71% high on average and 5.73% at worst** -- mean and maximum
together, so it is a uniform bias rather than a local artifact.

`5.7%` is larger than the 3% the delivery gate allows, and the gate passes at
`0.28%` and `0.58%` anyway.  Both are true: a uniform bias on every
neighbour mutual shifts L without redistributing the current much.  The
element is not thereby validated -- it gives a passing answer for a reason
unrelated to its accuracy.  Data:
`results/fin_partial_element_mutual_20260922.json`; the behaviour is pinned by
`tests/test_fin_peec_branch_mutual.py`.

A first reading of this attributed the error to the cross-section orientation,
since `PEECSegment` carries no frame and `MutualInductanceRectBar` builds one
from a global fallback axis.  That was wrong: the averaged kernel is not
reached at all, so its frame never applies.  The orientation sensitivity is
retained in the artifact as a secondary number -- `7.35%` -- because it
becomes the next question if the averaged path is ever adopted.

### Gate 3 opens: the two routes agree with each other and both miss (2026-09-22)

BEM-A and the surface PEEC agree on the delivery metrics to `0.58%`.  They also
both impose a Leontovich surface impedance, so that agreement tests their
discretisations and says nothing about the assumption they share.  An
interior-resolved A-V solve now sits beside them.

The reference is `radia.eddy_aphi`: a terminal-driven, time-harmonic A-V
formulation in which the conductor's interior is solved and the skin profile is
an output.  It reproduces the exact Bessel round-wire resistance to `2.3e-5`
over `a / delta` from `0.25` to `4`
(`results/eddy_aphi_round_wire_20260922.json`).

On the 48 mm fin at 5 kHz, where the skin depth is `0.9346 mm`, the reference is
mesh-converged -- `maxh` `0.9`, `0.7`, `0.5 mm` give a beak fraction spread of
`1.4e-3` and R within `0.11%`:

| route | R (uOhm) | beak loss | vs A-V | tip loss | vs A-V |
|---|---|---|---|---|---|
| A-V, interior resolved | **64.729** | 0.30157 | -- | 0.08724 | -- |
| BEM-A, SIBC | 41.582 | 0.31376 | +4.04% | 0.11638 | **+33.40%** |
| PEEC, SIBC | 40.582 | 0.31504 | +4.47% | 0.11846 | **+35.79%** |

The two surface-impedance routes differ from each other by `0.41%` on the beak
loss and land on the same side of the reference by `4%`, while the tip share --
the quantity most exposed to the assumption, since the nose radius is
`0.25 mm` against a `0.93 mm` skin -- is out by a third.  R is low by `36%`.
This is exactly the failure gate 3 exists to catch: two routes that share an
assumption agreeing with each other and being wrong together.  Data:
`results/beak_fin_three_route_20260922.json`.

Two things this does NOT say, both recorded in the artifact.  The fractions are
measured differently by construction -- a volume integral of `|J|^2 / 2 sigma`
against `|K|^2 Rs / 2` weighted by perimeter -- and those coincide exactly in
the thin-skin limit.  The divergence is that limit failing rather than a
mismatch to correct away, but it has to be stated.  And the delivery gate was
accepted at 150 kHz, where `delta / thickness` is `0.043` against `0.23` here:
this locates a failure at one end and does not locate the boundary.

Reaching 150 kHz is blocked, not merely expensive.  Resolving a `0.17 mm` skin
over a 48 mm fin isotropically is out of reach -- `maxh 0.5 mm` already gives
170561 elements, 975622 DOFs and ten minutes per frequency -- and the
anisotropic boundary layer that would do it cheaply cannot be built in Netgen
6.2.2606: `Mesh.BoundaryLayer` refuses with "Call syntax has changed", and the
replacement route raises `Need to register class
netgen::BoundaryLayerParameters for Archive using std::any`.  A different
route to the same end -- meshing the conductor separately in layers, or halving
the model on its symmetry plane -- is the next move, not a bigger machine.

One false lead worth not repeating: the first 48 mm sweep gave a non-monotonic
beak fraction (0.302, 0.331, 0.302 at 5, 15, 50 kHz), which looked like a
measurement artefact because the region masks cut through elements.  Varying
the integration order from 2 to 14 settled the fraction to `1e-4` by order 6,
so the measurement was fine and the 50 kHz row was simply under-resolved --
its skin is thinner than one element.

### 150 kHz reached by revolving the section, and the assumption fails there too (2026-09-22)

The three-dimensional route could not reach the frequency the delivery gate
was accepted at.  Revolving the section removes the third dimension without
removing the physics: the beak profile becomes the meridian of a ring, the
skin is resolved by a two-dimensional mesh, and a solve takes **0.3 s** at
43k degrees of freedom against 10 minutes at 976k.  `radia.eddy_axisym_ring`
drives `radia.axifem`'s already-validated stiffness and sigma-mass operators
with a loop voltage and recovers the ring current.

The solver carries the same anchor as the three-dimensional one.  Against the
exact Bessel round-wire resistance it converges to `-0.61%` at `a / delta` of
both 6 and 12, and the residual shrinks as `1/R0` -- `-2.70%` at `R0/a = 10`,
`-0.50%` at 30 -- so it is the ring's curvature rather than the solve.

On the beak section the sweep is converged in both directions: refining
`maxh` from `0.08` to `0.05 mm` moves the 150 kHz beak share by `0.05%`, and
doubling the radius from 200 to 400 mm moves it by `0.14%`.

| frequency | skin depth | beak loss | tip loss |
|---|---|---|---|
| 5 kHz | 0.935 mm | 0.31570 | 0.09119 |
| 15 kHz | 0.540 mm | 0.34770 | 0.16105 |
| 50 kHz | 0.296 mm | 0.32490 | 0.22959 |
| **150 kHz** | **0.171 mm** | **0.30753** | **0.24504** |

Against the 2-D surface-impedance reference at 150 kHz -- the same section, the
same cuts, converged to 1024 perimeter samples, differing only in that it
replaces the interior by a Leontovich impedance:

| | beak loss | tip loss |
|---|---|---|
| interior resolved | 0.30753 | 0.24504 |
| 2-D SIBC, n=1024 | 0.34167 | 0.20231 |
| **SIBC relative** | **+11.10%** | **-17.44%** |

The reason is geometric and can be read off the fixture.  A Leontovich
impedance is the leading term of an expansion in `delta / R`, the skin depth
over the local radius of curvature.  At 150 kHz the skin is `0.171 mm` and the
tip radius is `0.25 mm` (`TIP_RADIUS` in the generator), so **`delta / R` is
`0.68`: the first curvature correction is of order one, not of order a
percent** -- and the beak root is a pair of sharp corners, where `R` is zero
and no term of the expansion applies at all.  The assumption has to fail at
those features first, and they are what a beak fin is for.  That is where the
error is largest and it changes sign against the beak share, which is what a
redistribution looks like rather than a scale error.

*(Corrected 2026-09-22: an earlier revision of this section and of the commit
that introduced it quoted the tip radius as `0.125 mm` and concluded the skin
was thicker than the tip.  The fixture's radius is `0.25 mm`; at 150 kHz the
skin is thinner than the tip radius, not thicker.  The measured shares and the
`+11.10% / -17.44%` split are unaffected -- only the stated cause was wrong,
and it is the cause that motivates the correction below.  The 5 kHz row, where
the skin is `0.935 mm`, does sit beyond the tip radius.)*

This does not say the delivery gate's *pairwise* acceptance was wrongly
computed -- BEM-A and the PEEC do agree with each other to `0.58%`, and that
was never in doubt.  It says the quantity they agree on is not the one the
conductor produces.  Data: `results/beak_section_axisym_sweep_20260922.json`.

Scope: a ring is not a straight fin, and end effects are absent by
construction, so this is the mid-span comparison and not a terminal one.  The
curvature is measured rather than assumed away, but it is not zero.

## Required next gates

1. A synthetic straight beak-fin STEP and reproducible generator now live in
   `tests/coil_from_cad/fixtures/beak_fin_straight.step` and
   `validation_test/induction_heating/make_beak_fin_step.py`. The dedicated
   straight-prism section route preserves its single solid, 0.25 mm rounded
   tip, constant section, and direct perimeter (64 lanes in the fixture test).
   This is not representative of an as-built coil. **Done on the synthetic
   fixture (2026-09-22)**: station alignment, branch geometry, and convergence
   under perimeter and axial refinement are measured above, and missing tips
   and zero-area cells are rejected per level. **Still open**: representative
   as-built CAD, which the lab has to supply; nothing here establishes that a
   real coil's section behaves like this prism.
2. Preserve the experimental path's SIBC-consistent surface-band resistance
   and transverse MNA topology when promoting it. Replace or validate the
   rectangular Ruehli self term as a nonorthogonal surface partial element;
   do not route it through production's isolated-wire Dowell/Bessel allocation
   or proximity correction. **Measured, not passed (2026-09-22)**: the branch
   mutual is the bare filamentary formula and runs `5.7%` high on touching
   bands. Replacing it means reaching the cross-section-averaged kernel, which
   today is only accessible to sub-filaments of a shared parent.
3. Compare local current, tip/root loss, terminal impedance, and field at the
   workpiece against an independent 3-D A-phi/HCurl reference over frequency,
   conductivity, geometry, and mesh sweeps. KCL alone is not an accuracy gate.
   **Reference built and validated, gate NOT passed (2026-09-22)**: at 5 kHz
   the three-dimensional A-V solve puts the surface-impedance routes `4%` out
   on the beak share, a third out on the tip share and `36%` low on R, and
   they agree with each other throughout. Revolving the section then reached
   150 kHz, where the surface-impedance description is `+11.1%` on the beak
   share and `-17.4%` on the tip share, for the geometric reason that
   `delta / R` at the tip is `0.68` and infinite at the beak root corners.
   Mesh and frequency sweeps are done;
   conductivity and geometry sweeps are not, and the terminal quantities on
   the straight fin at 150 kHz remain out of reach in three dimensions.
4. Extend the workpiece weak and strong coupling APIs to accept per-branch
   currents and per-branch induced EMFs. The present K-by-K bundle reduction
   and one-current-per-filament field projection cannot represent transverse
   branches or longitudinal redistribution.
5. Only after those gates, add an explicit `fin-surface` IH backend and
   document its contract through the owning Radia MCP manual. Existing
   rectangle/circle production behavior must remain unchanged.

## Automatic fin detection from the STEP section (2026-09-19, later)

The beak/tip regions and the workpiece-side probe are no longer fixture
constants (`x >= 1.6 mm`, `x >= 3.5 mm`, probe at `x = 5 mm`).
`radia.fin_section.analyze_section` derives them from the dense CAD outline:

- local thickness by inward-normal ray casting; fin = contiguous thin run
  (`t < 0.5 * median(t)`), merged across the rounded tip and walked back to
  the geometric root (`t < 0.9 * t_body` stops at the shoulder faces);
- root chord, fin axis, extremity; tip radius by least-squares circle fit
  (chord-polygonisation insensitive); tip arc by circle membership;
- `fin_mask` = on/beyond the root plane, `tip_mask` = within `2 r_tip` of the
  extremity along the axis, `probe_points` = perpendicular segment at
  `max(5 delta, 4 r_tip)` beyond the tip spanning the body half width.

On `beak_fin_straight.step` this reproduces the former constants
(root 1.608 mm, tip threshold 3.4995 mm, probe centre 5.0008 mm, r_tip
0.2500 mm) and is rigid-motion invariant; rectangles and circles yield no
fin.  `build_hybrid_surface_topology_from_straight_prism_step` returns the
analysis in its metadata and accepts `lane_grading="auto"`, which grades the
lanes geometrically toward each fin tip with `tip_lanes` cells on the tip
arc (2-D check: graded 128 lanes reproduce the uniform-256 R to 0.03 %,
beak fractions to 1 %).  Default lane placement is unchanged
(`"uniform"`), so previous goldens are untouched.  `tests/test_fin_section.py`
covers detection, invariance, density stability, two-fin separation and
lane grading without CAD; `sibc2d_beak_reference.py` gained an
auto-detection golden and a graded-vs-uniform consistency test.

Verification in the build123d/native-PEEC environment passed all 25 focused
tests. The CAD sampler's fifth argument is the requested point count, so the
2048-point outline path is API-correct. A graded 128-lane STEP solve produced
eight lanes on the fitted CAD tip arc and completed the 896-branch native C++
PEEC assembly/solve. During verification, centre-only region masks exposed a
3.27% lane-phase jump when one graded dual cell straddled the detected root.
Integrated metrics now use `feature_panel_weights`, the exact overlap fraction
of each perimeter dual cell with the automatically detected fin/tip interval;
the declared graded-vs-uniform tolerances then pass without relaxation.
`tip_mean_absK_over_mean` is consequently an arc-length-weighted tip mean,
not a mean over lane indices; its existing 2.0 +/- 5% golden passes under this
definition. Arc coordinates use continuous nearest-segment projection, so two
lanes cannot silently collapse onto one dense-outline sample; duplicate or
unordered lanes fail loudly. Unit tests cover dual-cell arc-length conservation,
less than 0.5% change under a half-cell lane-phase shift, and duplicate-lane
rejection. A `FinFeature.fin_mask` is a single-feature half-space convenience;
for multiple fins in the same half-space, connected-interval weights are the
required integration route.

Not covered yet: fins whose section changes along the sweep (needs
station-wise analysis and lane transition rings), separate brazed fin
solids, and fins with filleted roots (the root then lands at the end of
the thin walk, not at a corner).

## STEP -> fin PEEC for general sweeps (2026-09-19, later)

`radia.fin_sweep` generalises the straight-prism fixture route:

- `station_outlines_from_step` returns dense per-station outlines in the
  station frame using the existing CAD extractors.  `route="auto"` takes
  a single z-extruded prism through direct z-sectioning and everything
  else through the spine sectioner; the sectioner's front half was
  factored out of `_filaments_from_section_planes` as
  `coil_from_cad._section_faces_from_solid` (behaviour of the filament
  builder unchanged) so no PEEC solver is built for outline extraction.
- `align_station_origins` re-registers each station's arc-length origin to
  its predecessor by FFT cross-correlation.  The CAD sampler starts at the
  boundary point nearest +u, which jumps when a fin appears; without this
  step lane k would change material line and the longitudinal branches
  would cross the section.
- `build_fin_graph` analyses every station (`fin_section.analyze_section`),
  grades lanes toward the fin tip where a fin exists and keeps uniform
  lanes elsewhere, lifts the lanes to 3-D with the station frame and puts
  circumferential rings at every interior station.  `graded_lane_arclengths`
  gained `max_ratio` (default 8) so a fin tapering to nothing along the
  sweep never demands more lanes than the budget; `max_ratio=None` keeps
  the strict fixture behaviour.
- `validation_test/induction_heating/fin_peec_from_step.py` runs STEP ->
  graph -> experimental surface PEEC and reports R, external L and, per fin
  station, beak current/loss fractions, axial centroid, tip lanes and the
  3-D probe field (`peec_proximity._biot_savart_H`) at
  `max(5 delta, 4 r_tip)` beyond the tip.

`tests/test_fin_sweep.py` covers the frames, origin registration, a
straight sweep whose fin tapers to zero (fin detected while present, absent
after, longitudinal branches stay within 1.6x the station spacing after
registration) and a 90 degree planar sweep (rings lie in the (u, v) planes,
branch count).  A build123d-backed regression also checks that direct STEP
sections retain connected outer-wire order, SI scale and fin detection.

Native verification now covers all three CAD cases.  The 6 mm fixture gives
29.1938 micro-ohm (the previous independent result was 29.19 micro-ohm), while
the 60 mm fixture gives 322.10 micro-ohm = 5.37 milliohm/m; the earlier
28.8 micro-ohm expectation belongs to the short fixture, not the 60 mm file.
A z-directed loft with 2.5% changing cap area runs through eight sections and
the complete 3-D graph/PEEC solve (160.38 micro-ohm).  Auto routing accepts
such a tapered axial loft without requiring equal cap areas, but uses an
axial-aspect guard so a flat planar coil extruded in z is not misclassified as
a z-directed conductor.  The curved `rect_torus_lofted_united.step` therefore
uses `step_section_planes`; eight CAD sections, 416 branches and the native
solve complete (1.263 milliohm, 141.51 nH; geometry-path smoke evidence, not a
fin accuracy reference).  The runner reports requested `tip_lanes` separately
from `tip_overlap_majority_panels`, since conservative dual-cell overlap can
make the latter larger than the placement target.

Still outside the automation: separate brazed fin solids (unite in CAD),
filleted roots (root lands at the end of the thin walk), the area-outlier
filter in `_section_faces_from_solid` (drops stations whose section area
differs >30 % from the median, i.e. very large fins appearing mid-sweep),
and the production hook (`calc_inductance.py` / IH operator assembler
still use the series bundle).

## Production boundary: `--coil-solver fin-surface` (2026-09-19, later)

Design in `docs/induction_heating/FIN_PEEC_PRODUCTION_BOUNDARY.md`.  The IH
assembler has no electromagnetic solve of its own; it calls
`calc_inductance.run_inductance` and consumes `qsurf.sol`.  Hence the fin
PEEC enters production as a third coil solver:

- `radia.fin_surface_solver.solve_fin_surface` is the single implementation
  (STEP -> `fin_sweep` graph -> surface PEEC -> branch currents -> per-station
  fin metrics); `fin_peec_from_step.py` now calls it, so the runner and the
  production backend report identical numbers.
- `calc_inductance._solve_coil_fin_surface` returns the existing
  `coil_data` contract with `source_type="filament"`, one two-point polyline
  per surface branch (axial and circumferential) and its complex current.
  Downstream A/B/phi_inc/msh consumers integrate paths independently, so
  weak BEM-SIBC coupling works unchanged.  Strong coupling is rejected
  (the coupled PEEC solver needs the K x K lane bundle; per-branch EMF is
  gate 4).  New CLI: `--fin-n-lanes --fin-n-stations --fin-n-outline
  --fin-lane-grading --fin-tip-lanes --fin-route`.
- `IHOperatorAssemblyOptions.coil_step_solver = "peec" | "fin-surface"` plus
  `fin_*` options; `_unit_current_argv` builds the CLI (unit tested).  The
  native config keeps `eddy_solver="peec"` (MEX contract) and records the
  variant in `eddy_method` and `geometry.coil_backend`; `unit_current`
  carries `fin_metrics` / `fin_sweep`.  The C++ runtime is untouched.

Boundary acceptance (interface, not accuracy) per the design document:
(1) `--coil-solver fin-surface --coil-only` on the straight fixture equals
`fin_peec_from_step.py`; (2) weak coupling with a workpiece `.vol` produces a
finite, non-negative `qsurf.sol`; (3) the assembler's `native_ih.json` with
`coil_step_solver="fin-surface"` passes `validateIHNativeConfig`; (4) record
the corner-region loss share and centroid of `qsurf` for `peec` vs
`fin-surface` on the same workpiece.  (1)-(4) were not executed in the
authoring environment (no NGSolve / native kernel); they were executed
on 2026-09-21 and the result is
`results/fin_surface_boundary_acceptance_20260921.json`.  (1) and (2)
pass, (3) satisfies the contract with the MATLAB validator itself not
run, and (4) is recorded for `fin-surface` only.

Two findings came out of running them.  The straight beak-fin fixture
is one open prism, so its field has no scalar potential on a nearby
surface and the surface-Poisson gate refuses it; the residual is
mesh-converged at 19.3%, so this is the open current path and not
discretisation, and criterion (2) is therefore recorded on the tracked
closed conductor.  On that case the fin-surface source passes the
shared SIBC reciprocity gate at 8.2% while the series-filament bundle
is refused at 13.2%, flat under both workpiece and perimeter
refinement and unchanged by the proximity iteration, although the two
surface losses agree to 1.3%.  That is an interface observation, not a
fin accuracy result.
