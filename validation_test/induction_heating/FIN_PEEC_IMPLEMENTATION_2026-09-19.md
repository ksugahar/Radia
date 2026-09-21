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

## Required next gates

1. A synthetic straight beak-fin STEP and reproducible generator now live in
   `tests/coil_from_cad/fixtures/beak_fin_straight.step` and
   `validation_test/induction_heating/make_beak_fin_step.py`. The dedicated
   straight-prism section route preserves its single solid, 0.25 mm rounded
   tip, constant section, and direct perimeter (64 lanes in the fixture test).
   This is not representative of an as-built coil. Obtain representative CAD,
   then check station alignment, branch geometry, and convergence under
   perimeter/axial refinement. Reject missing tips and zero-area cells.
2. Preserve the experimental path's SIBC-consistent surface-band resistance
   and transverse MNA topology when promoting it. Replace or validate the
   rectangular Ruehli self term as a nonorthogonal surface partial element;
   do not route it through production's isolated-wire Dowell/Bessel allocation
   or proximity correction.
3. Compare local current, tip/root loss, terminal impedance, and field at the
   workpiece against an independent 3-D A-phi/HCurl reference over frequency,
   conductivity, geometry, and mesh sweeps. KCL alone is not an accuracy gate.
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
