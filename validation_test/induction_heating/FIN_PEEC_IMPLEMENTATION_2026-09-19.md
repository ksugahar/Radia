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

## Required next gates

1. A synthetic straight beak-fin STEP and reproducible generator now live in
   `tests/coil_from_cad/fixtures/beak_fin_straight.step` and
   `validation_test/induction_heating/make_beak_fin_step.py`. The dedicated
   straight-prism section route preserves its single solid, 0.25 mm rounded
   tip, constant section, and direct perimeter (64 lanes in the fixture test).
   This is not representative of an as-built coil. Obtain representative CAD,
   then check station alignment, branch geometry, and convergence under
   perimeter/axial refinement. Reject missing tips and zero-area cells.
2. Implement a consistent surface-current discretization with passive SIBC
   Gram matrix and validated nonorthogonal partial elements. Avoid adding the
   existing isolated-wire Dowell/Bessel or proximity correction on top.
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
