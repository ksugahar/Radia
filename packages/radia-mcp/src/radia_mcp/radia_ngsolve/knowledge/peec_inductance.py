"""PEEC-inductance knowledge (mcp-server-radia-ngsolve, public tool).

Panel mode: PEEC inductance (coil only, STEP).  Computes L_coil / R_coil
from a coil solid STEP file via filament-based Biot-Savart + loop-bundle
PEEC solve.  No workpiece, no BEM, no FEM mesh -- this is the lightest
path in the IH panel family and the one a user reaches for when they
want just a quick "what is the inductance of this coil at this
frequency?" answer.

Filaments are placed on the cross-section **perimeter only** (thin-skin
regime, d/δ >= 3).  Use ``n_peri`` filaments around the arc-length
perimeter of each cross-section; no interior volume grid.

Source: src/radia/panels/calc_inductance.py,
        src/radia/coil_from_cad.py (STEP path -- the only input path
        since 4.13.0; .jou explicit-centerline input was retired
        per CLAUDE.md No-Fallbacks),
        src/radia/ih_design.py + the Induction Heating Simulink block
        -- Method = "PEEC inductance (coil only, STEP)".
"""


PEEC_IND_OVERVIEW = """
# PEEC inductance (coil only) — overview

## What it solves

Given a coil solid STEP file + conductivity + frequency,
returns:
  * L_coil_nH  — inductance in nanohenry (includes internal-inductance
    drop at high frequency via the Bessel `Z_cyl(ω)` formulation)
  * R_coil_mOhm — full round-wire **AC resistance** (Bessel
    `Z_cyl(ω) · L_filament`).  Fixed 2026-05-19: the panel now
    computes per-filament `Zs_fil_k = n_peri · (Z_cyl(ω) − R_DC/m) ·
    L_k` and passes it to `solve_loop_bundle`, so R recovers R_DC at
    low frequency AND the SIBC `(1+j)/(σδ·P)` asymptote at high
    frequency.  Cross-section radius is read from the topology
    (`cross_section_radius_m_mean` or derived from cell_wh).  Falls
    back to R_DC-only with a WARNING if the radius cannot be inferred
    (e.g. unusual rectangular bars without cross_section_area_m2_mean).
    See docs/esim/R_MISMATCH_PEEC_VS_BEMA.md.
  * Z_coil = R + j·ω·L (full complex impedance from the Bessel
    formulation)

No workpiece, no BEM, no FEM mesh.  The fastest IH panel mode; use
when the question is "what is the coil's L / R at this freq?" and
you do NOT need heating, losses in a workpiece, or field maps.

## When to pick this mode vs the others

| Question                                   | Use mode                             |
|--------------------------------------------|--------------------------------------|
| L, R of coil alone (no workpiece)          | **PEEC-inductance** (this page)      |
| P_wp (workpiece heating), 1-way            | PEEC+BEM                             |
| L + P_wp + P_coil exact, small mesh        | FEM A-V (calc_fem_coilmesh.py)       |
| Nonlinear steel workpiece, Karl iteration  | FEM + ESIM                           |

## Filament placement: perimeter only

Filaments sit on the cross-section OUTER BOUNDARY only — no interior
volume grid.  Parameter ``n_peri`` (default 16) chooses how many
filaments to spread around the perimeter, arc-length-equidistant.

Physical justification: at the target IH operating regime (d/δ >= 3,
1 kHz – 1 MHz for Cu / Al / brass), current crowds to the surface.
Volume-grid filament schemes (FastHenry nwinc/nhinc) spend DOFs on
interior cells that carry no current at these frequencies.  Perimeter
placement spends them all where the current actually flows.

Per-filament SIBC is now wired (2026-05-19): `_solve_coil_peec`
computes `Zs_fil_k = n_peri · (cylinder_ac_impedance(a_eq, σ, ω) −
R_DC_per_m) · L_filament_k` and passes it as the `Zs_fil` argument
of `solve_loop_bundle`.  The result is the full round-wire Bessel
AC impedance bundle (R_DC at ω → 0, `(1+j)/(σδ·2πa)` per unit length
at high ω).

Dowell for rectangular bars and ESIM for nonlinear steel are NOT
wired into this dispatcher yet — round / quasi-round cross sections
get the proper Bessel R; high-aspect-ratio ribbons fall back to the
equivalent-circle approximation (within ~10 % for square-like, errs
for thin ribbons).

## Input dispatch

``calc_peec_inductance.py --peec-step <path>`` accepts:
  * ``.step`` / ``.stp``  — geometry file; centerline auto-extracted from
                            solid topology (see ``centerline_extraction``).
  * other                 — raises ValueError.

Coordinates in METRES (CLAUDE.md "Unit System Policy: Radia always uses
meters").  Cubit set-up: ``set unit-system mks`` BEFORE generating the
coil geometry.

Since 4.13.0 the legacy ``.jou`` explicit-centerline input was retired
(CLAUDE.md No-Fallbacks): the STEP B-Rep is the canonical PEEC source
and centerline + cross-section are auto-extracted from it.  The retired
.jou path opened a silent inconsistency window when the .jou contained
``volume all scale K`` -- parse-only PEEC and Cubit-derived paths
disagreed on geometry size.

## Filament placement evolution (4.13 → 4.15)

  * v4.13.0: STEP-only centerline + constant equivalent-circle filaments
  * v4.14.0: per-station face sampling for NON-united multi-loft
             (variable cross-section: tapered, varying-width pancake)
  * v4.15.0: UV-map sampling on single-piece lateral surface (most
             general; arbitrary cross-section + arbitrary spine)

See the ``filament_dispatch`` topic for the 3-tier order, triggers,
and the recipe for keeping coil STEPs in the most general form
(do not ``unite`` before export).
"""


PEEC_IND_CENTERLINE = """
# PEEC inductance — STEP centerline extraction

**v4.48.1 (2026-05-15)**: ``extract_centerline_from_step`` in
``src/radia/coil_from_cad.py`` uses **classification-based single
dispatch** (CLAUDE.md "STEP-Only Centerline: Auto-Detect or Fail" +
"No Fallbacks - Fail Fast, Fail Loud").  Five positive-match
predicates route to exactly ONE extractor; the chosen extractor's
failures propagate as ``ValueError`` -- there is NO ``try/except``
cascade across paths, NO caller-supplied centerline JSON
(``path_points_m`` was removed), and NO silent "this path's output
looks wrong, try the next one" behaviour.

If auto-detect cannot recover a covering centerline, the panel /
CLI **raises with a HINT pointing at CAD regeneration or BEM-A
switch** -- the user fixes the geometry instead of papering over it.

## Predicate 1. Loft of profiles (multi-turn coil, tight pancake spiral)

**Trigger**: solid has >= 5 planar faces of consistent area
(``_collect_loft_cross_sections`` filters planar faces, checks
median-area +- 50% band).

**Algorithm**:
  1. Enumerate planar cross-section faces (``GeomType.PLANE``)
  2. Dedupe near-duplicate centroids (shared end-cap between
     adjacent loft volumes appears twice, once per owning volume --
     merge if < 10% of equivalent radius apart)
  3. Compute centroid of each unique cross-section
  4. Chain centroids via **nearest-neighbor with tangent continuity**:
     endpoint = point with largest (2nd-NN / 1st-NN) distance ratio;
     walk from endpoint, prefer forward-aligned next neighbor
  5. Return ordered polyline + equivalent-circle radius from mean area

**Validation** (Kubota's ``3turncoil.stp``, 4 MB, 382 loft volumes,
re-verified 2026-04-28 after .jou retirement):
  * STEP-only:  L = 426.265 nH (golden 426.25 nH, error 0.004%)
  * topology extraction: ~9 s

**Why this is faster than the old spine method**: no ``section()``
calls -- the cross-section faces ARE the sections.  Old code called
``section(solid, Plane)`` once per midpoint (100 calls = 100 OCC
boolean cuts = minutes on 10 MB STEPs).

## Predicate 2. United multi-turn pancake (CIRCLE-edge stations)

**Trigger** (``_collect_circle_edge_centers``): solid is a united
multi-turn pancake where boolean ``unite`` consumed the planar
end-caps but cross-section CIRCLE edges remain (split into 2
semicircles per cross-section by Cubit's unite operation).

**Algorithm** (``_centerline_from_circle_edge_centers``): group
CIRCLE edges by arc_center, dedupe the semicircle pairs into one
centre per cross-section, chain via NN + tangent continuity, return
polyline + median CIRCLE radius.  Handles Kubota's united
``3turncoil.stp`` class when end-caps have been consumed by
``unite``.

## Predicate 3. Single-loop revolution-sweep coil (gapped torus, rect-profile bend)

**Trigger**: solid has one or more revolution-type surfaces
(TORUS / CYLINDER / CONE / REVOLUTION) sharing a common axis +
at least one PLANE end-cap face.

**Algorithm** (``_centerline_from_revolution_sweep``):
  1. Enumerate ``build123d`` faces, map each to OCP BRepAdaptor_Surface
  2. Collect revolution-type surfaces sharing a common axis (location
     + direction match within 1e-6)
  3. Sweep angle = max(U) − min(U) across the union of U parameter
     intervals.  Merging splits a single surface into multiple faces
     each covering part of [u_min, u_max]; the union recovers the
     full sweep.
  4. Pick the smallest-area PLANE face as the cross-section cap
  5. Spine radius = distance from cap centroid to the sweep axis
  6. Cross-section area = cap.area → equivalent-square side for
     downstream ``filaments_from_polyline``

**Topology → faces seen** (verified 2026-04-21):

| Profile            | Lateral faces           | Cap faces            |
|--------------------|-------------------------|----------------------|
| Circle             | 4× TORUS                | 2× PLANE (disks)     |
| Rect 6×4mm         | 4× CYLINDER + 2× PLANE* | 2× PLANE (rect)      |
| Polygon N-sided    | N × CYLINDER            | 2× PLANE (poly)      |

\\* Rect profile produces 2 extra PLANE faces for top / bottom flat
sides of the sweep.  These are LARGER than the end caps, so sorting
by area and picking the smallest correctly identifies the cap.

**Validation** (all coil topologies STEP-only, re-verified 2026-04-28):

| Coil                           | Expected L | Got L     | Error    |
|--------------------------------|-----------|-----------|----------|
| Circular 355° torus            | 88.55 nH  | 85.100 nH | -3.9 %   |
| Rect 6×4 mm 355° torus         | ~88 nH    | 88.15 nH  | reference|
| 3-turn loft (Kubota 3turncoil) | 426.25 nH | 426.265 nH| +0.004 % |

## Predicate 4. OPEN coil with leads (longest open lateral rim)

**Trigger**: predicates 1-3 did not match AND
``coil_topology.extract_coil_topology(solid).is_open == True`` (the
solid has 2 cap faces detected -- BSPLINE-lofted "arc + leads" or
free-form helix).

**Algorithm** (``_centerline_from_open_spine``): find the LONGEST
OPEN lateral rim edge (``is_closed == False``, EXCLUDES cross-
section CIRCLE edges); sample it at ``n_segments + 1``
arc-length-equidistant points; section at midpoints to extract
cross-section centroids and areas.

**Why this path is critical** (keiko 2026-05-15 incident):
BSPLINE-lofted coils whose lateral is split into 2 equal halves at
the z-equator fail Predicate 1 (dominance ~0.5 < 0.8), have no
multi-station planar faces (Predicate 1 of filaments_from_step has
2 BSPLINE half-faces total), have no CIRCLE-edge stations
(Predicate 2 here -- the only CIRCLE edges sit at the caps, not
along the spine), and have no revolution-type lateral
(Predicate 3 here).  They DO have 2 PLANE cap faces, so
``extract_coil_topology`` returns ``is_open=True`` and routes here.
The long lateral rim correctly traces the LEAD extensions even
when the caps sit far from the round arc (e.g. caps at
``(+-13, +50, 0) mm`` with the arc at ``R ~ 30 mm`` -- the rim
follows ``arc + lead`` faithfully).

## Predicate 5. CLOSED full-revolution coil (no caps)

**Trigger** (``_centerline_from_topology_spine``): predicates 1-4 did
not match AND ``coil_topology.extract_coil_topology(solid).is_open ==
False`` (no cap faces -- typical for CLOSED full 360deg loop, e.g.
``ih_closed_torus_coil.step``).

**Algorithm**: ``coil_topology.generate_spine`` produces the full
360deg arc at the bbox-derived ``R_spine``; section once at the
midpoint to recover cross-section area; equivalent-square side from
the section area.

**CLOSED-only guard**: ``_centerline_from_topology_spine`` raises if
called with ``topo.is_open == True`` (programming-error indicator,
NOT a soft-fallback signal).  This guard catches keiko-class OPEN-
with-leads geometries that would otherwise get a bbox-radius planar
arc bypassing their leads.

## Units: METRES ONLY (4.12.0, 2026-04-28)

Per CLAUDE.md "Unit System Policy: Radia always uses meters", the
PEEC inductance pipeline accepts METRE inputs only.  No auto-
detection (the previous bbox heuristic was a CLAUDE.md No-Fallbacks
violation that masked Cubit STEP-writer unit-declaration bugs and
silently produced 10x / 1000x wrong inductances).

``build123d.import_step`` passes raw numeric coordinates through
without applying the STEP ``CONVERSION_BASED_UNIT`` declaration.
Cubit's STEP writer often declares ``MILLIMETRE`` while writing
metre-valued coordinates; that is OK because we ignore the
declaration entirely.  The contract is simple:

  * STEP coordinate values must be metres
  * .jou ``move Surface`` x/y/z must be metres
  * .jou ``create surface circle radius`` must be metres

``filaments_from_step(...)`` and ``parse_jou_centerline(...)``
default to ``cad_units_per_meter=1.0``.  Override only if you
KNOW your input is in another unit (e.g. legacy mm-mode Cubit
journals).

Cubit set-up to produce metre-unit output:
  ``set unit-system mks`` BEFORE building the coil geometry.
"""


PEEC_IND_FILAMENT_DISPATCH = """
# PEEC inductance — filament placement dispatch (n_peri perimeter mode)

Since **v4.48.1** (2026-05-15), ``filaments_from_step`` dispatches
filament placement through **classification-based single dispatch**
(CLAUDE.md "STEP-Only Centerline: Auto-Detect or Fail" + "No
Fallbacks - Fail Fast, Fail Loud").  Each predicate is a positive
match for a specific solid topology; the chosen path's failures
PROPAGATE as ValueError -- there is NO try/except cascade across
paths, NO ``path_points_m`` user override, and NO silent fall-through
to a different sampler.

Order in ``src/radia/coil_from_cad.py``:

## Path 1 — UV-map sampling (most general single-piece lateral)

**Predicate** (``_find_lateral_surface``): all of
  * solid has ONE dominant lateral surface (BSPLINE / CYLINDER /
    TORUS / REVOLUTION / EXTRUSION) with
    ``largest_area / total_area >= 0.8``
  * ``len(candidates) <= 3`` (defensive multi-piece guard)
  * the lateral has at least one closed UV axis
    (``IsUClosed or IsVClosed``)

When the predicate returns a face, downstream UV sampling MUST
succeed -- the UV-closure check is part of the predicate, so Path 1
does NOT carry a try/except.

**Algorithm** (``_sample_lateral_surface_uv`` →
``_filaments_from_lateral_surface_uv``):
  1. Sample on (n_stations x n_peri) parameter grid; the closed
     axis is the perimeter direction (cell-centred ``(k+0.5)/n_peri``),
     the open axis is the spine direction (``0..n_stations-1``).
  2. Build ``n_peri`` filaments from the (n_stations x n_peri x 3) grid.
  3. Per-segment cell area via 3D-shoelace on the perimeter polygon
     at each station -- **variable cross-section is automatic**, no
     equivalent-circle assumption.

## Path 2 — Per-station face sampling (NON-united multi-loft)

**Predicate** (``_try_extract_loft_with_profile``): solid is a
NON-united multi-loft with planar end-cap faces surviving in the
STEP (loft volumes preserved, no boolean ``unite``).

**Algorithm** (``_filaments_from_per_station_faces``): sample each
cross-section face's outer wire at ``n_peri`` arc-length-equispaced
points; project to the parallel-transport frame; filament k
connects sample k of station i to sample k of station i+1.
Handles RECTANGULAR / POLYGON / arbitrary cross-section + any
spine-direction shape variation.

## Path 2b — Circle-edge stations (united multi-turn pancake)

**Predicate** (``_collect_circle_edge_centers``): solid is a united
multi-turn pancake where the boolean unite consumed planar end-caps
but the cross-section CIRCLE edges remain (split into 2 semicircles
per cross-section by Cubit's unite).

**Algorithm**: group CIRCLE edges by arc_center, dedupe the
semicircle pairs into one centre per cross-section, chain via NN +
tangent continuity.  Handles Kubota's 3turncoil united.stp class.

## Path 2c — Section-planes (united multi-loft, non-circular)

**Gate** (predicate-style): ``n_plane >= 4 AND n_revolution_like == 0``.
**Trigger** (``_filaments_from_section_planes``): united multi-loft
with NON-circular cross-section (rect / polygon / shape transition).
Slow (one OCC section per station) but the only path that handles
rect-torus and similar united-rect-cross-section coils.

## Path 3 — Constant equivalent-circle (legacy catch-all)

**Trigger**: none of the above predicates matched.  The centerline
comes from ``extract_centerline_from_step``'s classification
dispatch (see "centerline" topic for the 5-predicate breakdown),
then ``filaments_from_polyline`` places ``n_peri`` filaments on the
equivalent-circle perimeter at the mean cross-section area.

After v4.48.1 the centerline coming back from
``extract_centerline_from_step`` either (a) covers the conductor
bbox -- success, or (b) fails the
``_check_filaments_cover_solid_bbox`` sanity check -- raises with a
diagnostic and a HINT pointing at CAD regeneration or BEM-A switch.
**No silent wrong-spine output** (the regression that caught keiko's
``1turn_coil_loft_outsideline.step`` on 2026-05-15).

## Verification (2026-04-28 on 100号機)

| Coil topology                        | Tier hit  | L (nH)    | vs golden  |
|--------------------------------------|-----------|-----------|------------|
| shipped torus (4 TORUS, gapped)      | 3 (sweep) | 85.100    | 85.10 OK   |
| synth multi-loft 1-loop, 20 stations | 2         | 79.582    | -6 % vs 88.55 (20-pt 355° polygon under-approx; same as 4.13.0 baseline) |
| 3-turn pancake (Kubota 3turncoil)    | 3 (loft)  | 426.265   | 426.25 (+0.004 %) |

UV-map (Tier 1) falls through on 4-TORUS gapped torus because both U
and V are closed AND ``len(candidates) > 3`` — the existing
``_centerline_from_torus_sweep`` analytical path takes over.  No
regression for shipped golden cases.

## Why the 3-tier order matters

  * Tier 1 (UV-map) is the SINGLE most general sampler but requires
    a single-piece lateral.  Multi-loft and united multi-turn pancake
    STEPs trip the fragmentation guard and fall through.
  * Tier 2 (per-station faces) needs the loft volumes preserved
    (no ``unite``).  Cubit + Cubit's STEP writer can produce either
    form depending on the export script; if you ``unite`` before
    export, Tier 2 is unreachable.
  * Tier 3 is the universal fallback; equivalent-circle assumption
    breaks for severely non-circular cross-sections.

## Recipe: keep coil STEPs in their MOST GENERAL form

  * Loft of profiles → DO NOT ``unite`` before STEP export.  Tier 2
    will sample each end-cap face faithfully.
  * Single sweep / extrusion → leave the lateral surface intact (no
    ``regularize`` that fragments it).  Tier 1 will UV-sample.
  * If you must ``unite`` (e.g. for downstream Cubit meshing), expect
    Tier 3 fallback and accept the equivalent-circle approximation.

## Unsupported topology (DO NOT silently trust the result)

The 3-tier dispatch covers the **smooth-coil regime**: continuous
spine + continuous cross-section + 2-port + non-magnetic conductor +
``d/δ ≥ 3``.  Outside that regime, Tier 3 (equivalent-circle)
"runs but is not physical".  The known unsupported / inaccurate
cases:

| Case                                              | Behaviour                  | Fix path |
|---------------------------------------------------|----------------------------|----------|
| Cross-section **discontinuity** along the spine (circle→rect, sudden step in area) | Tier 2 connects sample-k-station-i to sample-k-station-i+1 across the discontinuity → non-physical filament path; Tier 3 averages the area → loses detail | Path 2c partial coverage (smooth shape transitions OK; sharp area steps still produce non-physical filament paths) |
| **United** multi-loft with non-circular cross-section (rect / polygon / shape-transition) | **v4.24.0+: Path 2c (`_filaments_from_section_planes`)**: sections the solid at n_stations along a rotation-axis-derived spine; gives 138 nH on the synth rect torus golden (was NaN before v4.24.0).  Pre-v4.24.0 fell to Tier 3 (longest-edge equivalent-circle, returned NaN for rect cross-section). | Shipped v4.24.0; see `project_peec_phase_c_heavy_shipped_v4_24_0.md` |
| **T-junction / Y-branch** bus bar                  | single-spine polyline assumption breaks; ``filaments_from_polyline`` cannot represent topology | Multi-spine PEEC re-architecture (no concrete plan) |
| **>2-port** coil (multi-tap)                       | ``add_port(n1, n2)`` is single-pair; cannot express 3+ port impedance matrix | Multi-port PEEC topology extension |
| **Closed loop** with no external terminals         | port assignment fails; cannot define the source current | Out of scope; PEEC inherently needs at least 1 driving port |
| **Litz wire** (N strands per cross-section)        | perimeter is sampled as a single outer wire; bundle inter-strand mutual L is not represented | Bundle-PEEC (multiple parallel filament sets per spine station) |
| **Magnetic conductor** (μ_r > 1, e.g. iron busbar) | SIBC ``Z_s = (1+j) / (σ δ)`` is for non-magnetic; need ``Z_s = (1+j) sqrt(ω μ / (2 σ))`` with the correct μ | One-line patch in ``rad_peec_surface_impedance.cpp`` after deciding the σ + μ_r input UX |
| **Thick-skin** ``d/δ < 3``                         | perimeter-only filaments miss the interior current distribution | Out of scope for this panel; switch to FEM A-V (volumetric coil mesh) |
| **Self-intersecting** STEP                         | OCC topology load undefined behaviour | Reject at parse time |
| **Non-manifold** STEP                              | ``build123d.import_step`` partial-parses or silently drops shells | Add a manifold-check before dispatch |

**Action when the user asks for any of the above**: do NOT silently
fall back to Path 3. Either (a) refuse + explain the limitation
(if the inaccuracy is unacceptable), or (b) accept Path 3 with a
WARN log line that names the specific limitation (e.g. ``WARN:
united multi-loft falling back to equivalent-circle, expect -1 to
-5 % on L``).
"""


PEEC_IND_STEP_AUTHORING = """
# PEEC inductance -- authoring auto-detect-friendly STEPs

The v4.48.1 STEP-only centerline policy auto-detects the spine from
the solid's BREP topology without any user-supplied JSON.  This
section is the **author's recipe**: which Cubit / build123d /
external CAD operations PRODUCE a STEP that one of the 5 predicates
(Loft / Circle-edge / Revolution+Plane / OPEN / CLOSED) cleanly
matches, and which operations TRIP one of them.

## Quick decision table

| Author tool | Operation | Predicate hit | Notes |
|-------------|-----------|---------------|-------|
| Cubit | ``sweep surface ... axis ... angle A < 360`` (circle profile) | 3 (revolution-sweep) | TORUS lateral + 2 PLANE caps; gapped torus |
| Cubit | ``sweep surface ... axis ... angle 360``               | 5 (CLOSED) | No caps; ``coil_topology`` returns ``is_open=False`` |
| Cubit | ``sweep surface ... axis ... angle A < 360`` (rect profile) | 3 (revolution-sweep) | CYLINDER lateral + PLANE caps (2 top/bottom + 2 ends); cap area filter picks the smallest |
| Cubit | multi-station ``create surface ...`` + ``create volume loft surface 1 2 ... N``, NO ``unite`` | 1 (loft cross-sections) | N planar end-caps survive; needs N >= 5 |
| Cubit | multi-station loft + ``unite``, circular cross-section | 2 (CIRCLE-edge stations) | Cubit unite consumes end-caps but keeps cross-section CIRCLE arcs |
| Cubit | multi-station loft + ``unite``, NON-circular cross-section | 2c (section-planes) | Slow (one OCC section per station) but works |
| Cubit | single multi-section loft ``create volume loft surface 1 2 ... N`` | 1 or single BSPLINE lateral -> Path 1 UV | Smoother than chain-of-pairwise-lofts; preferred for mesh quality |
| build123d | ``revolve(profile, axis, angle=A)``                       | 3 (revolution-sweep) | Same as Cubit sweep |
| build123d | ``sweep(profile, path)`` with curved path                 | 4 (OPEN) | Single BSPLINE lateral if profile is closed; if profile vertex aligns across stations, lateral is dominant -> Path 1 |
| build123d | ``loft([profile_1, ..., profile_n])`` with consistent vertex count | 1 or Path 1 UV | Lofts preserve cross-section faces unless followed by ``Compound`` flattening |
| FreeCAD Part Design | ``Pad`` / ``Revolution`` | 3 or 5 | Same topology as Cubit sweep |
| Solidworks export AP214 | revolved sweep | 3 or 5 | Cap area heuristic picks the smallest; works |

## Concrete Cubit recipes (gapped torus, circle profile)

```
reset
set unit-system mks                 # METERS only (CLAUDE.md unit policy)
brick x 0.001 y 0.003 z 0.0001     # 1 mm x 3 mm rect profile -- NO, use circle below
delete volume all
create surface circle radius 0.0015 zplane    # 1.5 mm wire radius, in z-plane at origin
move surface 1 x 0.030                          # 30 mm arc radius
sweep surface 1 axis 0 0 0 0 0 1 angle 355      # 355 deg around z
# Result: 4x TORUS lateral + 2x PLANE caps -> Predicate 3 hit
export step "gapped_torus_circle.step" overwrite
```

## Concrete Cubit recipes (multi-turn pancake spiral, circle profile)

```
reset
set unit-system mks
# author N cross-section circles at successive (R, theta, z) stations
#define a Python loop in the .jou: for i in range(N): create surface circle ...
# then:
create volume loft surface 1 2 3 ... N            # NON-united N-station loft
# NO ``unite`` -> N planar end-caps survive -> Predicate 1 hit
export step "pancake_N_loft.step" overwrite
```

If you must ``unite`` (e.g. for downstream Cubit meshing), the
end-caps are consumed but cross-section CIRCLE arcs remain split
into 2 semicircles per cross-section -> Predicate 2 (CIRCLE-edge
stations) hits.

## Concrete build123d recipe (sweep along a curved path)

```python
import build123d as bd

# circular wire profile
profile = bd.Plane.XY * bd.Circle(0.0029)        # 2.9 mm radius

# curved spine in 3D (arc + lead)
path = bd.Plane.XY * (
    bd.Line((-0.012, 0.050, 0), (-0.012, 0.025, 0))   # left lead down
    + bd.ThreePointArc(
        (-0.012, 0.025, 0), (0, -0.030, 0), (0.012, 0.025, 0))
    + bd.Line((0.012, 0.025, 0), (0.012, 0.050, 0))   # right lead up
)

coil = bd.sweep(profile, path)
bd.export_step(coil, "coil_clean.step")
```

The resulting STEP has ONE dominant BSPLINE lateral + 2 PLANE caps
-> Predicate 1 hits (single-piece lateral, IsUClosed=True for the
circular cross-section direction).

## Anti-patterns (PRODUCE STEPs that the auto-detect rejects)

  * **Lateral split into 2 equal halves at the z-equator**.  Symptom:
    keiko 2026-05-15 ``1turn_coil_loft_outsideline.step`` -- the loft
    cross-section had no consistent vertex alignment, the resulting
    lateral surface split into ``BSPLINE_top`` and ``BSPLINE_bottom``
    of nearly equal area at z=0.  ``_find_lateral_surface`` rejects
    (dominance ~0.50 < 0.80), Predicate 1 misses; the solid has no
    multi-station planar faces, Predicate 1 of the centerline misses;
    no revolution surfaces, Predicate 3 misses; ``is_open=True``
    routes correctly to Predicate 4 (OPEN longest open-edge) which
    handles it.  **But**: if the lateral split is the only issue,
    re-author with consistent vertex alignment (e.g. specify
    ``vertex_side='inner'`` or ``'outer'`` in
    ``gen_1turn_coil_loft.py``) to get a single-piece BSPLINE
    lateral -> Predicate 1 cleanly hits (smoother sampling).

  * **Pairwise loft chain** (``loft surface 1 2; loft surface 2 3;
    ...; unite``).  Symptom: chain of N-1 C^0 creases at every
    cross-section -> tetmesh slivers (radia panel-samples ``3turnCoil_work.jou``
    measured 2026-05-08: scaled_jac min=0.219 with chain vs 0.349
    with single multi-section loft).  Prefer single multi-section
    loft ``create volume loft surface 1 2 ... N`` for both mesh
    quality AND clean Predicate 1 detection.

  * **Hardcoded entity IDs in .jou**.  CLAUDE.md "Journal File
    Portability Policy": humans pick IDs visually in the GUI, but
    LLM-authored .jou must identify entities by geometric properties
    (area, volume, centroid).  Hardcoded IDs change between Cubit
    versions and after imprint/merge.

  * **Non-manifold STEP** (open shells, dangling faces).
    ``build123d.import_step`` partial-parses; predicates miss
    silently.  Add a manifold check OR fail at import.

  * **Self-intersecting STEP** (rare; usually the result of swapping
    a degenerate loft profile).  OCC topology load is undefined
    behaviour; reject at parse.

## How to verify your STEP is auto-detect-friendly BEFORE running the panel

A 10-line build123d probe (also exposed by MCP tool
``peec_inductance(topic="centerline")`` -> "Predicate N hits when..."):

```python
import build123d as bd
from radia.coil_from_cad import extract_centerline_from_step
from radia.coil_topology import extract_coil_topology

step = "my_coil.step"
solid = bd.import_step(step)

bb = solid.bounding_box()
print(f"bbox: x=[{bb.min.X*1e3:+.1f},{bb.max.X*1e3:+.1f}] mm")
print(f"      y=[{bb.min.Y*1e3:+.1f},{bb.max.Y*1e3:+.1f}] mm")
print(f"      z=[{bb.min.Z*1e3:+.1f},{bb.max.Z*1e3:+.1f}] mm")

face_kinds = {}
for f in solid.faces():
    k = str(f.geom_type).split(".")[-1]
    face_kinds[k] = face_kinds.get(k, 0) + 1
print("face kinds:", face_kinds)

topo = extract_coil_topology(solid)
print(f"is_open={topo.is_open}  sweep_deg={topo.sweep_deg:.1f}")

# Final answer: which predicate will hit?
path_m, w, h = extract_centerline_from_step(step, 100, 1.0)
print(f"centerline ok: {len(path_m)} pts, equiv R = "
      f"{(w*h).mean()**0.5*1e3/3.14**0.5:.2f} mm")
print(f"bbox match:  filament y in "
      f"[{path_m[:,1].min()*1e3:+.1f}, {path_m[:,1].max()*1e3:+.1f}] mm")
```

If the print of the filament y range covers the solid bbox y range
(within ~ wire-radius slack), auto-detect is healthy.  If not, either
the STEP is anti-pattern (re-author per the recipes above) or you've
found a new geometry class -- file a bug report with the .step.
"""


PEEC_IND_JAPANESE_PATH = """
# PEEC inductance — Japanese / Unicode path support

Verified 2026-04-21 on LAB + 100号機 with paths like
``C:\\temp\\日本語テスト\\3turncoil.stp``:

| Layer                              | Status | Mechanism                               |
|------------------------------------|--------|------------------------------------------|
| Python argparse / os.path          | ✅     | Python 3.12 Unicode-aware                |
| build123d ``import_step``          | ✅     | pythonocc-core wide API on Windows      |
| netgen.occ ``OCCGeometry``         | ✅     | OCC wide API                             |
| Cubit plugin write (.vol / .json)  | ✅     | ``utf8_path.hpp`` → UTF-8 → wide API    |
| Panel → subprocess (PySide6)       | ✅     | QProcess uses CreateProcessW             |

## Caveats

  * **``cubit -batch <japanese.jou>`` from command-line**: fails.
    Cubit parses argv via narrow cp932 on Japanese Windows, the
    .jou path gets mojibaked before Cubit can open it.  This is a
    Cubit core limitation — we do not patch it.
  * **Workaround for CI scripts**: ASCII-path wrapper.jou containing
    ``playback "<utf8 japanese path>"``.  The ``playback`` command
    uses wide API internally and reads the Japanese path correctly.
  * SSH-based pwsh heredoc with Japanese path **arguments**: cp932
    mangling.  Use heredoc script bodies (not argv) — scp the script
    with Japanese chars inside, then ``ssh host 'pwsh -File script'``.

## End-user workflow is safe

Kubota's actual workflow:
  1. Cubit GUI → File → Open .jou  (Qt wide API → works)
  2. Cubit GUI → Solve → Radia-NGSolve → PEEC inductance panel
  3. Browse to .step  (QFileDialog wide API → works)
  4. Run → QProcess spawns calc_peec_inductance.py  (CreateProcessW
     → Japanese args preserved)
  5. ``calc_peec_inductance.py`` handles Japanese path internally,
     sibling-jou detection works on UTF-8 dir listing

All verified on 100号機 with L_coil = 426.245 nH on Japanese path.
"""


def get_peec_inductance_documentation(topic="all"):
    """Return PEEC-inductance panel documentation by topic.

    Topics:
        overview          - What it solves, when to use, perimeter placement
        centerline        - STEP centerline extraction (5 classification predicates)
        filament_dispatch - Filament placement paths (single-dispatch, no fallback)
        step_authoring    - Cubit / build123d recipes for auto-detect-friendly STEPs
        japanese_path     - Unicode / Japanese path support
    """
    topics = {
        "overview": PEEC_IND_OVERVIEW,
        "centerline": PEEC_IND_CENTERLINE,
        "centerline_extraction": PEEC_IND_CENTERLINE,
        "filament_dispatch": PEEC_IND_FILAMENT_DISPATCH,
        "filament_placement": PEEC_IND_FILAMENT_DISPATCH,
        "dispatch": PEEC_IND_FILAMENT_DISPATCH,
        "n_peri": PEEC_IND_FILAMENT_DISPATCH,
        "uv_map": PEEC_IND_FILAMENT_DISPATCH,
        "per_station": PEEC_IND_FILAMENT_DISPATCH,
        "step_authoring": PEEC_IND_STEP_AUTHORING,
        "authoring": PEEC_IND_STEP_AUTHORING,
        "cad_recipes": PEEC_IND_STEP_AUTHORING,
        "cubit_recipe": PEEC_IND_STEP_AUTHORING,
        "build123d_recipe": PEEC_IND_STEP_AUTHORING,
        "anti_patterns": PEEC_IND_STEP_AUTHORING,
        "japanese_path": PEEC_IND_JAPANESE_PATH,
        "unicode_path": PEEC_IND_JAPANESE_PATH,
        "japanese": PEEC_IND_JAPANESE_PATH,
    }
    unique_docs = [PEEC_IND_OVERVIEW, PEEC_IND_CENTERLINE,
                   PEEC_IND_FILAMENT_DISPATCH, PEEC_IND_STEP_AUTHORING,
                   PEEC_IND_JAPANESE_PATH]
    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(unique_docs)
    if topic in topics:
        return topics[topic]
    return (
        f"Unknown topic: '{topic}'. "
        f"Available: all, overview, centerline, filament_dispatch, japanese_path."
    )
