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

Source: src/radia/panels/calc_peec_inductance.py,
        src/radia/coil_from_cad.py (STEP path -- the only input path
        since 4.13.0; .jou explicit-centerline input was retired
        per CLAUDE.md No-Fallbacks),
        src/radia/radia_ih.py -- IHWindow with Method = "PEEC inductance
        (coil only, STEP)" (merged 2026-04-26 from the previously
        standalone radia_peec_inductance.py wrapper; the wrapper added
        no behaviour beyond auto-fill, which now lives on IHWindow).
"""


PEEC_IND_OVERVIEW = """
# PEEC inductance (coil only) — overview

## What it solves

Given a coil solid STEP file + conductivity + frequency,
returns:
  * L_coil_nH  — inductance in nanohenry
  * R_coil_mOhm — AC resistance including skin effect via SIBC
  * Z_coil = R + j·ω·L (complex impedance)

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

Each filament carries the per-filament SIBC impedance Z_s = (1+j)/(σδ)
for circular wire; Dowell for rectangular; ESIM for nonlinear steel.

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

``extract_centerline_from_step`` in ``src/radia/coil_from_cad.py``
dispatches on solid topology:

## 1. Loft of profiles (multi-turn coil, tight pancake spiral)

**Trigger**: solid has >= 5 planar faces of consistent area
(``_collect_loft_cross_sections`` filters planar faces, checks
median-area ± 50% band).

**Algorithm**:
  1. Enumerate planar cross-section faces (``GeomType.PLANE``)
  2. Dedupe near-duplicate centroids (shared end-cap between
     adjacent loft volumes appears twice, once per owning volume —
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
calls — the cross-section faces ARE the sections.  Old code called
``section(solid, Plane)`` once per midpoint (100 calls = 100 OCC
boolean cuts = minutes on 10 MB STEPs).

## 2. Single-loop revolution-sweep coil (gapped torus, rect-profile bend)

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

## 3. Generic swept coil fallback (non-revolution)

**Trigger**: falls through when both (1) and (2) do not apply (the
solid has no revolution-axis structure — e.g. irregular bend,
free-form helix extrusion along a non-axial path).

**Algorithm**: find the LONGEST open boundary edge
(``is_closed == False``) as the spine, sample it at
``n_segments + 1`` arc-length-equidistant points, section at
midpoints to extract cross-section centroids and areas.

This path rarely fires for panel-generated coils (Cubit sweep /
loft always produce revolution-type surfaces).  Reserved for
externally-authored STEP files with exotic topology.

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

Since v4.15.0 (2026-04-28), ``filaments_from_step`` dispatches filament
placement through a **3-tier fallback chain** to handle increasingly
general coil topologies.  Order in
``src/radia/coil_from_cad.py``:

## Tier 1 — UV-map sampling (NEW v4.15.0, most general)

**Trigger** (`_find_lateral_surface`): solid has ONE dominant lateral
surface (BSPLINE / CYLINDER / TORUS / REVOLUTION / EXTRUSION) by area:
  * ``largest_area / total_area >= 0.8`` (single-piece guard)
  * ``len(candidates) <= 3`` (defensive bound; multi-piece trips this)

**Algorithm** (`_sample_lateral_surface_uv` →
`_filaments_from_lateral_surface_uv`):
  1. Auto-detect closed UV axis (perimeter) vs open (spine) via
     ``Geom_Surface.IsUClosed`` / ``IsVClosed``
  2. Sample at:
       * spine: ``0..n_stations-1`` endpoints (inclusive)
       * perimeter: cell-centred ``(k+0.5)/n_peri``
  3. Build ``n_peri`` filaments from the (n_stations × n_peri × 3) grid
  4. Per-segment cell area via 3D-shoelace on the perimeter polygon at
     each station — **variable cross-section is automatic**, no
     equivalent-circle assumption.

Falls through to Tier 2 when the lateral surface is fragmented.

## Tier 2 — Per-station face sampling (v4.14.0)

**Trigger**: solid is a NON-united multi-loft with planar end-cap
faces surviving in the STEP (loft volumes preserved, no boolean
``unite``).

**Algorithm** (`_filaments_from_per_station_faces` →
`_sample_face_perimeter_in_pt_frame`):
  1. Enumerate planar cross-section faces along the spine
  2. Sample EACH cross-section face's outer wire at ``n_peri``
     **arc-length-equispaced** points
  3. Project each sample to the parallel-transport ``(u_hat, v_hat)``
     frame at that station
  4. Filament k connects sample k of station i to sample k of
     station i+1

Handles RECTANGULAR / POLYGON / arbitrary cross-section + any
spine-direction shape variation (tapered conductor, varying-width
pancake, twisted bus bar).

Falls through to Tier 3 when the loft is gapped or united.

## Tier 3 — Constant equivalent-circle (legacy, ≤ v4.13.0)

**Trigger**: catch-all for non-loft solids (gapped torus, open-spine
helix, fully united multi-turn pancake).

**Algorithm**: ``filaments_from_polyline`` with constant equivalent-
circle radius derived from mean cross-section area.

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
| Cross-section **discontinuity** along the spine (circle→rect, sudden step in area) | Tier 2 connects sample-k-station-i to sample-k-station-i+1 across the discontinuity → non-physical filament path; Tier 3 averages the area → loses detail | **Phase C-heavy** (deferred since v4.17.0; see ``project_peec_phase_c_heavy_deferred.md``) |
| **United** multi-loft with non-circular cross-section (rect / polygon / shape-transition) | Tier 1 UV-map trips because the lateral surface is fragmented across the ``unite()`` boolean; Tier 3 falls back to equivalent-circle (Kubota's 3turncoil case: 421.837 nH vs 426.25 nH golden = -1.0 %) | Phase C-heavy (BSPLINE-adjacency graph traversal + Hamiltonian path) |
| **T-junction / Y-branch** bus bar                  | single-spine polyline assumption breaks; ``filaments_from_polyline`` cannot represent topology | Multi-spine PEEC re-architecture (no concrete plan) |
| **>2-port** coil (multi-tap)                       | ``add_port(n1, n2)`` is single-pair; cannot express 3+ port impedance matrix | Multi-port PEEC topology extension |
| **Closed loop** with no external terminals         | port assignment fails; cannot define the source current | Out of scope; PEEC inherently needs at least 1 driving port |
| **Litz wire** (N strands per cross-section)        | perimeter is sampled as a single outer wire; bundle inter-strand mutual L is not represented | Bundle-PEEC (multiple parallel filament sets per spine station) |
| **Magnetic conductor** (μ_r > 1, e.g. iron busbar) | SIBC ``Z_s = (1+j) / (σ δ)`` is for non-magnetic; need ``Z_s = (1+j) sqrt(ω μ / (2 σ))`` with the correct μ | One-line patch in ``rad_peec_surface_impedance.cpp`` after deciding the σ + μ_r input UX |
| **Thick-skin** ``d/δ < 3``                         | perimeter-only filaments miss the interior current distribution | Out of scope for this panel; switch to FEM A-V (volumetric coil mesh) |
| **Self-intersecting** STEP                         | OCC topology load undefined behaviour | Reject at parse time |
| **Non-manifold** STEP                              | ``build123d.import_step`` partial-parses or silently drops shells | Add a manifold-check before dispatch |

**Action when the user asks for any of the above**: do NOT silently
fall back to Tier 3. Either (a) refuse + explain the limitation
(if the inaccuracy is unacceptable), or (b) accept Tier 3 with a
WARN log line that names the specific limitation (e.g. ``WARN:
united multi-loft falling back to equivalent-circle, expect -1 to
-5 % on L``).
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
        centerline        - STEP centerline extraction (loft / open-spine)
        filament_dispatch - 3-tier filament placement (UV-map / per-station / equiv-circle)
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
        "japanese_path": PEEC_IND_JAPANESE_PATH,
        "unicode_path": PEEC_IND_JAPANESE_PATH,
        "japanese": PEEC_IND_JAPANESE_PATH,
    }
    unique_docs = [PEEC_IND_OVERVIEW, PEEC_IND_CENTERLINE,
                   PEEC_IND_FILAMENT_DISPATCH, PEEC_IND_JAPANESE_PATH]
    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(unique_docs)
    if topic in topics:
        return topics[topic]
    return (
        f"Unknown topic: '{topic}'. "
        f"Available: all, overview, centerline, filament_dispatch, japanese_path."
    )
