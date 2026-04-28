"""PEEC-inductance knowledge (mcp-server-radia-ngsolve, public tool).

Panel mode: PEEC inductance (coil only, STEP).  Computes L_coil / R_coil
from a coil solid (STEP or Cubit .jou) via filament-based Biot-Savart +
loop-bundle PEEC solve.  No workpiece, no BEM, no FEM mesh — this is
the lightest path in the IH panel family and the one a user reaches
for when they want just a quick "what is the inductance of this coil
at this frequency?" answer.

Filaments are placed on the cross-section **perimeter only** (thin-skin
regime, d/δ >= 3).  Use ``n_peri`` filaments around the arc-length
perimeter of each cross-section; no interior volume grid.

Source: src/radia/panels/calc_peec_inductance.py,
        src/radia/coil_from_cad.py (STEP path),
        src/radia/coil_from_jou.py (.jou path),
        src/radia/radia_ih.py -- IHWindow with Method = "PEEC inductance
        (coil only, STEP)" (merged 2026-04-26 from the previously
        standalone radia_peec_inductance.py wrapper; the wrapper added
        no behaviour beyond auto-fill, which now lives on IHWindow).
"""


PEEC_IND_OVERVIEW = """
# PEEC inductance (coil only) — overview

## What it solves

Given a coil solid (STEP or Cubit .jou) + conductivity + frequency,
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
  * ``.jou``          — Cubit journal with explicit ``move Surface N x Y y Y z Z``
                        centerline (3turnCoil.jou pattern).  Fastest,
                        most accurate: parser reads the literal
                        coordinates, no CAD heuristic needed.
  * ``.step`` / ``.stp`` — geometry file; centerline auto-extracted from
                        solid topology (see `centerline_extraction`).
  * other             — raises ValueError.

Auto-preference: if the input is ``.step`` and a sibling ``.jou`` with the
same stem (case-insensitive) exists in the same directory, the panel
**automatically switches to the .jou path** and logs
``PEEC:found sibling .jou, preferring it: <name>``.  This matches the
Cubit workflow: the panel's ``ensure_jou_path()`` saves .jou before every
STEP export, so they naturally come as a pair.
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

**Validation** (Kubota's ``3turncoil.stp``, 10 MB, 382 loft volumes):
  * STEP-only:  L = 430.86 nH,  topology = 12.9 s
  * .jou path:  L = 426.25 nH,  topology =  3.5 s
  * Agreement: +1.1 %, well inside PEEC discretisation noise

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

**Validation** (all 4 coil topologies give correct L within a few %):

| Coil                           | Expected L | Got L    | Error   |
|--------------------------------|-----------|----------|---------|
| Circular 355° torus            | 88.55 nH  | 85.10 nH | -3.9 %  |
| Rect 6×4 mm 355° torus         | ~88 nH    | 88.15 nH | reference|
| 3-turn loft (sibling .jou)     | 426.25 nH | 426.25 nH| 0 %    |
| 3-turn loft (STEP only)        | 426.25 nH | 430.86 nH| +1.1 % |

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


PEEC_IND_JOU_PARSER = """
# PEEC inductance — .jou explicit centerline parser

``src/radia/coil_from_jou.py`` parses Cubit journal files that
explicitly define the centerline via ``move Surface`` commands.

## Canonical pattern

```
create surface circle radius 3.15 zplane           # cross-section
rotate Surface 1 angle ... direction ... include_merged
move Surface 1 x 129.394895 y -12.500000 z 10.000000 include_merged
create surface circle radius 3.15 zplane
rotate Surface 2 angle ... direction ... include_merged
move Surface 2 x 129.428679 y -12.500000 z 10.000050 include_merged
...
create volume loft surface 1 2
create volume loft surface 2 3
...
```

N ``move Surface`` commands → N centerline points (one per profile
station).  The first ``create surface circle radius R`` gives the
cross-section radius.

## Regex

```python
RE_CIRCLE = r'^\\s*create\\s+surface\\s+circle\\s+radius\\s+([\\d.eE+-]+)'
RE_MOVE   = r'^\\s*move\\s+Surface\\s+\\d+\\s+x\\s+([-\\d.eE+]+)'
            r'\\s+y\\s+([-\\d.eE+]+)\\s+z\\s+([-\\d.eE+]+)'
```

## Why prefer .jou over STEP

| Metric                              | STEP (cross-section method) | .jou (explicit) |
|-------------------------------------|-----------------------------|-----------------|
| Topology time (Kubota 3turncoil)    | 12.9 s                      | 3.5 s           |
| Cross-section radius accuracy       | face.area → sqrt(A/π)       | literal         |
| Centerline ordering                 | NN chain heuristic          | file order      |
| Works on non-loft coils (helix, swept) | ✅                      | only if .jou has move pattern |
| Robustness to CAD merge quirks      | needs dedupe tolerance      | irrelevant      |

.jou is the reference ground truth when available; the STEP path is
the fallback when only a STEP is in hand.

## Limitations

  * Parser reads only ``move Surface`` positions.  Rotation +
    non-circular profiles are NOT reconstructed — the cross-section
    is assumed to match ``RE_CIRCLE`` radius (circle).
  * Units: METRES (``cad_units_per_meter=1.0`` default since 4.12.0).
    Cubit set-up: ``set unit-system mks`` BEFORE generating the .jou.
    For legacy mm .jou files, pass ``cad_units_per_meter=1000.0``
    explicitly via the Python API.
  * Open path only (no closure).  The PEEC solver adds
    port_plus / port_minus at first / last centerline points.
"""


PEEC_IND_SIBLING_JOU = """
# PEEC inductance — sibling-jou auto-preference

When the user passes a ``.step`` to ``calc_peec_inductance.py`` and
a sibling ``.jou`` with matching stem (case-insensitive) is in the
same directory, the calc script **automatically switches to the
.jou path**:

```
PEEC:found sibling .jou, preferring it: 3turnCoil.jou
PEEC:JOU -> explicit centerline (n_peri=16)
PEEC:L_coil=426.245 nH, R_coil=0.3946 mOhm
```

## Why this is safe

  1. Cubit workflow: the IH panel's ``ensure_jou_path()`` saves the
     journal file before every ``radia_export / export step``.
     So the .jou is the source, the .step is a derived artifact.
     Preferring the .jou = preferring the source.
  2. Case-insensitive match: ``3turncoil.stp`` + ``3turnCoil.jou``
     (Kubota's pattern, capital-case difference from Windows/
     OS-X file creation) matches correctly.
  3. No silent fallback between differently-named files — matching
     requires the SAME STEM.  Unrelated .jou in the same dir is
     ignored.

## When the user wants STEP-only behaviour

  * Remove the sibling .jou, OR
  * Put the .step in a dedicated directory without a .jou, OR
  * Explicitly name the STEP with a stem no .jou matches

For testing STEP-only centerline extraction (the multi-turn loft
cross-section method), copy the .stp alone to ``C:\\temp\\stp_only\\``.

## Diagnostics in panel log

The panel's subprocess stderr includes ``PEEC:found sibling .jou,
preferring it: <name>`` when auto-switch happens.  This is logged
to ``C:\\radia_panel_log.txt`` so users can confirm which path was
actually taken.
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
| sibling-jou detection              | ✅     | ``os.listdir`` Unicode normalisation     |
| .jou parser (``open(utf-8)``)      | ✅     | utf-8 native                             |
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
        overview        - What it solves, when to use, perimeter placement
        centerline      - STEP centerline extraction (loft / open-spine)
        jou             - .jou explicit centerline parser
        sibling_jou     - Auto-prefer sibling .jou when co-located
        japanese_path   - Unicode / Japanese path support
    """
    topics = {
        "overview": PEEC_IND_OVERVIEW,
        "centerline": PEEC_IND_CENTERLINE,
        "centerline_extraction": PEEC_IND_CENTERLINE,
        "jou": PEEC_IND_JOU_PARSER,
        "jou_parser": PEEC_IND_JOU_PARSER,
        "sibling_jou": PEEC_IND_SIBLING_JOU,
        "sibling": PEEC_IND_SIBLING_JOU,
        "japanese_path": PEEC_IND_JAPANESE_PATH,
        "unicode_path": PEEC_IND_JAPANESE_PATH,
        "japanese": PEEC_IND_JAPANESE_PATH,
    }
    unique_docs = [PEEC_IND_OVERVIEW, PEEC_IND_CENTERLINE,
                   PEEC_IND_JOU_PARSER, PEEC_IND_SIBLING_JOU,
                   PEEC_IND_JAPANESE_PATH]
    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(unique_docs)
    if topic in topics:
        return topics[topic]
    return (
        f"Unknown topic: '{topic}'. "
        f"Available: all, overview, centerline, jou, sibling_jou, "
        f"japanese_path."
    )
