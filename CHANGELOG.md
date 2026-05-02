# Changelog

All notable changes to the `radia` package.  Format: each release lists
**what shipped** + **why** in compact form.  Packaged wheels on PyPI.

## 4.25.0 — BEM-A coil + unified inductance CLI + 6-method IH panel

Released 2026-05-02.  Surface-current coil solver lands as a peer of
PEEC; three scoped CLIs collapse into one.

### BEM-A coil solver (Phase C.1-C.5 shipped)

`radia.bem.efie_rwg` (~1 kLOC) — Weggler stabilized EFIE saddle on
HDivSurface RWG (= Lucy decomposition implicit), in-tree Sauter-Schwab
Duffy 4-D Galerkin assembler, source/sink port-driven self-inductance,
DC + AC SIBC closure for R extraction.  Tested against ngsolve.bem
oracle to 0.025-0.004 % on standard fixtures (5 test files, 13 cases).

Cross-method validation against PEEC golden on
`rect_torus_lofted_united.step` (8x6 mm rect, 50 kHz Cu):
- BEM-A converged 153.5 nH vs PEEC asymptote ~149 nH (n_peri → ∞)
- ~3.5 % residual gap is real modeling difference (BEM-A surface RWG
  resolves rect-corner current crowding; PEEC perimeter filaments do
  not).  Both methods are SIBC; both converge correctly within their
  discretisation classes.

### Unified `calc_inductance.py` (Refactor A1+B+C)

Single CLI replaces three scoped predecessors:

| OLD | NEW |
|-----|-----|
| `calc_peec_inductance.py` | `calc_inductance.py --coil-solver peec` |
| `calc_peec_bem.py` | `calc_inductance.py --coil-solver peec --vol …` |
| `calc_coil_bem_a_workpiece.py` | `calc_inductance.py --coil-solver bem-a --vol …` |

Dispatch on `--coil-solver {peec, bem-a}` × `--vol {present, absent}`
gives 4 modes (vacuum or weak-coupled, PEEC or BEM-A coil).  All
share the workpiece scalar BEM-SIBC block.

**Coupling terminology corrected**: previous "1-way forward" was
misleading because Telegen φ·(n·B) ΔL IS computed, capturing port-level
back-reaction even though coil J is fixed.  Renamed to **weak coupling**
throughout (CLI flag, panel labels, docstrings).  Strict one-way
(no ΔL) and strong coupling (coil J recomputed iteratively, FEM A-V)
remain distinct.

### IH panel: 6-method dropdown

`radia_ih.py` adds two BEM-A counterparts to the existing PEEC modes:

| Method | Coil | Workpiece | Notes |
|--------|------|-----------|-------|
| PEEC inductance (vacuum) | PEEC filament | — | existing |
| **BEM-A inductance (vacuum)** | BEM-A surface RWG | — | new |
| PEEC + BEM weak coupling | PEEC | scalar BEM-SIBC | renamed (was "1-way") |
| **BEM-A + BEM weak coupling** | BEM-A | scalar BEM-SIBC | new |
| PEEC coil + FEM Kelvin | PEEC | volumetric FEM-SIBC | unchanged |
| FEM A-V full | volumetric coil | volumetric FEM-SIBC | unchanged |

`coil_maxh` widget added (visible only for BEM-A modes).  19/19
`tests/panels/test_ih_panel_qt.py` PASS; 36/36 combined panel +
inductance + coil_topology golden tests PASS.

### PEEC unified OPEN/CLOSED cap-aware spine

New `radia.coil_topology` module (~280 LOC) — single source of truth
for OPEN vs CLOSED coil classification (cap detection, spine arc
parameters).  Consumed by:
- Path 2c (`_filaments_from_section_planes`) — `rect_torus_lofted_united`
  PEEC L 138.16 → 145.30 nH (filament was clipping ~14° off the 355°
  arc due to a `linspace(0, 2π)` fallback).
- Path 2b (`_filaments_from_circle_edges_per_station`) — single-turn
  coils now use topology spine; multi-turn helix guard
  (z_extent > 2 · median_r) keeps 3turncoil 422 nH on legacy NN-chain.
- Path 3 fallback — new `_centerline_from_topology_spine` helper
  inserted between `_centerline_from_torus_sweep` and
  `_centerline_from_open_spine`: closed-torus filament sweep
  178° → 358° (was tracing one half-arc seam edge).

Golden update: `peec_inductance_rect_united_50kHz_Cu.json`
138.16 → 145.30 nH, `regression_guard.min_L_nH` 120 → 130.

### Files

| Action | Path |
|--------|------|
| New | `src/radia/bem/efie_rwg.py` |
| New | `src/radia/coil_topology.py` |
| New | `src/radia/panels/calc_inductance.py` |
| Removed | `src/radia/panels/calc_peec_inductance.py` |
| Removed | `src/radia/panels/calc_peec_bem.py` |
| Removed | `src/radia/panels/calc_coil_bem_a_workpiece.py` |
| New | `tests/panels/test_inductance_golden.py` (6 tests) |
| New | `tests/test_coil_topology.py` (11 tests) |
| New | `tests/bem/test_coil_bem_a_*.py` (5 files, 13 tests) |
| Removed | `tests/panels/test_peec_inductance_golden.py` |
| Removed | `tests/panels/test_peec_bem_golden.py` |
| New | `examples/coil_bem_a/` (build123d demo + GMSH viz) |
| Updated | `src/radia/radia_ih.py` (6-method dropdown + dispatch) |
| Updated | `src/radia/panels/sync_registry.py` + `panel_registry.json` |

## 4.10.0 — Panel UX overhaul: font baseline, base-class unification, CoilBuilder wizards

Released 2026-04-26.  All-Python release; no C++ rebuild needed.
Focus is the panel runtime that EM, IH, and PCB share.

### Font baseline (readable on 2K)

The Qt OS default of 9pt Segoe UI was unreadable on the lab's 2K
displays (per `feedback_panel_vertical_space.md`).  After three turns
of bracketing with the user (9pt unreadable -> 11pt stingy -> 13pt
overshoot -> **12pt chosen**), `apply_panel_base_font(app)` is now
called from `run_app()` and from the test fixture; it inherits the OS
default family (Segoe UI on Windows, system sans-serif on Linux/macOS)
and bumps only the point size.

Hardcoded smaller fonts removed:
  - output text area: 9pt -> 12pt (Consolas)
  - status label (IH): "font-size: 11px" deleted -- de-emphasis is
    now colour-only (#888), font size inherits the baseline.

Regression guard: new `panel_qa.check_font_size_min` walks every
visible widget and fails if any renders below 10pt.  Width / height
thresholds rebudgeted (1400 RED / 1350 RED) for the larger glyphs.

### Panel base-class unification

EM and IH carried two identical 9-line copies of the section-header
helper.  One source of truth in `ModePanel`:

  - `_add_section(title, key=None)`           hoisted to base
  - `add_status_label(key=None)`              hoisted to base
  - `_method_combo` attribute convention      adopted by EM
                                              (was `_form_combo`)
  - `add_browse_action(key, label, callback)` new helper -- attach
                                              extra buttons to an
                                              existing browse row

PCB combo extraction `"0 (LU)".split()[0]` replaced with an explicit
`PCB_SOLVERS = {"LU": 0, "BiCGSTAB": 1, "HACApK": 2}` map.

### PEEC-inductance Window merged into IH

`radia_peec_inductance.py` (95-line wrapper) deleted.  The same
analysis is reached by selecting Method = "PEEC inductance (coil
only, STEP)" in the IH window.  IH's `__init__` auto-fills the STEP
field from the newest `*.step` / `*.stp` / `*.jou` in cwd when the
field is empty -- behaviour previously specific to the wrapper
window now applies to any IH launch in PEEC-inductance mode.

The Cubit launcher (`RadiaComp.cpp`) discovers panels by globbing
`radia_*.py` so the merged entry disappears from the menu without
any .ccl change.

### "New..." coil wizards (EM + IH)

Both panels expose a [New...] button next to the coil-source field:

  - **EM** -> writes a self-contained `.py` template (CoilBuilder
    racetrack with EDIT BLOCK at the top).  The user picks the save
    path, edits the 6 numbers (NI / START / WIDTH / HEIGHT /
    STRAIGHT / ARC_R), and points the panel's Coil script at it.
  - **IH** -> writes the same `.py` AND immediately runs it in-process
    to materialise a sibling `.step` (the format the IH PEEC modes
    consume).  User can later edit the `.py` and re-run
    `python coil.py` to refresh the `.step`.

Single COIL_TEMPLATE constant in `radia_gui_base.py` powers both.
Verified end-to-end: 84-segment closed racetrack, 2000 A, gap=1e-17 m,
.step = 125 KB OCC swept solid.

## 4.9.0 — EM canonical trio + 1/2 Kelvin Benchmark + MSC silent breakage fix

Released 2026-04-26.

### MSC silent breakage fix (production bug)

`calc_accel_msc.py` registered `add_material_args(include_custom=False)`,
so the EM panel's MSC mode (which sends `--material custom --mu-r <user>`
for the "mu_r (Linear)" material option) was rejected at argparse with
`error: argument --material: invalid choice: 'custom'`.  Clicking Run in
MSC mode produced no useful error in the GUI and no field result.  Fixed
by flipping to `include_custom=True` (matches the FEM Omega/A-Phi
convention in `calc_accel_magnet.py`).  Regression guard:
`tests/panels/test_em_msc_smoke.py` (6 sub-second static checks of the
MSC panel command vs argparse).

### Kelvin Benchmark — 1/2 sample joins 1/4

The EM panel's Kelvin Benchmark mode now ships **two** verification
samples (was 1/4 only):

| sample (frac) | reduction | error @ p=2 |
|---|---|---|
| `kelvin_benchmark_sphere_1_2.vol` | `{"y":"bn=0"}` | +1.07% |
| `kelvin_benchmark_sphere_1_4.vol` | `{"x":"bn=0","y":"bn=0"}` | +0.71% |

Build scripts refactored: parameterized core
`kelvin_benchmark_sphere_build.py` with `--frac 1_2|1_4|1_8` plus three
thin wrappers.  1/8 sphere build script ships for research only — the
.vol is **not** shipped because the magnetic-sphere-in-uniform-Hz BVP
fundamentally lacks 1/8 symmetry (z=0 mirror reverses Hzẑ source).
Realistic 1/8 EM workloads (C-yokes, dipoles) use the existing
`em/em_1-8_eighth.jou` C-yoke sample with the ELF `-x-y+z` convention.

Surprise (rho_min sweep): for compact-geometry Kelvin samples
(offset = 3R), the Periodic + sym BCs do most of the open-boundary
work — the `(R/r')^2` reluctivity is a small correction.  Capping
`Mu = mu_0` everywhere still gives 1/2 +0.34% / 1/4 -0.02% error.

### EM panel canonical trio (coil + .vol + BH)

Added two missing artifacts so the EM panel runs end-to-end on a
fresh checkout / pip install without requiring local Cubit
regeneration of every sample:

  - `src/radia/panels/samples/em_sample.vol` — half-z C-type dipole
    with auto-Kelvin (456k tets / 83k nodes).
  - `src/radia/panels/samples/em_sample_bh.txt` — 100-point steel
    BH curve (CEFC 2020 reference, 0..318 kA/m / 0..2.61 T).
    Full-precision source of the rounded built-in `STEEL_BH` table;
    functionally equivalent (relative diff < 2e-5).

Coil was already shipped as `em_sample_coil.py`.  Sub-second smoke
tests in `tests/panels/test_em_sample_artifacts.py` (8 tests) lock
ship-presence + parse-correctness for the trio.

### IH panel test refresh (4 stale test files repaired)

The 2026-04-19 IH panel restructure (4-method combo, separate
coil/wp material sections) left the test suite stale.  Rewrites:

  - `tests/panels/test_ih_panel_qt.py` — was 14/15 failing, now 17/17
    passing (method combo, per-method widget visibility, solver
    items, SIBC vs ESIM toggle, build_command roundtrip).
  - `tests/panels/test_panel_state_restore.py` — 6 tests rewritten
    from the obsolete 2-method labels to the canonical METHOD_*
    constants.  Round-trip widgets updated from retired
    `workpiece_mode` to surviving `impedance_model`.
  - `tests/panels/test_calc_inductance.py` — was failing at
    collection (ModuleNotFoundError on retired `bem_inductance`).
    `bem_inductance.py` now lives at
    `examples/induction_heating/bem_reference/`; tests `importorskip`
    cleanly when that dir is missing.
  - `tests/panels/test_ih_solvers.py` — same fix pattern as
    test_calc_inductance for the BEM precondition guard.

Also: `tests/panels/golden/em_eighth_mu1000.json` topology lock
refreshed for the post-2026-04-25 deterministic-anchor mesh
(ne 56289 → 56369, ndof 11695 → 11708; physics unchanged).
Slow physics regressions (`test_peec_bem_matches_2d_ref_coarse` 174s,
`test_fem_coilmesh_matches_2d_ref_gapped` 510s) both pass on the
shipped IH samples.

### Sample matrix documentation

Added `src/radia/panels/samples/README.md` indexing the canonical
samples per panel/method, identifying stale/research-only artifacts
that don't ship in CI-built wheels (gitignored .vol files).

### New panel logo

`src/radia/resources/2026_04_26.png` replaces the deleted
`radia_icon.png`.  `_icon_path()` now falls back to the most-recent
.ico/.png in `resources/` so future logo refreshes drop in without
a code change.

### radia-mcp 0.33.2

Updates `radia_ngsolve.kelvin_knowledge.benchmark_panel` topic to
document the 1/2 + 1/4 sphere samples + the rho_min sweep insight.
No API change.

## 4.7.0 — PEEC-inductance full geometry coverage + Japanese path + daemon speedup

Released 2026-04-22.

### PEEC-inductance — 4 coil topologies all STEP-only

The `radia_peec_inductance` panel mode ("PEEC inductance (coil only,
STEP)") can now extract a correct centerline from every coil topology
the lab uses, without requiring the user to supply an explicit `.jou`:

| Coil topology | Before 4.7.0 | After 4.7.0 |
|---|---|---|
| Circular gapped torus (single loop) | walker spine broken, L off 100× | ✅ TORUS analytical sweep, L = 85.10 nH (analytical 88.55, -3.9 %) |
| Rect 6 × 4 mm gapped torus | walker failed | ✅ CYLINDER analytical sweep, L = 88.15 nH |
| Multi-turn pancake loft (e.g. Kubota `3turncoil.stp`) | walker hung > 5 min, native crash | ✅ cross-section centroid NN chain with tangent continuity, L = 430.86 nH (matches `.jou` 426.25 nH to +1.1 %) in 12.9 s |
| `.jou` explicit centerline | | ✅ unchanged, fastest path |

New support code in `src/radia/coil_from_cad.py`:

- `_centerline_from_revolution_sweep()` — OCP `BRepAdaptor_Surface`
  extraction of TORUS / CYLINDER / CONE / REVOLUTION parameters,
  union of U-intervals = sweep angle, centroid of smallest PLANE
  face = spine radius.
- `_centerline_from_cross_sections()` — `GeomType.PLANE` face
  enumeration + dedupe near-duplicates + NN chain with tangent-
  continuity bias.
- `_auto_detect_cad_units()` — bbox-diagonal heuristic; handles
  Cubit STEPs that are sometimes mm-unit (values 10 – 1000) and
  sometimes m-unit (values 0.01 – 1) despite declaring
  `CONVERSION_BASED_UNIT('MILLIMETRE')`.

### Sibling `.jou` auto-preference (with format guard)

When `calc_peec_inductance.py --peec-step foo.step` is called and a
`foo.jou` sibling exists in the same directory, check whether the
`.jou` contains explicit `move Surface N x Y y Y z Z` lines before
switching to it.  A generic Cubit .jou (rect-sweep, brick setup)
no longer crashes the parser.

### Golden regression tests for PEEC-inductance

`tests/panels/test_peec_inductance_golden.py` — 3 tests locking in
the 3 canonical paths (TORUS analytical / cross-section centroid /
sibling .jou), each with a hard band + golden tolerance.

### IH panel workflows

PEEC+BEM (`calc_peec_bem.py`) and FEM A-V (`calc_fem_coilmesh.py`)
unchanged; verified on LAB + 100号機 that no regression.

### Packaging / deployment / documentation

- Cubit MCP daemon speedup (via `radia-mcp` 0.32.0 — see its CHANGELOG):
  cold start 30 – 60 s → 3 s.  VSCode restart 6 s → 0.01 s attach.
- `CLAUDE.md`: new Panel Samples Quality Policy + MCP Knowledge
  Placement Policy sections.
- `tests/panels/panel_qa.py`: 7 checks × every panel mode (deploy
  skill enforces).
- Deploy skill: per-mode Panel Mode Matrix + golden-range numeric
  assertions replace `status: "ok"` as the definition of "works".
- docs/research/IGA_CLN_DUAL_REDUCTION.md (LAB-local research plan).

### Known issues / non-goals

- Daemon auto-start across user logoff: Phase 2 deferred — Cubit
  GUI is desktop-bound so logoff-survival has no practical value.
- MCP server Python module hot-reload: still requires VSCode
  restart (OS-level limit, not addressable by daemon decoupling).

## 4.6.x — Panel QA + PEEC-inductance initial drop

(See git log before v4.7.0.  Adds the PEEC-inductance panel as a
new mode of the Radia-NGSolve launcher, plus Panel QA automated
layout checks.)
