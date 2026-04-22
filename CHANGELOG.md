# Changelog

All notable changes to the `radia` package.  Format: each release lists
**what shipped** + **why** in compact form.  Packaged wheels on PyPI.

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
