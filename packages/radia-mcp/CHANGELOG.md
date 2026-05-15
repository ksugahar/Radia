# Changelog

All notable changes to `radia-mcp`. Format: each release lists **what
shipped** + **why** in compact form. Older releases (≤ 0.4) are
omitted; the 0.5 → 0.6 jump is when the standalone `radia-mcp` wheel
crystallized as its own package.

## 0.48.4 — peec_inductance knowledge updated for v4.48.1 STEP-only centerline

Released 2026-05-16.  Pairs with radia 4.48.1 which replaced the
spine-extractor try/except cascades in `coil_from_cad.py` with
classification-based single dispatch + removed the `path_points_m`
parameter ("STEP-Only Centerline: Auto-Detect or Fail" policy in
CLAUDE.md).  The 0.48.3 bump was a version-coordination only; this
0.48.4 ships the actual knowledge updates:

### What shipped

- **`PEEC_IND_FILAMENT_DISPATCH`** (topic `filament_dispatch`): rewritten
  from "3-tier fallback chain" language to **classification-based
  single dispatch**.  Documents Path 1 (UV-map; predicate now
  includes the UV-closure check so downstream sampling MUST succeed
  -- no try/except in Path 1), Path 2 (per-station faces), Path 2b
  (CIRCLE-edge stations), Path 2c (section-planes), Path 3
  (equivalent-circle catch-all with the new fail-fast sanity check).
- **`PEEC_IND_CENTERLINE`** (topic `centerline`): expanded from 3
  paths to **5 classification predicates** (Loft / Circle-edge /
  Revolution-sweep / OPEN longest-edge / CLOSED full-revolution).
  Documents the CLOSED-only guard in `_centerline_from_topology_spine`
  and the keiko `1turn_coil_loft_outsideline.step` lesson (OPEN
  geometries with leads must route to Predicate 4, not 5).
- **New `PEEC_IND_STEP_AUTHORING`** (topic `step_authoring` +
  aliases `cubit_recipe`, `build123d_recipe`, `anti_patterns`):
  concrete recipes for authoring auto-detect-friendly STEPs.
  Quick-decision table mapping Cubit/build123d operations to
  predicate hits, full Cubit `.jou` recipes for gapped torus and
  multi-turn pancake, build123d `sweep()` recipe for curved
  spine + circular profile, anti-patterns (lateral split into 2
  halves, pairwise loft chain, hardcoded IDs, non-manifold,
  self-intersecting), and a 10-line build123d probe script for
  verifying a STEP is auto-detect-friendly BEFORE running the panel.

### Why

radia-mcp 0.48.3 (released 2026-05-15 alongside radia 4.48.1) only
bumped versions for release-triple coordination -- the knowledge
documents still described the obsolete try/except cascade.  Users
asking the `peec_inductance(topic=...)` MCP tool got stale guidance.
0.48.4 reconciles the knowledge layer with the v4.48.1 dispatcher.

## 0.40.0 — 3D CLN (Tanimoto-Kameari) knowledge module

New `radia_ngsolve.knowledge.cln_3d` module captures Tanimoto's 3D
Cauer Ladder Network (CLN) methods from W:/00_CAE/NGSolve/谷本/
master's thesis + production code (~25 notebooks). Covers:

  - **A-T**, **T-Ω**, **A-Φ** formulations (mathematical foundation,
    iteration pseudocode, common boilerplate)
  - **Constraint variants**: penalty stabilization, explicit Coulomb gauge
  - **Solver variants**: SparseSolvPy ICCG, accICCG, NGSolve CG, direct
  - **Validation**: cylindrical TM-mode analytical R/L, Schmidt drift
    diagnostic, bonus_intorder=8 critical setting
  - **Open research note**: Kameari + Kelvin combination remains
    unsolved (3D HCurl A-formulation gives ~25× discrepancy with
    mpmath BEM Foster target due to A_ext gauge unboundedness)

Five canonical notebooks embedded as `cln_notebooks/*.py` resources:
  - `CLN_AT.py` (primary 修論 reference, 7.4 KB)
  - `CLN_T_Omega.py` (T-Ω formulation, 7.6 KB)
  - `CLN_APhi.py` (A-Φ formulation, 8.6 KB)
  - `CLN_2D.py` (2D scalar reference, 2.7 KB)
  - `A_ICCG_production.py` (latest 2024-09-17 production, 6.9 KB)

New MCP tools:
  - `cln_3d(topic="all"|"overview"|"notebooks"|"formulas")`:
    structured documentation
  - `cln_3d_notebook(name="list"|"AT"|"T_Omega"|"APhi"|"2D"|"production")`:
    raw Python code retrieval


## 0.33.5 — Sync with radia 4.10.0 (PEEC-inductance Window merged into IH)

`radia_ngsolve.peec_inductance_knowledge` Source list updated: the
standalone `radia_peec_inductance.py` wrapper was merged into IHWindow
in radia 4.10.0; the analysis is now reached via Method dropdown.
Knowledge text re-points new users at the IHWindow path so MCP
suggestions stay accurate.

No behavioural changes to any MCP tool.

## 0.33.4 — Kelvin knowledge maturity pass (republished)

Same content as 0.33.3 but with a shortened pyproject `description`
field (PyPI's 512-char `summary` limit rejected 0.33.3's metadata
upload at 596 chars, so the wheel never made it to PyPI).  No
behavioural / knowledge changes vs the unreleased 0.33.3; see below
for the actual changes.

## 0.33.3 — Kelvin knowledge maturity pass

Knowledge-only release across 3 subpackages, capturing the
2026-04-26 1/2 + 1/4 Kelvin Benchmark debug session and clarifying
why the 1/8 case has two completely different answers depending on
which panel mode is asking.

- **`radia_ngsolve.kelvin_transformation` (`benchmark_panel` topic)**:
  - Why 1/8 is unsupported for the magnetic-sphere-in-uniform-Hz BVP
    (the source `H0 z_hat` reverses sign under z=0 mirror -- a
    physical limitation, not a Cubit/NGSolve bug).
  - **rho_min sweep diagnostic**: setting rho_min = R collapses
    Mu = mu_0 *(R/rho')^2 to uniform mu_0; if the answer becomes
    correct, the bug is in the Mu coefficient; if still wrong,
    the bug is in BCs / Periodic / mesh.  One solve isolates the
    layer.
  - Surprise: for compact geometry (Kelvin offset = 3*R), even
    Mu = mu_0 in the Kelvin region gives 1/2 +0.34% / 1/4 -0.02% --
    Periodic + sym BCs do most of the open-boundary work.
  - **Cubit-meshed Kelvin needs `-specialcf.normal`** in the
    reduced-Omega Neumann correction term (Cubit assigns surface
    normals with opposite sign to NGSolve's WorkPlane OCC; sign-
    flip A/B test takes 30 seconds and catches it).

- **`cubit` (new `kelvin_reduction_traps` topic)**:
  - Trap 1: `subtract A from B keep` is a silent no-op in Cubit
    2025.3 -- workaround is to drop `keep` and re-create A as a
    fresh primitive.
  - Trap 2: 1/8 octant copy-mesh anchor curve picking is non-
    deterministic (3 equal-length quarter-arcs); fix is
    `min(curves, key=(centroid_z, y, x))` -- 143/143 pairs at
    machine precision.
  - Trap 3: surface normal sign convention differs between Cubit
    and OCC (cross-ref to `radia_ngsolve.kelvin_transformation`).

- **`electromagnet` (new `symmetry_reductions` topic)**:
  - Two distinct Kelvin panel paths -- "Kelvin Benchmark" sphere
    (1/2, 1/4 only) vs "EM panel FEM/MSC" C-yoke (1/1, 1/2, 1/4,
    1/8).  Don't conflate.
  - C-yoke 1/8 sample paths and ELF CEFC 2020 convention
    `ht=0_x, ht=0_y, bn=0_z`.
  - "Don't add a 1/8 sphere benchmark" -- multi-hour debug trail
    capture so the next session doesn't re-investigate.

## 0.32.0 — PEEC-inductance public topic + Cubit daemon speedup

- **`peec_inductance` tool** in `mcp-server-radia-ngsolve`: 5 sub-topics
  (overview / centerline / jou / sibling_jou / japanese_path) promoted
  from LAB-private `mcp-server-ih` after the feature stabilised.
- **Cubit daemon license warmup**: `cubit_license_warmup.py` mirrors
  `coreform_cubit.ps1` renewals cache logic (3-day cache + 7-day
  expiry).  Cold daemon start 30 – 60 s → 3 s.
- **Cubit daemon Phase 1 attach**: per-user stable drop-dir
  (`%LOCALAPPDATA%\radia-mcp\cubit-session\`) + `pid.lock` discovery.
  VSCode restart → new MCP server attaches to living Cubit in
  **0.01 s** instead of re-spawning (6 s cold).
- `open_in_cubit`: same license warmup applied so one-shot GUI
  launches from VSCode also get the speedup.
- `cubit_session_status` reports `mode = owned | attached`.
- New MCP knowledge placement policy in `CLAUDE.md`: stable /
  general → public `radia-mcp` (PyPI), research-stage / lab-only →
  `S:\mcp-server\mcp-server-ih`.

## 0.23.x — YouTube + training pack + GitHub `.jou` search

- **0.23.1** (planned, docs-only): full README rewrite with badges /
  multi-server table / quickstart / lab stance / acknowledgments;
  CHANGELOG.md + CONTRIBUTING.md added. (You're reading it.)
- **0.23.0**: YouTube tutorial transcript scraping for
  `cubit_youtube` / `build123d_youtube` / `gmsh_youtube`
  sub-sources (`youtube-transcript-api` extra). Coreform training
  `examples_only.zip` (24 MB / 30 .jou) auto-folded into
  `cubit_local`. PAT-gated `gmsh_post_jou_github` GitHub-wide `.jou`
  code search. New optional extra `radia-mcp[youtube]`.

## 0.22.x — Universal CAD-MCP mesh backend

- **0.22.4**: lab stance refinement — FreeCAD marked `friendly` /
  `compat — Sugawara Lab respects the FreeCAD community`; build123d
  + Cubit explicitly tagged `主力 (push)` in `lab_policy` topic.
- **0.22.3**: Sugawara Lab primary-pair stance reflected in
  `lab_policy` KB topic + `list_cad_mcp_interop` payload (`lab`,
  `primary_pair` fields) + memory.
- **0.22.2**: build123d marked `PREFERRED` in adapter list, others
  flagged `compat`; `note` clarifies "new lab work should be
  authored in build123d".
- **0.22.1**: expanded CAD detection — `_find_openscad` /
  `_find_freecad` walk Windows `Program Files\FreeCAD*\bin\` and
  macOS `/Applications` so installed-but-not-on-PATH FreeCAD is
  auto-discovered.
- **0.22.0**: new server `mcp-server-radia-interop` —
  `any_step_to_cubit_hex` (universal STEP receiver) +
  `openscad_to_cubit_hex` (CLI) + `freecad_to_cubit_hex`
  (FreeCADCmd subprocess) + `list_cad_mcp_interop`. Position:
  "the mesh backend any CAD MCP can dispatch to."

## 0.21.0 — gmsh community scrape

- New `gmsh_examples(query)` + `gmsh_examples_refresh` MCP tools.
- Sub-sources `gmsh_issues` (gitlab.onelab.info, 3000+ tickets)
  and `gmsh_stackoverflow` (StackOverflow + SciComp.SE `[gmsh]`).
- FAMILIES["gmsh"] union for ranked retrieval.

## 0.20.0 — gmsh post-processing forged

- mcp-server-gmsh-post: bundled auto-generated **gmsh API
  reference** (651 entries across `model` / `view` / `option` /
  `fltk` / …, 2 008 lines, via `_gen_api_reference.py`).
- New cookbooks: `view_data_cookbook`
  (`$NodeData`/`$ElementData`/`$ElementNodeData` decision tree)
  and `physical_groups_cookbook` (dim/tag, downstream solver
  conventions).
- New tools: `gmsh_post_api` (focused tf-idf), `gmsh_post_quality`
  (min Jacobian / skew histogram), `gmsh_post_extract_physical`,
  `gmsh_post_boundary`, `gmsh_post_add_view_from_csv` (most-frequent
  post workflow).

## 0.19.0 — build123d depth gaps closed

- Bundled auto-generated **build123d API reference** (142 classes /
  65 functions / 1 673 lines, via `_gen_api_reference.py`).
- New cookbooks: `plane_axis_location_cookbook` (the 3 most-
  confused classes, 20+ worked recipes) and
  `builder_vs_algebra_rosetta` (side-by-side conversion table).
- New tool `build123d_api(query)` for API-focused tf-idf.

## 0.18.0 — Radia-specific build123d templates + STEP gating

- 7 new templates in `generate_build123d_script`: `magnet_ring`,
  `halbach_array`, `c_core`, `e_core`, `pole_piece`,
  `stator_lamination`, `racetrack_coil`.
- `build123d_inspect_step(path)` — OCCT validity / bbox /
  micro-edge ratio / labels report; gates external STEPs before
  Cubit.
- `build123d_heal(step_in, step_out)` — `OCP.ShapeFix_Shape`
  auto-repair (small edges / face orientation / degenerate fixes).

## 0.17.0 — build123d parity with Cubit

- `lint_build123d_script` + `lint_build123d_directory` (7 rules:
  `missing-buildpart-context`, `sweep-no-path`,
  `polyline-not-closed`, `buildsketch-ambiguous-arg`,
  `missing-export`, `cadquery-in-build123d`,
  `micro-fillet-radius`).
- `build123d_suggest_next(goal, script)` — state-aware (5 goals).
- `generate_build123d_script(pattern)` — 6 starter templates
  (helix_coil, l_bracket, cae_block, gear_bd_warehouse,
  fastener_assembly, sweep_square_path).
- `build123d_try(script)` — fresh subprocess; OCCT segfault
  containment + clean namespace.
- `build123d_to_cubit_hex(script, target_size)` — one-call
  pipeline (build123d → STEP → cubit_mesh_auto → live GUI replay).
- 3 new KB topics: `joints_and_mates`, `assemblies_and_compounds`,
  `cae_workflow_tips`.

## 0.16.0 — Unified search + safety gate

- GitHub PAT auto-discovery (`GITHUB_TOKEN` / `GH_TOKEN` /
  `gh auth token`); 60 → 5000 req/h on GitHub API + GraphQL access.
- Threaded Coreform forum walk (300 topics, ~30 s on 8 threads).
- `build123d_github_discussions` via GraphQL (PAT-gated, 50
  discussions).
- `cubit_ask` / `build123d_ask` unified retrieval across
  bundled KB + scraped examples + optional live web (`include_web`).
- Pre-flight check: `cubit_exec` / `execute_build123d` scan
  failure log for similar inputs (token Jaccard ≥ 0.6) and
  surface the past hint non-blockingly.
- `cubit_mesh_auto` geometry-split rung — auto-detects compound
  bodies (`vol ≤ 3 ∧ surf/vol ≥ 7`) and `webcut volume all with
  cylinder axis z` before retrying scheme auto.

## 0.15.0 — build123d community scrape

- `build123d_discussions(query)` — `gumyr/build123d` GitHub Issues
  + comments (anonymous REST, 60 issues default).

## 0.14.x — gmsh-post lab v4.1 standardization

- **0.14.1**: lint rule `gmsh-v22-deprecated` (HIGH) — flags
  `export mesh "...msh"` without `mesh_version 4.1`. Lab policy
  is v4.1 only; `.vol` (NETGEN native) is the sole exception for
  HO curved meshes.
- **0.14.0**: new server `mcp-server-gmsh-post` —
  `gmsh_post_inspect`, `gmsh_post_validate`, `gmsh_post_convert`
  (lifts any older .msh to v4.1), `gmsh_post_write_node_data` /
  `_element_data` (append `$NodeData` / `$ElementData` blocks
  while keeping the file v4.1-compliant), `gmsh_post_spec`.
  `cubit_exec_safely` — auto-checkpoint to `.cub5`, batch dry-run
  on the snapshot, replay on live GUI on success; silent-error
  detection via `cubit.get_error_count()` delta.

## 0.13.0 — CadQuery interop

- `execute_cadquery(script)` (sibling OCCT lib) +
  `cadquery_to_cubit_hex(script)` one-call pipeline.
- `radia-mcp[cadquery]` extra; integration with cadquery-mcp
  community.

## 0.12.0 — Multi-source example unions

- FAMILY mapping: `cubit` = `[cubit, cubit_local]`; `build123d`
  = `[build123d, bd_warehouse]`.
- `cubit_local` indexer walks `S:\CoreformCubit` (lab archive of
  ~145 .jou) + `S:\Radia\01_GitHub\examples` (~400 files); 753
  files indexed.
- `bd_warehouse` (15 modules: gear, bearing, fastener, flange,
  pipe, …).
- Forum seed queries 5 → 15.

## 0.11.0 — Scraped example libraries

- `build123d_examples(query)` — `gumyr/build123d/examples` (65
  curated scripts).
- `cubit_examples(query)` — Coreform forum (Discourse search.json,
  triple-backtick code-fence extract).

## 0.10.0 — Batch ladder safety pattern

- `cubit_batch_try(commands)` — disposable headless Cubit.
- `cubit_mesh_auto(step_path)` — scheme ladder
  (auto → sweep → polyhedron → tetmesh) batch-validated, winning
  recipe replayed in live GUI. 4-turn spiral coil yielded 1668
  hex on first run.

## 0.9.0 — Failure log + tf-idf retrieval + live web docs

- Persistent jsonl failure log per kind (`cubit` / `build123d`),
  fed into every `*_lookup`.
- tf-idf retrieval with heading boost replaces substring counter.
- `cubit_web_docs` (Discourse JSON for forum.coreform.com) +
  `build123d_web_docs` (readthedocs).

## 0.8.0 — Standalone wheel crystallized

- Plan A established (Cubit GUI + PyQt5 QTimer + file-drop IPC).
- `cubit_session.py` dual-mode (gui / batch) + auto-restart on
  RPC failure.
- `cubit_checkpoint(label)` / `cubit_restore(label)` — `.cub5`
  snapshot undo.
- `cubit_mesh_diagnose` (per-volume scheme alternatives),
  `cubit_suggest_next(goal)` (state-aware), `cubit_lookup(query)`
  (heading-chunk retrieval over 8000-line knowledge).
- 4-turn coil + KEIKO 6-letter text both produced pure hex
  meshes via the build123d → Cubit pipeline.

## 0.5 / 0.6 / 0.7 — Initial wheel

- Standalone `radia-mcp` package extracted from the `radia` core
  repo (Option Y restructure).
- `mcp-server-cubit` and `mcp-server-build123d` as the first two
  entry points.
- OCP CAD Viewer retired in favor of the persistent Cubit GUI.
