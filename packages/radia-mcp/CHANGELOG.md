# Changelog

All notable changes to `radia-mcp`. Format: each release lists **what
shipped** + **why** in compact form. Older releases (≤ 0.4) are
omitted; the 0.5 → 0.6 jump is when the standalone `radia-mcp` wheel
crystallized as its own package.

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
