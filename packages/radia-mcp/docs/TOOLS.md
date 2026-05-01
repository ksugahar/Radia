# radia-mcp Tools Inventory

Auto-generated from each server's `mcp.list_tools()` via `scripts/gen_tools_doc.py`. **Do not edit by hand** — regenerate after adding/renaming tools.

Total: **105 tools** across 9 MCP servers.

| Server (console-script) | Subpackage | Tools |
|---|---|---:|
| [`mcp-server-cubit`](#mcp-server-cubit) | `radia_mcp.cubit` | 42 |
| [`mcp-server-build123d`](#mcp-server-build123d) | `radia_mcp.build123d` | 28 |
| [`mcp-server-radia-ngsolve`](#mcp-server-radia-ngsolve) | `radia_mcp.radia_ngsolve` | 19 |
| [`mcp-server-gmsh`](#mcp-server-gmsh) | `radia_mcp.gmsh` | 6 |
| [`mcp-server-elf`](#mcp-server-elf) | `radia_mcp.elf` | 1 |
| [`mcp-server-electromagnet`](#mcp-server-electromagnet) | `radia_mcp.electromagnet` | 1 |
| [`mcp-server-ih`](#mcp-server-ih) | `radia_mcp.ih` | 2 |
| [`mcp-server-peec`](#mcp-server-peec) | `radia_mcp.peec` | 1 |
| [`mcp-server-radia-interop`](#mcp-server-radia-interop) | `radia_mcp.interop` | 5 |

## `mcp-server-cubit`

_Cubit hex-mesh export, Netgen/NGSolve curving, scripting + API reference_

Module: `radia_mcp.cubit.server`

| Tool | Description |
|---|---|
| `cubit_ask` | One-shot search across every Cubit knowledge surface we have. |
| `cubit_batch_try` | Dry-run a recipe in a fresh headless Cubit subprocess. |
| `cubit_checkpoint` | Save the current Cubit session state as a named checkpoint. |
| `cubit_curate_learned_recipes` | **Lab maintainer tool**: read accumulated `learned_recipes.jsonl`, |
| `cubit_diagnostics_guide` | Get the foundational mesh-diagnostics + cleanup + quality playbook. |
| `cubit_docs` | Get Cubit documentation: export formats, scripting guide, and API reference. |
| `cubit_examples` | Search Cubit journal examples from **multiple unioned sources**. |
| `cubit_examples_refresh` | Force-refresh every Cubit example sub-source. |
| `cubit_exec` | Send arbitrary Cubit commands to the persistent viewer session. |
| `cubit_exec_safely` | Execute commands against the live GUI **with a pre-save and a |
| `cubit_forum_tips` | Get practical Cubit meshing tips sourced from the Coreform forum. |
| `cubit_generate_dialog` | Generate a single complete PySide6 dialog script for a Cubit toolbar |
| `cubit_list_checkpoints` | List all saved Cubit checkpoints with size + mtime. |
| `cubit_lookup` | Retrieve relevant sections from the bundled Cubit knowledge base. |
| `cubit_mesh_apply_choice` | Apply the human's chosen variant from a previous |
| `cubit_mesh_auto` | Find a working mesh recipe by trying a scheme ladder in batch, |
| `cubit_mesh_diagnose` | Per-volume meshing diagnostic. |
| `cubit_mesh_race` | Race N recipes in parallel batch Cubits, replay the first/best |
| `cubit_mesh_race_review` | Race N AI variants + observe live human, **wait for all** to |
| `cubit_mesh_race_review_async` | **Background launch** of a race review — returns immediately |
| `cubit_mesh_race_smart` | **The flagship workflow**: AI inspects the live Cubit state, |
| `cubit_mesh_race_smart_async` | **Background variant of `cubit_mesh_race_smart`** — AI inspects |
| `cubit_mesh_race_status` | Check the status of a background race launched by |
| `cubit_mesh_race_with_human` | **The radia-mcp signature workflow.** |
| `cubit_probe` | Query the Cubit session for geometry/mesh statistics. |
| `cubit_recent_failures` | Return the last N failed Cubit invocations (from persistent log). |
| `cubit_restore` | Restore a previously-saved Cubit checkpoint by label. |
| `cubit_scaffold_toolbar` | Generate a complete Coreform-Cubit custom-toolbar skeleton on disk. |
| `cubit_session_shutdown` | Stop the persistent Cubit daemon. Next `cubit_show` relaunches. |
| `cubit_session_status` | Return diagnostic info about the Cubit session (pid, alive, bin_dir). |
| `cubit_show` | Load a file into the **persistent Cubit viewer** and optionally run |
| `cubit_snapshot` | Hardcopy the current Cubit view to a PNG file. |
| `cubit_suggest_next` | Suggest concrete next Cubit commands based on the current state. |
| `cubit_toolbar_guide` | Get documentation on building custom in-Cubit PySide6 toolbars. |
| `cubit_web_docs` | Fetch live Cubit documentation and grep for `query`. |
| `generate_cubit_script` | Generate a template Cubit Python script for common workflows. |
| `get_lint_rules` | List all available Cubit export lint rules with descriptions. |
| `lint_cubit_directory` | Lint all Python scripts in a directory for Cubit export convention violations. |
| `lint_cubit_script` | Lint a Python script for Cubit mesh export convention violations. |
| `netgen_code_example` | Get a ready-to-run Netgen export code example. |
| `netgen_workflow_guide` | Get step-by-step Netgen/NGSolve workflow documentation. |
| `open_in_cubit` | Open a file, or execute a list of commands, in **Cubit GUI**. |

## `mcp-server-build123d`

_build123d (Pythonic OCCT) + STEP/XCAF labels + Cubit pipeline interop_

Module: `radia_mcp.build123d.server`

| Tool | Description |
|---|---|
| `build123d_api` | Search the bundled build123d API reference (auto-generated from |
| `build123d_ask` | One-shot search across every build123d knowledge surface. |
| `build123d_discussions` | Search the build123d GitHub Issues archive (de-facto forum). |
| `build123d_examples` | Search build123d + **bd_warehouse** + **GitHub Issues** (unioned). |
| `build123d_examples_refresh` | Force-refresh all build123d sources from GitHub + YouTube. |
| `build123d_heal` | Run OCCT `ShapeFix_Shape` on a STEP file, write a healed copy. |
| `build123d_inspect_step` | Inspect an external STEP file via the build123d / OCCT importer |
| `build123d_lookup` | Retrieve relevant sections from the bundled build123d knowledge + |
| `build123d_recent_failures` | Return the last N failed `execute_build123d` invocations (from log). |
| `build123d_suggest_next` | Suggest concrete next build123d steps toward a goal. |
| `build123d_to_cubit_hex` | End-to-end: build123d script → STEP → `cubit_mesh_auto` (batch- |
| `build123d_try` | Dry-run a build123d script in a **fresh Python subprocess**. |
| `build123d_try_race` | Race N build123d script variants in parallel subprocesses, |
| `build123d_usage` | Get build123d CAD modeling documentation for CAE workflows. |
| `build123d_web_docs` | Fetch live build123d documentation (readthedocs) and grep for `query`. |
| `cadquery_to_cubit_hex` | End-to-end: cadquery script → STEP → Cubit `cubit_mesh_auto` |
| `execute_build123d` | Execute a build123d Python script and return geometry information. |
| `execute_cadquery` | Execute a **CadQuery** Python script (sibling OCCT CAD lib). |
| `generate_build123d_script` | Return a ready-to-run build123d boilerplate for a common pattern. |
| `generate_helix_coil` | Generate a helical coil with optional cross-section taper. |
| `inspect_geometry` | Inspect a STEP or BREP file for CAE quality. |
| `lint_build123d_directory` | Lint every `.py` file under `directory` (non-recursive default). |
| `lint_build123d_script` | Static-analysis lint for a single build123d `.py` script. |
| `preview_boolean` | Perform a boolean (union/subtract/intersect) of two build123d |
| `preview_extrude` | Extrude a build123d sketch by `amount` and show in Cubit. |
| `preview_shape_in_cubit` | Run a build123d script and show the resulting Shape in the |
| `preview_text` | Generate 3D extruded text and send it to the live Cubit viewer. |
| `section_along_path` | Section a STEP/BREP coil solid along a discrete path and extract |

## `mcp-server-radia-ngsolve`

_Radia + NGSolve coupled magnetostatics, Kelvin transformation, sparse solver_

Module: `radia_mcp.radia_ngsolve.server`

| Tool | Description |
|---|---|
| `analytical_formulas` | Get documentation for radia.analytical_formulas (closed-form reference layer). |
| `esim` | Get ESIM (Effective Surface Impedance Method) general documentation. |
| `get_radia_lint_rules` | List all available NGSolve lint rules with descriptions. |
| `gmsh_post_spec` | GMSH post-processing specification for Radia panels. |
| `install_deploy` | Radia install / deploy policy and recipes — 3-tier configuration |
| `kelvin_transformation` | Get Kelvin transformation documentation for open boundary FEM problems. |
| `lint_radia_directory` | Lint all Python scripts in a directory for NGSolve convention violations. |
| `lint_radia_script` | Lint a single Python script for Radia + NGSolve convention violations. |
| `md2html_usage` | Get md2html converter documentation (MathJax, reference links, styled HTML). |
| `ngsbem_inductance` | Get ngsolve.bem boundary element method documentation for inductance extraction. |
| `ngsolve_usage` | Get NGSolve finite element library usage documentation. |
| `panel_add_param` | Plan where to add a new parameter to a Radia-NGSolve panel. |
| `panel_describe_jp` | 現在のパネルソースを AST 解析して日本語で詳細に説明する。 |
| `panel_gui_pitfalls` | Pitfalls and lessons learned from Radia GUI / Cubit panel development. |
| `panel_schema` | Show Radia-NGSolve panel definitions with Japanese labels and physics. |
| `panel_widget_locations` | Return file:line locations for everything that touches a widget. |
| `peec_inductance` | Get documentation for the Radia PEEC-inductance (coil only, STEP) panel mode. |
| `radia_usage` | Get Radia C++ library usage documentation. |
| `sparsesolv` | Get sparsesolv documentation and code examples (now in ngsolve.la). |

## `mcp-server-gmsh`

_Gmsh script linting + post-processing spec helpers_

Module: `radia_mcp.gmsh.server`

| Tool | Description |
|---|---|
| `get_gmsh_lint_rules` | List all available GMSH lint rules with descriptions. |
| `gmsh_examples` | Get GMSH tutorial and example documentation. |
| `gmsh_reference` | Get GMSH technical reference (options, algorithms, fields, formats). |
| `gmsh_usage` | Get GMSH documentation for visualization and post-processing. |
| `lint_gmsh_directory` | Lint all Python scripts in a directory for GMSH policy violations. |
| `lint_gmsh_script` | Lint a single Python script for GMSH policy violations. |

## `mcp-server-elf`

_ELF (Electromagnetic Loss/Field) postprocessing knowledge_

Module: `radia_mcp.elf.server`

| Tool | Description |
|---|---|
| `elf_usage` | Get ELF600 electromagnetic field analysis documentation. |

## `mcp-server-electromagnet`

_Electromagnet design (symmetry reductions, BC choices)_

Module: `radia_mcp.electromagnet.server`

| Tool | Description |
|---|---|
| `electromagnet_usage` | Get accelerator electromagnet analysis documentation. |

## `mcp-server-ih`

_IH (induction-heating) coil + load workflow_

Module: `radia_mcp.ih.server`

| Tool | Description |
|---|---|
| `ih_sibc` | Get IH solver architecture and SIBC documentation. |
| `induction_heating` | Get induction heating simulation documentation. |

## `mcp-server-peec`

_PEEC (partial element equivalent circuit) inductance modeling_

Module: `radia_mcp.peec.server`

| Tool | Description |
|---|---|
| `peec_usage` | Get PEEC (Partial Element Equivalent Circuit) documentation. |

## `mcp-server-radia-interop`

_Cross-tool interop (CadQuery / build123d / Cubit STEP boundary)_

Module: `radia_mcp.interop.server`

| Tool | Description |
|---|---|
| `any_step_to_cubit_hex` | Universal CAD-MCP mesh backend: accept ANY STEP file and run it |
| `freecad_exec_safely` | Cubit-style safety pattern for FreeCAD: snapshot → batch dry-run |
| `freecad_to_cubit_hex` | Execute a FreeCAD script in a FreeCADCmd subprocess, export the |
| `list_cad_mcp_interop` | List registered CAD-MCP interop adapters + their availability. |
| `openscad_to_cubit_hex` | Execute OpenSCAD code, export STEP, run through `cubit_mesh_auto`. |

