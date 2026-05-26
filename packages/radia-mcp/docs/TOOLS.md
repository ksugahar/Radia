# radia-mcp Tools Inventory

Auto-generated from each server's `mcp.list_tools()` via `scripts/gen_tools_doc.py`. **Do not edit by hand** — regenerate after adding/renaming tools.

Total: **417 tools** across 40 MCP servers.

| Server (console-script) | Subpackage | Tools |
|---|---|---:|
| [`mcp-server-cubit`](#mcp-server-cubit) | `radia_mcp.cubit` | 44 |
| [`mcp-server-build123d`](#mcp-server-build123d) | `radia_mcp.build123d` | 29 |
| [`mcp-server-radia-interop`](#mcp-server-radia-interop) | `radia_mcp.interop` | 8 |
| [`mcp-server-gmsh`](#mcp-server-gmsh) | `radia_mcp.gmsh` | 7 |
| [`mcp-server-radia-ngsolve`](#mcp-server-radia-ngsolve) | `radia_mcp.radia_ngsolve` | 32 |
| [`mcp-server-fem`](#mcp-server-fem) | `radia_mcp.fem` | 12 |
| [`mcp-server-bem`](#mcp-server-bem) | `radia_mcp.bem` | 9 |
| [`mcp-server-matrix-solvers`](#mcp-server-matrix-solvers) | `radia_mcp.matrix_solvers` | 6 |
| [`mcp-server-mor`](#mcp-server-mor) | `radia_mcp.mor` | 9 |
| [`mcp-server-ih`](#mcp-server-ih) | `radia_mcp.ih` | 4 |
| [`mcp-server-peec`](#mcp-server-peec) | `radia_mcp.peec` | 4 |
| [`mcp-server-electromagnet`](#mcp-server-electromagnet) | `radia_mcp.electromagnet` | 3 |
| [`mcp-server-motor`](#mcp-server-motor) | `radia_mcp.motor` | 11 |
| [`mcp-server-accelerator`](#mcp-server-accelerator) | `radia_mcp.accelerator` | 3 |
| [`mcp-server-fusion`](#mcp-server-fusion) | `radia_mcp.fusion` | 3 |
| [`mcp-server-magnetic-materials`](#mcp-server-magnetic-materials) | `radia_mcp.magnetic_materials` | 7 |
| [`mcp-server-litz-transmission`](#mcp-server-litz-transmission) | `radia_mcp.litz_transmission` | 3 |
| [`mcp-server-rna-mec`](#mcp-server-rna-mec) | `radia_mcp.rna_mec` | 3 |
| [`mcp-server-topology-optimization`](#mcp-server-topology-optimization) | `radia_mcp.topology_optimization` | 4 |
| [`mcp-server-optuna`](#mcp-server-optuna) | `radia_mcp.optuna` | 6 |
| [`mcp-server-bayesian-opt`](#mcp-server-bayesian-opt) | `radia_mcp.bayesian_opt` | 3 |
| [`mcp-server-evolutionary`](#mcp-server-evolutionary) | `radia_mcp.evolutionary` | 3 |
| [`mcp-server-data-assimilation`](#mcp-server-data-assimilation) | `radia_mcp.data_assimilation` | 3 |
| [`mcp-server-gnn`](#mcp-server-gnn) | `radia_mcp.gnn` | 3 |
| [`mcp-server-pinn`](#mcp-server-pinn) | `radia_mcp.pinn` | 3 |
| [`mcp-server-wpt`](#mcp-server-wpt) | `radia_mcp.wpt` | 7 |
| [`mcp-server-ndt`](#mcp-server-ndt) | `radia_mcp.ndt` | 4 |
| [`mcp-server-metamaterial`](#mcp-server-metamaterial) | `radia_mcp.metamaterial` | 4 |
| [`mcp-server-nmr-mri`](#mcp-server-nmr-mri) | `radia_mcp.nmr_mri` | 2 |
| [`mcp-server-maglev-linear`](#mcp-server-maglev-linear) | `radia_mcp.maglev_linear` | 3 |
| [`mcp-server-team-benchmark`](#mcp-server-team-benchmark) | `radia_mcp.team_benchmark` | 7 |
| [`mcp-server-differential-forms`](#mcp-server-differential-forms) | `radia_mcp.differential_forms` | 13 |
| [`mcp-server-mathematica`](#mcp-server-mathematica) | `radia_mcp.mathematica` | 11 |
| [`mcp-server-md2html`](#mcp-server-md2html) | `radia_mcp.md2html` | 2 |
| [`mcp-server-paper-writing`](#mcp-server-paper-writing) | `radia_mcp.paper_writing` | 92 |
| [`mcp-server-graph`](#mcp-server-graph) | `radia_mcp.graph` | 8 |
| [`mcp-server-chart2d`](#mcp-server-chart2d) | `radia_mcp.chart2d` | 24 |
| [`mcp-server-literature-index`](#mcp-server-literature-index) | `radia_mcp.literature_index` | 9 |
| [`mcp-server-radia-meta`](#mcp-server-radia-meta) | `radia_mcp.meta` | 6 |
| [`mcp-server-panel-review`](#mcp-server-panel-review) | `radia_mcp.panel_review` | 3 |

## `mcp-server-cubit`

_Cubit mesh scripting, hex/tet workflow, export formats_

Module: `radia_mcp.cubit.server`

| Tool | Description |
|---|---|
| `cubit_ask` | One-shot search across every Cubit knowledge surface we have. |
| `cubit_batch_try` | Dry-run a recipe in a fresh headless Cubit subprocess. |
| `cubit_checkpoint` | Save the current Cubit session state as a named checkpoint. |
| `cubit_cpp_sdk_guide` | Get documentation on building Cubit C++ SDK plugins. |
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
| `cubit_status` | (no description) |
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

_build123d STEP authoring (CAD-as-code) + Cubit interop_

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
| `build123d_status` | (no description) |
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

## `mcp-server-radia-interop`

_Cross-CAD interop (STEP/IGES/CadQuery <-> Cubit/Netgen)_

Module: `radia_mcp.interop.server`

| Tool | Description |
|---|---|
| `any_step_to_cubit_hex` | Universal CAD-MCP mesh backend: accept ANY STEP file and run it |
| `freecad_exec_safely` | Cubit-style safety pattern for FreeCAD: snapshot → batch dry-run |
| `freecad_to_cubit_hex` | Execute a FreeCAD script in a FreeCADCmd subprocess, export the |
| `interop_comsol_lab_tips` | Sugahara Lab (Kindai University) COMSOL practical tips compendium. |
| `interop_comsol_livelink` | COMSOL LiveLink (Java + MATLAB + MPh Python) knowledge. |
| `list_cad_mcp_interop` | List registered CAD-MCP interop adapters + their availability. |
| `openscad_to_cubit_hex` | Execute OpenSCAD code, export STEP, run through `cubit_mesh_auto`. |
| `radia_interop_status` | (no description) |

## `mcp-server-gmsh`

_GMSH MSH v4.1 inspect/validate/convert/write_node_data_

Module: `radia_mcp.gmsh.server`

| Tool | Description |
|---|---|
| `get_gmsh_lint_rules` | List all available GMSH lint rules with descriptions. |
| `gmsh_examples` | Get GMSH tutorial and example documentation. |
| `gmsh_reference` | Get GMSH technical reference (options, algorithms, fields, formats). |
| `gmsh_status` | (no description) |
| `gmsh_usage` | Get GMSH documentation for visualization and post-processing. |
| `lint_gmsh_directory` | Lint all Python scripts in a directory for GMSH policy violations. |
| `lint_gmsh_script` | Lint a single Python script for GMSH policy violations. |

## `mcp-server-radia-ngsolve`

_Radia + NGSolve: Kelvin / sparsesolv / CLN / PEEC / analytical formulas / lint_

Module: `radia_mcp.radia_ngsolve.server`

| Tool | Description |
|---|---|
| `analytical_formulas` | Get documentation for radia.analytical_formulas (closed-form reference layer). |
| `axifemm_documentation` | Get radia-axifemm documentation: Henrotte axisymmetric Q-element FE |
| `basis_functions` | Finite-element basis function library — Mathematica-canonical |
| `bem_cln` | Get BEM-CLN (per-element multipole CLN with Schur-F termination) |
| `cln_3d` | Get 3D Cauer Ladder Network (CLN) / Kameari-Tanimoto iteration |
| `cln_3d_notebook` | Retrieve Tanimoto's raw 3D CLN notebook Python code. |
| `cln_sibc_orthogonal` | Get CLN expansion-point + SIBC orthogonal-residual theory documentation. |
| `cln_sphere_dd_pipeline` | Get the Sphere DD (double-double, ~32 digit) VIM Cauer Ladder Network |
| `esim` | Get ESIM (Effective Surface Impedance Method) general documentation. |
| `get_radia_lint_rules` | List all available NGSolve lint rules with descriptions. |
| `gmsh_post_spec` | GMSH post-processing specification for Radia panels. |
| `install_deploy` | Radia install / deploy policy and recipes — 3-tier configuration |
| `kelvin_identify_post_hoc` | Add Kelvin Periodic Identifications to an existing NGSolve mesh |
| `kelvin_transformation` | Get Kelvin transformation documentation for open boundary FEM problems. |
| `lint_radia_directory` | Lint all Python scripts in a directory for NGSolve convention violations. |
| `lint_radia_script` | Lint a single Python script for Radia + NGSolve convention violations. |
| `md2html_usage` | Get md2html converter documentation (MathJax, reference links, styled HTML). |
| `mmm_core` | MMM (Magnetic Moment Method) core theory + Radia heritage. |
| `ngsbem_inductance` | Get ngsolve.bem boundary element method documentation for inductance extraction. |
| `ngsolve_usage` | Get NGSolve finite element library usage documentation. |
| `panel_add_param` | Plan where to add a new parameter to a Radia-NGSolve panel. |
| `panel_describe_jp` | 現在のパネルソースを AST 解析して日本語で詳細に説明する。 |
| `panel_gui_pitfalls` | Pitfalls and lessons learned from Radia GUI / Cubit panel development. |
| `panel_schema` | Show Radia-NGSolve panel definitions with Japanese labels and physics. |
| `panel_widget_locations` | Return file:line locations for everything that touches a widget. |
| `peec_inductance` | Get documentation for the Radia PEEC-inductance (coil only, STEP) panel mode. |
| `radia_ngsolve_status` | (no description) |
| `radia_usage` | Get Radia C++ library usage documentation. |
| `release_workflow` | Triple-package release workflow for the Radia monorepo |
| `sparsesolv` | Get sparsesolv documentation and code examples. |
| `standalone_panels` | Cubit-bypass standalone launch of the four Radia-NGSolve panels — |
| `topology_optimization` | Topology optimization knowledge mirror for NGSolve users. |

## `mcp-server-fem`

_FEM formulations theory layer (A-Omega / T-Omega / H / Reduced / Darwin, edge / HO / XFEM / IGA / DG, gauging + Kelvin, MSFEM, Schur circuit coupling, NGSolve hierarchical)_

Module: `radia_mcp.fem.server`

| Tool | Description |
|---|---|
| `fem_elements` | Element technology: edge (Nedelec), high-order, XFEM, isogeometric, DG. |
| `fem_equivalence_source` | Equivalence-theorem near-field source (CST Near-Field Source equivalent). |
| `fem_gauge_open_boundary` | Gauging + open boundary techniques. |
| `fem_large_scale_special` | Large-scale, error theory, multi-scale (Hollaus MSFEM), misc techniques. |
| `fem_ngsolve_hierarchy` | NGSolve hierarchical H(curl) bases - Zaglmayr / nograds / tree-cotree. |
| `fem_nonconforming_mesh_coupling` | Non-conforming mesh coupling: mortar / Nitsche / FETI-DP / BDDC / DG / |
| `fem_overview` | FEM landscape: lab stack, decision tree, genealogy. |
| `fem_potential_formulations` | Potential formulations: A-Omega, T-Omega, H, Reduced, Darwin. |
| `fem_status` | (no description) |
| `fem_time_domain_axisym` | Time-domain, axisymmetric (Henrotte), harmonic balance, HF, circuit coupling. |
| `fem_xfem_comsol` | XFEM in COMSOL Multiphysics (Jafari-Broumand-Vahab-Khalili 2021). |
| `fem_xfem_em_hiruma` | EM-XFEM (Hiruma 2023): electromagnetic XFEM for eddy-current |

## `mcp-server-bem`

_MoM/BEM theory: RWG, EFIE/MFIE/CFIE/PMCHWT, Loop-Star, Calderon, Radia MMM/MSC, HACApK, FEM-BEM_

Module: `radia_mcp.bem.server`

| Tool | Description |
|---|---|
| `bem_fem_bem_hybrid` | FEM-BEM hybrid methods for open-boundary EM. |
| `bem_h_matrix` | H-matrix / ACA acceleration for BEM. |
| `bem_low_freq` | Low-frequency BEM stabilization. |
| `bem_mmm_msc` | Magnetic Moment Method (MMM) and Surface Charge (MSC) -- Radia's core. |
| `bem_mom_foundations` | MoM foundations: Harrington 1968, RWG 1982, wire-grid (NEC). |
| `bem_overview` | BEM/MoM landscape: lab stack, decision tree, genealogy. |
| `bem_sommerfeld_layered` | Sommerfeld integrals / layered-medium Green's function theory. |
| `bem_status` | (no description) |
| `bem_surface_ie` | Surface Integral Equations: EFIE, MFIE, CFIE, PMCHWT. |

## `mcp-server-matrix-solvers`

_Sparse solver theory + decision tree: Krylov (CG/BiCGSTAB/GMRES/COCG/COCR/IDR), preconditioners (AMG, Hiptmair-Xu AMS), Biro-Preis A-V, tree-cotree_

Module: `radia_mcp.matrix_solvers.server`

| Tool | Description |
|---|---|
| `matrix_solvers_direct` | Direct sparse solvers: LU, PARDISO, MUMPS, SuperLU. |
| `matrix_solvers_em_specific` | EM-specific solver / formulation choices. |
| `matrix_solvers_krylov` | Krylov subspace methods: CG, BiCGSTAB, GMRES, COCG, COCR, IDR(s). |
| `matrix_solvers_overview` | Solver landscape: lab stack, decision tree, genealogy. |
| `matrix_solvers_preconditioners` | Preconditioner catalog: classical, AMG, AMS (Hiptmair-Xu). |
| `matrix_solvers_status` | (no description) |

## `mcp-server-mor`

_Model Order Reduction: PRIMA, Cauer Ladder Network, hyperreduction (DEIM)_

Module: `radia_mcp.mor.server`

| Tool | Description |
|---|---|
| `mor_bibliography` | Search the MOR bibliography catalog (87 papers in lab library). |
| `mor_cln` | Cauer Ladder Network (CLN) -- the lab-specialty MOR for eddy-current FEM. |
| `mor_cln_advanced` | CLN advanced physics extensions -- HF, nonlinear, circuit, |
| `mor_cln_collab` | CLN external-collaborator work + Cauer-form conversion + |
| `mor_cln_multiport` | CLN multi-port + multi-expansion-point + 3D extensions. |
| `mor_cln_practice` | CLN MATLAB+COMSOL practice corpus -- foundations + 2020_11_04 lab |
| `mor_cln_specialty` | CLN lab-signature techniques -- termination, Hiruma method, |
| `mor_status` | (no description) |
| `mor_systematic` | Systematic MOR knowledge -- distilled from the deGruyter 3-volume |

## `mcp-server-ih`

_Induction heating: SIBC, ESIM, Karl iteration, workpiece coupling_

Module: `radia_mcp.ih.server`

| Tool | Description |
|---|---|
| `ih_esim` | Induction-heating ESIM (Effective Surface Impedance Method) usage. |
| `ih_sibc` | Get IH solver architecture and SIBC documentation. |
| `ih_status` | (no description) |
| `induction_heating` | Get induction heating simulation documentation. |

## `mcp-server-peec`

_PEEC filament/panel, FastHenry parser, HOIBC, Carstensen AC copper loss_

Module: `radia_mcp.peec.server`

| Tool | Description |
|---|---|
| `peec_carstensen_ac_loss` | Carstensen-Dowell analytical AC copper-loss formulas for stranded |
| `peec_hoibc` | HOIBC (Higher Order Impedance Boundary Conditions) — extension of |
| `peec_status` | (no description) |
| `peec_usage` | Get PEEC (Partial Element Equivalent Circuit) documentation. |

## `mcp-server-electromagnet`

_Accelerator electromagnet: CoilBuilder, Hantila, Play/Energy hysteresis_

Module: `radia_mcp.electromagnet.server`

| Tool | Description |
|---|---|
| `electromagnet_status` | (no description) |
| `electromagnet_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 13 topics. |
| `electromagnet_usage` | Get accelerator electromagnet analysis documentation. |

## `mcp-server-motor`

_Motor analysis: ONELAB transient, Hollaus effective material (lamination), Wakao autoencoder topology, Kaimori-Mifune Darwin TD_

Module: `radia_mcp.motor.server`

| Tool | Description |
|---|---|
| `motor_bibliography` | Search the motor analysis bibliography catalog. |
| `motor_darwin_model` | Darwin-model time-domain formulation (capacitive + inductive coupling). |
| `motor_em_force_extras` | Forward to `differential_forms_em_force_extras` -- advanced EM force |
| `motor_em_force_recipe` | Practical NGSolve EM-force recipe for motor analysis. |
| `motor_femm_transient` | FEMM newbuild transient solver — Lange-Henrotte-Hameyer 2009 |
| `motor_henrotte_lineage` | The Henrotte–Hameyer–RWTH research arc (energy-consistent E&M FE). |
| `motor_hollaus_eddy` | Karl Hollaus / TU Wien MSFEM for laminated-iron eddy currents. |
| `motor_hollaus_genealogy` | Visualize the Karl Hollaus / TU Wien MSFEM research genealogy |
| `motor_onelab` | ONELAB/GetDP electric-machine reference template knowledge. |
| `motor_status` | (no description) |
| `motor_topology_optimization` | SynRM topology optimization (Wakao 2025 autoencoder + level-set). |

## `mcp-server-accelerator`

_Accelerator physics: beam optics, dipole/quad/sext magnets, undulator/wiggler_

Module: `radia_mcp.accelerator.server`

| Tool | Description |
|---|---|
| `accelerator` | Accelerator magnet design with Radia + radia-mcp. |
| `accelerator_status` | (no description) |
| `accelerator_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 4 topics. |

## `mcp-server-fusion`

_Fusion reactor magnets: tokamak ITER + stellarator LHD/W7-X/heliotron lineage_

Module: `radia_mcp.fusion.server`

| Tool | Description |
|---|---|
| `fusion` | Fusion reactor magnet knowledge. |
| `fusion_status` | (no description) |
| `fusion_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 12 topics. |

## `mcp-server-magnetic-materials`

_Magnetic materials: hysteresis (Play/Energy lab core), iron loss (Bertotti/Steinmetz/iGSE), JIS silicon steel, PM datasheets, Osborn demag factor_

Module: `radia_mcp.magnetic_materials.server`

| Tool | Description |
|---|---|
| `magnetic_materials_demagnetization` | Demagnetization factor N (反磁場係数): Osborn 1945 closed-form |
| `magnetic_materials_hysteresis` | Hysteresis model catalog & decision tree. |
| `magnetic_materials_iron_loss` | Iron loss models: Steinmetz family, Bertotti 3-term, Carstensen, |
| `magnetic_materials_permanent_magnet` | Permanent magnet datasheets: NdFeB, SmCo, Ferrite, AlNiCo |
| `magnetic_materials_radia_status` | Radia magnetic material implementation status (Mat classes). |
| `magnetic_materials_silicon_steel` | JIS silicon steel grade database + processing/handling notes. |
| `magnetic_materials_status` | (no description) |

## `mcp-server-litz-transmission`

_Litz wire AC loss (Dowell, homogenization, magnetic-plated wire) + multiconductor transmission line theory_

Module: `radia_mcp.litz_transmission.server`

| Tool | Description |
|---|---|
| `litz_transmission` | Litz wire + transmission line knowledge. |
| `litz_transmission_status` | (no description) |
| `litz_transmission_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 14 topics. |

## `mcp-server-rna-mec`

_RNA / Magnetic Equivalent Circuit. ★ Lab specialty: dynamic hysteresis MEC (Play + Cauer)_

Module: `radia_mcp.rna_mec.server`

| Tool | Description |
|---|---|
| `rna_mec` | Reluctance Network Analysis / Magnetic Equivalent Circuit. |
| `rna_mec_status` | (no description) |
| `rna_mec_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 13 topics. |

## `mcp-server-topology-optimization`

_Topology optimization: SIMP, level set, ON/OFF, MMA, Wakao autoencoder+LS SynRM_

Module: `radia_mcp.topology_optimization.server`

| Tool | Description |
|---|---|
| `topology_opt_applications` | Practical applications. |
| `topology_opt_shape_optimization` | Shape optimization for nonlinear magnetostatics. |
| `topology_opt_topology_derivative` | Topological derivative for changing topology (adding/removing material). |
| `topology_optimization_status` | (no description) |

## `mcp-server-optuna`

_Optuna black-box optimization (Sano-Akiba-Imamura 2023 textbook)_

Module: `radia_mcp.optuna.server`

| Tool | Description |
|---|---|
| `optuna_algorithm` | Optuna algorithm internals: samplers, MO, constraints, pruning. |
| `optuna_kanamori2016_textbook` | Kanamori et al. (2016) continuous-optimization textbook companion. |
| `optuna_lab_applications` | Lab applications: how Optuna plugs into Radia / NGSolve work. |
| `optuna_recipes_advanced` | Advanced lab BBO recipes that wire Optuna onto a Stage-2 calc_*.py. |
| `optuna_status` | (no description) |
| `optuna_usage` | Optuna usage: basics, storage, visualization. |

## `mcp-server-bayesian-opt`

_BO + GP regression + FMQA + surrogate models (57 lab files; ARD kernel, PI-GP, multi-fidelity)_

Module: `radia_mcp.bayesian_opt.server`

| Tool | Description |
|---|---|
| `bayesian_opt` | Bayesian optimization, GP regression, FMQA, surrogate models. |
| `bayesian_opt_status` | (no description) |
| `bayesian_opt_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 8 topics. |

## `mcp-server-evolutionary`

_GA / DE / PSO / CMA-ES / Immune / NSGA-II for EM_

Module: `radia_mcp.evolutionary.server`

| Tool | Description |
|---|---|
| `evolutionary` | Evolutionary computation algorithms for EM optimization. |
| `evolutionary_status` | (no description) |
| `evolutionary_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 7 topics. |

## `mcp-server-data-assimilation`

_Kalman / EnKF / 4D-Var for EM state estimation + sensor fusion_

Module: `radia_mcp.data_assimilation.server`

| Tool | Description |
|---|---|
| `data_assimilation` | Data assimilation for EM state estimation + sensor fusion. |
| `data_assimilation_status` | (no description) |
| `data_assimilation_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 4 topics. |

## `mcp-server-gnn`

_Graph Neural Networks for PDE/EM. Physics-Embedded GNN, E(n)-GNN / NequIP / MACE_

Module: `radia_mcp.gnn.server`

| Tool | Description |
|---|---|
| `gnn` | Graph Neural Networks for PDE / EM problems. |
| `gnn_status` | (no description) |
| `gnn_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 4 topics. |

## `mcp-server-pinn`

_Physics-Informed Neural Networks + Gaussian Processes for EM_

Module: `radia_mcp.pinn.server`

| Tool | Description |
|---|---|
| `pinn` | Physics-Informed Neural Networks (PINN) + Gaussian Processes (PI-GP) |
| `pinn_status` | (no description) |
| `pinn_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 5 topics. |

## `mcp-server-wpt`

_Wireless Power Transfer: coil + compensation (SS/LCC/LCL), efficiency, IEC 61980 / SAE J2954, FOD, dynamic EV / robot / bearingless motor, capacitive / microwave / metamaterial_

Module: `radia_mcp.wpt.server`

| Tool | Description |
|---|---|
| `wpt_alternatives` | Alternative WPT: capacitive, microwave/rectenna, metamaterial. |
| `wpt_applications` | WPT applications: dynamic EV, robot, bearingless motor. |
| `wpt_coil_compensation` | Coil design + compensation topology + resonance matching. |
| `wpt_efficiency_safety` | Efficiency (Q, k, kQ) + safety + IEC/SAE standards. |
| `wpt_fod` | ★ Foreign Object Detection (FOD) — lab core research. |
| `wpt_overview` | WPT landscape: regimes, decision tree, lab focus. |
| `wpt_status` | (no description) |

## `mcp-server-ndt`

_Non-destructive testing: eddy current testing, magnetic flux leakage, MFL signal analysis_

Module: `radia_mcp.ndt.server`

| Tool | Description |
|---|---|
| `ndt` | Electromagnetic non-destructive testing (NDT / NDE) knowledge. |
| `ndt_bibliography` | Search the NDT/NDE bibliography catalog. |
| `ndt_status` | (no description) |
| `ndt_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 14 topics. |

## `mcp-server-metamaterial`

_Metamaterials: homogenization, effective medium, periodic structures_

Module: `radia_mcp.metamaterial.server`

| Tool | Description |
|---|---|
| `metamaterial` | Electromagnetic metamaterial knowledge. |
| `metamaterial_bibliography` | Search the metamaterial bibliography catalog of cited PDFs. |
| `metamaterial_status` | (no description) |
| `metamaterial_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 12 topics. |

## `mcp-server-nmr-mri`

_NMR/MRI: gradient coils, B0 shimming, RF coils, field uniformity_

Module: `radia_mcp.nmr_mri.server`

| Tool | Description |
|---|---|
| `nmr_mri_bibliography` | Search the NMR/MRI bibliography catalog. |
| `nmr_mri_status` | (no description) |

## `mcp-server-maglev-linear`

_Maglev (EMS/EDS/SCMaglev/Halbach/bearingless ★) + linear drives (LIM/LSM). Lab specialty: bearingless + WPT_

Module: `radia_mcp.maglev_linear.server`

| Tool | Description |
|---|---|
| `maglev_linear` | Magnetic levitation + linear drive knowledge. |
| `maglev_linear_status` | (no description) |
| `maglev_linear_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 10 topics. |

## `mcp-server-team-benchmark`

_TEAM Workshop benchmark problems reference layer (30 problems × physics class). ★ Lab core: 13, 20, 23, 32, 33b_

Module: `radia_mcp.team_benchmark.server`

| Tool | Description |
|---|---|
| `team_benchmark_status` | (no description) |
| `team_catalog` | TEAM Workshop benchmark catalog. |
| `team_eddy_current` | TEAM eddy current problems: 1a/1b, 2, 3, 4, 5, 7, 9, 21, 24. |
| `team_force_motion` | TEAM force / motion / levitation: 17, 20, 23, 28, 33b ★ lab core. |
| `team_magnetostatic` | TEAM magnetostatic problems: 6 (sphere), 13 (nonlinear yoke). |
| `team_ndt_inverse` | TEAM NDT and inverse / optimization problems. |
| `team_special` | TEAM special problems: hysteresis (32), motors (30b/34), HF (18/19/29). |

## `mcp-server-differential-forms`

_Differential forms / exterior calculus for EM: de Rham complex, cohomology, EM forces theory_

Module: `radia_mcp.differential_forms.server`

| Tool | Description |
|---|---|
| `differential_forms_basics` | Basic differential-form machinery. |
| `differential_forms_bibliography` | Bibliography of the source PDFs distilled into this server. |
| `differential_forms_de_rham` | de Rham complex, Sobolev spaces of forms, Maxwell's house. |
| `differential_forms_em_force_extras` | Additional EM force topics beyond the 7-method catalog (2026-05-22). |
| `differential_forms_em_force_recipe` | PRACTICAL NGSolve recipe for EM force computation (2026-05-22). |
| `differential_forms_feec` | Finite Element Exterior Calculus (Arnold-Falk-Winther 2006). |
| `differential_forms_forces` | Electromagnetic forces in differential-form language. |
| `differential_forms_homology` | Chain complex, homology, Betti numbers, tree-cotree gauge. |
| `differential_forms_kelvin_lab_studies` | Sugahara Lab practical Kelvin-transform case studies (2020-2023). |
| `differential_forms_mathematica_recipes` | Wolfram Language recipes for symbolic verification, pairing |
| `differential_forms_maxwell` | Maxwell's equations in differential-form language. |
| `differential_forms_status` | (no description) |
| `differential_forms_whitney` | Whitney elements: discrete differential forms on a simplicial mesh. |

## `mcp-server-mathematica`

_Mathematica recipes: vector calc, Kelvin transform, symbolic Maxwell, evaluation pipeline_

Module: `radia_mcp.mathematica.server`

| Tool | Description |
|---|---|
| `mathematica_check_identity` | 式 LHS == RHS が常に成り立つかを FullSimplify で判定。 |
| `mathematica_differentiate` | 微分。 |
| `mathematica_evaluate` | Evaluate Wolfram Language code via wolframscript and return result. |
| `mathematica_integrate` | 積分 (定積分 / 不定積分)。 |
| `mathematica_server_status` | (no description) |
| `mathematica_simplify` | Wolfram FullSimplify[expression, assumptions] を実行。 |
| `mathematica_solve` | 方程式 (系) を Solve で解く。 |
| `mathematica_status` | Diagnostic: check wolframscript availability + version + license. |
| `mathematica_to_tex` | Wolfram の TeXForm を文字列で取得 (paper / 数式 DB 登録用)。 |
| `mathematica_unit_convert` | 物理単位変換 (Wolfram Quantity / UnitConvert)。 |
| `mathematica_vector_calc` | ベクトル解析: Curl / Div / Grad / Laplacian / Cross / Dot。 |

## `mcp-server-md2html`

_Markdown -> self-contained HTML with MathJax v3 + tables + fenced code + codehilite + base64 image embed + [N] reference linking + UTF-8 / cp932 fallback for legacy Japanese files.  Promoted from mcp-server-document._

Module: `radia_mcp.md2html.server`

| Tool | Description |
|---|---|
| `md2html_convert` | Convert a Markdown file to a self-contained HTML file. |
| `md2html_status` | (no description) |

## `mcp-server-paper-writing`

_Journal paper writing skill suite: 67 IMRaD / abstract / citation / figure / equation / hedge-counting / passive-voice / paragraph / PDF-edge-overflow / page-limit lint tools + PDF download (IEEE Xplore / ScienceDirect / Emerald with cookies) + Plan-B composite-score modules (reproducibility, statistical reporting, IMRaD discussion, journal fit).  v0.88.0 adds 5 image-based PDF layout verification tools (pymupdf-render thumbnail strip + whitespace anomaly detection + per-page PNG + float-far-from-reference check) and the tex_figure_placement knowledge module (htbp / placeins / widths / anti-patterns / per-journal profiles).  Promoted from mcp-server-document._

Module: `radia_mcp.paper_writing.server`

| Tool | Description |
|---|---|
| `paper_writing_abstract_strength` | Abstract の強度を 4 要素 (problem / method / result-with-number / impact) |
| `paper_writing_acronym_usage_audit` | 略語の使用頻度監査 (grant_writing 実装の re-export)。 |
| `paper_writing_adaptive_health_report` | paper T8 health_report の severity 判定を context で adjust。 |
| `paper_writing_analyze_sentences` | 文長分析 (和文)。journal では長文を避けて読みやすさ重視。 |
| `paper_writing_arxiv_extract_equations` | Extract all displayed equations from a LaTeX source. |
| `paper_writing_arxiv_fetch_latex_source` | Fetch the LaTeX source of an arXiv preprint. |
| `paper_writing_arxiv_search` | Search arXiv via the official Atom XML API. |
| `paper_writing_check_abstract_background_ratio` | Abstract 内で background (導入文) が占める割合を推定。 |
| `paper_writing_check_citation_usage` | TeX 本文中の \cite{} キーと bib file の entries を突合。 |
| `paper_writing_check_english_redflags` | 英文論文の典型的 red flag を検出 (冠詞、時制、自動詞/他動詞 の混同)。 |
| `paper_writing_check_equation_numbering` | 方程式番号 (1), (2), ... の連番欠落 / 重複をチェック。 |
| `paper_writing_check_figure_caption_showing` | Figure caption が showing (describe) 形か telling (claim) 形か判定。 |
| `paper_writing_check_figure_forward_reference` | 図/表の \label と \ref の整合チェック (孤立ラベル / 未解決参照)。 |
| `paper_writing_check_floats_far_from_reference` | Detect figures whose \ref{} appears far from the actual float. |
| `paper_writing_check_imrad_balance` | IMRAD 各セクションの字数バランスを検証する。 |
| `paper_writing_check_kanji_ratio` | 漢字比率 check (本多『日本語の作文技術』第四章 re-export)。 |
| `paper_writing_check_misuse_japanese` | 『問題な日本語』由来の現代誤用検出 (re-export)。 |
| `paper_writing_check_notation_variants` | 和文表記ゆれ検出 (grant_writing 実装の re-export)。 |
| `paper_writing_check_overfull_hbox` | LaTeX ログ中の Overfull \hbox 警告をカウント。journal では許容ゼロ。 |
| `paper_writing_check_paragraph_length` | 段落の字数/語数が適正範囲内か検出。 |
| `paper_writing_check_paragraph_opener` | Introduction/Abstract の段落冒頭が禁断フレーズで始まるか検出。 |
| `paper_writing_check_passive_voice_ratio` | 英文の受動態比率を推定。Wallwork §8: Methods は passive 可、 |
| `paper_writing_check_pdf_advanced_anomalies` | 論文 PDF を実際に読んで、reviewer-2 が刺してくる高頻度の体裁 |
| `paper_writing_check_pdf_edge_overflow` | 論文 PDF を実際に読んで、本文/数式/図表が紙面の端からはみ出して |
| `paper_writing_check_pdf_obvious_errors` | 論文 PDF を実際に開いて「明らかな体裁エラー」を 6 種類スキャン。 |
| `paper_writing_check_prose_density` | Detect compression-induced prose-density anti-patterns. |
| `paper_writing_check_self_citation_ratio` | 自己引用率を算出。Wallwork: <20% 推奨。 |
| `paper_writing_check_sentence_ending_variety` | 文末表現の単調さ / 連続を検出 (中島・塚本 §1.3.2)。 |
| `paper_writing_check_strong_adjective_budget` | 強調副詞/形容詞の過剰使用を検出。 |
| `paper_writing_check_subject_predicate_distance` | 主述直結原則 (本多 p.22) — 「は、」「が、」 主題マーカーの |
| `paper_writing_check_tense_consistency` | Discussion の 3-part tense (hypothesis=現在, result=過去, background=現在完了) 混在検出。 |
| `paper_writing_check_typography_hacks` | Detect typography hacks used to cram content past a page limit. |
| `paper_writing_check_undefined_variables` | Detect math symbols that are used but never defined. |
| `paper_writing_check_word_repetition` | 同一単語の近接障害を検出する (中島・塚本『知的な科学・技術文章の書き方』§1.3.5)。 |
| `paper_writing_citation_health_4_axes` | 論文 reference の 4 軸 health 診断。 |
| `paper_writing_citation_workflow_recipe` | Return the lab's mandatory citation-verification workflow recipe. |
| `paper_writing_claim_quantification` | Unquantified hype claims を検出。 |
| `paper_writing_classify_reviewer_comment` | Reviewer コメントを佐藤 Q40 の 4-tier に自動分類。 |
| `paper_writing_contribution_clarity_score` | Introduction 末尾の Contribution リスト/段落の明確度診断。 |
| `paper_writing_count_underlines` | TeX/LaTeX 内の下線コマンドを実測。論文では原則ゼロを目指す。 |
| `paper_writing_count_weak_expressions` | 弱気修飾語の出現数。journal 論文では conclusion / contribution で |
| `paper_writing_detect_overlapping_text_blocks` | Detect text-on-text overlap (block-vs-block IoU on every page). |
| `paper_writing_detect_page_whitespace_anomalies` | Flag pages that are mostly whitespace (sign of bad float placement). |
| `paper_writing_detect_text_image_overlap` | Detect text blocks overlapping with images on every page. |
| `paper_writing_detect_text_overflow_page` | Detect text blocks extending past the page CropBox / MediaBox. |
| `paper_writing_discussion_structure_4_elements` | Discussion section の 4 要素 (interpretation / limitations / |
| `paper_writing_doi_to_bibtex` | Generate a BibTeX entry for a DOI via Crossref metadata. |
| `paper_writing_em_paper_style` | EM-domain-specific style/notation/convention knowledge. |
| `paper_writing_em_submission_gate` | One-shot EM-paper pre-submission gate. |
| `paper_writing_emerald_download_pdf` | Download an Emerald Publishing PDF (e.g. COMPEL journal). |
| `paper_writing_external_sources_recipe` | Return the GitHub-survey writeup + decision tree + credits. |
| `paper_writing_extract_abstract` | Extract the abstract body from a .tex file. |
| `paper_writing_fetch_and_cite` | One-shot: download IEEE PDF + generate BibTeX entry from a DOI. |
| `paper_writing_figure_referencing_coverage` | Every \label{fig:X}/\label{tab:X} が本文で \ref される回数を集計。 |
| `paper_writing_find_undefined_acronyms` | Latin 略語の初出定義 check (grant_writing 実装の re-export)。 |
| `paper_writing_generate_cover_letter` | 投稿 cover letter の skeleton を生成。 |
| `paper_writing_generate_response_letter` | Reviewer コメントへの応答文テンプレ (佐藤 Q40 "agree first, pivot")。 |
| `paper_writing_given_new_ordering` | Wallwork §3.4-3.6: Given-New 情報配置ルールの heuristic 評価。 |
| `paper_writing_health_report` | paper_writing Plan B の全 T1-T7 を束ねた統合レポート。 |
| `paper_writing_ieee_doi_to_arnumber` | Resolve an IEEE DOI to its IEEE Xplore arnumber. |
| `paper_writing_ieee_download_pdf` | Download an IEEE Xplore PDF via cookie-seeded curl-like session. |
| `paper_writing_journal_fit_assessment` | target journal の aims & scope と論文の fit を診断。 |
| `paper_writing_layout_thumbnail_strip` | Render every page as a thumbnail tile into ONE composite PNG. |
| `paper_writing_layout_visual_recipe` | Return the lab recipe for image-based PDF layout verification. |
| `paper_writing_limitation_statement_presence` | Discussion 内で limitation/caveat 段落の有無・位置・充実度を検査。 |
| `paper_writing_lint_bedrock` | 木下 10 原則 + 本多 + 知的 の bedrock 診断 (re-export)。 |
| `paper_writing_lint_reference_format` | .bib の reference エントリーの完全性と形式を検証。 |
| `paper_writing_next_5_actions` | paper health_report の priority_issues を impact / effort で再 sort、 |
| `paper_writing_pdf_overlap_recipe` | Return the recipe for PDF overlap/overflow detection. |
| `paper_writing_related_work_density` | Introduction 内の \cite 密度・自己引用比率・年度分布を診断。 |
| `paper_writing_render_pages_to_png` | Render PDF pages to PNG files for visual inspection. |
| `paper_writing_reproducibility_open_science_check` | Reproducibility & Open Science の 6 軸を一括診断。 |
| `paper_writing_resolve_doi` | Look up a DOI's metadata via the Crossref public API. |
| `paper_writing_resolve_input_chain` | MCP tool wrapper for resolve_input_chain(). |
| `paper_writing_reviewer_2_trigger_summary` | 悪意ある reviewer-2 が最も突いてくるポイントを weighted union で列挙。 |
| `paper_writing_rewrite_suggest` | paper 用 rewrite candidate 生成 (11 target × 3-5 candidate)。 |
| `paper_writing_root_cause_diagnosis` | health_report 結果を横断 pattern matching し、論文の根本原因を診断。 |
| `paper_writing_run_full_workflow` | paper 用 Phase 1-5 を 1 コール chain 実行。 |
| `paper_writing_sciencedirect_download_pdf` | Download an Elsevier ScienceDirect PDF. |
| `paper_writing_semantic_scholar_citations` | List the papers CITING a given paper. |
| `paper_writing_semantic_scholar_lookup` | Look up a paper via Semantic Scholar API. |
| `paper_writing_semantic_scholar_references` | List the references CITED BY a given paper. |
| `paper_writing_statistical_reporting_compliance` | 各 p 値の周辺で effect size / CI / sample size の有無を診断。 |
| `paper_writing_status` | (no description) |
| `paper_writing_suggest_concept_drops` | Suggest specific concepts to drop when prose is over the |
| `paper_writing_suggest_redundancy_fixes` | 冗長表現 25 パターンの置換候補提示 (re-export)。 |
| `paper_writing_tex_figure_placement` | LaTeX figure placement knowledge: float specifiers, placeins, |
| `paper_writing_title_abstract_conclusion_triangle` | Title / Abstract / Conclusion の三角形整合性を診断。 |
| `paper_writing_usage` | Journal 論文 (IEEE / IEEJ / APS / Elsevier 等) の作文技術ガイド全体。 |
| `paper_writing_validate_abstract_length` | Abstract 字数 / 語数が制限内か検証。言語を自動判定。 |
| `paper_writing_validate_pdf_pages` | PDF のページ数が投稿制限内か検証。pymupdf が必要。 |
| `paper_writing_verify_citation` | Verify a citation BEFORE inserting it into the paper. |

## `mcp-server-graph`

_Sugahara Lab publication-figure style guide: IEEE / IEEJ font/size profiles, MATLAB + Matplotlib snippets, lab style rules (units in parentheses, no in-figure title, Times New Roman serif).  Promoted from mcp-server-document._

Module: `radia_mcp.graph.server`

| Tool | Description |
|---|---|
| `flux_line_recipe` | Lab flux-line tracing + visualization recipes (FEM post-processing). |
| `graph_size_for_target` | Recommend output figure size + font settings for a target embedding. |
| `graph_status` | (no description) |
| `graph_style_guide` | Return the lab-standard graph style guide. |
| `paper_figure_8cm_recipe` | Generate a Python recipe for the 8-cm-column lab anchor (1 or 2 panels). |
| `paper_figure_profiles` | List paper-quality figure profiles + their exact journal geometry. |
| `paper_figure_quality_rules` | Why paper-quality figures need a margin-efficiency gate. |
| `paper_figure_recipe` | Generate a self-contained Python recipe for a paper-quality figure. |

## `mcp-server-chart2d`

_22 paper-quality 2D charts as MCP tools: line / loglog / semilog / step / errorbar / fill_between / bode / histogram / bar / box / violin / ecdf / contour / contourf / pcolormesh / quiver / streamplot / imshow / polar / scatter / phase (Nyquist).  Each accepts return_mode='recipe' (Python text) | 'image' (MCP Image inline) | 'both'.  Inherits radia_mcp.graph profile + gate stack._

Module: `radia_mcp.chart2d.server`

| Tool | Description |
|---|---|
| `chart2d_bar` | Vertical bar chart. |
| `chart2d_bode` | Bode pair: magnitude (top) + phase (bottom), shared frequency axis. |
| `chart2d_box` | Box-and-whisker plot. |
| `chart2d_catalog` | List the 22 2D chart types and their groupings. |
| `chart2d_contour` | Contour lines of Z(X, Y). |
| `chart2d_contourf` | Filled contours.  Default viridis = perceptually uniform + CVD safe. |
| `chart2d_ecdf` | Empirical cumulative distribution. |
| `chart2d_errorbar` | Line + error bars for measured data with uncertainty. |
| `chart2d_fill_between` | Filled band between y and y2 (defaults to 0).  Confidence / hysteresis. |
| `chart2d_histogram` | 1D histogram. |
| `chart2d_imshow` | 2D array as image / heatmap. |
| `chart2d_line` | Standard line plot. y can be 1D or 2D (multiple traces). |
| `chart2d_loglog` | Log-log plot.  Straight line of slope alpha = power-law y=k*x^alpha. |
| `chart2d_pcolormesh` | Irregular-grid heatmap. |
| `chart2d_phase` | Complex-plane (Re vs Im) -- Nyquist / impedance locus / root locus. |
| `chart2d_polar` | Polar / radial plot. |
| `chart2d_quiver` | 2D vector field as arrows. |
| `chart2d_scatter` | Scatter plot of (x, y) points. |
| `chart2d_semilogx` | Linear y, log10 x.  Bode-magnitude / frequency-domain default. |
| `chart2d_semilogy` | Log10 y.  Decaying / growing quantities (relaxation, noise tails). |
| `chart2d_status` | (no description) |
| `chart2d_step` | Stair-step plot.  PWM / discrete-time / step histogram. |
| `chart2d_streamplot` | Field lines of a 2D vector field on a REGULAR grid. |
| `chart2d_violin` | Violin plot (KDE-as-fill). |

## `mcp-server-literature-index`

_★ Meta-MCP: full-text search across 3,889 lab literature files in W:/00_電磁界解析_

Module: `radia_mcp.literature_index.server`

| Tool | Description |
|---|---|
| `literature_build_vector_index` | Build / extend the ChromaDB vector index from PDFs. |
| `literature_by_folder` | List papers in a specific top-level / nested folder. |
| `literature_folder_tree` | List all top-level folders with file count + size. |
| `literature_index_cancel` | Cancel the running indexing job (cooperative; checks every PDF). |
| `literature_index_job_status` | Status of the most recent indexing job (idle / running / |
| `literature_index_status` | (no description) |
| `literature_search` | Search the lab literature corpus by filename keywords. |
| `literature_semantic_search` | Semantic search over indexed PDF text via ChromaDB + sentence- |
| `literature_stats` | Index statistics + cache info. |

## `mcp-server-radia-meta`

_★ RECOMMENDED FIRST CALL. Cross-server catalog of all radia_mcp.* servers — answers "which tool covers concept X?" without trial-and-error._

Module: `radia_mcp.meta.server`

| Tool | Description |
|---|---|
| `radia_mcp_by_tag` | Servers tagged with `tag`. |
| `radia_mcp_get` | Look up one server by short name (e.g. 'optuna', 'ih', 'kelvin'). |
| `radia_mcp_health` | Probe importability of every radia_mcp.* subpackage. |
| `radia_mcp_overview` | Authoritative catalog of all radia_mcp.* servers. |
| `radia_mcp_related` | Servers that pair well with `name` (e.g. radia_mcp_related('optuna') |
| `radia_meta_status` | (no description) |

## `mcp-server-panel-review`

_Radia GUI panel review skill chain (panel-cli-diff / panel-review / panel-qt-test / panel-preview / panel-smoke) + bug catalogue. Surfaces the 13-check list and known panel bug patterns through one queryable tool._

Module: `radia_mcp.panel_review.server`

| Tool | Description |
|---|---|
| `panel_review` | Get Radia GUI panel review skill-chain documentation and bug catalogue. |
| `panel_review_status` | (no description) |
| `panel_review_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 10 topics. |

