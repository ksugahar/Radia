# radia-mcp Tools Inventory

Auto-generated from each server's `mcp.list_tools()` via `scripts/gen_tools_doc.py`. **Do not edit by hand** — regenerate after adding/renaming tools.

Total: **285 tools** across 37 MCP servers.

| Server (console-script) | Subpackage | Tools |
|---|---|---:|
| [`mcp-server-cubit`](#mcp-server-cubit) | `radia_mcp.cubit` | 44 |
| [`mcp-server-build123d`](#mcp-server-build123d) | `radia_mcp.build123d` | 29 |
| [`mcp-server-radia-interop`](#mcp-server-radia-interop) | `radia_mcp.interop` | 8 |
| [`mcp-server-gmsh`](#mcp-server-gmsh) | `radia_mcp.gmsh` | 7 |
| [`mcp-server-radia-ngsolve`](#mcp-server-radia-ngsolve) | `radia_mcp.radia_ngsolve` | 31 |
| [`mcp-server-fem`](#mcp-server-fem) | `radia_mcp.fem` | 9 |
| [`mcp-server-bem`](#mcp-server-bem) | `radia_mcp.bem` | 8 |
| [`mcp-server-matrix-solvers`](#mcp-server-matrix-solvers) | `radia_mcp.matrix_solvers` | 6 |
| [`mcp-server-mor`](#mcp-server-mor) | `radia_mcp.mor` | 4 |
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
| [`mcp-server-optuna`](#mcp-server-optuna) | `radia_mcp.optuna` | 4 |
| [`mcp-server-bayesian-opt`](#mcp-server-bayesian-opt) | `radia_mcp.bayesian_opt` | 3 |
| [`mcp-server-evolutionary`](#mcp-server-evolutionary) | `radia_mcp.evolutionary` | 3 |
| [`mcp-server-mcmc`](#mcp-server-mcmc) | `radia_mcp.mcmc` | 6 |
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

## `mcp-server-fem`

_FEM formulations theory layer (A-Omega / T-Omega / H / Reduced / Darwin, edge / HO / XFEM / IGA / DG, gauging + Kelvin, MSFEM, Schur circuit coupling, NGSolve hierarchical)_

Module: `radia_mcp.fem.server`

| Tool | Description |
|---|---|
| `fem_elements` | Element technology: edge (Nedelec), high-order, XFEM, isogeometric, DG. |
| `fem_gauge_open_boundary` | Gauging + open boundary techniques. |
| `fem_large_scale_special` | Large-scale, error theory, multi-scale (Hollaus MSFEM), misc techniques. |
| `fem_ngsolve_hierarchy` | NGSolve hierarchical H(curl) bases - Zaglmayr / nograds / tree-cotree. |
| `fem_overview` | FEM landscape: lab stack, decision tree, genealogy. |
| `fem_potential_formulations` | Potential formulations: A-Omega, T-Omega, H, Reduced, Darwin. |
| `fem_status` | (no description) |
| `fem_time_domain_axisym` | Time-domain, axisymmetric (Henrotte), harmonic balance, HF, circuit coupling. |
| `fem_xfem_comsol` | XFEM in COMSOL Multiphysics (Jafari-Broumand-Vahab-Khalili 2021). |

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
| `optuna_lab_applications` | Lab applications: how Optuna plugs into Radia / NGSolve work. |
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

## `mcp-server-mcmc`

_MCMC + MCTS + SPM for EM. ★ Hokkaido Sato 2023 MCTS PM motor + Yin 2024 inductor + Saotome 1995 SPM_

Module: `radia_mcp.mcmc.server`

| Tool | Description |
|---|---|
| `mcmc_algorithms` | MCMC algorithm theory: MH / Gibbs / HMC / NUTS / PT, diagnostics. |
| `mcmc_em_applications` | MCMC applied to EM engineering problems. |
| `mcmc_libraries` | MCMC libraries: PyMC / Stan / emcee / NumPyro / BlackJAX. |
| `mcmc_mcts_lab` | Monte Carlo Tree Search for EM design -- lab lineage. |
| `mcmc_spm` | Sampled Pattern Matching (Saotome 1995) -- cosine-similarity |
| `mcmc_status` | (no description) |

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
| `radia_mcp_get` | Look up one server by short name (e.g. 'mcmc', 'ih', 'kelvin'). |
| `radia_mcp_health` | Probe importability of every radia_mcp.* subpackage. |
| `radia_mcp_overview` | Authoritative catalog of all radia_mcp.* servers. |
| `radia_mcp_related` | Servers that pair well with `name` (e.g. radia_mcp_related('mcmc') |
| `radia_meta_status` | (no description) |

## `mcp-server-panel-review`

_Radia GUI panel review skill chain (panel-cli-diff / panel-review / panel-qt-test / panel-preview / panel-smoke) + bug catalogue. Surfaces the 13-check list and known panel bug patterns through one queryable tool._

Module: `radia_mcp.panel_review.server`

| Tool | Description |
|---|---|
| `panel_review` | Get Radia GUI panel review skill-chain documentation and bug catalogue. |
| `panel_review_status` | (no description) |
| `panel_review_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 10 topics. |

