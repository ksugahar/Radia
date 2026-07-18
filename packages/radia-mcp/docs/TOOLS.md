# radia-mcp Tools Inventory

Auto-generated from each server's `mcp.list_tools()` via `scripts/gen_tools_doc.py`. **Do not edit by hand** — regenerate after adding/renaming tools.

Total: **787 tools** across 42 MCP servers.

| Server (console-script) | Subpackage | Tools |
|---|---|---:|
| [`mcp-server-cubit`](#mcp-server-cubit) | `radia_mcp.cubit` | 82 |
| [`mcp-server-build123d`](#mcp-server-build123d) | `radia_mcp.build123d` | 74 |
| [`mcp-server-gmsh`](#mcp-server-gmsh) | `radia_mcp.gmsh` | 12 |
| [`mcp-server-radia-ngsolve`](#mcp-server-radia-ngsolve) | `radia_mcp.radia_ngsolve` | 140 |
| [`mcp-server-radia-matlab`](#mcp-server-radia-matlab) | `radia_mcp.matlab` | 5 |
| [`mcp-server-radia-streamfunction`](#mcp-server-radia-streamfunction) | `radia_mcp.streamfunction` | 3 |
| [`mcp-server-fem`](#mcp-server-fem) | `radia_mcp.fem` | 11 |
| [`mcp-server-bem`](#mcp-server-bem) | `radia_mcp.bem` | 7 |
| [`mcp-server-matrix-solvers`](#mcp-server-matrix-solvers) | `radia_mcp.matrix_solvers` | 6 |
| [`mcp-server-mor`](#mcp-server-mor) | `radia_mcp.mor` | 9 |
| [`mcp-server-ih`](#mcp-server-ih) | `radia_mcp.ih` | 4 |
| [`mcp-server-peec`](#mcp-server-peec) | `radia_mcp.peec` | 4 |
| [`mcp-server-electromagnet`](#mcp-server-electromagnet) | `radia_mcp.electromagnet` | 3 |
| [`mcp-server-motor`](#mcp-server-motor) | `radia_mcp.motor` | 44 |
| [`mcp-server-accelerator`](#mcp-server-accelerator) | `radia_mcp.accelerator` | 4 |
| [`mcp-server-fusion-reactor`](#mcp-server-fusion-reactor) | `radia_mcp.fusion_reactor` | 3 |
| [`mcp-server-magnetic-materials`](#mcp-server-magnetic-materials) | `radia_mcp.magnetic_materials` | 8 |
| [`mcp-server-litz-transmission`](#mcp-server-litz-transmission) | `radia_mcp.litz_transmission` | 4 |
| [`mcp-server-rna-mec`](#mcp-server-rna-mec) | `radia_mcp.rna_mec` | 3 |
| [`mcp-server-topology-optimization`](#mcp-server-topology-optimization) | `radia_mcp.topology_optimization` | 7 |
| [`mcp-server-bayesian-opt`](#mcp-server-bayesian-opt) | `radia_mcp.bayesian_opt` | 3 |
| [`mcp-server-evolutionary`](#mcp-server-evolutionary) | `radia_mcp.evolutionary` | 3 |
| [`mcp-server-data-assimilation`](#mcp-server-data-assimilation) | `radia_mcp.data_assimilation` | 3 |
| [`mcp-server-gnn`](#mcp-server-gnn) | `radia_mcp.gnn` | 3 |
| [`mcp-server-pinn`](#mcp-server-pinn) | `radia_mcp.pinn` | 3 |
| [`mcp-server-pcb`](#mcp-server-pcb) | `radia_mcp.pcb` | 7 |
| [`mcp-server-ndt`](#mcp-server-ndt) | `radia_mcp.ndt` | 4 |
| [`mcp-server-metamaterial`](#mcp-server-metamaterial) | `radia_mcp.metamaterial` | 4 |
| [`mcp-server-nmr-mri`](#mcp-server-nmr-mri) | `radia_mcp.nmr_mri` | 2 |
| [`mcp-server-maglev`](#mcp-server-maglev) | `radia_mcp.maglev` | 4 |
| [`mcp-server-team-benchmark`](#mcp-server-team-benchmark) | `radia_mcp.team_benchmark` | 7 |
| [`mcp-server-differential-forms`](#mcp-server-differential-forms) | `radia_mcp.differential_forms` | 13 |
| [`mcp-server-mathematica`](#mcp-server-mathematica) | `radia_mcp.mathematica` | 11 |
| [`mcp-server-md2html`](#mcp-server-md2html) | `radia_mcp.md2html` | 2 |
| [`mcp-server-chart2d`](#mcp-server-chart2d) | `radia_mcp.chart2d` | 24 |
| [`mcp-server-paper-writing`](#mcp-server-paper-writing) | `radia_mcp.paper_writing` | 179 |
| [`mcp-server-grant-writing`](#mcp-server-grant-writing) | `radia_mcp.grant_writing` | 18 |
| [`mcp-server-poster`](#mcp-server-poster) | `radia_mcp.poster` | 32 |
| [`mcp-server-literature-index`](#mcp-server-literature-index) | `radia_mcp.literature_index` | 9 |
| [`mcp-server-document-meta`](#mcp-server-document-meta) | `radia_mcp.document_meta` | 11 |
| [`mcp-server-radia-meta`](#mcp-server-radia-meta) | `radia_mcp.meta` | 9 |
| [`mcp-server-panel-review`](#mcp-server-panel-review) | `radia_mcp.panel_review` | 3 |

## `mcp-server-cubit`

_Cubit mesh scripting, hex/tet workflow, export formats_

Module: `radia_mcp.cubit.server`

| Tool | Description |
|---|---|
| `cubit_ask` | One-shot search across every Cubit knowledge surface we have. |
| `cubit_ato_levelset_sculpt_source_replay_gate` | Gate official ATO provenance, MBG migration, and headless replay. |
| `cubit_audit_summary` | Return a machine-readable Cubit export-lint audit summary. |
| `cubit_batch_try` | Dry-run a recipe in a fresh headless Cubit subprocess. |
| `cubit_boundary_layer_candidate_gate` | Select a non-inverted boundary-layer sweep candidate with export closure. |
| `cubit_boundary_layer_journal_recovery_gate` | Gate three-parameter, pairwise, headless recovery of a failed journal. |
| `cubit_checkpoint` | Save the current Cubit session state as a named checkpoint. |
| `cubit_conformal_hex_pyramid_tet_interface_gate` | Gate a conformal hex-pyramid-tet interface and independent volume sum. |
| `cubit_cpp_sdk_guide` | Get documentation on building Cubit C++ SDK plugins. |
| `cubit_curate_learned_recipes` | **Lab maintainer tool**: read accumulated `learned_recipes.jsonl`, |
| `cubit_diagnostics_guide` | Get the foundational mesh-diagnostics + cleanup + quality playbook. |
| `cubit_docs` | Get Cubit documentation: export formats, scripting guide, and API reference. |
| `cubit_embedded_pipe_source_recovery_gate` | Gate source-journal replay and semantically classified version recovery. |
| `cubit_embedded_region_mixed_transition_gate` | Gate hex-led tet/pyramid recovery, quality, interfaces, and Gmsh 4.1. |
| `cubit_examples` | Search Cubit journal examples from **multiple unioned sources**. |
| `cubit_examples_refresh` | Force-refresh every Cubit example sub-source. |
| `cubit_exec` | Send arbitrary Cubit commands to the persistent viewer session. |
| `cubit_exec_safely` | Execute commands against the live GUI **with a pre-save and a |
| `cubit_forum_tips` | Get practical Cubit meshing tips sourced from the Coreform forum. |
| `cubit_generate_dialog` | Generate a single complete PySide6 dialog script for a Cubit toolbar |
| `cubit_gmsh_v41_inventory` | Parse inline ASCII Gmsh 4.1 by entity blocks and validate connectivity. |
| `cubit_gmsh_v41_mixed_order_gate` | Gate Gmsh 4.1 mixed topology/order while retaining .vol label authority. |
| `cubit_headless_netgen_export_gate` | Gate migration from a GUI plugin export command to native headless Netgen export. |
| `cubit_helical_conductor_source_gate` | Gate a helical-conductor source replay and classified tet fallback. |
| `cubit_helical_partition_mesh_gate` | Gate a many-body helical mesh against quality and parsed .vol inventory. |
| `cubit_hex_refinement_geometry_gate` | Detect curved-geometry error plateaus in an all-hex refinement series. |
| `cubit_journal_reproducibility_gate` | Compare two Cubit journals without inventing a script root cause. |
| `cubit_levelset_sculpt_hex_validation_gate` | Gate coarse/fine Sculpt all-hex quality and Gmsh volume closure. |
| `cubit_list_checkpoints` | List all saved Cubit checkpoints with size + mtime. |
| `cubit_live_mixed_mesh_gate` | Gate a source-journal hex+pyramid+tet replay from headless Cubit Python. |
| `cubit_loft_high_order_vol_series_gate` | Gate all-hex loft topology, curved payload, sidecars, and quality for orders 1-5. |
| `cubit_lookup` | Retrieve relevant sections from the bundled Cubit knowledge base. |
| `cubit_mapped_boundary_layer_shell_gate` | Gate mapped all-hex boundary layers by nodal shells, quality, and scale. |
| `cubit_mesh_apply_choice` | Apply the human's chosen variant from a previous |
| `cubit_mesh_auto` | Find a working mesh recipe by trying a scheme ladder in batch, |
| `cubit_mesh_carrying_straight_sweep_gate` | Gate a straight include_mesh sweep, topology lattice, and Gmsh export. |
| `cubit_mesh_carrying_straight_sweep_source_replay_gate` | Gate official Help provenance, headless replay, and no-mesh control. |
| `cubit_mesh_diagnose` | Per-volume meshing diagnostic. |
| `cubit_mesh_race` | Race N recipes in parallel batch Cubits, replay the first/best |
| `cubit_mesh_race_review` | Race N AI variants + observe live human, **wait for all** to |
| `cubit_mesh_race_review_async` | **Background launch** of a race review — returns immediately |
| `cubit_mesh_race_smart` | **The flagship workflow**: AI inspects the live Cubit state, |
| `cubit_mesh_race_smart_async` | **Background variant of `cubit_mesh_race_smart`** — AI inspects |
| `cubit_mesh_race_status` | Check the status of a background race launched by |
| `cubit_mesh_race_with_human` | **The radia-mcp signature workflow.** |
| `cubit_mixed_order_series_gate` | Validate mixed-mesh topology and routing across export orders. |
| `cubit_mixed_transition_source_gate` | Gate source commands, headless diagnostics, and quality API recovery. |
| `cubit_partial_volume_hex_diagnosis_gate` | Gate a truthful partial-volume/low-quality hex rejection. |
| `cubit_partitioned_sweep_compatibility_gate` | Gate a legacy webcut/partition journal promoted to an all-hex sweep. |
| `cubit_power_tools_cleanup_source_replay_gate` | Gate an official Power Tools cleanup trace and console diagnosis. |
| `cubit_probe` | Query the Cubit session for geometry/mesh statistics. |
| `cubit_pyramid_degenerate_hex_export_gate` | Gate CPYRAM versus nopyramid decks, including order-2 linearization. |
| `cubit_pyramid_mixed_export_gate` | Gate explicit hex/pyramid/tet preservation in Gmsh and Nastran. |
| `cubit_pyramid_source_plugin_replay_gate` | Gate legacy source migration through executable-owned plugin startup. |
| `cubit_recent_failures` | Return the last N failed Cubit invocations (from persistent log). |
| `cubit_region_owned_mixed_mesh_gate` | Gate region-owned conductor hex and air tet/pyramid topology. |
| `cubit_restore` | Restore a previously-saved Cubit checkpoint by label. |
| `cubit_scaffold_toolbar` | Generate a complete Coreform-Cubit custom-toolbar skeleton on disk. |
| `cubit_session_shutdown` | Stop the persistent Cubit daemon. Next `cubit_show` relaunches. |
| `cubit_session_status` | Return diagnostic info about the Cubit session (pid, alive, bin_dir). |
| `cubit_show` | Load a file into the **persistent Cubit viewer** and optionally run |
| `cubit_snapshot` | Hardcopy the current Cubit view to a PNG file. |
| `cubit_source_journal_replay_gate` | Gate synchronous, headless replay and expected mesh disposition. |
| `cubit_status` | (no description) |
| `cubit_structured_hex_lattice_gate` | Gate structured all-hex counts, quality, and Gmsh volume closure. |
| `cubit_structured_hex_source_replay_gate` | Gate source commands, license diagnostics, and headless exit semantics. |
| `cubit_suggest_next` | Suggest concrete next Cubit commands based on the current state. |
| `cubit_sweep_along_curve_gate` | Gate an all-hex mesh-carrying curve sweep and headless launcher evidence. |
| `cubit_symmetric_swept_mixed_mesh_gate` | Gate symmetric CAD, hex/pyramid/tet topology, quality, and Gmsh closure. |
| `cubit_symmetric_swept_source_replay_gate` | Gate source-journal headless replay and public mixed-mesh closure. |
| `cubit_toolbar_guide` | Get documentation on building custom in-Cubit PySide6 toolbars. |
| `cubit_vol_inventory` | Return semantic element inventory for a Netgen `.vol` export. |
| `cubit_web_docs` | Fetch live Cubit documentation and grep for `query`. |
| `cubit_webcut_conformal_hex_gate` | Gate webcut volume drift, partition balance, interfaces, and hex quality. |
| `cubit_webcut_journal_execution_gate` | Gate source-journal operation order and headless process evidence. |
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
| `build123d_brep_mass_topology_roundtrip_gate` | Gate CAD roundtrips by mass properties, centroid semantics and B-rep Euler topology. |
| `build123d_cad_handoff_manifest` | Run the final build123d CAD handoff manifest gate from JSON inputs. |
| `build123d_cad_route_source_contract` | Gate a build123d CAD package before Cubit hex/mixed route promotion. |
| `build123d_cross_kernel_mass_topology_diagnosis_gate` | Diagnose STEP portability while separating evidence quality from acceptance. |
| `build123d_cubit_quality_ledger_handoff` | Bind build123d CAD rows to a Cubit mesh-quality ledger identity gate. |
| `build123d_cubit_solver_route_handoff` | Bind build123d CAD rows to a Cubit mixed solver-route manifest gate. |
| `build123d_curved_shell_step_semantics_gate` | Diagnose topology-preserving curved STEP mass loss across CAD kernels. |
| `build123d_curved_step_topology_crosscheck_gate` | Gate curved STEP mass and exact topology across independent imports. |
| `build123d_discussions` | Search the build123d GitHub Issues archive (de-facto forum). |
| `build123d_drafted_housing_cross_kernel_gate` | Gate drafted housing mass/topology across B-rep, STEP, Cubit, and Gmsh. |
| `build123d_drafted_housing_source_replay_gate` | Gate tagged draft/fillet/hole source and headless mesh-companion replay. |
| `build123d_dual_api_perforated_board_gate` | Gate equivalent Builder/Algebra perforated boards through two CAD imports. |
| `build123d_dual_api_prismatic_pattern_gate` | Gate native dual-API parity separately from external STEP-kernel bias. |
| `build123d_dual_api_source_replay_gate` | Gate immutable upstream dual-API execution and headless CAD replay. |
| `build123d_examples` | Search build123d + **bd_warehouse** + **GitHub Issues** (unioned). |
| `build123d_examples_refresh` | Force-refresh all build123d sources from GitHub + YouTube. |
| `build123d_external_cad_mass_topology_gate` | Crosscheck two CAD kernels without confusing entity centers with mass centroids. |
| `build123d_external_cad_volume_evidence_package` | Bundle dual-source CAD volume evidence before reuse. |
| `build123d_faceted_edit_portability_gate` | Separate faceted CAD portability from downstream mesh readiness. |
| `build123d_faceted_source_replay_gate` | Gate tagged source, dependent STL, viewer stub, and external replay. |
| `build123d_heal` | Run OCCT `ShapeFix_Shape` on a STEP file, write a healed copy. |
| `build123d_heat_exchanger_source_recovery_gate` | Gate the upstream heat-exchanger replay and rotation recovery. |
| `build123d_inspect_step` | Inspect an external STEP file via the build123d / OCCT importer |
| `build123d_jointed_assembly_heal_invariance_gate` | Verify that STEP solid closure loss persists across heal/noheal imports. |
| `build123d_jointed_assembly_source_replay_gate` | Gate immutable source, joint graph, and headless external-CAD evidence. |
| `build123d_jointed_assembly_step_closure_gate` | Diagnose a component-level solid closure loss in a jointed STEP assembly. |
| `build123d_loft_example_source_replay_gate` | Gate the immutable upstream loft source and headless CAD replay. |
| `build123d_lofted_shell_handoff_gate` | Gate a bounded lofted-shell CAD handoff without solver overclaim. |
| `build123d_lookup` | Retrieve relevant sections from the bundled build123d knowledge + |
| `build123d_mass_property_crosscheck` | Compare build123d volume/area/bbox rows with one or more CAD sources. |
| `build123d_motor_housing_thermal_reference` | Return analytic volume/area/body-count data for a finned motor housing. |
| `build123d_nested_assembly_volume_gate` | Distinguish a zero parent Compound from an empty CAD handoff. |
| `build123d_path_sweep_handoff_gate` | Gate a curved build123d sweep through analytic path and STEP/CAD checks. |
| `build123d_path_sweep_source_contract_gate` | Gate the source-native build123d sweep idiom and ``is_valid`` API form. |
| `build123d_patterned_compound_translation_gate` | Diagnose dominant curved-body STEP bias without solver-ready overclaim. |
| `build123d_perforated_prism_roundtrip_gate` | Check STEP roundtrip volume and through-hole boundary topology. |
| `build123d_platonic_solid_family_gate` | Gate all five Platonic solids by topology, analytic volume and CAD replay. |
| `build123d_recent_failures` | Return the last N failed `execute_build123d` invocations (from log). |
| `build123d_reflection_rotation_handoff_gate` | Gate reflection failures and a proper-rotation two-body STEP handoff. |
| `build123d_repeated_cavity_dual_api_gate` | Gate dual APIs and four STEP imports for a repeated-feature cavity solid. |
| `build123d_repeated_cavity_source_replay_gate` | Gate immutable dual sources, STEP identities, and headless CAD replay. |
| `build123d_status` | (no description) |
| `build123d_step_portability_diagnosis_gate` | Diagnose whether STEP mass loss occurs in export or external import. |
| `build123d_stud_wall_source_replay_gate` | Gate exact stud-wall source, RigidJoints, and headless CAD replay. |
| `build123d_suggest_next` | Suggest concrete next build123d steps toward a goal. |
| `build123d_tea_cup_source_contract_gate` | Gate the upstream tea-cup source and headless portability diagnosis. |
| `build123d_to_cubit_hex` | End-to-end: build123d script → STEP → `cubit_mesh_auto` (batch- |
| `build123d_try` | Dry-run a build123d script in a **fresh Python subprocess**. |
| `build123d_try_race` | Race N build123d script variants in parallel subprocesses, |
| `build123d_upstream_example_roundtrip_gate` | Gate source identity and STEP self-roundtrip for an upstream build123d example. |
| `build123d_upstream_source_external_cad_contract_gate` | Gate immutable upstream execution and an explicit external-CAD decision. |
| `build123d_usage` | Get build123d CAD modeling documentation for CAE workflows. |
| `build123d_vase_external_solid_contract_gate` | Gate an exact vase replay and reject zero-volume external solids. |
| `build123d_volume_crosscheck` | Compare build123d reference volumes with Cubit or external-CAD volumes. |
| `build123d_volume_crosscheck_source_coverage_gate` | Require Cubit/external CAD source coverage after a volume crosscheck. |
| `build123d_volume_crosscheck_source_identity_gate` | Require source identity metadata after an external CAD volume crosscheck. |
| `build123d_volume_crosscheck_with_units` | Normalize explicit cubic units before comparing CAD volumes. |
| `build123d_web_docs` | Fetch live build123d documentation (readthedocs) and grep for `query`. |
| `build123d_wrap_faces_rotational_source_replay_gate` | Gate immutable wrap_faces, thicken, and rotational-pattern replay. |
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

## `mcp-server-gmsh`

_GMSH MSH v4.1 inspect/validate/convert/write_node_data_

Module: `radia_mcp.gmsh.server`

| Tool | Description |
|---|---|
| `get_gmsh_lint_rules` | List all available GMSH lint rules with descriptions. |
| `gmsh_audit_summary` | Return a machine-readable GMSH lint audit summary. |
| `gmsh_examples` | Get GMSH tutorial and example documentation. |
| `gmsh_mesh_generation_remediation_plan` | List scripts that still use GMSH as a mesh generator. |
| `gmsh_numsubedges_remediation_plan` | List scripts that need high-order GMSH display settings. |
| `gmsh_post_display_contract` | Return the shared .geo/.geo.opt/.msh.opt contract for Gmsh post artifacts. |
| `gmsh_post_display_gate` | Validate a Radia/Gypsilab Gmsh post-display manifest. |
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
| `acoustic_duct_band_gap_gate` | Gate a confined acoustic band gap against empty and free-space controls. |
| `acoustic_fembem_cross_learnings` | Method-selection and validation cross-learnings for radia-acoustic, distilled |
| `adjoint_gradient_scaling_gate` | Gate reverse-mode solve scaling, FD agreement and ascent direction. |
| `airgap_motor_workflow` | Get AGE rotating machine workflow documentation -- nonlinear iron + AGE coupling. |
| `alternate_eddy_loss_formulation_gate` | Gate volume-resolved and surface-impedance losses as non-additive alternatives. |
| `analytical_formulas` | Get documentation for radia.analytical_formulas (closed-form reference layer). |
| `autodiff_harmonic_balance_convergence_gate` | Gate AD harmonic balance without mean-only false convergence. |
| `axifem_documentation` | Get radia-core axifem documentation: Henrotte axisymmetric FE |
| `balanced_mcp_learning_profile` | Return the ten-stage equal public/source MCP learning contract. |
| `basis_functions` | Finite-element basis function library — Mathematica-canonical |
| `bem_cln` | Get BEM-CLN (per-element multipole CLN with Schur-F termination) |
| `cln_3d` | Get 3D Cauer Ladder Network (CLN) / Kameari-Tanimoto iteration |
| `cln_3d_notebook` | Retrieve Tanimoto's raw 3D CLN notebook Python code. |
| `cln_sibc_orthogonal` | Get CLN expansion-point + SIBC orthogonal-residual theory documentation. |
| `cln_sphere_dd_pipeline` | Get the Sphere DD (double-double, ~32 digit) VIM Cauer Ladder Network |
| `cogging_torque_periodicity_gate` | Gate a zero-current torque sweep over one slot/pole LCM period. |
| `coil_self_resonance_sweep_gate` | Gate complex coil impedance, self-resonance, and sweep replay. |
| `complex_vector_field_maximum_gate` | Gate complex vector-field magnitudes and per-material maxima. |
| `conductive_network_resistance_monotonicity_gate` | Gate Rayleigh resistance monotonicity for conductive contact networks. |
| `coupled_cq_refinement_gate` | Gate coupled FEM/BEM CQ symbols, contour balance, and refinement. |
| `cq_response_reality_gate` | Gate a coupled CQ solve, including its real time-domain reconstruction. |
| `cq_scattering_arrival_gate` | Gate CQ scattered-field causality against a geometric ray arrival. |
| `cyclic_terminal_phasor_balance_gate` | Gate cyclic voltage/current triplets and all-terminal phasor KCL. |
| `cyclic_terminal_source_sweep_gate` | Gate cyclic terminal charges without assuming formulations are identical. |
| `cylindrical_conductor_skin_bessel_gate` | Gate cylindrical skin-effect identities and exact Bessel structure. |
| `dtn_coarse_mesh` | Why open-boundary methods stay accurate on COARSE meshes -- a spectral |
| `dual_formulation_force_error_convergence_gate` | Gate force-error convergence envelopes across two or more formulations. |
| `dual_formulation_symmetric_field_profile_gate` | Gate full-profile agreement and symmetry for two field formulations. |
| `energy_budgeted_trace_kkt_gate` | Gate KKT closure for an energy-budgeted FEM/BEM trace fit. |
| `esim` | Get ESIM (Effective Surface Impedance Method) general documentation. |
| `fem_bem_capstone_suite_gate` | Gate a ten-case first-order FEM/BEM reference capstone suite. |
| `fem_bem_schur` | Get FEM-BEM Schur coupling documentation -- exact open boundary for interior FEM. |
| `femm_parity_documentation` | Get FEMM-parity documentation: which FEMM (Finite Element Method Magnetics, |
| `finite_solenoid_surface_current_gate` | Gate a finite-solenoid surface-current profile and signed linearity. |
| `force_coenergy_displacement_gate` | Gate direct force against the central derivative of magnetic coenergy. |
| `force_position_profile_gate` | Gate a force-position sweep without assuming it is monotonic. |
| `force_validation` | EM force extraction in NGSolve + independent <-> NGSolve cross-validation. |
| `fsi_scattering_invariants_gate` | Gate lossless FSI reciprocity, energy closure, and exterior-method agreement. |
| `get_radia_lint_rules` | List all available NGSolve lint rules with descriptions. |
| `global_local_optimization_replay_gate` | Gate a stochastic global-search to derivative-checked local-polish replay. |
| `gmsh_post_spec` | GMSH post-processing specification for Radia panels. |
| `grounded_sphere_capacitance_convergence_gate` | Gate grounded-sphere image-series convergence and mixed-boundary energy. |
| `hall_effect_transverse_voltage_gate` | Gate Hall voltage by coefficient, drive, field, and replay controls. |
| `harmonic_current_port_power_energy_identity_gate` | Gate peak-phasor port, loss, energy, flux, and profile identities. |
| `harmonic_magnetic_force_triplet_closure_gate` | Gate harmonic body-force methods and source/body action-reaction closure. |
| `harmonic_zero_net_circuit_gate` | Gate zero-net harmonic phasors, Faraday sign, loss, and force metadata. |
| `hartmann_profile_gate` | Gate a Hartmann-number sweep against an independent channel profile. |
| `hdiv_vim` | HDiv-type VIM (Volume Integral Method) demag operator -- the lab's FEEC H(div) RT |
| `helmholtz_double_layer_low_frequency_gate` | Gate the quadratic low-frequency correction of a Helmholtz double layer. |
| `helmholtz_dual_formulation_axis_gate` | Gate Helmholtz-coil axis symmetry, flatness, and formulation agreement. |
| `heterogeneous_current_flow_p1_reintegration_gate` | Gate heterogeneous current-flow P1 reintegration and sign covariance. |
| `heterogeneous_part_mesh_replay_gate` | Diagnose deterministic heterogeneous part-mesh replay drift. |
| `hmatrix_compression_scaling_gate` | Gate H-matrix accuracy, bounded rank, and subquadratic storage scaling. |
| `homogenized_bundle_impedance_comparison_gate` | Gate a stranded-bundle approximation against an explicit reference. |
| `hysteresis_minor_loop_replay_gate` | Gate history, knot normalization, signed loss, and exact loop replay. |
| `inductance_energy_mutual_gate` | Gate L=2W/I^2 and an analytic one-direction mutual inductance. |
| `inductance_matrix_family_gate` | Gate two-winding matrices, identities, replay, and turn scaling. |
| `install_deploy` | Radia install / deploy policy and recipes — 2-tier configuration |
| `kelvin_identify_post_hoc` | Add Kelvin Periodic Identifications to an existing NGSolve mesh |
| `kelvin_transformation` | Get Kelvin transformation documentation for open boundary FEM problems. |
| `leakage_inductance_closure_gate` | Gate compensated-energy and unit-current-matrix leakage inductance. |
| `linear_axisymmetric_circuit_energy_gate` | Gate current, flux, field, and energy identities on one fixed mesh. |
| `linear_eddy_levitation_force_gate` | Gate linear harmonic levitation force by dual extraction and I-squared laws. |
| `linear_induction_frequency_sweep_gate` | Gate a linear-induction frequency sweep by thrust, loss, and phase balance. |
| `linear_magnetization_scaling_gate` | Gate source scaling plus an independent refined P1 FEM reference. |
| `linear_sphere_geometry_convergence_gate` | Gate first-order sphere tri/tet geometry convergence and replay. |
| `linked_study_silent_noop_gate` | Verify a linked native run that returned without creating solver results. |
| `lint_radia_directory` | Lint all Python scripts in a directory for NGSolve convention violations. |
| `lint_radia_script` | Lint a single Python script for Radia + NGSolve convention violations. |
| `loop_learning` | Public-safe CAE loop learning rules distilled from repeated validation |
| `loss_temperature_coupling_gate` | Gate an electromagnetic-loss to transient-temperature handoff. |
| `lossy_dielectric_complex_power_refinement_gate` | Gate lossy-dielectric constitutive, energy, complex-power, and mesh closure. |
| `magnetic_conductive_shield_frequency_gate` | Gate low-frequency magnetic loading and high-frequency eddy shielding. |
| `magnetic_force_method_profile_gate` | Gate magnetic-force profiles with explicit body/surface selection scope. |
| `magnetostatic_open_boundary_equivalence_gate` | Gate gauge-invariant equivalence of two magnetostatic open-boundary solutions. |
| `manual_auto_mixed_mesh_preservation_gate` | Gate exact manual-region preservation and bounded automatic remeshing. |
| `material_contrast_force_gate` | Gate null, attraction, and increasing-repulsion material-force cases. |
| `matlab_acoustic_fembem_agent_guide` | Agent guide for MATLAB acoustic FEM-BEM / Gypsilab-style workflows. |
| `md2html_usage` | Get md2html converter documentation (MathJax, reference links, styled HTML). |
| `motion_coupled_eddy_levitation_transient_gate` | Gate motion-coupled lift while detecting aliased force output times. |
| `moving_conductor_eddy_brake_gate` | Gate motion, Lorentz-force, and Joule-loss table identities. |
| `multiconductor_capacitance_cross_formulation_gate` | Gate N-conductor Maxwell matrices across volume and boundary formulations. |
| `multiport_impedance_sweep_gate` | Gate common-grid, positive-real, nontrivial complex impedance sweeps. |
| `ngsbem_inductance` | Get ngsolve.bem boundary element method documentation for inductance extraction. |
| `ngsolve_usage` | Get NGSolve finite element library usage documentation. |
| `nonlinear_actuator_saturation_knee_gate` | Gate an axisymmetric nonlinear actuator by a shared L/F saturation knee. |
| `nonlinear_bh_piecewise_material_gate` | Gate secant and left-interval differential permeability from B-H rows. |
| `nonlinear_inductance_sweep_gate` | Gate nonlinear apparent/incremental matrices, duality, and replay. |
| `one_port_power_balance_sweep_gate` | Gate passive one-port accepted power against S11 and reference impedance. |
| `one_port_vi_s_impedance_gate` | Gate one-port S, V/I, impedance-transform, and power identities. |
| `opposed_busbar_skin_force_gate` | Gate AC skin/proximity, phasor identities, and Lorentz action-reaction. |
| `panel_add_param` | Plan where to add a new parameter to a Radia-NGSolve panel. |
| `panel_describe_jp` | 現在のパネルソースを AST 解析して日本語で詳細に説明する。 |
| `panel_gui_pitfalls` | Pitfalls and lessons learned from Radia GUI / Cubit panel development. |
| `panel_schema` | Show Radia-NGSolve panel definitions with Japanese labels and physics. |
| `panel_widget_locations` | Return file:line locations for everything that touches a widget. |
| `parallel_wire_force_refinement_gate` | Gate a reciprocal two-wire force refinement sweep without requiring monotone error. |
| `passive_axial_bearing_stiffness_gate` | Gate signed force, action-reaction, axial stability, and sweep replay. |
| `peec_inductance` | Get documentation for the Radia PEEC-inductance (coil only, STEP) panel mode. |
| `periodic_unwrapped_pm_machine_replay_gate` | Gate topology-aware PM-machine field symmetry and replay stability. |
| `permanent_magnet_recoil_state_gate` | Gate nonlinear, open-circuit, and partial-recoil PM field states. |
| `physics_result_preflight_gate` | Gate physics namespace, selection, solution, and license metadata before result evaluation. |
| `pwm_controlled_motor_loss_gate` | Gate PWM current-control and aggregate/harmonic loss-table identities. |
| `radar_range_angle_localization_gate` | Gate wideband range-angle localization of multiple targets. |
| `radar_range_rcs_profile_gate` | Gate wideband range-RCS localization, method agreement, and analytic amplitude. |
| `radia_ngsolve_status` | (no description) |
| `radia_usage` | Get Radia C++ library usage documentation. |
| `radial_bearing_force_symmetry_gate` | Gate magnetic-body force with equal and mirrored excitation controls. |
| `reciprocal_two_port_power_sweep_gate` | Gate complex two-port reciprocity, symmetry, passivity, and power closure. |
| `regularized_trace_inverse_path_gate` | Gate a P1 trace Tikhonov path, L-curve, Morozov, and replay. |
| `release_workflow` | Release-QUD workflow for the Radia monorepo |
| `rf_sweep_artifact_summary_gate` | Gate a solved two-port sweep artifact and its process-neutral metadata. |
| `rotating_conductor_transient_gate` | Gate moving-axis migration, rotational kinematics, and loss partition. |
| `rotational_eddy_brake_energy_gate` | Gate free rotational braking with angular impulse and field energy. |
| `rotational_kinematics_time_axis_gate` | Gate a result-table time axis using angle/speed kinematics. |
| `rwg_hcurl_trace_consistency_gate` | Gate RWG/HCurl trace topology, de Rham closure, and reference matrices. |
| `single_loop_source_normalized_field_gate` | Gate a single-loop field transfer across two port formulations. |
| `skin_effect_adaptive_energy_loss_gate` | Gate current-port identities and adaptive skin-effect loss convergence. |
| `source_free_static_null_solution_gate` | Gate a source-free static Maxwell solve against the exact zero solution. |
| `source_off_linear_relaxation_gate` | Gate a linear source-off RL relaxation using total current and field decay. |
| `sparsesolv` | Get sparsesolv documentation and code examples. |
| `standalone_panels` | Retired standalone PySide panel topic.  The canonical Radia panel surface |
| `static_field_shim_family_gate` | Gate static-field scale, ROI uniformity, shim sensitivity, and map quality. |
| `symmetric_axial_field_profile_gate` | Gate an origin-centered axial profile by analytic value and symmetry. |
| `symmetric_complex_field_curve_gate` | Gate an even- or odd-sampled complex field curve by mirror symmetry. |
| `taskmanager` | NGSolve TaskManager parallelism — usage, MKL interaction, audit, C++. |
| `technical_reports` | Return implementation-oriented knowledge distilled from IEEJ reports. |
| `thermal_robin_boundary_balance_gate` | Gate signed Robin heat balance, mesh plateau, replay, and reflection. |
| `three_phase_winding_power_balance_gate` | Gate three-phase balance, STAR KCL, and coupled-winding copper power. |
| `transient_conductor_replay_identity_gate` | Gate full transient conductor histories, identities, and independent replay. |
| `transient_coupled_coil_response_gate` | Gate a passive shorted-secondary transient induced-current history. |
| `twin_conductor_skin_effect_frequency_gate` | Gate passive twin-conductor R/L and impedance trends over frequency. |
| `two_body_force_magnitude_replay_gate` | Gate unsigned two-body force balance and two fresh solver replays. |
| `two_conductor_capacitance_identity_gate` | Gate two-conductor capacitance using terminal charge and field energy. |
| `two_conductor_capacitance_matrix_gate` | Gate reciprocal Maxwell and mutual capacitance matrix representations. |
| `two_terminal_dc_conduction_power_gate` | Gate current closure, Joule power, adaptive convergence, and replay. |
| `two_winding_frequency_faraday_gate` | Gate two-winding complex response against linked-flux Faraday identity. |
| `urn` | Universal Relaxation Network (URN): causal/passive rational fitting of a |
| `urn_fit` | Fit a complex frequency response with a Universal Relaxation Network and |
| `voice_coil_force_flux_sweep_gate` | Gate a PM voice-coil current sweep by force, flux, symmetry, and mesh evidence. |

## `mcp-server-radia-matlab`

_Official MATLAB MCP composition and generic ML/RL gates_

Module: `radia_mcp.matlab.server`

| Tool | Description |
|---|---|
| `matlab_agent_guide` | (no description) |
| `matlab_extension_contract` | (no description) |
| `matlab_official_server_config` | (no description) |
| `matlab_radia_acoustic_interface_contract` | (no description) |
| `radia_matlab_status` | (no description) |

## `mcp-server-radia-streamfunction`

_Stream-function coil design: (ACA+)+TSVD least-norm, FE-direct psi, regularisation / folded-Tikhonov Pareto front, single-stroke chain, sheet-metal levers_

Module: `radia_mcp.streamfunction.server`

| Tool | Description |
|---|---|
| `radia_streamfunction_status` | (no description) |
| `radia_streamfunction_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 25 topics. |
| `streamfunction` | Get Stream-Function (SF) coil-design documentation. |

## `mcp-server-fem`

_FEM formulations theory layer (A-Omega / T-Omega / H / Reduced / Darwin, edge / HO / XFEM / IGA / DG, gauging + Kelvin, MSFEM, Schur circuit coupling, NGSolve hierarchical)_

Module: `radia_mcp.fem.server`

| Tool | Description |
|---|---|
| `fem_elements` | Element technology: edge (Nedelec), high-order, XFEM, isogeometric, DG. |
| `fem_equivalence_source` | Equivalence-theorem near-field source (Schelkunoff/Love -- Stratton-Chu). |
| `fem_gauge_open_boundary` | Gauging + open boundary techniques. |
| `fem_large_scale_special` | Large-scale, error theory, multi-scale (Hollaus MSFEM), misc techniques. |
| `fem_ngsolve_hierarchy` | NGSolve hierarchical H(curl) bases - Zaglmayr / nograds / tree-cotree. |
| `fem_nonconforming_mesh_coupling` | Non-conforming mesh coupling: mortar / Nitsche / FETI-DP / BDDC / DG / |
| `fem_overview` | FEM landscape: lab stack, decision tree, genealogy. |
| `fem_potential_formulations` | Potential formulations: A-Omega, T-Omega, H, Reduced, Darwin. |
| `fem_status` | (no description) |
| `fem_time_domain_axisym` | Time-domain, axisymmetric (Henrotte), harmonic balance, HF, circuit coupling. |
| `fem_xfem_em_hiruma` | EM-XFEM (Hiruma 2023): electromagnetic XFEM for eddy-current |

## `mcp-server-bem`

_MoM/BEM theory: RWG, EFIE/MFIE/CFIE/PMCHWT, Loop-Star, Calderon, Radia HDiv-VIM, HACApK, FEM-BEM_

Module: `radia_mcp.bem.server`

| Tool | Description |
|---|---|
| `bem_fem_bem_hybrid` | FEM-BEM hybrid methods for open-boundary EM. |
| `bem_h_matrix` | H-matrix / ACA acceleration for BEM. |
| `bem_low_freq` | Low-frequency BEM stabilization. |
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
| `electromagnet_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 14 topics. |
| `electromagnet_usage` | Get accelerator electromagnet analysis documentation. |

## `mcp-server-motor`

_Motor analysis: ONELAB transient, Hollaus effective material (lamination), Wakao autoencoder topology, Kaimori-Mifune Darwin TD_

Module: `radia_mcp.motor.server`

| Tool | Description |
|---|---|
| `motor_age_quality` | NGSolve AGE quality gates for radia-motor. |
| `motor_age_validation_plan` | Route a motor prompt to the required NGSolve AGE quality gates. |
| `motor_angle_periodic_rom` | HCurl Eddy Bubble + HDiv-MMM angle-periodic motor ROM knowledge. |
| `motor_bibliography` | Search the motor analysis bibliography catalog. |
| `motor_darwin_model` | Darwin-model time-domain formulation (capacitive + inductive coupling). |
| `motor_deck_bridge` | Public-safe motor deck corpus bridge for radia-motor. |
| `motor_dual_lane_training_catalog` | Return the public-safe wide motor learning catalog. |
| `motor_dual_lane_training_gate` | Check that the public motor catalog is complete and provenance-scrubbed. |
| `motor_dual_lane_training_route` | Route a motor prompt to one catalog case and both radia-motor lanes. |
| `motor_dual_torque_method_curve_gate` | Gate two independently evaluated static-torque curves. |
| `motor_electrothermal_result_chain_gate` | Gate a four-stage motor electrothermal result handoff. |
| `motor_em_force_extras` | Forward to `differential_forms_em_force_extras` -- advanced EM force |
| `motor_em_force_recipe` | Practical NGSolve EM-force recipe for motor analysis. |
| `motor_femm_transient` | FEMM newbuild transient solver — Lange-Henrotte-Hameyer 2009 |
| `motor_field_quick_check` | First-order 2D magnetic-circuit/BEM-like motor quick check. |
| `motor_force_report_method_metadata_gate` | Gate a force report using independent methods and action-reaction. |
| `motor_force_rotation_covariance_gate` | Check that a planar force vector follows a rotated excitation/geometry. |
| `motor_henrotte_lineage` | The Henrotte–Hameyer–RWTH research arc (energy-consistent E&M FE). |
| `motor_hollaus_eddy` | Karl Hollaus / TU Wien MSFEM for laminated-iron eddy currents. |
| `motor_hollaus_genealogy` | Visualize the Karl Hollaus / TU Wien MSFEM research genealogy |
| `motor_ipm_two_run_ldlq_gate` | Gate same-angle PM-only/current-on runs and extract ``Ld``/``Lq``. |
| `motor_magnet_model_handoff_gate` | Gate a converged source result and two-file downstream magnet model. |
| `motor_mirror_symmetric_three_magnet_handoff_gate` | Gate grouped magnetization vectors, mirror symmetry, and fresh replay. |
| `motor_motion_table_coordinate_gate` | Validate independent 3D translation and rotation motion tables. |
| `motor_onelab` | ONELAB/GetDP electric-machine reference template knowledge. |
| `motor_periodic_torque_sampling_gate` | Validate periodic torque sampling and FFT endpoint ownership. |
| `motor_permanent_magnet_demagnetization_history_gate` | Gate irreversible PM state across one history or a replayed case family. |
| `motor_permanent_magnet_force_pair_gate` | Gate attraction/repulsion reversal for a facing permanent-magnet pair. |
| `motor_phase_flux_park_alignment_gate` | Gate a PM-only three-phase flux sweep in the rotating d/q frame. |
| `motor_planar_coupling` | 2D PLANAR machine modelling in radia: HDiv-VIM soft-iron demag + the shared |
| `motor_rotating_circuit_transient_gate` | Gate rotating-circuit identities and endpoint state before FFT use. |
| `motor_status` | (no description) |
| `motor_thermal_handoff_gate` | Validate one motor-loss table for both LPTN and 3D all-hex thermal paths. |
| `motor_topology_optimization` | SynRM topology optimization (Wakao 2025 autoencoder + level-set). |
| `motor_transient_no_load_load_cycle_gate` | Gate paired no-load and loaded three-phase transient cycles. |
| `motor_triple_check_artifact_gate` | Validate a combined AGE and HDiv-MMM/HCurl eddy-bubble motor artifact. |
| `motor_triple_check_plan` | Plan the standard radia-motor comparison. |
| `motor_tritool_cross_reference` | Tri-tool cross-reference: FEMM / JMAG / radia-ngsolve (相互学習). |
| `motor_validation_artifact_gate` | Check whether a motor cross-validation artifact can train radia-motor. |
| `motor_validation_lane_template` | Return the JSON artifact template for a motor validation lane. |
| `motor_validation_lanes` | Cross-validation lane policy for radia-motor. |
| `motor_validation_router` | Route a motor prompt to a public deck, field quick check, and NGSolve AGE validation. |
| `motor_variable_magnet_material_gate` | Gate variable-PM material parameters and their authoritative source. |
| `motor_virtual_work_width_ladder_gate` | Select a coenergy-difference angle width against independent torque. |

## `mcp-server-accelerator`

_Accelerator physics: beam optics, dipole/quad/sext magnets, undulator/wiggler_

Module: `radia_mcp.accelerator.server`

| Tool | Description |
|---|---|
| `accelerator` | Accelerator magnet design with Radia + radia-mcp. |
| `accelerator_magnetic_trajectory_pair_gate` | Gate paired charged-particle trajectories with magnetic field off/on. |
| `accelerator_status` | (no description) |
| `accelerator_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 17 topics. |

## `mcp-server-fusion-reactor`

_Fusion reactor magnets: tokamak ITER + stellarator LHD/W7-X/heliotron lineage_

Module: `radia_mcp.fusion_reactor.server`

| Tool | Description |
|---|---|
| `fusion_reactor` | Fusion reactor magnet knowledge. |
| `fusion_reactor_status` | (no description) |
| `fusion_reactor_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 12 topics. |

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
| `periodic_hysteresis_loss_energy_gate` | Gate periodic hysteresis power by cycle energy and loss closure. |

## `mcp-server-litz-transmission`

_Litz wire AC loss (Dowell, homogenization, magnetic-plated wire) + multiconductor transmission line theory_

Module: `radia_mcp.litz_transmission.server`

| Tool | Description |
|---|---|
| `litz_proximity_approximation_pair_gate` | Validate a reduced proximity-effect bundle against an explicit model. |
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
| `topology_opt_cae_ai_artifact_gate` | Gate CAE-AI artifacts before they are promoted as engineering results. |
| `topology_opt_nonlinear_lsq_multistart_gate` | Gate nonlinear least-squares multistart, Jacobian, and KKT evidence. |
| `topology_opt_shape_optimization` | Shape optimization for nonlinear magnetostatics. |
| `topology_opt_simplex_stationarity_audit_gate` | Audit derivative-free convergence using independent stationarity checks. |
| `topology_opt_topology_derivative` | Topological derivative for changing topology (adding/removing material). |
| `topology_optimization_status` | (no description) |

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

## `mcp-server-pcb`

_Wireless Power Transfer: coil + compensation (SS/LCC/LCL), efficiency, IEC 61980 / SAE J2954, FOD, dynamic EV / robot / bearingless motor, capacitive / microwave / metamaterial_

Module: `radia_mcp.pcb.server`

| Tool | Description |
|---|---|
| `pcb_alternatives` | Alternative WPT: capacitive, microwave/rectenna, metamaterial. |
| `pcb_applications` | WPT applications: dynamic EV, robot, bearingless motor. |
| `pcb_coil_compensation` | Coil design + compensation topology + resonance matching. |
| `pcb_efficiency_safety` | Efficiency (Q, k, kQ) + safety + IEC/SAE standards. |
| `pcb_fod` | ★ Foreign Object Detection (FOD) — lab core research. |
| `pcb_overview` | WPT landscape: regimes, decision tree, lab focus. |
| `pcb_status` | (no description) |

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
| `metamaterial_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 13 topics. |

## `mcp-server-nmr-mri`

_NMR/MRI: gradient coils, B0 shimming, RF coils, field uniformity_

Module: `radia_mcp.nmr_mri.server`

| Tool | Description |
|---|---|
| `nmr_mri_bibliography` | Search the NMR/MRI bibliography catalog. |
| `nmr_mri_status` | (no description) |

## `mcp-server-maglev`

_Magnetic levitation, UNIFIED: maglev systems (EMS/EDS/PM/SC/Halbach) + levitation FORCE physics (induction/EML/AMB/superconducting/diamagnetic/Earnshaw/force-computation). Lab research: Radia IEM<->FEM weak coupling + Cauer Ladder Network MOR for control-coupled maglev (Yano, CAE-AI)._

Module: `radia_mcp.maglev.server`

| Tool | Description |
|---|---|
| `maglev` | Magnetic levitation knowledge -- maglev systems + levitation force physics. |
| `maglev_status` | (no description) |
| `maglev_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 19 topics. |
| `rotating_conductor_periodic_settling_gate` | Gate full-turn convergence of a rotating-conductor eddy response. |

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
| `differential_forms_gauge_invariance_gate` | Gate physical B/loss invariance without treating A as invariant. |
| `differential_forms_homology` | Chain complex, homology, Betti numbers, tree-cotree gauge. |
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

## `mcp-server-chart2d`

_22 paper-quality 2D charts as MCP tools: line / loglog / semilog / step / errorbar / fill_between / bode / histogram / bar / box / violin / ecdf / contour / contourf / pcolormesh / quiver / streamplot / imshow / polar / scatter / phase (Nyquist).  Each accepts return_mode='recipe' (Python text) | 'image' (MCP Image inline) | 'both'.  Inherits radia_mcp.figure profile + gate stack._

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

## `mcp-server-paper-writing`

_Journal paper / digest writing helpers: IMRaD, abstract, citation, figure, equation, PDF layout, and reviewer-trigger lints. Also serves the merged presentation_* slide lint + PPTX tools (2026-07-17) and the merged figure_* / paper_figure_* publication-figure tools (2026-07-18: the standalone presentation and figure servers were retired -- AI cannot yet author slide decks end-to-end, and figure was a shared middle layer now unified here)._

Module: `radia_mcp.paper_writing.server`

| Tool | Description |
|---|---|
| `figure_audit_embeds` | Lint every \includegraphics in a LaTeX file for figure embeds that |
| `figure_design_principles` | The figure-MAKING (作図, *sakuzu*) DESIGN canon, distilled from the |
| `figure_diagram_recipes` | Flowchart + conceptual/schematic DIAGRAM recipes -- the diagram-DRAWING skill |
| `figure_everyday_recipe` | Matplotlib recipe for the lab's EVERYDAY analysis figure. |
| `figure_matlab2tikz_recipe` | Generate a MATLAB recipe that exports the current figure to TikZ |
| `figure_office_export_recipe` | MATLAB recipe to export the current figure for Word / PowerPoint |
| `figure_size_for_target` | Recommend output figure size + font settings for a target embedding. |
| `figure_style_guide` | Return the lab-standard graph style guide. |
| `paper_figure_profiles` | List paper-quality figure profiles + their exact journal geometry. |
| `paper_figure_quality_rules` | Why paper-quality figures need a margin-efficiency gate. |
| `paper_figure_recipe` | Generate a self-contained Python recipe for a paper-quality figure. |
| `paper_writing_abstract_strength` | Abstract の強度を 4 要素 (problem / method / result-with-number / impact) |
| `paper_writing_acronym_usage_audit` | 略語の使用頻度監査 (grant_writing 実装の re-export)。 |
| `paper_writing_adaptive_health_report` | paper T8 health_report の severity 判定を context で adjust。 |
| `paper_writing_analyze_sentences` | 文長分析 (和文)。journal では長文を避けて読みやすさ重視。 |
| `paper_writing_arxiv_extract_equations` | Extract all displayed equations from a LaTeX source. |
| `paper_writing_arxiv_fetch_latex_source` | Fetch the LaTeX source of an arXiv preprint. |
| `paper_writing_arxiv_search` | Search arXiv via the official Atom XML API. |
| `paper_writing_check_abstract_background_ratio` | Abstract 内で background (導入文) が占める割合を推定。 |
| `paper_writing_check_abstract_no_math_no_citation` | Abstract 内に数式 (math), citation, domain acronym が混入していないか検出。 |
| `paper_writing_check_citation_usage` | TeX 本文中の \cite{} キーと bib file の entries を突合。 |
| `paper_writing_check_digest_human_review_triggers` | Detect one-page digest issues learned from Sugahara human review. |
| `paper_writing_check_english_redflags` | 英文論文の典型的 red flag を検出 (冠詞、時制、自動詞/他動詞 の混同)。 |
| `paper_writing_check_equation_numbering` | 方程式番号 (1), (2), ... の連番欠落 / 重複をチェック。 |
| `paper_writing_check_figure_caption_showing` | Figure caption が showing (describe) 形か telling (claim) 形か判定。 |
| `paper_writing_check_figure_forward_reference` | 図/表の \label と \ref の整合チェック (孤立ラベル / 未解決参照)。 |
| `paper_writing_check_figure_uses_pdf` | `\includegraphics` paths must be vector (.pdf / .eps), not raster. |
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
| `paper_writing_normalize_terminology` | Normalize known terminology variants in paper text. |
| `paper_writing_normalize_terminology_file` | Normalize known terminology variants in a TeX/text file. |
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
| `presentation_acronym_usage_audit` | 略語使用頻度監査 (re-export)。 |
| `presentation_adaptive_health_report` | pptx health_report の severity を venue で adjust。 |
| `presentation_add_citation_footer` | Add a small citation footnote textbox along the BOTTOM of one |
| `presentation_analyze_sentences` | 文長分析。スライドは短文指向。 |
| `presentation_arrow_usage` | 矢印 shape (line connector with arrow) の過剰使用検出。 |
| `presentation_chart_simplification_check` | Chart 簡素化 (Cole Knaflic style) の 5 軸診断。 |
| `presentation_check_bullet_count_per_slide` | 1 slide の bullet 数が上限超過を検出 (Miller 7±2). |
| `presentation_check_bullet_ending_style` | bullet 末尾の「。」有無が統一されているか. |
| `presentation_check_color_accessibility` | R+G 近接色ペアを検出 (protanopia/deuteranopia で区別困難). |
| `presentation_check_color_count_per_slide` | 宮野『研究発表のためのスライドデザイン』S12: 3 色使い原則の検査 (v0.9.0)。 |
| `presentation_check_hedge_on_key_slides` | pptx で Result / Conclusion / Summary スライドに弱気修飾語が |
| `presentation_check_image_text_ratio` | 1 slide の image 面積比が min 未満を検出 (Zen style). |
| `presentation_check_kanji_ratio` | スライド台本の漢字比率 check (re-export)。 |
| `presentation_check_logo_on_every_slide` | 全スライドに同じロゴ画像が繰り返し配置されているかを検出。 |
| `presentation_check_misuse_japanese` | 台本の現代誤用検出 (re-export)。 |
| `presentation_check_notation_variants` | スライドテキストの表記ゆれ検出 (re-export)。 |
| `presentation_check_over_politeness` | 学会発表で過剰に丁寧な言い回しを検出。木下 p.235。 |
| `presentation_check_overfull_hbox` | beamer ログ中の Overfull \hbox をカウント。スライドでは致命的。 |
| `presentation_check_pie_3d_charts` | pptx 内の chart shape を走査し、pie / doughnut / 3D chart を NG 検出。 |
| `presentation_check_pptx_font_size` | pptx font size < 下限を検出。 |
| `presentation_check_progress_indicator` | outline / section-header slides for progress indication を検出。 |
| `presentation_check_qa_backup_slides` | pptx に Q&A backup slide (hidden or named) が N 枚以上あるか確認. |
| `presentation_check_script_paragraph_length` | 発表原稿の 1 パラグラフが 200-300 字目安から大きく外れていないか。 |
| `presentation_check_slide_density` | 1 スライドあたりの文字密度チェック (テキストを直接渡す)。 |
| `presentation_check_slide_line_count` | pptx の各 slide で text 行数が木下推奨の範囲内か検証。 |
| `presentation_check_slide_title_verb` | 各 slide の title が claim 形式か名詞句止まりか。 |
| `presentation_check_takehome_slide` | pptx 最終 3 枚以内に Take-home / Summary / まとめ slide があるか確認。 |
| `presentation_check_time_13_rule` | 木下 1/3 則 — 前半で全員わかる話、中盤で大半が分かった気、 |
| `presentation_check_time_14_rule` | 木下 1/4 則 — 10 分講演を 4 等分 (intro/method/result/discussion) した |
| `presentation_check_underline_in_pptx` | pptx runs の font.underline を走査して下線密度を診断する。 |
| `presentation_citation_audit` | Check numeric ``[N]`` citations on the slides against the |
| `presentation_cite_format` | Format ONE reference into the styles used on talk slides. |
| `presentation_count_slides` | スライド数を count。beamer (.tex) の \begin{frame} か、 |
| `presentation_count_underlines` | beamer ソース内の下線コマンドを実測。 |
| `presentation_count_weak_expressions` | 弱気修飾語の出現。presentation では key slide 上で使うと信頼感低下。 |
| `presentation_embed_tts_audio_in_pptx` | Generate per-slide TTS MP3 and embed it into a PowerPoint deck. |
| `presentation_equation_slide_compliance` | 数式 slide の理系プレゼン compliance を診断。 |
| `presentation_estimate_per_slide_time` | 原稿を slide 境界で分割し、各 slide の発表時間を推定. |
| `presentation_estimate_speaking_time` | 原稿テキストから発表時間を推定。 |
| `presentation_extract_pptx_text` | pptx の各 slide のテキストを抽出。密度チェックや文字起こしに。 |
| `presentation_figure_slide_compliance` | Figure / chart slide の理系プレゼン compliance を診断。 |
| `presentation_find_undefined_acronyms` | スライド内略語の初出定義 check (re-export)。 |
| `presentation_font_consistency` | deck 内で使用されているフォントファミリーの数を集計。 |
| `presentation_health_report` | presentation Plan B の全 T1-T11 を束ねた統合レポート。 |
| `presentation_lint_bedrock` | 台本・スライド注釈の bedrock lint (木下 10 原則、re-export)。 |
| `presentation_mini_imrad_structure_check` | 理系プレゼンの mini-IMRAD 構造 (7 phases) 充足度を診断。 |
| `presentation_next_5_actions` | (no description) |
| `presentation_opening_hook_strength` | First 2 slides' text hook 強度診断。 |
| `presentation_qa_anticipation_list` | slide 内容から予想 Q&A 質問 list を生成。 |
| `presentation_qa_from_history` | ★ Q&A REHEARSAL: surface the real (or anticipated) questions asked |
| `presentation_references_slide` | Build a "References" slide from a list of full citation lines. |
| `presentation_results_slide_statistical_evidence` | Results slide の統計報告 4 要素 compliance を診断 (paper T12 の plot 版)。 |
| `presentation_rewrite_suggest` | (no description) |
| `presentation_rikei_minimalism_score` | 理系プレゼンの minimalism 5 軸を per slide 診断。 |
| `presentation_root_cause_diagnosis` | pptx の health_report を横断 pattern matching、根本原因診断。 |
| `presentation_run_full_workflow` | pptx を 1 コール chain 実行。 |
| `presentation_script_vs_slide_coverage` | 台本 (speaker_note) が slide 内容を網羅しているかを per slide 診断。 |
| `presentation_single_message_per_slide_semantic` | 1 主張 vs 複数テーマ slide を意味単位で診断。 |
| `presentation_slide_density_balance` | deck 内の char-count 分布の不均衡を Gini-like 指標で評価。 |
| `presentation_slide_titles_outline_coherence` | 全 slide のタイトルだけ並べて outline 化し、論理整合を診断。 |
| `presentation_speaker_note_ratio` | slide.notes_slide の speaker note 長 vs slide text 長の比率を検査。 |
| `presentation_speaking_pace_estimate` | speaker_note の文字数から WPM (日本語は文字数/分) で発表時間を推定。 |
| `presentation_suggest_redundancy_fixes` | 冗長表現 25 パターンの置換候補 (re-export)。 |
| `presentation_takehome_strength` | 最終スライド or last-3-slides の Take-home 品質診断。 |
| `presentation_talk_feedback_lookup` | Query the learned conference-talk field-note catalog. |
| `presentation_talk_feedback_stats` | Counts of the conference-talk field-note catalog (by venue / status |
| `presentation_text_density_per_slide_western_style` | 欧米式 text-heavy slide を検出し、日本理系向けに修正提案。 |
| `presentation_title_body_alignment_check` | Title が body の主張を要約しているかを per slide 診断。 |
| `presentation_usage` | 学会発表スライド (IEEJ SA / IEEE conference / セミナー) の作文技術ガイド全体。 |
| `presentation_validate_pdf_pages` | スライド PDF のページ数を実測。発表時間 / slot との整合を検証。 |
| `presentation_visual_text_ratio_score` | per-slide visual/text ratio の distribution を score 化。 |

## `mcp-server-grant-writing`

_Grant proposal helpers: Japanese technical-prose lint, section coverage, budget alignment, recommendation-letter template, and KDDI Digital Innovation social-implementation checks._

Module: `radia_mcp.grant_writing.server`

| Tool | Description |
|---|---|
| `grant_writing_acronym_usage_audit` | 略語の使用頻度と初出形式を監査し、3 段階の推奨を返す。 |
| `grant_writing_analyze_sentences` | Analyze Japanese sentence length for grant proposals. |
| `grant_writing_budget_alignment_check` | Check that budget items are tied to verification and implementation. |
| `grant_writing_check_kanji_ratio` | 漢字比率の偏りを検出。本多『日本語の作文技術』第四章に基づく。 |
| `grant_writing_check_misuse_japanese` | 『問題な日本語』由来の現代誤用 15 パターン検出。 |
| `grant_writing_check_notation_variants` | 同一テキスト内で同じ概念が複数の表記で書かれていないかを検出。 |
| `grant_writing_check_subject_predicate_distance` | 主述の直結原則 (本多 p.22): 主語と述語の間の距離が遠い文を検出。 |
| `grant_writing_count_weak_expressions` | Count hedges and grant-specific non-commitment phrases. |
| `grant_writing_find_undefined_acronyms` | Latin 略語の初出で定義 (〜 or 〜の略) が近くに無いものを検出。 |
| `grant_writing_health_report` | Integrated grant-writing health report. |
| `grant_writing_kddi_digital_check` | KDDI Foundation Digital Innovation / social implementation check. |
| `grant_writing_kddi_power_electronics_focus_check` | Check the current KDDI power-electronics-board CAE-AI framing. |
| `grant_writing_lint_bedrock` | 木下 10 原則 + 本多 + 知的 による和文技術文章 bedrock 診断。 |
| `grant_writing_recommendation_letter_template` | Return a one-page recommendation-letter draft template. |
| `grant_writing_section_presence` | Check whether a proposal draft contains the expected review axes. |
| `grant_writing_status` | (no description) |
| `grant_writing_suggest_redundancy_fixes` | 和文の典型的冗長表現 25 パターンを検出し置換候補を示す。 |
| `grant_writing_usage` | Return the grant-writing guide. |

## `mcp-server-poster`

_Conference poster generation and lint: templates, viewing-distance font size, color contrast, zone balance, QR audit, and print readiness._

Module: `radia_mcp.poster.server`

| Tool | Description |
|---|---|
| `poster_adaptive_health_report` | Health report tuned for a target conference's review culture. |
| `poster_betterposter_billboard_lint` | Lint the central billboard text for plain-language compliance. |
| `poster_caption_self_contained` | Score each ``\caption`` / ``\captionof{figure}`` for self-sufficiency. |
| `poster_color_contrast_wcag` | Check WCAG 2.1 AA contrast for *actually-paired* text/bg combinations. |
| `poster_color_count_321` | Count unique colors and warn if the palette violates the 3-color rule. |
| `poster_colorblind_hint` | Simulate deuteranopia and flag pairs that collapse to similar colors. |
| `poster_compile` | Compile a poster .tex to PDF. |
| `poster_elevator_pitch_generate` | Render a 3-minute speakable script from the poster source. |
| `poster_figures_audit` | Check that every ``\includegraphics{path}`` resolves to a file. |
| `poster_font_embed_check` | Run ``pdffonts`` and report any non-embedded font. |
| `poster_fontsize_by_distance` | Verify that ``\fontsize{X}`` values are large enough for their role. |
| `poster_from_paper_tex` | Convert a paper .tex into a Kelvin-style poster skeleton. |
| `poster_from_pptx` | Convert a PowerPoint poster draft to a Kelvin-style A1 .tex. |
| `poster_health_report` | Run Tier 1-2 lints and produce a weighted health score. |
| `poster_jp_font_check` | Inspect Japanese font family declarations in a poster .tex. |
| `poster_line_length` | Flag sentences that are too long for poster reading. |
| `poster_lint` | Lint a poster .tex against poster-specific (not slide) criteria. |
| `poster_next_5_actions` | Return the top 5 actions ranked by (impact - 0.5*effort). |
| `poster_print_readiness_audit` | Audit a poster for print-readiness: paper size + figure DPI. |
| `poster_qa_anticipation_list` | Anticipate likely poster Q&A and tag each with reviewer-type motivation. |
| `poster_qr_audit` | Audit QR code(s) in a poster for prominence + labeling + URL reachability. |
| `poster_qr_inject` | Inject a labeled QR code into a poster .tex. |
| `poster_rewrite_suggest` | Return 3-4 candidate phrasings for a target poster element. |
| `poster_root_cause_diagnosis` | Diagnose which of the 5 typical poster failure patterns apply. |
| `poster_run_full_workflow` | Chain the Intelligence Layer phases into one call. |
| `poster_skill_doc` | Return the poster sub-skill's ``skill.md`` documentation as text. |
| `poster_status` | (no description) |
| `poster_template_betterposter` | Return (or write) the A0-landscape #betterposter template. |
| `poster_template_kelvin` | Return (or write) the A1-portrait Japanese poster template. |
| `poster_typography_lints` | Run cheap regex-based typography hygiene checks. |
| `poster_word_count` | Lint a poster's word budget against Purrington's ≤1000-word target. |
| `poster_zone_balance_check` | Check that minipage column widths approximate 0.25 / 0.50 / 0.25. |

## `mcp-server-literature-index`

_★ Meta-MCP: full-text search across 3,889 lab literature files in public-safe curated corpus_

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
| `literature_semantic_search` | Semantic search over indexed text via ChromaDB + sentence-transformers |
| `literature_stats` | Index statistics + cache info. |

## `mcp-server-document-meta`

_Cross-cutting document/repo helpers: deadline, version diff, templates, lint-all, result-saving notebook audits, examples->docs/validation_test promotion audits, and root-level panels migration impact checks._

Module: `radia_mcp.document_meta.server`

| Tool | Description |
|---|---|
| `document_meta_deadline_countdown` | 任意の締切までの日数と推奨アクションを返す。 |
| `document_meta_diff_versions` | 2 つのテキスト file の unified diff を返す (作文 version 比較)。 |
| `document_meta_examples_migration_policy` | Return the current Radia examples/ migration policy. |
| `document_meta_examples_notebook_audit` | Audit examples -> docs/ipynb or validation_test promotion state. |
| `document_meta_lint_all` | Run every applicable radia-mcp lint over one text / TeX file. |
| `document_meta_notebook_result_audit` | Audit docs notebooks for saved results and synchronized result JSON. |
| `document_meta_panel_layout_audit` | Audit impact of moving panel surfaces toward repo-root ``panels/``. |
| `document_meta_status` | (no description) |
| `document_meta_template_loader` | 学術 document の定型 skeleton を返す。 |
| `document_meta_write_docs_notebook_result_jsons` | Batch-write synchronized result JSON sidecars for executed docs notebooks. |
| `document_meta_write_notebook_result_json` | Write a durable JSON sidecar summarising a saved-result notebook. |

## `mcp-server-radia-meta`

_★ RECOMMENDED FIRST CALL. Cross-server catalog of all radia_mcp.* servers — answers "which tool covers concept X?" without trial-and-error._

Module: `radia_mcp.meta.server`

| Tool | Description |
|---|---|
| `bug_patterns_lookup` | Query the learned bug-pattern catalog. |
| `bug_patterns_stats` | Counts of catalogued bug patterns by severity + topic. |
| `radia_mcp_by_tag` | Servers tagged with `tag`. |
| `radia_mcp_get` | Look up one server by short name (e.g. 'bayesian-opt', 'ih', 'kelvin'). |
| `radia_mcp_golden_gate` | Machine-readable golden-quality gate for the radia-mcp server fleet. |
| `radia_mcp_health` | Probe importability of every radia_mcp.* subpackage. |
| `radia_mcp_overview` | Authoritative catalog of all radia_mcp.* servers. |
| `radia_mcp_related` | Servers that pair well with `name` (e.g. radia_mcp_related('bayesian-opt') |
| `radia_meta_status` | (no description) |

## `mcp-server-panel-review`

_Radia notebook panel review and construction contract (DesignSpec / Workbench / result artifacts / validation_test / no-PySide gate), including the cubit_panels migration route._

Module: `radia_mcp.panel_review.server`

| Tool | Description |
|---|---|
| `panel_review` | Get Radia notebook panel review / construction documentation. |
| `panel_review_status` | (no description) |
| `panel_review_topics` | Authoritative list of topics accepted by this server's dispatcher tool. Returns 13 topics. |
