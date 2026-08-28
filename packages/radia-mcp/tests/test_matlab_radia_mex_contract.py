import asyncio
import json
from pathlib import Path

from radia_mcp.matlab import (
    matlab_cad_topology_build,
    matlab_optimize_build,
    matlab_optimize_resume,
    matlab_optuna_benchmark_plan,
    matlab_optuna_compatibility_contract,
    matlab_optuna_health,
    matlab_optuna_mcp_route,
    matlab_optuna_oracle_audit,
    matlab_optuna_oracle_plan,
    matlab_optuna_release_gate,
    matlab_optuna_simulink_contract,
    matlab_radia_mex_contract,
    matlab_sheet_metal_topology_build,
)


def test_radia_mex_contract_reads_the_cpp_command_inventory():
    contract = matlab_radia_mex_contract("mex")

    assert contract["status"] == "ready"
    assert contract["command_count"] == 364
    assert not any(
        command.startswith("optuna.")
        for command in contract["command_names"]
    )
    assert contract["matlab_wrapper_count"] >= 133
    assert contract["matlab_optuna_distribution_health"]["ok"] is True
    assert contract["matlab_optuna_file_count"] == (
        contract["matlab_optuna_expected_file_count"]
    )
    root = Path(__file__).resolve().parents[3]
    coverage = json.loads(
        (root / "matlab" / "optuna49_api_coverage.json").read_text(
            encoding="utf-8"
        )
    )
    public_count = len(coverage["entries"])
    assert contract["matlab_optuna_public_api"] == {
        "entry_count": public_count,
        "present_count": public_count,
        "missing_count": 0,
        "verified_count": public_count,
        "partial_count": 0,
        "unmapped_count": 0,
        "complete": True,
    }
    assert contract["optuna_mex_command_count"] == 20
    assert contract["matlab_optuna_class_count"] >= 92
    assert contract["matlab_optuna_function_count"] >= 85
    assert {
        "TPESampler",
        "MOTPESampler",
        "CmaEsSampler",
        "GPSampler",
        "NSGAIISampler",
        "NSGAIIISampler",
        "BruteForceSampler",
        "LiveMonitor",
    }.issubset(
        contract["matlab_optuna_classes"]
    )
    assert "ngsolve.matrix_dump" in contract["command_names"]
    assert "ngsolve.mesh.set_deformation" in contract["command_names"]
    assert "ngsolve.mesh.trafo_quality" in contract["command_names"]
    assert "simulink.state_space.output" in contract["command_names"]
    assert "simulink.state_space.update" in contract["command_names"]
    assert "simulink.state_space.snapshot" in contract["command_names"]
    assert "simulink.state_space.restore" in contract["command_names"]
    optuna_commands = {
        "optuna.pareto.rank_crowding",
        "optuna.parzen.log_pdf_numerical",
        "optuna.parzen.log_pdf_categorical",
        "optuna.tpe.best_numerical",
        "optuna.tpe.best_joint",
        "optuna.tpe.best_numerical_observations",
        "optuna.tpe.best_joint_observations",
        "optuna.tpe.history.reset",
        "optuna.tpe.history.append_complete",
        "optuna.tpe.best_grouped_history",
        "optuna.random_state.create",
        "optuna.random_state.rand",
        "optuna.random_state.randn",
        "optuna.random_state.randi",
        "optuna.random_state.randperm",
        "optuna.random_state.snapshot",
        "optuna.random_state.restore",
        "optuna.random_state.destroy",
    }
    assert optuna_commands.issubset(contract["optuna_mex_command_names"])
    assert not any(name.startswith("optuna.") for name in contract["command_names"])
    assert "hacapk.charge_gram.reduce_configured_candidate_directional_schur" in contract["command_names"]
    assert {
        "ih.eddy.create",
        "ih.eddy.output",
        "ih.eddy.destroy",
        "ih.thermal.create",
        "ih.thermal.output",
        "ih.thermal.update",
        "ih.thermal.reset",
        "ih.thermal.destroy",
    }.issubset(contract["command_names"])
    assert "hdiv.field_evaluator.from_cloud" in contract["command_names"]
    assert "hacapk.charge_gram.configure_charge_map" in contract["command_names"]
    assert "hacapk.charge_gram.set_image_rotations" in contract["command_names"]
    assert "hacapk.charge_gram.configured_linear_material_element_blocks" in contract["command_names"]
    assert "hacapk.charge_gram.configured_linear_material_candidate_clusters" in contract["command_names"]
    assert "hlu.set_trunc_tol" in contract["command_names"]
    # The native Lie-map pipeline and 3D reference-orbit tracker remain on the
    # regular Radia gateway after the Optuna commands are split out.
    assert {
        "beam.lie.map_tensors_spoly",
        "beam.lie.dragt_finn_factorize",
        "beam.lie.apply_dragt_finn_batch",
        "beam.orbit.track_reference_3d",
        "beam.orbit.track_reference_to_plane",
    }.issubset(contract["command_names"])
    assert {
        "topopt.abe_element_fill_plan",
        "beam.orbit.track_reference_to_plane",
        "hacapk.charge_gram.configured_field_values_shape_derivative",
        "hacapk.charge_gram.configured_active_hmatrix_stats",
        "hdiv.field_evaluator.field_gradient",
    }.issubset(contract["command_names"])
    assert contract["command_groups"]["radia-core"] >= 70
    assert contract["pybind_public_count"] == 99
    assert contract["pybind_covered_count"] == 99
    assert contract["pybind_missing"] == []
    assert contract["pybind_internal_numerical_count"] == 28
    assert contract["pybind_internal_missing"] == []
    assert contract["pybind_internal_unclassified"] == []
    assert contract["command_groups"]["acoustic"] == 7
    assert "acoustic.elastic_sphere" in contract["command_names"]
    assert contract["command_groups"]["axifem"] == 2
    assert "axifem.q1_magnetic_element_matrices" in contract["command_names"]
    assert "axifem.q2_magnetic_element_matrices" in contract["command_names"]
    # _ChargeGramHMatrix.charge_sigma (the sigma-normalization diagnostic
    # from the roundoff-amplification fix) is EXCLUDED with a reason, and
    # exclusions leave the relevant surface; cyclic image setup expands the
    # covered stateful surface to 126 entries, including field-gradient and
    # configured directional-Schur and shape-derivative bindings.
    assert contract["pybind_class_surface_count"] == 126
    assert ("_ChargeGramHMatrix.charge_sigma"
            in contract["pybind_class_exclusions"])
    assert (
        "_ChargeGramHMatrix._reduce_configured_candidate_directional_schur"
        in contract["pybind_class_exclusions"]
    )
    assert contract["pybind_class_covered_count"] == contract["pybind_class_surface_count"]
    assert contract["pybind_class_missing_commands"] == []
    assert contract["pybind_class_unmapped"] == []
    assert contract["parity_status"] == "complete_for_radia_pybind_numerics"
    assert contract["retired_unsafe_constructor_leaks"] == []
    assert contract["retired_unsafe_c_abi_leaks"] == []
    for name in contract["retired_unsafe_constructors"]:
        assert name not in contract["command_names"]
    assert "ngsolve.matrix_dump" in contract["ngsolve_boundary"]["project_bridge_commands"]
    assert "radia.ngsolve.space_info" in contract["ngsolve_boundary"]["canonical_matlab_names"]
    assert "radia.ngsolve.matrix_dump" in contract["ngsolve_boundary"]["canonical_matlab_names"]
    assert contract["verified_contract"]["python_numerical_parity_gate"] == (
        "runtests('tests/matlab/test_radia_ngsolve_parity.m')"
    )
    assert contract["verified_contract"]["acoustic_python_mex_parity_gate"] == (
        "runtests('tests/matlab/test_acoustic_mex.m')"
    )
    assert contract["verified_contract"]["axifem_python_mex_parity_gate"] == (
        "runtests('tests/matlab/test_axifem_mex.m')"
    )
    assert contract["verified_contract"]["hcurl_topology_python_mex_parity_gate"] == (
        "runtests('tests/matlab/test_hcurl_topology_optimization.m')"
    )
    assert contract["verified_contract"]["topology_two_level_gate"] == (
        "runtests('tests/matlab/test_topology_optimization.m')"
    )
    assert "testOptunaParzenLogPdfKernels" in contract["verified_contract"][
        "optuna_native_kernel_gate"
    ]
    assert contract["verified_contract"]["optuna_native_kernel_benchmark"].endswith(
        "results_matlab_optuna_mex_benchmark_20260806.json"
    )
    assert contract["verified_contract"]["optuna49_performance_benchmark"].endswith(
        "results_matlab_optuna49_performance_20260825.json"
    )
    assert contract["verified_contract"]["native_motor_family_artifact"].endswith(
        "native_motor_angle_family.json"
    )
    assert "libiomp5md.dll" in contract["verified_contract"][
        "openmp_runtime_policy"
    ]
    all_contract = matlab_radia_mex_contract("all")
    family = all_contract["topic_data"]["sections"]["simulink"][
        "native_state_space_overloads"
    ]["periodic_motor_family"]
    assert "mechanical_angle" in family
    assert "torque" in family


def test_optuna_simulink_contract_is_table_backed():
    contract = matlab_optuna_simulink_contract()

    assert contract["status"] == "ready"
    assert contract["package"] == "radia.optuna"
    assert contract["distribution"] == "radia-optuna"
    assert contract["upstream_oracle_version"] == "optuna==4.9.0"
    assert contract["upstream_oracle"]["oracle_owner"].startswith(
        "optuna==4.9.0"
    )
    assert contract["mcp_ownership"]["routes"]["shared"]["owner"] == (
        "optuna/optuna-mcp"
    )
    assert "TrialTable" in contract["tables"]
    assert "ObjectiveTable" in contract["tables"]
    assert "ConstraintTable" in contract["tables"]
    assert "SamplerStateTable" in contract["tables"]
    assert contract["schema"].endswith("/v3")
    assert contract["upstream_oracle"]["ok"] is True
    assert contract["upstream_oracle"]["oracle_versions"]["optuna"] == "4.9.0"
    assert contract["native_acceleration"]["upstream_python_gp_python_per_trial"] is True
    assert contract["native_acceleration"]["full_optimizer_in_cpp"] is False
    assert contract["native_acceleration"]["gateway"] == "optuna_mex"
    assert contract["native_acceleration"]["command_count"] == (
        contract["distribution_health"]["native"]["command_count"]
    )
    assert contract["distribution_health"]["ok"] is True
    assert contract["native_acceleration"]["required"] is True
    assert contract["native_acceleration"]["missing_mex_fallback"] is False
    assert contract["cae_trial_contract"]["success_schema"] == (
        "radia.optuna.cae-trial.v1"
    )
    assert contract["cae_trial_contract"]["failure_schema"] == (
        "radia.optuna.cae-failure.v1"
    )
    assert "SimulinkRunner" in contract["classes"]
    assert "SheetMetalRunner" in contract["classes"]
    assert contract["multi_objective"]["selection"].startswith("bestTrial")
    assert contract["multi_objective"]["samplers"] == [
        "RandomSampler", "MOTPESampler", "NSGAIISampler"
    ]
    assert "parsim" in contract["parallel_trials"]["simulink"]
    assert "parfeval" in contract["parallel_trials"]["ltspice"]
    assert "complex" in contract["ltspice_integrated_workflow"]["raw"]
    assert len(contract["simulink_blocks"]) == 9
    assert "distributed-field kernels" in contract["simulink_blocks"][0]
    assert "LUT and lumped IH builders are removed" in contract["simulink_blocks"][2]
    assert contract["team28"]["frequency_hz"] == 50
    team28 = contract["team28"]
    assert team28["validated_dynamic_scope"] == (
        "cycle_averaged_mechanical_motion"
    )
    assert team28["electromagnetic_model_class"] == (
        "fixed_frequency_cycle_averaged_force_height_lut"
    )
    assert team28["height_coupling"] == "quasi_steady_interpolation"
    assert team28["electromagnetic_state_transient_included"] is False
    assert team28["motional_emf_included"] is False
    assert team28["damping_identified_from_measurement"] is False
    assert "full_electromagnetic_transient" in team28["unsupported_claims"]
    assert team28["artifact_gate"].endswith(
        "team28_cycle_averaged_motion_gate"
    )
    assert contract["hcurl_eddy_cln"]["mex_kernel"] == "hybrid_vim.solve"
    assert contract["hcurl_eddy_cln"]["moving_family"].startswith("ExportHCurlEddyCLNFamilyJSON")
    native_family = contract["hcurl_eddy_cln"]["native_motor_angle_family"]
    assert native_family["matlab_factory"] == "radia.simulink.makeMotorAngleFamily"
    assert native_family["simulink_builder"] == (
        "radia.simulink.buildMotorAngleFamilyModel"
    )
    assert native_family["verified_tests"] == 74
    assert contract["reinforcement_learning_workflow"]
    topology = contract["cad_topology_optimization"]
    assert topology["sensitivity_policy"].startswith("No cell-wise finite differences")
    assert "Cubit Sculpt" in topology["cad_reconstruction_route"]
    assert topology["cubit_validation_gates"] == [
        "cubit_ato_levelset_sculpt_source_replay_gate",
        "cubit_levelset_sculpt_hex_validation_gate",
    ]
    assert "existing hex mesh" in topology["boundary"]
    adjoint = topology["adjoint_optimization"]
    assert adjoint["status"] == "ready"
    assert "radia.topopt.optimizeAdjoint" in adjoint["matlab_api"]
    assert "radia.topopt.optimizeHCurlActivationAdjoint" in adjoint["matlab_api"]
    assert "never" in adjoint["finite_difference_policy"]
    assert topology["sheet_metal"]["two_level_loop"]["inner"].startswith("5-20")
    assert topology["sheet_metal"]["optuna_runner"] == "radia.optuna.SheetMetalRunner"
    assert topology["sheet_metal"]["simulink_block"] == (
        "Optimization/Sheet Metal Optimization"
    )
    assert contract["sampler_quality"]["python_parity_claim"].startswith(
        "Only behavior mapped"
    )
    hcurl = topology["sheet_metal"]["hcurl_eddy_bubble"]
    assert hcurl["status"] == "native-mex-ready"
    assert hcurl["python_boundary"] == "none in the MATLAB optimization loop"
    assert "radia.topopt.optimizeAdjoint" in hcurl["matlab_api"]
    assert "radia.topopt.optimizeHCurlActivationAdjoint" in hcurl["matlab_api"]
    assert {
        "hcurl.topopt.operator.*",
        "hcurl.topopt.resistance_shape_tangents",
        "hcurl.topopt.cell_curl_grams",
        "hcurl.topopt.multifrequency_joule",
        "hcurl.topopt.activation_multifrequency_joule",
    }.issubset(hcurl["mex_commands"])


def test_root_readme_publishes_native_topology_mex_parity():
    root = Path(__file__).resolve().parents[3]
    readme = " ".join((root / "README.md").read_text(encoding="utf-8").split())
    assert "HCurl-based reduced and topology workflows" in readme
    assert "HCurl multifrequency topology gradients" in readme

    matlab_readme = " ".join(
        (root / "matlab" / "README.md").read_text(encoding="utf-8").split()
    )
    assert "126 stateful class members" in matlab_readme
    assert "All 253 entries are covered by the current 364-command gateway" in matlab_readme
    assert "20-command `optuna_mex`" in matlab_readme

    parity_doc = (root / "docs" / "api" / "MATLAB_MEX_NGSOLVE_PARITY.md").read_text(
        encoding="utf-8"
    )
    assert "| Stateful pybind11 class surface | 126 / 126 covered |" in parity_doc
    assert "| Radia MEX gateway commands | 364 |" in parity_doc
    assert "| Optuna MEX gateway commands | 20 |" in parity_doc


def test_server_registers_bridge_tools():
    from radia_mcp.matlab.server import matlab_radia_mex_contract as mex_tool
    from radia_mcp.matlab.server import mcp

    tool_names = {item.name for item in asyncio.run(mcp.list_tools())}
    assert "matlab_radia_mex_contract" in tool_names
    assert "matlab_optuna_simulink_contract" in tool_names
    assert "matlab_optuna_mcp_route" in tool_names
    assert "matlab_optuna_health" in tool_names
    assert "matlab_optuna_oracle_plan" in tool_names
    assert "matlab_optuna_benchmark_plan" in tool_names
    assert "matlab_optuna_release_gate" in tool_names
    assert "matlab_optuna_compatibility_contract" in tool_names
    assert "matlab_optuna_oracle_audit" in tool_names
    assert "matlab_optimize_build" in tool_names
    assert "matlab_optimize_resume" in tool_names
    assert "matlab_cad_topology_build" in tool_names
    assert "matlab_sheet_metal_topology_build" in tool_names

    payload = json.loads(mex_tool("ngsolve"))
    assert payload["topic"] == "ngsolve"
    assert payload["topic_data"]["owner"] == "NGSolve"


def test_optuna_mcp_route_keeps_shared_tools_upstream_and_matlab_differences_local():
    contract = matlab_optuna_mcp_route()
    shared = contract["routes"]["shared"]
    matlab = contract["routes"]["matlab"]
    differential = contract["routes"]["differential"]
    stewardship = contract["routes"]["stewardship"]

    assert contract["policy"] == (
        "upstream for shared behavior; radia-mcp for MATLAB differences"
    )
    assert shared["owner"] == "optuna/optuna-mcp"
    assert "live MCP tools/list" in shared["authority"]
    assert shared["verified_snapshot"]["sampler_seed_exposed"] is False
    assert matlab["owner"] == "radia-mcp/radia-matlab"
    assert matlab["distribution"] == "radia-optuna"
    assert "matlab_optuna_health" in matlab["tools"]
    assert "matlab_optuna_oracle_plan" in matlab["tools"]
    assert "matlab_optuna_benchmark_plan" in matlab["tools"]
    assert "matlab_optuna_release_gate" in matlab["tools"]
    assert "table/MAT progress persistence and resume code generation" in (
        matlab["capabilities"]
    )
    assert "a second Optuna MCP server or optuna-mcp proxy" in (
        matlab["does_not_own"]
    )
    assert differential["behavioral_oracle"] == "optuna==4.9.0"
    assert "does not expose a seed" in differential["seeded_numeric_route"]
    assert stewardship["upstream_runtime_bundled"] is False
    assert stewardship["validation_operation"]["shared_or_production_storage"] is False
    assert stewardship["validation_operation"]["dashboard_in_automated_tests"] is False
    assert (
        stewardship["trademark_attribution"]
        == "Optuna, the Optuna logo and any related marks are trademarks of "
        "Preferred Networks, Inc."
    )
    assert {item["license"] for item in stewardship["upstream_licenses"]} == {"MIT"}


def test_optuna_quality_helpers_are_exported_from_the_matlab_package():
    health = matlab_optuna_health()
    assert health["ok"] is True
    assert health["distribution"]["matlab_file_count"] == (
        health["distribution"]["expected_matlab_file_count"]
    )
    assert matlab_optuna_oracle_plan()["status"] == "ready"
    assert matlab_optuna_benchmark_plan()["status"] == "ready"
    assert callable(matlab_optuna_release_gate)


def test_radia_matlab_tools_do_not_shadow_verified_upstream_optuna_mcp_tools():
    from radia_mcp.matlab.server import mcp

    root = Path(__file__).resolve().parents[3]
    fixture = json.loads(
        (root / "tests" / "matlab" / "fixtures" /
         "optuna49_mcp_oracle.json").read_text(encoding="utf-8")
    )
    upstream_tools = set(fixture["tools"])
    radia_tools = {item.name for item in asyncio.run(mcp.list_tools())}

    assert fixture["optuna_version"] == "4.9.0"
    assert fixture["optuna_mcp_version"] == "0.2.0"
    assert upstream_tools.isdisjoint(radia_tools)


def test_optuna_compatibility_and_oracle_audit_are_checked():
    contract = matlab_optuna_compatibility_contract()
    assert contract["ok"] is True
    assert contract["transport"]["public_mcp_contract"] == "stdio"
    assert contract["transport"]["mcp_sampler_seed_supported"] is False
    closure = contract["public_api_closure"]
    assert closure["surface_entry_count"] == 816
    assert closure["surface_present_count"] == 816
    assert closure["surface_missing_count"] == 0
    assert closure["oracle_verified_count"] == 816
    assert closure["oracle_partial_count"] == 0
    assert closure["oracle_unmapped_count"] == 0
    assert closure["full_compatibility_complete"] is True

    audit = matlab_optuna_oracle_audit()
    assert audit["ok"] is True
    assert audit["test_function_count"] == contract["test_counts"]["total"]
    assert audit["manifest_entry_count"] == contract["test_counts"]["total"]
    assert audit["policy_identical"] is True
    assert audit["missing_manifest_entries"] == []
    assert audit["stale_manifest_entries"] == []


def test_optimize_server_builds_multiobjective_ltspice_code():
    payload = matlab_optimize_build({
        "name": "loss-ripple",
        "directions": ["minimize", "minimize"],
        "sampler": "motpe",
        "n_trials": 24,
        "parallel": True,
        "runner": {
            "kind": "ltspice",
            "netlist": r"C:\temp\drive.cir",
            "configure_fcn": "configureDriveTrial",
            "score_fcn": "scoreDriveTrial",
        },
    })
    assert payload["runtime_owner"] == "MathWorks MATLAB MCP Server"
    assert "radia.optuna.MOTPESampler" in payload["matlab_code"]
    assert "runner.optimizeParallel(study,24)" in payload["matlab_code"]
    assert "pareto=study.paretoFront()" in payload["matlab_code"]

    resume = matlab_optimize_resume(r"C:\temp\loss-ripple.mat", 10, parallel=True)
    assert "runner.optimizeParallel(study,10)" in resume["matlab_code"]


def test_optimize_server_builds_cae_aware_native_simulink_code():
    payload = matlab_optimize_build({
        "name": "thermal-design",
        "directions": ["minimize"],
        "n_trials": 16,
        "parallel": True,
        "runner": {
            "kind": "simulink",
            "model": "radia_ih_design",
            "configure_fcn": "configureIHTrial",
            "score_fcn": "scoreIHTrial",
            "constraint_fcn": "constrainIHTrial",
            "validation_fcn": "validateIHTrial",
            "result_fcn": "collectIHArtifacts",
            "failure_classifier_fcn": "classifyIHFailure",
            "use_fast_restart": True,
            "continue_on_error": True,
            "batch_size": 3,
            "context": {"geometry": "workpiece", "mesh": "team36.vol"},
        },
    })
    code = payload["matlab_code"]
    assert payload["schema"].endswith("/v3")
    assert "radia.optuna.SimulinkRunner" in code
    assert "ConstraintFcn=@constrainIHTrial" in code
    assert "ValidationFcn=@validateIHTrial" in code
    assert "ResultFcn=@collectIHArtifacts" in code
    assert "FailureClassifierFcn=@classifyIHFailure" in code
    assert "UseFastRestart=true" in code
    assert "BatchSize=3,ContinueOnError=true" in code
    assert "Context=jsondecode" in code
    assert payload["result_contract"]["cae_success"] == (
        "radia.optuna.cae-trial.v1 in trial user attributes"
    )


def test_optimize_builder_v3_covers_seeded_sampler_surface():
    cases = {
        "random": "RandomSampler(41)",
        "tpe": "TPESampler(Seed=41",
        "cmaes": "CmaEsSampler(Seed=41",
        "gp": "GPSampler(Seed=41",
        "nsgaii": "NSGAIISampler(Seed=41",
        "nsgaiii": "NSGAIIISampler(Seed=41",
        "qmc": "QMCSampler(QMCType=\"sobol\",Scramble=true,Seed=41)",
        "bruteforce": "BruteForceSampler(Seed=41",
    }
    for name, expected in cases.items():
        sampler = {"name": name, "seed": 41}
        if name == "qmc":
            sampler["scramble"] = True
        payload = matlab_optimize_build({
            "directions": ["minimize", "minimize"] if name in {"nsgaii", "nsgaiii"} else ["minimize"],
            "sampler": sampler,
            "n_trials": 3,
            "live_monitor": False,
            "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
        })
        assert expected in payload["matlab_code"]
        assert payload["oracle"]["explicit_seed"] == 41
        assert payload["oracle"]["classification"] == "upstream-python"

    grid = matlab_optimize_build({
        "sampler": {"name": "grid", "seed": 7, "search_space": {"x": [1, 2]}},
        "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
    })
    assert "GridSampler(jsondecode" in grid["matlab_code"]

    fixed = matlab_optimize_build({
        "sampler": {
            "name": "partial_fixed", "seed": 9,
            "fixed_params": {"mode": "A"},
            "base_sampler": {"name": "random", "seed": 9},
        },
        "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
    })
    assert "PartialFixedSampler(jsondecode" in fixed["matlab_code"]
    assert "RandomSampler(9)" in fixed["matlab_code"]

    advanced_tpe = matlab_optimize_build({
        "sampler": {
            "name": "tpe", "seed": 13,
            "gamma_fcn": "customGamma", "weights_fcn": "customWeights",
            "multivariate": True, "group": True,
            "warn_independent_sampling": True,
            "categorical_distance_fcn": {"mode": "modeDistance"},
        },
        "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
    })
    assert "GammaFcn=@customGamma" in advanced_tpe["matlab_code"]
    assert "WeightsFcn=@customWeights" in advanced_tpe["matlab_code"]
    assert "Multivariate=true,Group=true" in advanced_tpe["matlab_code"]
    assert "WarnIndependentSampling=true" in advanced_tpe["matlab_code"]
    assert (
        "CategoricalDistanceFcn=containers.Map({'mode'},{@modeDistance})"
        in advanced_tpe["matlab_code"]
    )

    cma_independent = matlab_optimize_build({
        "sampler": {
            "name": "cmaes", "seed": 31,
            "independent_sampler": {"name": "random", "seed": 211},
            "warn_independent_sampling": False,
        },
        "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
    })
    assert "IndependentSampler=radia.optuna.RandomSampler(211)" in (
        cma_independent["matlab_code"]
    )
    assert "WarnIndependentSampling=false" in cma_independent["matlab_code"]


def test_optimize_builder_classifies_parallel_and_rejects_invalid_sampler_contracts():
    parallel = matlab_optimize_build({
        "sampler": {"name": "tpe", "seed": 17, "multivariate": True},
        "parallel": True,
        "runner": {
            "kind": "ltspice", "netlist": r"C:\temp\a.cir",
            "configure_fcn": "configureTrial", "score_fcn": "scoreTrial",
        },
    })
    assert parallel["oracle"]["classification"] == "matlab-integration"
    assert "Multivariate=true" in parallel["matlab_code"]

    import pytest

    with pytest.raises(ValueError, match="CmaEsSampler supports only one objective"):
        matlab_optimize_build({
            "directions": ["minimize", "maximize"], "sampler": "cmaes",
            "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
        })
    with pytest.raises(ValueError, match="requires non-empty sampler.search_space"):
        matlab_optimize_build({
            "sampler": "grid",
            "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
        })
    with pytest.raises(ValueError, match="requires sampler.multivariate=true"):
        matlab_optimize_build({
            "sampler": {"name": "tpe", "group": True},
            "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
        })


def test_optimize_server_builds_cubit_vim_lp_code():
    payload=matlab_cad_topology_build({
        "design_data":r"C:\temp\design.mat",
        "linearize_fcn":"linearizeMagnetVIM",
        "objective_weights":[1.0,-0.25],
        "volume_fraction":0.4,
        "move_limit":0.1,
        "output_journal":r"C:\temp\density.jou",
    })
    assert payload["status"]=="ready"
    assert "radia.topopt.optimizeVIMLP" in payload["matlab_code"]
    assert "radia.topopt.writeCubitJournal" in payload["matlab_code"]
    assert payload["gradient_policy"].startswith("no cell-wise finite differences")


def test_optimize_server_builds_sheet_metal_mesh_routing_code():
    payload=matlab_sheet_metal_topology_build({
        "design_data":r"C:\temp\sheet.mat","mesh_path":r"C:\temp\sheet.vol",
        "linearize_fcn":"linearizeSheetVIM","deformation_fcn":"makeSheetDeformation",
        "objective_fcn":"evaluateSheetObjective","rebuild_hmatrix_fcn":"rebuildSheetHMatrix",
        "cubit_rebuild_fcn":"rebuildSheetWithCubit","inner_iterations":10,
        "activation_remove_threshold":0.3,"activation_restore_threshold":0.7,
        "cubit_batch_interval":4,"cubit_batch_fraction":0.2,
    })
    assert payload["status"]=="ready"
    assert "radia.topopt.optimizeHexSheetTopology" in payload["matlab_code"]
    assert "InnerIterations=10" in payload["matlab_code"]
    assert payload["inner_iteration_range"]==[5,20]
    assert payload["activation_hysteresis"]=={
        "remove_threshold":0.3,"restore_threshold":0.7}
    assert payload["cubit_batching"]=={
        "maximum_pending_iterations":4,"pending_fraction":0.2}
    assert "CubitBatchInterval=4" in payload["matlab_code"]
    assert payload["hmatrix_rebuild_policy"].startswith("exactly once")
    assert payload["mesh_routes"]==["ngsolve_deform","ngsolve_refine","cubit_rebuild"]
