import asyncio
import json
import re
from pathlib import Path

from radia_mcp.matlab import (
    matlab_optuna_benchmark_plan,
    matlab_optuna_compatibility_contract,
    matlab_optuna_health,
    matlab_optuna_oracle_audit,
    matlab_optuna_oracle_plan,
    matlab_optuna_release_gate,
    matlab_optuna_simulink_contract,
    matlab_radia_mex_contract,
)


def _expected_optuna_health_errors() -> list[str]:
    """Errors the Optuna health gate must report for the checked ledger.

    Health reports the evidence backlog honestly. While the required public
    API is not completely evidence-mapped, that must be its only error, so no
    other distribution regression can hide behind it; once the ledger closes,
    no error is expected.
    """
    root = Path(__file__).resolve().parents[2]
    coverage = json.loads(
        (root / "packages/radia-optuna/matlab" / "optuna50_api_coverage.json").read_text(encoding="utf-8")
    )
    if coverage["full_compatibility_complete"]:
        return []
    return ["required public API scope is not completely evidence-mapped"]


def test_radia_mex_contract_routes_only_to_supported_gpu_owner():
    contract = matlab_radia_mex_contract("all")
    assert contract["schema"] == "radia-mcp.matlab-radia-mex/v2"
    assert contract["cpu_mex_supported"] is False
    assert contract["runtime_checked"] is False
    assert contract["command_names"] == []
    assert contract["mex_entrypoint"] == "radia_cuda_csr_pcg"
    assert contract["optuna_mex_command_count"] == 22
    sections = contract["topic_data"]["sections"]
    assert sections["mex"]["owner"] == "radia-simulink"
    assert sections["ngsolve"]["cpu_execution"] == "Python/pybind11/NGSolve"
    assert sections["simulink"]["owner"] == "radia-simulink"
    assert "No CPU Radia MEX" in sections["limitations"]["items"][0]


def test_acoustic_validation_references_stay_out_of_native_gateways(monkeypatch):
    root = Path(__file__).resolve().parents[2]
    monkeypatch.setenv("RADIA_REPO_ROOT", str(root))
    contract = matlab_radia_mex_contract("mex")

    assert not any(
        command.startswith("acoustic.") for command in contract["command_names"]
    )
    assert contract["cpu_mex_supported"] is False
    assert contract["runtime_checked"] is False
    assert "_Acoustic" not in (
        root / "src" / "lib" / "radia_pybind.cpp"
    ).read_text(encoding="utf-8")
    for path in sorted((root / "packages/radia-simulink/matlab" / "+radia" / "+acoustic").glob("*.m")):
        source = path.read_text(encoding="utf-8")
        assert "callMex" not in source
        assert "radia_mex" not in source


def test_optuna_simulink_contract_is_table_backed():
    contract = matlab_optuna_simulink_contract()

    assert contract["status"] == "ready"
    assert contract["package"] == "radia.optuna"
    assert contract["distribution"] == "radia-optuna"
    assert contract["upstream_oracle_version"] == "optuna==5.0.0"
    assert contract["upstream_oracle"]["oracle_owner"].startswith(
        "optuna==5.0.0"
    )
    assert contract["mcp_ownership"]["routes"]["shared"]["owner"] == (
        "optuna/optuna-mcp"
    )
    assert "TrialTable" in contract["tables"]
    assert "ObjectiveTable" in contract["tables"]
    assert "ConstraintTable" in contract["tables"]
    assert "SamplerStateTable" in contract["tables"]
    assert contract["schema"].endswith("/v5")
    assert contract["upstream_oracle"]["ok"] is True
    assert contract["upstream_oracle"]["oracle_versions"]["optuna"] == "5.0.0"
    assert contract["native_acceleration"]["upstream_python_gp_python_per_trial"] is True
    assert contract["native_acceleration"]["full_optimizer_in_cpp"] is False
    assert contract["native_acceleration"]["gateway"] == "optuna_mex"
    assert contract["native_acceleration"]["command_count"] == (
        contract["distribution_health"]["native"]["command_count"]
    )
    expected_errors = _expected_optuna_health_errors()
    assert contract["distribution_health"]["errors"] == expected_errors
    assert contract["distribution_health"]["ok"] is (not expected_errors)
    assert contract["native_acceleration"]["required"] is True
    assert contract["native_acceleration"]["missing_mex_fallback"] is False
    assert contract["cae_trial_contract"]["success_schema"] == (
        "radia.optuna.cae-trial.v1"
    )
    assert contract["cae_trial_contract"]["failure_schema"] == (
        "radia.optuna.cae-failure.v1"
    )
    assert "SimulinkRunner" in contract["classes"]
    assert "OptimizationSession" in contract["classes"]
    assert "OptimizationParameter" in contract["classes"]
    assert "SheetMetalRunner" in contract["classes"]
    assert contract["multi_objective"]["selection"].startswith("bestTrial")
    assert contract["multi_objective"]["samplers"] == [
        "RandomSampler", "TPESampler", "NSGAIISampler"
    ]
    assert "parsim" in contract["parallel_trials"]["simulink"]
    assert "parfeval" in contract["parallel_trials"]["ltspice"]
    assert "complex" in contract["ltspice_integrated_workflow"]["raw"]
    assert len(contract["simulink_blocks"]) == 3
    assert contract["cpu_radia_mex_supported"] is False
    assert contract["runtime_checked"] is False
    historical_blocks = contract["legacy_domain_models"]["descriptions"]
    assert "radia_optuna_teaching.slx" in contract["simulink_blocks"][-1]
    assert "distributed-field kernels" in historical_blocks[0]
    assert "IHMonitorBusV1" in historical_blocks[0]
    optuna_block = next(
        item for item in contract["simulink_blocks"]
        if "OptunaMonitorBusV1" in item
    )
    assert "six-input/eighteen-output" in optuna_block
    assert "LUT and lumped IH builders are removed" in historical_blocks[2]
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
    foster = contract["hcurl_eddy_foster"]
    assert foster["moving_family"].startswith("ExportHCurlEddyFosterFamilyJSON")
    assert "no MEX kernel" in foster["harmonic_solve"]
    assert "instantaneous K*z(t)*i(t)" in foster["force"]
    native_family = foster["native_motor_angle_family"]
    assert native_family["matlab_factory"] == "radia.simulink.makeMotorAngleFamily"
    assert native_family["simulink_builder"] == (
        "radia.simulink.buildMotorAngleFamilyModel"
    )
    assert "verified_tests" not in native_family
    artifact = json.loads(
        (Path(__file__).resolve().parents[2] / native_family["validation_artifact"]).read_text(
            encoding="utf-8"
        )
    )
    assert native_family["recorded_tests"] == artifact["test_count"]
    assert artifact["passed_count"] == artifact["test_count"]
    assert native_family["evidence_status"].startswith("recorded historical")
    assert native_family["status"] == "unsupported-cpu-mex"
    assert "do not certify current CUDA support" in native_family["evidence_status"]
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
    assert adjoint["status"] == "source-only; CUDA acceptance pending"
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
    assert hcurl["status"] == "unsupported-cpu-mex"
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


def test_root_readme_publishes_separate_cpu_and_gpu_ownership():
    root = Path(__file__).resolve().parents[2]
    readme = (root / "README.md").read_text(encoding="utf-8")
    assert "CPU execution belongs to `radia`" in readme
    assert "GPU kernels belong to `radia-simulink`" in readme
    assert "explicitly selected CPU route compiled with mexcuda" in readme
    assert "creates no CUDA context" in readme
    assert "no implicit fallback" in readme
    assert "Build.ps1 -MatlabMexOnly" not in readme


def test_optuna_quality_helpers_are_exported_from_the_matlab_package():
    expected_errors = _expected_optuna_health_errors()
    expected_status = "error" if expected_errors else "ready"
    health = matlab_optuna_health()
    assert health["errors"] == expected_errors
    assert health["ok"] is (not expected_errors)
    assert health["distribution"]["matlab_file_count"] == (
        health["distribution"]["expected_matlab_file_count"]
    )
    assert matlab_optuna_oracle_plan()["status"] == expected_status
    assert matlab_optuna_benchmark_plan()["status"] == expected_status
    assert callable(matlab_optuna_release_gate)


def test_radia_matlab_tools_do_not_shadow_verified_upstream_optuna_mcp_tools():
    from radia_mcp.matlab.server import mcp

    root = Path(__file__).resolve().parents[2]
    fixture = json.loads(
        (root / "packages" / "radia-optuna" / "oracles" /
         "optuna50_mcp_oracle.json").read_text(encoding="utf-8")
    )
    upstream_tools = set(fixture["tools"])
    radia_tools = {item.name for item in asyncio.run(mcp.list_tools())}

    assert fixture["optuna_version"] == "5.0.0"
    assert fixture["optuna_mcp_version"] == "0.2.0"
    assert upstream_tools.isdisjoint(radia_tools)


def test_optuna_compatibility_and_oracle_audit_are_checked():
    contract = matlab_optuna_compatibility_contract()
    assert contract["ok"] is True
    assert contract["transport"]["public_mcp_contract"] == "stdio"
    assert contract["transport"]["mcp_sampler_seed_supported"] is False
    closure = contract["public_api_closure"]
    assert closure["surface_entry_count"] == 812
    assert closure["surface_present_count"] == 812
    assert closure["surface_missing_count"] == 0
    coverage = json.loads((Path(__file__).resolve().parents[2] /
        "packages/radia-optuna/matlab/optuna50_api_coverage.json").read_text(encoding="utf-8"))
    entries = coverage["entries"]
    required = [entry for entry in entries if entry["scope"] == "required"]
    assert closure["oracle_verified_count"] == sum(e["oracle_status"] == "verified" for e in entries)
    assert closure["oracle_asserted_count"] == sum(e["oracle_status"] == "asserted" for e in entries)
    assert closure["oracle_partial_count"] == 0
    assert closure["oracle_unmapped_count"] == 0
    assert closure["required_entry_count"] == 401
    assert closure["required_oracle_mapped_count"] == sum(e["oracle_status"] == "verified" for e in required)
    assert closure["required_oracle_asserted_count"] == sum(e["oracle_status"] == "asserted" for e in required)
    assert closure["required_oracle_unmapped_count"] == sum(e["oracle_status"] != "verified" for e in required)
    assert closure["full_compatibility_complete"] is all(e["oracle_status"] == "verified" for e in required)

    audit = matlab_optuna_oracle_audit()
    assert audit["ok"] is True
    assert audit["test_function_count"] == contract["test_counts"]["total"]
    assert audit["manifest_entry_count"] == contract["test_counts"]["total"]
    assert audit["policy_identical"] is True
    assert audit["missing_manifest_entries"] == []
    assert audit["stale_manifest_entries"] == []
