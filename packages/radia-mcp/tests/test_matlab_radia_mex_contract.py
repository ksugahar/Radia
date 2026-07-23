import asyncio
import json
from pathlib import Path

from radia_mcp.matlab import (
    matlab_cad_topology_build,
    matlab_optimize_build,
    matlab_optimize_resume,
    matlab_sheet_metal_topology_build,
    matlab_optuna_simulink_contract,
    matlab_radia_mex_contract,
)


def test_radia_mex_contract_reads_the_cpp_command_inventory():
    contract = matlab_radia_mex_contract("mex")

    assert contract["status"] == "ready"
    assert contract["command_count"] == 311
    assert contract["matlab_wrapper_count"] >= 133
    assert contract["matlab_optuna_class_count"] == 12
    assert {"TPESampler", "MOTPESampler", "CmaEsSampler", "NSGAIISampler", "LiveMonitor"}.issubset(
        contract["matlab_optuna_classes"]
    )
    assert "ngsolve.matrix_dump" in contract["command_names"]
    assert "ngsolve.mesh.set_deformation" in contract["command_names"]
    assert "ngsolve.mesh.trafo_quality" in contract["command_names"]
    assert "hdiv.field_evaluator.from_cloud" in contract["command_names"]
    assert "hacapk.charge_gram.configure_charge_map" in contract["command_names"]
    assert "hlu.set_trunc_tol" in contract["command_names"]
    assert contract["command_groups"]["radia-core"] >= 70
    assert contract["pybind_public_count"] == 94
    assert contract["pybind_covered_count"] == 94
    assert contract["pybind_missing"] == []
    assert contract["pybind_internal_numerical_count"] == 27
    assert contract["pybind_internal_missing"] == []
    assert contract["pybind_internal_unclassified"] == []
    assert contract["command_groups"]["acoustic"] == 7
    assert "acoustic.elastic_sphere" in contract["command_names"]
    assert contract["command_groups"]["axifem"] == 1
    assert "axifem.q1_magnetic_element_matrices" in contract["command_names"]
    assert contract["pybind_class_surface_count"] == 111
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


def test_optuna_simulink_contract_is_table_backed():
    contract = matlab_optuna_simulink_contract()

    assert contract["status"] == "ready"
    assert contract["package"] == "radia.optuna"
    assert "TrialTable" in contract["tables"]
    assert "ObjectiveTable" in contract["tables"]
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
    assert contract["hcurl_eddy_cln"]["mex_kernel"] == "hybrid_vim.solve"
    assert contract["hcurl_eddy_cln"]["moving_family"].startswith("ExportHCurlEddyCLNFamilyJSON")
    assert contract["reinforcement_learning_workflow"]
    topology = contract["cad_topology_optimization"]
    assert topology["sensitivity_policy"].startswith("No cell-wise finite differences")
    assert "Cubit Sculpt" in topology["cad_reconstruction_route"]
    assert topology["cubit_validation_gates"] == [
        "cubit_ato_levelset_sculpt_source_replay_gate",
        "cubit_levelset_sculpt_hex_validation_gate",
    ]
    assert "existing hex mesh" in topology["boundary"]
    assert topology["sheet_metal"]["two_level_loop"]["inner"].startswith("5-20")
    assert topology["sheet_metal"]["optuna_runner"] == "radia.optuna.SheetMetalRunner"
    assert topology["sheet_metal"]["simulink_block"] == (
        "Optimization/Sheet Metal Optimization"
    )
    assert contract["sampler_quality"]["python_parity_claim"].startswith("none")
    hcurl = topology["sheet_metal"]["hcurl_eddy_bubble"]
    assert hcurl["status"] == "native-mex-ready"
    assert hcurl["python_boundary"] == "none in the MATLAB optimization loop"
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
    assert "111 stateful class members" in matlab_readme
    assert "All 232 entries are covered by the current 311-command gateway" in matlab_readme

    parity_doc = (root / "docs" / "api" / "MATLAB_MEX_NGSOLVE_PARITY.md").read_text(
        encoding="utf-8"
    )
    assert "| Stateful pybind11 class surface | 111 / 111 covered |" in parity_doc
    assert "| MEX gateway commands | 311 |" in parity_doc


def test_server_registers_bridge_tools():
    from radia_mcp.matlab.server import mcp, matlab_radia_mex_contract as mex_tool

    tool_names = {item.name for item in asyncio.run(mcp.list_tools())}
    assert "matlab_radia_mex_contract" in tool_names
    assert "matlab_optuna_simulink_contract" in tool_names
    assert "matlab_optimize_build" in tool_names
    assert "matlab_optimize_resume" in tool_names
    assert "matlab_cad_topology_build" in tool_names
    assert "matlab_sheet_metal_topology_build" in tool_names

    payload = json.loads(mex_tool("ngsolve"))
    assert payload["topic"] == "ngsolve"
    assert payload["topic_data"]["owner"] == "NGSolve"


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
