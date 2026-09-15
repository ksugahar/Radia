import asyncio
import json

from radia_mcp.matlab import (
    matlab_cad_topology_build,
    matlab_optimize_build,
    matlab_optimize_resume,
    matlab_optuna_mcp_route,
    matlab_sheet_metal_topology_build,
)


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
    assert "matlab_validation_catalog" in tool_names
    assert "matlab_validation_run" in tool_names
    catalog = mcp._tool_manager._tools["matlab_validation_catalog"].fn()
    assert "matlab_optuna_release_gate" in {
        item["name"] for item in catalog["operations"]
    }
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
    assert differential["behavioral_oracle"] == "optuna==5.0.0"
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


def test_optimize_server_builds_multiobjective_ltspice_code():
    payload = matlab_optimize_build({
        "name": "loss-ripple",
        "directions": ["minimize", "minimize"],
        "sampler": "tpe",
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
    assert "radia.optuna.TPESampler" in payload["matlab_code"]
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
            "constant_liar": True,
        },
        "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
    })
    assert "GammaFcn=@customGamma" in advanced_tpe["matlab_code"]
    assert "WeightsFcn=@customWeights" in advanced_tpe["matlab_code"]
    assert "Multivariate=true,Group=true" in advanced_tpe["matlab_code"]
    assert "WarnIndependentSampling=true" in advanced_tpe["matlab_code"]
    assert "ConstantLiar=true" in advanced_tpe["matlab_code"]

    mutated_nsga = matlab_optimize_build({
        "directions": ["minimize", "minimize"],
        "sampler": {
            "name": "nsgaii", "seed": 13,
            "mutation": {"name": "polynomial", "eta": 15},
        },
        "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
    })
    assert "Mutation=radia.optuna.nsgaii.PolynomialMutation(Eta=15)" in (
        mutated_nsga["matlab_code"]
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
    automatic_grouped = matlab_optimize_build({
        "sampler": {"name": "tpe", "group": True},
        "runner": {"kind": "objective", "objective_fcn": "objectiveFcn"},
    })
    assert "Group=true" in automatic_grouped["matlab_code"]
    assert "Multivariate=" not in automatic_grouped["matlab_code"]
    with pytest.raises(ValueError, match="requires sampler.multivariate=true"):
        matlab_optimize_build({
            "sampler": {
                "name": "tpe", "group": True, "multivariate": False,
            },
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
