"""Executable-manual helpers for MATLAB/Simulink/LTspice optimization."""

from __future__ import annotations

import json
import re
from typing import Any, Mapping

from .optuna_boundary import matlab_optuna_mcp_route


_FUNCTION = re.compile(r"^[A-Za-z]\w*(?:\.[A-Za-z]\w*)*$")
_SAMPLERS = {
    "random": "radia.optuna.RandomSampler(0)",
    "tpe": "radia.optuna.TPESampler(Seed=0)",
    "cmaes": "radia.optuna.CmaEsSampler(Seed=0)",
    "motpe": "radia.optuna.MOTPESampler(Seed=0)",
    "nsgaii": "radia.optuna.NSGAIISampler(Seed=0)",
}


def matlab_optimize_build(spec: Mapping[str, Any] | str) -> dict[str, Any]:
    """Validate an optimization spec and generate official-MCP-ready MATLAB."""
    if isinstance(spec, str):
        spec = json.loads(spec)
    if not isinstance(spec, Mapping):
        raise ValueError("optimization spec must be a JSON object")
    name = str(spec.get("name", "radia-study"))
    directions = [str(item) for item in spec.get("directions", ["minimize"])]
    if not directions or any(item not in {"minimize", "maximize"} for item in directions):
        raise ValueError("directions must contain minimize or maximize")
    sampler = str(spec.get("sampler", "motpe" if len(directions) > 1 else "tpe")).lower()
    if sampler not in _SAMPLERS:
        raise ValueError(f"unknown sampler {sampler!r}; expected {sorted(_SAMPLERS)}")
    if len(directions) > 1 and sampler in {"tpe", "cmaes"}:
        raise ValueError("multi-objective studies require random, motpe, or nsgaii")
    n_trials = int(spec.get("n_trials", 20))
    if n_trials <= 0:
        raise ValueError("n_trials must be positive")
    runner = dict(spec.get("runner", {"kind": "objective"}))
    kind = str(runner.get("kind", "objective")).lower()
    if kind not in {"objective", "simulink", "ltspice"}:
        raise ValueError("runner.kind must be objective, simulink, or ltspice")
    parallel = bool(spec.get("parallel", kind in {"simulink", "ltspice"}))
    storage = str(spec.get("storage", ""))
    live = bool(spec.get("live_monitor", True))
    code = _matlab_code(name, directions, sampler, n_trials, runner, parallel, storage, live)
    return {
        "schema": "radia-mcp.matlab-optimize-build/v2",
        "status": "ready",
        "runtime_owner": "MathWorks MATLAB MCP Server",
        "domain_owner": "radia-mcp.matlab.optimize",
        "mcp_ownership": matlab_optuna_mcp_route("composition")["route"],
        "execute_with": "official MATLAB MCP evaluate_matlab_code",
        "spec": {
            "name": name,
            "directions": directions,
            "sampler": sampler,
            "n_trials": n_trials,
            "runner": runner,
            "parallel": parallel,
            "storage": storage,
            "live_monitor": live,
        },
        "matlab_code": code,
        "result_contract": {
            "trials": "study.TrialTable",
            "objectives": "study.ObjectiveTable",
            "constraints": "study.ConstraintTable with c <= 0 feasibility",
            "cae_success": "radia.optuna.cae-trial.v1 in trial user attributes",
            "cae_failure": "radia.optuna.cae-failure.v1 in trial user attributes",
            "pareto": "study.paretoFront() for multi-objective studies",
        },
    }


def matlab_optimize_resume(storage_path: str, n_trials: int, *, parallel: bool = False) -> dict[str, Any]:
    """Generate code that resumes a MAT-file-backed Study."""
    if not storage_path.strip():
        raise ValueError("storage_path is required")
    if int(n_trials) <= 0:
        raise ValueError("n_trials must be positive")
    path = _quote(storage_path)
    return {
        "schema": "radia-mcp.matlab-optimize-resume/v1",
        "status": "ready",
        "execute_with": "official MATLAB MCP evaluate_matlab_code",
        "mcp_ownership": matlab_optuna_mcp_route("composition")["route"],
        "matlab_code": (
            f"study=radia.optuna.Study(StoragePath='{path}');\n"
            "% Recreate the same objective or runner, then execute one of:\n"
            f"% study.optimize(@objectiveFcn,{int(n_trials)});\n"
            f"% runner.optimize{'Parallel' if parallel else ''}(study,{int(n_trials)});\n"
            "pareto=study.paretoFront();\n"
        ),
    }


def matlab_cad_topology_build(spec: Mapping[str, Any] | str) -> dict[str, Any]:
    """Generate a Cubit + linearized-Radia-VIM + LP optimization workflow."""
    if isinstance(spec, str):
        spec=json.loads(spec)
    if not isinstance(spec,Mapping):
        raise ValueError("CAD topology spec must be a JSON object")
    design_data=str(spec.get("design_data", ""))
    if not design_data:
        raise ValueError("design_data MAT-file is required")
    linearize=_function(spec.get("linearize_fcn"),"linearize_fcn")
    weights=spec.get("objective_weights")
    if not isinstance(weights,list) or not weights:
        raise ValueError("objective_weights must be a non-empty numeric array")
    weights=[float(value) for value in weights]
    volume_fraction=float(spec.get("volume_fraction",0.5))
    if not 0<volume_fraction<=1:
        raise ValueError("volume_fraction must be in (0,1]")
    move_limit=float(spec.get("move_limit",0.2)); iterations=int(spec.get("max_iterations",30))
    journal=str(spec.get("output_journal","C:/temp/radia_topopt_density.jou"))
    threshold=float(spec.get("threshold",0.5))
    weight_expr="["+";".join(f"{value:.17g}" for value in weights)+"]"
    code=(
        f"data=load('{_quote(design_data)}');\n"
        f"result=radia.topopt.optimizeVIMLP(data.initialDensity,data.cellVolumes,{volume_fraction:.17g},@{linearize},"
        f"ObjectiveWeights={weight_expr},MoveLimit={move_limit:.17g},MaxIterations={iterations});\n"
        f"cubit=radia.topopt.writeCubitJournal('{_quote(journal)}',data.elementIds,result.density,Threshold={threshold:.17g});\n"
    )
    return {
        "schema":"radia-mcp.matlab-cad-topology-build/v1",
        "status":"ready",
        "runtime_owner":"MathWorks MATLAB MCP Server",
        "domain_owner":"radia-mcp.matlab.optimize",
        "method":"Cubit design cells -> analytic Radia-VIM linearization -> bounded LP -> Cubit material blocks",
        "gradient_policy":"no cell-wise finite differences; solve A*dm_i=db_i-dA_i*m for all design cells",
        "matlab_code":code,
        "design_data_contract":{
            "initialDensity":"n-by-1 values in [0,1]",
            "cellVolumes":"n-by-1 positive CAD/mesh cell volumes",
            "elementIds":"n-by-1 unique positive Cubit hex IDs",
        },
        "execute_with":"official MATLAB MCP evaluate_matlab_code",
    }


def matlab_sheet_metal_topology_build(spec: Mapping[str, Any] | str) -> dict[str, Any]:
    """Generate the two-level GetTrafo/Cubit sheet-topology workflow."""
    if isinstance(spec,str): spec=json.loads(spec)
    if not isinstance(spec,Mapping): raise ValueError("sheet-metal topology spec must be a JSON object")
    design_data=str(spec.get("design_data",""))
    if not design_data: raise ValueError("design_data MAT-file is required")
    mesh_path=str(spec.get("mesh_path",""))
    if not mesh_path: raise ValueError("mesh_path VOL-file is required")
    linearize=_function(spec.get("linearize_fcn"),"linearize_fcn")
    deformation=_function(spec.get("deformation_fcn"),"deformation_fcn")
    objective=_function(spec.get("objective_fcn"),"objective_fcn")
    rebuild_hmatrix=_function(spec.get("rebuild_hmatrix_fcn"),"rebuild_hmatrix_fcn")
    cubit_rebuild=_function(spec.get("cubit_rebuild_fcn"),"cubit_rebuild_fcn")
    inner_iterations=int(spec.get("inner_iterations",10))
    max_outer_iterations=int(spec.get("max_outer_iterations",10))
    activation_remove_threshold=float(spec.get("activation_remove_threshold",0.35))
    activation_restore_threshold=float(spec.get("activation_restore_threshold",0.65))
    cubit_batch_interval=int(spec.get("cubit_batch_interval",5))
    cubit_batch_fraction=float(spec.get("cubit_batch_fraction",0.05))
    if not 5<=inner_iterations<=20: raise ValueError("inner_iterations must be between 5 and 20")
    if max_outer_iterations<1: raise ValueError("max_outer_iterations must be positive")
    if not 0<=activation_remove_threshold<activation_restore_threshold<=1:
        raise ValueError("activation hysteresis must satisfy 0 <= remove < restore <= 1")
    if cubit_batch_interval<1 or not 0<cubit_batch_fraction<=1:
        raise ValueError("invalid Cubit batching controls")
    work_directory=str(spec.get("work_directory",r"C:\temp\radia_hex_topopt"))
    code=(
        f"data=load('{_quote(design_data)}');\n"
        f"mesh=radia.ngsolve.Mesh.create('{_quote(mesh_path)}');\n"
        "initialState=struct('mesh',mesh,'model',data.model,"
        "'normal_displacement',data.normalDisplacement,'thickness',data.thickness,"
        "'activation',data.activation,'objective',data.objective);\n"
        "result=radia.topopt.optimizeHexSheetTopology(initialState,"
        f"@{linearize},@{deformation},@{objective},@{rebuild_hmatrix},@{cubit_rebuild},"
        f"data.elementSizes,InnerIterations={inner_iterations},"
        f"MaxOuterIterations={max_outer_iterations},"
        f"ActivationRemoveThreshold={activation_remove_threshold:.17g},"
        f"ActivationRestoreThreshold={activation_restore_threshold:.17g},"
        f"CubitBatchInterval={cubit_batch_interval},"
        f"CubitBatchFraction={cubit_batch_fraction:.17g},"
        f"WorkDirectory='{_quote(work_directory)}');\n"
    )
    return {
        "schema":"radia-mcp.matlab-sheet-metal-topology-build/v3","status":"ready",
        "runtime_owner":"MathWorks MATLAB MCP Server","domain_owner":"radia-mcp.matlab.optimize",
        "method":"5-20 continuous-activation GetTrafo iterations -> conditional Cubit topology commit -> one H-matrix rebuild",
        "matlab_code":code,"execute_with":"official MATLAB MCP evaluate_matlab_code",
        "mesh_routes":["ngsolve_deform","ngsolve_refine","cubit_rebuild"],
        "inner_iteration_range":[5,20],
        "activation_hysteresis":{
            "remove_threshold":activation_remove_threshold,
            "restore_threshold":activation_restore_threshold,
        },
        "cubit_batching":{
            "maximum_pending_iterations":cubit_batch_interval,
            "pending_fraction":cubit_batch_fraction,
        },
        "hmatrix_rebuild_policy":"exactly once after a successful Cubit rebuild; never inside the GetTrafo inner loop",
        "design_data_contract":{"model":"current VIM/H-matrix model","normalDisplacement":"n-by-1","thickness":"n-by-1","activation":"n-by-1 continuous values in [0,1]","objective":"finite scalar","elementSizes":"n-by-1 positive values"},
    }


def _matlab_code(name, directions, sampler, n_trials, runner, parallel, storage, live):
    direction_expr = "[" + ",".join(f'\"{item}\"' for item in directions) + "]"
    progress = "monitor=radia.optuna.LiveMonitor();\nprogress=@monitor.update;" if live else "progress=[];"
    storage_arg = f",StoragePath='{_quote(storage)}'" if storage else ""
    lines = [
        progress,
        f"study=radia.optuna.createStudy(study_name=\"{_quote(name)}\",directions={direction_expr},"
        f"sampler={_SAMPLERS[sampler]},ProgressFcn=progress{storage_arg});",
    ]
    kind = str(runner.get("kind", "objective")).lower()
    if kind == "objective":
        objective = _function(runner.get("objective_fcn"), "runner.objective_fcn")
        lines.append(f"results=study.optimize(@{objective},{n_trials});")
    else:
        score = _function(runner.get("score_fcn"), "runner.score_fcn")
        if kind == "simulink":
            model = _quote(str(runner.get("model", "")))
            if not model:
                raise ValueError("runner.model is required")
            constructor = [f"ScoreFcn=@{score}"]
            for field, option in (
                ("configure_fcn", "ConfigureFcn"),
                ("constraint_fcn", "ConstraintFcn"),
                ("validation_fcn", "ValidationFcn"),
                ("result_fcn", "ResultFcn"),
                ("failure_classifier_fcn", "FailureClassifierFcn"),
            ):
                function = _optional_function(runner.get(field), f"runner.{field}")
                if function:
                    constructor.append(f"{option}=@{function}")
            stop_time = str(runner.get("stop_time", ""))
            if stop_time:
                constructor.append(f"StopTime=\"{_quote(stop_time)}\"")
            use_fast_restart = bool(runner.get("use_fast_restart", True))
            constructor.append(f"UseFastRestart={str(use_fast_restart).lower()}")
            context = runner.get("context", {})
            if not isinstance(context, Mapping):
                raise ValueError("runner.context must be a JSON object")
            context_json = _quote(json.dumps(context, ensure_ascii=False, allow_nan=False))
            constructor.append(f"Context=jsondecode('{context_json}')")
            lines.append(
                f"runner=radia.optuna.SimulinkRunner('{model}'," +
                ",".join(constructor) + ");"
            )
        else:
            configure = _function(runner.get("configure_fcn"), "runner.configure_fcn")
            netlist = _quote(str(runner.get("netlist", "")))
            if not netlist:
                raise ValueError("runner.netlist is required")
            lines.append(
                f"runner=radia.optuna.LTspiceRunner('{netlist}',ConfigureFcn=@{configure},ScoreFcn=@{score});"
            )
        method = "optimizeParallel" if parallel else "optimize"
        if kind == "simulink":
            continue_on_error = str(bool(runner.get("continue_on_error", True))).lower()
            method_options = [f"ContinueOnError={continue_on_error}"]
            if parallel:
                batch_size = int(runner.get("batch_size", 4))
                if batch_size <= 0:
                    raise ValueError("runner.batch_size must be positive")
                method_options.insert(0, f"BatchSize={batch_size}")
            lines.append(
                f"results=runner.{method}(study,{n_trials}," +
                ",".join(method_options) + ");"
            )
        else:
            lines.append(f"results=runner.{method}(study,{n_trials});")
    lines.append("if numel(study.Directions)>1, pareto=study.paretoFront(); else, best=study.bestTrial(); end")
    return "\n".join(lines) + "\n"


def _function(value, field):
    text = str(value or "")
    if not _FUNCTION.fullmatch(text):
        raise ValueError(f"{field} must be a qualified MATLAB function name")
    return text


def _optional_function(value, field):
    text = str(value or "")
    if text and not _FUNCTION.fullmatch(text):
        raise ValueError(f"{field} must be a qualified MATLAB function name")
    return text


def _quote(value):
    return str(value).replace("'", "''")


__all__ = ["matlab_optimize_build", "matlab_optimize_resume", "matlab_cad_topology_build", "matlab_sheet_metal_topology_build"]
