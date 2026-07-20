"""Executable-manual helpers for MATLAB/Simulink/LTspice optimization."""

from __future__ import annotations

import json
import re
from typing import Any, Mapping


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
        "schema": "radia-mcp.matlab-optimize-build/v1",
        "status": "ready",
        "runtime_owner": "MathWorks MATLAB MCP Server",
        "domain_owner": "radia-mcp.matlab.optimize",
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
    """Generate one Radia-VIM/LP sheet update and adaptive mesh-routing step."""
    if isinstance(spec,str): spec=json.loads(spec)
    if not isinstance(spec,Mapping): raise ValueError("sheet-metal topology spec must be a JSON object")
    design_data=str(spec.get("design_data",""))
    if not design_data: raise ValueError("design_data MAT-file is required")
    linearize=_function(spec.get("linearize_fcn"),"linearize_fcn")
    volume_max=float(spec.get("volume_max",0)); displacement_move=float(spec.get("displacement_move",0))
    thickness_move=float(spec.get("thickness_move",0)); bounds=spec.get("thickness_bounds")
    if volume_max<=0 or displacement_move<=0 or thickness_move<=0: raise ValueError("positive volume_max and move limits are required")
    if not isinstance(bounds,list) or len(bounds)!=2 or not 0<float(bounds[0])<=float(bounds[1]): raise ValueError("thickness_bounds must be [positive_min,max]")
    code=(
        f"data=load('{_quote(design_data)}');\n"
        f"model={linearize}(data.normalDisplacement,data.thickness,data.activation);\n"
        "update=radia.topopt.solveSheetMetalLP(data.normalDisplacement,data.thickness,data.activation,model.objective_gradient,data.cellAreas,"
        f"VolumeMax={volume_max:.17g},DisplacementMove={displacement_move:.17g},ThicknessMove={thickness_move:.17g},"
        f"ThicknessBounds=[{float(bounds[0]):.17g},{float(bounds[1]):.17g}],Laplacian=data.laplacian,CurvatureLimit=data.curvatureLimit);\n"
        "acceptance=radia.topopt.acceptTrafoStep(model.quality_fcn,model.relative_displacements);\n"
        "meshRoute=acceptance.route;\n"
    )
    return {
        "schema":"radia-mcp.matlab-sheet-metal-topology-build/v1","status":"ready",
        "runtime_owner":"MathWorks MATLAB MCP Server","domain_owner":"radia-mcp.matlab.optimize",
        "method":"Radia-VIM linearization -> local trust-region LP -> Trafo quality backtracking -> NGSolve deform/refine or Cubit rebuild",
        "matlab_code":code,"execute_with":"official MATLAB MCP evaluate_matlab_code",
        "mesh_routes":["ngsolve_deform","ngsolve_refine","cubit_rebuild"],
        "design_data_contract":{"normalDisplacement":"n-by-1","thickness":"n-by-1","activation":"n-by-1", "cellAreas":"n-by-1","laplacian":"k-by-n","curvatureLimit":"scalar or k-by-1"},
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
        configure = _function(runner.get("configure_fcn"), "runner.configure_fcn")
        score = _function(runner.get("score_fcn"), "runner.score_fcn")
        if kind == "simulink":
            model = _quote(str(runner.get("model", "")))
            if not model:
                raise ValueError("runner.model is required")
            lines.append(
                f"runner=radia.optuna.SimulinkRunner('{model}',ConfigureFcn=@{configure},ScoreFcn=@{score});"
            )
        else:
            netlist = _quote(str(runner.get("netlist", "")))
            if not netlist:
                raise ValueError("runner.netlist is required")
            lines.append(
                f"runner=radia.optuna.LTspiceRunner('{netlist}',ConfigureFcn=@{configure},ScoreFcn=@{score});"
            )
        method = "optimizeParallel" if parallel else "optimize"
        lines.append(f"results=runner.{method}(study,{n_trials});")
    lines.append("if numel(study.Directions)>1, pareto=study.paretoFront(); else, best=study.bestTrial(); end")
    return "\n".join(lines) + "\n"


def _function(value, field):
    text = str(value or "")
    if not _FUNCTION.fullmatch(text):
        raise ValueError(f"{field} must be a qualified MATLAB function name")
    return text


def _quote(value):
    return str(value).replace("'", "''")


__all__ = ["matlab_optimize_build", "matlab_optimize_resume", "matlab_cad_topology_build", "matlab_sheet_metal_topology_build"]
