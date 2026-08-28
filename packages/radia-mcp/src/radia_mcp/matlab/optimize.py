"""Executable-manual helpers for MATLAB/Simulink/LTspice optimization."""

from __future__ import annotations

import json
import re
from typing import Any, Mapping

from .optuna_boundary import matlab_optuna_mcp_route


_FUNCTION = re.compile(r"^[A-Za-z]\w*(?:\.[A-Za-z]\w*)*$")
_SAMPLER_NAMES = {
    "bruteforce",
    "cmaes",
    "gp",
    "grid",
    "motpe",
    "nsgaii",
    "nsgaiii",
    "partial_fixed",
    "qmc",
    "random",
    "tpe",
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
    sampler = _sampler_spec(
        spec.get("sampler", "motpe" if len(directions) > 1 else "tpe"),
        directions,
    )
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
    code = _matlab_code(
        name, directions, sampler["matlab_code"], n_trials, runner,
        parallel, storage, live,
    )
    oracle_classification = "matlab-integration" if parallel else "upstream-python"
    oracle_reason = (
        "Parallel scheduling can change proposal consumption and is a MATLAB "
        "integration contract, even when the sampler's sequential seeded path is oracled."
        if parallel
        else "Sequential execution with the same explicit seed is eligible for the Optuna 4.9.0 differential oracle."
    )
    return {
        "schema": "radia-mcp.matlab-optimize-build/v3",
        "status": "ready",
        "runtime_owner": "MathWorks MATLAB MCP Server",
        "domain_owner": "radia-mcp.matlab.optimize",
        "mcp_ownership": matlab_optuna_mcp_route("composition")["route"],
        "execute_with": "official MATLAB MCP evaluate_matlab_code",
        "spec": {
            "name": name,
            "directions": directions,
            "sampler": sampler["normalized"],
            "n_trials": n_trials,
            "runner": runner,
            "parallel": parallel,
            "storage": storage,
            "live_monitor": live,
        },
        "matlab_code": code,
        "oracle": {
            "version": "4.9.0",
            "classification": oracle_classification,
            "reason": oracle_reason,
            "explicit_seed": sampler["seed"],
            "supported_surface": sampler["supported_surface"],
        },
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


def _matlab_code(name, directions, sampler_code, n_trials, runner, parallel, storage, live):
    direction_expr = "[" + ",".join(f'\"{item}\"' for item in directions) + "]"
    progress = "monitor=radia.optuna.LiveMonitor();\nprogress=@monitor.update;" if live else "progress=[];"
    storage_arg = f",StoragePath='{_quote(storage)}'" if storage else ""
    lines = [
        progress,
        f"study=radia.optuna.createStudy(study_name=\"{_quote_double(name)}\",directions={direction_expr},"
        f"sampler={sampler_code},ProgressFcn=progress{storage_arg});",
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
                constructor.append(f"StopTime=\"{_quote_double(stop_time)}\"")
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


def _sampler_spec(value: Any, directions: list[str], *, nested: bool = False) -> dict[str, Any]:
    if isinstance(value, str):
        options: dict[str, Any] = {"name": value}
    elif isinstance(value, Mapping):
        options = dict(value)
    else:
        raise ValueError("sampler must be a name or JSON object")
    name = str(options.pop("name", "")).lower()
    if name not in _SAMPLER_NAMES:
        raise ValueError(f"unknown sampler {name!r}; expected {sorted(_SAMPLER_NAMES)}")
    if nested and name == "partial_fixed":
        raise ValueError("nested partial_fixed samplers are not supported")
    if len(directions) > 1 and name == "cmaes":
        raise ValueError("CmaEsSampler supports only one objective")
    if len(directions) == 1 and name == "motpe":
        raise ValueError("MOTPESampler requires multiple objectives")

    seed = _int_option(options, "seed", 0, minimum=0)
    args: list[str] = []
    supported_surface = "seeded sequential proposal sequence"
    if name == "random":
        code = f"radia.optuna.RandomSampler({seed})"
    elif name in {"tpe", "motpe"}:
        class_name = "TPESampler" if name == "tpe" else "MOTPESampler"
        args = [f"Seed={seed}"]
        _append_int(args, options, "n_startup_trials", "NStartupTrials", 10, 0)
        _append_int(args, options, "n_ei_candidates", "NumberOfEIChoices", 24, 1)
        _append_float(args, options, "gamma", "Gamma", 0.1, 0.0, 1.0)
        _append_int(args, options, "max_good_trials", "MaxGoodTrials", 25, 1)
        _append_function(args, options, "gamma_fcn", "GammaFcn")
        _append_function(args, options, "weights_fcn", "WeightsFcn")
        _append_float(args, options, "prior_weight", "PriorWeight", 1.0, 0.0)
        _append_bool(args, options, "consider_magic_clip", "ConsiderMagicClip", True)
        _append_bool(args, options, "consider_endpoints", "ConsiderEndpoints", False)
        if name == "tpe":
            multivariate = _bool_option(options, "multivariate", False)
            group = _bool_option(options, "group", False)
            if group and not multivariate:
                raise ValueError("sampler.group requires sampler.multivariate=true")
            args.append(f"Multivariate={_matlab_bool(multivariate)}")
            args.append(f"Group={_matlab_bool(group)}")
            _append_bool(
                args, options, "warn_independent_sampling",
                "WarnIndependentSampling", False,
            )
            _append_bool(args, options, "constant_liar", "ConstantLiar", False)
        _append_function(args, options, "constraints_fcn", "ConstraintsFcn")
        _append_categorical_distances(args, options)
        code = f"radia.optuna.{class_name}({','.join(args)})"
        supported_surface = (
            "seeded scalar/mixed/grouped-multivariate/constrained TPE with callable gamma/weights and categorical distances"
            if name == "tpe"
            else "seeded multi-objective constrained TPE with callable gamma/weights and categorical distances"
        )
    elif name == "cmaes":
        args = [f"Seed={seed}"]
        _append_int(args, options, "n_startup_trials", "NStartupTrials", 1, 0)
        _append_int(args, options, "population_size", "PopulationSize", 0, 0)
        _append_float(args, options, "sigma0", "Sigma0", None, 0.0)
        _append_bool(args, options, "consider_pruned_trials", "ConsiderPrunedTrials", False)
        _append_bool(args, options, "warn_independent_sampling", "WarnIndependentSampling", True)
        if "independent_sampler" in options:
            independent = _sampler_spec(
                options.pop("independent_sampler"), directions, nested=True,
            )
            args.append(f"IndependentSampler={independent['matlab_code']}")
        if "x0" in options:
            x0 = options.pop("x0")
            if not isinstance(x0, Mapping):
                raise ValueError("sampler.x0 must be a JSON object")
            args.append(f"X0={_jsondecode(x0)}")
        code = f"radia.optuna.CmaEsSampler({','.join(args)})"
        supported_surface = "seeded numeric CMA-ES with independent-sampler controls"
    elif name == "gp":
        args = [f"Seed={seed}"]
        _append_int(args, options, "n_startup_trials", "NStartupTrials", 10, 0)
        _append_bool(args, options, "deterministic_objective", "DeterministicObjective", False)
        backend = str(options.pop("backend", "upstream-python"))
        if backend not in {"upstream-python", "matlab-native"}:
            raise ValueError("sampler.backend must be upstream-python or matlab-native")
        args.append(f'Backend="{backend}"')
        # Bounds mirror GPSampler.m: CandidateCount is mustBePositive and
        # LocalSearchCount is mustBeNonnegative.  Tighter limits here rejected
        # values MATLAB itself accepts.
        _append_int(args, options, "candidate_count", "CandidateCount", None, 1)
        _append_int(args, options, "local_search_count", "LocalSearchCount", None, 0)
        _append_int(args, options, "monte_carlo_samples", "MonteCarloSamples", None, 1)
        _append_function(args, options, "constraints_fcn", "ConstraintsFcn")
        code = f"radia.optuna.GPSampler({','.join(args)})"
        supported_surface = (
            "exact pinned Optuna 4.9.0 GP including post-startup acquisition"
            if backend == "upstream-python"
            else "MATLAB-native GP integration only"
        )
    elif name in {"nsgaii", "nsgaiii"}:
        class_name = "NSGAIISampler" if name == "nsgaii" else "NSGAIIISampler"
        args = [f"Seed={seed}"]
        _append_int(args, options, "population_size", "PopulationSize", 50, 2)
        _append_float(args, options, "mutation_probability", "MutationProbability", None, 0.0, 1.0)
        _append_float(args, options, "crossover_probability", "CrossoverProbability", 0.9, 0.0, 1.0)
        _append_float(args, options, "swapping_probability", "SwappingProbability", 0.5, 0.0, 1.0)
        _append_function(args, options, "constraints_fcn", "ConstraintsFcn")
        if name == "nsgaiii":
            _append_int(args, options, "dividing_parameter", "DividingParameter", 3, 1)
            if "reference_points" in options:
                points = options.pop("reference_points")
                if not isinstance(points, list):
                    raise ValueError("sampler.reference_points must be a JSON array")
                args.append(f"ReferencePoints={_matlab_matrix(points)}")
        if "crossover" in options:
            args.append(f"Crossover={_crossover(options.pop('crossover'))}")
        code = f"radia.optuna.{class_name}({','.join(args)})"
        supported_surface = "seeded NSGA-II/III population and crossover behavior"
    elif name == "qmc":
        qmc_type = str(options.pop("qmc_type", "sobol")).lower()
        if qmc_type not in {"sobol", "halton"}:
            raise ValueError("sampler.qmc_type must be sobol or halton")
        scramble = _bool_option(options, "scramble", False)
        code = (
            "radia.optuna.QMCSampler("
            f'QMCType="{qmc_type}",Scramble={_matlab_bool(scramble)},Seed={seed})'
        )
        supported_surface = "seeded scrambled or deterministic unscrambled QMC"
    elif name == "grid":
        search_space = options.pop("search_space", None)
        if not isinstance(search_space, Mapping) or not search_space:
            raise ValueError("GridSampler requires non-empty sampler.search_space")
        code = (
            "radia.optuna.GridSampler("
            f"{_jsondecode(search_space)},Seed={seed})"
        )
        supported_surface = "finite Cartesian grid exhaustion"
    elif name == "bruteforce":
        avoid = _bool_option(options, "avoid_premature_stop", False)
        code = (
            "radia.optuna.BruteForceSampler("
            f"Seed={seed},AvoidPrematureStop={_matlab_bool(avoid)})"
        )
        supported_surface = "fixed and conditional define-by-run exhaustion"
    else:
        fixed = options.pop("fixed_params", None)
        if not isinstance(fixed, Mapping) or not fixed:
            raise ValueError("PartialFixedSampler requires non-empty sampler.fixed_params")
        base = _sampler_spec(options.pop("base_sampler", "tpe"), directions, nested=True)
        code = (
            "radia.optuna.PartialFixedSampler("
            f"{_jsondecode(fixed)},{base['matlab_code']})"
        )
        supported_surface = "fixed parameters plus an oracled base sampler"

    if options:
        raise ValueError(f"unsupported {name} sampler options: {sorted(options)}")
    normalized = {"name": name, "seed": seed, "matlab_code": code}
    return {
        "name": name,
        "seed": seed,
        "matlab_code": code,
        "normalized": normalized,
        "supported_surface": supported_surface,
    }


def _int_option(options, key, default, *, minimum=None, maximum=None):
    raw = options.pop(key, default)
    if raw is None:
        return None
    # A non-numeric value used to surface as int()'s own "invalid literal"
    # message, which names neither the sampler option nor the expected type.
    try:
        rejected = isinstance(raw, bool) or int(raw) != raw
    except (TypeError, ValueError):
        rejected = True
    if rejected:
        raise ValueError(f"sampler.{key} must be an integer, got {raw!r}")
    value = int(raw)
    if minimum is not None and value < minimum:
        raise ValueError(f"sampler.{key} must be >= {minimum}")
    if maximum is not None and value > maximum:
        raise ValueError(f"sampler.{key} must be <= {maximum}")
    return value


def _float_option(options, key, default, *, minimum=None, maximum=None):
    raw = options.pop(key, default)
    if raw is None:
        return None
    value = float(raw)
    if minimum is not None and value < minimum:
        raise ValueError(f"sampler.{key} must be >= {minimum}")
    if maximum is not None and value > maximum:
        raise ValueError(f"sampler.{key} must be <= {maximum}")
    return value


def _bool_option(options, key, default):
    raw = options.pop(key, default)
    if not isinstance(raw, bool):
        raise ValueError(f"sampler.{key} must be boolean")
    return raw


def _append_int(args, options, key, matlab_name, default, minimum=0, maximum=None):
    value = _int_option(options, key, default, minimum=minimum, maximum=maximum)
    if value is not None:
        args.append(f"{matlab_name}={value}")


def _append_float(args, options, key, matlab_name, default, minimum=None, maximum=None):
    value = _float_option(options, key, default, minimum=minimum, maximum=maximum)
    if value is not None:
        args.append(f"{matlab_name}={value:.17g}")


def _append_bool(args, options, key, matlab_name, default):
    value = _bool_option(options, key, default)
    args.append(f"{matlab_name}={_matlab_bool(value)}")


def _append_function(args, options, key, matlab_name):
    if key not in options:
        return
    value = _optional_function(options.pop(key), f"sampler.{key}")
    if value:
        args.append(f"{matlab_name}=@{value}")


def _append_categorical_distances(args, options):
    key = "categorical_distance_fcn"
    if key not in options:
        return
    mapping = options.pop(key)
    if not isinstance(mapping, Mapping):
        raise ValueError(f"sampler.{key} must be a JSON object")
    if not mapping:
        return
    parameter_names: list[str] = []
    functions: list[str] = []
    for parameter_name in sorted(mapping, key=str):
        name = str(parameter_name)
        function = _function(mapping[parameter_name], f"sampler.{key}.{name}")
        parameter_names.append(f"'{_quote(name)}'")
        functions.append(f"@{function}")
    args.append(
        "CategoricalDistanceFcn=containers.Map(" +
        "{" + ",".join(parameter_names) + "}," +
        "{" + ",".join(functions) + "})"
    )


def _matlab_bool(value):
    return "true" if value else "false"


def _jsondecode(value):
    encoded = json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    return f"jsondecode('{_quote(encoded)}')"


def _matlab_matrix(value):
    if not value or not all(isinstance(row, list) and row for row in value):
        raise ValueError("reference_points must be a non-empty rectangular matrix")
    width = len(value[0])
    if any(len(row) != width for row in value):
        raise ValueError("reference_points must be rectangular")
    return "[" + ";".join(
        " ".join(f"{float(item):.17g}" for item in row) for row in value
    ) + "]"


def _crossover(value):
    name = str(value).lower()
    classes = {
        "uniform": "UniformCrossover",
        "blxalpha": "BLXAlphaCrossover",
        "sbx": "SBXCrossover",
        "vsbx": "VSBXCrossover",
        "spx": "SPXCrossover",
        "undx": "UNDXCrossover",
    }
    if name not in classes:
        raise ValueError(f"unknown NSGA-II crossover {name!r}; expected {sorted(classes)}")
    return f"radia.optuna.nsgaii.{classes[name]}()"


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
    """Escape for a MATLAB single-quoted char array ('...')."""
    return str(value).replace("'", "''")


def _quote_double(value):
    """Escape for a MATLAB double-quoted string ("...").

    MATLAB doubles the delimiter to escape it, and the delimiter differs
    between the two literal forms.  Using the single-quote escaper inside a
    double-quoted string left `"` unescaped (breaking the generated code) and
    turned `'` into `''` (silently corrupting the value).
    """
    return str(value).replace('"', '""')


__all__ = ["matlab_optimize_build", "matlab_optimize_resume", "matlab_cad_topology_build", "matlab_sheet_metal_topology_build"]
