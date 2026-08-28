"""Compare the Optuna 4.9.0 public inventory with the MATLAB package surface."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
INVENTORY_PATH = Path(__file__).with_name("optuna49_public_api.json")
ORACLE_PATH = Path(__file__).with_name("optuna49_oracle.json")
MATLAB_DIRECTORY = ROOT / "matlab" / "+radia" / "+optuna"
DESTINATION = ROOT / "matlab" / "optuna49_api_coverage.json"
FUNCTION_PATTERN = re.compile(
    r"(?m)^\s*function\s+(?:\[[^]]*\]|\w+)\s*=\s*(\w+)\s*\(|"
    r"^\s*function\s+(\w+)\s*\("
)
ABSTRACT_METHOD_PATTERN = re.compile(
    r"(?m)^\s*(?:\[[^]]*\]|\w+)\s*=\s*(\w+)\s*\(|"
    r"^\s*(\w+)\s*\("
)


def _class_blocks(source: str, keyword: str) -> list[tuple[str, str]]:
    pattern = re.compile(
        rf"(?ms)^(?P<indent>[ \t]*){keyword}(?P<qualifier>[^\r\n]*)\r?\n"
        rf"(?P<body>.*?)(?=^(?P=indent)end[ \t]*$)"
    )
    return [
        (match.group("qualifier"), match.group("body"))
        for match in pattern.finditer(source)
    ]


def _public_get_access(qualifier: str) -> bool:
    compact = re.sub(r"\s+", "", qualifier).lower()
    return (
        re.search(
            r"(?:^|[(,])(?:access|getaccess)=(?:private|protected)(?:[,)]|$)",
            compact,
        )
        is None
    )
MODULE_EQUIVALENTS = {
    "artifacts": True,
    "distributions": True,
    "exceptions": True,
    "importance": True,
    "integration": True,
    "logging": True,
    "matplotlib": True,
    "nsgaii": True,
    "pruners": True,
    "samplers": True,
    "search_space": True,
    "storages": True,
    "study": True,
    "trial": True,
    "version": True,
    "visualization": True,
}
VERIFIED_FUNCTIONS = {
    "copy_study",
    "check_distribution_compatibility",
    "create_study",
    "create_trial",
    "delete_study",
    "get_all_study_names",
    "get_all_study_summaries",
    "get_param_importances",
    "intersection_search_space",
    "distribution_to_json",
    "fail_stale_trials",
    "json_to_distribution",
    "load_study",
    "is_available",
    "plot_contour",
    "plot_edf",
    "plot_hypervolume_history",
    "plot_intermediate_values",
    "plot_optimization_history",
    "plot_parallel_coordinate",
    "plot_param_importances",
    "plot_pareto_front",
    "plot_rank",
    "plot_slice",
    "plot_terminator_improvement",
    "plot_timeline",
    "run_grpc_proxy_server",
}
CLASS_EQUIVALENTS = {
    "_CachedStorage": "CachedStorage",
    "_GroupDecomposedSearchSpace": "GroupDecomposedSearchSpace",
    "_SearchSpaceGroup": "SearchSpaceGroup",
}
SYMBOL_EQUIVALENTS = {
    "__version__": "version",
}
INTEGRATION_EXPORTS = {
    "AllenNLPExecutor",
    "AllenNLPPruningCallback",
    "BoTorchSampler",
    "CatBoostPruningCallback",
    "ChainerMNStudy",
    "ChainerPruningExtension",
    "DaskStorage",
    "FastAIPruningCallback",
    "FastAIV2PruningCallback",
    "KerasPruningCallback",
    "LightGBMPruningCallback",
    "LightGBMTuner",
    "LightGBMTunerCV",
    "MLflowCallback",
    "MXNetPruningCallback",
    "OptunaSearchCV",
    "PyCmaSampler",
    "PyTorchIgnitePruningHandler",
    "PyTorchLightningPruningCallback",
    "ShapleyImportanceEvaluator",
    "SkorchPruningCallback",
    "TFKerasPruningCallback",
    "TensorBoardCallback",
    "TensorFlowPruningHook",
    "TorchDistributedTrial",
    "WeightsAndBiasesCallback",
    "XGBoostPruningCallback",
}
VERIFIED_SYMBOLS = VERIFIED_FUNCTIONS | INTEGRATION_EXPORTS | {
    "CategoricalChoiceType",
    "DISTRIBUTION_CLASSES",
    "CRITICAL",
    "DEBUG",
    "ERROR",
    "FATAL",
    "INFO",
    "WARN",
    "WARNING",
    "create_default_formatter",
    "disable_default_handler",
    "disable_propagation",
    "enable_default_handler",
    "enable_propagation",
    "get_logger",
    "get_verbosity",
    "set_verbosity",
    "report_cross_validation_scores",
    "download_artifact",
    "get_all_artifact_meta",
    "upload_artifact",
    "__version__",
    "CLIUsageError",
    "DuplicatedStudyError",
    "ExperimentalWarning",
    "OptunaError",
    "StorageInternalError",
    "TrialPruned",
    "UpdateFinishedTrialError",
    "BaseErrorEvaluator",
    "BaseImprovementEvaluator",
    "BaseTerminator",
    "CrossValidationErrorEvaluator",
    "EMMREvaluator",
    "MedianErrorEvaluator",
    "RegretBoundEvaluator",
    "StaticErrorEvaluator",
    "ArtifactMeta",
    "Backoff",
    "Boto3ArtifactStore",
    "FileSystemArtifactStore",
    "GCSArtifactStore",
    "GrpcStorageProxy",
    "BaseGASampler",
    "BaseStorage",
    "BaseJournalBackend",
    "BaseJournalLogStorage",
    "InMemoryStorage",
    "JournalFileBackend",
    "JournalFileOpenLock",
    "JournalFileStorage",
    "JournalFileSymlinkLock",
    "JournalRedisBackend",
    "JournalRedisStorage",
    "JournalStorage",
    "RDBStorage",
    "RetryFailedTrialCallback",
    "RetryHeartbeatStaleTrialCallback",
    "_CachedStorage",
}
SAMPLER_PUBLIC_MEMBERS = {
    "after_trial",
    "before_trial",
    "infer_relative_search_space",
    "sample_independent",
    "sample_relative",
}
GA_PUBLIC_MEMBERS = {
    "get_parent_population",
    "get_population",
    "get_trial_generation",
    "population_size",
    "select_parent",
}
STORAGE_PUBLIC_MEMBERS = {
    "check_trial_is_updatable",
    "create_new_study",
    "create_new_trial",
    "delete_study",
    "get_all_studies",
    "get_all_trials",
    "get_best_trial",
    "get_n_trials",
    "get_study_directions",
    "get_study_id_from_name",
    "get_study_name_from_id",
    "get_study_system_attrs",
    "get_study_user_attrs",
    "get_trial",
    "get_trial_id_from_study_id_trial_number",
    "get_trial_number_from_id",
    "get_trial_param",
    "get_trial_params",
    "get_trial_system_attrs",
    "get_trial_user_attrs",
    "remove_session",
    "set_study_system_attr",
    "set_study_user_attr",
    "set_trial_intermediate_value",
    "set_trial_param",
    "set_trial_state_values",
    "set_trial_system_attr",
    "set_trial_user_attr",
}
CROSSOVER_PUBLIC_MEMBERS = {"crossover", "n_parents"}
EXCEPTION_PUBLIC_MEMBERS = {"add_note", "args", "with_traceback"}
VERIFIED_MEMBERS = {
    "_GroupDecomposedSearchSpace": {"calculate"},
    "_SearchSpaceGroup": {"add_distributions", "search_spaces"},
    "BaseImportanceEvaluator": {"evaluate"},
    "BaseGASampler": SAMPLER_PUBLIC_MEMBERS
    | GA_PUBLIC_MEMBERS
    | {"reseed_rng"},
    "BaseStorage": STORAGE_PUBLIC_MEMBERS,
    "BaseJournalBackend": {"append_logs", "read_logs"},
    "BaseJournalLogStorage": {"append_logs", "read_logs"},
    "InMemoryStorage": STORAGE_PUBLIC_MEMBERS,
    "JournalFileBackend": {"append_logs", "read_logs"},
    "JournalFileOpenLock": {"acquire", "release"},
    "JournalFileStorage": {"append_logs", "read_logs"},
    "JournalFileSymlinkLock": {"acquire", "release"},
    "JournalRedisBackend": {
        "append_logs",
        "load_snapshot",
        "read_logs",
        "save_snapshot",
    },
    "JournalRedisStorage": {
        "append_logs",
        "load_snapshot",
        "read_logs",
        "save_snapshot",
    },
    "JournalStorage": STORAGE_PUBLIC_MEMBERS | {"restore_replay_result"},
    "RDBStorage": STORAGE_PUBLIC_MEMBERS
    | {
        "get_all_versions",
        "get_current_version",
        "get_failed_trial_callback",
        "get_head_version",
        "get_heartbeat_interval",
        "get_heartbeat_stale_trial_callback",
        "record_heartbeat",
        "upgrade",
    },
    "_CachedStorage": STORAGE_PUBLIC_MEMBERS
    | {
        "get_failed_trial_callback",
        "get_heartbeat_interval",
        "get_heartbeat_stale_trial_callback",
        "record_heartbeat",
    },
    "RetryFailedTrialCallback": {
        "retried_trial_number",
        "retry_history",
    },
    "RetryHeartbeatStaleTrialCallback": {
        "retried_trial_number",
        "retry_history",
    },
    "ArtifactMeta": {"artifact_id", "encoding", "filename", "mimetype"},
    "Backoff": {"open_reader", "remove", "write"},
    "Boto3ArtifactStore": {"open_reader", "remove", "write"},
    "BaseErrorEvaluator": {"evaluate"},
    "BaseImprovementEvaluator": {"evaluate"},
    "BaseTerminator": {"should_terminate"},
    "CLIUsageError": EXCEPTION_PUBLIC_MEMBERS,
    "DuplicatedStudyError": EXCEPTION_PUBLIC_MEMBERS,
    "ExperimentalWarning": EXCEPTION_PUBLIC_MEMBERS,
    "CrossValidationErrorEvaluator": {"evaluate"},
    "EMMREvaluator": {"evaluate"},
    "FanovaImportanceEvaluator": {"evaluate"},
    "FileSystemArtifactStore": {"open_reader", "remove", "write"},
    "GCSArtifactStore": {"open_reader", "remove", "write"},
    "GrpcStorageProxy": STORAGE_PUBLIC_MEMBERS
    | {"close", "wait_server_ready"},
    "MeanDecreaseImpurityImportanceEvaluator": {"evaluate"},
    "MedianErrorEvaluator": {"evaluate"},
    "RegretBoundEvaluator": {"evaluate"},
    "OptunaError": EXCEPTION_PUBLIC_MEMBERS,
    "PedAnovaImportanceEvaluator": {"evaluate"},
    "StorageInternalError": EXCEPTION_PUBLIC_MEMBERS,
    "StaticErrorEvaluator": {"evaluate"},
    "TrialPruned": EXCEPTION_PUBLIC_MEMBERS,
    "UpdateFinishedTrialError": EXCEPTION_PUBLIC_MEMBERS,
    "BaseDistribution": {
        "single",
        "to_external_repr",
        "to_internal_repr",
    },
    "CategoricalDistribution": {
        "single",
        "to_external_repr",
        "to_internal_repr",
    },
    "FloatDistribution": {
        "single",
        "to_external_repr",
        "to_internal_repr",
    },
    "IntDistribution": {
        "single",
        "to_external_repr",
        "to_internal_repr",
    },
    "UniformDistribution": {
        "single",
        "to_external_repr",
        "to_internal_repr",
    },
    "LogUniformDistribution": {
        "single",
        "to_external_repr",
        "to_internal_repr",
    },
    "DiscreteUniformDistribution": {
        "q",
        "single",
        "to_external_repr",
        "to_internal_repr",
    },
    "IntUniformDistribution": {
        "single",
        "to_external_repr",
        "to_internal_repr",
    },
    "IntLogUniformDistribution": {
        "single",
        "to_external_repr",
        "to_internal_repr",
    },
    "BasePruner": {"prune"},
    "BaseSampler": SAMPLER_PUBLIC_MEMBERS | {"reseed_rng"},
    "BaseTrial": {
        "datetime_start",
        "distributions",
        "number",
        "params",
        "report",
        "set_system_attr",
        "set_user_attr",
        "should_prune",
        "suggest_categorical",
        "suggest_discrete_uniform",
        "suggest_float",
        "suggest_int",
        "suggest_loguniform",
        "suggest_uniform",
        "system_attrs",
        "user_attrs",
    },
    "Study": {
        "add_trial",
        "add_trials",
        "ask",
        "best_params",
        "best_trial",
        "best_trials",
        "best_value",
        "direction",
        "directions",
        "enqueue_trial",
        "get_trials",
        "metric_names",
        "optimize",
        "set_metric_names",
        "set_system_attr",
        "set_user_attr",
        "stop",
        "system_attrs",
        "tell",
        "trials",
        "trials_dataframe",
        "user_attrs",
    },
    "StudySummary": {"direction", "directions", "system_attrs"},
    "Trial": {
        "datetime_start",
        "distributions",
        "number",
        "params",
        "report",
        "relative_params",
        "set_system_attr",
        "set_user_attr",
        "should_prune",
        "suggest_categorical",
        "suggest_discrete_uniform",
        "suggest_float",
        "suggest_int",
        "suggest_loguniform",
        "suggest_uniform",
        "system_attrs",
        "user_attrs",
    },
    "FixedTrial": {
        "datetime_start",
        "distributions",
        "number",
        "params",
        "report",
        "set_system_attr",
        "set_user_attr",
        "should_prune",
        "suggest_categorical",
        "suggest_discrete_uniform",
        "suggest_float",
        "suggest_int",
        "suggest_loguniform",
        "suggest_uniform",
        "system_attrs",
        "user_attrs",
    },
    "FrozenTrial": {
        "datetime_complete",
        "datetime_start",
        "distributions",
        "duration",
        "intermediate_values",
        "last_step",
        "number",
        "params",
        "report",
        "set_system_attr",
        "set_user_attr",
        "should_prune",
        "state",
        "suggest_categorical",
        "suggest_discrete_uniform",
        "suggest_float",
        "suggest_int",
        "suggest_loguniform",
        "suggest_uniform",
        "system_attrs",
        "user_attrs",
        "value",
        "values",
    },
    "BestValueStagnationEvaluator": {"evaluate"},
    "Terminator": {"should_terminate"},
    "IntersectionSearchSpace": {"calculate"},
    "BruteForceSampler": SAMPLER_PUBLIC_MEMBERS | {"reseed_rng"},
    "CmaEsSampler": SAMPLER_PUBLIC_MEMBERS | {"reseed_rng"},
    "GPSampler": SAMPLER_PUBLIC_MEMBERS | {"reseed_rng"},
    "GridSampler": SAMPLER_PUBLIC_MEMBERS | {"is_exhausted", "reseed_rng"},
    "NSGAIIISampler": SAMPLER_PUBLIC_MEMBERS | GA_PUBLIC_MEMBERS | {"reseed_rng"},
    "NSGAIISampler": SAMPLER_PUBLIC_MEMBERS | GA_PUBLIC_MEMBERS | {"reseed_rng"},
    "PartialFixedSampler": SAMPLER_PUBLIC_MEMBERS | {"reseed_rng"},
    "QMCSampler": SAMPLER_PUBLIC_MEMBERS | {"reseed_rng"},
    "RandomSampler": SAMPLER_PUBLIC_MEMBERS | {"reseed_rng"},
    "TPESampler": SAMPLER_PUBLIC_MEMBERS | {"hyperopt_parameters", "reseed_rng"},
    "BLXAlphaCrossover": CROSSOVER_PUBLIC_MEMBERS,
    "BaseCrossover": CROSSOVER_PUBLIC_MEMBERS,
    "SBXCrossover": CROSSOVER_PUBLIC_MEMBERS,
    "SPXCrossover": CROSSOVER_PUBLIC_MEMBERS,
    "UNDXCrossover": CROSSOVER_PUBLIC_MEMBERS,
    "UniformCrossover": CROSSOVER_PUBLIC_MEMBERS,
    "VSBXCrossover": CROSSOVER_PUBLIC_MEMBERS,
    "HyperbandPruner": {"prune"},
    "MedianPruner": {"prune"},
    "NopPruner": {"prune"},
    "PatientPruner": {"prune"},
    "PercentilePruner": {"prune"},
    "SuccessiveHalvingPruner": {"prune"},
    "ThresholdPruner": {"prune"},
    "WilcoxonPruner": {"prune"},
    "StudyDirection": {
        "NOT_SET",
        "MINIMIZE",
        "MAXIMIZE",
        "as_integer_ratio",
        "bit_count",
        "bit_length",
        "conjugate",
        "denominator",
        "from_bytes",
        "imag",
        "is_integer",
        "name",
        "numerator",
        "real",
        "to_bytes",
        "value",
    },
    "TrialState": {
        "RUNNING",
        "COMPLETE",
        "PRUNED",
        "FAIL",
        "WAITING",
        "is_finished",
        "as_integer_ratio",
        "bit_count",
        "bit_length",
        "conjugate",
        "denominator",
        "from_bytes",
        "imag",
        "is_integer",
        "name",
        "numerator",
        "real",
        "to_bytes",
        "value",
    },
}
CLASS_ORACLE_SECTIONS = {
    "_GroupDecomposedSearchSpace": ("search_space",),
    "_SearchSpaceGroup": ("search_space",),
    "BaseImportanceEvaluator": ("base_components", "importance"),
    "FanovaImportanceEvaluator": ("importance",),
    "MeanDecreaseImpurityImportanceEvaluator": ("importance",),
    "PedAnovaImportanceEvaluator": ("importance",),
    "BaseDistribution": ("distributions", "distribution_public_members"),
    "CategoricalDistribution": ("distributions", "distribution_public_members"),
    "DiscreteUniformDistribution": ("distributions", "distribution_public_members"),
    "FloatDistribution": ("distributions", "distribution_public_members"),
    "IntDistribution": ("distributions", "distribution_public_members"),
    "IntLogUniformDistribution": ("distributions", "distribution_public_members"),
    "IntUniformDistribution": ("distributions", "distribution_public_members"),
    "LogUniformDistribution": ("distributions", "distribution_public_members"),
    "UniformDistribution": ("distributions", "distribution_public_members"),
    "BasePruner": ("base_components", "pruners"),
    "HyperbandPruner": ("pruners",),
    "MedianPruner": ("pruners",),
    "NopPruner": ("pruners",),
    "PatientPruner": ("pruners",),
    "PercentilePruner": ("pruners",),
    "SuccessiveHalvingPruner": ("pruners",),
    "ThresholdPruner": ("pruners",),
    "WilcoxonPruner": ("pruners",),
    "BaseSampler": ("base_components", "sampler_public_members", "sampler_reseed"),
    "BruteForceSampler": (
        "brute_force_sampler_seed_29",
        "conditional_brute_force_sampler_seed_79",
        "sampler_public_members",
        "sampler_reseed",
    ),
    "CmaEsSampler": (
        "cmaes_independent_sampler_seed_31",
        "cmaes_sampler_seed_31",
        "cmaes_advanced",
        "sampler_public_members",
        "sampler_reseed",
    ),
    "GPSampler": (
        "gp_constraints_sampler_seed_89",
        "gp_sampler_seed_53",
        "sampler_public_members",
        "sampler_reseed",
    ),
    "GridSampler": ("grid_sampler_seed_17", "sampler_public_members", "sampler_reseed"),
    "NSGAIISampler": (
        "nsgaii_sampler_seed_19",
        "sampler_public_members",
        "sampler_reseed",
    ),
    "NSGAIIISampler": (
        "nsgaiii_sampler_seed_23",
        "sampler_public_members",
        "sampler_reseed",
    ),
    "PartialFixedSampler": (
        "partial_fixed_sampler_seed_61",
        "sampler_public_members",
        "sampler_reseed",
    ),
    "QMCSampler": (
        "native_sobol_high_dimension",
        "qmc_warnings",
        "sampler_public_members",
        "sampler_reseed",
        "scrambled_qmc_sampler_seed_47",
        "unscrambled_qmc_sampler_seed_71",
    ),
    "RandomSampler": ("random_sampler_seed_123", "sampler_public_members", "sampler_reseed"),
    "TPESampler": (
        "custom_tpe_sampler_gamma_weights",
        "mixed_tpe_sampler_seed_43",
        "multiobjective_tpe_sampler_seed_41",
        "multivariate_tpe_sampler_seed_67",
        "sampler_public_members",
        "sampler_reseed",
        "tpe_categorical_distance",
        "tpe_group",
        "tpe_sampler_seed_37",
        "tpe_constant_liar_seed_127",
    ),
    "BaseCrossover": ("base_components", "nsgaii_crossovers_seed_73"),
    "BLXAlphaCrossover": ("nsgaii_crossovers_seed_73",),
    "SBXCrossover": ("nsgaii_crossovers_seed_73",),
    "SPXCrossover": ("nsgaii_crossovers_seed_73",),
    "UNDXCrossover": ("nsgaii_crossovers_seed_73",),
    "UniformCrossover": ("nsgaii_crossovers_seed_73",),
    "VSBXCrossover": ("nsgaii_crossovers_seed_73",),
    "IntersectionSearchSpace": ("search_space",),
    "BestValueStagnationEvaluator": ("terminator",),
    "MaxTrialsCallback": ("terminator",),
    "Terminator": ("terminator",),
    "TerminatorCallback": ("terminator",),
    "Study": ("core_api", "study_management", "tell", "trials_dataframe"),
    "StudyDirection": ("enums",),
    "StudySummary": ("study_management",),
    "BaseTrial": ("base_components", "base_trial"),
    "FixedTrial": ("fixed_trial",),
    "FrozenTrial": ("frozen_trial",),
    "Trial": ("base_trial", "core_api", "tell"),
    "TrialState": ("enums",),
}


def _sha256(path: Path) -> str:
    content = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(content).hexdigest().upper()


def _normalized(name: str) -> str:
    return re.sub(r"_", "", name).lower()


def _matlab_surface() -> tuple[set[str], dict[str, dict[str, set[str]]]]:
    names: set[str] = set()
    members: dict[str, dict[str, set[str]]] = {}
    parents: dict[str, str] = {}
    for path in sorted(MATLAB_DIRECTORY.rglob("*.m")):
        if "+internal" in path.parts:
            continue
        names.add(path.stem)
        source = path.read_text(encoding="utf-8")
        class_match = re.search(
            r"(?m)^\s*classdef(?:\s*\([^)]*\))?\s+(\w+)"
            r"(?:\s*<\s*([\w.]+))?",
            source,
        )
        if class_match and class_match.group(2):
            parents[class_match.group(1)] = class_match.group(2).split(".")[-1]
        found: set[str] = set()
        for qualifier, block in _class_blocks(source, "methods"):
            if not _public_get_access(qualifier):
                continue
            found.update(
                candidate
                for groups in FUNCTION_PATTERN.findall(block)
                for candidate in groups
                if candidate
            )
            if "abstract" in qualifier.lower():
                found.update(
                    candidate
                    for groups in ABSTRACT_METHOD_PATTERN.findall(block)
                    for candidate in groups
                    if candidate
                )
        if "classdef" not in source:
            found.update(
                candidate
                for groups in FUNCTION_PATTERN.findall(source)
                for candidate in groups
                if candidate
            )
        class_members = members.setdefault(path.stem, {})
        for name in found:
            class_members.setdefault(_normalized(name), set()).add(name)
        for qualifier, block in _class_blocks(source, "properties"):
            if not _public_get_access(qualifier):
                continue
            for line in block.splitlines():
                match = re.match(r"\s*([A-Za-z]\w*)\b", line)
                if match and not line.lstrip().startswith("%"):
                    name = match.group(1)
                    class_members.setdefault(_normalized(name), set()).add(name)
        for qualifier, block in _class_blocks(source, "enumeration"):
            if not _public_get_access(qualifier):
                continue
            for line in block.splitlines():
                match = re.match(r"\s*([A-Za-z]\w*)\b", line)
                if match and not line.lstrip().startswith("%"):
                    name = match.group(1)
                    class_members.setdefault(_normalized(name), set()).add(name)
    changed = True
    while changed:
        changed = False
        for child, parent in parents.items():
            if child not in members or parent not in members:
                continue
            for normalized, candidates in members[parent].items():
                destination = members[child].setdefault(normalized, set())
                before = len(destination)
                destination.update(candidates)
                changed = changed or len(destination) != before
    return names, members


def _entry(
    upstream: str,
    kind: str,
    present: bool,
    matlab_name: str | None,
    oracle_status: str = "not-mapped",
) -> dict[str, object]:
    return {
        "kind": kind,
        "matlab_name": matlab_name,
        "oracle_status": oracle_status,
        "surface_status": "present" if present else "missing",
        "upstream": upstream,
    }


def build_coverage() -> dict[str, Any]:
    inventory = json.loads(INVENTORY_PATH.read_text(encoding="utf-8"))
    if inventory.get("optuna_version") != "4.9.0":
        raise RuntimeError("The public API inventory is not pinned to Optuna 4.9.0.")
    oracle = json.loads(ORACLE_PATH.read_text(encoding="utf-8"))
    if oracle.get("optuna_version") != "4.9.0":
        raise RuntimeError("The differential oracle is not pinned to Optuna 4.9.0.")
    required_sections = {
        section
        for sections in CLASS_ORACLE_SECTIONS.values()
        for section in sections
    }
    missing_sections = sorted(required_sections.difference(oracle))
    if missing_sections:
        raise RuntimeError(
            "Class oracle mappings reference missing fixture sections: "
            + ", ".join(missing_sections)
        )
    names, members = _matlab_surface()
    entries: list[dict[str, object]] = []
    for module in inventory["modules"]:
        module_name = str(module["module"])
        for symbol in module["symbols"]:
            name = str(symbol["name"])
            upstream = f"{module_name}.{name}"
            kind = str(symbol["kind"])
            if kind == "module":
                present = MODULE_EQUIVALENTS.get(name, False)
                matlab_name = f"radia.optuna ({name} namespace)" if present else None
            else:
                surface_name = CLASS_EQUIVALENTS.get(
                    name, SYMBOL_EQUIVALENTS.get(name, name)
                )
                present = surface_name in names
                matlab_name = f"radia.optuna.{surface_name}" if present else None
            surface_name = CLASS_EQUIVALENTS.get(name, name)
            class_members_complete = kind == "class" and all(
                surface_name in members
                and _normalized(str(member["name"])) in members[surface_name]
                and str(member["name"]) in VERIFIED_MEMBERS.get(name, set())
                for member in symbol.get("members", [])
            )
            if present and kind == "module" or present and name in VERIFIED_SYMBOLS:
                oracle_status = "verified"
            elif present and kind == "class" and name in CLASS_ORACLE_SECTIONS:
                oracle_status = "verified" if class_members_complete else "partial"
            else:
                oracle_status = "not-mapped"
            entries.append(
                _entry(upstream, kind, present, matlab_name, oracle_status)
            )
            if kind != "class":
                continue
            for member in symbol["members"]:
                member_name = str(member["name"])
                normalized_member = _normalized(member_name)
                member_present = (
                    surface_name in members
                    and normalized_member in members[surface_name]
                )
                if member_present:
                    candidates = members[surface_name][normalized_member]
                    actual_member = (
                        member_name
                        if member_name in candidates
                        else min(candidates, key=lambda value: (value.lower(), value))
                    )
                else:
                    actual_member = None
                member_oracle = (
                    "verified"
                    if member_present
                    and member_name in VERIFIED_MEMBERS.get(name, set())
                    else "not-mapped"
                )
                entries.append(
                    _entry(
                        f"{upstream}.{member_name}",
                        f"class-{member['kind']}",
                        member_present,
                        f"radia.optuna.{surface_name}.{actual_member}"
                        if member_present
                        else None,
                        member_oracle,
                    )
                )
    present_count = sum(entry["surface_status"] == "present" for entry in entries)
    missing_count = len(entries) - present_count
    verified_count = sum(entry["oracle_status"] == "verified" for entry in entries)
    partial_count = sum(entry["oracle_status"] == "partial" for entry in entries)
    return {
        "schema": "radia.optuna49-api-coverage.v1",
        "upstream_version": "4.9.0",
        "upstream_inventory": "tests/matlab/fixtures/optuna49_public_api.json",
        "upstream_inventory_sha256": _sha256(INVENTORY_PATH),
        "upstream_oracle": "tests/matlab/fixtures/optuna49_oracle.json",
        "upstream_oracle_sha256": _sha256(ORACLE_PATH),
        "class_oracle_sections": CLASS_ORACLE_SECTIONS,
        "closure_rule": (
            "full_compatibility_complete is true only when every required public "
            "symbol/member is present and has an upstream differential-oracle mapping; "
            "the two documented MATLAB extensions do not waive shared behavior"
        ),
        "allowed_matlab_extensions": [
            "parallel execution and scheduling",
            "MATLAB table and MAT-file storage",
        ],
        "surface_entry_count": len(entries),
        "surface_present_count": present_count,
        "surface_missing_count": missing_count,
        "oracle_verified_count": verified_count,
        "oracle_partial_count": partial_count,
        "oracle_unmapped_count": len(entries)-verified_count-partial_count,
        "full_compatibility_complete": missing_count == 0 and verified_count == len(entries),
        "entries": entries,
    }


def main() -> None:
    DESTINATION.write_text(
        json.dumps(build_coverage(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(DESTINATION)


if __name__ == "__main__":
    main()
