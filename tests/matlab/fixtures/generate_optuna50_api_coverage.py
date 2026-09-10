"""Compare the Optuna 5.0.0 public inventory with the MATLAB package surface."""

from __future__ import annotations

import ast
import hashlib
import json
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
INVENTORY_PATH = Path(__file__).with_name("optuna50_public_api.json")
ORACLE_PATH = Path(__file__).with_name("optuna50_oracle.json")
MATLAB_DIRECTORY = ROOT / "matlab" / "+radia" / "+optuna"
DESTINATION = ROOT / "matlab" / "optuna50_api_coverage.json"
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
        "constraints",
        "datetime_start",
        "distributions",
        "number",
        "params",
        "report",
        "set_constraint",
        "set_user_attr",
        "should_prune",
        "suggest_categorical",
        "suggest_discrete_uniform",
        "suggest_float",
        "suggest_int",
        "suggest_loguniform",
        "suggest_uniform",
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
        "set_user_attr",
        "stop",
        "tell",
        "trials",
        "trials_dataframe",
        "user_attrs",
    },
    "StudySummary": {"direction", "directions"},
    "Trial": {
        "constraints",
        "datetime_start",
        "distributions",
        "number",
        "params",
        "report",
        "relative_params",
        "set_constraint",
        "set_user_attr",
        "should_prune",
        "suggest_categorical",
        "suggest_discrete_uniform",
        "suggest_float",
        "suggest_int",
        "suggest_loguniform",
        "suggest_uniform",
        "user_attrs",
    },
    "FixedTrial": {
        "constraints",
        "datetime_start",
        "distributions",
        "number",
        "params",
        "report",
        "set_constraint",
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
        "constraints",
        "datetime_complete",
        "datetime_start",
        "distributions",
        "duration",
        "intermediate_values",
        "last_step",
        "number",
        "params",
        "report",
        "set_constraint",
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
    "BaseMutation": {"mutation"},
    "PolynomialMutation": {"mutation"},
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
        "tpe_group",
        "tpe_sampler_seed_37",
        "tpe_constant_liar_seed_127",
    ),
    "BaseCrossover": ("base_components", "nsgaii_crossovers_seed_73"),
    "BaseMutation": ("nsgaii_mutation",),
    "BLXAlphaCrossover": ("nsgaii_crossovers_seed_73",),
    "SBXCrossover": ("nsgaii_crossovers_seed_73",),
    "SPXCrossover": ("nsgaii_crossovers_seed_73",),
    "UNDXCrossover": ("nsgaii_crossovers_seed_73",),
    "UniformCrossover": ("nsgaii_crossovers_seed_73",),
    "VSBXCrossover": ("nsgaii_crossovers_seed_73",),
    "PolynomialMutation": ("nsgaii_mutation",),
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
        if "+internal" in path.parts or "private" in path.parts:
            continue
        class_folders = [part[1:] for part in path.parts if part.startswith("@")]
        owner_name = class_folders[-1] if class_folders else path.stem
        names.add(owner_name)
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
        class_members = members.setdefault(owner_name, {})
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


# Why a module is not held to the "port it" bar. Anything absent from this
# table is required: it must be present and oracle-mapped before compatibility
# closure can be true. A downgrade must name what discharges it.
SCOPE_RULES = {
    "optuna.storages": {
        "scope": "bridged",
        "reason": (
            "MATLAB keeps history in normalized tables and MAT files. "
            "Interoperability is discharged by the explicit storage bridge."
        ),
        "discharged_by": [
            "radia.optuna.export_study",
            "radia.optuna.import_study",
            "radia_optuna.bridge",
        ],
    },
    "optuna.visualization": {
        "scope": "replaced",
        "reason": (
            "MATLAB plotting functions consume the same normalized study "
            "history; upstream Python visualizations remain available after "
            "an explicit storage handoff."
        ),
        "discharged_by": [
            "radia.optuna.plot_*",
            "radia_optuna.bridge",
        ],
    },
    "optuna.integration": {
        "scope": "replaced",
        "reason": (
            "MATLAB-native integrations and explicit Python ecosystem "
            "adapters replace framework-specific Python callback objects."
        ),
        "discharged_by": [
            "radia.optuna integration classes",
            "radia_optuna.bridge",
        ],
    },
    "optuna.artifacts": {
        "scope": "replaced",
        "reason": (
            "MATLAB artifact stores and Radia result contracts replace "
            "Python artifact-store objects."
        ),
        "discharged_by": ["radia.optuna artifact-store classes"],
    },
    "optuna.logging": {
        "scope": "replaced",
        "reason": (
            "MATLAB reports through warning/error identifiers and display "
            "settings rather than a Python logging hierarchy."
        ),
        "discharged_by": ["MATLAB warning and error identifiers"],
    },
    "optuna.exceptions": {
        "scope": "replaced",
        "reason": (
            "MATLAB signals failures through named exception classes and "
            "radia:optuna:* identifiers."
        ),
        "discharged_by": [
            "radia.optuna exception classes",
            "radia:optuna:* error identifiers",
        ],
    },
}


# The inventory walks __mro__, so an IntEnum drags in int methods and an
# exception class drags in BaseException methods. Those are Python-language
# surface, not Optuna API.
PYTHON_LANGUAGE_MEMBERS = frozenset(
    name
    for base in (int, BaseException)
    for name in dir(base)
    if not name.startswith("_")
)


def _scope_for(upstream: str) -> str:
    leaf = upstream.rsplit(".", 1)[-1]
    if leaf in PYTHON_LANGUAGE_MEMBERS:
        return "python-language"
    if leaf.startswith("_"):
        return "out-of-scope"
    for prefix, rule in SCOPE_RULES.items():
        if upstream == prefix or upstream.startswith(prefix + "."):
            return str(rule["scope"])
    return "required"


# --- constructor-default audit ------------------------------------------
#
# The inventory records classes and their public members but not __init__, so
# a changed constructor default was invisible here while changing every seeded
# result.  Optuna 5.0 is that case: TPESampler moved multivariate False -> None
# and constant_liar False -> True; only the second had a named default test.
#
# The oracle records every upstream default literally.  This pairs them with
# the MATLAB `arguments` block and requires each parameter to either match, or
# be declared below with a reason.  A declaration that no longer applies is an
# error too, so the tables cannot rot.

MATLAB_ARGUMENT = re.compile(
    r"^\s*options\.(?P<name>\w+)(?P<declaration>[^=\n]*?)=\s*(?P<default>.+?)\s*$"
)

# Upstream parameter -> MATLAB argument, where stripping underscores and
# lowercasing does not already pair them.
CONSTRUCTOR_PARAMETER_ALIASES: dict[str, str] = {
    "constraints_func": "ConstraintsFcn",
    "crossover_prob": "CrossoverProbability",
    "mutation_prob": "MutationProbability",
    "swapping_prob": "SwappingProbability",
    "n_ei_candidates": "NumberOfEIChoices",
    "n_min_trials": "MinCompletedTrials",
    "popsize": "PopulationSize",
    "weights": "WeightsFcn",
    "gamma": "GammaFcn",
}

UNSET = (
    "MATLAB types this argument as a scalar double, which cannot hold [], so "
    "NaN is the unset sentinel for upstream None."
)

# (class, upstream parameter) -> why the MATLAB default differs but agrees.
# Every entry states what upstream's literal default resolves to; where that
# needed measuring rather than reading, the measured value is quoted.
CONSTRUCTOR_DEFAULT_EQUIVALENCES: dict[tuple[str, str], str] = {
    ("CmaEsSampler", "popsize"): (
        "MATLAB PopulationSize=0 is the unset sentinel; the constructor "
        "rejects anything between zero and two, so zero cannot be a real "
        "population."
    ),
    ("CmaEsSampler", "sigma0"): UNSET,
    ("CmaEsSampler", "x0"): (
        "MATLAB X0 is typed struct, so an empty struct is the unset value for "
        "upstream None."
    ),
    ("EMMREvaluator", "seed"): UNSET,
    ("FloatDistribution", "step"): UNSET,
    ("FrozenTrial", "values"): UNSET,
    ("MaxTrialsCallback", "states"): (
        "Upstream carries a TrialState tuple; MATLAB carries the equivalent "
        "state name as a string, matching how every other MATLAB entry point "
        "spells a trial state."
    ),
    ("NSGAIIISampler", "mutation_prob"): UNSET,
    ("NSGAIIISampler", "reference_points"): (
        "MATLAB ReferencePoints is a numeric matrix, so a 0x0 matrix is the "
        "unset value for upstream None."
    ),
    ("NSGAIISampler", "mutation_prob"): UNSET,
    ("RDBStorage", "engine_kwargs"): (
        "MATLAB engine_kwargs is typed struct, so an empty struct is the "
        "unset value for upstream None."
    ),
    ("RDBStorage", "grace_period"): UNSET,
    ("RDBStorage", "heartbeat_interval"): UNSET,
    ("RegretBoundEvaluator", "seed"): UNSET,
    ("RetryFailedTrialCallback", "max_retry"): UNSET,
    ("RetryHeartbeatStaleTrialCallback", "max_retry"): UNSET,
    ("SBXCrossover", "eta"): UNSET,
    ("SPXCrossover", "epsilon"): UNSET,
    ("TPESampler", "consider_endpoints"): (
        "Optuna 5.0 types this bool | None, where None resolves internally. "
        "Measured on the pinned build: TPESampler() resolves it to False, "
        "which is the MATLAB default."
    ),
    ("TPESampler", "consider_magic_clip"): (
        "Optuna 5.0 types this bool | None. Measured on the pinned build: "
        "TPESampler() resolves it to True, which is the MATLAB default."
    ),
    ("TPESampler", "prior_weight"): (
        "Optuna 5.0 types this float | None. Measured on the pinned build: "
        "TPESampler() resolves it to 1.0, which is the MATLAB default."
    ),
    ("TPESampler", "warn_independent_sampling"): (
        "Optuna 5.0 types this bool | None. Measured on the pinned build: "
        "TPESampler() resolves it to False, which is the MATLAB default."
    ),
    ("ThresholdPruner", "lower"): UNSET,
    ("ThresholdPruner", "upper"): UNSET,
    ("UNDXCrossover", "sigma_eta"): UNSET,
    ("VSBXCrossover", "eta"): UNSET,
}

# (class, upstream parameter) -> a default that genuinely disagrees with
# upstream.  These are recorded limitations, NOT equivalences: the audit
# counts them separately and the tests pin the set, so a new one cannot
# appear unnoticed.  Do not move an entry here to make a run pass.
CONSTRUCTOR_DEFAULT_DIVERGENCES: dict[tuple[str, str], str] = {
    ("Terminator", "improvement_evaluator"): (
        "Upstream Terminator() resolves None to RegretBoundEvaluator "
        "(measured on the pinned build); MATLAB defaults to "
        "BestValueStagnationEvaluator, so an unconfigured Terminator stops on "
        "a different criterion than upstream."
    ),
    ("FanovaImportanceEvaluator", "seed"): (
        "Upstream seed=None means fresh entropy per instance; MATLAB defaults "
        "to a fixed 0, so two unseeded MATLAB evaluators agree with each "
        "other where two upstream ones do not."
    ),
    ("MeanDecreaseImpurityImportanceEvaluator", "seed"): (
        "Same divergence as FanovaImportanceEvaluator.seed: upstream None is "
        "fresh entropy, MATLAB pins 0."
    ),
}

# (class, upstream parameter) -> why MATLAB carries no counterpart.
CONSTRUCTOR_PARAMETERS_NOT_IMPLEMENTED: dict[tuple[str, str], str] = {
    ("Boto3ArtifactStore", "client"): (
        "Artifact stores are scoped out-of-scope; a boto3 client object has "
        "no MATLAB counterpart."
    ),
    ("GPSampler", "independent_sampler"): (
        "MATLAB GPSampler does not expose an independent-sampler override; "
        "recorded as a limitation rather than a silently different knob."
    ),
    ("GPSampler", "warn_independent_sampling"): (
        "Follows independent_sampler: with no override there is no "
        "independent-sampling warning to configure."
    ),
    ("NSGAIIISampler", "elite_population_selection_strategy"): (
        "MATLAB NSGA-III uses the built-in reference-point elite selection "
        "and exposes no strategy override."
    ),
    ("SBXCrossover", "uniform_crossover_prob"): (
        "MATLAB implements the SBX crossover proper; the per-gene uniform "
        "mixing knobs are not exposed."
    ),
    ("SBXCrossover", "use_child_gene_prob"): (
        "Follows uniform_crossover_prob."
    ),
    ("VSBXCrossover", "uniform_crossover_prob"): (
        "Follows SBXCrossover.uniform_crossover_prob."
    ),
    ("VSBXCrossover", "use_child_gene_prob"): (
        "Follows SBXCrossover.use_child_gene_prob."
    ),
}

_MATLAB_ABSENT = {
    "[]",
    "{}",
    "double.empty(1,0)",
    "double.empty",
    "string.empty(1,0)",
    "string.empty",
    "cell(1,0)",
}


def _canonical_upstream_default(literal: str) -> tuple:
    text = literal.strip()
    if text == "None":
        return ("absent",)
    if text in ("True", "False"):
        return ("bool", text == "True")
    try:
        return ("number", float(text))
    except ValueError:
        pass
    if len(text) >= 2 and text[0] in "'\"" and text[-1] == text[0]:
        return ("string", text[1:-1])
    return ("opaque", text)


def _canonical_matlab_default(literal: str) -> tuple:
    text = literal.strip().rstrip(";").strip()
    if text in _MATLAB_ABSENT:
        return ("absent",)
    if text in ("true", "false"):
        return ("bool", text == "true")
    try:
        return ("number", float(text))
    except ValueError:
        pass
    if len(text) >= 2 and text[0] in "'\"" and text[-1] == text[0]:
        return ("string", text[1:-1])
    return ("opaque", text)


def _matlab_constructor_defaults() -> dict[str, dict[str, str]]:
    """Defaulted `arguments` entries of every public MATLAB constructor."""
    constructors: dict[str, dict[str, str]] = {}
    for path in sorted(MATLAB_DIRECTORY.rglob("*.m")):
        # Judge folders relative to the package root: +internal and MATLAB
        # private/ helpers are not public, and a checkout path that merely
        # contains such a folder name must not hide the whole surface.
        relative = path.relative_to(MATLAB_DIRECTORY)
        if "+internal" in relative.parts or "private" in relative.parts:
            continue
        class_folders = [
            part[1:] for part in relative.parts[:-1] if part.startswith("@")
        ]
        if class_folders and path.stem != class_folders[-1]:
            continue  # External method file; the constructor is in the class file.
        stem = path.stem
        source = path.read_text(encoding="utf-8")
        match = re.search(
            rf"function\s+\w+\s*=\s*{re.escape(stem)}\s*\([^)]*\)"
            rf"\s*\n\s*arguments\s*\n(?P<block>.*?)\n\s*end",
            source,
            re.S,
        )
        if not match:
            continue
        # MATLAB continues a line with "..."; join first, or every continued
        # argument is silently dropped and the audit under-reports.
        block = re.sub(r"\.\.\.[^\n]*\n\s*", " ", match.group("block"))
        defaults: dict[str, str] = {}
        for line in block.splitlines():
            line = line.split("%", 1)[0]
            hit = MATLAB_ARGUMENT.match(line)
            if hit:
                defaults[hit.group("name")] = hit.group("default").strip()
        if defaults:
            constructors[stem] = defaults
    return constructors


def _audit_constructor_defaults(oracle: dict[str, Any]) -> dict[str, object]:
    record = oracle.get("constructor_defaults")
    if not isinstance(record, dict) or "modules" not in record:
        raise RuntimeError(
            "The oracle carries no constructor_defaults section; regenerate it "
            "with the pinned Optuna before building the coverage ledger."
        )
    matlab = _matlab_constructor_defaults()
    matched: list[str] = []
    equivalent: list[str] = []
    not_implemented: list[str] = []
    divergent: list[str] = []
    undeclared: list[str] = []
    used_aliases: set[str] = set()
    used_equivalences: set[tuple[str, str]] = set()
    used_absences: set[tuple[str, str]] = set()
    used_divergences: set[tuple[str, str]] = set()

    for module_name, classes in sorted(record["modules"].items()):
        for class_name, parameters in sorted(classes.items()):
            if class_name not in matlab:
                # Not implemented in MATLAB at all; the surface ledger above
                # already accounts for that, so there is no default to audit.
                continue
            arguments = matlab[class_name]
            normalized = {
                name.replace("_", "").lower(): name for name in arguments
            }
            for parameter, detail in sorted(parameters.items()):
                key = (class_name, parameter)
                alias = CONSTRUCTOR_PARAMETER_ALIASES.get(parameter)
                if alias and alias in arguments:
                    used_aliases.add(parameter)
                    argument = alias
                else:
                    argument = normalized.get(parameter.replace("_", "").lower())
                label = f"{module_name}.{class_name}.{parameter}"
                if argument is None:
                    if key in CONSTRUCTOR_PARAMETERS_NOT_IMPLEMENTED:
                        used_absences.add(key)
                        not_implemented.append(label)
                    else:
                        undeclared.append(
                            f"{label}: no MATLAB argument, and no declared reason"
                        )
                    continue
                upstream_value = _canonical_upstream_default(detail["default_repr"])
                matlab_value = _canonical_matlab_default(arguments[argument])
                if upstream_value == matlab_value:
                    matched.append(label)
                elif key in CONSTRUCTOR_DEFAULT_EQUIVALENCES:
                    used_equivalences.add(key)
                    equivalent.append(label)
                elif key in CONSTRUCTOR_DEFAULT_DIVERGENCES:
                    used_divergences.add(key)
                    divergent.append(label)
                else:
                    undeclared.append(
                        f"{label}: upstream {detail['default_repr']} vs MATLAB "
                        f"{argument} = {arguments[argument]}"
                    )

    if undeclared:
        raise RuntimeError(
            "Constructor defaults disagree with the pinned oracle and carry no "
            "declared reason:\n  " + "\n  ".join(sorted(undeclared))
        )
    stale_aliases = sorted(set(CONSTRUCTOR_PARAMETER_ALIASES) - used_aliases)
    stale_equivalences = sorted(
        f"{cls}.{param}"
        for cls, param in set(CONSTRUCTOR_DEFAULT_EQUIVALENCES) - used_equivalences
    )
    stale_absences = sorted(
        f"{cls}.{param}"
        for cls, param in set(CONSTRUCTOR_PARAMETERS_NOT_IMPLEMENTED) - used_absences
    )
    stale_divergences = sorted(
        f"{cls}.{param}"
        for cls, param in set(CONSTRUCTOR_DEFAULT_DIVERGENCES) - used_divergences
    )
    stale = (
        stale_aliases + stale_equivalences + stale_absences + stale_divergences
    )
    if stale:
        raise RuntimeError(
            "Constructor-default declarations no longer apply; remove them: "
            + ", ".join(stale)
        )

    return {
        "audited_class_count": len(
            {
                label.rsplit(".", 1)[0]
                for label in matched + equivalent + not_implemented + divergent
            }
        ),
        "matched_count": len(matched),
        "declared_equivalent_count": len(equivalent),
        "declared_not_implemented_count": len(not_implemented),
        "declared_divergent_count": len(divergent),
        "declared_divergences": {
            f"{cls}.{param}": reason
            for (cls, param), reason in sorted(
                CONSTRUCTOR_DEFAULT_DIVERGENCES.items()
            )
        },
        "declared_equivalences": {
            f"{cls}.{param}": reason
            for (cls, param), reason in sorted(
                CONSTRUCTOR_DEFAULT_EQUIVALENCES.items()
            )
        },
        "declared_not_implemented": {
            f"{cls}.{param}": reason
            for (cls, param), reason in sorted(
                CONSTRUCTOR_PARAMETERS_NOT_IMPLEMENTED.items()
            )
        },
        "parameter_aliases": dict(sorted(CONSTRUCTOR_PARAMETER_ALIASES.items())),
    }


def _qualified_optuna_references(tree: ast.AST) -> set[str]:
    """Conservative identity evidence: never infer an owner from a suffix.

    Instance aliases/dynamic dispatch need explicit owner-aware contracts.
    Until those exist they remain asserted, not verified. Strings, local
    variable names and unrelated receivers are not Optuna API identities.
    """
    def qualified(node: ast.AST) -> str | None:
        if isinstance(node, ast.Name) and node.id == "optuna":
            return "optuna"
        if isinstance(node, ast.Attribute):
            owner = qualified(node.value)
            return f"{owner}.{node.attr}" if owner else None
        return None

    return {
        name
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute)
        if (name := qualified(node)) is not None
    }


def _oracle_generator_sections() -> dict[str, set[str]]:
    """Map qualified Optuna references to their fixture producer sections.

    Presence in VERIFIED_SYMBOLS or VERIFIED_MEMBERS is an assertion, not
    evidence. The mapping below is derived from the pinned upstream fixture
    generators, so an entry without a generating section cannot be counted as
    verified.
    """
    fixtures = Path(__file__).resolve().parent
    oracle_source = fixtures / "generate_optuna50_oracle.py"
    mcp_source = fixtures / "generate_optuna50_mcp_oracle.py"
    tree = ast.parse(oracle_source.read_text(encoding="utf-8"))
    functions = {
        node.name: node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef)
    }
    if "build_oracle" not in functions:
        raise RuntimeError(
            f"{oracle_source.name} no longer defines build_oracle(); "
            "oracle-section derivation cannot be trusted."
        )

    section_producers: dict[str, str] = {}
    returned = next(
        (
            node.value
            for node in functions["build_oracle"].body
            if isinstance(node, ast.Return) and isinstance(node.value, ast.Dict)
        ),
        None,
    )
    if returned is None:
        raise RuntimeError("build_oracle() must directly return its fixture dictionary.")
    for key, value in zip(returned.keys, returned.values):
        if (
            isinstance(key, ast.Constant)
            and isinstance(key.value, str)
            and isinstance(value, ast.Call)
            and isinstance(value.func, ast.Name)
        ):
            section_producers[key.value] = value.func.id

    fixture_sections = set(
        json.loads(
            (fixtures / "optuna50_oracle.json").read_text(encoding="utf-8")
        )
    )
    unknown = sorted(set(section_producers) - fixture_sections)
    if unknown:
        raise RuntimeError(
            "build_oracle() names sections absent from optuna50_oracle.json: "
            + ", ".join(unknown)
        )

    def referenced(
        function_name: str, depth: int = 2, seen: set[str] | None = None
    ) -> set[str]:
        seen = seen or set()
        found: set[str] = set()
        if function_name in seen or function_name not in functions:
            return found
        seen.add(function_name)
        found.update(_qualified_optuna_references(functions[function_name]))
        for node in ast.walk(functions[function_name]):
            if (
                depth > 0
                and isinstance(node, ast.Call)
                and isinstance(node.func, ast.Name)
            ):
                found |= referenced(node.func.id, depth - 1, seen)
        return found

    sections: dict[str, set[str]] = {}
    for section, producer in section_producers.items():
        for name in referenced(producer):
            sections.setdefault(name, set()).add(section)

    mcp_tree = ast.parse(mcp_source.read_text(encoding="utf-8"))
    for name in _qualified_optuna_references(mcp_tree):
        sections.setdefault(name, set()).add("mcp:" + mcp_source.stem)
    return sections


ORACLE_SECTIONS_BY_NAME = _oracle_generator_sections()


def _oracle_sections_for(upstream: str) -> list[str]:
    return sorted(ORACLE_SECTIONS_BY_NAME.get(upstream, ()))


def _entry(
    upstream: str,
    kind: str,
    present: bool,
    matlab_name: str | None,
    oracle_status: str = "not-mapped",
) -> dict[str, object]:
    sections = _oracle_sections_for(upstream)
    if oracle_status in ("verified", "partial") and not sections:
        oracle_status = "asserted"
    return {
        "kind": kind,
        "matlab_name": matlab_name,
        "oracle_sections": sections,
        "oracle_status": oracle_status,
        "scope": _scope_for(upstream),
        "surface_status": "present" if present else "missing",
        "upstream": upstream,
    }


def _matlab_qualified_names() -> dict[str, str]:
    """Derive resolvable public names from MATLAB package directories."""
    result: dict[str, str] = {}
    for path in sorted(MATLAB_DIRECTORY.rglob("*.m")):
        relative = path.relative_to(MATLAB_DIRECTORY)
        if "+internal" in relative.parts or "private" in relative.parts:
            continue
        folders = relative.parts[:-1]
        class_folders = [part[1:] for part in folders if part.startswith("@")]
        if class_folders and path.stem != class_folders[-1]:
            continue  # External method: identity belongs to its owning class.
        packages = [part for part in folders if not part.startswith("@")]
        if any(not part.startswith("+") for part in packages):
            raise RuntimeError(f"Non-package public MATLAB path: {relative}")
        qualified = ".".join(
            ["radia", "optuna", *(part[1:] for part in packages), path.stem]
        )
        if path.stem in result:
            raise RuntimeError(f"Ambiguous MATLAB public basename: {path.stem}")
        result[path.stem] = qualified
    return result


def build_coverage() -> dict[str, Any]:
    inventory = json.loads(INVENTORY_PATH.read_text(encoding="utf-8"))
    if inventory.get("optuna_version") != "5.0.0":
        raise RuntimeError("The public API inventory is not pinned to Optuna 5.0.0.")
    oracle = json.loads(ORACLE_PATH.read_text(encoding="utf-8"))
    if oracle.get("optuna_version") != "5.0.0":
        raise RuntimeError("The differential oracle is not pinned to Optuna 5.0.0.")
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
    qualified_names = _matlab_qualified_names()
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
                if present and (MATLAB_DIRECTORY / f"+{name}").is_dir():
                    matlab_name = f"radia.optuna.{name}"
            else:
                surface_name = CLASS_EQUIVALENTS.get(
                    name, SYMBOL_EQUIVALENTS.get(name, name)
                )
                present = surface_name in names
                matlab_name = qualified_names[surface_name] if present else None
            surface_name = CLASS_EQUIVALENTS.get(name, name)
            class_members_complete = kind == "class" and all(
                surface_name in members
                and _normalized(str(member["name"])) in members[surface_name]
                and str(member["name"]) in VERIFIED_MEMBERS.get(name, set())
                for member in symbol.get("members", [])
            )
            if present and (kind == "module" or name in VERIFIED_SYMBOLS):
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
                        f"{qualified_names[surface_name]}.{actual_member}"
                        if member_present
                        else None,
                        member_oracle,
                    )
                )
    present_count = sum(entry["surface_status"] == "present" for entry in entries)
    missing_count = len(entries) - present_count
    verified_count = sum(entry["oracle_status"] == "verified" for entry in entries)
    partial_count = sum(entry["oracle_status"] == "partial" for entry in entries)
    asserted_count = sum(
        entry["oracle_status"] == "asserted" for entry in entries
    )
    required = [entry for entry in entries if entry["scope"] == "required"]
    required_present = [
        entry for entry in required if entry["surface_status"] == "present"
    ]
    required_mapped = [
        entry
        for entry in required_present
        if entry["oracle_status"] not in ("not-mapped", "asserted")
    ]
    required_asserted = [
        entry
        for entry in required_present
        if entry["oracle_status"] == "asserted"
    ]
    scope_counts: dict[str, int] = {}
    for entry in entries:
        scope = str(entry["scope"])
        scope_counts[scope] = scope_counts.get(scope, 0) + 1
    complete = (
        len(required_present) == len(required)
        and len(required_mapped) == len(required)
        and not required_asserted
    )
    return {
        "schema": "radia.optuna50-api-coverage.v1",
        "upstream_version": "5.0.0",
        "upstream_inventory": "tests/matlab/fixtures/optuna50_public_api.json",
        "upstream_inventory_sha256": _sha256(INVENTORY_PATH),
        "upstream_oracle": "tests/matlab/fixtures/optuna50_oracle.json",
        "upstream_oracle_sha256": _sha256(ORACLE_PATH),
        "class_oracle_sections": CLASS_ORACLE_SECTIONS,
        "closure_rule": (
            "full_compatibility_complete is true only when every entry whose scope "
            "is 'required' is present and has evidence derived from an upstream "
            "differential-oracle section; assertion-only mappings do not pass"
        ),
        "scope_rules": SCOPE_RULES,
        "constructor_default_audit": _audit_constructor_defaults(oracle),
        "oracle_asserted_count": asserted_count,
        "required_oracle_asserted_count": len(required_asserted),
        "scope_counts": scope_counts,
        "required_entry_count": len(required),
        "required_present_count": len(required_present),
        "required_missing_count": len(required) - len(required_present),
        "required_oracle_mapped_count": len(required_mapped),
        "required_oracle_unmapped_count": len(required) - len(required_mapped),
        "allowed_matlab_extensions": [
            "parallel execution and scheduling",
            "MATLAB table and MAT-file storage",
        ],
        "surface_entry_count": len(entries),
        "surface_present_count": present_count,
        "surface_missing_count": missing_count,
        "oracle_verified_count": verified_count,
        "oracle_partial_count": partial_count,
        "oracle_unmapped_count": len(entries)
        - verified_count
        - partial_count
        - asserted_count,
        "full_compatibility_complete": complete,
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
