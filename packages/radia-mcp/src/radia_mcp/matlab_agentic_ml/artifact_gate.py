"""Small, solver-neutral contract gate for MATLAB ML/RL result artifacts."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping


_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_TASKS = {"regression", "classification", "reinforcement_learning"}


def _integer(value: object, *, minimum: int) -> bool:
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and value >= minimum
    )


def validate_matlab_ml_rl_artifact(value: Mapping[str, object]) -> dict[str, object]:
    """Check provenance, split/evaluation separation, and task-specific metrics.

    This is an artifact gate, not a trainer. It deliberately refuses to infer
    validation from a training score or from an unpinned MATLAB session.
    """
    if not isinstance(value, Mapping):
        raise TypeError("artifact must be a mapping")
    task = value.get("task")
    checks: dict[str, bool] = {
        "schema_is_pinned": value.get("schema") == "cae-ai-lab.matlab-ml-rl-result.v1",
        "task_is_supported": task in _TASKS,
        "matlab_release_is_pinned": bool(str(value.get("matlab_release", "")).strip()),
        "session_owner_is_recorded": str(value.get("session_owner", "")).startswith("matlab:"),
        "seed_is_recorded": _integer(value.get("random_seed"), minimum=0),
        "result_digest_is_valid": bool(_DIGEST.fullmatch(str(value.get("result_sha256", "")))),
        "accepted_digest_matches": value.get("accepted_result_sha256") == value.get("result_sha256"),
        "training_and_evaluation_are_disjoint": value.get("training_evaluation_disjoint") is True,
    }
    if task in {"regression", "classification"}:
        folds = value.get("cross_validation_folds")
        checks["cross_validation_is_explicit"] = _integer(folds, minimum=2)
        checks["holdout_or_cv_score_is_finite"] = math.isfinite(float(value.get("validation_score"))) if _number(value.get("validation_score")) else False
    elif task == "reinforcement_learning":
        checks["environment_is_pinned"] = bool(str(value.get("environment_id", "")).strip())
        checks["training_episodes_are_positive"] = _integer(
            value.get("training_episodes"), minimum=1
        )
        checks["evaluation_episodes_are_positive"] = _integer(
            value.get("evaluation_episodes"), minimum=1
        )
        checks["evaluation_return_is_finite"] = math.isfinite(float(value.get("evaluation_mean_return"))) if _number(value.get("evaluation_mean_return")) else False
    else:
        checks["task_specific_contract_is_present"] = False
    return {
        "policy": "matlab_ml_rl_artifact_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def validate_matlab_ml_rl_v44_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate the v44 split/evaluation identity carried by a solver summary.

    The gate is deliberately independent of a trainer.  It checks the two
    durable promises that an MCP can verify after execution: held-out data are
    disjoint for supervised learning, and RL evaluation uses the same pinned
    environment with a fresh evaluation seed and no exploration.
    """
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v44_identity")
    if not isinstance(identity, Mapping):
        return None

    checks: dict[str, bool] = {}
    supervised = identity.get("supervised")
    rl = identity.get("reinforcement_learning")
    checks["supervised_record_is_mapping"] = isinstance(supervised, Mapping)
    checks["rl_record_is_mapping"] = isinstance(rl, Mapping)
    if isinstance(supervised, Mapping):
        checks.update(_supervised_v44_checks(supervised))
    if isinstance(rl, Mapping):
        checks.update(_rl_v44_checks(rl))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v44",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def validate_matlab_ml_rl_v45_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate v45 official ML/RL and Agentic Toolkit replay bindings."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v45_identity")
    if not isinstance(identity, Mapping):
        return None
    records = {
        "supervised": identity.get("supervised"),
        "reinforcement_learning": identity.get("reinforcement_learning"),
        "agentic_toolkit": identity.get("agentic_toolkit"),
        "mlrl_checkpoint": identity.get("mlrl_checkpoint"),
    }
    checks = {
        "v45_supervised_record_is_mapping": isinstance(records["supervised"], Mapping),
        "v45_rl_record_is_mapping": isinstance(records["reinforcement_learning"], Mapping),
        "v45_agentic_record_is_mapping": isinstance(records["agentic_toolkit"], Mapping),
        "v45_checkpoint_record_is_mapping": isinstance(records["mlrl_checkpoint"], Mapping),
    }
    if isinstance(records["supervised"], Mapping):
        checks.update(_supervised_v45_checks(records["supervised"]))
    if isinstance(records["reinforcement_learning"], Mapping):
        checks.update(_rl_v45_checks(records["reinforcement_learning"]))
    if isinstance(records["agentic_toolkit"], Mapping):
        checks.update(_agentic_v45_checks(records["agentic_toolkit"]))
    if isinstance(records["mlrl_checkpoint"], Mapping):
        checks.update(_checkpoint_v45_checks(records["mlrl_checkpoint"]))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v45",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def validate_matlab_ml_rl_v46_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate official ML/RL replay finiteness, seed, restart, and agent-call bindings."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v46_identity")
    if not isinstance(identity, Mapping):
        return None
    records = {name: identity.get(name) for name in ("supervised", "reinforcement_learning", "agentic_toolkit", "mlrl_checkpoint")}
    checks = {f"v46_{name}_record_is_mapping": isinstance(value, Mapping) for name, value in records.items()}
    if isinstance(records["supervised"], Mapping):
        checks.update(_supervised_v46_checks(records["supervised"]))
    if isinstance(records["reinforcement_learning"], Mapping):
        checks.update(_rl_v46_checks(records["reinforcement_learning"]))
    if isinstance(records["agentic_toolkit"], Mapping):
        checks.update(_agentic_v46_checks(records["agentic_toolkit"]))
    if isinstance(records["mlrl_checkpoint"], Mapping):
        checks.update(_checkpoint_v46_checks(records["mlrl_checkpoint"]))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v46",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def _v46_release_and_digest(value: Mapping[str, object]) -> dict[str, bool]:
    return {
        "v46_release_is_bound": _same_release(value),
        "v46_digest_is_bound": _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256"),
    }


def _supervised_v46_checks(value: Mapping[str, object]) -> dict[str, bool]:
    checks = {
        "v46_supervised_nonfinite_policy_is_bound": value.get("nonfinite_policy") == value.get("result_nonfinite_policy") == "drop_with_count" and value.get("nonfinite_input_count") == value.get("result_nonfinite_input_count") == 0,
        "v46_supervised_split_and_normalization_are_bound": bool(str(value.get("split_generation", "")).strip()) and value.get("result_split_generation") == value.get("split_generation") and value.get("normalization_fit_scope") == value.get("result_normalization_fit_scope") == "training_only",
        "v46_supervised_worker_seed_is_bound": isinstance(value.get("worker_seed"), int) and value.get("worker_seed") >= 0 and value.get("result_worker_seed") == value.get("worker_seed"),
        "v46_supervised_restart_state_is_fresh": value.get("restart_state") == value.get("result_restart_state") == "fresh_training",
    }
    checks.update(_v46_release_and_digest(value))
    return checks


def _rl_v46_checks(value: Mapping[str, object]) -> dict[str, bool]:
    checks = {
        "v46_rl_episode_timeout_is_bound": isinstance(value.get("episode_timeout_steps"), int) and value.get("episode_timeout_steps") > 0 and value.get("result_episode_timeout_steps") == value.get("episode_timeout_steps"),
        "v46_rl_termination_is_environment_defined": value.get("termination_semantics") == value.get("result_termination_semantics") == "environment_defined",
        "v46_rl_exploration_is_training_only": value.get("exploration_mode") == value.get("result_exploration_mode") == "training_only" and value.get("evaluation_mode") == value.get("result_evaluation_mode") == "greedy_no_exploration",
        "v46_rl_checkpoint_is_bound": bool(str(value.get("checkpoint_generation", "")).strip()) and value.get("result_checkpoint_generation") == value.get("checkpoint_generation"),
    }
    checks.update(_v46_release_and_digest(value))
    return checks


def _agentic_v46_checks(value: Mapping[str, object]) -> dict[str, bool]:
    checks = {
        "v46_agentic_argument_shape_is_json_object": value.get("argument_schema") == value.get("result_argument_schema") == "json_object" and value.get("argument_shape_valid") == value.get("result_argument_shape_valid") is True,
        "v46_agentic_existing_session_and_timeout_are_bound": value.get("session_detection") == value.get("result_session_detection") == "existing_shared_matlab" and isinstance(value.get("timeout_s"), (int, float)) and float(value.get("timeout_s")) > 0.0 and value.get("result_timeout_s") == value.get("timeout_s"),
        "v46_agentic_error_class_is_bound": value.get("error_class") == value.get("result_error_class") == "none",
        "v46_agentic_tool_arguments_are_bound": _valid_digest(value.get("tool_arguments_sha256")) and value.get("result_tool_arguments_sha256") == value.get("tool_arguments_sha256") and _same_release(value) and str(value.get("owner", "")).startswith("matlab:") and value.get("result_owner") == value.get("owner"),
    }
    return checks


def _checkpoint_v46_checks(value: Mapping[str, object]) -> dict[str, bool]:
    order = value.get("checkpoint_order")
    checks = {
        "v46_checkpoint_worker_seed_is_bound": isinstance(value.get("worker_seed"), int) and value.get("worker_seed") >= 0 and value.get("result_worker_seed") == value.get("worker_seed"),
        "v46_checkpoint_order_is_monotone": isinstance(order, list) and order == value.get("result_checkpoint_order") and order == sorted(order),
        "v46_checkpoint_optimizer_state_is_bound": value.get("optimizer_state") == value.get("result_optimizer_state") == "adam_ready",
        "v46_checkpoint_generation_and_owner_are_bound": bool(str(value.get("checkpoint_generation", "")).strip()) and value.get("result_checkpoint_generation") == value.get("checkpoint_generation") and _same_release(value) and str(value.get("owner", "")).startswith("matlab:") and value.get("result_owner") == value.get("owner"),
    }
    checks.update({"v46_checkpoint_digest_is_bound": _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256")})
    return checks


def _same_release(value: Mapping[str, object]) -> bool:
    release = str(value.get("release_id", "")).strip()
    return bool(release) and value.get("result_release_id") == release


def _supervised_v45_checks(value: Mapping[str, object]) -> dict[str, bool]:
    return {
        "v45_supervised_schema_is_pinned": value.get("schema") == "cae-ai-lab.matlab-ml-rl-result.v3",
        "v45_datastore_is_pinned": bool(str(value.get("datastore_id", "")).strip()),
        "v45_split_generation_is_bound": bool(str(value.get("split_generation", "")).strip()) and value.get("result_split_generation") == value.get("split_generation"),
        "v45_normalization_is_fit_on_training_only": value.get("normalization_fit_scope") == "training_only" and value.get("result_normalization_fit_scope") == "training_only",
        "v45_hyperparameter_selection_uses_training_cv": value.get("hyperparameter_selection_source") == "training_cross_validation" and value.get("result_hyperparameter_selection_source") == "training_cross_validation",
        "v45_holdout_is_separate": bool(str(value.get("holdout_id", "")).strip()) and value.get("holdout_id") not in set(_string_list(value.get("training_ids"))),
        "v45_model_card_matches_release": value.get("model_card_release") == value.get("release_id"),
        "v45_release_is_bound": _same_release(value),
        "v45_digest_is_bound": _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256"),
    }


def _rl_v45_checks(value: Mapping[str, object]) -> dict[str, bool]:
    discount = value.get("discount_factor")
    return {
        "v45_rl_replay_buffer_is_bound": bool(str(value.get("replay_buffer_generation", "")).strip()) and value.get("result_replay_buffer_generation") == value.get("replay_buffer_generation"),
        "v45_rl_termination_semantics_are_environment_defined": value.get("termination_semantics") == "environment_defined" and value.get("result_termination_semantics") == "environment_defined",
        "v45_rl_discount_is_valid_and_bound": _number(discount) and 0.0 < float(discount) <= 1.0 and value.get("result_discount_factor") == discount,
        "v45_rl_policy_seed_is_recorded": isinstance(value.get("policy_seed"), int) and value.get("policy_seed") >= 0,
        "v45_rl_evaluation_is_separate": value.get("evaluation_mode") == "greedy_no_exploration" and value.get("result_evaluation_mode") == "greedy_no_exploration",
        "v45_rl_release_is_bound": _same_release(value),
        "v45_rl_digest_is_bound": _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256"),
    }


def _agentic_v45_checks(value: Mapping[str, object]) -> dict[str, bool]:
    return {
        "v45_agentic_capability_route_is_pinned": value.get("capability_route") == "matlab_tool" and value.get("result_capability_route") == "matlab_tool",
        "v45_agentic_consent_is_recorded": value.get("consent_recorded") is True and value.get("result_consent_recorded") is True,
        "v45_agentic_tool_arguments_digest_is_bound": _valid_digest(value.get("tool_arguments_sha256")) and value.get("result_tool_arguments_sha256") == value.get("tool_arguments_sha256"),
        "v45_agentic_existing_session_is_explicit": value.get("session_detection") == "existing_shared_matlab" and value.get("result_session_detection") == "existing_shared_matlab",
        "v45_agentic_release_and_owner_are_bound": _same_release(value) and str(value.get("owner", "")).startswith("matlab:") and value.get("result_owner") == value.get("owner"),
    }


def _checkpoint_v45_checks(value: Mapping[str, object]) -> dict[str, bool]:
    return {
        "v45_checkpoint_generation_is_bound": bool(str(value.get("checkpoint_generation", "")).strip()) and value.get("result_checkpoint_generation") == value.get("checkpoint_generation"),
        "v45_datastore_state_is_replayed": value.get("datastore_state") == "replayed" and value.get("result_datastore_state") == "replayed",
        "v45_optimizer_is_recorded": bool(str(value.get("optimizer", "")).strip()) and value.get("result_optimizer") == value.get("optimizer"),
        "v45_discount_is_bound": _number(value.get("discount_factor")) and value.get("result_discount_factor") == value.get("discount_factor"),
        "v45_evaluation_id_is_separate": bool(str(value.get("evaluation_id", "")).strip()) and value.get("evaluation_id") != value.get("checkpoint_generation"),
        "v45_checkpoint_release_and_owner_are_bound": _same_release(value) and str(value.get("owner", "")).startswith("matlab:") and value.get("result_owner") == value.get("owner"),
        "v45_checkpoint_digest_is_bound": _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256"),
    }


def _supervised_v44_checks(value: Mapping[str, object]) -> dict[str, bool]:
    training_ids = _string_list(value.get("training_ids"))
    evaluation_ids = _string_list(value.get("evaluation_ids"))
    validation_metric = value.get("validation_metric")
    metric_name = validation_metric.get("name") if isinstance(validation_metric, Mapping) else None
    metric_value = validation_metric.get("value") if isinstance(validation_metric, Mapping) else None
    return {
        "supervised_schema_is_pinned": value.get("schema") == "cae-ai-lab.matlab-ml-rl-result.v2",
        "supervised_task_is_regression_or_classification": value.get("task") in {"regression", "classification"},
        "supervised_release_is_pinned": bool(str(value.get("matlab_release", "")).strip()),
        "supervised_session_owner_is_recorded": str(value.get("session_owner", "")).startswith("matlab:"),
        "supervised_seed_is_recorded": isinstance(value.get("random_seed"), int) and value.get("random_seed") >= 0,
        "supervised_split_is_named": bool(str(value.get("split_id", "")).strip()),
        "supervised_training_ids_are_nonempty": bool(training_ids),
        "supervised_evaluation_ids_are_nonempty": bool(evaluation_ids),
        "supervised_training_and_evaluation_ids_are_disjoint": bool(training_ids) and bool(evaluation_ids) and set(training_ids).isdisjoint(evaluation_ids),
        "supervised_training_and_evaluation_are_disjoint": value.get("training_evaluation_disjoint") is True,
        "supervised_validation_metric_is_finite": _number(metric_value),
        "supervised_validation_metric_is_not_training_score": bool(str(metric_name).strip()) and not str(metric_name).lower().startswith("training"),
        "supervised_result_digest_is_valid": _valid_digest(value.get("result_sha256")),
        "supervised_accepted_digest_matches": value.get("accepted_result_sha256") == value.get("result_sha256"),
        "supervised_timing_breakdown_is_bounded": _valid_timings(value.get("timing_breakdown_s")),
    }


def _rl_v44_checks(value: Mapping[str, object]) -> dict[str, bool]:
    evaluation_return = value.get("evaluation_mean_return")
    evaluation_std = value.get("evaluation_std_return")
    return {
        "rl_schema_is_pinned": value.get("schema") == "cae-ai-lab.matlab-ml-rl-result.v2",
        "rl_task_is_supported": value.get("task") == "reinforcement_learning",
        "rl_release_is_pinned": bool(str(value.get("matlab_release", "")).strip()),
        "rl_session_owner_is_recorded": str(value.get("session_owner", "")).startswith("matlab:"),
        "rl_seed_is_recorded": isinstance(value.get("random_seed"), int) and value.get("random_seed") >= 0,
        "rl_environment_is_pinned": bool(str(value.get("environment_id", "")).strip()),
        "rl_evaluation_environment_matches": value.get("evaluation_environment_id") == value.get("environment_id") and bool(str(value.get("evaluation_environment_id", "")).strip()),
        "rl_training_episodes_are_positive": isinstance(value.get("training_episodes"), int) and value.get("training_episodes") > 0,
        "rl_evaluation_episodes_are_positive": isinstance(value.get("evaluation_episodes"), int) and value.get("evaluation_episodes") > 0,
        "rl_evaluation_seed_is_fresh": isinstance(value.get("evaluation_seed"), int) and value.get("evaluation_seed") >= 0 and value.get("evaluation_seed") != value.get("random_seed"),
        "rl_evaluation_has_no_exploration": value.get("exploration_during_evaluation") is False,
        "rl_training_and_evaluation_are_disjoint": value.get("training_evaluation_disjoint") is True,
        "rl_evaluation_return_is_finite": _number(evaluation_return),
        "rl_evaluation_std_is_finite_nonnegative": _number(evaluation_std) and float(evaluation_std) >= 0.0,
        "rl_result_digest_is_valid": _valid_digest(value.get("result_sha256")),
        "rl_accepted_digest_matches": value.get("accepted_result_sha256") == value.get("result_sha256"),
        "rl_evaluation_digest_is_valid": _valid_digest(value.get("evaluation_result_sha256")),
        "rl_accepted_evaluation_digest_matches": value.get("accepted_evaluation_result_sha256") == value.get("evaluation_result_sha256"),
        "rl_timing_breakdown_is_bounded": _valid_timings(value.get("timing_breakdown_s")),
    }


def _string_list(value: object) -> list[str]:
    if not isinstance(value, (list, tuple)):
        return []
    return [item for item in value if isinstance(item, str) and item.strip()]


def _valid_digest(value: object) -> bool:
    return bool(_DIGEST.fullmatch(str(value)))


def _valid_timings(value: object) -> bool:
    if not isinstance(value, Mapping) or not 1 <= len(value) <= 4:
        return False
    return all(_number(item) and float(item) >= 0.0 for item in value.values())


def _number(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False
