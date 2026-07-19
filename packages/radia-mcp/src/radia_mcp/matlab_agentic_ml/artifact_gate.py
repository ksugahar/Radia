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
