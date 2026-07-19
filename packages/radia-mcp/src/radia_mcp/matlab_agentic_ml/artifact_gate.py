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


def _number(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False
