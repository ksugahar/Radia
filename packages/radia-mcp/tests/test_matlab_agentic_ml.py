from __future__ import annotations

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_artifact


def _artifact(task: str) -> dict[str, object]:
    value: dict[str, object] = {
        "schema": "cae-ai-lab.matlab-ml-rl-result.v1",
        "task": task,
        "matlab_release": "R2026a",
        "session_owner": "matlab:shared-session",
        "random_seed": 843,
        "result_sha256": "a" * 64,
        "accepted_result_sha256": "a" * 64,
        "training_evaluation_disjoint": True,
    }
    if task in {"regression", "classification"}:
        value.update({"cross_validation_folds": 5, "validation_score": 0.1})
    else:
        value.update({"environment_id": "cae:design-control-v1", "training_episodes": 100, "evaluation_episodes": 20, "evaluation_mean_return": 1.5})
    return value


def test_supervised_artifact_requires_explicit_cross_validation():
    assert validate_matlab_ml_rl_artifact(_artifact("regression"))["status"] == "ok"


def test_rl_artifact_requires_disjoint_evaluation_and_environment():
    assert validate_matlab_ml_rl_artifact(_artifact("reinforcement_learning"))["status"] == "ok"


def test_training_only_result_is_rejected():
    artifact = _artifact("reinforcement_learning")
    artifact["training_evaluation_disjoint"] = False
    assert validate_matlab_ml_rl_artifact(artifact)["status"] == "needs_attention"
