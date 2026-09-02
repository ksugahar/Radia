from __future__ import annotations

import asyncio
import json

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


def test_boolean_integer_fields_are_rejected():
    artifact = _artifact("reinforcement_learning")
    artifact["random_seed"] = True
    artifact["training_episodes"] = True
    artifact["evaluation_episodes"] = True
    result = validate_matlab_ml_rl_artifact(artifact)
    assert result["status"] == "needs_attention"
    assert result["checks"]["seed_is_recorded"] is False
    assert result["checks"]["training_episodes_are_positive"] is False
    assert result["checks"]["evaluation_episodes_are_positive"] is False


def test_matlab_mcp_exposes_ml_guide_and_artifact_gate():
    from radia_mcp.matlab.server import matlab_ml_rl_artifact_gate, mcp

    tool_names = {item.name for item in asyncio.run(mcp.list_tools())}
    assert {"matlab_agentic_ml_guide", "matlab_validation_catalog"} <= tool_names
    catalog = mcp._tool_manager._tools["matlab_validation_catalog"].fn()
    assert "matlab_ml_rl_artifact_gate" in {
        item["name"] for item in catalog["operations"]
    }
    result = json.loads(matlab_ml_rl_artifact_gate(json.dumps(_artifact("regression"))))
    assert result["status"] == "ok"
