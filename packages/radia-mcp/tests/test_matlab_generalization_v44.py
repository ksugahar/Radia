from __future__ import annotations

import json
from copy import deepcopy

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v44_identity
from radia_mcp.radia_ngsolve.server import regularized_trace_inverse_path_gate
from test_matlab_generalization_v42 import _summary_v42


_CASE_IDS = (
    "v44_public_supervised_split_validation_metric_result_identity_mismatch",
    "v44_public_rl_environment_evaluation_episode_seed_exploration_result_identity_mismatch",
    "v44_source_official_mcp_ml_result_release_session_split_owner_mismatch",
    "v44_source_official_mcp_rl_evaluation_environment_seed_episode_result_mismatch",
)


def _identity() -> dict[str, object]:
    return {
        "supervised": {
            "schema": "cae-ai-lab.matlab-ml-rl-result.v2",
            "task": "regression",
            "matlab_release": "R2026a",
            "session_owner": "matlab:shared-session-v44",
            "random_seed": 844,
            "split_id": "holdout:ml-rl-v44",
            "training_ids": ["train-01", "train-02"],
            "evaluation_ids": ["eval-01"],
            "training_evaluation_disjoint": True,
            "validation_metric": {"name": "rmse", "value": 0.1, "units": "normalized"},
            "result_sha256": "5" * 64,
            "accepted_result_sha256": "5" * 64,
            "timing_breakdown_s": {"train": 0.2, "evaluate": 0.1},
        },
        "reinforcement_learning": {
            "schema": "cae-ai-lab.matlab-ml-rl-result.v2",
            "task": "reinforcement_learning",
            "matlab_release": "R2026a",
            "session_owner": "matlab:shared-session-v44",
            "random_seed": 845,
            "environment_id": "cae:rl-design-control-v44",
            "evaluation_environment_id": "cae:rl-design-control-v44",
            "training_episodes": 100,
            "evaluation_episodes": 20,
            "evaluation_seed": 1845,
            "exploration_during_evaluation": False,
            "training_evaluation_disjoint": True,
            "evaluation_mean_return": 1.5,
            "evaluation_std_return": 0.1,
            "result_sha256": "6" * 64,
            "accepted_result_sha256": "6" * 64,
            "evaluation_result_sha256": "7" * 64,
            "accepted_evaluation_result_sha256": "7" * 64,
            "timing_breakdown_s": {"train": 0.3, "evaluate": 0.2},
        },
    }


def _summary() -> dict[str, object]:
    value = deepcopy(_summary_v42())
    value["matlab_ml_rl_v44_identity"] = _identity()
    return value


def test_v44_public_ml_rl_identity_positive_is_accepted() -> None:
    result = validate_matlab_ml_rl_v44_identity(_summary())
    assert result is not None
    assert result["status"] == "ok"
    assert len(_CASE_IDS) == 4


def test_v44_public_supervised_rejects_split_and_validation_contamination() -> None:
    value = _summary()
    identity = value["matlab_ml_rl_v44_identity"]["supervised"]
    identity.update({
        "evaluation_ids": ["train-02"],
        "training_evaluation_disjoint": False,
        "validation_metric": {"name": "training_rmse", "value": 0.0},
        "accepted_result_sha256": "8" * 64,
    })
    result = regularized_trace_inverse_path_gate(json.dumps(value))
    assert json.loads(result)["status"] == "needs_attention"


def test_v44_public_rl_rejects_environment_seed_episode_and_exploration_mismatch() -> None:
    value = _summary()
    identity = value["matlab_ml_rl_v44_identity"]["reinforcement_learning"]
    identity.update({
        "evaluation_environment_id": "cae:rl-design-control-old",
        "evaluation_seed": 845,
        "evaluation_episodes": 0,
        "exploration_during_evaluation": True,
        "accepted_evaluation_result_sha256": "9" * 64,
    })
    result = validate_matlab_ml_rl_v44_identity(value)
    assert result is not None
    assert result["status"] == "needs_attention"
