from __future__ import annotations

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v45_identity


_PROMOTED_CASE_IDS = (
    "v45_public_supervised_datastore_normalization_hyperparameter_holdout_modelcard_release_mismatch",
    "v45_public_reinforcementlearning_replay_termination_discount_policy_seed_evaluation_release_mismatch",
    "v45_source_agentictoolkit_capability_routing_consent_tool_arguments_session_release_owner_mismatch",
    "v45_source_mlrl_checkpoint_datastore_seed_optimizer_discount_evaluation_release_owner_mismatch",
)


def _identity() -> dict[str, object]:
    return {
        "matlab_ml_rl_v45_identity": {
            "supervised": {"schema": "cae-ai-lab.matlab-ml-rl-result.v3", "datastore_id": "ds:845", "split_generation": "split:845", "result_split_generation": "split:845", "normalization_fit_scope": "training_only", "result_normalization_fit_scope": "training_only", "hyperparameter_selection_source": "training_cross_validation", "result_hyperparameter_selection_source": "training_cross_validation", "holdout_id": "holdout:845", "training_ids": ["train:1"], "model_card_release": "R2026a", "release_id": "R2026a", "result_release_id": "R2026a", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64},
            "reinforcement_learning": {"replay_buffer_generation": "replay:845", "result_replay_buffer_generation": "replay:845", "termination_semantics": "environment_defined", "result_termination_semantics": "environment_defined", "discount_factor": 0.99, "result_discount_factor": 0.99, "policy_seed": 845, "evaluation_mode": "greedy_no_exploration", "result_evaluation_mode": "greedy_no_exploration", "release_id": "R2026a", "result_release_id": "R2026a", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64},
            "agentic_toolkit": {"capability_route": "matlab_tool", "result_capability_route": "matlab_tool", "consent_recorded": True, "result_consent_recorded": True, "tool_arguments_sha256": "c" * 64, "result_tool_arguments_sha256": "c" * 64, "session_detection": "existing_shared_matlab", "result_session_detection": "existing_shared_matlab", "release_id": "R2026a", "result_release_id": "R2026a", "owner": "matlab:shared", "result_owner": "matlab:shared"},
            "mlrl_checkpoint": {"checkpoint_generation": "checkpoint:845", "result_checkpoint_generation": "checkpoint:845", "datastore_state": "replayed", "result_datastore_state": "replayed", "optimizer": "adam", "result_optimizer": "adam", "discount_factor": 0.99, "result_discount_factor": 0.99, "evaluation_id": "eval:845", "release_id": "R2026a", "result_release_id": "R2026a", "owner": "matlab:shared", "result_owner": "matlab:shared", "result_sha256": "d" * 64, "accepted_result_sha256": "d" * 64},
        }
    }


def test_v45_matlab_ml_rl_identity_positive() -> None:
    result = validate_matlab_ml_rl_v45_identity(_identity())
    assert result is not None
    assert result["status"] == "ok"


def test_v45_matlab_identity_rejects_holdout_and_agentic_consent_mutations() -> None:
    identity = _identity()
    identity["matlab_ml_rl_v45_identity"]["supervised"]["holdout_id"] = "train:1"
    identity["matlab_ml_rl_v45_identity"]["agentic_toolkit"]["result_consent_recorded"] = False
    result = validate_matlab_ml_rl_v45_identity(identity)
    assert result is not None
    assert result["status"] == "needs_attention"
