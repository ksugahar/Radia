from __future__ import annotations

from radia_mcp.matlab_agentic_ml.artifact_gate import validate_matlab_ml_rl_v46_identity


_PROMOTED_CASE_IDS = (
    "v46_public_supervised_ml_nan_inf_split_normalization_worker_seed_restart_mismatch",
    "v46_public_reinforcementlearning_episode_timeout_termination_exploration_checkpoint_mismatch",
    "v46_source_tool_agentic_toolkit_argument_shape_error_session_attach_timeout_mismatch",
    "v46_source_tool_mlrl_parallel_worker_seed_checkpoint_order_optimizer_state_mismatch",
)


def _summary():
    return {"matlab_ml_rl_v46_identity": {
        "supervised": {"nonfinite_policy": "drop_with_count", "result_nonfinite_policy": "drop_with_count", "nonfinite_input_count": 0, "result_nonfinite_input_count": 0, "split_generation": "split:test", "result_split_generation": "split:test", "normalization_fit_scope": "training_only", "result_normalization_fit_scope": "training_only", "worker_seed": 846, "result_worker_seed": 846, "restart_state": "fresh_training", "result_restart_state": "fresh_training", "release_id": "R2026a", "result_release_id": "R2026a", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64},
        "reinforcement_learning": {"episode_timeout_steps": 200, "result_episode_timeout_steps": 200, "termination_semantics": "environment_defined", "result_termination_semantics": "environment_defined", "exploration_mode": "training_only", "result_exploration_mode": "training_only", "checkpoint_generation": "checkpoint:test", "result_checkpoint_generation": "checkpoint:test", "evaluation_mode": "greedy_no_exploration", "result_evaluation_mode": "greedy_no_exploration", "release_id": "R2026a", "result_release_id": "R2026a", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64},
        "agentic_toolkit": {"argument_schema": "json_object", "result_argument_schema": "json_object", "argument_shape_valid": True, "result_argument_shape_valid": True, "session_detection": "existing_shared_matlab", "result_session_detection": "existing_shared_matlab", "timeout_s": 90.0, "result_timeout_s": 90.0, "error_class": "none", "result_error_class": "none", "tool_arguments_sha256": "c" * 64, "result_tool_arguments_sha256": "c" * 64, "release_id": "R2026a", "result_release_id": "R2026a", "owner": "matlab:shared", "result_owner": "matlab:shared"},
        "mlrl_checkpoint": {"worker_seed": 846, "result_worker_seed": 846, "checkpoint_order": [0, 1, 2], "result_checkpoint_order": [0, 1, 2], "optimizer_state": "adam_ready", "result_optimizer_state": "adam_ready", "checkpoint_generation": "checkpoint:test", "result_checkpoint_generation": "checkpoint:test", "release_id": "R2026a", "result_release_id": "R2026a", "owner": "matlab:shared", "result_owner": "matlab:shared", "result_sha256": "d" * 64, "accepted_result_sha256": "d" * 64},
    }}


def test_v46_matlab_ml_rl_identity_accepts_closed_artifact():
    result = validate_matlab_ml_rl_v46_identity(_summary())
    assert result and result["status"] == "ok" and all(result["checks"].values())


def test_v46_matlab_ml_rl_identity_rejects_nan_seed_episode_and_argument_mutations():
    summary = _summary()
    identity = summary["matlab_ml_rl_v46_identity"]
    identity["supervised"]["result_nonfinite_policy"] = "keep_nan"
    identity["reinforcement_learning"]["result_episode_timeout_steps"] = 1
    identity["agentic_toolkit"]["result_argument_shape_valid"] = False
    identity["mlrl_checkpoint"]["result_checkpoint_order"] = [0, 2, 1]
    result = validate_matlab_ml_rl_v46_identity(summary)
    assert result and result["status"] == "needs_attention"
