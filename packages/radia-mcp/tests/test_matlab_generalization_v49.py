from __future__ import annotations

from copy import deepcopy

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v49_identity


PROMOTED_CASE_IDS = {
    "v49_public_rl_replay_buffer_observation_action_reward_terminal_seed_checkpoint_owner_mismatch",
    "v49_public_ml_normalization_fit_scope_class_encoding_split_fold_metric_owner_mismatch",
}


def _summary() -> dict[str, object]:
    rl = "rl-v49"; ml = "ml-v49"; agentic = "agentic-v49"; parallel = "parallel-v49"
    rows = ["transition:0", "transition:1", "transition:2"]
    observations = [[0.0, 1.0], [0.5, 0.8], [1.0, 0.2]]
    splits = {"train": ["sample:0", "sample:1", "sample:2", "sample:3"], "validation": ["sample:4", "sample:5"], "test": ["sample:6", "sample:7"]}
    workers = [1, 2, 3]; seeds = [7101, 7102, 7103]; streams = ["Threefry:1", "Threefry:2", "Threefry:3"]
    return {"matlab_ml_rl_v49_identity": {
        "rl_replay": {
            "generation": rl, **{key: rl for key in ("buffer_generation", "observation_generation", "action_generation", "reward_generation", "terminal_generation", "seed_generation", "checkpoint_generation", "policy_generation", "result_generation")},
            "replay_row_keys": rows, "result_replay_row_keys": rows, "observations": observations, "result_observations": observations,
            "actions": [0, 1, 0], "result_actions": [0, 1, 0], "rewards": [0.1, 1.0, -0.2], "result_rewards": [0.1, 1.0, -0.2],
            "terminals": [False, False, True], "result_terminals": [False, False, True], "episode_seeds": [4101, 4102, 4103], "result_episode_seeds": [4101, 4102, 4103],
            "checkpoint_owner": "checkpoint:rl", "result_checkpoint_owner": "checkpoint:rl", "policy_owner": "policy:rl", "result_policy_owner": "policy:rl",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        "ml_normalization": {
            "generation": ml, **{key: ml for key in ("normalization_generation", "class_generation", "split_generation", "fold_generation", "metric_generation", "model_generation", "result_generation")},
            "normalization_fit_scope": "training_partition_only", "result_normalization_fit_scope": "training_partition_only",
            "normalization_fit_rows": splits["train"], "result_normalization_fit_rows": splits["train"], "class_encoding": {"cold": 0, "hot": 1}, "result_class_encoding": {"cold": 0, "hot": 1},
            "split_row_keys": splits, "result_split_row_keys": splits, "training_fold_ids": [0, 1, 0, 1], "result_training_fold_ids": [0, 1, 0, 1],
            "metric_row_keys": ["fold:0", "fold:1"], "result_metric_row_keys": ["fold:0", "fold:1"], "metric_values": [0.8, 0.85], "result_metric_values": [0.8, 0.85],
            "model_owner": "model:classifier", "result_model_owner": "model:classifier", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
        "agentic_cancel": {
            "generation": agentic, **{key: agentic for key in ("timeout_generation", "cancel_generation", "partial_output_generation", "cleanup_generation", "tool_call_generation", "result_generation")},
            "timeout_s": 120.0, "result_timeout_s": 120.0, "timed_out": True, "result_timed_out": True,
            "cancel_requested": True, "result_cancel_requested": True, "cancel_completed": True, "result_cancel_completed": True,
            "partial_output_policy": "discard", "result_partial_output_policy": "discard", "session_cleanup": "released", "result_session_cleanup": "released",
            "tool_call_owner": "tool-call:cancel", "result_tool_call_owner": "tool-call:cancel", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
        },
        "parallel_resume": {
            "generation": parallel, **{key: parallel for key in ("worker_generation", "seed_generation", "stream_generation", "trial_generation", "resume_generation", "checkpoint_generation", "experiment_generation", "result_generation")},
            "worker_ids": workers, "result_worker_ids": workers, "worker_seeds": seeds, "result_worker_seeds": seeds, "random_streams": streams, "result_random_streams": streams,
            "trial_state": "completed:12", "result_trial_state": "completed:12", "resume_state": "resumed_from:8", "result_resume_state": "resumed_from:8",
            "checkpoint_owner": "checkpoint:parallel", "result_checkpoint_owner": "checkpoint:parallel", "experiment_owner": "experiment:parallel", "result_experiment_owner": "experiment:parallel",
            "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        },
    }}


def test_v49_positive_ml_rl_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v49_identity(_summary())
    assert result and result["status"] == "ok"


def test_v49_rl_replay_mutation_is_rejected() -> None:
    summary = _summary(); row = summary["matlab_ml_rl_v49_identity"]["rl_replay"]
    row.update({"result_replay_row_keys": ["transition:1", "transition:0", "transition:2"], "result_terminals": [False, True, False], "result_episode_seeds": [4102, 4101, 4103], "result_policy_owner": "policy:other"})
    result = validate_matlab_ml_rl_v49_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v49_ml_normalization_mutation_is_rejected() -> None:
    summary = _summary(); row = summary["matlab_ml_rl_v49_identity"]["ml_normalization"]
    row.update({"result_normalization_fit_scope": "all_samples", "result_class_encoding": {"cold": 1, "hot": 0}, "result_training_fold_ids": [1, 0, 1, 0], "result_model_owner": "model:other"})
    result = validate_matlab_ml_rl_v49_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v49_self_consistent_data_leakage_and_shared_streams_are_rejected() -> None:
    summary = deepcopy(_summary()); identity = summary["matlab_ml_rl_v49_identity"]
    ml = identity["ml_normalization"]; ml["normalization_fit_scope"] = ml["result_normalization_fit_scope"] = "all_samples"
    parallel = identity["parallel_resume"]; parallel["random_streams"] = parallel["result_random_streams"] = ["global", "global", "global"]
    result = validate_matlab_ml_rl_v49_identity(summary)
    assert result and result["status"] == "needs_attention"
