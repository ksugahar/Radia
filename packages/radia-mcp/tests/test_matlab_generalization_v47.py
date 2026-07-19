from __future__ import annotations

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v47_identity


PROMOTED_CASE_IDS = {
    "v47_public_ml_datastore_label_order_preprocess_hash_partition_owner_mismatch",
    "v47_public_rl_observation_action_spec_environment_reset_policy_owner_mismatch",
}


def _summary() -> dict[str, object]:
    ml_generation = "ml-v47"
    rl_generation = "rl-v47"
    agentic_generation = "agentic-v47"
    trial_generation = "trial-v47"
    partition = {"train": ["sample-1", "sample-2"], "validation": ["sample-3"]}
    return {
        "matlab_ml_rl_v47_identity": {
            "ml_datastore": {
                "generation": ml_generation,
                **{key: ml_generation for key in ("datastore_generation", "label_generation", "preprocess_generation", "partition_generation", "owner_generation", "result_generation")},
                "datastore_order": ["sample-1", "sample-2", "sample-3"], "result_datastore_order": ["sample-1", "sample-2", "sample-3"],
                "labels": ["cat", "dog", "cat"], "result_labels": ["cat", "dog", "cat"],
                "preprocess_sha256": "1" * 64, "result_preprocess_sha256": "1" * 64,
                "partition": partition, "result_partition": partition,
                "partition_owner": "partition:test", "result_partition_owner": "partition:test",
                "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
            },
            "reinforcement_learning": {
                "generation": rl_generation,
                **{key: rl_generation for key in ("observation_generation", "action_generation", "environment_generation", "reset_generation", "policy_generation", "result_generation")},
                "observation_spec_sha256": "3" * 64, "result_observation_spec_sha256": "3" * 64,
                "action_spec_sha256": "4" * 64, "result_action_spec_sha256": "4" * 64,
                "environment_id": "environment:test", "result_environment_id": "environment:test",
                "reset_policy": "deterministic_seeded", "result_reset_policy": "deterministic_seeded",
                "policy_owner": "policy:test", "result_policy_owner": "policy:test",
                "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
            },
            "agentic_toolkit": {
                "generation": agentic_generation,
                **{key: agentic_generation for key in ("call_generation", "correlation_generation", "session_generation", "workspace_generation", "release_generation", "result_generation")},
                "tool_call_id": "call:test", "result_tool_call_id": "call:test",
                "correlation_id": "correlation:test", "result_correlation_id": "correlation:test",
                "session_identity": "matlab:shared", "result_session_identity": "matlab:shared",
                "workspace_sha256": "6" * 64, "result_workspace_sha256": "6" * 64,
                "release_id": "R2026a", "result_release_id": "R2026a",
                "result_sha256": "7" * 64, "accepted_result_sha256": "7" * 64,
            },
            "experiment_trials": {
                "generation": trial_generation,
                **{key: trial_generation for key in ("experiment_generation", "trial_generation", "cv_generation", "result_index_generation", "result_generation")},
                "trial_row_keys": ["trial=0", "trial=1", "trial=2"], "result_trial_row_keys": ["trial=0", "trial=1", "trial=2"],
                "cv_partition_ids": ["fold=0", "fold=1", "fold=2"], "result_cv_partition_ids": ["fold=0", "fold=1", "fold=2"],
                "result_indices": [0, 1, 2], "replayed_result_indices": [0, 1, 2],
                "experiment_owner": "experiment:test", "result_experiment_owner": "experiment:test",
                "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
            },
        }
    }


def test_v47_positive_ml_rl_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v47_identity(_summary())
    assert result and result["status"] == "ok"


def test_v47_ml_datastore_partition_mutation_is_rejected() -> None:
    summary = _summary()
    row = summary["matlab_ml_rl_v47_identity"]["ml_datastore"]
    row["result_datastore_order"] = ["sample-2", "sample-1", "sample-3"]
    row["result_labels"] = ["dog", "cat", "cat"]
    row["result_preprocess_sha256"] = "a" * 64
    row["result_partition_owner"] = "partition:other"
    result = validate_matlab_ml_rl_v47_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v47_rl_spec_environment_policy_mutation_is_rejected() -> None:
    summary = _summary()
    row = summary["matlab_ml_rl_v47_identity"]["reinforcement_learning"]
    row["result_observation_spec_sha256"] = "b" * 64
    row["result_action_spec_sha256"] = "c" * 64
    row["result_environment_id"] = "environment:other"
    row["result_reset_policy"] = "random_unseeded"
    row["result_policy_owner"] = "policy:other"
    result = validate_matlab_ml_rl_v47_identity(summary)
    assert result and result["status"] == "needs_attention"
