from copy import deepcopy

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v51_identity


PROMOTED_CASE_IDS = {
    "v51_public_rl_replay_transition_terminal_truncation_discount_priority_target_owner_mismatch",
    "v51_public_classification_classorder_labelencoding_split_confusionmatrix_model_owner_mismatch",
}


def _summary() -> dict[str, object]:
    generation = "matlab-public-v51"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    transitions = ["transition:1", "transition:2", "transition:3", "transition:4"]
    classes = ["normal", "warning", "fault"]
    identity = {
        "rl_replay": {
            "generation": generation, **{name: generation for name in ("transition_generation", "terminal_generation", "discount_generation", "priority_generation", "target_generation", "owner_generation", "result_generation")},
            "transition_ids": transitions, "result_transition_ids": transitions, "terminal_flags": [False, False, False, True],
            "result_terminal_flags": [False, False, False, True], "truncation_flags": [False] * 4,
            "result_truncation_flags": [False] * 4, "discounts": [0.99, 0.99, 0.99, 0.0],
            "result_discounts": [0.99, 0.99, 0.99, 0.0], "priorities": [0.2, 0.5, 0.4, 1.0],
            "result_priorities": [0.2, 0.5, 0.4, 1.0], "target_network_sha256": "1" * 64,
            "result_target_network_sha256": "1" * 64, "buffer_owner": "buffer:prioritized-v51",
            "result_buffer_owner": "buffer:prioritized-v51", **result,
        },
        "classification": {
            "generation": generation, **{name: generation for name in ("class_generation", "label_generation", "split_generation", "confusion_generation", "owner_generation", "result_generation")},
            "class_order": classes, "result_class_order": classes, "label_encoding": {"normal": 0, "warning": 1, "fault": 2},
            "result_label_encoding": {"normal": 0, "warning": 1, "fault": 2},
            "data_split_sha256": {"train": "2" * 64, "validation": "3" * 64, "test": "4" * 64},
            "result_data_split_sha256": {"train": "2" * 64, "validation": "3" * 64, "test": "4" * 64},
            "confusion_matrix": [[18, 1, 0], [2, 15, 1], [0, 1, 12]],
            "result_confusion_matrix": [[18, 1, 0], [2, 15, 1], [0, 1, 12]],
            "confusion_matrix_axes": {"rows": classes, "columns": classes},
            "result_confusion_matrix_axes": {"rows": classes, "columns": classes},
            "model_owner": "model:classifier-v51", "result_model_owner": "model:classifier-v51", **result,
        },
        "agentic_edit": {
            "generation": generation, **{name: generation for name in ("precondition_generation", "patch_generation", "rollback_generation", "tool_call_generation", "owner_generation", "result_generation")},
            "target_path": "workspace:project/src/agent.m", "result_target_path": "workspace:project/src/agent.m",
            "file_precondition_sha256": "5" * 64, "result_file_precondition_sha256": "5" * 64,
            "patch_sha256": "6" * 64, "result_patch_sha256": "6" * 64, "rollback_state": "available",
            "result_rollback_state": "available", "rollback_sha256": "5" * 64, "result_rollback_sha256": "5" * 64,
            "tool_call_owner": "tool-call:edit-v51", "result_tool_call_owner": "tool-call:edit-v51", **result,
        },
        "parallel_rng": {
            "generation": generation, **{name: generation for name in ("rng_generation", "substream_generation", "worker_generation", "reduction_generation", "checkpoint_generation", "owner_generation", "result_generation")},
            "rng_algorithm": "Threefry", "result_rng_algorithm": "Threefry", "rng_seed": 510901, "result_rng_seed": 510901,
            "substream_ids": [1, 2, 3, 4], "result_substream_ids": [1, 2, 3, 4], "worker_map": [1, 2, 1, 2],
            "result_worker_map": [1, 2, 1, 2], "reduction_order": [1, 3, 2, 4], "result_reduction_order": [1, 3, 2, 4],
            "checkpoint_sha256": "7" * 64, "result_checkpoint_sha256": "7" * 64,
            "pool_owner": "pool:threads-v51", "result_pool_owner": "pool:threads-v51", **result,
        },
    }
    return {"matlab_ml_rl_v51_identity": identity}


def test_v51_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v51_identity(_summary())
    assert result and result["status"] == "ok"


def test_v51_frozen_public_counterfactuals_are_rejected() -> None:
    summary = deepcopy(_summary())
    summary["matlab_ml_rl_v51_identity"]["rl_replay"].update({"result_transition_ids": ["transition:4"], "result_buffer_owner": "buffer:stale"})
    summary["matlab_ml_rl_v51_identity"]["classification"].update({"result_class_order": ["fault", "warning", "normal"], "result_model_owner": "model:stale"})
    result = validate_matlab_ml_rl_v51_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v51_self_consistent_wrong_semantics_are_rejected() -> None:
    summary = deepcopy(_summary())
    rl = summary["matlab_ml_rl_v51_identity"]["rl_replay"]
    rl["discounts"] = rl["result_discounts"] = [0.99] * 4
    classification = summary["matlab_ml_rl_v51_identity"]["classification"]
    classification["class_order"] = classification["result_class_order"] = ["fault", "warning", "normal"]
    result = validate_matlab_ml_rl_v51_identity(summary)
    assert result and result["status"] == "needs_attention"
