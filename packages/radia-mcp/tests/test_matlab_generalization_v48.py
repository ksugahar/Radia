from __future__ import annotations

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v48_identity


PROMOTED_CASE_IDS = {
    "v48_public_autodiff_parameter_order_gradient_tape_objective_fd_spotcheck_owner_mismatch",
    "v48_public_sequence_model_padding_mask_length_shuffle_minibatch_checkpoint_owner_mismatch",
}


def _summary() -> dict[str, object]:
    autodiff = "autodiff-v48"
    sequence = "sequence-v48"
    agentic = "agentic-v48"
    experiment = "experiment-v48"
    return {"matlab_ml_rl_v48_identity": {
        "autodiff": {
            "generation": autodiff, **{key: autodiff for key in ("parameter_generation", "tape_generation", "objective_generation", "gradient_generation", "fd_generation", "checkpoint_generation", "result_generation")},
            "parameter_order": ["w1", "w2"], "result_parameter_order": ["w1", "w2"],
            "gradient_tape_sha256": "1" * 64, "result_gradient_tape_sha256": "1" * 64,
            "objective_id": "objective:validation-loss", "result_objective_id": "objective:validation-loss",
            "gradient": [1.5, -0.25], "result_gradient": [1.5, -0.25],
            "fd_spotcheck_indices": [0, 1], "result_fd_spotcheck_indices": [0, 1],
            "fd_gradient": [1.5001, -0.2499], "result_fd_gradient": [1.5001, -0.2499],
            "checkpoint_owner": "checkpoint:autodiff", "result_checkpoint_owner": "checkpoint:autodiff",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
        "sequence_model": {
            "generation": sequence, **{key: sequence for key in ("padding_generation", "mask_generation", "length_generation", "shuffle_generation", "minibatch_generation", "checkpoint_generation", "result_generation")},
            "padding_policy": "right_zero", "result_padding_policy": "right_zero",
            "padding_mask": [[1, 1, 1], [1, 1, 0]], "result_padding_mask": [[1, 1, 1], [1, 1, 0]],
            "sequence_lengths": [3, 2], "result_sequence_lengths": [3, 2],
            "shuffle_order": [1, 0], "result_shuffle_order": [1, 0],
            "minibatch_row_keys": ["sequence:1", "sequence:0"], "result_minibatch_row_keys": ["sequence:1", "sequence:0"],
            "checkpoint_owner": "checkpoint:sequence", "result_checkpoint_owner": "checkpoint:sequence",
            "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
        },
        "agentic_workspace": {
            "generation": agentic, **{key: agentic for key in ("mutation_generation", "diff_generation", "approval_generation", "rollback_generation", "tool_call_generation", "result_generation")},
            "workspace_mutations": ["models/train.m", "tests/test_train.m"], "result_workspace_mutations": ["models/train.m", "tests/test_train.m"],
            "workspace_diff_sha256": "4" * 64, "result_workspace_diff_sha256": "4" * 64,
            "approval_scope": "workspace_write:approved_paths", "result_approval_scope": "workspace_write:approved_paths",
            "rollback_state": "not_required", "result_rollback_state": "not_required",
            "tool_call_owner": "tool-call:test", "result_tool_call_owner": "tool-call:test",
            "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
        },
        "experiment_selection": {
            "generation": experiment, **{key: experiment for key in ("metric_generation", "direction_generation", "encoding_generation", "trial_generation", "selection_generation", "result_generation")},
            "metric_name": "validation_loss", "result_metric_name": "validation_loss",
            "metric_direction": "minimize", "result_metric_direction": "minimize",
            "categorical_encoding": {"optimizer": ["adam", "sgdm"]}, "result_categorical_encoding": {"optimizer": ["adam", "sgdm"]},
            "trial_row_keys": ["trial=0", "trial=1", "trial=2"], "result_trial_row_keys": ["trial=0", "trial=1", "trial=2"],
            "metric_values": [0.3, 0.2, 0.25], "result_metric_values": [0.3, 0.2, 0.25],
            "best_trial_index": 1, "result_best_trial_index": 1,
            "best_result_row": "trial=1", "result_best_result_row": "trial=1",
            "experiment_owner": "experiment:test", "result_experiment_owner": "experiment:test",
            "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64,
        },
    }}


def test_v48_positive_ml_rl_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v48_identity(_summary())
    assert result and result["status"] == "ok"


def test_v48_autodiff_mutation_is_rejected() -> None:
    summary = _summary()
    row = summary["matlab_ml_rl_v48_identity"]["autodiff"]
    row.update({"result_parameter_order": ["w2", "w1"], "result_gradient_tape_sha256": "a" * 64, "result_objective_id": "objective:training-loss", "result_gradient": [-0.25, 1.5], "result_fd_gradient": [1.2, -0.1], "result_checkpoint_owner": "checkpoint:old"})
    result = validate_matlab_ml_rl_v48_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v48_sequence_mutation_is_rejected() -> None:
    summary = _summary()
    row = summary["matlab_ml_rl_v48_identity"]["sequence_model"]
    row.update({"result_padding_policy": "left_zero", "result_padding_mask": [[1, 1, 1], [0, 1, 1]], "result_shuffle_order": [0, 1], "result_checkpoint_owner": "checkpoint:old"})
    result = validate_matlab_ml_rl_v48_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v48_self_consistent_invalid_conventions_are_rejected() -> None:
    summary = _summary()
    identity = summary["matlab_ml_rl_v48_identity"]
    identity["sequence_model"]["padding_policy"] = identity["sequence_model"]["result_padding_policy"] = "left_zero"
    identity["agentic_workspace"]["approval_scope"] = identity["agentic_workspace"]["result_approval_scope"] = "workspace_write:any"
    result = validate_matlab_ml_rl_v48_identity(summary)
    assert result and result["status"] == "needs_attention"
