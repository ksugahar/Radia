from copy import deepcopy

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v52_identity


PROMOTED_CASE_IDS = {
    "v52_public_ml_gradientclip_optimizerstate_schedule_checkpoint_model_owner_mismatch",
    "v52_public_rl_recurrent_hiddenstate_episode_reset_sequencemask_policy_owner_mismatch",
}


def _summary() -> dict[str, object]:
    generation = "matlab-v52-test"
    generations = lambda names: {name: generation for name in names}
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    hidden = ["1" * 64, "2" * 64, "3" * 64, "4" * 64]
    hierarchy = ["plant/velocity", "controller/command", "sensor/current"]
    identity = {
        "ml_training": {
            "generation": generation,
            **generations(("gradient_generation", "optimizer_generation", "schedule_generation", "checkpoint_generation", "owner_generation", "result_generation")),
            "gradient_clip_threshold": 1.0, "result_gradient_clip_threshold": 1.0,
            "gradient_clip_method": "global-l2norm", "result_gradient_clip_method": "global-l2norm",
            "optimizer_name": "adam", "result_optimizer_name": "adam",
            "optimizer_step": 240, "result_optimizer_step": 240,
            "optimizer_state_sha256": "4" * 64, "result_optimizer_state_sha256": "4" * 64,
            "learning_rate_schedule": [0.001, 0.0005, 0.0001], "result_learning_rate_schedule": [0.001, 0.0005, 0.0001],
            "checkpoint_sha256": "5" * 64, "result_checkpoint_sha256": "5" * 64,
            "model_owner": "model:training-v52", "result_model_owner": "model:training-v52", **result,
        },
        "recurrent_rl": {
            "generation": generation,
            **generations(("hidden_generation", "reset_generation", "mask_generation", "policy_generation", "owner_generation", "result_generation")),
            "hidden_state_sha256": hidden, "result_hidden_state_sha256": hidden,
            "episode_reset": [True, False, False, False], "result_episode_reset": [True, False, False, False],
            "sequence_mask": [True, True, True, False], "result_sequence_mask": [True, True, True, False],
            "policy_sha256": "6" * 64, "result_policy_sha256": "6" * 64,
            "policy_owner": "policy:recurrent-v52", "result_policy_owner": "policy:recurrent-v52", **result,
        },
        "codegen": {
            "generation": generation,
            **generations(("target_generation", "numeric_generation", "config_generation", "build_generation", "owner_generation", "result_generation")),
            "target_hardware": "ARM-Cortex-M7", "result_target_hardware": "ARM-Cortex-M7",
            "numeric_type": "single", "result_numeric_type": "single",
            "configuration_sha256": "7" * 64, "result_configuration_sha256": "7" * 64,
            "build_sha256": "8" * 64, "result_build_sha256": "8" * 64,
            "build_owner": "build:embedded-v52", "result_build_owner": "build:embedded-v52", **result,
        },
        "simulink_data_inspector": {
            "generation": generation,
            **generations(("run_generation", "signal_generation", "unit_generation", "interpolation_generation", "owner_generation", "result_generation")),
            "run_id": "run:sdi-v52", "result_run_id": "run:sdi-v52",
            "signal_hierarchy": hierarchy, "result_signal_hierarchy": hierarchy,
            "signal_units": ["m/s", "N", "A"], "result_signal_units": ["m/s", "N", "A"],
            "interpolation": ["linear", "zoh", "zoh"], "result_interpolation": ["linear", "zoh", "zoh"],
            "session_owner": "session:sdi-v52", "result_session_owner": "session:sdi-v52", **result,
        },
    }
    return {"matlab_ml_rl_v52_identity": identity}


def test_v52_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v52_identity(_summary())
    assert result and result["status"] == "ok"


def test_v52_frozen_public_counterfactuals_are_rejected() -> None:
    summary = deepcopy(_summary())
    training = summary["matlab_ml_rl_v52_identity"]["ml_training"]
    training.update({"result_optimizer_name": "sgdm", "result_checkpoint_sha256": "0" * 64, "result_model_owner": "model:stale"})
    recurrent = summary["matlab_ml_rl_v52_identity"]["recurrent_rl"]
    recurrent.update({"result_hidden_state_sha256": list(reversed(recurrent["hidden_state_sha256"])), "result_episode_reset": [False, False, True, False], "result_sequence_mask": [True, False, True, True], "result_policy_owner": "policy:stale"})
    result = validate_matlab_ml_rl_v52_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v52_self_consistent_invalid_training_and_sequence_semantics_are_rejected() -> None:
    summary = deepcopy(_summary())
    training = summary["matlab_ml_rl_v52_identity"]["ml_training"]
    training["learning_rate_schedule"] = training["result_learning_rate_schedule"] = [0.001, 0.002, 0.0001]
    recurrent = summary["matlab_ml_rl_v52_identity"]["recurrent_rl"]
    recurrent["episode_reset"] = recurrent["result_episode_reset"] = [False, False, True, False]
    recurrent["sequence_mask"] = recurrent["result_sequence_mask"] = [True, False, True, True]
    result = validate_matlab_ml_rl_v52_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v52_source_contracts_reject_unsupported_or_unowned_replay() -> None:
    summary = deepcopy(_summary())
    codegen = summary["matlab_ml_rl_v52_identity"]["codegen"]
    codegen["numeric_type"] = codegen["result_numeric_type"] = "half"
    sdi = summary["matlab_ml_rl_v52_identity"]["simulink_data_inspector"]
    sdi["interpolation"] = sdi["result_interpolation"] = ["spline", "zoh", "zoh"]
    result = validate_matlab_ml_rl_v52_identity(summary)
    assert result and result["status"] == "needs_attention"
