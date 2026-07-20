from copy import deepcopy
import math

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v53_identity


PROMOTED_CASE_IDS = {
    "v53_public_offline_rl_behaviorpolicy_importanceweight_support_dataset_owner_mismatch",
    "v53_public_ml_probability_calibration_classorder_temperature_model_owner_mismatch",
}


def _softmax(row: list[float], temperature: float) -> list[float]:
    scaled = [item / temperature for item in row]
    maximum = max(scaled)
    values = [math.exp(item - maximum) for item in scaled]
    total = sum(values)
    return [item / total for item in values]


def _summary() -> dict[str, object]:
    offline_generation = "offline-rl-v53-test"
    calibration_generation = "calibration-v53-test"
    behavior = [0.40, 0.25, 0.20, 0.15]
    target = [0.35, 0.30, 0.20, 0.15]
    weights = [target_item / behavior_item for target_item, behavior_item in zip(target, behavior)]
    logits = [[2.0, 0.5, -1.0], [0.1, 1.2, -0.2]]
    temperature = 1.4
    probabilities = [_softmax(row, temperature) for row in logits]
    identity = {
        "offline_rl": {
            "generation": offline_generation,
            **{name: offline_generation for name in ("behavior_generation", "weight_generation", "support_generation", "dataset_generation", "owner_generation", "result_generation")},
            "behavior_policy_sha256": "1" * 64, "result_behavior_policy_sha256": "1" * 64,
            "target_policy_sha256": "2" * 64, "result_target_policy_sha256": "2" * 64,
            "behavior_action_probability": behavior, "result_behavior_action_probability": behavior,
            "target_action_probability": target, "result_target_action_probability": target,
            "importance_weight": weights, "result_importance_weight": weights,
            "support_coverage_fraction": 1.0, "result_support_coverage_fraction": 1.0,
            "dataset_owner": "dataset:offline-v53", "result_dataset_owner": "dataset:offline-v53",
            "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
        },
        "probability_calibration": {
            "generation": calibration_generation,
            **{name: calibration_generation for name in ("logit_generation", "probability_generation", "class_generation", "temperature_generation", "owner_generation", "result_generation")},
            "logits": logits, "result_logits": logits,
            "calibrated_probability": probabilities, "result_calibrated_probability": probabilities,
            "class_order": ["class:normal", "class:warning", "class:fault"],
            "result_class_order": ["class:normal", "class:warning", "class:fault"],
            "temperature": temperature, "result_temperature": temperature,
            "model_owner": "model:calibration-v53", "result_model_owner": "model:calibration-v53",
            "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        },
    }
    return {"matlab_ml_rl_v53_identity": identity}


def test_v53_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v53_identity(_summary())
    assert result and result["status"] == "ok"


def test_v53_frozen_public_counterfactuals_are_rejected() -> None:
    summary = deepcopy(_summary())
    offline = summary["matlab_ml_rl_v53_identity"]["offline_rl"]
    offline.update({"result_behavior_policy_sha256": "8" * 64, "result_importance_weight": [10.0] * 4, "result_support_coverage_fraction": 0.25, "result_dataset_owner": "dataset:stale"})
    calibration = summary["matlab_ml_rl_v53_identity"]["probability_calibration"]
    calibration.update({"result_calibrated_probability": [[0.9, 0.9, 0.9]], "result_class_order": ["class:fault", "class:normal"], "result_temperature": 0.1, "result_model_owner": "model:stale"})
    result = validate_matlab_ml_rl_v53_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v53_self_consistent_invalid_offline_semantics_are_rejected() -> None:
    summary = deepcopy(_summary())
    offline = summary["matlab_ml_rl_v53_identity"]["offline_rl"]
    offline["importance_weight"] = offline["result_importance_weight"] = [1.0] * 4
    offline["support_coverage_fraction"] = offline["result_support_coverage_fraction"] = 0.5
    result = validate_matlab_ml_rl_v53_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v53_self_consistent_invalid_calibration_semantics_are_rejected() -> None:
    summary = deepcopy(_summary())
    calibration = summary["matlab_ml_rl_v53_identity"]["probability_calibration"]
    calibration["calibrated_probability"] = calibration["result_calibrated_probability"] = [[0.8, 0.1, 0.1], [0.8, 0.1, 0.1]]
    calibration["class_order"] = calibration["result_class_order"] = ["class:normal", "class:normal", "class:fault"]
    calibration["temperature"] = calibration["result_temperature"] = -1.0
    result = validate_matlab_ml_rl_v53_identity(summary)
    assert result and result["status"] == "needs_attention"
