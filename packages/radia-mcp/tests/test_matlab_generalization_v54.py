from copy import deepcopy
import math

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v54_identity


PROMOTED_CASE_IDS = {
    "v54_public_rl_prioritizedreplay_priority_beta_seed_policy_buffer_owner_mismatch",
    "v54_public_ml_featurepreprocess_normalization_categoryorder_split_model_owner_mismatch",
}


def _weights(priorities: list[float], alpha: float, beta: float) -> tuple[list[float], list[float]]:
    scaled = [value**alpha for value in priorities]
    probabilities = [value / sum(scaled) for value in scaled]
    raw = [(len(priorities) * probability) ** (-beta) for probability in probabilities]
    return probabilities, [value / max(raw) for value in raw]


def _summary() -> dict[str, object]:
    replay_generation = "replay-v54-test"
    preprocess_generation = "preprocess-v54-test"
    priorities = [1.0, 2.0, 4.0, 8.0]
    probabilities, weights = _weights(priorities, 0.6, 0.4)
    categories = {"material": ["category:steel", "category:copper"], "state": ["category:normal", "category:fault"]}
    split = {"train": [0, 1, 2, 3, 4, 5], "validation": [6, 7], "test": [8, 9]}
    identity = {
        "prioritized_replay": {
            "generation": replay_generation,
            **{name: replay_generation for name in ("priority_generation", "beta_generation", "seed_generation", "policy_generation", "owner_generation", "result_generation")},
            "priorities": priorities, "result_priorities": priorities,
            "priority_alpha": 0.6, "result_priority_alpha": 0.6,
            "sampling_probability": probabilities, "result_sampling_probability": probabilities,
            "importance_beta": 0.4, "result_importance_beta": 0.4,
            "importance_weight": weights, "result_importance_weight": weights,
            "rng_seed": 8675309, "result_rng_seed": 8675309,
            "policy_checkpoint_sha256": "1" * 64, "result_policy_checkpoint_sha256": "1" * 64,
            "buffer_owner": "buffer:v54", "result_buffer_owner": "buffer:v54",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
        "feature_preprocess": {
            "generation": preprocess_generation,
            **{name: preprocess_generation for name in ("normalization_generation", "category_generation", "split_generation", "model_generation", "owner_generation", "result_generation")},
            "feature_mean": [10.0, 2.0], "result_feature_mean": [10.0, 2.0],
            "feature_std": [2.0, 0.5], "result_feature_std": [2.0, 0.5],
            "category_order": categories, "result_category_order": categories,
            "data_split_indices": split, "result_data_split_indices": split,
            "model_checkpoint_sha256": "3" * 64, "result_model_checkpoint_sha256": "3" * 64,
            "preprocessing_owner": "preprocess:v54", "result_preprocessing_owner": "preprocess:v54",
            "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        },
    }
    return {"matlab_ml_rl_v54_identity": identity}


def test_v54_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v54_identity(_summary())
    assert result and result["status"] == "ok"


def test_v54_frozen_public_counterfactuals_are_rejected() -> None:
    summary = deepcopy(_summary())
    summary["matlab_ml_rl_v54_identity"]["prioritized_replay"]["result_rng_seed"] = 42
    summary["matlab_ml_rl_v54_identity"]["feature_preprocess"]["result_model_checkpoint_sha256"] = "9" * 64
    result = validate_matlab_ml_rl_v54_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v54_self_consistent_invalid_semantics_are_rejected() -> None:
    summary = deepcopy(_summary())
    replay = summary["matlab_ml_rl_v54_identity"]["prioritized_replay"]
    replay["importance_weight"] = replay["result_importance_weight"] = [1.0] * 4
    preprocess = summary["matlab_ml_rl_v54_identity"]["feature_preprocess"]
    preprocess["data_split_indices"] = preprocess["result_data_split_indices"] = {"train": [0, 1], "validation": [1], "test": [2]}
    result = validate_matlab_ml_rl_v54_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v54_malformed_values_reject_without_raising() -> None:
    summary = deepcopy(_summary())
    summary["matlab_ml_rl_v54_identity"]["prioritized_replay"]["importance_beta"] = [0.4]
    summary["matlab_ml_rl_v54_identity"]["feature_preprocess"]["feature_std"] = [[2.0], 0.5]
    result = validate_matlab_ml_rl_v54_identity(summary)
    assert result and result["status"] == "needs_attention"
