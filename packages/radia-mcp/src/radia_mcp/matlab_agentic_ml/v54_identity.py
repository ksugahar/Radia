"""Replay-identity gates for prioritized replay and ML preprocessing."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence


_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_TOL = 1.0e-10


def validate_matlab_ml_rl_v54_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate v54 prioritized-replay and feature-preprocessing identities."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v54_identity")
    if not isinstance(identity, Mapping):
        return None
    replay = identity.get("prioritized_replay")
    preprocess = identity.get("feature_preprocess")
    checks = {
        "v54_prioritized_replay_record_is_mapping": isinstance(replay, Mapping),
        "v54_feature_preprocess_record_is_mapping": isinstance(preprocess, Mapping),
    }
    if isinstance(replay, Mapping):
        checks.update(_prioritized_checks(replay))
    if isinstance(preprocess, Mapping):
        checks.update(_preprocess_checks(preprocess))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v54",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def _prioritized_checks(value: Mapping[str, object]) -> dict[str, bool]:
    priorities = _float_vector(value.get("priorities"))
    probabilities = _float_vector(value.get("sampling_probability"))
    weights = _float_vector(value.get("importance_weight"))
    alpha = value.get("priority_alpha")
    beta = value.get("importance_beta")
    parameters_ok = (
        bool(priorities)
        and all(item > 0.0 for item in priorities)
        and _finite(alpha) and 0.0 < float(alpha) <= 1.0
        and _finite(beta) and 0.0 <= float(beta) <= 1.0
    )
    if parameters_ok:
        scaled = [item ** float(alpha) for item in priorities]
        total = sum(scaled)
        expected_probabilities = [item / total for item in scaled]
        raw_weights = [(len(priorities) * item) ** (-float(beta)) for item in expected_probabilities]
        maximum = max(raw_weights)
        expected_weights = [item / maximum for item in raw_weights]
    else:
        expected_probabilities = []
        expected_weights = []
    probability_ok = len(probabilities) == len(expected_probabilities) and all(
        _close(actual, expected) for actual, expected in zip(probabilities, expected_probabilities)
    )
    weights_ok = len(weights) == len(expected_weights) and all(
        _close(actual, expected) for actual, expected in zip(weights, expected_weights)
    )
    seed = value.get("rng_seed")
    return {
        "v54_replay_generation_is_closed": _generation_closed(value, ("priority_generation", "beta_generation", "seed_generation", "policy_generation", "owner_generation", "result_generation")),
        "v54_replay_priorities_and_alpha_are_bound": parameters_ok and value.get("result_priorities") == value.get("priorities") and value.get("result_priority_alpha") == alpha,
        "v54_replay_sampling_probability_is_recomputed": probability_ok and value.get("result_sampling_probability") == value.get("sampling_probability"),
        "v54_replay_beta_and_importance_weight_are_recomputed": weights_ok and value.get("result_importance_beta") == beta and value.get("result_importance_weight") == value.get("importance_weight"),
        "v54_replay_seed_is_nonnegative_and_bound": isinstance(seed, int) and not isinstance(seed, bool) and seed >= 0 and value.get("result_rng_seed") == seed,
        "v54_replay_policy_checkpoint_is_bound": _digest_equal(value, "policy_checkpoint_sha256", "result_policy_checkpoint_sha256"),
        "v54_replay_buffer_owner_is_bound": _prefixed_equal(value, "buffer_owner", "result_buffer_owner", "buffer:"),
        "v54_replay_result_digest_is_bound": _result_digest_bound(value),
    }


def _preprocess_checks(value: Mapping[str, object]) -> dict[str, bool]:
    means = _float_vector(value.get("feature_mean"))
    standard_deviations = _float_vector(value.get("feature_std"))
    normalization_ok = bool(means) and len(means) == len(standard_deviations) and all(item > 0.0 for item in standard_deviations)
    category_order = value.get("category_order")
    categories_ok = isinstance(category_order, Mapping) and bool(category_order)
    if categories_ok:
        for feature, categories in category_order.items():
            if not (
                isinstance(feature, str) and feature
                and isinstance(categories, Sequence) and not isinstance(categories, (str, bytes)) and bool(categories)
                and len(set(categories)) == len(categories)
                and all(isinstance(category, str) and category.startswith("category:") and len(category) > len("category:") for category in categories)
            ):
                categories_ok = False
                break
    split = value.get("data_split_indices")
    split_ok = isinstance(split, Mapping) and set(split) == {"train", "validation", "test"}
    partitions: list[list[int]] = []
    if split_ok:
        for name in ("train", "validation", "test"):
            indices = split[name]
            if not (
                isinstance(indices, Sequence) and not isinstance(indices, (str, bytes)) and bool(indices)
                and all(isinstance(index, int) and not isinstance(index, bool) and index >= 0 for index in indices)
                and len(set(indices)) == len(indices)
            ):
                split_ok = False
                break
            partitions.append(list(indices))
        if split_ok:
            combined = [index for partition in partitions for index in partition]
            split_ok = len(set(combined)) == len(combined) and sorted(combined) == list(range(len(combined)))
    return {
        "v54_preprocess_generation_is_closed": _generation_closed(value, ("normalization_generation", "category_generation", "split_generation", "model_generation", "owner_generation", "result_generation")),
        "v54_preprocess_normalization_is_valid_and_bound": normalization_ok and value.get("result_feature_mean") == value.get("feature_mean") and value.get("result_feature_std") == value.get("feature_std"),
        "v54_preprocess_category_order_is_unique_and_bound": categories_ok and value.get("result_category_order") == category_order,
        "v54_preprocess_split_is_disjoint_complete_and_bound": split_ok and value.get("result_data_split_indices") == split,
        "v54_preprocess_model_checkpoint_is_bound": _digest_equal(value, "model_checkpoint_sha256", "result_model_checkpoint_sha256"),
        "v54_preprocess_owner_is_bound": _prefixed_equal(value, "preprocessing_owner", "result_preprocessing_owner", "preprocess:"),
        "v54_preprocess_result_digest_is_bound": _result_digest_bound(value),
    }


def _float_vector(value: object) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value or not all(_finite(item) for item in value):
        return []
    return [float(item) for item in value]


def _generation_closed(value: Mapping[str, object], fields: Sequence[str]) -> bool:
    generation = str(value.get("generation", "")).strip()
    return bool(generation) and all(value.get(field) == generation for field in fields)


def _prefixed_equal(value: Mapping[str, object], left: str, right: str, prefix: str) -> bool:
    return str(value.get(left, "")).startswith(prefix) and value.get(left) == value.get(right)


def _digest_equal(value: Mapping[str, object], left: str, right: str) -> bool:
    return _valid_digest(value.get(left)) and value.get(right) == value.get(left)


def _result_digest_bound(value: Mapping[str, object]) -> bool:
    return _digest_equal(value, "result_sha256", "accepted_result_sha256")


def _valid_digest(value: object) -> bool:
    return isinstance(value, str) and bool(_DIGEST.fullmatch(value))


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=_TOL, abs_tol=_TOL)
