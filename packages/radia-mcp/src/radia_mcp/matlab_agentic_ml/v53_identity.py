"""Replay-identity gates for offline RL and probability calibration artifacts."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence


_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_TOL = 1.0e-10


def validate_matlab_ml_rl_v53_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate v53 offline-RL and calibration replay identities."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v53_identity")
    if not isinstance(identity, Mapping):
        return None
    offline = identity.get("offline_rl")
    calibration = identity.get("probability_calibration")
    checks = {
        "v53_offline_rl_record_is_mapping": isinstance(offline, Mapping),
        "v53_probability_calibration_record_is_mapping": isinstance(calibration, Mapping),
    }
    if isinstance(offline, Mapping):
        checks.update(_offline_checks(offline))
    if isinstance(calibration, Mapping):
        checks.update(_calibration_checks(calibration))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v53",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def _offline_checks(value: Mapping[str, object]) -> dict[str, bool]:
    behavior = _float_vector(value.get("behavior_action_probability"))
    target = _float_vector(value.get("target_action_probability"))
    weights = _float_vector(value.get("importance_weight"))
    probability_ok = (
        bool(behavior)
        and len(behavior) == len(target)
        and all(item > 0.0 for item in behavior)
        and all(item >= 0.0 for item in target)
        and _close(sum(behavior), 1.0)
        and _close(sum(target), 1.0)
    )
    weight_ok = (
        probability_ok
        and len(weights) == len(behavior)
        and all(item > 0.0 for item in weights)
        and all(_close(weight, target_item / behavior_item) for weight, target_item, behavior_item in zip(weights, target, behavior))
    )
    coverage = value.get("support_coverage_fraction")
    expected_coverage = sum(item > 0.0 for item in behavior) / len(behavior) if behavior else math.nan
    coverage_ok = _finite(coverage) and 0.0 <= float(coverage) <= 1.0 and _close(float(coverage), expected_coverage)
    return {
        "v53_offline_generation_is_closed": _generation_closed(value, ("behavior_generation", "weight_generation", "support_generation", "dataset_generation", "owner_generation", "result_generation")),
        "v53_offline_policy_digests_are_bound": _digest_equal(value, "behavior_policy_sha256", "result_behavior_policy_sha256") and _digest_equal(value, "target_policy_sha256", "result_target_policy_sha256"),
        "v53_offline_probabilities_have_shared_support": probability_ok and value.get("result_behavior_action_probability") == value.get("behavior_action_probability") and value.get("result_target_action_probability") == value.get("target_action_probability"),
        "v53_offline_importance_weights_are_recomputed": weight_ok and value.get("result_importance_weight") == value.get("importance_weight"),
        "v53_offline_support_coverage_is_recomputed": coverage_ok and value.get("result_support_coverage_fraction") == coverage,
        "v53_offline_dataset_owner_is_bound": _prefixed_equal(value, "dataset_owner", "result_dataset_owner", "dataset:"),
        "v53_offline_result_digest_is_bound": _result_digest_bound(value),
    }


def _calibration_checks(value: Mapping[str, object]) -> dict[str, bool]:
    logits = _float_matrix(value.get("logits"))
    probabilities = _float_matrix(value.get("calibrated_probability"))
    class_order = value.get("class_order")
    temperature = value.get("temperature")
    width = len(logits[0]) if logits else 0
    logits_ok = bool(logits) and width > 0 and all(len(row) == width for row in logits)
    classes_ok = (
        isinstance(class_order, list)
        and len(class_order) == width
        and len(set(class_order)) == width
        and all(isinstance(item, str) and item.startswith("class:") and len(item) > 6 for item in class_order)
    )
    temperature_ok = _finite(temperature) and float(temperature) > 0.0
    probabilities_ok = (
        logits_ok
        and temperature_ok
        and len(probabilities) == len(logits)
        and all(len(row) == width and all(0.0 <= item <= 1.0 for item in row) and _close(sum(row), 1.0) for row in probabilities)
    )
    recomputed = [_softmax(row, float(temperature)) for row in logits] if logits_ok and temperature_ok else []
    recomputed_ok = probabilities_ok and all(
        all(_close(actual, expected) for actual, expected in zip(actual_row, expected_row))
        for actual_row, expected_row in zip(probabilities, recomputed)
    )
    return {
        "v53_calibration_generation_is_closed": _generation_closed(value, ("logit_generation", "probability_generation", "class_generation", "temperature_generation", "owner_generation", "result_generation")),
        "v53_calibration_logits_are_bound": logits_ok and value.get("result_logits") == value.get("logits"),
        "v53_calibration_class_order_is_unique_and_bound": classes_ok and value.get("result_class_order") == class_order,
        "v53_calibration_temperature_is_positive_and_bound": temperature_ok and value.get("result_temperature") == temperature,
        "v53_calibration_probabilities_are_recomputed": recomputed_ok and value.get("result_calibrated_probability") == value.get("calibrated_probability"),
        "v53_calibration_model_owner_is_bound": _prefixed_equal(value, "model_owner", "result_model_owner", "model:"),
        "v53_calibration_result_digest_is_bound": _result_digest_bound(value),
    }


def _softmax(row: list[float], temperature: float) -> list[float]:
    scaled = [item / temperature for item in row]
    maximum = max(scaled)
    exponentials = [math.exp(item - maximum) for item in scaled]
    total = sum(exponentials)
    return [item / total for item in exponentials]


def _float_vector(value: object) -> list[float]:
    if not isinstance(value, list) or not all(_finite(item) for item in value):
        return []
    return [float(item) for item in value]


def _float_matrix(value: object) -> list[list[float]]:
    if not isinstance(value, list) or not value:
        return []
    rows = [_float_vector(row) for row in value]
    return rows if all(rows) else []


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
    return bool(_DIGEST.fullmatch(str(value)))


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=_TOL, abs_tol=_TOL)
