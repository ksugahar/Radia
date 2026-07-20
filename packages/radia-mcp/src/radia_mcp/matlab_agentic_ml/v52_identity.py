"""Replay-identity gates for MATLAB training, recurrent RL, codegen, and SDI."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping


_DIGEST = re.compile(r"^[0-9a-f]{64}$")


def validate_matlab_ml_rl_v52_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate v52 ML/RL and agent-tool replay identities."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v52_identity")
    if not isinstance(identity, Mapping):
        return None
    records = {
        name: identity.get(name)
        for name in ("ml_training", "recurrent_rl", "codegen", "simulink_data_inspector")
    }
    checks = {
        f"v52_{name}_record_is_mapping": isinstance(value, Mapping)
        for name, value in records.items()
    }
    if isinstance(records["ml_training"], Mapping):
        checks.update(_training_checks(records["ml_training"]))
    if isinstance(records["recurrent_rl"], Mapping):
        checks.update(_recurrent_checks(records["recurrent_rl"]))
    if isinstance(records["codegen"], Mapping):
        checks.update(_codegen_checks(records["codegen"]))
    if isinstance(records["simulink_data_inspector"], Mapping):
        checks.update(_sdi_checks(records["simulink_data_inspector"]))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v52",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def _training_checks(value: Mapping[str, object]) -> dict[str, bool]:
    threshold = value.get("gradient_clip_threshold")
    step = value.get("optimizer_step")
    schedule = value.get("learning_rate_schedule")
    schedule_ok = (
        isinstance(schedule, list)
        and bool(schedule)
        and all(_finite(item) and float(item) > 0.0 for item in schedule)
        and all(float(left) >= float(right) for left, right in zip(schedule, schedule[1:]))
    )
    return {
        "v52_training_generation_is_closed": _generation_closed(value, ("gradient_generation", "optimizer_generation", "schedule_generation", "checkpoint_generation", "owner_generation", "result_generation")),
        "v52_training_gradient_clip_is_bound": _finite(threshold) and float(threshold) > 0.0 and value.get("gradient_clip_method") in {"global-l2norm", "global-l1norm", "elementwise-value"} and value.get("result_gradient_clip_threshold") == threshold and value.get("result_gradient_clip_method") == value.get("gradient_clip_method"),
        "v52_training_optimizer_state_is_bound": value.get("optimizer_name") in {"adam", "sgdm", "rmsprop"} and _integer(step) and step >= 0 and value.get("result_optimizer_name") == value.get("optimizer_name") and value.get("result_optimizer_step") == step and _digest_equal(value, "optimizer_state_sha256", "result_optimizer_state_sha256"),
        "v52_training_schedule_is_monotone_and_bound": schedule_ok and value.get("result_learning_rate_schedule") == schedule,
        "v52_training_checkpoint_is_bound": _digest_equal(value, "checkpoint_sha256", "result_checkpoint_sha256"),
        "v52_training_model_owner_is_bound": _prefixed_equal(value, "model_owner", "result_model_owner", "model:"),
        "v52_training_result_digest_is_bound": _result_digest_bound(value),
    }


def _recurrent_checks(value: Mapping[str, object]) -> dict[str, bool]:
    hidden = value.get("hidden_state_sha256")
    reset = value.get("episode_reset")
    mask = value.get("sequence_mask")
    hidden_ok = isinstance(hidden, list) and bool(hidden) and all(_valid_digest(item) for item in hidden)
    count = len(hidden) if hidden_ok else 0
    reset_ok = isinstance(reset, list) and len(reset) == count and all(isinstance(item, bool) for item in reset) and reset[0] and sum(reset) == 1
    mask_ok = isinstance(mask, list) and len(mask) == count and all(isinstance(item, bool) for item in mask) and mask[0] and all(left or not right for left, right in zip(mask, mask[1:]))
    return {
        "v52_recurrent_generation_is_closed": _generation_closed(value, ("hidden_generation", "reset_generation", "mask_generation", "policy_generation", "owner_generation", "result_generation")),
        "v52_recurrent_hidden_state_is_bound": hidden_ok and value.get("result_hidden_state_sha256") == hidden,
        "v52_recurrent_episode_reset_is_initial_and_bound": reset_ok and value.get("result_episode_reset") == reset,
        "v52_recurrent_sequence_mask_is_prefix_and_bound": mask_ok and value.get("result_sequence_mask") == mask,
        "v52_recurrent_policy_is_bound": _digest_equal(value, "policy_sha256", "result_policy_sha256"),
        "v52_recurrent_policy_owner_is_bound": _prefixed_equal(value, "policy_owner", "result_policy_owner", "policy:"),
        "v52_recurrent_result_digest_is_bound": _result_digest_bound(value),
    }


def _codegen_checks(value: Mapping[str, object]) -> dict[str, bool]:
    target = value.get("target_hardware")
    return {
        "v52_codegen_generation_is_closed": _generation_closed(value, ("target_generation", "numeric_generation", "config_generation", "build_generation", "owner_generation", "result_generation")),
        "v52_codegen_target_is_bound": isinstance(target, str) and bool(target.strip()) and value.get("result_target_hardware") == target,
        "v52_codegen_numeric_type_is_bound": value.get("numeric_type") in {"double", "single", "fixed-point"} and value.get("result_numeric_type") == value.get("numeric_type"),
        "v52_codegen_configuration_is_bound": _digest_equal(value, "configuration_sha256", "result_configuration_sha256"),
        "v52_codegen_build_is_bound": _digest_equal(value, "build_sha256", "result_build_sha256"),
        "v52_codegen_owner_is_bound": _prefixed_equal(value, "build_owner", "result_build_owner", "build:"),
        "v52_codegen_result_digest_is_bound": _result_digest_bound(value),
    }


def _sdi_checks(value: Mapping[str, object]) -> dict[str, bool]:
    hierarchy = value.get("signal_hierarchy")
    units = value.get("signal_units")
    interpolation = value.get("interpolation")
    hierarchy_ok = isinstance(hierarchy, list) and bool(hierarchy) and len(hierarchy) == len(set(hierarchy)) and all(isinstance(item, str) and "/" in item and not item.startswith("/") for item in hierarchy)
    count = len(hierarchy) if hierarchy_ok else 0
    units_ok = isinstance(units, list) and len(units) == count and all(isinstance(item, str) and bool(item.strip()) for item in units)
    interpolation_ok = isinstance(interpolation, list) and len(interpolation) == count and all(item in {"linear", "zoh"} for item in interpolation)
    return {
        "v52_sdi_generation_is_closed": _generation_closed(value, ("run_generation", "signal_generation", "unit_generation", "interpolation_generation", "owner_generation", "result_generation")),
        "v52_sdi_run_is_bound": _prefixed_equal(value, "run_id", "result_run_id", "run:"),
        "v52_sdi_hierarchy_is_unique_and_bound": hierarchy_ok and value.get("result_signal_hierarchy") == hierarchy,
        "v52_sdi_units_are_bound": units_ok and value.get("result_signal_units") == units,
        "v52_sdi_interpolation_is_supported_and_bound": interpolation_ok and value.get("result_interpolation") == interpolation,
        "v52_sdi_owner_is_bound": _prefixed_equal(value, "session_owner", "result_session_owner", "session:"),
        "v52_sdi_result_digest_is_bound": _result_digest_bound(value),
    }


def _generation_closed(value: Mapping[str, object], fields: tuple[str, ...]) -> bool:
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


def _integer(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))
