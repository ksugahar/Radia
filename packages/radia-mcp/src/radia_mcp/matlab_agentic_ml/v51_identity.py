"""Identity gates for RL replay, classification, edits, and parallel RNG."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence


_DIGEST = re.compile(r"^[0-9a-f]{64}$")


def validate_matlab_ml_rl_v51_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v51_identity")
    if not isinstance(identity, Mapping):
        return None
    records = {name: identity.get(name) for name in ("rl_replay", "classification", "agentic_edit", "parallel_rng")}
    checks = {f"v51_{name}_record_is_mapping": isinstance(value, Mapping) for name, value in records.items()}
    if isinstance(records["rl_replay"], Mapping):
        checks.update(_rl_checks(records["rl_replay"]))
    if isinstance(records["classification"], Mapping):
        checks.update(_classification_checks(records["classification"]))
    if isinstance(records["agentic_edit"], Mapping):
        checks.update(_edit_checks(records["agentic_edit"]))
    if isinstance(records["parallel_rng"], Mapping):
        checks.update(_parallel_checks(records["parallel_rng"]))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v51",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def _rl_checks(value: Mapping[str, object]) -> dict[str, bool]:
    transitions = value.get("transition_ids")
    terminal = value.get("terminal_flags")
    truncation = value.get("truncation_flags")
    discounts = value.get("discounts")
    priorities = value.get("priorities")
    transition_ok = isinstance(transitions, list) and bool(transitions) and len(transitions) == len(set(transitions)) and all(isinstance(item, str) and item.startswith("transition:") for item in transitions)
    count = len(transitions) if transition_ok else 0
    flags_ok = all(isinstance(items, list) and len(items) == count and all(isinstance(item, bool) for item in items) for items in (terminal, truncation))
    discount_ok = isinstance(discounts, list) and len(discounts) == count and all(_finite(item) and 0.0 <= float(item) <= 1.0 for item in discounts)
    priority_ok = isinstance(priorities, list) and len(priorities) == count and all(_finite(item) and float(item) > 0.0 for item in priorities)
    semantics_ok = flags_ok and discount_ok and all(not (done and truncated) and (not done or math.isclose(float(discount), 0.0, abs_tol=0.0)) for done, truncated, discount in zip(terminal, truncation, discounts))
    return {
        "v51_rl_generation_is_closed": _generation_closed(value, ("transition_generation", "terminal_generation", "discount_generation", "priority_generation", "target_generation", "owner_generation", "result_generation")),
        "v51_rl_transitions_are_ordered_and_bound": transition_ok and value.get("result_transition_ids") == transitions,
        "v51_rl_terminal_truncation_discount_semantics_are_bound": semantics_ok and value.get("result_terminal_flags") == terminal and value.get("result_truncation_flags") == truncation and value.get("result_discounts") == discounts,
        "v51_rl_priorities_are_positive_and_bound": priority_ok and value.get("result_priorities") == priorities,
        "v51_rl_target_network_is_bound": _valid_digest(value.get("target_network_sha256")) and value.get("result_target_network_sha256") == value.get("target_network_sha256"),
        "v51_rl_buffer_owner_is_bound": _prefixed_equal(value, "buffer_owner", "result_buffer_owner", "buffer:"),
        "v51_rl_result_digest_is_bound": _digest_bound(value),
    }


def _classification_checks(value: Mapping[str, object]) -> dict[str, bool]:
    classes = value.get("class_order")
    encoding = value.get("label_encoding")
    splits = value.get("data_split_sha256")
    matrix = value.get("confusion_matrix")
    axes = value.get("confusion_matrix_axes")
    classes_ok = isinstance(classes, list) and len(classes) >= 2 and len(classes) == len(set(classes)) and all(isinstance(item, str) and bool(item) for item in classes)
    encoding_ok = isinstance(encoding, Mapping) and classes_ok and list(encoding) == classes and list(encoding.values()) == list(range(len(classes)))
    splits_ok = isinstance(splits, Mapping) and set(splits) == {"train", "validation", "test"} and all(_valid_digest(item) for item in splits.values()) and len(set(splits.values())) == 3
    matrix_ok = isinstance(matrix, list) and classes_ok and len(matrix) == len(classes) and all(isinstance(row, list) and len(row) == len(classes) and all(_integer(item) and item >= 0 for item in row) for row in matrix)
    axes_ok = isinstance(axes, Mapping) and axes.get("rows") == classes and axes.get("columns") == classes
    return {
        "v51_classification_generation_is_closed": _generation_closed(value, ("class_generation", "label_generation", "split_generation", "confusion_generation", "owner_generation", "result_generation")),
        "v51_classification_class_order_is_bound": classes_ok and value.get("result_class_order") == classes,
        "v51_classification_label_encoding_is_bound": encoding_ok and value.get("result_label_encoding") == encoding,
        "v51_classification_split_is_bound": splits_ok and value.get("result_data_split_sha256") == splits,
        "v51_classification_confusion_axes_are_bound": matrix_ok and axes_ok and value.get("result_confusion_matrix") == matrix and value.get("result_confusion_matrix_axes") == axes,
        "v51_classification_model_owner_is_bound": _prefixed_equal(value, "model_owner", "result_model_owner", "model:"),
        "v51_classification_result_digest_is_bound": _digest_bound(value),
    }


def _edit_checks(value: Mapping[str, object]) -> dict[str, bool]:
    path = value.get("target_path")
    precondition = value.get("file_precondition_sha256")
    rollback = value.get("rollback_sha256")
    path_ok = isinstance(path, str) and path.startswith("workspace:project/") and ".." not in path.replace("\\", "/").split("/")
    return {
        "v51_edit_generation_is_closed": _generation_closed(value, ("precondition_generation", "patch_generation", "rollback_generation", "tool_call_generation", "owner_generation", "result_generation")),
        "v51_edit_target_is_safe_and_bound": path_ok and value.get("result_target_path") == path,
        "v51_edit_precondition_is_bound": _valid_digest(precondition) and value.get("result_file_precondition_sha256") == precondition,
        "v51_edit_patch_is_bound": _valid_digest(value.get("patch_sha256")) and value.get("result_patch_sha256") == value.get("patch_sha256"),
        "v51_edit_rollback_restores_precondition": value.get("rollback_state") == value.get("result_rollback_state") == "available" and _valid_digest(rollback) and rollback == precondition and value.get("result_rollback_sha256") == rollback,
        "v51_edit_tool_call_owner_is_bound": _prefixed_equal(value, "tool_call_owner", "result_tool_call_owner", "tool-call:"),
        "v51_edit_result_digest_is_bound": _digest_bound(value),
    }


def _parallel_checks(value: Mapping[str, object]) -> dict[str, bool]:
    substreams = value.get("substream_ids")
    workers = value.get("worker_map")
    reduction = value.get("reduction_order")
    substreams_ok = isinstance(substreams, list) and bool(substreams) and len(substreams) == len(set(substreams)) and all(_integer(item) and item > 0 for item in substreams)
    workers_ok = isinstance(workers, list) and substreams_ok and len(workers) == len(substreams) and all(_integer(item) and item > 0 for item in workers)
    reduction_ok = isinstance(reduction, list) and substreams_ok and len(reduction) == len(substreams) and set(reduction) == set(substreams)
    seed = value.get("rng_seed")
    return {
        "v51_parallel_generation_is_closed": _generation_closed(value, ("rng_generation", "substream_generation", "worker_generation", "reduction_generation", "checkpoint_generation", "owner_generation", "result_generation")),
        "v51_parallel_rng_is_bound": value.get("rng_algorithm") in {"Threefry", "Philox"} and value.get("result_rng_algorithm") == value.get("rng_algorithm") and _integer(seed) and seed >= 0 and value.get("result_rng_seed") == seed,
        "v51_parallel_substreams_are_bound": substreams_ok and value.get("result_substream_ids") == substreams,
        "v51_parallel_worker_map_is_bound": workers_ok and value.get("result_worker_map") == workers,
        "v51_parallel_reduction_order_is_bound": reduction_ok and value.get("result_reduction_order") == reduction,
        "v51_parallel_checkpoint_is_bound": _valid_digest(value.get("checkpoint_sha256")) and value.get("result_checkpoint_sha256") == value.get("checkpoint_sha256"),
        "v51_parallel_pool_owner_is_bound": _prefixed_equal(value, "pool_owner", "result_pool_owner", "pool:"),
        "v51_parallel_result_digest_is_bound": _digest_bound(value),
    }


def _generation_closed(value: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(value.get("generation", "")).strip()
    return bool(generation) and all(value.get(field) == generation for field in fields)


def _prefixed_equal(value: Mapping[str, object], left: str, right: str, prefix: str) -> bool:
    return str(value.get(left, "")).startswith(prefix) and value.get(left) == value.get(right)


def _valid_digest(value: object) -> bool:
    return bool(_DIGEST.fullmatch(str(value)))


def _digest_bound(value: Mapping[str, object]) -> bool:
    return _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256")


def _integer(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))
