"""Identity gates for MATLAB optimization, training, and agentic artifacts."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping


_DIGEST = re.compile(r"^[0-9a-f]{64}$")


def validate_matlab_ml_rl_v50_identity(
    summary: Mapping[str, object],
) -> dict[str, object] | None:
    """Bind optimization/training results to the state that produced them."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v50_identity")
    if not isinstance(identity, Mapping):
        return None
    records = {
        name: identity.get(name)
        for name in (
            "bayesopt",
            "mixed_precision",
            "agentic_sandbox",
            "tall_datastore",
        )
    }
    checks = {
        f"v50_{name}_record_is_mapping": isinstance(value, Mapping)
        for name, value in records.items()
    }
    if isinstance(records["bayesopt"], Mapping):
        checks.update(_bayesopt_checks(records["bayesopt"]))
    if isinstance(records["mixed_precision"], Mapping):
        checks.update(_mixed_precision_checks(records["mixed_precision"]))
    if isinstance(records["agentic_sandbox"], Mapping):
        checks.update(_agentic_sandbox_checks(records["agentic_sandbox"]))
    if isinstance(records["tall_datastore"], Mapping):
        checks.update(_tall_datastore_checks(records["tall_datastore"]))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v50",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def _bayesopt_checks(value: Mapping[str, object]) -> dict[str, bool]:
    constraints = value.get("constraints")
    acquisition = value.get("acquisition_state")
    seed = value.get("rng_seed")
    noise = value.get("objective_noise_sigma")
    return {
        "v50_bayesopt_generation_is_closed": _generation_closed(
            value,
            (
                "objective_generation",
                "noise_generation",
                "constraint_generation",
                "seed_generation",
                "acquisition_generation",
                "model_generation",
                "result_generation",
            ),
        ),
        "v50_bayesopt_objective_is_bound": _nonempty_equal(
            value, "objective_name", "result_objective_name"
        ),
        "v50_bayesopt_noise_is_nonnegative_and_bound": _finite_number(noise)
        and float(noise) >= 0.0
        and value.get("result_objective_noise_sigma") == noise,
        "v50_bayesopt_constraints_are_explicit_and_bound": _valid_constraints(
            constraints
        )
        and value.get("result_constraints") == constraints,
        "v50_bayesopt_seed_is_bound": _integer(seed)
        and int(seed) >= 0
        and value.get("result_rng_seed") == seed,
        "v50_bayesopt_acquisition_state_is_bound": _valid_acquisition(acquisition)
        and value.get("result_acquisition_state") == acquisition,
        "v50_bayesopt_model_owner_is_bound": _prefixed_equal(
            value, "model_owner", "result_model_owner", "model:"
        ),
        "v50_bayesopt_result_digest_is_bound": _digest_bound(value),
    }


def _mixed_precision_checks(value: Mapping[str, object]) -> dict[str, bool]:
    policy = value.get("precision_policy")
    loss_scale = value.get("loss_scale")
    optimizer = value.get("optimizer_state")
    checkpoint = value.get("checkpoint_sha256")
    return {
        "v50_mixed_precision_generation_is_closed": _generation_closed(
            value,
            (
                "precision_generation",
                "loss_scale_generation",
                "optimizer_generation",
                "checkpoint_generation",
                "network_generation",
                "result_generation",
            ),
        ),
        "v50_mixed_precision_policy_is_bound": policy
        in {"mixed-fp16", "mixed-bfloat16"}
        and value.get("result_precision_policy") == policy,
        "v50_mixed_precision_loss_scale_is_positive_and_bound": _finite_number(
            loss_scale
        )
        and float(loss_scale) > 0.0
        and value.get("result_loss_scale") == loss_scale,
        "v50_mixed_precision_optimizer_state_is_bound": _valid_optimizer(optimizer)
        and value.get("result_optimizer_state") == optimizer,
        "v50_mixed_precision_checkpoint_is_bound": _valid_digest(checkpoint)
        and value.get("result_checkpoint_sha256") == checkpoint,
        "v50_mixed_precision_network_owner_is_bound": _prefixed_equal(
            value, "network_owner", "result_network_owner", "network:"
        ),
        "v50_mixed_precision_result_digest_is_bound": _digest_bound(value),
    }


def _agentic_sandbox_checks(value: Mapping[str, object]) -> dict[str, bool]:
    roots = value.get("allowed_roots")
    requested = value.get("requested_path")
    resolved = value.get("resolved_path")
    roots_ok = (
        isinstance(roots, list)
        and bool(roots)
        and len(roots) == len(set(roots))
        and all(isinstance(root, str) and root.startswith("workspace:") for root in roots)
    )
    paths_ok = roots_ok and all(
        isinstance(path, str)
        and ".." not in path.replace("\\", "/").split("/")
        and any(path == root or path.startswith(root + "/") for root in roots)
        for path in (requested, resolved)
    )
    return {
        "v50_agentic_sandbox_generation_is_closed": _generation_closed(
            value,
            (
                "root_generation",
                "resolution_generation",
                "traversal_generation",
                "approval_generation",
                "tool_call_generation",
                "result_generation",
            ),
        ),
        "v50_agentic_allowed_roots_are_bound": roots_ok
        and value.get("result_allowed_roots") == roots,
        "v50_agentic_paths_stay_inside_allowed_roots": paths_ok
        and value.get("result_requested_path") == requested
        and value.get("result_resolved_path") == resolved,
        "v50_agentic_symlink_and_traversal_are_safe": value.get(
            "symlink_resolution"
        )
        == value.get("result_symlink_resolution")
        == "inside_allowed_root"
        and value.get("traversal_detected") is False
        and value.get("result_traversal_detected") is False,
        "v50_agentic_approval_is_bound": _prefixed_equal(
            value, "approval_id", "result_approval_id", "approval:"
        ),
        "v50_agentic_tool_call_owner_is_bound": _prefixed_equal(
            value, "tool_call_owner", "result_tool_call_owner", "tool-call:"
        ),
        "v50_agentic_result_digest_is_bound": _digest_bound(value),
    }


def _tall_datastore_checks(value: Mapping[str, object]) -> dict[str, bool]:
    partitions = value.get("partition_ids")
    workers = value.get("worker_order")
    cursor = value.get("resume_cursor")
    checkpoint = value.get("checkpoint_sha256")
    partitions_ok = (
        isinstance(partitions, list)
        and bool(partitions)
        and len(partitions) == len(set(partitions))
        and all(isinstance(item, str) and item.startswith("partition:") for item in partitions)
    )
    workers_ok = (
        isinstance(workers, list)
        and partitions_ok
        and len(workers) == len(partitions)
        and all(_integer(worker) and int(worker) > 0 for worker in workers)
    )
    cursor_ok = isinstance(cursor, str) and bool(
        re.fullmatch(r"partition:(\d+)/row:(\d+)", cursor)
    )
    if cursor_ok:
        match = re.fullmatch(r"partition:(\d+)/row:(\d+)", cursor)
        cursor_ok = bool(match and f"partition:{match.group(1)}" in partitions)
    return {
        "v50_tall_generation_is_closed": _generation_closed(
            value,
            (
                "partition_generation",
                "worker_generation",
                "checkpoint_generation",
                "resume_generation",
                "datastore_generation",
                "result_generation",
            ),
        ),
        "v50_tall_partitions_are_unique_and_bound": partitions_ok
        and value.get("result_partition_ids") == partitions,
        "v50_tall_worker_order_is_bound": workers_ok
        and value.get("result_worker_order") == workers,
        "v50_tall_checkpoint_is_bound": _valid_digest(checkpoint)
        and value.get("result_checkpoint_sha256") == checkpoint,
        "v50_tall_resume_cursor_is_bound": cursor_ok
        and value.get("result_resume_cursor") == cursor,
        "v50_tall_datastore_owner_is_bound": _prefixed_equal(
            value, "datastore_owner", "result_datastore_owner", "datastore:"
        ),
        "v50_tall_result_digest_is_bound": _digest_bound(value),
    }


def _generation_closed(value: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(value.get("generation", "")).strip()
    return bool(generation) and all(value.get(field) == generation for field in fields)


def _valid_constraints(value: object) -> bool:
    return isinstance(value, list) and bool(value) and all(
        isinstance(item, Mapping)
        and bool(str(item.get("name", "")).strip())
        and item.get("sense") in {"<=", ">="}
        and _finite_number(item.get("limit"))
        for item in value
    )


def _valid_acquisition(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and value.get("name")
        in {"expected-improvement", "expected-improvement-plus", "lower-confidence-bound"}
        and _integer(value.get("iteration"))
        and int(value["iteration"]) > 0
        and str(value.get("incumbent", "")).startswith("trial:")
    )


def _valid_optimizer(value: object) -> bool:
    if not isinstance(value, Mapping):
        return False
    common = (
        value.get("name") in {"adam", "sgdm"}
        and _integer(value.get("iteration"))
        and int(value["iteration"]) >= 0
        and _finite_number(value.get("learn_rate"))
        and float(value["learn_rate"]) > 0.0
    )
    if value.get("name") == "adam":
        return common and _valid_digest(value.get("moment1_sha256")) and _valid_digest(
            value.get("moment2_sha256")
        )
    return common


def _nonempty_equal(value: Mapping[str, object], left: str, right: str) -> bool:
    return bool(str(value.get(left, "")).strip()) and value.get(left) == value.get(right)


def _prefixed_equal(
    value: Mapping[str, object], left: str, right: str, prefix: str
) -> bool:
    return str(value.get(left, "")).startswith(prefix) and value.get(left) == value.get(right)


def _valid_digest(value: object) -> bool:
    return isinstance(value, str) and bool(_DIGEST.fullmatch(value))


def _digest_bound(value: Mapping[str, object]) -> bool:
    return _valid_digest(value.get("result_sha256")) and value.get(
        "accepted_result_sha256"
    ) == value.get("result_sha256")


def _integer(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _finite_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(
        float(value)
    )
