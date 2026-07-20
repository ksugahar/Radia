"""Small, solver-neutral contract gate for MATLAB ML/RL result artifacts."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping


_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_TASKS = {"regression", "classification", "reinforcement_learning"}


def _integer(value: object, *, minimum: int) -> bool:
    return (
        isinstance(value, int)
        and not isinstance(value, bool)
        and value >= minimum
    )


def validate_matlab_ml_rl_artifact(value: Mapping[str, object]) -> dict[str, object]:
    """Check provenance, split/evaluation separation, and task-specific metrics.

    This is an artifact gate, not a trainer. It deliberately refuses to infer
    validation from a training score or from an unpinned MATLAB session.
    """
    if not isinstance(value, Mapping):
        raise TypeError("artifact must be a mapping")
    task = value.get("task")
    checks: dict[str, bool] = {
        "schema_is_pinned": value.get("schema") == "cae-ai-lab.matlab-ml-rl-result.v1",
        "task_is_supported": task in _TASKS,
        "matlab_release_is_pinned": bool(str(value.get("matlab_release", "")).strip()),
        "session_owner_is_recorded": str(value.get("session_owner", "")).startswith("matlab:"),
        "seed_is_recorded": _integer(value.get("random_seed"), minimum=0),
        "result_digest_is_valid": bool(_DIGEST.fullmatch(str(value.get("result_sha256", "")))),
        "accepted_digest_matches": value.get("accepted_result_sha256") == value.get("result_sha256"),
        "training_and_evaluation_are_disjoint": value.get("training_evaluation_disjoint") is True,
    }
    if task in {"regression", "classification"}:
        folds = value.get("cross_validation_folds")
        checks["cross_validation_is_explicit"] = _integer(folds, minimum=2)
        checks["holdout_or_cv_score_is_finite"] = math.isfinite(float(value.get("validation_score"))) if _number(value.get("validation_score")) else False
    elif task == "reinforcement_learning":
        checks["environment_is_pinned"] = bool(str(value.get("environment_id", "")).strip())
        checks["training_episodes_are_positive"] = _integer(
            value.get("training_episodes"), minimum=1
        )
        checks["evaluation_episodes_are_positive"] = _integer(
            value.get("evaluation_episodes"), minimum=1
        )
        checks["evaluation_return_is_finite"] = math.isfinite(float(value.get("evaluation_mean_return"))) if _number(value.get("evaluation_mean_return")) else False
    else:
        checks["task_specific_contract_is_present"] = False
    return {
        "policy": "matlab_ml_rl_artifact_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def validate_matlab_ml_rl_v44_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate the v44 split/evaluation identity carried by a solver summary.

    The gate is deliberately independent of a trainer.  It checks the two
    durable promises that an MCP can verify after execution: held-out data are
    disjoint for supervised learning, and RL evaluation uses the same pinned
    environment with a fresh evaluation seed and no exploration.
    """
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v44_identity")
    if not isinstance(identity, Mapping):
        return None

    checks: dict[str, bool] = {}
    supervised = identity.get("supervised")
    rl = identity.get("reinforcement_learning")
    checks["supervised_record_is_mapping"] = isinstance(supervised, Mapping)
    checks["rl_record_is_mapping"] = isinstance(rl, Mapping)
    if isinstance(supervised, Mapping):
        checks.update(_supervised_v44_checks(supervised))
    if isinstance(rl, Mapping):
        checks.update(_rl_v44_checks(rl))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v44",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def validate_matlab_ml_rl_v45_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate v45 official ML/RL and Agentic Toolkit replay bindings."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v45_identity")
    if not isinstance(identity, Mapping):
        return None
    records = {
        "supervised": identity.get("supervised"),
        "reinforcement_learning": identity.get("reinforcement_learning"),
        "agentic_toolkit": identity.get("agentic_toolkit"),
        "mlrl_checkpoint": identity.get("mlrl_checkpoint"),
    }
    checks = {
        "v45_supervised_record_is_mapping": isinstance(records["supervised"], Mapping),
        "v45_rl_record_is_mapping": isinstance(records["reinforcement_learning"], Mapping),
        "v45_agentic_record_is_mapping": isinstance(records["agentic_toolkit"], Mapping),
        "v45_checkpoint_record_is_mapping": isinstance(records["mlrl_checkpoint"], Mapping),
    }
    if isinstance(records["supervised"], Mapping):
        checks.update(_supervised_v45_checks(records["supervised"]))
    if isinstance(records["reinforcement_learning"], Mapping):
        checks.update(_rl_v45_checks(records["reinforcement_learning"]))
    if isinstance(records["agentic_toolkit"], Mapping):
        checks.update(_agentic_v45_checks(records["agentic_toolkit"]))
    if isinstance(records["mlrl_checkpoint"], Mapping):
        checks.update(_checkpoint_v45_checks(records["mlrl_checkpoint"]))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v45",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def validate_matlab_ml_rl_v46_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate official ML/RL replay finiteness, seed, restart, and agent-call bindings."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v46_identity")
    if not isinstance(identity, Mapping):
        return None
    records = {name: identity.get(name) for name in ("supervised", "reinforcement_learning", "agentic_toolkit", "mlrl_checkpoint")}
    checks = {f"v46_{name}_record_is_mapping": isinstance(value, Mapping) for name, value in records.items()}
    if isinstance(records["supervised"], Mapping):
        checks.update(_supervised_v46_checks(records["supervised"]))
    if isinstance(records["reinforcement_learning"], Mapping):
        checks.update(_rl_v46_checks(records["reinforcement_learning"]))
    if isinstance(records["agentic_toolkit"], Mapping):
        checks.update(_agentic_v46_checks(records["agentic_toolkit"]))
    if isinstance(records["mlrl_checkpoint"], Mapping):
        checks.update(_checkpoint_v46_checks(records["mlrl_checkpoint"]))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v46",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def validate_matlab_ml_rl_v47_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate ML/RL data, environment, agent-call, and experiment lineage."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v47_identity")
    if not isinstance(identity, Mapping):
        return None
    records = {
        name: identity.get(name)
        for name in ("ml_datastore", "reinforcement_learning", "agentic_toolkit", "experiment_trials")
    }
    checks = {f"v47_{name}_record_is_mapping": isinstance(value, Mapping) for name, value in records.items()}
    if isinstance(records["ml_datastore"], Mapping):
        checks.update(_ml_datastore_v47_checks(records["ml_datastore"]))
    if isinstance(records["reinforcement_learning"], Mapping):
        checks.update(_rl_v47_checks(records["reinforcement_learning"]))
    if isinstance(records["agentic_toolkit"], Mapping):
        checks.update(_agentic_v47_checks(records["agentic_toolkit"]))
    if isinstance(records["experiment_trials"], Mapping):
        checks.update(_experiment_v47_checks(records["experiment_trials"]))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v47",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def validate_matlab_ml_rl_v48_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate autodiff, sequence, agentic-workspace, and experiment selection lineage."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v48_identity")
    if not isinstance(identity, Mapping):
        return None
    records = {name: identity.get(name) for name in ("autodiff", "sequence_model", "agentic_workspace", "experiment_selection")}
    checks = {f"v48_{name}_record_is_mapping": isinstance(value, Mapping) for name, value in records.items()}
    if isinstance(records["autodiff"], Mapping):
        checks.update(_autodiff_v48_checks(records["autodiff"]))
    if isinstance(records["sequence_model"], Mapping):
        checks.update(_sequence_v48_checks(records["sequence_model"]))
    if isinstance(records["agentic_workspace"], Mapping):
        checks.update(_agentic_workspace_v48_checks(records["agentic_workspace"]))
    if isinstance(records["experiment_selection"], Mapping):
        checks.update(_experiment_selection_v48_checks(records["experiment_selection"]))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v48",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def validate_matlab_ml_rl_v49_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate replay-buffer, normalization, cancellation, and parallel-resume identity."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v49_identity")
    if not isinstance(identity, Mapping):
        return None
    records = {
        name: identity.get(name)
        for name in ("rl_replay", "ml_normalization", "agentic_cancel", "parallel_resume")
    }
    checks = {f"v49_{name}_record_is_mapping": isinstance(value, Mapping) for name, value in records.items()}
    if isinstance(records["rl_replay"], Mapping):
        checks.update(_rl_replay_v49_checks(records["rl_replay"]))
    if isinstance(records["ml_normalization"], Mapping):
        checks.update(_ml_normalization_v49_checks(records["ml_normalization"]))
    if isinstance(records["agentic_cancel"], Mapping):
        checks.update(_agentic_cancel_v49_checks(records["agentic_cancel"]))
    if isinstance(records["parallel_resume"], Mapping):
        checks.update(_parallel_resume_v49_checks(records["parallel_resume"]))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v49",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def _rl_replay_v49_checks(value: Mapping[str, object]) -> dict[str, bool]:
    rows = value.get("replay_row_keys")
    observations = value.get("observations")
    actions = value.get("actions")
    rewards = value.get("rewards")
    terminals = value.get("terminals")
    seeds = value.get("episode_seeds")
    count = len(rows) if isinstance(rows, list) else 0
    observations_ok = (
        isinstance(observations, list)
        and len(observations) == count
        and bool(observations)
        and all(_finite_list(row) and len(row) == len(observations[0]) for row in observations)
    )
    return {
        "v49_rl_replay_generation_is_closed": _generation_v47(value, ("buffer_generation", "observation_generation", "action_generation", "reward_generation", "terminal_generation", "seed_generation", "checkpoint_generation", "policy_generation", "result_generation")),
        "v49_rl_replay_rows_are_unique_and_bound": isinstance(rows, list) and bool(rows) and len(rows) == len(set(rows)) and value.get("result_replay_row_keys") == rows,
        "v49_rl_observation_action_reward_rows_are_bound": observations_ok and value.get("result_observations") == observations and _finite_list(actions, count) and value.get("result_actions") == actions and _finite_list(rewards, count) and value.get("result_rewards") == rewards,
        "v49_rl_terminal_rows_are_bound": isinstance(terminals, list) and len(terminals) == count and all(isinstance(item, bool) for item in terminals) and any(terminals) and value.get("result_terminals") == terminals,
        "v49_rl_episode_seeds_are_unique_and_bound": isinstance(seeds, list) and len(seeds) == count and len(seeds) == len(set(seeds)) and all(isinstance(seed, int) and not isinstance(seed, bool) and seed >= 0 for seed in seeds) and value.get("result_episode_seeds") == seeds,
        "v49_rl_checkpoint_policy_owners_are_bound": str(value.get("checkpoint_owner", "")).startswith("checkpoint:") and value.get("result_checkpoint_owner") == value.get("checkpoint_owner") and str(value.get("policy_owner", "")).startswith("policy:") and value.get("result_policy_owner") == value.get("policy_owner"),
        "v49_rl_result_digest_is_bound": _digest_v47(value),
    }


def _ml_normalization_v49_checks(value: Mapping[str, object]) -> dict[str, bool]:
    splits = value.get("split_row_keys")
    train = splits.get("train") if isinstance(splits, Mapping) else None
    validation = splits.get("validation") if isinstance(splits, Mapping) else None
    test = splits.get("test") if isinstance(splits, Mapping) else None
    split_lists = (train, validation, test)
    all_rows = [row for group in split_lists if isinstance(group, list) for row in group]
    encoding = value.get("class_encoding")
    folds = value.get("training_fold_ids")
    metric_rows = value.get("metric_row_keys")
    metric_values = value.get("metric_values")
    encoding_values = list(encoding.values()) if isinstance(encoding, Mapping) else []
    return {
        "v49_ml_normalization_generation_is_closed": _generation_v47(value, ("normalization_generation", "class_generation", "split_generation", "fold_generation", "metric_generation", "model_generation", "result_generation")),
        "v49_ml_split_is_disjoint_and_bound": isinstance(splits, Mapping) and set(splits) == {"train", "validation", "test"} and all(isinstance(group, list) and bool(group) for group in split_lists) and len(all_rows) == len(set(all_rows)) and value.get("result_split_row_keys") == splits,
        "v49_ml_normalization_fits_training_only": value.get("normalization_fit_scope") == value.get("result_normalization_fit_scope") == "training_partition_only" and value.get("normalization_fit_rows") == train and value.get("result_normalization_fit_rows") == train,
        "v49_ml_class_encoding_is_canonical_and_bound": isinstance(encoding, Mapping) and bool(encoding) and all(isinstance(name, str) and name for name in encoding) and all(isinstance(index, int) and not isinstance(index, bool) and index >= 0 for index in encoding_values) and sorted(encoding_values) == list(range(len(encoding_values))) and value.get("result_class_encoding") == encoding,
        "v49_ml_fold_identity_is_bound": isinstance(folds, list) and isinstance(train, list) and len(folds) == len(train) and len(set(folds)) >= 2 and all(isinstance(fold, int) and not isinstance(fold, bool) and fold >= 0 for fold in folds) and value.get("result_training_fold_ids") == folds,
        "v49_ml_metric_rows_are_unique_finite_and_bound": isinstance(metric_rows, list) and bool(metric_rows) and len(metric_rows) == len(set(metric_rows)) and _finite_list(metric_values, len(metric_rows)) and value.get("result_metric_row_keys") == metric_rows and value.get("result_metric_values") == metric_values,
        "v49_ml_model_owner_is_bound": str(value.get("model_owner", "")).startswith("model:") and value.get("result_model_owner") == value.get("model_owner"),
        "v49_ml_result_digest_is_bound": _digest_v47(value),
    }


def _agentic_cancel_v49_checks(value: Mapping[str, object]) -> dict[str, bool]:
    timeout = value.get("timeout_s")
    return {
        "v49_agentic_cancel_generation_is_closed": _generation_v47(value, ("timeout_generation", "cancel_generation", "partial_output_generation", "cleanup_generation", "tool_call_generation", "result_generation")),
        "v49_agentic_timeout_is_positive_and_bound": _number(timeout) and float(timeout) > 0.0 and value.get("result_timeout_s") == timeout,
        "v49_agentic_timeout_cancel_is_complete": value.get("timed_out") is value.get("result_timed_out") is True and value.get("cancel_requested") is value.get("result_cancel_requested") is True and value.get("cancel_completed") is value.get("result_cancel_completed") is True,
        "v49_agentic_partial_output_is_discarded": value.get("partial_output_policy") == value.get("result_partial_output_policy") == "discard",
        "v49_agentic_session_is_released": value.get("session_cleanup") == value.get("result_session_cleanup") == "released",
        "v49_agentic_tool_call_owner_is_bound": str(value.get("tool_call_owner", "")).startswith("tool-call:") and value.get("result_tool_call_owner") == value.get("tool_call_owner"),
        "v49_agentic_result_digest_is_bound": _digest_v47(value),
    }


def _parallel_resume_v49_checks(value: Mapping[str, object]) -> dict[str, bool]:
    workers = value.get("worker_ids")
    seeds = value.get("worker_seeds")
    streams = value.get("random_streams")
    count = len(workers) if isinstance(workers, list) else 0
    return {
        "v49_parallel_resume_generation_is_closed": _generation_v47(value, ("worker_generation", "seed_generation", "stream_generation", "trial_generation", "resume_generation", "checkpoint_generation", "experiment_generation", "result_generation")),
        "v49_parallel_workers_are_unique_and_bound": isinstance(workers, list) and bool(workers) and len(workers) == len(set(workers)) and all(isinstance(worker, int) and not isinstance(worker, bool) and worker > 0 for worker in workers) and value.get("result_worker_ids") == workers,
        "v49_parallel_seeds_are_unique_and_bound": isinstance(seeds, list) and len(seeds) == count and len(seeds) == len(set(seeds)) and all(isinstance(seed, int) and not isinstance(seed, bool) and seed >= 0 for seed in seeds) and value.get("result_worker_seeds") == seeds,
        "v49_parallel_streams_are_independent_and_bound": isinstance(streams, list) and len(streams) == count and len(streams) == len(set(streams)) and all(isinstance(stream, str) and stream.startswith("Threefry:") for stream in streams) and value.get("result_random_streams") == streams,
        "v49_parallel_trial_resume_state_is_bound": str(value.get("trial_state", "")).startswith("completed:") and value.get("result_trial_state") == value.get("trial_state") and str(value.get("resume_state", "")).startswith("resumed_from:") and value.get("result_resume_state") == value.get("resume_state"),
        "v49_parallel_checkpoint_experiment_owners_are_bound": str(value.get("checkpoint_owner", "")).startswith("checkpoint:") and value.get("result_checkpoint_owner") == value.get("checkpoint_owner") and str(value.get("experiment_owner", "")).startswith("experiment:") and value.get("result_experiment_owner") == value.get("experiment_owner"),
        "v49_parallel_result_digest_is_bound": _digest_v47(value),
    }


def _generation_v47(value: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(value.get("generation", "")).strip()
    return bool(generation) and all(value.get(field) == generation for field in fields)


def _digest_v47(value: Mapping[str, object]) -> bool:
    return _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256")


def _ml_datastore_v47_checks(value: Mapping[str, object]) -> dict[str, bool]:
    order = value.get("datastore_order")
    labels = value.get("labels")
    partition = value.get("partition")
    train = partition.get("train") if isinstance(partition, Mapping) else None
    validation = partition.get("validation") if isinstance(partition, Mapping) else None
    return {
        "v47_ml_generation_is_closed": _generation_v47(
            value,
            ("datastore_generation", "label_generation", "preprocess_generation", "partition_generation", "owner_generation", "result_generation"),
        ),
        "v47_ml_datastore_labels_are_bound": isinstance(order, list) and bool(order) and len(order) == len(set(order)) and order == value.get("result_datastore_order") and isinstance(labels, list) and len(labels) == len(order) and labels == value.get("result_labels"),
        "v47_ml_preprocess_is_bound": _valid_digest(value.get("preprocess_sha256")) and value.get("result_preprocess_sha256") == value.get("preprocess_sha256"),
        "v47_ml_partition_is_disjoint_and_bound": isinstance(train, list) and bool(train) and isinstance(validation, list) and bool(validation) and set(train).isdisjoint(validation) and value.get("result_partition") == partition,
        "v47_ml_partition_owner_is_bound": str(value.get("partition_owner", "")).startswith("partition:") and value.get("result_partition_owner") == value.get("partition_owner"),
        "v47_ml_result_digest_is_bound": _digest_v47(value),
    }


def _rl_v47_checks(value: Mapping[str, object]) -> dict[str, bool]:
    return {
        "v47_rl_generation_is_closed": _generation_v47(
            value,
            ("observation_generation", "action_generation", "environment_generation", "reset_generation", "policy_generation", "result_generation"),
        ),
        "v47_rl_observation_spec_is_bound": _valid_digest(value.get("observation_spec_sha256")) and value.get("result_observation_spec_sha256") == value.get("observation_spec_sha256"),
        "v47_rl_action_spec_is_bound": _valid_digest(value.get("action_spec_sha256")) and value.get("result_action_spec_sha256") == value.get("action_spec_sha256"),
        "v47_rl_environment_reset_is_bound": str(value.get("environment_id", "")).startswith("environment:") and value.get("result_environment_id") == value.get("environment_id") and value.get("reset_policy") == value.get("result_reset_policy") == "deterministic_seeded",
        "v47_rl_policy_owner_is_bound": str(value.get("policy_owner", "")).startswith("policy:") and value.get("result_policy_owner") == value.get("policy_owner"),
        "v47_rl_result_digest_is_bound": _digest_v47(value),
    }


def _agentic_v47_checks(value: Mapping[str, object]) -> dict[str, bool]:
    return {
        "v47_agentic_generation_is_closed": _generation_v47(
            value,
            ("call_generation", "correlation_generation", "session_generation", "workspace_generation", "release_generation", "result_generation"),
        ),
        "v47_agentic_call_correlation_is_bound": str(value.get("tool_call_id", "")).startswith("call:") and value.get("result_tool_call_id") == value.get("tool_call_id") and str(value.get("correlation_id", "")).startswith("correlation:") and value.get("result_correlation_id") == value.get("correlation_id"),
        "v47_agentic_session_workspace_release_is_bound": str(value.get("session_identity", "")).startswith("matlab:") and value.get("result_session_identity") == value.get("session_identity") and _valid_digest(value.get("workspace_sha256")) and value.get("result_workspace_sha256") == value.get("workspace_sha256") and _same_release(value),
        "v47_agentic_result_digest_is_bound": _digest_v47(value),
    }


def _experiment_v47_checks(value: Mapping[str, object]) -> dict[str, bool]:
    trials = value.get("trial_row_keys")
    folds = value.get("cv_partition_ids")
    indices = value.get("result_indices")
    return {
        "v47_experiment_generation_is_closed": _generation_v47(
            value,
            ("experiment_generation", "trial_generation", "cv_generation", "result_index_generation", "result_generation"),
        ),
        "v47_experiment_trial_rows_are_unique_and_bound": isinstance(trials, list) and bool(trials) and len(trials) == len(set(trials)) and value.get("result_trial_row_keys") == trials,
        "v47_experiment_cv_partitions_are_unique_and_bound": isinstance(folds, list) and len(folds) == len(trials or []) and len(folds) == len(set(folds)) and value.get("result_cv_partition_ids") == folds,
        "v47_experiment_result_indices_are_unique_ordered_and_bound": isinstance(indices, list) and len(indices) == len(trials or []) and all(isinstance(index, int) and index >= 0 for index in indices) and indices == sorted(indices) and len(indices) == len(set(indices)) and value.get("replayed_result_indices") == indices,
        "v47_experiment_owner_is_bound": str(value.get("experiment_owner", "")).startswith("experiment:") and value.get("result_experiment_owner") == value.get("experiment_owner"),
        "v47_experiment_result_digest_is_bound": _digest_v47(value),
    }


def _finite_list(value: object, length: int | None = None) -> bool:
    return isinstance(value, list) and (length is None or len(value) == length) and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)


def _autodiff_v48_checks(value: Mapping[str, object]) -> dict[str, bool]:
    parameters = value.get("parameter_order")
    gradient = value.get("gradient")
    indices = value.get("fd_spotcheck_indices")
    fd_gradient = value.get("fd_gradient")
    spotcheck_ok = (
        isinstance(indices, list)
        and _finite_list(fd_gradient, len(indices))
        and isinstance(gradient, list)
        and all(isinstance(index, int) and 0 <= index < len(gradient) for index in indices)
        and len(indices) == len(set(indices))
        and all(abs(float(gradient[index]) - float(fd)) <= 1.0e-3 * max(1.0, abs(float(fd))) for index, fd in zip(indices, fd_gradient))
    )
    return {
        "v48_autodiff_generation_is_closed": _generation_v47(value, ("parameter_generation", "tape_generation", "objective_generation", "gradient_generation", "fd_generation", "checkpoint_generation", "result_generation")),
        "v48_autodiff_parameters_and_gradient_are_bound": isinstance(parameters, list) and bool(parameters) and len(parameters) == len(set(parameters)) and value.get("result_parameter_order") == parameters and _finite_list(gradient, len(parameters)) and value.get("result_gradient") == gradient,
        "v48_autodiff_tape_and_objective_are_bound": _valid_digest(value.get("gradient_tape_sha256")) and value.get("result_gradient_tape_sha256") == value.get("gradient_tape_sha256") and str(value.get("objective_id", "")).startswith("objective:") and value.get("result_objective_id") == value.get("objective_id"),
        "v48_autodiff_fd_spotcheck_is_bound": spotcheck_ok and value.get("result_fd_spotcheck_indices") == indices and value.get("result_fd_gradient") == fd_gradient,
        "v48_autodiff_checkpoint_owner_is_bound": str(value.get("checkpoint_owner", "")).startswith("checkpoint:") and value.get("result_checkpoint_owner") == value.get("checkpoint_owner"),
        "v48_autodiff_result_digest_is_bound": _digest_v47(value),
    }


def _sequence_v48_checks(value: Mapping[str, object]) -> dict[str, bool]:
    masks = value.get("padding_mask")
    lengths = value.get("sequence_lengths")
    shuffle = value.get("shuffle_order")
    rows = value.get("minibatch_row_keys")
    mask_ok = isinstance(masks, list) and isinstance(lengths, list) and len(masks) == len(lengths) > 0
    if mask_ok:
        width = max(lengths) if all(isinstance(length, int) and length > 0 for length in lengths) else 0
        mask_ok = width > 0 and all(mask == [1] * length + [0] * (width - length) for mask, length in zip(masks, lengths))
    return {
        "v48_sequence_generation_is_closed": _generation_v47(value, ("padding_generation", "mask_generation", "length_generation", "shuffle_generation", "minibatch_generation", "checkpoint_generation", "result_generation")),
        "v48_sequence_padding_mask_lengths_are_bound": value.get("padding_policy") == value.get("result_padding_policy") == "right_zero" and mask_ok and value.get("result_padding_mask") == masks and value.get("result_sequence_lengths") == lengths,
        "v48_sequence_shuffle_minibatch_rows_are_bound": isinstance(shuffle, list) and sorted(shuffle) == list(range(len(lengths or []))) and value.get("result_shuffle_order") == shuffle and isinstance(rows, list) and rows == [f"sequence:{index}" for index in shuffle] and value.get("result_minibatch_row_keys") == rows,
        "v48_sequence_checkpoint_owner_is_bound": str(value.get("checkpoint_owner", "")).startswith("checkpoint:") and value.get("result_checkpoint_owner") == value.get("checkpoint_owner"),
        "v48_sequence_result_digest_is_bound": _digest_v47(value),
    }


def _agentic_workspace_v48_checks(value: Mapping[str, object]) -> dict[str, bool]:
    mutations = value.get("workspace_mutations")
    return {
        "v48_agentic_workspace_generation_is_closed": _generation_v47(value, ("mutation_generation", "diff_generation", "approval_generation", "rollback_generation", "tool_call_generation", "result_generation")),
        "v48_agentic_workspace_mutations_are_bound": isinstance(mutations, list) and bool(mutations) and len(mutations) == len(set(mutations)) and all(isinstance(path, str) and path and ".." not in path for path in mutations) and value.get("result_workspace_mutations") == mutations,
        "v48_agentic_workspace_diff_is_bound": _valid_digest(value.get("workspace_diff_sha256")) and value.get("result_workspace_diff_sha256") == value.get("workspace_diff_sha256"),
        "v48_agentic_approval_and_rollback_are_bound": value.get("approval_scope") == value.get("result_approval_scope") == "workspace_write:approved_paths" and value.get("rollback_state") == value.get("result_rollback_state") == "not_required",
        "v48_agentic_tool_call_owner_is_bound": str(value.get("tool_call_owner", "")).startswith("tool-call:") and value.get("result_tool_call_owner") == value.get("tool_call_owner"),
        "v48_agentic_result_digest_is_bound": _digest_v47(value),
    }


def _experiment_selection_v48_checks(value: Mapping[str, object]) -> dict[str, bool]:
    direction = value.get("metric_direction")
    encoding = value.get("categorical_encoding")
    rows = value.get("trial_row_keys")
    metrics = value.get("metric_values")
    best = value.get("best_trial_index")
    expected = None
    if isinstance(metrics, list) and metrics and _finite_list(metrics):
        expected = min(range(len(metrics)), key=metrics.__getitem__) if direction == "minimize" else max(range(len(metrics)), key=metrics.__getitem__) if direction == "maximize" else None
    encoding_ok = isinstance(encoding, Mapping) and bool(encoding) and all(isinstance(options, list) and bool(options) and len(options) == len(set(options)) for options in encoding.values())
    return {
        "v48_experiment_generation_is_closed": _generation_v47(value, ("metric_generation", "direction_generation", "encoding_generation", "trial_generation", "selection_generation", "result_generation")),
        "v48_experiment_metric_direction_is_bound": bool(str(value.get("metric_name", ""))) and value.get("result_metric_name") == value.get("metric_name") and direction in {"minimize", "maximize"} and value.get("result_metric_direction") == direction,
        "v48_experiment_categorical_encoding_is_bound": encoding_ok and value.get("result_categorical_encoding") == encoding,
        "v48_experiment_trials_and_metrics_are_bound": isinstance(rows, list) and bool(rows) and len(rows) == len(set(rows)) and value.get("result_trial_row_keys") == rows and _finite_list(metrics, len(rows)) and value.get("result_metric_values") == metrics,
        "v48_experiment_best_trial_is_bound": isinstance(best, int) and best == expected and value.get("result_best_trial_index") == best and value.get("best_result_row") == rows[best] and value.get("result_best_result_row") == rows[best],
        "v48_experiment_owner_is_bound": str(value.get("experiment_owner", "")).startswith("experiment:") and value.get("result_experiment_owner") == value.get("experiment_owner"),
        "v48_experiment_result_digest_is_bound": _digest_v47(value),
    }


def _v46_release_and_digest(value: Mapping[str, object]) -> dict[str, bool]:
    return {
        "v46_release_is_bound": _same_release(value),
        "v46_digest_is_bound": _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256"),
    }


def _supervised_v46_checks(value: Mapping[str, object]) -> dict[str, bool]:
    checks = {
        "v46_supervised_nonfinite_policy_is_bound": value.get("nonfinite_policy") == value.get("result_nonfinite_policy") == "drop_with_count" and value.get("nonfinite_input_count") == value.get("result_nonfinite_input_count") == 0,
        "v46_supervised_split_and_normalization_are_bound": bool(str(value.get("split_generation", "")).strip()) and value.get("result_split_generation") == value.get("split_generation") and value.get("normalization_fit_scope") == value.get("result_normalization_fit_scope") == "training_only",
        "v46_supervised_worker_seed_is_bound": isinstance(value.get("worker_seed"), int) and value.get("worker_seed") >= 0 and value.get("result_worker_seed") == value.get("worker_seed"),
        "v46_supervised_restart_state_is_fresh": value.get("restart_state") == value.get("result_restart_state") == "fresh_training",
    }
    checks.update(_v46_release_and_digest(value))
    return checks


def _rl_v46_checks(value: Mapping[str, object]) -> dict[str, bool]:
    checks = {
        "v46_rl_episode_timeout_is_bound": isinstance(value.get("episode_timeout_steps"), int) and value.get("episode_timeout_steps") > 0 and value.get("result_episode_timeout_steps") == value.get("episode_timeout_steps"),
        "v46_rl_termination_is_environment_defined": value.get("termination_semantics") == value.get("result_termination_semantics") == "environment_defined",
        "v46_rl_exploration_is_training_only": value.get("exploration_mode") == value.get("result_exploration_mode") == "training_only" and value.get("evaluation_mode") == value.get("result_evaluation_mode") == "greedy_no_exploration",
        "v46_rl_checkpoint_is_bound": bool(str(value.get("checkpoint_generation", "")).strip()) and value.get("result_checkpoint_generation") == value.get("checkpoint_generation"),
    }
    checks.update(_v46_release_and_digest(value))
    return checks


def _agentic_v46_checks(value: Mapping[str, object]) -> dict[str, bool]:
    checks = {
        "v46_agentic_argument_shape_is_json_object": value.get("argument_schema") == value.get("result_argument_schema") == "json_object" and value.get("argument_shape_valid") == value.get("result_argument_shape_valid") is True,
        "v46_agentic_existing_session_and_timeout_are_bound": value.get("session_detection") == value.get("result_session_detection") == "existing_shared_matlab" and isinstance(value.get("timeout_s"), (int, float)) and float(value.get("timeout_s")) > 0.0 and value.get("result_timeout_s") == value.get("timeout_s"),
        "v46_agentic_error_class_is_bound": value.get("error_class") == value.get("result_error_class") == "none",
        "v46_agentic_tool_arguments_are_bound": _valid_digest(value.get("tool_arguments_sha256")) and value.get("result_tool_arguments_sha256") == value.get("tool_arguments_sha256") and _same_release(value) and str(value.get("owner", "")).startswith("matlab:") and value.get("result_owner") == value.get("owner"),
    }
    return checks


def _checkpoint_v46_checks(value: Mapping[str, object]) -> dict[str, bool]:
    order = value.get("checkpoint_order")
    checks = {
        "v46_checkpoint_worker_seed_is_bound": isinstance(value.get("worker_seed"), int) and value.get("worker_seed") >= 0 and value.get("result_worker_seed") == value.get("worker_seed"),
        "v46_checkpoint_order_is_monotone": isinstance(order, list) and order == value.get("result_checkpoint_order") and order == sorted(order),
        "v46_checkpoint_optimizer_state_is_bound": value.get("optimizer_state") == value.get("result_optimizer_state") == "adam_ready",
        "v46_checkpoint_generation_and_owner_are_bound": bool(str(value.get("checkpoint_generation", "")).strip()) and value.get("result_checkpoint_generation") == value.get("checkpoint_generation") and _same_release(value) and str(value.get("owner", "")).startswith("matlab:") and value.get("result_owner") == value.get("owner"),
    }
    checks.update({"v46_checkpoint_digest_is_bound": _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256")})
    return checks


def _same_release(value: Mapping[str, object]) -> bool:
    release = str(value.get("release_id", "")).strip()
    return bool(release) and value.get("result_release_id") == release


def _supervised_v45_checks(value: Mapping[str, object]) -> dict[str, bool]:
    return {
        "v45_supervised_schema_is_pinned": value.get("schema") == "cae-ai-lab.matlab-ml-rl-result.v3",
        "v45_datastore_is_pinned": bool(str(value.get("datastore_id", "")).strip()),
        "v45_split_generation_is_bound": bool(str(value.get("split_generation", "")).strip()) and value.get("result_split_generation") == value.get("split_generation"),
        "v45_normalization_is_fit_on_training_only": value.get("normalization_fit_scope") == "training_only" and value.get("result_normalization_fit_scope") == "training_only",
        "v45_hyperparameter_selection_uses_training_cv": value.get("hyperparameter_selection_source") == "training_cross_validation" and value.get("result_hyperparameter_selection_source") == "training_cross_validation",
        "v45_holdout_is_separate": bool(str(value.get("holdout_id", "")).strip()) and value.get("holdout_id") not in set(_string_list(value.get("training_ids"))),
        "v45_model_card_matches_release": value.get("model_card_release") == value.get("release_id"),
        "v45_release_is_bound": _same_release(value),
        "v45_digest_is_bound": _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256"),
    }


def _rl_v45_checks(value: Mapping[str, object]) -> dict[str, bool]:
    discount = value.get("discount_factor")
    return {
        "v45_rl_replay_buffer_is_bound": bool(str(value.get("replay_buffer_generation", "")).strip()) and value.get("result_replay_buffer_generation") == value.get("replay_buffer_generation"),
        "v45_rl_termination_semantics_are_environment_defined": value.get("termination_semantics") == "environment_defined" and value.get("result_termination_semantics") == "environment_defined",
        "v45_rl_discount_is_valid_and_bound": _number(discount) and 0.0 < float(discount) <= 1.0 and value.get("result_discount_factor") == discount,
        "v45_rl_policy_seed_is_recorded": isinstance(value.get("policy_seed"), int) and value.get("policy_seed") >= 0,
        "v45_rl_evaluation_is_separate": value.get("evaluation_mode") == "greedy_no_exploration" and value.get("result_evaluation_mode") == "greedy_no_exploration",
        "v45_rl_release_is_bound": _same_release(value),
        "v45_rl_digest_is_bound": _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256"),
    }


def _agentic_v45_checks(value: Mapping[str, object]) -> dict[str, bool]:
    return {
        "v45_agentic_capability_route_is_pinned": value.get("capability_route") == "matlab_tool" and value.get("result_capability_route") == "matlab_tool",
        "v45_agentic_consent_is_recorded": value.get("consent_recorded") is True and value.get("result_consent_recorded") is True,
        "v45_agentic_tool_arguments_digest_is_bound": _valid_digest(value.get("tool_arguments_sha256")) and value.get("result_tool_arguments_sha256") == value.get("tool_arguments_sha256"),
        "v45_agentic_existing_session_is_explicit": value.get("session_detection") == "existing_shared_matlab" and value.get("result_session_detection") == "existing_shared_matlab",
        "v45_agentic_release_and_owner_are_bound": _same_release(value) and str(value.get("owner", "")).startswith("matlab:") and value.get("result_owner") == value.get("owner"),
    }


def _checkpoint_v45_checks(value: Mapping[str, object]) -> dict[str, bool]:
    return {
        "v45_checkpoint_generation_is_bound": bool(str(value.get("checkpoint_generation", "")).strip()) and value.get("result_checkpoint_generation") == value.get("checkpoint_generation"),
        "v45_datastore_state_is_replayed": value.get("datastore_state") == "replayed" and value.get("result_datastore_state") == "replayed",
        "v45_optimizer_is_recorded": bool(str(value.get("optimizer", "")).strip()) and value.get("result_optimizer") == value.get("optimizer"),
        "v45_discount_is_bound": _number(value.get("discount_factor")) and value.get("result_discount_factor") == value.get("discount_factor"),
        "v45_evaluation_id_is_separate": bool(str(value.get("evaluation_id", "")).strip()) and value.get("evaluation_id") != value.get("checkpoint_generation"),
        "v45_checkpoint_release_and_owner_are_bound": _same_release(value) and str(value.get("owner", "")).startswith("matlab:") and value.get("result_owner") == value.get("owner"),
        "v45_checkpoint_digest_is_bound": _valid_digest(value.get("result_sha256")) and value.get("accepted_result_sha256") == value.get("result_sha256"),
    }


def _supervised_v44_checks(value: Mapping[str, object]) -> dict[str, bool]:
    training_ids = _string_list(value.get("training_ids"))
    evaluation_ids = _string_list(value.get("evaluation_ids"))
    validation_metric = value.get("validation_metric")
    metric_name = validation_metric.get("name") if isinstance(validation_metric, Mapping) else None
    metric_value = validation_metric.get("value") if isinstance(validation_metric, Mapping) else None
    return {
        "supervised_schema_is_pinned": value.get("schema") == "cae-ai-lab.matlab-ml-rl-result.v2",
        "supervised_task_is_regression_or_classification": value.get("task") in {"regression", "classification"},
        "supervised_release_is_pinned": bool(str(value.get("matlab_release", "")).strip()),
        "supervised_session_owner_is_recorded": str(value.get("session_owner", "")).startswith("matlab:"),
        "supervised_seed_is_recorded": isinstance(value.get("random_seed"), int) and value.get("random_seed") >= 0,
        "supervised_split_is_named": bool(str(value.get("split_id", "")).strip()),
        "supervised_training_ids_are_nonempty": bool(training_ids),
        "supervised_evaluation_ids_are_nonempty": bool(evaluation_ids),
        "supervised_training_and_evaluation_ids_are_disjoint": bool(training_ids) and bool(evaluation_ids) and set(training_ids).isdisjoint(evaluation_ids),
        "supervised_training_and_evaluation_are_disjoint": value.get("training_evaluation_disjoint") is True,
        "supervised_validation_metric_is_finite": _number(metric_value),
        "supervised_validation_metric_is_not_training_score": bool(str(metric_name).strip()) and not str(metric_name).lower().startswith("training"),
        "supervised_result_digest_is_valid": _valid_digest(value.get("result_sha256")),
        "supervised_accepted_digest_matches": value.get("accepted_result_sha256") == value.get("result_sha256"),
        "supervised_timing_breakdown_is_bounded": _valid_timings(value.get("timing_breakdown_s")),
    }


def _rl_v44_checks(value: Mapping[str, object]) -> dict[str, bool]:
    evaluation_return = value.get("evaluation_mean_return")
    evaluation_std = value.get("evaluation_std_return")
    return {
        "rl_schema_is_pinned": value.get("schema") == "cae-ai-lab.matlab-ml-rl-result.v2",
        "rl_task_is_supported": value.get("task") == "reinforcement_learning",
        "rl_release_is_pinned": bool(str(value.get("matlab_release", "")).strip()),
        "rl_session_owner_is_recorded": str(value.get("session_owner", "")).startswith("matlab:"),
        "rl_seed_is_recorded": isinstance(value.get("random_seed"), int) and value.get("random_seed") >= 0,
        "rl_environment_is_pinned": bool(str(value.get("environment_id", "")).strip()),
        "rl_evaluation_environment_matches": value.get("evaluation_environment_id") == value.get("environment_id") and bool(str(value.get("evaluation_environment_id", "")).strip()),
        "rl_training_episodes_are_positive": isinstance(value.get("training_episodes"), int) and value.get("training_episodes") > 0,
        "rl_evaluation_episodes_are_positive": isinstance(value.get("evaluation_episodes"), int) and value.get("evaluation_episodes") > 0,
        "rl_evaluation_seed_is_fresh": isinstance(value.get("evaluation_seed"), int) and value.get("evaluation_seed") >= 0 and value.get("evaluation_seed") != value.get("random_seed"),
        "rl_evaluation_has_no_exploration": value.get("exploration_during_evaluation") is False,
        "rl_training_and_evaluation_are_disjoint": value.get("training_evaluation_disjoint") is True,
        "rl_evaluation_return_is_finite": _number(evaluation_return),
        "rl_evaluation_std_is_finite_nonnegative": _number(evaluation_std) and float(evaluation_std) >= 0.0,
        "rl_result_digest_is_valid": _valid_digest(value.get("result_sha256")),
        "rl_accepted_digest_matches": value.get("accepted_result_sha256") == value.get("result_sha256"),
        "rl_evaluation_digest_is_valid": _valid_digest(value.get("evaluation_result_sha256")),
        "rl_accepted_evaluation_digest_matches": value.get("accepted_evaluation_result_sha256") == value.get("evaluation_result_sha256"),
        "rl_timing_breakdown_is_bounded": _valid_timings(value.get("timing_breakdown_s")),
    }


def _string_list(value: object) -> list[str]:
    if not isinstance(value, (list, tuple)):
        return []
    return [item for item in value if isinstance(item, str) and item.strip()]


def _valid_digest(value: object) -> bool:
    return isinstance(value, str) and bool(_DIGEST.fullmatch(value))


def _valid_timings(value: object) -> bool:
    if not isinstance(value, Mapping) or not 1 <= len(value) <= 4:
        return False
    return all(_number(item) and float(item) >= 0.0 for item in value.values())


def _number(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False
