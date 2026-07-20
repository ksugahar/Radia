"""Identity gates for n-step returns and cross-validation evidence."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence


_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_TOL = 1.0e-10


def validate_matlab_ml_rl_v55_identity(summary: Mapping[str, object]) -> dict[str, object] | None:
    """Validate v55 n-step-return and cross-validation identities."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity = summary.get("matlab_ml_rl_v54_identity")
    if not isinstance(identity, Mapping):
        return None
    nstep = identity.get("nstep_return")
    cross_validation = identity.get("cross_validation")
    if nstep is None and cross_validation is None:
        return None
    checks = {
        "v55_nstep_return_record_is_mapping": isinstance(nstep, Mapping),
        "v55_cross_validation_record_is_mapping": isinstance(cross_validation, Mapping),
    }
    if isinstance(nstep, Mapping):
        checks.update(_nstep_checks(nstep))
    if isinstance(cross_validation, Mapping):
        checks.update(_cross_validation_checks(cross_validation))
    return {
        "policy": "matlab_ml_rl_artifact_gate_v55",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
    }


def _nstep_checks(value: Mapping[str, object]) -> dict[str, bool]:
    gamma = value.get("gamma")
    rewards = _float_vector(value.get("rewards"))
    terminal = value.get("terminal")
    bootstrap = value.get("bootstrap_value")
    reported_return = value.get("n_step_return")
    parameters_ok = (
        _finite(gamma)
        and 0.0 <= float(gamma) <= 1.0
        and bool(rewards)
        and isinstance(terminal, bool)
        and _finite(bootstrap)
        and _finite(reported_return)
    )
    return_ok = False
    if parameters_ok:
        expected = sum(float(gamma) ** index * reward for index, reward in enumerate(rewards))
        if not terminal:
            expected += float(gamma) ** len(rewards) * float(bootstrap)
        return_ok = _close(float(reported_return), expected)
    return {
        "v55_nstep_generation_is_closed": _generation_closed(value, ("reward_generation", "gamma_generation", "terminal_generation", "bootstrap_generation", "return_generation", "trajectory_generation", "owner_generation", "result_generation")),
        "v55_nstep_parameters_are_valid_and_bound": parameters_ok and value.get("result_gamma") == gamma and value.get("result_rewards") == value.get("rewards") and value.get("result_terminal") is terminal and value.get("result_bootstrap_value") == bootstrap,
        "v55_nstep_return_is_recomputed_and_bound": return_ok and value.get("result_n_step_return") == reported_return,
        "v55_nstep_trajectory_is_bound": _prefixed_equal(value, "trajectory_id", "result_trajectory_id", "trajectory:"),
        "v55_nstep_policy_owner_is_bound": _prefixed_equal(value, "policy_owner", "result_policy_owner", "policy:"),
        "v55_nstep_result_digest_is_bound": _result_digest_bound(value),
    }


def _cross_validation_checks(value: Mapping[str, object]) -> dict[str, bool]:
    fold_ids = _integer_vector(value.get("fold_id_per_sample"))
    labels = value.get("class_labels")
    folds_ok = bool(fold_ids) and min(fold_ids) == 0 and sorted(set(fold_ids)) == list(range(max(fold_ids) + 1)) and len(set(fold_ids)) >= 2
    labels_ok = (
        isinstance(labels, Sequence)
        and not isinstance(labels, (str, bytes))
        and len(labels) == len(fold_ids)
        and all(isinstance(label, str) and label.startswith("class:") and len(label) > len("class:") for label in labels)
    )
    stratified_ok = folds_ok and labels_ok
    if stratified_ok:
        classes = set(labels)
        stratified_ok = all({labels[index] for index, fold in enumerate(fold_ids) if fold == fold_id} == classes for fold_id in set(fold_ids))
    fit_rows = value.get("preprocess_fit_rows")
    fit_scope_ok = folds_ok and isinstance(fit_rows, Mapping) and set(fit_rows) == {str(fold_id) for fold_id in set(fold_ids)}
    if fit_scope_ok:
        all_rows = set(range(len(fold_ids)))
        for fold_id in set(fold_ids):
            rows = _integer_vector(fit_rows[str(fold_id)])
            expected = all_rows - {index for index, fold in enumerate(fold_ids) if fold == fold_id}
            if len(rows) != len(set(rows)) or set(rows) != expected:
                fit_scope_ok = False
                break
    seed = value.get("rng_seed")
    metrics = _float_vector(value.get("fold_metrics"))
    aggregate = value.get("aggregate_metric")
    metrics_ok = folds_ok and len(metrics) == len(set(fold_ids)) and _finite(aggregate) and _close(float(aggregate), sum(metrics) / len(metrics))
    return {
        "v55_cross_validation_generation_is_closed": _generation_closed(value, ("fold_generation", "stratification_generation", "preprocess_generation", "seed_generation", "metric_generation", "owner_generation", "result_generation")),
        "v55_cross_validation_folds_are_complete_and_bound": folds_ok and value.get("result_fold_id_per_sample") == value.get("fold_id_per_sample"),
        "v55_cross_validation_is_stratified_and_bound": stratified_ok and value.get("result_class_labels") == labels,
        "v55_cross_validation_preprocess_scope_is_train_only_and_bound": fit_scope_ok and value.get("result_preprocess_fit_rows") == fit_rows,
        "v55_cross_validation_seed_is_nonnegative_and_bound": isinstance(seed, int) and not isinstance(seed, bool) and seed >= 0 and value.get("result_rng_seed") == seed,
        "v55_cross_validation_metrics_are_aggregated_and_bound": metrics_ok and value.get("result_fold_metrics") == value.get("fold_metrics") and value.get("result_aggregate_metric") == aggregate,
        "v55_cross_validation_model_owner_is_bound": _prefixed_equal(value, "model_owner", "result_model_owner", "model:"),
        "v55_cross_validation_result_digest_is_bound": _result_digest_bound(value),
    }


def _float_vector(value: object) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value or not all(_finite(item) for item in value):
        return []
    return [float(item) for item in value]


def _integer_vector(value: object) -> list[int]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value or not all(isinstance(item, int) and not isinstance(item, bool) and item >= 0 for item in value):
        return []
    return list(value)


def _generation_closed(value: Mapping[str, object], fields: Sequence[str]) -> bool:
    generation = str(value.get("generation", "")).strip()
    return bool(generation) and all(value.get(field) == generation for field in fields)


def _prefixed_equal(value: Mapping[str, object], left: str, right: str, prefix: str) -> bool:
    return str(value.get(left, "")).startswith(prefix) and value.get(left) == value.get(right)


def _result_digest_bound(value: Mapping[str, object]) -> bool:
    digest = value.get("result_sha256")
    return (
        isinstance(digest, str)
        and bool(_DIGEST.fullmatch(digest))
        and value.get("accepted_result_sha256") == digest
    )


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=_TOL, abs_tol=_TOL)
