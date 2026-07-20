from copy import deepcopy

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v55_identity


PROMOTED_CASE_IDS = {
    "v55_public_rl_nstepreturn_gamma_terminal_reward_bootstrap_trajectory_owner_mismatch",
    "v55_public_ml_crossvalidation_fold_stratification_preprocess_seed_metric_owner_mismatch",
}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _summary() -> dict[str, object]:
    gamma = 0.9
    rewards = [1.0, 0.5, -0.2]
    nstep_return = sum(gamma**index * reward for index, reward in enumerate(rewards)) + gamma**len(rewards) * 2.0
    folds = [0, 0, 1, 1, 2, 2]
    labels = ["class:A", "class:B"] * 3
    fit_rows = {"0": [2, 3, 4, 5], "1": [0, 1, 4, 5], "2": [0, 1, 2, 3]}
    metrics = [0.8, 0.7, 0.9]
    return {
        "matlab_ml_rl_v54_identity": {
            "nstep_return": {
                **_generations("nstep-v55-test", ("reward_generation", "gamma_generation", "terminal_generation", "bootstrap_generation", "return_generation", "trajectory_generation", "owner_generation", "result_generation")),
                "gamma": gamma, "result_gamma": gamma,
                "rewards": rewards, "result_rewards": rewards,
                "terminal": False, "result_terminal": False,
                "bootstrap_value": 2.0, "result_bootstrap_value": 2.0,
                "n_step_return": nstep_return, "result_n_step_return": nstep_return,
                "trajectory_id": "trajectory:v55", "result_trajectory_id": "trajectory:v55",
                "policy_owner": "policy:nstep-v55", "result_policy_owner": "policy:nstep-v55",
                "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
            },
            "cross_validation": {
                **_generations("crossval-v55-test", ("fold_generation", "stratification_generation", "preprocess_generation", "seed_generation", "metric_generation", "owner_generation", "result_generation")),
                "fold_id_per_sample": folds, "result_fold_id_per_sample": folds,
                "class_labels": labels, "result_class_labels": labels,
                "preprocess_fit_rows": fit_rows, "result_preprocess_fit_rows": fit_rows,
                "rng_seed": 1729, "result_rng_seed": 1729,
                "fold_metrics": metrics, "result_fold_metrics": metrics,
                "aggregate_metric": sum(metrics) / len(metrics), "result_aggregate_metric": sum(metrics) / len(metrics),
                "model_owner": "model:crossval-v55", "result_model_owner": "model:crossval-v55",
                "result_sha256": "9" * 64, "accepted_result_sha256": "9" * 64,
            },
        }
    }


def test_v55_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v55_identity(_summary())
    assert result and result["status"] == "ok"


def test_v55_frozen_public_counterfactuals_are_rejected() -> None:
    summary = deepcopy(_summary())
    nstep = summary["matlab_ml_rl_v54_identity"]["nstep_return"]
    nstep.update({"result_gamma": 0.5, "result_rewards": [99.0], "result_terminal": True, "result_bootstrap_value": 0.0, "result_n_step_return": -1.0, "result_trajectory_id": "trajectory:stale", "result_policy_owner": "policy:stale"})
    cross_validation = summary["matlab_ml_rl_v54_identity"]["cross_validation"]
    cross_validation.update({"result_fold_id_per_sample": [0] * 6, "result_class_labels": ["class:A"] * 6, "result_preprocess_fit_rows": {"0": [0, 1]}, "result_rng_seed": 42, "result_fold_metrics": [1.0], "result_aggregate_metric": 1.0, "result_model_owner": "model:stale"})
    result = validate_matlab_ml_rl_v55_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v55_self_consistent_invalid_return_and_terminal_semantics_are_rejected() -> None:
    summary = deepcopy(_summary())
    nstep = summary["matlab_ml_rl_v54_identity"]["nstep_return"]
    nstep["terminal"] = nstep["result_terminal"] = True
    result = validate_matlab_ml_rl_v55_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v55_self_consistent_cv_leakage_and_unstratified_folds_are_rejected() -> None:
    summary = deepcopy(_summary())
    cross_validation = summary["matlab_ml_rl_v54_identity"]["cross_validation"]
    labels = ["class:A", "class:A", "class:B", "class:B", "class:A", "class:B"]
    fit_rows = {"0": [0, 1, 2, 3, 4, 5], "1": [0, 1, 4, 5], "2": [0, 1, 2, 3]}
    cross_validation["class_labels"] = cross_validation["result_class_labels"] = labels
    cross_validation["preprocess_fit_rows"] = cross_validation["result_preprocess_fit_rows"] = fit_rows
    result = validate_matlab_ml_rl_v55_identity(summary)
    assert result and result["status"] == "needs_attention"
