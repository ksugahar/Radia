from __future__ import annotations

from copy import deepcopy

from radia_mcp.matlab_agentic_ml import validate_matlab_ml_rl_v50_identity


PROMOTED_CASE_IDS = {
    "v50_public_bayesopt_objective_noise_constraint_seed_acquisition_model_owner_mismatch",
    "v50_public_deep_learning_mixed_precision_loss_scale_optimizer_state_checkpoint_owner_mismatch",
}


def _summary() -> dict[str, object]:
    bayes = "bayes-v50"
    mixed = "mixed-v50"
    sandbox = "sandbox-v50"
    tall = "tall-v50"
    constraints = [{"name": "temperature_c", "sense": "<=", "limit": 120.0}]
    acquisition = {
        "name": "expected-improvement-plus",
        "iteration": 18,
        "incumbent": "trial:12",
    }
    optimizer = {
        "name": "adam",
        "iteration": 240,
        "learn_rate": 0.001,
        "moment1_sha256": "1" * 64,
        "moment2_sha256": "2" * 64,
    }
    roots = ["workspace:project", "workspace:artifacts"]
    partitions = ["partition:0", "partition:1", "partition:2", "partition:3"]
    workers = [1, 2, 1, 2]
    return {"matlab_ml_rl_v50_identity": {
        "bayesopt": {
            "generation": bayes,
            **{key: bayes for key in ("objective_generation", "noise_generation", "constraint_generation", "seed_generation", "acquisition_generation", "model_generation", "result_generation")},
            "objective_name": "validation_loss", "result_objective_name": "validation_loss",
            "objective_noise_sigma": 0.02, "result_objective_noise_sigma": 0.02,
            "constraints": constraints, "result_constraints": constraints,
            "rng_seed": 50123, "result_rng_seed": 50123,
            "acquisition_state": acquisition, "result_acquisition_state": acquisition,
            "model_owner": "model:bayes-v50", "result_model_owner": "model:bayes-v50",
            "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
        },
        "mixed_precision": {
            "generation": mixed,
            **{key: mixed for key in ("precision_generation", "loss_scale_generation", "optimizer_generation", "checkpoint_generation", "network_generation", "result_generation")},
            "precision_policy": "mixed-fp16", "result_precision_policy": "mixed-fp16",
            "loss_scale": 1024.0, "result_loss_scale": 1024.0,
            "optimizer_state": optimizer, "result_optimizer_state": optimizer,
            "checkpoint_sha256": "4" * 64, "result_checkpoint_sha256": "4" * 64,
            "network_owner": "network:mixed-v50", "result_network_owner": "network:mixed-v50",
            "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
        },
        "agentic_sandbox": {
            "generation": sandbox,
            **{key: sandbox for key in ("root_generation", "resolution_generation", "traversal_generation", "approval_generation", "tool_call_generation", "result_generation")},
            "allowed_roots": roots, "result_allowed_roots": roots,
            "requested_path": "workspace:project/src/train.m", "result_requested_path": "workspace:project/src/train.m",
            "resolved_path": "workspace:project/src/train.m", "result_resolved_path": "workspace:project/src/train.m",
            "symlink_resolution": "inside_allowed_root", "result_symlink_resolution": "inside_allowed_root",
            "traversal_detected": False, "result_traversal_detected": False,
            "approval_id": "approval:v50", "result_approval_id": "approval:v50",
            "tool_call_owner": "tool-call:v50", "result_tool_call_owner": "tool-call:v50",
            "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64,
        },
        "tall_datastore": {
            "generation": tall,
            **{key: tall for key in ("partition_generation", "worker_generation", "checkpoint_generation", "resume_generation", "datastore_generation", "result_generation")},
            "partition_ids": partitions, "result_partition_ids": partitions,
            "worker_order": workers, "result_worker_order": workers,
            "checkpoint_sha256": "7" * 64, "result_checkpoint_sha256": "7" * 64,
            "resume_cursor": "partition:2/row:4096", "result_resume_cursor": "partition:2/row:4096",
            "datastore_owner": "datastore:v50", "result_datastore_owner": "datastore:v50",
            "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
        },
    }}


def test_v50_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v50_identity(_summary())
    assert result and result["status"] == "ok"


def test_v50_bayesopt_mutation_is_rejected() -> None:
    summary = deepcopy(_summary())
    row = summary["matlab_ml_rl_v50_identity"]["bayesopt"]
    row.update({"result_objective_name": "training_loss", "result_constraints": [], "result_rng_seed": 999, "result_model_owner": "model:foreign"})
    result = validate_matlab_ml_rl_v50_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v50_mixed_precision_mutation_is_rejected() -> None:
    summary = deepcopy(_summary())
    row = summary["matlab_ml_rl_v50_identity"]["mixed_precision"]
    row.update({"result_precision_policy": "fp32", "result_loss_scale": 1.0, "result_checkpoint_sha256": "9" * 64, "result_network_owner": "network:foreign"})
    result = validate_matlab_ml_rl_v50_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v50_self_consistent_unsafe_states_are_rejected() -> None:
    summary = deepcopy(_summary())
    mixed = summary["matlab_ml_rl_v50_identity"]["mixed_precision"]
    mixed["precision_policy"] = mixed["result_precision_policy"] = "fp32"
    sandbox = summary["matlab_ml_rl_v50_identity"]["agentic_sandbox"]
    sandbox["requested_path"] = sandbox["result_requested_path"] = "workspace:project/../private/key.txt"
    sandbox["traversal_detected"] = sandbox["result_traversal_detected"] = True
    result = validate_matlab_ml_rl_v50_identity(summary)
    assert result and result["status"] == "needs_attention"
