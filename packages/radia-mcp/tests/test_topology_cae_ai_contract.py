import json

import pytest

from radia_mcp.topology_optimization.cae_ai_contract import evaluate_cae_ai_artifact
from radia_mcp.topology_optimization.server import topology_opt_cae_ai_artifact_gate


COMMON = {
    "problem_id": "motor-slot-shape-001",
    "seed": 17,
    "input_schema": "design-parameters/v1",
    "output_schema": "candidate-geometry/v1",
    "training_data_provenance": "public-safe generated training manifest sha256:abc",
    "model_or_algorithm_version": "model-0.3.0",
    "run_date_utc": "2026-07-11T00:00:00Z",
    "metrics": {
        "objective": {"value": 0.82, "unit": "dimensionless"},
        "volume": {"value": 1.4e-5, "unit": "m^3"},
    },
    "forward_solver_verification": {
        "solver": "radia-ngsolve",
        "case_id": "forward-001",
        "observables": ["torque_Nm", "volume_m3"],
        "tolerances": {"torque_relative": 0.02, "volume_relative": 1e-9},
        "status": "verified",
    },
}


def test_diffusion_artifact_requires_forward_physics_closure():
    artifact = {
        **COMMON,
        "noise_schedule": "cosine",
        "sampling_steps": 50,
        "conditioning": {"target_torque_Nm": 8.0},
        "generated_candidate_count": 16,
    }
    result = evaluate_cae_ai_artifact("diffusion", artifact)
    assert result["status"] == "ok"
    assert all(result["checks"].values())

    wrapped = json.loads(topology_opt_cae_ai_artifact_gate("diffusion", json.dumps(artifact)))
    assert wrapped["status"] == "ok"


@pytest.mark.parametrize(
    ("family", "extra"),
    [
        ("normalizing_flow", {"jacobian_logdet_check": "passed", "likelihood_metric": "nll", "invertibility_check": "passed"}),
        ("reinforcement_learning", {"state_schema": "state/v1", "action_schema": "action/v1", "reward_definition": "-loss", "termination_rule": "100 steps", "evaluation_episodes": 20}),
        ("pseudoinverse", {"matrix_shape": [20, 8], "rank": 8, "singular_values": [4.0, 2.0, 1.0], "regularization": "tsvd", "residual_norm": 0.01, "solution_norm": 1.2}),
    ],
)
def test_other_cae_ai_families_have_explicit_contracts(family, extra):
    assert evaluate_cae_ai_artifact(family, {**COMMON, **extra})["status"] == "ok"


def test_artifact_without_forward_verification_is_rejected():
    artifact = {
        **COMMON,
        "noise_schedule": "cosine",
        "sampling_steps": 50,
        "conditioning": {"target_torque_Nm": 8.0},
        "generated_candidate_count": 16,
        "forward_solver_verification": {"status": "candidate"},
    }
    result = evaluate_cae_ai_artifact("diffusion", artifact)
    assert result["status"] == "needs_attention"
    assert result["checks"]["forward_solver_verified"] is False


def test_metrics_without_units_are_rejected():
    artifact = {
        **COMMON,
        "noise_schedule": "cosine",
        "sampling_steps": 50,
        "conditioning": {"target_torque_Nm": 8.0},
        "generated_candidate_count": 16,
        "metrics": {"objective": {"value": 0.82}},
    }
    result = evaluate_cae_ai_artifact("diffusion", artifact)
    assert result["status"] == "needs_attention"
    assert result["checks"]["metrics_have_values_and_units"] is False


def test_invalid_family_and_json_are_rejected_at_boundary():
    with pytest.raises(ValueError, match="method_family"):
        evaluate_cae_ai_artifact("gan", COMMON)
    with pytest.raises(ValueError, match="valid JSON"):
        topology_opt_cae_ai_artifact_gate("diffusion", "not-json")
