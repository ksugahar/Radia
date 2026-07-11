"""Artifact contract for reproducible CAE-AI optimization and inverse design."""
from __future__ import annotations

import json


FAMILY_REQUIREMENTS = {
    "diffusion": (
        "noise_schedule",
        "sampling_steps",
        "conditioning",
        "generated_candidate_count",
    ),
    "normalizing_flow": (
        "jacobian_logdet_check",
        "likelihood_metric",
        "invertibility_check",
    ),
    "reinforcement_learning": (
        "state_schema",
        "action_schema",
        "reward_definition",
        "termination_rule",
        "evaluation_episodes",
    ),
    "pseudoinverse": (
        "matrix_shape",
        "rank",
        "singular_values",
        "regularization",
        "residual_norm",
        "solution_norm",
    ),
}

COMMON_REQUIREMENTS = (
    "problem_id",
    "seed",
    "input_schema",
    "output_schema",
    "training_data_provenance",
    "model_or_algorithm_version",
    "run_date_utc",
    "metrics",
    "forward_solver_verification",
)


def _present(value: object) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (list, dict, tuple)):
        return bool(value)
    return True


def evaluate_cae_ai_artifact(method_family: str, artifact: dict) -> dict:
    """Validate reproducibility and forward-physics closure for a CAE-AI run."""
    family = method_family.strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "flow": "normalizing_flow",
        "normalising_flow": "normalizing_flow",
        "rl": "reinforcement_learning",
        "pseudo_inverse": "pseudoinverse",
        "pinv": "pseudoinverse",
    }
    family = aliases.get(family, family)
    if family not in FAMILY_REQUIREMENTS:
        raise ValueError(
            "method_family must be diffusion, normalizing_flow, "
            "reinforcement_learning, or pseudoinverse"
        )

    required = COMMON_REQUIREMENTS + FAMILY_REQUIREMENTS[family]
    missing = [name for name in required if not _present(artifact.get(name))]
    errors = [f"missing required field: {name}" for name in missing]

    verification = artifact.get("forward_solver_verification")
    if not isinstance(verification, dict):
        errors.append("forward_solver_verification must be an object")
        verification = {}
    verification_required = ("solver", "case_id", "observables", "tolerances", "status")
    verification_missing = [
        name for name in verification_required if not _present(verification.get(name))
    ]
    errors.extend(
        f"forward_solver_verification missing: {name}" for name in verification_missing
    )
    if verification.get("status") != "verified":
        errors.append("forward_solver_verification.status must be 'verified'")

    metrics = artifact.get("metrics")
    if not isinstance(metrics, dict):
        errors.append("metrics must be an object with named values and units")
        metrics = {}
    metrics_have_units = bool(metrics) and all(
        isinstance(row, dict)
        and _present(row.get("value"))
        and _present(row.get("unit"))
        for row in metrics.values()
    )
    if not metrics_have_units:
        errors.append("every metric must contain value and unit")

    if family == "diffusion":
        try:
            if int(artifact.get("sampling_steps", 0)) <= 0:
                errors.append("sampling_steps must be positive")
            if int(artifact.get("generated_candidate_count", 0)) <= 0:
                errors.append("generated_candidate_count must be positive")
        except (TypeError, ValueError):
            errors.append("diffusion counts must be positive integers")
    elif family == "reinforcement_learning":
        try:
            if int(artifact.get("evaluation_episodes", 0)) <= 0:
                errors.append("evaluation_episodes must be positive")
        except (TypeError, ValueError):
            errors.append("evaluation_episodes must be a positive integer")
    elif family == "pseudoinverse":
        try:
            residual = float(artifact.get("residual_norm"))
            solution = float(artifact.get("solution_norm"))
            if residual < 0 or solution < 0:
                errors.append("residual_norm and solution_norm must be non-negative")
        except (TypeError, ValueError):
            errors.append("residual_norm and solution_norm must be numeric")

    return {
        "schema": "cae-ai-lab.artifact-gate/v1",
        "policy": "ai_candidate_requires_forward_physics_verification",
        "method_family": family,
        "status": "ok" if not errors else "needs_attention",
        "required_fields": list(required),
        "missing_fields": missing,
        "checks": {
            "common_metadata_complete": not any(name in missing for name in COMMON_REQUIREMENTS),
            "family_metadata_complete": not any(
                name in missing for name in FAMILY_REQUIREMENTS[family]
            ),
            "metrics_have_values_and_units": metrics_have_units,
            "forward_solver_verified": not verification_missing
            and verification.get("status") == "verified",
        },
        "errors": errors,
        "teaching_note": (
            "A generated, learned, or regularized candidate is not a CAE result until "
            "the owning forward solver reproduces named observables within saved tolerances."
        ),
    }


def cae_ai_artifact_gate(method_family: str, artifact_json: str) -> str:
    """JSON wrapper for the topology-optimization MCP server."""
    try:
        artifact = json.loads(artifact_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"artifact_json must be valid JSON: {exc.msg}") from exc
    if not isinstance(artifact, dict):
        raise ValueError("artifact_json must decode to an object")
    return json.dumps(
        evaluate_cae_ai_artifact(method_family, artifact), indent=2, sort_keys=True
    )
