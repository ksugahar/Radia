from copy import deepcopy

from radia_mcp.radia_ngsolve.transform_normalization_v48 import validate_public_v48_identity


PROMOTED_CASE_IDS = {
    "v48_public_ale_reference_current_configuration_force_quadrature_owner_mismatch",
    "v48_public_segregated_solver_variable_scaling_residual_norm_iteration_solution_owner_mismatch",
}


def _records() -> dict[str, object]:
    ale_generation = "ale-force-v48"
    segregated_generation = "segregated-v48"
    groups = ["magnetic_vector_potential", "temperature", "displacement"]
    scaling = {"magnetic_vector_potential": 1.0, "temperature": 300.0, "displacement": 1.0e-3}
    iterations = ["iter=1|group=magnetic_vector_potential", "iter=1|group=temperature", "iter=1|group=displacement"]
    return {
        "ale_reference_current_force_quadrature_owner_identity": {
            "generation": ale_generation,
            "reference_mesh_generation": ale_generation,
            "current_mesh_generation": ale_generation,
            "quadrature_generation": ale_generation,
            "normal_generation": ale_generation,
            "result_generation": ale_generation,
            "reference_configuration_id": "ale/reference-v48",
            "result_reference_configuration_id": "ale/reference-v48",
            "current_configuration_id": "ale/current-v48",
            "result_current_configuration_id": "ale/current-v48",
            "quadrature_rule": "gauss-surface-order-4",
            "result_quadrature_rule": "gauss-surface-order-4",
            "normal_orientation_sha256": "1" * 64,
            "result_normal_orientation_sha256": "1" * 64,
            "body_owner": "body:moving-domain-v48",
            "result_body_owner": "body:moving-domain-v48",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
        "segregated_variable_scaling_residual_iteration_solution_identity": {
            "generation": segregated_generation,
            "variable_group_generation": segregated_generation,
            "scaling_generation": segregated_generation,
            "residual_generation": segregated_generation,
            "iteration_generation": segregated_generation,
            "solution_generation": segregated_generation,
            "result_generation": segregated_generation,
            "variable_groups": groups,
            "result_variable_groups": groups,
            "variable_scaling": scaling,
            "result_variable_scaling": scaling,
            "residual_norm": "scaled_l2",
            "result_residual_norm": "scaled_l2",
            "iteration_rows": iterations,
            "result_iteration_rows": iterations,
            "solution_owner": "solution:segregated-v48",
            "result_solution_owner": "solution:segregated-v48",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def test_v48_positive_replays_are_accepted() -> None:
    result = validate_public_v48_identity(_records())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v48_mixed_ale_configuration_is_rejected() -> None:
    records = deepcopy(_records())
    row = records["ale_reference_current_force_quadrature_owner_identity"]
    row["result_current_configuration_id"] = "ale/current-old"
    row["result_quadrature_rule"] = "gauss-surface-order-2"
    row["result_normal_orientation_sha256"] = "a" * 64
    row["result_body_owner"] = "body:fixed-old"
    assert validate_public_v48_identity(records)["status"] == "needs_attention"


def test_v48_mixed_segregated_solver_rows_are_rejected() -> None:
    records = deepcopy(_records())
    row = records["segregated_variable_scaling_residual_iteration_solution_identity"]
    row["result_variable_groups"] = list(reversed(row["variable_groups"]))
    row["result_variable_scaling"]["temperature"] = 1.0
    row["result_residual_norm"] = "unscaled_linf"
    row["result_iteration_rows"] = row["iteration_rows"][:-1]
    row["result_solution_owner"] = "solution:segregated-old"
    assert validate_public_v48_identity(records)["status"] == "needs_attention"
