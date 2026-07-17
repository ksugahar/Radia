from __future__ import annotations

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v24 import _summary_v24


_PROMOTED_CASE_IDS = (
    "v25_public_hmatrix_block_tree_admissibility_permutation_tolerance_kernel_mesh_mismatch",
    "v25_public_ad_parameter_tape_material_operator_mesh_objective_gradient_mismatch",
)


def _summary_v25():
    summary = _summary_v24()
    summary[
        "hmatrix_block_tree_admissibility_permutation_tolerance_kernel_mesh_generation_identity"
    ] = {
        "hmatrix_generation": "hmatrix-201",
        "block_tree_hmatrix_generation": "hmatrix-201",
        "admissibility_hmatrix_generation": "hmatrix-201",
        "permutation_hmatrix_generation": "hmatrix-201",
        "tolerance_hmatrix_generation": "hmatrix-201",
        "kernel_hmatrix_generation": "hmatrix-201",
        "mesh_hmatrix_generation": "hmatrix-201",
        "result_hmatrix_generation": "hmatrix-201",
        "matrix_shape": [4, 4],
        "result_matrix_shape": [4, 4],
        "block_tree_sha256": "1" * 64,
        "result_block_tree_sha256": "1" * 64,
        "admissibility_rule": "diameter_le_eta_distance",
        "result_admissibility_rule": "diameter_le_eta_distance",
        "admissibility_eta": 1.5,
        "result_admissibility_eta": 1.5,
        "row_permutation": [2, 0, 3, 1],
        "result_row_permutation": [2, 0, 3, 1],
        "column_permutation": [1, 3, 0, 2],
        "result_column_permutation": [1, 3, 0, 2],
        "relative_tolerance": 1.0e-6,
        "result_relative_tolerance": 1.0e-6,
        "kernel_id": "helmholtz_single_layer_p1",
        "result_kernel_id": "helmholtz_single_layer_p1",
        "boundary_mesh_sha256": "2" * 64,
        "result_boundary_mesh_sha256": "2" * 64,
        "hmatrix_result_sha256": "3" * 64,
        "reported_hmatrix_result_sha256": "3" * 64,
    }
    summary[
        "ad_parameter_tape_material_operator_mesh_objective_primal_gradient_generation_identity"
    ] = {
        "ad_generation": "ad-gradient-201",
        "parameter_tape_ad_generation": "ad-gradient-201",
        "material_law_ad_generation": "ad-gradient-201",
        "operator_ad_generation": "ad-gradient-201",
        "mesh_ad_generation": "ad-gradient-201",
        "objective_ad_generation": "ad-gradient-201",
        "primal_ad_generation": "ad-gradient-201",
        "gradient_ad_generation": "ad-gradient-201",
        "result_ad_generation": "ad-gradient-201",
        "parameter_names": ["density", "bulk_modulus"],
        "result_parameter_names": ["density", "bulk_modulus"],
        "parameter_values": [1.2, 142000.0],
        "result_parameter_values": [1.2, 142000.0],
        "material_law_id": "linear_acoustic_fluid",
        "result_material_law_id": "linear_acoustic_fluid",
        "objective_id": "receiver_pressure_l2",
        "result_objective_id": "receiver_pressure_l2",
        "parameter_tape_sha256": "4" * 64,
        "result_parameter_tape_sha256": "4" * 64,
        "assembled_operator_sha256": "5" * 64,
        "result_assembled_operator_sha256": "5" * 64,
        "mesh_sha256": "6" * 64,
        "result_mesh_sha256": "6" * 64,
        "primal_solution_sha256": "7" * 64,
        "result_primal_solution_sha256": "7" * 64,
        "ad_gradient": [0.25, -0.004],
        "finite_difference_gradient": [0.25000001, -0.0040000001],
        "maximum_gradient_relative_error": 4.0e-8,
        "gradient_relative_tolerance": 1.0e-5,
        "gradient_result_sha256": "8" * 64,
        "reported_gradient_result_sha256": "8" * 64,
    }
    return summary


def test_v25_public_positive_hmatrix_and_ad_identity() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v25())["status"] == "ok"


def test_v25_public_hmatrix_identity_mismatch() -> None:
    summary = _summary_v25()
    identity = summary[
        "hmatrix_block_tree_admissibility_permutation_tolerance_kernel_mesh_generation_identity"
    ]
    identity.update(
        {
            "block_tree_hmatrix_generation": "hmatrix-200",
            "admissibility_hmatrix_generation": "hmatrix-199",
            "permutation_hmatrix_generation": "hmatrix-198",
            "tolerance_hmatrix_generation": "hmatrix-197",
            "kernel_hmatrix_generation": "hmatrix-196",
            "mesh_hmatrix_generation": "hmatrix-195",
            "result_matrix_shape": [5, 4],
            "result_block_tree_sha256": "f" * 64,
            "result_admissibility_rule": "always_admissible",
            "result_admissibility_eta": 3.0,
            "result_row_permutation": [2, 0, 3, 3],
            "result_relative_tolerance": 1.0e-2,
            "result_kernel_id": "laplace_single_layer_p0",
            "result_boundary_mesh_sha256": "0" * 64,
            "reported_hmatrix_result_sha256": "1" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hmatrix_uses_current_block_tree_admissibility_permutations_tolerance_kernel_and_mesh"
    ]


def test_v25_public_ad_gradient_identity_mismatch() -> None:
    summary = _summary_v25()
    identity = summary[
        "ad_parameter_tape_material_operator_mesh_objective_primal_gradient_generation_identity"
    ]
    identity.update(
        {
            "parameter_tape_ad_generation": "ad-gradient-200",
            "material_law_ad_generation": "ad-gradient-199",
            "operator_ad_generation": "ad-gradient-198",
            "mesh_ad_generation": "ad-gradient-197",
            "objective_ad_generation": "ad-gradient-196",
            "primal_ad_generation": "ad-gradient-195",
            "result_parameter_names": ["bulk_modulus", "density"],
            "result_parameter_values": [142000.0, 1.2],
            "result_material_law_id": "nonlinear_fluid_previous",
            "result_objective_id": "source_power_previous",
            "result_parameter_tape_sha256": "2" * 64,
            "result_assembled_operator_sha256": "3" * 64,
            "result_mesh_sha256": "4" * 64,
            "result_primal_solution_sha256": "5" * 64,
            "finite_difference_gradient": [-0.25, 0.04],
            "maximum_gradient_relative_error": 2.0,
            "reported_gradient_result_sha256": "6" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "ad_gradient_uses_current_tape_material_operator_mesh_objective_and_primal"
    ]
