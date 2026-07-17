from __future__ import annotations

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v26 import _summary_v26


_PROMOTED_CASE_IDS = (
    "v27_public_cq_adaptive_contour_quadrature_order_startup_correction_error_estimator_restart_mismatch",
    "v27_public_p1_fembem_boundary_orientation_quadrature_singular_treatment_trace_matrix_mismatch",
)


def _summary_v27():
    summary = _summary_v26()
    generation = "adaptive-cq-311"
    summary[
        "cq_adaptive_contour_quadrature_order_startup_correction_error_estimator_restart_generation_identity"
    ] = {
        "cq_generation": generation,
        "contour_cq_generation": generation,
        "quadrature_order_cq_generation": generation,
        "startup_correction_cq_generation": generation,
        "error_estimator_cq_generation": generation,
        "restart_cq_generation": generation,
        "result_cq_generation": generation,
        "contour_family": "lubich_bdf2_circle",
        "result_contour_family": "lubich_bdf2_circle",
        "contour_radii": [0.82, 0.9, 0.95],
        "result_contour_radii": [0.82, 0.9, 0.95],
        "quadrature_orders": [16, 32, 64],
        "result_quadrature_orders": [16, 32, 64],
        "startup_correction": "bdf2_consistent_two_step",
        "result_startup_correction": "bdf2_consistent_two_step",
        "error_estimator": "successive_contour_l2_relative",
        "result_error_estimator": "successive_contour_l2_relative",
        "relative_tolerance": 1.0e-5,
        "result_relative_tolerance": 1.0e-5,
        "estimated_relative_errors": [2.0e-3, 1.5e-4, 8.0e-6],
        "result_estimated_relative_errors": [2.0e-3, 1.5e-4, 8.0e-6],
        "restart_step": 40,
        "result_restart_step": 40,
        "restart_state_sha256": "1" * 64,
        "loaded_restart_state_sha256": "1" * 64,
        "time_history_sha256": "2" * 64,
        "accepted_time_history_sha256": "2" * 64,
    }
    generation = "p1-fembem-311"
    summary[
        "p1_fembem_boundary_orientation_quadrature_singular_treatment_trace_matrix_mesh_generation_identity"
    ] = {
        "coupling_generation": generation,
        "boundary_orientation_coupling_generation": generation,
        "quadrature_coupling_generation": generation,
        "singular_treatment_coupling_generation": generation,
        "trace_coupling_generation": generation,
        "matrix_coupling_generation": generation,
        "mesh_coupling_generation": generation,
        "result_coupling_generation": generation,
        "fem_basis_order": 1,
        "bem_basis_order": 1,
        "volume_element": "tet",
        "boundary_element": "tri",
        "boundary_orientation": "volume_outward",
        "result_boundary_orientation": "volume_outward",
        "regular_quadrature": "triangle_degree_4",
        "result_regular_quadrature": "triangle_degree_4",
        "singular_treatment": "duffy_p1_galerkin",
        "result_singular_treatment": "duffy_p1_galerkin",
        "trace_shape": [48, 120],
        "result_trace_shape": [48, 120],
        "trace_matrix_sha256": "3" * 64,
        "result_trace_matrix_sha256": "3" * 64,
        "fem_matrix_sha256": "4" * 64,
        "result_fem_matrix_sha256": "4" * 64,
        "bem_matrix_sha256": "5" * 64,
        "result_bem_matrix_sha256": "5" * 64,
        "volume_mesh_sha256": "6" * 64,
        "result_volume_mesh_sha256": "6" * 64,
        "boundary_mesh_sha256": "7" * 64,
        "result_boundary_mesh_sha256": "7" * 64,
        "coupled_result_sha256": "8" * 64,
        "accepted_coupled_result_sha256": "8" * 64,
    }
    return summary


def test_v27_public_positive_adaptive_cq_and_p1_fembem_identities() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v27())["status"] == "ok"


def test_v27_public_rejects_adaptive_cq_identity_mismatch() -> None:
    summary = _summary_v27()
    identity = summary[
        "cq_adaptive_contour_quadrature_order_startup_correction_error_estimator_restart_generation_identity"
    ]
    identity.update(
        {
            "contour_cq_generation": "adaptive-cq-310",
            "quadrature_order_cq_generation": "adaptive-cq-309",
            "restart_cq_generation": "adaptive-cq-308",
            "result_contour_family": "talbot_untracked",
            "result_contour_radii": [0.7],
            "result_quadrature_orders": [24, 48],
            "result_startup_correction": "none",
            "result_error_estimator": "absolute_peak",
            "result_relative_tolerance": 1.0e-2,
            "result_estimated_relative_errors": [1.0e-1],
            "result_restart_step": 20,
            "loaded_restart_state_sha256": "f" * 64,
            "accepted_time_history_sha256": "0" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "adaptive_cq_uses_current_contour_order_startup_error_estimator_restart_and_history"
    ]


def test_v27_public_rejects_p1_fembem_identity_mismatch() -> None:
    summary = _summary_v27()
    identity = summary[
        "p1_fembem_boundary_orientation_quadrature_singular_treatment_trace_matrix_mesh_generation_identity"
    ]
    identity.update(
        {
            "boundary_orientation_coupling_generation": "p1-fembem-310",
            "quadrature_coupling_generation": "p1-fembem-309",
            "trace_coupling_generation": "p1-fembem-308",
            "bem_basis_order": 0,
            "boundary_element": "quad",
            "result_boundary_orientation": "volume_inward",
            "result_regular_quadrature": "triangle_centroid",
            "result_singular_treatment": "diagonal_zeroed",
            "result_trace_shape": [47, 120],
            "result_trace_matrix_sha256": "1" * 64,
            "result_fem_matrix_sha256": "2" * 64,
            "result_bem_matrix_sha256": "3" * 64,
            "result_volume_mesh_sha256": "4" * 64,
            "result_boundary_mesh_sha256": "5" * 64,
            "accepted_coupled_result_sha256": "6" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "p1_fembem_uses_current_boundary_orientation_quadrature_singular_trace_matrices_and_mesh"
    ]
