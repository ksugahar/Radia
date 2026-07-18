from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v33 import _summary_v33


_PROMOTED_CASE_IDS = (
    "v34_public_fembem_reciprocity_radiation_power_trace_orientation_mesh_solution_mismatch",
    "v34_public_nonlinear_eigenvalue_contour_moment_rank_residual_biorthogonality_mismatch",
)


def _summary_v34():
    payload = deepcopy(_summary_v33())
    generation = "fembem-coupled-391"
    payload[
        "fembem_reciprocity_radiation_power_interior_energy_trace_orientation_boundary_volume_map_frequency_mesh_solution_identity"
    ] = {
        "fembem_generation": generation,
        **{
            key: generation
            for key in (
                "transfer_fembem_generation", "radiation_fembem_generation",
                "interior_fembem_generation", "trace_fembem_generation",
                "map_fembem_generation", "frequency_fembem_generation",
                "mesh_fembem_generation", "solution_fembem_generation",
                "result_fembem_generation",
            )
        },
        "frequency_hz": 1000.0, "result_frequency_hz": 1000.0,
        "transfer_ab_ri": [0.2, 0.1], "result_transfer_ab_ri": [0.2, 0.1],
        "transfer_ba_ri": [0.2, 0.1], "result_transfer_ba_ri": [0.2, 0.1],
        "reciprocity_tolerance": 1.0e-8, "result_reciprocity_tolerance": 1.0e-8,
        "radiated_power_w": 0.5, "result_radiated_power_w": 0.5,
        "boundary_flux_power_w": 0.5, "result_boundary_flux_power_w": 0.5,
        "interior_energy_j": 0.25, "result_interior_energy_j": 0.25,
        "trace_orientation": "outward_volume_to_boundary",
        "result_trace_orientation": "outward_volume_to_boundary",
        "boundary_volume_node_map": [1, 4, 7],
        "result_boundary_volume_node_map": [1, 4, 7],
        "trace_node_ids": [1, 4, 7], "result_trace_node_ids": [1, 4, 7],
        "mesh_owner": "fembem/mesh-391", "accepted_mesh_owner": "fembem/mesh-391",
        "mesh_sha256": "1" * 64, "accepted_mesh_sha256": "1" * 64,
        "solution_owner": "fembem/solution-391",
        "accepted_solution_owner": "fembem/solution-391",
        "solution_sha256": "2" * 64, "accepted_solution_sha256": "2" * 64,
    }
    generation = "nonlinear-eigen-contour-391"
    contour = [[1.5, -1.0], [2.5, 0.0], [1.5, 1.0], [0.5, 0.0]]
    gram = [[[1.0, 0.0], [0.0, 0.0]], [[0.0, 0.0], [1.0, 0.0]]]
    payload[
        "nonlinear_eigen_contour_orientation_quadrature_moment_rank_count_residual_biorthogonality_pole_result_identity"
    ] = {
        "nonlinear_eigen_generation": generation,
        **{
            key: generation
            for key in (
                "contour_eigen_generation", "quadrature_eigen_generation",
                "moment_eigen_generation", "rank_eigen_generation",
                "count_eigen_generation", "residual_eigen_generation",
                "biorthogonality_eigen_generation", "pole_eigen_generation",
                "result_eigen_generation",
            )
        },
        "contour_orientation": "counterclockwise",
        "result_contour_orientation": "counterclockwise",
        "contour_points_ri": contour, "result_contour_points_ri": deepcopy(contour),
        "quadrature_rule": "trapezoidal_periodic",
        "result_quadrature_rule": "trapezoidal_periodic",
        "moment_ranks": [2, 2], "result_moment_ranks": [2, 2],
        "numerical_rank": 2, "result_numerical_rank": 2,
        "enclosed_eigenvalue_count": 2, "result_enclosed_eigenvalue_count": 2,
        "eigenvalues_ri": [[1.0, 0.0], [2.0, 0.0]],
        "result_eigenvalues_ri": [[1.0, 0.0], [2.0, 0.0]],
        "residual_norms": [1.0e-10, 2.0e-10],
        "result_residual_norms": [1.0e-10, 2.0e-10],
        "biorthogonality_gram_ri": gram,
        "result_biorthogonality_gram_ri": deepcopy(gram),
        "pole_owner": "nonlinear-eigen/poles-391",
        "accepted_pole_owner": "nonlinear-eigen/poles-391",
        "pole_sha256": "3" * 64, "accepted_pole_sha256": "3" * 64,
        "result_owner": "nonlinear-eigen/result-391",
        "accepted_result_owner": "nonlinear-eigen/result-391",
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v34_public_positive_fembem_and_nonlinear_eigen_closure():
    assert regularized_trace_inverse_path_gate(_summary_v34())["status"] == "ok"


def test_v34_public_fembem_reciprocity_radiation_power_trace_orientation_mesh_solution_mismatch():
    payload = _summary_v34()
    identity = payload[
        "fembem_reciprocity_radiation_power_interior_energy_trace_orientation_boundary_volume_map_frequency_mesh_solution_identity"
    ]
    identity.update({
        "transfer_fembem_generation": "fembem-coupled-390",
        "map_fembem_generation": "fembem-coupled-389",
        "result_fembem_generation": "fembem-coupled-388",
        "result_frequency_hz": 900.0, "result_transfer_ab_ri": [0.3, 0.2],
        "result_transfer_ba_ri": [0.1, -0.2], "result_reciprocity_tolerance": 1.0e-3,
        "result_radiated_power_w": -0.2, "result_boundary_flux_power_w": 0.8,
        "result_interior_energy_j": -0.1,
        "result_trace_orientation": "inward_boundary_to_volume",
        "result_boundary_volume_node_map": [7, 4, 1],
        "result_trace_node_ids": [1, 2, 3],
        "accepted_mesh_owner": "fembem/old-mesh", "accepted_mesh_sha256": "b" * 64,
        "accepted_solution_owner": "fembem/old-solution",
        "accepted_solution_sha256": "c" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fembem_uses_current_reciprocal_transfer_radiation_power_interior_energy_trace_map_frequency_mesh_and_solution"
    ]


def test_v34_public_nonlinear_eigenvalue_contour_moment_rank_residual_biorthogonality_mismatch():
    payload = _summary_v34()
    identity = payload[
        "nonlinear_eigen_contour_orientation_quadrature_moment_rank_count_residual_biorthogonality_pole_result_identity"
    ]
    identity.update({
        "contour_eigen_generation": "nonlinear-eigen-contour-390",
        "rank_eigen_generation": "nonlinear-eigen-contour-389",
        "result_eigen_generation": "nonlinear-eigen-contour-388",
        "result_contour_orientation": "clockwise",
        "result_contour_points_ri": [[0.5, 0.0], [1.5, 1.0]],
        "result_quadrature_rule": "open_newton_cotes",
        "result_moment_ranks": [3, 1], "result_numerical_rank": 3,
        "result_enclosed_eigenvalue_count": 1,
        "result_eigenvalues_ri": [[3.0, 1.0]], "result_residual_norms": [0.2],
        "result_biorthogonality_gram_ri": [[[0.0, 1.0]]],
        "accepted_pole_owner": "nonlinear-eigen/old-poles",
        "accepted_pole_sha256": "d" * 64,
        "accepted_result_owner": "nonlinear-eigen/old-result",
        "accepted_result_sha256": "e" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_eigenpairs_use_current_contour_orientation_quadrature_moments_rank_count_residual_biorthogonality_poles_and_result"
    ]


def test_v34_public_self_consistent_nonreciprocal_transfer_is_rejected():
    payload = _summary_v34()
    identity = payload[
        "fembem_reciprocity_radiation_power_interior_energy_trace_orientation_boundary_volume_map_frequency_mesh_solution_identity"
    ]
    identity["transfer_ba_ri"] = [0.1, -0.2]
    identity["result_transfer_ba_ri"] = [0.1, -0.2]
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v34_public_self_consistent_clockwise_contour_is_rejected():
    payload = _summary_v34()
    identity = payload[
        "nonlinear_eigen_contour_orientation_quadrature_moment_rank_count_residual_biorthogonality_pole_result_identity"
    ]
    points = list(reversed(identity["contour_points_ri"]))
    identity["contour_points_ri"] = points
    identity["result_contour_points_ri"] = deepcopy(points)
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"
