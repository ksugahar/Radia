from __future__ import annotations

import math
from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v38 import _summary_v38


_PROMOTED_CASE_IDS = (
    "v39_public_nonlinear_fem_newton_consistent_tangent_linesearch_residual_energy_mesh_mismatch",
    "v39_public_cq_contour_frequency_interpolation_aliasing_passivity_error_timehistory_mismatch",
)
_NONLINEAR_KEY = (
    "nonlinear_fem_newton_residual_consistent_tangent_linesearch_step_energy_"
    "mesh_owner_result_identity"
)
_CQ_KEY = (
    "cq_contour_frequency_interpolation_aliasing_passivity_reconstruction_time_"
    "operator_result_identity"
)


def _summary_v39() -> dict:
    payload = deepcopy(_summary_v38())
    generation = "nonlinear-fem-715"
    residuals = [1.0, 0.2, 0.03, 1.0e-3, 1.0e-8]
    energy = [0.0, 0.8, 1.1, 1.18, 1.1801]
    mirrored = {
        "nonlinear_formulation": "total_lagrangian_hyperelastic",
        "residual_norm_history": residuals,
        "consistent_tangent_matrix": [[6.0, -2.0], [-2.0, 4.0]],
        "directional_tangent_product": [4.0, 0.0],
        "finite_difference_directional_derivative": [4.000000001, 0.0],
        "tangent_relative_tolerance": 1.0e-8,
        "newton_step_history": [
            [-0.2, 0.1],
            [-0.04, 0.02],
            [-0.006, 0.003],
            [-0.0002, 0.0001],
        ],
        "line_search_alpha_history": [1.0, 1.0, 1.0, 1.0],
        "line_search_trial_residual_norm": residuals[1:],
        "line_search_armijo_constant": 1.0e-4,
        "strain_energy_history_j": energy,
        "external_work_final_j": energy[-1],
        "energy_balance_residual_j": 0.0,
        "convergence_tolerance": 1.0e-7,
        "nonlinear_mesh_sha256": "1" * 64,
    }
    payload[_NONLINEAR_KEY] = {
        "nonlinear_generation": generation,
        **{
            key: generation
            for key in (
                "residual_generation",
                "tangent_generation",
                "step_generation",
                "linesearch_generation",
                "iteration_generation",
                "energy_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "nonlinear_owner": "fem/nonlinear-715",
        "accepted_nonlinear_owner": "fem/nonlinear-715",
        "nonlinear_result_sha256": "2" * 64,
        "accepted_nonlinear_result_sha256": "2" * 64,
    }

    generation = "cq-contour-715"
    count, radius = 8, 0.92
    contour = [
        [
            radius * math.cos(2.0 * math.pi * index / count),
            radius * math.sin(2.0 * math.pi * index / count),
        ]
        for index in range(count)
    ]
    history = [0.0, 1.0, 0.6, 0.3, 0.12, 0.04, 0.01, 0.0]
    mirrored = {
        "cq_method": "bdf2",
        "time_step_s": 1.0e-4,
        "time_step_count": count,
        "contour_radius": radius,
        "contour_nodes_complex": contour,
        "frequency_interpolation_relative_error": 2.0e-4,
        "maximum_frequency_interpolation_relative_error": 1.0e-3,
        "aliasing_error_bound": 5.0e-5,
        "maximum_aliasing_error": 1.0e-4,
        "minimum_transfer_passivity_eigenvalue": 0.02,
        "time_reconstruction_relative_error": 3.0e-4,
        "maximum_time_reconstruction_relative_error": 1.0e-3,
        "time_history": history,
        "reconstructed_time_history": history,
        "cq_operator_sha256": "3" * 64,
    }
    payload[_CQ_KEY] = {
        "cq_generation": generation,
        **{
            key: generation
            for key in (
                "contour_generation",
                "frequency_generation",
                "interpolation_generation",
                "aliasing_generation",
                "passivity_generation",
                "reconstruction_generation",
                "time_generation",
                "operator_generation",
                "result_generation",
            )
        },
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "operator_owner": "cq/operator-715",
        "accepted_operator_owner": "cq/operator-715",
        "cq_result_sha256": "4" * 64,
        "accepted_cq_result_sha256": "4" * 64,
    }
    return payload


def test_v39_public_positive_nonlinear_fem_and_cq_contour_closure() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v39())["status"] == "ok"


def test_v39_public_nonlinear_fem_newton_consistent_tangent_linesearch_residual_energy_mesh_mismatch() -> None:
    payload = _summary_v39()
    identity = payload[_NONLINEAR_KEY]
    identity.update(
        {
            "tangent_generation": "nonlinear-fem-714",
            "energy_generation": "nonlinear-fem-713",
            "result_generation": "nonlinear-fem-712",
            "result_residual_norm_history": [1.0, 2.0],
            "result_consistent_tangent_matrix": [[-1.0]],
            "result_directional_tangent_product": [9.0],
            "result_line_search_alpha_history": [-1.0],
            "result_line_search_trial_residual_norm": [2.0],
            "result_strain_energy_history_j": [1.0, -1.0],
            "result_energy_balance_residual_j": 9.0,
            "accepted_nonlinear_owner": "fem/old",
            "accepted_nonlinear_result_sha256": "a" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_fem_uses_current_residual_tangent_newton_linesearch_energy_mesh_owner_and_result"
    ]


def test_v39_public_cq_contour_frequency_interpolation_aliasing_passivity_error_timehistory_mismatch() -> None:
    payload = _summary_v39()
    identity = payload[_CQ_KEY]
    identity.update(
        {
            "contour_generation": "cq-contour-714",
            "passivity_generation": "cq-contour-713",
            "result_generation": "cq-contour-712",
            "result_contour_nodes_complex": [[2.0, 0.0]],
            "result_frequency_interpolation_relative_error": 1.0,
            "result_aliasing_error_bound": 1.0,
            "result_minimum_transfer_passivity_eigenvalue": -1.0,
            "result_time_reconstruction_relative_error": 1.0,
            "result_reconstructed_time_history": [9.0],
            "accepted_operator_owner": "cq/old",
            "accepted_cq_result_sha256": "b" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_contour_uses_current_nodes_interpolation_aliasing_passivity_reconstruction_time_operator_and_result"
    ]


def test_v39_public_rejects_self_consistent_increasing_newton_residual() -> None:
    payload = _summary_v39()
    identity = payload[_NONLINEAR_KEY]
    residuals = [1.0, 0.2, 0.3, 1.0e-3, 1.0e-8]
    identity["residual_norm_history"] = residuals
    identity["result_residual_norm_history"] = residuals
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_non_circular_cq_contour() -> None:
    payload = _summary_v39()
    identity = payload[_CQ_KEY]
    contour = deepcopy(identity["contour_nodes_complex"])
    contour[1][0] *= 0.5
    identity["contour_nodes_complex"] = contour
    identity["result_contour_nodes_complex"] = contour
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v39_public_accepts_general_symmetric_positive_definite_tangent() -> None:
    payload = _summary_v39()
    identity = payload[_NONLINEAR_KEY]
    tangent = [[5.0, 1.0, 0.0], [1.0, 4.0, 1.0], [0.0, 1.0, 3.0]]
    product = [5.0, 1.0, 0.0]
    steps = [[*row, 0.0] for row in identity["newton_step_history"]]
    identity["consistent_tangent_matrix"] = tangent
    identity["result_consistent_tangent_matrix"] = tangent
    identity["directional_tangent_product"] = product
    identity["result_directional_tangent_product"] = product
    identity["finite_difference_directional_derivative"] = product
    identity["result_finite_difference_directional_derivative"] = product
    identity["newton_step_history"] = steps
    identity["result_newton_step_history"] = steps
    assert regularized_trace_inverse_path_gate(payload)["status"] == "ok"
