from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v35 import _summary_v35

_PROMOTED_CASE_IDS = (
    "v36_public_cq_adaptive_timestep_contour_restart_interpolation_causality_energy_mismatch",
    "v36_public_fembem_shape_derivative_mesh_morph_normal_velocity_trace_jacobian_fd_mismatch",
)


def _summary_v36():
    payload = deepcopy(_summary_v35())
    generation = "cq-adaptive-412"
    steps = [1.0e-4, 5.0e-5, 5.0e-5, 1.0e-4]
    times = [0.0]
    for step in steps:
        times.append(times[-1] + step)
    payload[
        "cq_adaptive_timestep_contour_restart_interpolation_causality_energy_operator_result_identity"
    ] = {
        "adaptive_cq_generation": generation,
        **{
            key: generation
            for key in (
                "timestep_generation",
                "contour_generation",
                "restart_generation",
                "interpolation_generation",
                "causality_generation",
                "energy_generation",
                "operator_generation",
                "result_generation",
            )
        },
        "cq_method": "bdf2",
        "result_cq_method": "bdf2",
        "time_step_history_s": steps,
        "result_time_step_history_s": list(steps),
        "time_samples_s": times,
        "result_time_samples_s": list(times),
        "contour_radius": 0.8,
        "result_contour_radius": 0.8,
        "laplace_anchor_real_per_s": [0.22 / step for step in steps],
        "result_laplace_anchor_real_per_s": [0.22 / step for step in steps],
        "restart_step": 3,
        "result_restart_step": 3,
        "restart_history_state": [0.3, 0.2],
        "result_restart_history_state": [0.3, 0.2],
        "history_interpolation": "piecewise_linear_causal",
        "result_history_interpolation": "piecewise_linear_causal",
        "prehistory_max_abs": 0.0,
        "result_prehistory_max_abs": 0.0,
        "discrete_energy_j": [1.0, 0.8, 0.7, 0.6, 0.5],
        "result_discrete_energy_j": [1.0, 0.8, 0.7, 0.6, 0.5],
        "operator_owner": "cq/operator-412",
        "accepted_operator_owner": "cq/operator-412",
        "operator_sha256": "1" * 64,
        "accepted_operator_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "shape-derivative-412"
    step = 1.0e-3
    reference = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    velocity = [0.1, 0.2, 0.3]
    morphed = [
        [node[0], node[1], node[2] + step * normal_velocity]
        for node, normal_velocity in zip(reference, velocity)
    ]
    derivative = 0.6
    payload[
        "fembem_shape_derivative_morph_normal_velocity_trace_jacobian_objective_fd_mesh_owner_result_identity"
    ] = {
        "shape_generation": generation,
        **{
            key: generation
            for key in (
                "morph_generation",
                "normal_generation",
                "trace_generation",
                "jacobian_generation",
                "objective_generation",
                "fd_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "shape_step": step,
        "result_shape_step": step,
        "reference_nodes_m": reference,
        "result_reference_nodes_m": deepcopy(reference),
        "normal_velocity_m": velocity,
        "result_normal_velocity_m": list(velocity),
        "morphed_nodes_m": morphed,
        "result_morphed_nodes_m": deepcopy(morphed),
        "trace_node_map": [1, 2, 3],
        "result_trace_node_map": [1, 2, 3],
        "geometry_jacobian_determinant": [1.0, 1.0, 1.0],
        "result_geometry_jacobian_determinant": [1.0, 1.0, 1.0],
        "objective_directional_derivative": derivative,
        "result_objective_directional_derivative": derivative,
        "objective_minus": 2.0 - step * derivative,
        "result_objective_minus": 2.0 - step * derivative,
        "objective_plus": 2.0 + step * derivative,
        "result_objective_plus": 2.0 + step * derivative,
        "mesh_owner": "shape/mesh-412",
        "accepted_mesh_owner": "shape/mesh-412",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "shape_result_sha256": "4" * 64,
        "accepted_shape_result_sha256": "4" * 64,
    }
    return payload


def test_v36_public_positive_adaptive_cq_and_shape_derivative_closure():
    assert regularized_trace_inverse_path_gate(_summary_v36())["status"] == "ok"


def test_v36_public_adaptive_cq_manifest_mismatch_is_rejected():
    payload = _summary_v36()
    identity = payload[
        "cq_adaptive_timestep_contour_restart_interpolation_causality_energy_operator_result_identity"
    ]
    identity.update(
        {
            "timestep_generation": "cq-adaptive-411",
            "result_time_step_history_s": [-1.0],
            "result_time_samples_s": [1.0, 0.0],
            "result_restart_step": 99,
            "result_history_interpolation": "future_hold",
            "result_prehistory_max_abs": 1.0,
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "adaptive_cq_uses_current_timesteps_contour_restart_interpolation_causality_energy_operator_and_result"
    ]


def test_v36_public_shape_derivative_manifest_mismatch_is_rejected():
    payload = _summary_v36()
    identity = payload[
        "fembem_shape_derivative_morph_normal_velocity_trace_jacobian_objective_fd_mesh_owner_result_identity"
    ]
    identity.update(
        {
            "morph_generation": "shape-derivative-411",
            "result_shape_step": -1.0e-3,
            "result_morphed_nodes_m": [[9.0, 9.0, 9.0]],
            "result_geometry_jacobian_determinant": [-1.0],
            "result_objective_directional_derivative": -0.6,
            "accepted_shape_result_sha256": "c" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fembem_shape_derivative_uses_current_morph_normal_velocity_trace_jacobian_fd_mesh_and_result"
    ]


def test_v36_public_self_consistent_energy_growth_is_rejected():
    payload = _summary_v36()
    identity = payload[
        "cq_adaptive_timestep_contour_restart_interpolation_causality_energy_operator_result_identity"
    ]
    identity["discrete_energy_j"] = [1.0, 0.8, 0.9, 0.6, 0.5]
    identity["result_discrete_energy_j"] = list(identity["discrete_energy_j"])
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v36_public_self_consistent_incorrect_mesh_morph_is_rejected():
    payload = _summary_v36()
    identity = payload[
        "fembem_shape_derivative_morph_normal_velocity_trace_jacobian_objective_fd_mesh_owner_result_identity"
    ]
    identity["morphed_nodes_m"] = deepcopy(identity["reference_nodes_m"])
    identity["result_morphed_nodes_m"] = deepcopy(identity["reference_nodes_m"])
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"
