from __future__ import annotations

import math
from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v34 import _summary_v34

_PROMOTED_CASE_IDS = (
    "v35_public_cq_acoustic_laplace_contour_weight_passivity_trace_timestep_history_mismatch",
    "v35_public_fembem_autodiff_complex_gradient_wirtinger_shape_fd_objective_mismatch",
)


def _summary_v35():
    payload = deepcopy(_summary_v34())
    generation = "cq-acoustic-411"
    timestep = 1.0e-4
    radius = 0.8
    zeta = []
    laplace = []
    for index in range(4):
        angle = -2.0 * math.pi * index / 4
        point = radius * complex(math.cos(angle), math.sin(angle))
        transformed = (1.5 - 2.0 * point + 0.5 * point * point) / timestep
        zeta.append([point.real, point.imag])
        laplace.append([transformed.real, transformed.imag])
    payload[
        "cq_acoustic_laplace_contour_weight_passivity_trace_timestep_history_mesh_owner_result_identity"
    ] = {
        "cq_generation": generation,
        **{
            key: generation
            for key in (
                "contour_cq_generation",
                "weight_cq_generation",
                "passivity_cq_generation",
                "trace_cq_generation",
                "timestep_cq_generation",
                "history_cq_generation",
                "mesh_cq_generation",
                "owner_cq_generation",
                "result_cq_generation",
            )
        },
        "cq_method": "bdf2",
        "result_cq_method": "bdf2",
        "contour_radius": radius,
        "result_contour_radius": radius,
        "zeta_points_ri": zeta,
        "result_zeta_points_ri": deepcopy(zeta),
        "laplace_points_ri": laplace,
        "result_laplace_points_ri": deepcopy(laplace),
        "cq_weights": [0.4, 0.3, 0.2, 0.1],
        "result_cq_weights": [0.4, 0.3, 0.2, 0.1],
        "boundary_impedance_ri": [[2.0, 0.1], [2.2, 0.2], [2.4, 0.3], [2.6, 0.4]],
        "result_boundary_impedance_ri": [[2.0, 0.1], [2.2, 0.2], [2.4, 0.3], [2.6, 0.4]],
        "trace_orientation": "outward_volume_to_boundary",
        "result_trace_orientation": "outward_volume_to_boundary",
        "fem_trace_node_ids": [1, 4, 7],
        "result_fem_trace_node_ids": [1, 4, 7],
        "bem_trace_node_ids": [1, 4, 7],
        "result_bem_trace_node_ids": [1, 4, 7],
        "trace_sign": 1,
        "result_trace_sign": 1,
        "time_step_s": timestep,
        "result_time_step_s": timestep,
        "history_length": 4,
        "result_history_length": 4,
        "time_samples_s": [index * timestep for index in range(4)],
        "result_time_samples_s": [index * timestep for index in range(4)],
        "mesh_owner": "cq-acoustic/mesh-411",
        "accepted_mesh_owner": "cq-acoustic/mesh-411",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_owner": "cq-acoustic/result-411",
        "accepted_result_owner": "cq-acoustic/result-411",
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "fembem-autodiff-411"
    design = [1.2, -0.4]
    direction = [0.6, 0.8]
    objective = 0.5 * sum(item * item for item in design)
    directional = sum(item * tangent for item, tangent in zip(design, direction))
    payload["fembem_autodiff_wirtinger_objective_shape_fd_trace_mesh_owner_gradient_identity"] = {
        "autodiff_generation": generation,
        **{
            key: generation
            for key in (
                "complex_autodiff_generation",
                "wirtinger_autodiff_generation",
                "objective_autodiff_generation",
                "shape_autodiff_generation",
                "finite_difference_autodiff_generation",
                "trace_autodiff_generation",
                "mesh_autodiff_generation",
                "owner_autodiff_generation",
                "gradient_autodiff_generation",
                "result_autodiff_generation",
            )
        },
        "complex_design_ri": design,
        "result_complex_design_ri": list(design),
        "objective_scaling": "one_half_l2_squared",
        "result_objective_scaling": "one_half_l2_squared",
        "objective_value": objective,
        "result_objective_value": objective,
        "wirtinger_convention": "dJ_dconjugate_z",
        "result_wirtinger_convention": "dJ_dconjugate_z",
        "wirtinger_gradient_ri": [0.5 * item for item in design],
        "result_wirtinger_gradient_ri": [0.5 * item for item in design],
        "real_gradient_ri": design,
        "result_real_gradient_ri": list(design),
        "shape_direction_ri": direction,
        "result_shape_direction_ri": list(direction),
        "shape_step": 1.0e-6,
        "result_shape_step": 1.0e-6,
        "finite_difference_directional_derivative": directional,
        "result_finite_difference_directional_derivative": directional,
        "trace_node_map": [1, 4, 7],
        "result_trace_node_map": [1, 4, 7],
        "mesh_owner": "fembem-autodiff/mesh-411",
        "accepted_mesh_owner": "fembem-autodiff/mesh-411",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "gradient_owner": "fembem-autodiff/gradient-411",
        "accepted_gradient_owner": "fembem-autodiff/gradient-411",
        "gradient_sha256": "4" * 64,
        "accepted_gradient_sha256": "4" * 64,
    }
    return payload


def test_v35_public_positive_cq_and_fembem_autodiff_closure():
    assert regularized_trace_inverse_path_gate(_summary_v35())["status"] == "ok"


def test_v35_public_cq_manifest_mismatch_is_rejected():
    payload = _summary_v35()
    identity = payload[
        "cq_acoustic_laplace_contour_weight_passivity_trace_timestep_history_mesh_owner_result_identity"
    ]
    identity.update(
        {
            "contour_cq_generation": "cq-acoustic-410",
            "result_cq_method": "backward_euler",
            "result_laplace_points_ri": [[-1.0, 0.0]],
            "result_boundary_impedance_ri": [[-2.0, 0.0]],
            "result_trace_sign": -1,
            "result_time_step_s": -1.0e-4,
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_acoustics_use_current_bdf2_laplace_contour_weights_passivity_trace_timestep_history_mesh_and_result"
    ]


def test_v35_public_autodiff_manifest_mismatch_is_rejected():
    payload = _summary_v35()
    identity = payload[
        "fembem_autodiff_wirtinger_objective_shape_fd_trace_mesh_owner_gradient_identity"
    ]
    identity.update(
        {
            "complex_autodiff_generation": "fembem-autodiff-410",
            "result_objective_scaling": "l2_squared",
            "result_wirtinger_convention": "dJ_dz",
            "result_real_gradient_ri": [0.0, 0.0],
            "result_shape_step": -1.0e-6,
            "accepted_gradient_sha256": "b" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fembem_autodiff_uses_current_wirtinger_objective_shape_fd_trace_mesh_and_gradient"
    ]


def test_v35_public_self_consistent_wrong_cq_laplace_symbol_is_rejected():
    payload = _summary_v35()
    identity = payload[
        "cq_acoustic_laplace_contour_weight_passivity_trace_timestep_history_mesh_owner_result_identity"
    ]
    wrong = [[row[0] + 10.0, row[1]] for row in identity["laplace_points_ri"]]
    identity["laplace_points_ri"] = wrong
    identity["result_laplace_points_ri"] = deepcopy(wrong)
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v35_public_self_consistent_wrong_wirtinger_gradient_is_rejected():
    payload = _summary_v35()
    identity = payload[
        "fembem_autodiff_wirtinger_objective_shape_fd_trace_mesh_owner_gradient_identity"
    ]
    identity["wirtinger_gradient_ri"] = [1.2, -0.4]
    identity["result_wirtinger_gradient_ri"] = [1.2, -0.4]
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"
