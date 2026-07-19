from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v39 import _summary_v39


_PROMOTED_CASE_IDS = (
    "v40_public_johnson_nedelec_fembem_trace_normal_orientation_operator_energy_residual_mismatch",
    "v40_public_adjoint_hessian_design_gradient_vector_product_kkt_fd_owner_mismatch",
)
_FEMBEM_KEY = (
    "johnson_nedelec_volume_trace_normal_single_double_layer_sign_residual_"
    "energy_mesh_owner_result_identity"
)
_OPTIMIZATION_KEY = (
    "adjoint_hessian_design_objective_constraint_hvp_kkt_fd_model_owner_"
    "result_identity"
)


def _summary_v40() -> dict:
    payload = deepcopy(_summary_v39())
    generation = "johnson-nedelec-724"
    mirrored = {
        "volume_trace_matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        "boundary_normals": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        "outward_reference_vectors": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        "single_layer_matrix": [[2.0, 0.25], [0.25, 1.5]],
        "double_layer_matrix": [[0.1, -0.05], [0.02, 0.08]],
        "coupling_sign": -1.0,
        "interface_residual_vector": [1.0e-10, -1.0e-10],
        "interface_residual_tolerance": 1.0e-8,
        "interior_energy_flux_w": 2.5,
        "exterior_energy_flux_w": -2.5,
        "energy_flux_residual_w": 0.0,
        "fembem_mesh_sha256": "1" * 64,
    }
    payload[_FEMBEM_KEY] = {
        "fembem_generation": generation,
        **{key: generation for key in (
            "trace_generation", "normal_generation", "single_layer_generation",
            "double_layer_generation", "coupling_generation", "residual_generation",
            "energy_generation", "mesh_generation", "owner_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "fembem_owner": "acoustic/fembem-724",
        "accepted_fembem_owner": "acoustic/fembem-724",
        "fembem_result_sha256": "2" * 64,
        "accepted_fembem_result_sha256": "2" * 64,
    }

    generation = "adjoint-hessian-724"
    mirrored = {
        "design_variables": [0.4, 0.6],
        "objective_gradient": [-0.3, -0.3],
        "constraint_jacobian": [[1.0, 1.0]],
        "lagrange_multipliers": [0.3],
        "constraint_values": [0.0],
        "hessian_vector_direction": [1.0, -0.5],
        "adjoint_hessian_vector_product": [2.0, -1.0],
        "finite_difference_hessian_vector_product": [2.000000001, -1.0],
        "hessian_vector_relative_tolerance": 1.0e-8,
        "kkt_stationarity_residual": [0.0, 0.0],
        "kkt_residual_tolerance": 1.0e-8,
        "optimization_model_sha256": "3" * 64,
    }
    payload[_OPTIMIZATION_KEY] = {
        "optimization_generation": generation,
        **{key: generation for key in (
            "design_generation", "gradient_generation", "constraint_generation",
            "adjoint_generation", "hessian_generation", "kkt_generation",
            "finite_difference_generation", "model_generation", "owner_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "model_owner": "optimization/model-724",
        "accepted_model_owner": "optimization/model-724",
        "optimization_result_sha256": "4" * 64,
        "accepted_optimization_result_sha256": "4" * 64,
    }
    return payload


def test_v40_public_positive_johnson_nedelec_and_adjoint_hessian_closure() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v40())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v40_public_johnson_nedelec_fembem_trace_normal_orientation_operator_energy_residual_mismatch() -> None:
    payload = _summary_v40()
    payload[_FEMBEM_KEY].update(
        {
            "trace_generation": "johnson-nedelec-723",
            "result_boundary_normals": [[-1.0, 0.0, 0.0]],
            "result_single_layer_matrix": [[-1.0]],
            "result_coupling_sign": 1.0,
            "result_energy_flux_residual_w": 9.0,
            "accepted_fembem_owner": "acoustic/old",
            "accepted_fembem_result_sha256": "a" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "johnson_nedelec_uses_current_trace_normals_operators_sign_residual_energy_mesh_owner_and_result"
    ]


def test_v40_public_adjoint_hessian_design_gradient_vector_product_kkt_fd_owner_mismatch() -> None:
    payload = _summary_v40()
    payload[_OPTIMIZATION_KEY].update(
        {
            "design_generation": "adjoint-hessian-723",
            "result_design_variables": [9.0],
            "result_objective_gradient": [9.0],
            "result_adjoint_hessian_vector_product": [9.0],
            "result_kkt_stationarity_residual": [9.0],
            "accepted_model_owner": "optimization/old",
            "accepted_optimization_result_sha256": "b" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "adjoint_hessian_uses_current_design_gradients_constraints_hvp_kkt_fd_model_owner_and_result"
    ]


def test_v40_public_rejects_self_consistent_energy_flux_gap() -> None:
    payload = _summary_v40()
    payload[_FEMBEM_KEY]["energy_flux_residual_w"] = 1.0
    payload[_FEMBEM_KEY]["result_energy_flux_residual_w"] = 1.0
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_negative_hessian_curvature() -> None:
    payload = _summary_v40()
    identity = payload[_OPTIMIZATION_KEY]
    identity["adjoint_hessian_vector_product"] = [-2.0, 1.0]
    identity["result_adjoint_hessian_vector_product"] = [-2.0, 1.0]
    identity["finite_difference_hessian_vector_product"] = [-2.000000001, 1.0]
    identity["result_finite_difference_hessian_vector_product"] = [-2.000000001, 1.0]
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"
