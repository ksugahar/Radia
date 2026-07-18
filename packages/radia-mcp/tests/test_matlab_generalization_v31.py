from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v30 import _summary_v30

_PROMOTED_CASE_IDS = (
    "v31_public_complex_automatic_differentiation_wirtinger_branch_gradient_finite_difference_mismatch",
    "v31_public_pde_quadratic_curved_mesh_to_vol_midnode_boundary_face_orientation_mismatch",
)


def _summary_v31():
    payload = deepcopy(_summary_v30())
    generation = "complex-ad-351"
    payload["complex_ad_wirtinger_conjugation_branch_scaling_fd_mesh_result_identity"] = {
        "ad_generation": generation,
        "wirtinger_ad_generation": generation,
        "conjugation_ad_generation": generation,
        "branch_ad_generation": generation,
        "scaling_ad_generation": generation,
        "finite_difference_ad_generation": generation,
        "mesh_ad_generation": generation,
        "result_ad_generation": generation,
        "wirtinger_convention": "dJ_dconj_z",
        "result_wirtinger_convention": "dJ_dconj_z",
        "adjoint_conjugation": "conjugate_transpose",
        "result_adjoint_conjugation": "conjugate_transpose",
        "objective_branch": "real_objective",
        "result_objective_branch": "real_objective",
        "design_variable_scaling": [1.0, 0.1, 10.0],
        "result_design_variable_scaling": [1.0, 0.1, 10.0],
        "gradient_ri": [[0.2, -0.1], [0.05, 0.03], [-0.4, 0.2]],
        "finite_difference_gradient_ri": [
            [0.20000001, -0.10000001],
            [0.05000001, 0.02999999],
            [-0.39999999, 0.20000001],
        ],
        "finite_difference_relative_error": 5.0e-8,
        "finite_difference_tolerance": 1.0e-6,
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "pde-p2-vol-351"
    payload["pde_quadratic_curved_vol_midnode_tet_boundary_region_order_mesh_identity"] = {
        "mesh_generation": generation,
        "midnode_mesh_generation": generation,
        "tet_mesh_generation": generation,
        "boundary_mesh_generation": generation,
        "region_mesh_generation": generation,
        "order_mesh_generation": generation,
        "result_mesh_generation": generation,
        "geometry_order": 2,
        "result_geometry_order": 2,
        "tet_connectivity": [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]],
        "result_tet_connectivity": [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]],
        "curved_midnode_sha256": "3" * 64,
        "result_curved_midnode_sha256": "3" * 64,
        "boundary_tri_connectivity": [[1, 3, 2, 7, 6, 5]],
        "result_boundary_tri_connectivity": [[1, 3, 2, 7, 6, 5]],
        "boundary_orientation": [1],
        "result_boundary_orientation": [1],
        "tet_region_labels": [11],
        "result_tet_region_labels": [11],
        "boundary_region_labels": [21],
        "result_boundary_region_labels": [21],
        "mesh_sha256": "4" * 64,
        "result_mesh_sha256": "4" * 64,
    }
    return payload


def test_v31_public_positive_complex_ad_and_pde_p2_vol_identities() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v31())["status"] == "ok"


def test_v31_public_complex_automatic_differentiation_wirtinger_branch_gradient_finite_difference_mismatch() -> None:
    payload = _summary_v31()
    identity = payload[
        "complex_ad_wirtinger_conjugation_branch_scaling_fd_mesh_result_identity"
    ]
    identity.update(
        {
            "wirtinger_ad_generation": "complex-ad-350",
            "finite_difference_ad_generation": "complex-ad-349",
            "result_wirtinger_convention": "dJ_dz",
            "result_adjoint_conjugation": "transpose_without_conjugation",
            "result_objective_branch": "imaginary_branch",
            "result_design_variable_scaling": [1.0, 1.0, 1.0],
            "finite_difference_gradient_ri": [[-0.2, 0.1]],
            "finite_difference_relative_error": 0.5,
            "result_mesh_sha256": "9" * 64,
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "complex_ad_uses_current_wirtinger_conjugation_branch_scaling_fd_mesh_and_result"
    ]


def test_v31_public_pde_quadratic_curved_mesh_to_vol_midnode_boundary_face_orientation_mismatch() -> None:
    payload = _summary_v31()
    identity = payload[
        "pde_quadratic_curved_vol_midnode_tet_boundary_region_order_mesh_identity"
    ]
    identity.update(
        {
            "midnode_mesh_generation": "pde-p2-vol-350",
            "boundary_mesh_generation": "pde-p2-vol-349",
            "result_geometry_order": 1,
            "result_tet_connectivity": [[1, 3, 2, 4, 5, 7, 6, 8, 10, 9]],
            "result_curved_midnode_sha256": "b" * 64,
            "result_boundary_tri_connectivity": [[1, 2, 3, 5, 6, 7]],
            "result_boundary_orientation": [-1],
            "result_tet_region_labels": [12],
            "result_boundary_region_labels": [22],
            "result_mesh_sha256": "c" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pde_quadratic_vol_uses_current_midnodes_tets_boundary_orientation_regions_order_and_mesh"
    ]
