import copy
import json

import pytest

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from radia_mcp.radia_ngsolve.server import (
    regularized_trace_inverse_path_gate as mcp_gate,
)


def _summary():
    return {
        "schema": "regularized_trace_inverse_path/v1",
        "mesh": {
            "volume_element": "tetrahedron",
            "boundary_element": "triangle",
            "polynomial_order": 1,
            "volume_nodes": 5,
            "surface_nodes": 4,
            "tetrahedra": 1,
            "triangles": 4,
            "trace_rows": 4,
            "fem_unknowns": 5,
            "trace_nnz": 4,
        },
        "problem": {"noise_norm": 0.05},
        "path": {
            "alphas": [0.0, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0],
            "solution_norms": [10.0, 9.5, 8.0, 5.0, 2.0, 1.0],
            "trace_residual_norms": [0.0, 0.01, 0.03, 0.15, 0.8, 2.0],
            "weighted_trace_residuals": [0.0, 0.005, 0.02, 0.06, 0.3, 1.0],
            "normal_equation_residuals": [1.0e-13] * 6,
            "gradient_check_max_abs_errors": [2.0e-9] * 6,
        },
        "lcurve": {"selected_index": 3, "selected_alpha": 1.0e-2},
        "morozov": {"selected_index": 4, "selected_alpha": 1.0e-2},
        "crosscheck": {
            "reference_solver_count": 2,
            "max_solution_relative_error": 4.0e-14,
            "max_trace_relative_error": 5.0e-14,
            "max_regularized_objective_relative_error": 2.0e-14,
            "zero_alpha_objective_absolute_error": 3.0e-25,
        },
        "replay": {
            "count": 2,
            "selectors_identical": True,
            "max_relative_error": 0.0,
        },
    }


def test_regularized_trace_inverse_accepts_recomputed_choices_and_replay():
    result = regularized_trace_inverse_path_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert result["lcurve"]["selected_alpha"] == 1.0e-2
    assert result["morozov"]["selected_alpha"] == 1.0e-2


def test_regularized_trace_inverse_rejects_stale_choice_and_crosscheck():
    bad = copy.deepcopy(_summary())
    bad["lcurve"]["selected_index"] = 2
    bad["crosscheck"]["max_solution_relative_error"] = 1.0e-3
    bad["replay"]["max_relative_error"] = 1.0e-4
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["reported_lcurve_choice_matches"] is False
    assert result["checks"]["two_independent_linear_references_close"] is False
    assert result["checks"]["deterministic_replay_closes"] is False


def test_regularized_trace_inverse_rejects_nonmonotone_weighted_residual():
    bad = copy.deepcopy(_summary())
    bad["path"]["weighted_trace_residuals"][4] = 0.001
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["weighted_residual_increases_along_path"] is False
    assert result["checks"]["morozov_recomputation_passes"] is False


def test_regularized_trace_inverse_mcp_dispatches_and_rejects_bad_shape():
    result = json.loads(mcp_gate(json.dumps(_summary())))
    assert result["status"] == "ok"
    invalid = json.loads(mcp_gate('{"path": {}}'))
    assert invalid["status"] == "invalid_input"


def test_regularized_trace_inverse_rejects_gradient_and_replay_drift():
    bad = copy.deepcopy(_summary())
    bad["path"]["gradient_check_max_abs_errors"][3] = 1.0e-2
    bad["replay"]["max_relative_error"] = 0.1
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["finite_difference_gradients_close"] is False
    assert result["checks"]["deterministic_replay_closes"] is False


@pytest.mark.parametrize(
    "case_id",
    ["boundary_element", "alpha_order", "weighted_residual", "lcurve_choice", "reference_error"],
)
def test_counterfactual_curriculum90_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "boundary_element":
        bad["mesh"]["boundary_element"] = "quadrilateral"
    elif case_id == "alpha_order":
        bad["path"]["alphas"][3] = 5.0e-4
    elif case_id == "weighted_residual":
        bad["path"]["weighted_trace_residuals"][4] = 0.0
    elif case_id == "lcurve_choice":
        bad["lcurve"]["selected_index"] = 2
    else:
        bad["crosscheck"]["max_solution_relative_error"] = 1.0e-2
    assert regularized_trace_inverse_path_gate(bad)["status"] == "needs_attention"


def test_generalization_v3s_rejects_trace_sparsity_mismatch():
    bad = copy.deepcopy(_summary())
    bad["mesh"]["trace_nnz"] = 1
    assert regularized_trace_inverse_path_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v4_volume_element", "v4_polynomial_order", "v4_fem_unknown_count", "v4_normal_equation_residual", "v4_morozov_alpha"],
)
def test_counterfactual_curriculum90_v4_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "v4_volume_element":
        bad["mesh"]["volume_element"] = "hexahedron"
    elif case_id == "v4_polynomial_order":
        bad["mesh"]["polynomial_order"] = 2
    elif case_id == "v4_fem_unknown_count":
        bad["mesh"]["fem_unknowns"] = 4
    elif case_id == "v4_normal_equation_residual":
        bad["path"]["normal_equation_residuals"][2] = 1.0e-2
    else:
        bad["morozov"]["selected_alpha"] = 1.0
    assert regularized_trace_inverse_path_gate(bad)["status"] == "needs_attention"


def test_generalization_v5_rejects_surface_node_trace_mismatch():
    bad = copy.deepcopy(_summary())
    bad["mesh"]["surface_nodes"] = 3
    assert regularized_trace_inverse_path_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    ["v6_public_lcurve_index_alpha_mismatch", "v6_public_crosscheck_objective_drift"],
)
def test_generalization_v6_public(case_id):
    bad = copy.deepcopy(_summary())
    if case_id == "v6_public_lcurve_index_alpha_mismatch":
        bad["lcurve"]["selected_index"] += 1
    else:
        bad["crosscheck"]["max_regularized_objective_relative_error"] = 1.0e-2
    assert regularized_trace_inverse_path_gate(bad)["status"] == "needs_attention"


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_alpha_path_row_permutation",
        "v7_public_gradient_false_pass_tiny_step",
    ],
)
def test_generalization_v7_public(case_id):
    bad = copy.deepcopy(_summary())
    row_ids = [f"alpha-row-{index}" for index in range(6)]
    bad["path"].update(
        {
            "alpha_row_ids": row_ids,
            "solution_row_ids": row_ids.copy(),
            "residual_row_ids": row_ids.copy(),
            "gradient_check_step_sizes": [1.0e-6] * 6,
            "gradient_check_parameter_scales": [1.0] * 6,
            "gradient_check_objective_pair_deltas": [1.0e-6] * 6,
        }
    )
    if case_id == "v7_public_alpha_path_row_permutation":
        bad["path"]["solution_row_ids"][2:4] = reversed(
            bad["path"]["solution_row_ids"][2:4]
        )
        bad["path"]["residual_row_ids"][2:4] = reversed(
            bad["path"]["residual_row_ids"][2:4]
        )
    else:
        bad["path"]["gradient_check_step_sizes"][3] = 1.0e-320
        bad["path"]["gradient_check_objective_pair_deltas"][3] = 0.0
    assert regularized_trace_inverse_path_gate(bad)["status"] == "needs_attention"


def _with_v8_generations(summary):
    count = len(summary["path"]["alphas"])
    summary["path"].update(
        {
            "parameter_generation_ids": [
                f"parameter-{index}" for index in range(count)
            ],
            "gradient_parameter_generation_ids": [
                f"parameter-{index}" for index in range(count)
            ],
            "path_run_generation_ids": ["regularization-run-42"] * count,
            "solution_run_generation_ids": ["regularization-run-42"] * count,
        }
    )
    return summary


def _with_v9_bindings(summary):
    summary = _with_v8_generations(summary)
    count = len(summary["path"]["alphas"])
    summary["path"].update(
        {
            "objective_quadrature_generation_ids": [
                "boundary-quadrature-42" for _ in range(count)
            ],
            "gradient_quadrature_generation_ids": [
                "boundary-quadrature-42" for _ in range(count)
            ],
        }
    )
    summary["lcurve"]["curvature_parameterization"] = {
        "path_coordinate": "log10_alpha",
        "curvature_coordinate": "log10_alpha",
        "coordinate_transform_recorded": True,
    }
    return summary


def _with_v10_identity(summary):
    summary = _with_v9_bindings(summary)
    summary["design_variable_identity"] = {
        "design_variable_ids": ["thickness", "density", "damping"],
        "adjoint_gradient_design_variable_ids": [
            "thickness",
            "density",
            "damping",
        ],
        "finite_difference_design_variable_ids": [
            "thickness",
            "density",
            "damping",
        ],
        "design_generation": "design-42",
        "adjoint_gradient_design_generation": "design-42",
        "finite_difference_design_generation": "design-42",
    }
    summary["convolution_quadrature_identity"] = {
        "time_grid_step_s": 1.0e-4,
        "weight_generation_step_s": 1.0e-4,
        "time_grid_method": "BDF2",
        "weight_generation_method": "BDF2",
        "time_grid_generation": "time-grid-42",
        "weight_time_grid_generation": "time-grid-42",
    }
    return summary


def _with_v11_identity(summary):
    summary = _with_v10_identity(summary)
    summary["cq_inverse_laplace_contour_identity"] = {
        "contour_generation": "cq-contour-43",
        "transfer_sample_contour_generation": "cq-contour-43",
        "laplace_sample_ids": ["s-0", "s-1", "s-2", "s-3"],
        "transfer_sample_ids": ["s-0", "s-1", "s-2", "s-3"],
        "sqrt_branch_conventions": ["principal_outgoing"] * 4,
        "inverse_laplace_branch_convention": "principal_outgoing",
    }
    mesh_digest = "7" * 64
    summary["fembem_trace_orientation_identity"] = {
        "volume_mesh_sha256": mesh_digest,
        "trace_mesh_sha256": mesh_digest,
        "volume_mesh_generation": "vol-mesh-43",
        "trace_mesh_generation": "vol-mesh-43",
        "boundary_orientation_mesh_generation": "vol-mesh-43",
        "trace_boundary_triangle_digest": "trace-triangles-43",
        "oriented_boundary_triangle_digest": "trace-triangles-43",
        "outward_orientation_verified": True,
    }
    return summary


def test_accepts_v8_parameter_and_regularization_run_generations():
    result = regularized_trace_inverse_path_gate(_with_v8_generations(_summary()))
    assert result["status"] == "ok"


def test_v8_public_gradient_previous_parameter_generation():
    bad = _with_v8_generations(_summary())
    bad["path"]["gradient_parameter_generation_ids"][3] = "parameter-2"
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["gradient_uses_current_parameter_generation"] is False


def test_v8_public_regularization_restart_row_reuse():
    bad = _with_v8_generations(_summary())
    bad["path"]["solution_run_generation_ids"][3] = "regularization-run-41"
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["solution_rows_share_regularization_run_generation"] is False


def test_v9_public_gradient_quadrature_generation_mismatch():
    bad = _with_v9_bindings(_summary())
    bad["path"]["gradient_quadrature_generation_ids"][3] = (
        "boundary-quadrature-41"
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["gradient_uses_current_boundary_quadrature_generation"]
        is False
    )


def test_v9_public_regularization_curvature_parameterization_mismatch():
    bad = _with_v9_bindings(_summary())
    bad["lcurve"]["curvature_parameterization"].update(
        {
            "curvature_coordinate": "natural_log_alpha",
            "coordinate_transform_recorded": False,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["lcurve_curvature_uses_recorded_path_parameterization"]
        is False
    )


def test_v10_public_adjoint_gradient_design_order_mismatch():
    bad = _with_v10_identity(_summary())
    bad["design_variable_identity"]["finite_difference_design_variable_ids"] = [
        "density",
        "thickness",
        "damping",
    ]
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "adjoint_and_finite_difference_share_design_variable_order"
        ]
        is False
    )


def test_v10_public_cq_weights_time_step_method_mismatch():
    bad = _with_v10_identity(_summary())
    bad["convolution_quadrature_identity"].update(
        {
            "weight_generation_step_s": 2.0e-4,
            "weight_generation_method": "BDF1",
            "weight_time_grid_generation": "time-grid-41",
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["cq_weights_match_current_time_grid_and_method"] is False


def test_v11_public_cq_inverse_laplace_contour_branch_mismatch():
    bad = _with_v11_identity(_summary())
    bad["cq_inverse_laplace_contour_identity"]["sqrt_branch_conventions"][
        2
    ] = "opposite_incoming"
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "cq_transfer_samples_share_inverse_laplace_contour_branch"
        ]
        is False
    )


def test_v11_public_fembem_trace_orientation_mesh_mismatch():
    bad = _with_v11_identity(_summary())
    bad["fembem_trace_orientation_identity"].update(
        {
            "boundary_orientation_mesh_generation": "vol-mesh-42",
            "oriented_boundary_triangle_digest": "trace-triangles-42",
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "fembem_trace_orientation_matches_current_volume_mesh"
        ]
        is False
    )
