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
