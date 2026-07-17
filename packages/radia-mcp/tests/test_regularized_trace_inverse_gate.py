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


def _with_v12_identity(summary):
    summary = _with_v11_identity(summary)
    summary["cq_starting_weights_identity"] = {
        "multistep_method": "BDF2",
        "multistep_order": 2,
        "convolution_weight_method": "BDF2",
        "startup_correction_method": "BDF2_starting_correction",
        "startup_correction_order": 2,
        "weight_generation": "cq-weights-14",
        "startup_correction_weight_generation": "cq-weights-14",
    }
    summary["hmatrix_cluster_permutation_identity"] = {
        "boundary_mesh_sha256": "8" * 64,
        "boundary_triangle_ordering_sha256": "9" * 64,
        "cluster_permutation_triangle_ordering_sha256": "9" * 64,
        "boundary_mesh_generation": "boundary-mesh-14",
        "cluster_permutation_mesh_generation": "boundary-mesh-14",
        "cluster_permutation_generation": "cluster-permutation-14",
        "hmatrix_assembly_permutation_generation": "cluster-permutation-14",
    }
    return summary


def _with_v13_identity(summary):
    summary = _with_v12_identity(summary)
    summary["fembem_trace_surface_normal_orientation_identity"] = {
        "boundary_mesh_generation": "boundary-mesh-15",
        "fem_trace_boundary_mesh_generation": "boundary-mesh-15",
        "bem_surface_mesh_generation": "boundary-mesh-15",
        "fem_outward_normal_generation": "surface-normal-15",
        "bem_normal_generation": "surface-normal-15",
        "trace_operator_normal_generation": "surface-normal-15",
        "fem_normal_orientation": "outward",
        "bem_normal_orientation": "outward",
        "trace_normal_sign": 1,
        "fem_boundary_triangle_sha256": "b" * 64,
        "bem_boundary_triangle_sha256": "b" * 64,
    }
    summary["cq_inverse_z_transform_fft_normalization_identity"] = {
        "frequency_sample_count": 128,
        "time_sample_count": 128,
        "forward_fft_normalization": "unscaled",
        "inverse_fft_normalization": "one_over_n",
        "inverse_fft_scale": 1.0 / 128.0,
        "cq_frequency_generation": "cq-frequency-15",
        "inverse_transform_frequency_generation": "cq-frequency-15",
        "transform_convention": (
            "forward_negative_exponent_inverse_positive_exponent"
        ),
        "reconstruction_transform_convention": (
            "forward_negative_exponent_inverse_positive_exponent"
        ),
    }
    return summary


def _with_v14_identity(summary):
    summary = _with_v13_identity(summary)
    summary["fembem_sesquilinear_inner_product_identity"] = {
        "coupling_generation": "fembem-coupling-16",
        "fem_operator_generation": "fembem-coupling-16",
        "bem_operator_generation": "fembem-coupling-16",
        "trace_operator_generation": "fembem-coupling-16",
        "inner_product_convention": "conjugate_test_linear_trial",
        "fem_inner_product_convention": "conjugate_test_linear_trial",
        "bem_inner_product_convention": "conjugate_test_linear_trial",
        "trace_inner_product_convention": "conjugate_test_linear_trial",
        "test_argument_conjugated": True,
        "trial_argument_conjugated": False,
        "operator_metadata_sha256": "1" * 64,
        "assembled_operator_metadata_sha256": "1" * 64,
    }
    summary["hmatrix_low_rank_tolerance_norm_basis_identity"] = {
        "hmatrix_assembly_generation": "hmatrix-assembly-16",
        "low_rank_factor_generation": "hmatrix-assembly-16",
        "compression_tolerance": 1.0e-5,
        "measured_relative_error": 8.0e-6,
        "tolerance_norm_basis": "relative_frobenius",
        "measured_error_norm_basis": "relative_frobenius",
        "acceptance_norm_basis": "relative_frobenius",
        "tolerance_calibration_generation": "hmatrix-tolerance-16",
        "acceptance_tolerance_generation": "hmatrix-tolerance-16",
        "dense_block_sha256": "2" * 64,
        "accepted_block_source_sha256": "2" * 64,
    }
    summary["complex_operator_symmetry_residual_norm_identity"] = {
        "operator_generation": "complex-operator-17",
        "assembled_operator_generation": "complex-operator-17",
        "residual_operator_generation": "complex-operator-17",
        "operator_symmetry_class": "complex_symmetric",
        "assembly_symmetry_check": "transpose",
        "residual_symmetry_check": "transpose",
        "hermitian_symmetry_check_applied": False,
        "residual_norm_basis": "complex_euclidean",
        "reported_residual_norm_basis": "complex_euclidean",
        "residual_norm": 8.0e-8,
        "reported_residual_norm": 8.0e-8,
        "operator_metadata_sha256": "1" * 64,
        "residual_operator_metadata_sha256": "1" * 64,
    }
    summary["hmatrix_admissibility_cluster_diameter_metric_identity"] = {
        "cluster_tree_generation": "cluster-tree-17",
        "source_cluster_generation": "cluster-tree-17",
        "target_cluster_generation": "cluster-tree-17",
        "diameter_metric": "euclidean_l2",
        "admissibility_diameter_metric": "euclidean_l2",
        "threshold_calibration_diameter_metric": "euclidean_l2",
        "source_cluster_diameter": 0.2,
        "target_cluster_diameter": 0.1,
        "cluster_separation": 1.0,
        "admissibility_eta": 0.5,
        "admissible": True,
        "threshold_calibration_generation": "admissibility-17",
        "acceptance_threshold_generation": "admissibility-17",
        "cluster_geometry_sha256": "2" * 64,
        "admissibility_geometry_sha256": "2" * 64,
    }
    return summary


def _with_v16_identity(summary):
    summary = _with_v14_identity(summary)
    summary["cq_convolution_weight_ztransform_branch_order_identity"] = {
        "cq_scheme_generation": "cq-scheme-18",
        "convolution_weight_generation": "cq-scheme-18",
        "transfer_sample_generation": "cq-scheme-18",
        "z_transform_branch": "principal",
        "convolution_weight_z_transform_branch": "principal",
        "transfer_sample_z_transform_branch": "principal",
        "multistep_method": "bdf2",
        "convolution_weight_multistep_method": "bdf2",
        "transfer_sample_multistep_method": "bdf2",
        "multistep_order": 2,
        "convolution_weight_multistep_order": 2,
        "transfer_sample_multistep_order": 2,
        "scheme_metadata_sha256": "1" * 64,
        "convolution_weight_scheme_metadata_sha256": "1" * 64,
        "transfer_sample_scheme_metadata_sha256": "1" * 64,
    }
    summary["fembem_trace_matrix_boundary_node_generation_identity"] = {
        "volume_mesh_generation": "vol-mesh-18",
        "boundary_node_generation": "boundary-nodes-18",
        "fem_operator_boundary_node_generation": "boundary-nodes-18",
        "bem_operator_boundary_node_generation": "boundary-nodes-18",
        "trace_matrix_boundary_node_generation": "boundary-nodes-18",
        "boundary_node_count": 42,
        "trace_matrix_row_count": 42,
        "trace_matrix_nonzero_count": 42,
        "boundary_node_ids_sha256": "2" * 64,
        "fem_boundary_node_ids_sha256": "2" * 64,
        "bem_boundary_node_ids_sha256": "2" * 64,
        "trace_matrix_boundary_node_ids_sha256": "2" * 64,
    }
    return summary


def _with_v17_identity(summary):
    summary = _with_v16_identity(summary)
    summary["cq_frequency_contour_radius_damping_generation_identity"] = {
        "frequency_sampling_generation": "cq-frequency-sampling-19",
        "transfer_sample_frequency_generation": "cq-frequency-sampling-19",
        "contour_metadata_frequency_generation": "cq-frequency-sampling-19",
        "inverse_transform_frequency_generation": "cq-frequency-sampling-19",
        "time_step_s": 2.5e-4,
        "transfer_sample_time_step_s": 2.5e-4,
        "sample_count": 512,
        "transfer_sample_count": 512,
        "contour_radius": 0.985,
        "inverse_transform_contour_radius": 0.985,
        "damping_factor": 0.97,
        "inverse_transform_damping_factor": 0.97,
        "contour_metadata_sha256": "5" * 64,
        "inverse_transform_contour_metadata_sha256": "5" * 64,
    }
    summary["hmatrix_aca_pivot_cluster_permutation_generation_identity"] = {
        "cluster_tree_generation": "cluster-tree-19",
        "source_cluster_permutation_generation": "cluster-permutation-19",
        "target_cluster_permutation_generation": "cluster-permutation-19",
        "aca_source_cluster_permutation_generation": "cluster-permutation-19",
        "aca_target_cluster_permutation_generation": "cluster-permutation-19",
        "row_pivot_global_ids": [11, 17, 23],
        "aca_row_pivot_global_ids": [11, 17, 23],
        "column_pivot_global_ids": [5, 13, 29],
        "aca_column_pivot_global_ids": [5, 13, 29],
        "cluster_permutation_sha256": "6" * 64,
        "aca_cluster_permutation_sha256": "6" * 64,
    }
    return summary


def _with_v18_identity(summary):
    summary = _with_v17_identity(summary)
    summary[
        "cq_inverse_transform_conjugate_frequency_bin_ownership_identity"
    ] = {
        "transfer_sample_generation": "cq-transfer-20",
        "frequency_bin_transfer_sample_generation": "cq-transfer-20",
        "inverse_transform_transfer_sample_generation": "cq-transfer-20",
        "sample_count": 512,
        "inverse_transform_sample_count": 512,
        "positive_frequency_bin_indices": [1, 2, 3],
        "conjugate_frequency_bin_indices": [511, 510, 509],
        "inverse_transform_conjugate_frequency_bin_indices": [511, 510, 509],
        "dc_bin_index": 0,
        "nyquist_bin_index": 256,
        "real_response_conjugate_symmetry": True,
        "inverse_transform_conjugate_symmetry": True,
        "frequency_bin_map_sha256": "5" * 64,
        "inverse_transform_frequency_bin_map_sha256": "5" * 64,
    }
    summary[
        "hmatrix_block_cluster_bounding_box_mesh_scale_generation_identity"
    ] = {
        "mesh_generation": "mesh-20",
        "cluster_tree_mesh_generation": "mesh-20",
        "block_cluster_mesh_generation": "mesh-20",
        "mesh_scale_generation": "mesh-scale-20",
        "cluster_bounding_box_mesh_scale_generation": "mesh-scale-20",
        "block_admissibility_mesh_scale_generation": "mesh-scale-20",
        "mesh_length_scale_to_m": 1.0e-3,
        "cluster_bounding_box_length_scale_to_m": 1.0e-3,
        "source_cluster_bounding_box_m": [0.0, 0.0, 0.0, 0.01, 0.02, 0.03],
        "block_source_cluster_bounding_box_m": [
            0.0,
            0.0,
            0.0,
            0.01,
            0.02,
            0.03,
        ],
        "cluster_bounding_box_sha256": "6" * 64,
        "block_admissibility_bounding_box_sha256": "6" * 64,
    }
    return summary


def _with_v19_identity(summary):
    summary = _with_v18_identity(summary)
    contour_nodes = [[0.95, 0.0], [0.0, 0.95], [-0.95, 0.0], [0.0, -0.95]]
    summary["cq_laplace_contour_fft_scaling_timestep_generation_identity"] = {
        "cq_generation": "cq-21",
        "laplace_contour_cq_generation": "cq-21",
        "fft_scaling_cq_generation": "cq-21",
        "timestep_cq_generation": "cq-21",
        "inverse_transform_cq_generation": "cq-21",
        "sample_count": 512,
        "fft_sample_count": 512,
        "time_step_s": 1.0e-5,
        "inverse_transform_time_step_s": 1.0e-5,
        "fft_scaling": 1.0 / 512.0,
        "inverse_transform_fft_scaling": 1.0 / 512.0,
        "laplace_contour_nodes_re_im": contour_nodes,
        "inverse_transform_laplace_contour_nodes_re_im": contour_nodes,
        "contour_fft_timestep_sha256": "a" * 64,
        "inverse_transform_contour_fft_timestep_sha256": "a" * 64,
    }
    summary["hmatrix_aca_tolerance_norm_block_rank_generation_identity"] = {
        "block_generation": "hmatrix-block-21",
        "aca_tolerance_block_generation": "hmatrix-block-21",
        "norm_estimate_block_generation": "hmatrix-block-21",
        "retained_rank_block_generation": "hmatrix-block-21",
        "assembled_block_generation": "hmatrix-block-21",
        "block_ids": [101, 102, 103],
        "aca_block_ids": [101, 102, 103],
        "aca_tolerances": [1.0e-6, 1.0e-6, 1.0e-7],
        "applied_aca_tolerances": [1.0e-6, 1.0e-6, 1.0e-7],
        "block_norm_estimates": [12.0, 8.0, 3.0],
        "aca_block_norm_estimates": [12.0, 8.0, 3.0],
        "retained_ranks": [8, 6, 4],
        "assembled_block_ranks": [8, 6, 4],
        "aca_block_table_sha256": "b" * 64,
        "assembled_aca_block_table_sha256": "b" * 64,
    }
    return summary


def _with_v20_identity(summary):
    summary = _with_v19_identity(summary)
    summary["fembem_trace_orientation_normal_node_map_generation_identity"] = {
        "volume_mesh_generation": "vol-mesh-22",
        "boundary_trace_volume_mesh_generation": "vol-mesh-22",
        "boundary_mesh_generation": "boundary-mesh-22",
        "triangle_orientation_boundary_mesh_generation": "boundary-mesh-22",
        "normal_boundary_mesh_generation": "boundary-mesh-22",
        "node_map_boundary_mesh_generation": "boundary-mesh-22",
        "fem_trace_boundary_mesh_generation": "boundary-mesh-22",
        "bem_trace_boundary_mesh_generation": "boundary-mesh-22",
        "boundary_triangle_ids": [101, 102, 103],
        "trace_triangle_ids": [101, 102, 103],
        "triangle_orientations": [1, -1, 1],
        "applied_triangle_orientations": [1, -1, 1],
        "normal_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "applied_normal_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "volume_boundary_node_ids": [11, 12, 13, 14],
        "trace_volume_node_ids": [11, 12, 13, 14],
        "trace_map_sha256": "4" * 64,
        "assembled_trace_map_sha256": "4" * 64,
    }
    summary["cq_causality_conjugate_symmetry_contour_pair_generation_identity"] = {
        "cq_generation": "cq-22",
        "frequency_pair_cq_generation": "cq-22",
        "contour_index_cq_generation": "cq-22",
        "causality_window_cq_generation": "cq-22",
        "inverse_transform_cq_generation": "cq-22",
        "positive_frequency_ids": [1, 2, 3],
        "conjugate_frequency_ids": [511, 510, 509],
        "inverse_transform_conjugate_frequency_ids": [511, 510, 509],
        "contour_indices": [1, 2, 3, 509, 510, 511],
        "inverse_transform_contour_indices": [1, 2, 3, 509, 510, 511],
        "causality_window_samples": [0, 511],
        "inverse_transform_causality_window_samples": [0, 511],
        "real_time_response": True,
        "inverse_transform_conjugate_symmetry": True,
        "precausal_max_abs": 1.0e-12,
        "precausal_tolerance": 1.0e-10,
        "cq_pair_table_sha256": "5" * 64,
        "inverse_transform_pair_table_sha256": "5" * 64,
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


def test_v12_public_cq_starting_weights_multistep_startup_mismatch():
    bad = _with_v12_identity(_summary())
    bad["cq_starting_weights_identity"].update(
        {
            "startup_correction_method": "BDF1_starting_correction",
            "startup_correction_order": 1,
            "startup_correction_weight_generation": "cq-weights-13",
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["cq_starting_weights_match_multistep_startup_scheme"]
        is False
    )


def test_v12_public_hmatrix_cluster_permutation_boundary_mesh_mismatch():
    bad = _with_v12_identity(_summary())
    bad["hmatrix_cluster_permutation_identity"].update(
        {
            "cluster_permutation_triangle_ordering_sha256": "a" * 64,
            "cluster_permutation_mesh_generation": "boundary-mesh-13",
            "hmatrix_assembly_permutation_generation": "cluster-permutation-13",
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "hmatrix_cluster_permutation_matches_boundary_triangle_order"
        ]
        is False
    )


def test_v13_public_fembem_trace_surface_normal_orientation_mismatch():
    bad = _with_v13_identity(_summary())
    bad["fembem_trace_surface_normal_orientation_identity"].update(
        {
            "bem_normal_generation": "surface-normal-14",
            "bem_normal_orientation": "inward",
            "trace_normal_sign": -1,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "fembem_trace_and_surface_share_outward_normal_orientation"
        ]
        is False
    )


def test_v13_public_cq_inverse_z_transform_fft_normalization_mismatch():
    bad = _with_v13_identity(_summary())
    bad["cq_inverse_z_transform_fft_normalization_identity"].update(
        {
            "inverse_fft_normalization": "unscaled",
            "inverse_fft_scale": 1.0,
            "reconstruction_transform_convention": "legacy_unscaled_inverse",
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "cq_inverse_z_transform_uses_matching_fft_normalization"
        ]
        is False
    )


def test_accepts_v14_sesquilinear_and_hmatrix_norm_lineage():
    result = regularized_trace_inverse_path_gate(_with_v14_identity(_summary()))
    assert result["status"] == "ok"


def test_v14_public_fembem_sesquilinear_inner_product_conjugation_mismatch():
    bad = _with_v14_identity(_summary())
    bad["fembem_sesquilinear_inner_product_identity"].update(
        {
            "bem_operator_generation": "fembem-coupling-15",
            "bem_inner_product_convention": "linear_test_conjugate_trial",
            "test_argument_conjugated": False,
            "trial_argument_conjugated": True,
            "assembled_operator_metadata_sha256": "5" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "fembem_operators_share_sesquilinear_conjugation_convention"
        ]
        is False
    )


def test_v14_public_hmatrix_low_rank_tolerance_norm_basis_mismatch():
    bad = _with_v14_identity(_summary())
    bad["hmatrix_low_rank_tolerance_norm_basis_identity"].update(
        {
            "measured_error_norm_basis": "relative_spectral",
            "acceptance_norm_basis": "relative_spectral",
            "acceptance_tolerance_generation": "hmatrix-tolerance-15",
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "hmatrix_low_rank_acceptance_uses_calibrated_norm_basis"
        ]
        is False
    )


def test_v15_public_complex_symmetric_hermitian_residual_norm_mismatch():
    bad = _with_v14_identity(_summary())
    bad["complex_operator_symmetry_residual_norm_identity"].update(
        {
            "residual_symmetry_check": "conjugate_transpose",
            "hermitian_symmetry_check_applied": True,
            "residual_operator_metadata_sha256": "5" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "complex_operator_residual_uses_transpose_symmetry_class"
        ]
        is False
    )


def test_v15_public_hmatrix_admissibility_cluster_diameter_metric_mismatch():
    bad = _with_v14_identity(_summary())
    bad["hmatrix_admissibility_cluster_diameter_metric_identity"].update(
        {
            "admissibility_diameter_metric": "infinity_linf",
            "acceptance_threshold_generation": "admissibility-16",
            "admissibility_geometry_sha256": "5" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "hmatrix_admissibility_uses_one_cluster_diameter_metric"
        ]
        is False
    )


def test_accepts_v16_cq_weight_and_trace_boundary_node_lineage():
    result = regularized_trace_inverse_path_gate(_with_v16_identity(_summary()))
    assert result["status"] == "ok"


def test_v16_public_cq_convolution_weight_ztransform_branch_order_mismatch():
    bad = _with_v16_identity(_summary())
    bad["cq_convolution_weight_ztransform_branch_order_identity"].update(
        {
            "convolution_weight_generation": "cq-scheme-17",
            "convolution_weight_z_transform_branch": "negative_real_cut",
            "convolution_weight_multistep_method": "bdf1",
            "convolution_weight_multistep_order": 1,
            "convolution_weight_scheme_metadata_sha256": "5" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "cq_convolution_weights_share_ztransform_branch_and_multistep_order"
        ]
        is False
    )


def test_v16_public_fembem_trace_matrix_boundary_node_generation_mismatch():
    bad = _with_v16_identity(_summary())
    bad["fembem_trace_matrix_boundary_node_generation_identity"].update(
        {
            "trace_matrix_boundary_node_generation": "boundary-nodes-17",
            "trace_matrix_row_count": 41,
            "trace_matrix_nonzero_count": 41,
            "trace_matrix_boundary_node_ids_sha256": "5" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "fembem_trace_matrix_uses_current_boundary_node_generation"
        ]
        is False
    )


def test_accepts_v17_cq_contour_and_hmatrix_aca_pivot_lineage():
    result = regularized_trace_inverse_path_gate(_with_v17_identity(_summary()))
    assert result["status"] == "ok"


def test_v17_public_cq_frequency_contour_radius_damping_generation_mismatch():
    bad = _with_v17_identity(_summary())
    bad["cq_frequency_contour_radius_damping_generation_identity"].update(
        {
            "contour_metadata_frequency_generation": "cq-frequency-sampling-18",
            "inverse_transform_frequency_generation": "cq-frequency-sampling-18",
            "inverse_transform_contour_radius": 0.96,
            "inverse_transform_damping_factor": 0.92,
            "inverse_transform_contour_metadata_sha256": "9" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "cq_inverse_transform_uses_current_frequency_contour_metadata"
        ]
        is False
    )


def test_v17_public_hmatrix_aca_pivot_cluster_permutation_generation_mismatch():
    bad = _with_v17_identity(_summary())
    bad["hmatrix_aca_pivot_cluster_permutation_generation_identity"].update(
        {
            "aca_source_cluster_permutation_generation": "cluster-permutation-18",
            "aca_target_cluster_permutation_generation": "cluster-permutation-18",
            "aca_row_pivot_global_ids": [17, 11, 23],
            "aca_column_pivot_global_ids": [13, 5, 29],
            "aca_cluster_permutation_sha256": "9" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "hmatrix_aca_pivots_use_current_cluster_permutation"
        ]
        is False
    )


def test_accepts_v18_cq_conjugate_bins_and_hmatrix_mesh_scale_lineage():
    result = regularized_trace_inverse_path_gate(_with_v18_identity(_summary()))
    assert result["status"] == "ok"


def test_v18_public_cq_inverse_transform_conjugate_frequency_bin_ownership_generation_mismatch():
    bad = _with_v18_identity(_summary())
    bad[
        "cq_inverse_transform_conjugate_frequency_bin_ownership_identity"
    ].update(
        {
            "inverse_transform_transfer_sample_generation": "cq-transfer-19",
            "inverse_transform_sample_count": 510,
            "inverse_transform_conjugate_frequency_bin_indices": [509, 508, 507],
            "inverse_transform_conjugate_symmetry": False,
            "inverse_transform_frequency_bin_map_sha256": "9" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "cq_inverse_transform_uses_current_conjugate_frequency_bins"
        ]
        is False
    )


def test_v18_public_hmatrix_block_cluster_bounding_box_mesh_scale_generation_mismatch():
    bad = _with_v18_identity(_summary())
    bad[
        "hmatrix_block_cluster_bounding_box_mesh_scale_generation_identity"
    ].update(
        {
            "cluster_bounding_box_mesh_scale_generation": "mesh-scale-19",
            "block_admissibility_mesh_scale_generation": "mesh-scale-19",
            "cluster_bounding_box_length_scale_to_m": 1.0,
            "block_source_cluster_bounding_box_m": [
                0.0,
                0.0,
                0.0,
                10.0,
                20.0,
                30.0,
            ],
            "block_admissibility_bounding_box_sha256": "9" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "hmatrix_block_clusters_use_current_mesh_scale_bounding_boxes"
        ]
        is False
    )


def test_accepts_v19_cq_contour_fft_and_hmatrix_aca_lineage():
    result = regularized_trace_inverse_path_gate(_with_v19_identity(_summary()))
    assert result["status"] == "ok"


def test_v19_public_cq_laplace_contour_fft_scaling_timestep_generation_mismatch():
    bad = _with_v19_identity(_summary())
    bad["cq_laplace_contour_fft_scaling_timestep_generation_identity"].update(
        {
            "fft_scaling_cq_generation": "cq-20",
            "timestep_cq_generation": "cq-20",
            "inverse_transform_time_step_s": 2.0e-5,
            "inverse_transform_fft_scaling": 1.0 / 256.0,
            "inverse_transform_laplace_contour_nodes_re_im": [
                [0.9, 0.0],
                [0.0, 0.9],
                [-0.9, 0.0],
                [0.0, -0.9],
            ],
            "inverse_transform_contour_fft_timestep_sha256": "9" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "cq_inverse_transform_uses_current_contour_fft_scaling_and_timestep"
        ]
        is False
    )


def test_v19_public_hmatrix_aca_tolerance_norm_block_rank_generation_mismatch():
    bad = _with_v19_identity(_summary())
    bad["hmatrix_aca_tolerance_norm_block_rank_generation_identity"].update(
        {
            "aca_tolerance_block_generation": "hmatrix-block-20",
            "norm_estimate_block_generation": "hmatrix-block-20",
            "assembled_block_generation": "hmatrix-block-20",
            "applied_aca_tolerances": [1.0e-3, 1.0e-3, 1.0e-3],
            "aca_block_norm_estimates": [3.0, 8.0, 12.0],
            "assembled_block_ranks": [4, 6, 8],
            "assembled_aca_block_table_sha256": "9" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "hmatrix_aca_uses_current_block_tolerance_norm_and_rank"
        ]
        is False
    )


def test_accepts_v20_trace_map_and_cq_causality_lineage():
    result = regularized_trace_inverse_path_gate(_with_v20_identity(_summary()))
    assert result["status"] == "ok"


def test_v20_public_fembem_trace_orientation_normal_node_map_generation_mismatch():
    bad = _with_v20_identity(_summary())
    bad["fembem_trace_orientation_normal_node_map_generation_identity"].update(
        {
            "triangle_orientation_boundary_mesh_generation": "boundary-mesh-21",
            "normal_boundary_mesh_generation": "boundary-mesh-21",
            "trace_triangle_ids": [103, 102, 101],
            "applied_triangle_orientations": [-1, 1, -1],
            "applied_normal_sha256": ["3" * 64, "2" * 64, "1" * 64],
            "trace_volume_node_ids": [14, 13, 12, 11],
            "assembled_trace_map_sha256": "f" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fembem_trace_uses_current_orientation_normals_and_volume_node_map"
    ]


def test_v20_public_cq_causality_conjugate_symmetry_contour_pair_generation_mismatch():
    bad = _with_v20_identity(_summary())
    bad["cq_causality_conjugate_symmetry_contour_pair_generation_identity"].update(
        {
            "frequency_pair_cq_generation": "cq-21",
            "causality_window_cq_generation": "cq-21",
            "inverse_transform_conjugate_frequency_ids": [509, 510, 511],
            "inverse_transform_contour_indices": [1, 2, 3, 511, 510, 509],
            "inverse_transform_causality_window_samples": [16, 511],
            "inverse_transform_conjugate_symmetry": False,
            "precausal_max_abs": 1.0e-4,
            "inverse_transform_pair_table_sha256": "f" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(bad)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_inverse_transform_uses_current_causal_conjugate_contour_pairs"
    ]
