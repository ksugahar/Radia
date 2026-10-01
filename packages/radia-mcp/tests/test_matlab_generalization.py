from __future__ import annotations

import json
from copy import deepcopy

from radia_mcp.matlab_agentic_ml import (
    validate_matlab_ml_rl_v44_identity,
    validate_matlab_ml_rl_v45_identity,
    validate_matlab_ml_rl_v47_identity,
    validate_matlab_ml_rl_v48_identity,
    validate_matlab_ml_rl_v49_identity,
    validate_matlab_ml_rl_v50_identity,
    validate_matlab_ml_rl_v51_identity,
    validate_matlab_ml_rl_v52_identity,
    validate_matlab_ml_rl_v53_identity,
    validate_matlab_ml_rl_v54_identity,
    validate_matlab_ml_rl_v55_identity,
)
from radia_mcp.matlab_agentic_ml.artifact_gate import validate_matlab_ml_rl_v46_identity
from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from radia_mcp.radia_ngsolve.server import (
    regularized_trace_inverse_path_gate as server_regularized_trace_inverse_path_gate,
)

from _matlab_generalization_payloads import (
    _CQ_KEY_V39,
    _CQ_KEY_V41,
    _DUCT_KEY,
    _FEMBEM_KEY,
    _HMATRIX_KEY,
    _LOW_KEY,
    _NONLINEAR_KEY,
    _OPTIMIZATION_KEY,
    _identity_v45,
    _summary_v23,
    _summary_v24,
    _summary_v25,
    _summary_v26,
    _summary_v27,
    _summary_v28,
    _summary_v29,
    _summary_v30,
    _summary_v31,
    _summary_v32,
    _summary_v33,
    _summary_v34,
    _summary_v35,
    _summary_v36,
    _summary_v37,
    _summary_v38,
    _summary_v39,
    _summary_v40,
    _summary_v41,
    _summary_v42,
    _summary_v44,
    _summary_v46,
    _summary_v47,
    _summary_v48,
    _summary_v49,
    _summary_v50,
    _summary_v51,
    _summary_v52,
    _summary_v53,
    _summary_v54,
    _summary_v55,
)


def test_v23_public_parallel_pool_worker_path_device_rng_code_generation_mismatch():
    summary = _summary_v23()
    summary["parallel_pool_worker_path_device_rng_code_generation_identity"].update(
        {
            "worker_path_pool_generation": "parallel-pool-50",
            "device_pool_generation": "parallel-pool-49",
            "rng_pool_generation": "parallel-pool-48",
            "code_pool_generation": "parallel-pool-47",
            "result_pool_generation": "parallel-pool-46",
            "result_worker_ids": [1, 2, 4],
            "result_worker_code_paths": ["toolbox/a", "toolbox/old", "toolbox/a"],
            "result_device_assignments": ["cpu:0", "gpu:0", "cpu:3"],
            "result_random_stream_seeds": [101, 999, 404],
            "result_worker_code_sha256": "c" * 64,
            "assembled_parallel_result_sha256": "d" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "parallel_results_use_current_worker_paths_devices_rng_and_code"
    ]


def test_v23_public_autodiff_tape_variable_order_mesh_objective_generation_mismatch():
    summary = _summary_v23()
    summary["autodiff_tape_variable_order_mesh_objective_generation_identity"].update(
        {
            "variable_order_tape_generation": "autodiff-tape-50",
            "mesh_tape_generation": "autodiff-tape-49",
            "objective_scaling_tape_generation": "autodiff-tape-48",
            "primal_solve_tape_generation": "autodiff-tape-47",
            "gradient_result_tape_generation": "autodiff-tape-46",
            "gradient_variable_ids": ["thickness", "radius", "impedance"],
            "gradient_mesh_sha256": "e" * 64,
            "gradient_objective_id": "mass",
            "gradient_objective_scale": 1.0,
            "gradient_primal_state_sha256": "f" * 64,
            "reported_gradient_table_sha256": "0" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "autodiff_gradients_use_current_tape_variables_mesh_objective_and_primal"
    ]


def test_v24_public_fembem_trace_normal_node_order_unit_mismatch() -> None:
    summary = _summary_v24()
    summary["fembem_trace_normal_interface_node_order_unit_generation_identity"].update(
        {
            "trace_coupling_generation": "fembem-100",
            "normal_coupling_generation": "fembem-99",
            "node_order_coupling_generation": "fembem-98",
            "unit_coupling_generation": "fembem-97",
            "operator_coupling_generation": "fembem-96",
            "result_trace_orientation": "boundary_to_volume",
            "result_outward_normal_convention": "inward_to_volume",
            "result_interface_node_ids": [1, 3, 2, 4],
            "result_boundary_triangles": [[1, 3, 2], [1, 2, 4]],
            "result_outward_normals": [[0.0, 0.0, -1.0], [0.0, -1.0, 0.0]],
            "result_physical_units": {"pressure": "kPa", "normal_velocity": "mm/s"},
            "result_interface_mesh_sha256": "a" * 64,
            "result_coupled_operator_sha256": "b" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fembem_trace_uses_current_normals_nodes_units_mesh_and_operator"
    ]


def test_v24_public_cq_contour_weight_startup_causality_window_mismatch() -> None:
    summary = _summary_v24()
    summary[
        "cq_contour_weight_startup_causality_window_result_generation_identity"
    ].update(
        {
            "contour_cq_generation": "cq-time-100",
            "weight_cq_generation": "cq-time-99",
            "startup_cq_generation": "cq-time-98",
            "causality_window_cq_generation": "cq-time-97",
            "time_grid_cq_generation": "cq-time-96",
            "result_contour_points_ri": [[1.2, 0.0], [0.0, 0.7], [-0.7, 0.0], [0.0, -0.7]],
            "result_cq_weights_ri": [[1.0, 0.0], [-1.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            "result_startup_weights_ri": [[0.0, 0.0], [0.0, 0.0]],
            "result_time_samples_s": [0.001, 0.002, 0.003, 0.004],
            "result_causality_window_s": [-0.001, 0.003],
            "result_prehistory_norm": 0.2,
            "reported_cq_result_sha256": "c" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_time_history_uses_current_contour_weights_startup_and_causality_window"
    ]


def test_v25_public_hmatrix_identity_mismatch() -> None:
    summary = _summary_v25()
    identity = summary[
        "hmatrix_block_tree_admissibility_permutation_tolerance_kernel_mesh_generation_identity"
    ]
    identity.update(
        {
            "block_tree_hmatrix_generation": "hmatrix-200",
            "admissibility_hmatrix_generation": "hmatrix-199",
            "permutation_hmatrix_generation": "hmatrix-198",
            "tolerance_hmatrix_generation": "hmatrix-197",
            "kernel_hmatrix_generation": "hmatrix-196",
            "mesh_hmatrix_generation": "hmatrix-195",
            "result_matrix_shape": [5, 4],
            "result_block_tree_sha256": "f" * 64,
            "result_admissibility_rule": "always_admissible",
            "result_admissibility_eta": 3.0,
            "result_row_permutation": [2, 0, 3, 3],
            "result_relative_tolerance": 1.0e-2,
            "result_kernel_id": "laplace_single_layer_p0",
            "result_boundary_mesh_sha256": "0" * 64,
            "reported_hmatrix_result_sha256": "1" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hmatrix_uses_current_block_tree_admissibility_permutations_tolerance_kernel_and_mesh"
    ]


def test_v25_public_ad_gradient_identity_mismatch() -> None:
    summary = _summary_v25()
    identity = summary[
        "ad_parameter_tape_material_operator_mesh_objective_primal_gradient_generation_identity"
    ]
    identity.update(
        {
            "parameter_tape_ad_generation": "ad-gradient-200",
            "material_law_ad_generation": "ad-gradient-199",
            "operator_ad_generation": "ad-gradient-198",
            "mesh_ad_generation": "ad-gradient-197",
            "objective_ad_generation": "ad-gradient-196",
            "primal_ad_generation": "ad-gradient-195",
            "result_parameter_names": ["bulk_modulus", "density"],
            "result_parameter_values": [142000.0, 1.2],
            "result_material_law_id": "nonlinear_fluid_previous",
            "result_objective_id": "source_power_previous",
            "result_parameter_tape_sha256": "2" * 64,
            "result_assembled_operator_sha256": "3" * 64,
            "result_mesh_sha256": "4" * 64,
            "result_primal_solution_sha256": "5" * 64,
            "finite_difference_gradient": [-0.25, 0.04],
            "maximum_gradient_relative_error": 2.0,
            "reported_gradient_result_sha256": "6" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "ad_gradient_uses_current_tape_material_operator_mesh_objective_and_primal"
    ]


def test_v26_public_rejects_cq_transfer_identity_mismatch() -> None:
    summary = _summary_v26()
    identity = summary[
        "cq_contour_radius_timestep_laplace_branch_transfer_operator_inverse_transform_generation_identity"
    ]
    identity.update(
        {
            "contour_cq_generation": "cq-transfer-300",
            "result_contour_radius": 0.8,
            "result_time_step_s": 2.0e-4,
            "result_laplace_branch": "negative_sqrt_incoming",
            "result_laplace_points_ri": [[-100.0, 0.0]],
            "result_transfer_operator_id": "laplace_single_layer_p0",
            "result_transfer_operator_sha256": "1" * 63 + "2",
            "result_inverse_transform": "direct_dft_unsigned",
            "reported_time_history_sha256": "2" * 63 + "3",
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_history_uses_current_contour_timestep_branch_transfer_and_inverse_transform"
    ]


def test_v26_public_rejects_fembem_coupling_identity_mismatch() -> None:
    summary = _summary_v26()
    identity = summary[
        "fembem_trace_map_normal_material_wavenumber_coupling_matrix_mesh_generation_identity"
    ]
    identity.update(
        {
            "trace_map_coupling_generation": "fembem-coupling-300",
            "result_trace_map_sha256": "3" * 63 + "4",
            "result_normal_orientation": "volume_inward",
            "result_fluid_density_kg_m3": 1000.0,
            "result_sound_speed_m_s": 1480.0,
            "result_wavenumber_ri_m_inv": [4.2, -0.02],
            "result_coupling_matrix_sha256": "5" * 63 + "6",
            "result_volume_mesh_sha256": "6" * 63 + "7",
            "result_boundary_mesh_sha256": "7" * 63 + "8",
            "reported_coupled_result_sha256": "8" * 63 + "9",
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fembem_coupling_uses_current_trace_normals_material_wavenumber_matrices_and_mesh"
    ]


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


def test_v28_public_rejects_hmatrix_aca_identity_mismatch() -> None:
    summary = _summary_v28()
    identity = summary[
        "hmatrix_aca_cluster_permutation_admissibility_rank_tolerance_kernel_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "permutation_hmatrix_generation": "hmatrix-aca-320",
            "mesh_hmatrix_generation": "hmatrix-aca-319",
            "result_cluster_permutation": [1, 2, 3, 4],
            "result_admissibility_rule": "strong",
            "result_admissibility_eta": 0.5,
            "result_aca_rank": 3,
            "result_relative_tolerance": 1.0e-2,
            "result_kernel": "laplace-p0",
            "loaded_cluster_tree_sha256": "d" * 64,
            "result_mesh_sha256": "e" * 64,
            "accepted_result_sha256": "f" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hmatrix_aca_uses_current_clusters_permutation_admissibility_rank_tolerance_kernel_mesh_and_result"
    ]


def test_v28_public_rejects_calderon_cq_identity_mismatch() -> None:
    summary = _summary_v28()
    identity = summary[
        "calderon_cq_operator_v_k_trace_normal_frequency_grid_inverse_transform_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "v_calderon_generation": "calderon-cq-320",
            "normal_calderon_generation": "calderon-cq-319",
            "result_v_operator_sha256": "0" * 64,
            "result_k_operator_sha256": "1" * 64,
            "result_trace_basis": "p0-cell",
            "result_trace_shape": [47, 120],
            "result_boundary_normal": "volume-inward",
            "result_laplace_frequency_ri": [[0.0, 10.0]],
            "result_inverse_transform": "direct-real-ifft",
            "result_boundary_mesh_sha256": "2" * 64,
            "accepted_result_sha256": "3" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "calderon_cq_uses_current_v_k_trace_normals_frequency_inverse_mesh_and_result"
    ]


def test_v29_public_bem_near_singular_quadrature_distance_element_size_adaptive_order_reference_mismatch() -> None:
    summary = _summary_v29()
    identity = summary[
        "bem_near_singular_quadrature_distance_element_size_adaptive_order_reference_result_generation_identity"
    ]
    identity.update(
        {
            "target_quadrature_generation": "near-singular-330",
            "mesh_quadrature_generation": "near-singular-329",
            "result_target_distance_m": 1.0e-3,
            "result_element_size_m": 2.0e-3,
            "result_distance_size_ratio": 0.5,
            "result_adaptive_order": 4,
            "result_quadrature_rule": "triangle-centroid",
            "result_coordinate_map": "global-cartesian",
            "result_kernel": "laplace-p0",
            "computed_integral_ri": [0.2, 0.0],
            "relative_error": 0.5,
            "result_element_mesh_sha256": "a" * 64,
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "bem_near_singular_quadrature_uses_current_distance_size_order_map_kernel_reference_mesh_and_result"
    ]


def test_v29_public_fembem_energy_flux_reciprocity_interface_trace_orientation_frequency_mismatch() -> None:
    summary = _summary_v29()
    identity = summary[
        "fembem_energy_flux_reciprocity_interface_trace_orientation_frequency_incident_result_generation_identity"
    ]
    identity.update(
        {
            "trace_coupling_generation": "fembem-energy-330",
            "frequency_coupling_generation": "fembem-energy-329",
            "result_interface_trace_basis": "p0-cell",
            "result_interface_trace_shape": [47, 120],
            "result_normal_orientation": "volume-inward",
            "result_frequency_hz": 1000.0,
            "result_incident_field_sha256": "c" * 64,
            "result_reciprocity_pair_ids": ["source-b/receiver-a", "source-a/receiver-b"],
            "result_reciprocity_values_ri": [[0.4, 0.0], [0.7, 0.0]],
            "reciprocity_relative_error": 0.3,
            "bem_radiated_power_w": 0.5,
            "energy_flux_relative_error": 0.6,
            "accepted_coupled_result_sha256": "d" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fembem_energy_flux_and_reciprocity_use_current_trace_normal_frequency_incident_field_and_result"
    ]


def test_v30_public_hmatrix_recompression_svd_tolerance_norm_rank_permutation_operator_error_mismatch() -> None:
    summary = _summary_v30()
    identity = summary[
        "hmatrix_recompression_svd_tolerance_norm_rank_permutation_operator_mesh_result_identity"
    ]
    identity.update(
        {
            "svd_hmatrix_generation": "hmatrix-recompress-340",
            "rank_hmatrix_generation": "hmatrix-recompress-339",
            "result_svd_basis": "mass-weighted",
            "result_tolerance": 1.0e-2,
            "result_tolerance_norm": "frobenius-absolute",
            "result_block_ranks_before": [8, 10, 12],
            "result_block_ranks_after": [9, 11, 13],
            "result_row_permutation": [0, 1, 2],
            "result_column_permutation": [0, 1, 2],
            "result_operator_relative_error": 0.1,
            "result_mesh_sha256": "b" * 64,
            "accepted_result_sha256": "c" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hmatrix_recompression_uses_current_svd_tolerance_norm_ranks_permutations_operator_mesh_and_result"
    ]


def test_v30_public_cq_restart_block_history_startup_weights_time_index_sample_count_digest_mismatch() -> None:
    summary = _summary_v30()
    identity = summary[
        "cq_restart_block_history_startup_weight_time_index_sample_contour_operator_result_identity"
    ]
    identity.update(
        {
            "block_cq_generation": "cq-restart-340",
            "history_cq_generation": "cq-restart-339",
            "result_block_size": 8,
            "result_completed_block_ids": [0, 2],
            "result_history_sample_count": 40,
            "result_startup_weights_ri": [[0.0, 0.0]],
            "result_restart_time_index": 47,
            "result_total_sample_count": 64,
            "result_contour_owner_sha256": "d" * 64,
            "result_operator_owner_sha256": "e" * 64,
            "loaded_history_sha256": "f" * 64,
            "accepted_result_sha256": "0" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_block_restart_uses_current_blocks_history_startup_weights_time_samples_owners_and_result"
    ]


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


def test_v32_public_adaptive_cq_timestep_contour_rebuild_interpolation_history_error_restart_mismatch() -> None:
    payload = _summary_v32()
    identity = payload[
        "adaptive_cq_timestep_contour_history_interpolation_error_restart_operator_mesh_result_identity"
    ]
    identity.update(
        {
            "timestep_cq_generation": "adaptive-cq-360",
            "operator_cq_generation": "adaptive-cq-359",
            "result_cq_generation": "adaptive-cq-358",
            "result_timestep_schedule_s": [1.0e-5, 5.0e-6, 1.0e-5],
            "result_contour_rebuild_indices": [0],
            "result_history_interpolation": "linear_noncausal",
            "result_local_error_estimates": [0.2, 0.1, 0.05],
            "result_local_error_tolerance": 1.0e-4,
            "result_restart_index": 1,
            "result_operator_owner_sha256": "b" * 64,
            "loaded_history_sha256": "c" * 64,
            "result_mesh_sha256": "d" * 64,
            "accepted_result_sha256": "e" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "adaptive_cq_uses_current_timesteps_contour_rebuild_history_interpolation_error_restart_operator_mesh_and_result"
    ]


def test_v32_public_fembem_modal_transient_mass_damping_initial_condition_energy_balance_mismatch() -> None:
    payload = _summary_v32()
    identity = payload[
        "fembem_modal_transient_mass_damping_initial_projection_truncation_energy_mesh_history_result_identity"
    ]
    identity.update(
        {
            "mass_modal_generation": "modal-fembem-360",
            "history_modal_generation": "modal-fembem-359",
            "result_modal_generation": "modal-fembem-358",
            "result_mass_normalization": "euclidean",
            "result_damping_model": "none",
            "result_rayleigh_coefficients": [0.0, 0.0],
            "result_initial_displacement_projection": [0.0, 1.0],
            "result_initial_velocity_projection": [1.0, 0.0],
            "result_modal_count": 2,
            "result_truncation_frequency_hz": 600.0,
            "result_initial_energy_j": 0.4,
            "result_radiated_energy_j": 0.6,
            "result_dissipated_energy_j": -0.1,
            "result_final_energy_j": 0.3,
            "result_mesh_sha256": "f" * 64,
            "accepted_time_history_owner": "fembem/old-history",
            "accepted_result_sha256": "0" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "modal_fembem_transient_uses_current_mass_damping_initial_projection_truncation_energy_mesh_history_and_result"
    ]


def test_v33_public_calderon_projector_v_k_kt_w_mass_duality_normal_quadrature_mismatch() -> None:
    payload = _summary_v33()
    identity = payload[
        "calderon_projector_p1_v_k_kt_w_mass_duality_normal_quadrature_mesh_owner_result_identity"
    ]
    identity.update(
        {
            "operator_calderon_generation": "calderon-p1-380",
            "mesh_calderon_generation": "calderon-p1-379",
            "result_calderon_generation": "calderon-p1-378",
            "result_trial_space": "P0",
            "result_test_space": "P0",
            "result_projector_convention": "exterior_stale_normal",
            "result_block_order": ["neumann", "dirichlet"],
            "result_v_sign": 1,
            "result_k_sign": -1,
            "result_kt_sign": 1,
            "result_w_sign": 1,
            "result_mass_duality_residual": 0.2,
            "result_mass_duality_tolerance": 1.0e-4,
            "result_normal_orientation": "inward",
            "result_singular_quadrature": "centroid_regular",
            "result_projector_residual": 0.4,
            "result_projector_tolerance": 1.0e-4,
            "result_operator_sha256": "c" * 64,
            "result_mass_sha256": "d" * 64,
            "result_boundary_mesh_sha256": "e" * 64,
            "accepted_result_owner": "fembem/calderon/old",
            "accepted_result_sha256": "f" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "calderon_projector_uses_current_p1_spaces_v_k_kt_w_mass_duality_normals_quadrature_mesh_owner_and_result"
    ]


def test_v33_public_cq_symbol_contour_conjugate_symmetry_causal_ifft_parseval_passivity_mismatch() -> None:
    payload = _summary_v33()
    identity = payload[
        "cq_symbol_contour_transfer_conjugate_causal_ifft_parseval_passivity_timestep_operator_result_identity"
    ]
    identity.update(
        {
            "symbol_cq_generation": "cq-physical-380",
            "transfer_cq_generation": "cq-physical-379",
            "result_cq_generation": "cq-physical-378",
            "result_multistep_symbol": "BDF1",
            "result_symbol_coefficients": [1.0, -1.0],
            "result_contour_radius": 1.1,
            "result_transfer_samples_ri": [[1.0, 0.1], [0.8, 0.3]],
            "result_time_response": [0.4, -0.2],
            "result_negative_time_energy": 0.1,
            "result_time_domain_work": 0.1,
            "result_frequency_domain_work": 0.7,
            "result_parseval_tolerance": 1.0e-4,
            "result_minimum_real_transfer": -0.3,
            "result_passivity_sign": "negative_real_transfer",
            "result_timestep_s": -1.0e-5,
            "result_operator_family": "stale_operator",
            "result_operator_sha256": "0" * 64,
            "accepted_result_owner": "fembem/cq/old",
            "accepted_result_sha256": "1" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_uses_current_symbol_contour_conjugate_transfer_causal_ifft_parseval_passivity_timestep_operator_and_result"
    ]


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


def test_v37_public_hmatrix_dense_reference_error_memory_rank_complexity_mesh_owner_mismatch():
    payload = _summary_v37()
    identity = payload["hmatrix_dense_reference_error_tolerance_rank_memory_complexity_mesh_operator_benchmark_result_identity"]
    identity.update({
        "dense_generation": "hmatrix-benchmark-512", "memory_generation": "hmatrix-benchmark-511",
        "result_generation": "hmatrix-benchmark-510", "result_dense_reference_relative_error": [9.0],
        "result_relative_tolerance": -1.0, "result_maximum_block_rank": [999],
        "result_dense_memory_bytes": [1], "result_hmatrix_memory_bytes": [999999999],
        "result_memory_complexity_exponent": 3.0, "result_rank_complexity_exponent": 2.0,
        "result_boundary_mesh_sha256": "8" * 64, "accepted_operator_owner": "hmatrix/old",
        "accepted_benchmark_owner": "benchmark/old", "accepted_hmatrix_result_sha256": "9" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["hmatrix_benchmarks_use_current_dense_error_tolerance_rank_memory_complexity_mesh_owners_and_result"]


def test_v37_public_multifrequency_fembem_adjoint_gradient_quadrature_trace_fd_owner_mismatch():
    payload = _summary_v37()
    identity = payload["multifrequency_fembem_adjoint_weight_objective_quadrature_trace_gradient_fd_mesh_owner_result_identity"]
    identity.update({
        "weight_generation": "multifrequency-adjoint-512", "trace_generation": "multifrequency-adjoint-511",
        "result_generation": "multifrequency-adjoint-510", "result_frequency_weights": [2.0, -1.0],
        "result_objective_complex": [[99.0, 99.0]], "result_weighted_objective_complex": [99.0, 99.0],
        "result_quadrature_order": [0], "result_trace_node_map": [3, 2, 1],
        "result_frequency_gradient": [99.0], "result_accumulated_gradient": -99.0,
        "result_finite_difference_gradient": 99.0, "result_fembem_mesh_sha256": "a" * 64,
        "accepted_adjoint_owner": "fembem/old", "accepted_adjoint_result_sha256": "b" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["multifrequency_fembem_adjoints_use_current_weights_objective_quadrature_trace_gradient_fd_mesh_owner_and_result"]


def test_v37_public_rejects_self_consistent_dense_storage_as_hmatrix():
    payload = _summary_v37()
    identity = payload["hmatrix_dense_reference_error_tolerance_rank_memory_complexity_mesh_operator_benchmark_result_identity"]
    identity["hmatrix_memory_bytes"] = list(identity["dense_memory_bytes"])
    identity["result_hmatrix_memory_bytes"] = list(identity["dense_memory_bytes"])
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_wrong_gradient_accumulation():
    payload = _summary_v37()
    identity = payload["multifrequency_fembem_adjoint_weight_objective_quadrature_trace_gradient_fd_mesh_owner_result_identity"]
    identity["accumulated_gradient"] = 99.0
    identity["result_accumulated_gradient"] = 99.0
    identity["finite_difference_gradient"] = 99.0
    identity["result_finite_difference_gradient"] = 99.0
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v38_public_simp_topology_filter_projection_volume_compliance_adjoint_fd_kkt_mismatch():
    payload = _summary_v38()
    identity = payload["simp_topology_density_filter_projection_volume_compliance_adjoint_fd_kkt_mesh_owner_result_identity"]
    identity.update({
        "filter_generation": "simp-topology-613",
        "kkt_generation": "simp-topology-612",
        "result_generation": "simp-topology-611",
        "result_density_filter_matrix": [[2.0, -1.0]],
        "result_filtered_density": [2.0],
        "result_projection_beta_continuation": [8.0, 1.0],
        "result_projected_density": [-1.0, 2.0],
        "result_volume_fraction": 2.0,
        "result_compliance": -12.0,
        "result_adjoint_compliance_gradient": [9.0],
        "result_finite_difference_compliance_gradient": [-9.0],
        "result_kkt_stationarity_residual": 9.0,
        "accepted_topology_owner": "optimization/old",
        "accepted_topology_result_sha256": "a" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["simp_topology_uses_current_density_filter_projection_volume_compliance_adjoint_fd_kkt_mesh_owner_and_result"]


def test_v38_public_fembem_model_reduction_projection_stability_passivity_moment_error_owner_mismatch():
    payload = _summary_v38()
    identity = payload["fembem_model_reduction_projection_order_stability_passivity_moment_frequency_error_full_mesh_owner_result_identity"]
    identity.update({
        "projection_generation": "fembem-reduction-613",
        "passivity_generation": "fembem-reduction-612",
        "result_generation": "fembem-reduction-611",
        "result_reduced_order": 5,
        "result_trial_projection_basis": [[1.0]],
        "result_test_projection_basis": [[-1.0]],
        "result_biorthogonality_gram": [[-1.0]],
        "result_reduced_model_poles": [[1.0, 0.0]],
        "result_minimum_passivity_eigenvalue": -1.0,
        "result_matched_moments_reduced": [[9.0, 9.0]],
        "result_frequency_response_relative_error": [1.0],
        "accepted_full_model_owner": "fembem/old",
        "accepted_reduction_owner": "fembem/old",
        "accepted_reduction_result_sha256": "b" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["fembem_model_reduction_uses_current_projection_order_stability_passivity_moments_frequency_error_full_model_mesh_owners_and_result"]


def test_v38_public_rejects_self_consistent_nonstochastic_density_filter():
    payload = _summary_v38()
    identity = payload["simp_topology_density_filter_projection_volume_compliance_adjoint_fd_kkt_mesh_owner_result_identity"]
    identity["density_filter_matrix"] = [[1.5, -0.5, 0.0, 0.0]] * 4
    identity["result_density_filter_matrix"] = deepcopy(identity["density_filter_matrix"])
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_unstable_reduced_pole():
    payload = _summary_v38()
    identity = payload["fembem_model_reduction_projection_order_stability_passivity_moment_frequency_error_full_mesh_owner_result_identity"]
    identity["reduced_model_poles"] = [[1.0, 0.0], [-2.1, 0.0]]
    identity["result_reduced_model_poles"] = deepcopy(identity["reduced_model_poles"])
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


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
    identity = payload[_CQ_KEY_V39]
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
    identity = payload[_CQ_KEY_V39]
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


def test_v41_public_cq_acoustic_causality_passivity_timestep_ztransform_energy_history_mismatch() -> None:
    payload = _summary_v41()
    payload[_CQ_KEY_V41].update({
        "timestep_generation": "cq-acoustic-730",
        "result_time_step_s": -1.0,
        "result_z_transform_radius": 1.1,
        "result_pressure_history": [1.0, 9.0],
        "result_minimum_passivity_real_part": -1.0,
        "result_energy_balance_residual_j": 1.0,
        "accepted_cq_owner": "acoustic/old",
        "accepted_cq_result_sha256": "a" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_acoustic_history_uses_current_causality_passivity_timestep_ztransform_energy_mesh_owner_and_result"
    ]


def test_v41_public_hmatrix_admissibility_cluster_rank_tolerance_matvec_error_memory_mismatch() -> None:
    payload = _summary_v41()
    payload[_HMATRIX_KEY].update({
        "cluster_generation": "hmatrix-730",
        "result_cluster_permutation": [1, 1, 9],
        "result_admissibility_eta": -1.0,
        "result_numerical_ranks": [[-1]],
        "result_measured_matvec_relative_error": 1.0,
        "result_compressed_memory_bytes": 262144,
        "accepted_hmatrix_owner": "acoustic/old",
        "accepted_hmatrix_result_sha256": "b" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hmatrix_compression_uses_current_clusters_admissibility_ranks_tolerance_matvec_memory_mesh_owner_and_result"
    ]


def test_v41_public_rejects_self_consistent_noncausal_pressure_history() -> None:
    payload = _summary_v41()
    payload[_CQ_KEY_V41]["pressure_history"] = [1.0, 0.2, 0.3, 0.15, 0.04]
    payload[_CQ_KEY_V41]["result_pressure_history"] = [1.0, 0.2, 0.3, 0.15, 0.04]
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_hmatrix_error_or_memory_regression() -> None:
    payload = _summary_v41()
    payload[_HMATRIX_KEY]["measured_matvec_relative_error"] = 2.0e-5
    payload[_HMATRIX_KEY]["result_measured_matvec_relative_error"] = 2.0e-5
    payload[_HMATRIX_KEY]["compressed_memory_bytes"] = 131072
    payload[_HMATRIX_KEY]["result_compressed_memory_bytes"] = 131072
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v42_public_positive_lowfrequency_bem_and_duct_scattering_closure() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v42())["status"] == "ok"


def test_v42_public_lowfrequency_bem_stabilized_kernel_condition_staticlimit_charge_energy_mismatch() -> None:
    payload = _summary_v42()
    payload[_LOW_KEY].update({
        "kernel_generation": "low-frequency-bem-841",
        "energy_generation": "low-frequency-bem-840",
        "result_generation": "low-frequency-bem-839",
        "result_stabilized_kernel": "unstabilized",
        "result_static_limit_residual": [1.0],
        "result_condition_estimate": [-1.0],
        "result_boundary_charge_c": [1.0],
        "result_potential_energy_j": -1.0,
        "accepted_boundary_mesh_owner": "stale:mesh",
        "accepted_bem_result_sha256": "9" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "lowfrequency_bem_uses_current_stabilized_kernel_static_limit_condition_charge_energy_mesh_owner_and_result"
    ]


def test_v42_public_duct_scattering_transmission_reflection_loss_power_modal_balance_mismatch() -> None:
    payload = _summary_v42()
    payload[_DUCT_KEY].update({
        "mode_generation": "duct-scattering-841",
        "power_generation": "duct-scattering-840",
        "result_generation": "duct-scattering-839",
        "result_incident_mode_pressure": [0.0, 0.0],
        "result_reflected_mode_pressure": [2.0, 0.0],
        "result_transmitted_mode_pressure": [-1.0, 0.0],
        "result_transmission_loss_db": -10.0,
        "result_modal_power_balance_residual_w": 1.0,
        "accepted_mesh_owner": "stale:mesh",
        "accepted_duct_result_sha256": "a" * 64,
    })
    result = regularized_trace_inverse_path_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "duct_scattering_uses_current_modes_pressures_transmission_reflection_loss_power_mesh_owner_and_result"
    ]


def test_v42_public_rejects_self_consistent_non_neutral_lowfrequency_charge() -> None:
    payload = _summary_v42()
    payload[_LOW_KEY]["boundary_charge_c"] = [1.0e-9, 1.0e-9]
    payload[_LOW_KEY]["result_boundary_charge_c"] = [1.0e-9, 1.0e-9]
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_duct_power_imbalance() -> None:
    payload = _summary_v42()
    payload[_DUCT_KEY]["reflected_modal_power_w"] = 0.4
    payload[_DUCT_KEY]["result_reflected_modal_power_w"] = 0.4
    assert regularized_trace_inverse_path_gate(payload)["status"] == "needs_attention"


def test_v44_public_ml_rl_identity_positive_is_accepted() -> None:
    result = validate_matlab_ml_rl_v44_identity(_summary_v44())
    assert result is not None
    assert result["status"] == "ok"


def test_v44_public_supervised_rejects_split_and_validation_contamination() -> None:
    value = _summary_v44()
    identity = value["matlab_ml_rl_v44_identity"]["supervised"]
    identity.update({
        "evaluation_ids": ["train-02"],
        "training_evaluation_disjoint": False,
        "validation_metric": {"name": "training_rmse", "value": 0.0},
        "accepted_result_sha256": "8" * 64,
    })
    result = server_regularized_trace_inverse_path_gate(json.dumps(value))
    assert json.loads(result)["status"] == "needs_attention"


def test_v44_public_rl_rejects_environment_seed_episode_and_exploration_mismatch() -> None:
    value = _summary_v44()
    identity = value["matlab_ml_rl_v44_identity"]["reinforcement_learning"]
    identity.update({
        "evaluation_environment_id": "cae:rl-design-control-old",
        "evaluation_seed": 845,
        "evaluation_episodes": 0,
        "exploration_during_evaluation": True,
        "accepted_evaluation_result_sha256": "9" * 64,
    })
    result = validate_matlab_ml_rl_v44_identity(value)
    assert result is not None
    assert result["status"] == "needs_attention"


def test_v45_matlab_ml_rl_identity_positive() -> None:
    result = validate_matlab_ml_rl_v45_identity(_identity_v45())
    assert result is not None
    assert result["status"] == "ok"


def test_v45_matlab_identity_rejects_holdout_and_agentic_consent_mutations() -> None:
    identity = _identity_v45()
    identity["matlab_ml_rl_v45_identity"]["supervised"]["holdout_id"] = "train:1"
    identity["matlab_ml_rl_v45_identity"]["agentic_toolkit"]["result_consent_recorded"] = False
    result = validate_matlab_ml_rl_v45_identity(identity)
    assert result is not None
    assert result["status"] == "needs_attention"


def test_v46_matlab_ml_rl_identity_accepts_closed_artifact():
    result = validate_matlab_ml_rl_v46_identity(_summary_v46())
    assert result and result["status"] == "ok" and all(result["checks"].values())


def test_v46_matlab_ml_rl_identity_rejects_nan_seed_episode_and_argument_mutations():
    summary = _summary_v46()
    identity = summary["matlab_ml_rl_v46_identity"]
    identity["supervised"]["result_nonfinite_policy"] = "keep_nan"
    identity["reinforcement_learning"]["result_episode_timeout_steps"] = 1
    identity["agentic_toolkit"]["result_argument_shape_valid"] = False
    identity["mlrl_checkpoint"]["result_checkpoint_order"] = [0, 2, 1]
    result = validate_matlab_ml_rl_v46_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v47_positive_ml_rl_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v47_identity(_summary_v47())
    assert result and result["status"] == "ok"


def test_v47_ml_datastore_partition_mutation_is_rejected() -> None:
    summary = _summary_v47()
    row = summary["matlab_ml_rl_v47_identity"]["ml_datastore"]
    row["result_datastore_order"] = ["sample-2", "sample-1", "sample-3"]
    row["result_labels"] = ["dog", "cat", "cat"]
    row["result_preprocess_sha256"] = "a" * 64
    row["result_partition_owner"] = "partition:other"
    result = validate_matlab_ml_rl_v47_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v47_rl_spec_environment_policy_mutation_is_rejected() -> None:
    summary = _summary_v47()
    row = summary["matlab_ml_rl_v47_identity"]["reinforcement_learning"]
    row["result_observation_spec_sha256"] = "b" * 64
    row["result_action_spec_sha256"] = "c" * 64
    row["result_environment_id"] = "environment:other"
    row["result_reset_policy"] = "random_unseeded"
    row["result_policy_owner"] = "policy:other"
    result = validate_matlab_ml_rl_v47_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v48_positive_ml_rl_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v48_identity(_summary_v48())
    assert result and result["status"] == "ok"


def test_v48_autodiff_mutation_is_rejected() -> None:
    summary = _summary_v48()
    row = summary["matlab_ml_rl_v48_identity"]["autodiff"]
    row.update({"result_parameter_order": ["w2", "w1"], "result_gradient_tape_sha256": "a" * 64, "result_objective_id": "objective:training-loss", "result_gradient": [-0.25, 1.5], "result_fd_gradient": [1.2, -0.1], "result_checkpoint_owner": "checkpoint:old"})
    result = validate_matlab_ml_rl_v48_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v48_sequence_mutation_is_rejected() -> None:
    summary = _summary_v48()
    row = summary["matlab_ml_rl_v48_identity"]["sequence_model"]
    row.update({"result_padding_policy": "left_zero", "result_padding_mask": [[1, 1, 1], [0, 1, 1]], "result_shuffle_order": [0, 1], "result_checkpoint_owner": "checkpoint:old"})
    result = validate_matlab_ml_rl_v48_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v48_self_consistent_invalid_conventions_are_rejected() -> None:
    summary = _summary_v48()
    identity = summary["matlab_ml_rl_v48_identity"]
    identity["sequence_model"]["padding_policy"] = identity["sequence_model"]["result_padding_policy"] = "left_zero"
    identity["agentic_workspace"]["approval_scope"] = identity["agentic_workspace"]["result_approval_scope"] = "workspace_write:any"
    result = validate_matlab_ml_rl_v48_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v49_positive_ml_rl_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v49_identity(_summary_v49())
    assert result and result["status"] == "ok"


def test_v49_rl_replay_mutation_is_rejected() -> None:
    summary = _summary_v49(); row = summary["matlab_ml_rl_v49_identity"]["rl_replay"]
    row.update({"result_replay_row_keys": ["transition:1", "transition:0", "transition:2"], "result_terminals": [False, True, False], "result_episode_seeds": [4102, 4101, 4103], "result_policy_owner": "policy:other"})
    result = validate_matlab_ml_rl_v49_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v49_ml_normalization_mutation_is_rejected() -> None:
    summary = _summary_v49(); row = summary["matlab_ml_rl_v49_identity"]["ml_normalization"]
    row.update({"result_normalization_fit_scope": "all_samples", "result_class_encoding": {"cold": 1, "hot": 0}, "result_training_fold_ids": [1, 0, 1, 0], "result_model_owner": "model:other"})
    result = validate_matlab_ml_rl_v49_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v49_self_consistent_data_leakage_and_shared_streams_are_rejected() -> None:
    summary = deepcopy(_summary_v49()); identity = summary["matlab_ml_rl_v49_identity"]
    ml = identity["ml_normalization"]; ml["normalization_fit_scope"] = ml["result_normalization_fit_scope"] = "all_samples"
    parallel = identity["parallel_resume"]; parallel["random_streams"] = parallel["result_random_streams"] = ["global", "global", "global"]
    result = validate_matlab_ml_rl_v49_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v50_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v50_identity(_summary_v50())
    assert result and result["status"] == "ok"


def test_v50_bayesopt_mutation_is_rejected() -> None:
    summary = deepcopy(_summary_v50())
    row = summary["matlab_ml_rl_v50_identity"]["bayesopt"]
    row.update({"result_objective_name": "training_loss", "result_constraints": [], "result_rng_seed": 999, "result_model_owner": "model:foreign"})
    result = validate_matlab_ml_rl_v50_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v50_mixed_precision_mutation_is_rejected() -> None:
    summary = deepcopy(_summary_v50())
    row = summary["matlab_ml_rl_v50_identity"]["mixed_precision"]
    row.update({"result_precision_policy": "fp32", "result_loss_scale": 1.0, "result_checkpoint_sha256": "9" * 64, "result_network_owner": "network:foreign"})
    result = validate_matlab_ml_rl_v50_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v50_self_consistent_unsafe_states_are_rejected() -> None:
    summary = deepcopy(_summary_v50())
    mixed = summary["matlab_ml_rl_v50_identity"]["mixed_precision"]
    mixed["precision_policy"] = mixed["result_precision_policy"] = "fp32"
    sandbox = summary["matlab_ml_rl_v50_identity"]["agentic_sandbox"]
    sandbox["requested_path"] = sandbox["result_requested_path"] = "workspace:project/../private/key.txt"
    sandbox["traversal_detected"] = sandbox["result_traversal_detected"] = True
    result = validate_matlab_ml_rl_v50_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v51_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v51_identity(_summary_v51())
    assert result and result["status"] == "ok"


def test_v51_frozen_public_counterfactuals_are_rejected() -> None:
    summary = deepcopy(_summary_v51())
    summary["matlab_ml_rl_v51_identity"]["rl_replay"].update({"result_transition_ids": ["transition:4"], "result_buffer_owner": "buffer:stale"})
    summary["matlab_ml_rl_v51_identity"]["classification"].update({"result_class_order": ["fault", "warning", "normal"], "result_model_owner": "model:stale"})
    result = validate_matlab_ml_rl_v51_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v51_self_consistent_wrong_semantics_are_rejected() -> None:
    summary = deepcopy(_summary_v51())
    rl = summary["matlab_ml_rl_v51_identity"]["rl_replay"]
    rl["discounts"] = rl["result_discounts"] = [0.99] * 4
    classification = summary["matlab_ml_rl_v51_identity"]["classification"]
    classification["class_order"] = classification["result_class_order"] = ["fault", "warning", "normal"]
    result = validate_matlab_ml_rl_v51_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v52_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v52_identity(_summary_v52())
    assert result and result["status"] == "ok"


def test_v52_frozen_public_counterfactuals_are_rejected() -> None:
    summary = deepcopy(_summary_v52())
    training = summary["matlab_ml_rl_v52_identity"]["ml_training"]
    training.update({"result_optimizer_name": "sgdm", "result_checkpoint_sha256": "0" * 64, "result_model_owner": "model:stale"})
    recurrent = summary["matlab_ml_rl_v52_identity"]["recurrent_rl"]
    recurrent.update({"result_hidden_state_sha256": list(reversed(recurrent["hidden_state_sha256"])), "result_episode_reset": [False, False, True, False], "result_sequence_mask": [True, False, True, True], "result_policy_owner": "policy:stale"})
    result = validate_matlab_ml_rl_v52_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v52_self_consistent_invalid_training_and_sequence_semantics_are_rejected() -> None:
    summary = deepcopy(_summary_v52())
    training = summary["matlab_ml_rl_v52_identity"]["ml_training"]
    training["learning_rate_schedule"] = training["result_learning_rate_schedule"] = [0.001, 0.002, 0.0001]
    recurrent = summary["matlab_ml_rl_v52_identity"]["recurrent_rl"]
    recurrent["episode_reset"] = recurrent["result_episode_reset"] = [False, False, True, False]
    recurrent["sequence_mask"] = recurrent["result_sequence_mask"] = [True, False, True, True]
    result = validate_matlab_ml_rl_v52_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v52_source_contracts_reject_unsupported_or_unowned_replay() -> None:
    summary = deepcopy(_summary_v52())
    codegen = summary["matlab_ml_rl_v52_identity"]["codegen"]
    codegen["numeric_type"] = codegen["result_numeric_type"] = "half"
    sdi = summary["matlab_ml_rl_v52_identity"]["simulink_data_inspector"]
    sdi["interpolation"] = sdi["result_interpolation"] = ["spline", "zoh", "zoh"]
    result = validate_matlab_ml_rl_v52_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v53_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v53_identity(_summary_v53())
    assert result and result["status"] == "ok"


def test_v53_frozen_public_counterfactuals_are_rejected() -> None:
    summary = deepcopy(_summary_v53())
    offline = summary["matlab_ml_rl_v53_identity"]["offline_rl"]
    offline.update({"result_behavior_policy_sha256": "8" * 64, "result_importance_weight": [10.0] * 4, "result_support_coverage_fraction": 0.25, "result_dataset_owner": "dataset:stale"})
    calibration = summary["matlab_ml_rl_v53_identity"]["probability_calibration"]
    calibration.update({"result_calibrated_probability": [[0.9, 0.9, 0.9]], "result_class_order": ["class:fault", "class:normal"], "result_temperature": 0.1, "result_model_owner": "model:stale"})
    result = validate_matlab_ml_rl_v53_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v53_self_consistent_invalid_offline_semantics_are_rejected() -> None:
    summary = deepcopy(_summary_v53())
    offline = summary["matlab_ml_rl_v53_identity"]["offline_rl"]
    offline["importance_weight"] = offline["result_importance_weight"] = [1.0] * 4
    offline["support_coverage_fraction"] = offline["result_support_coverage_fraction"] = 0.5
    result = validate_matlab_ml_rl_v53_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v53_self_consistent_invalid_calibration_semantics_are_rejected() -> None:
    summary = deepcopy(_summary_v53())
    calibration = summary["matlab_ml_rl_v53_identity"]["probability_calibration"]
    calibration["calibrated_probability"] = calibration["result_calibrated_probability"] = [[0.8, 0.1, 0.1], [0.8, 0.1, 0.1]]
    calibration["class_order"] = calibration["result_class_order"] = ["class:normal", "class:normal", "class:fault"]
    calibration["temperature"] = calibration["result_temperature"] = -1.0
    result = validate_matlab_ml_rl_v53_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v54_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v54_identity(_summary_v54())
    assert result and result["status"] == "ok"


def test_v54_frozen_public_counterfactuals_are_rejected() -> None:
    summary = deepcopy(_summary_v54())
    summary["matlab_ml_rl_v54_identity"]["prioritized_replay"]["result_rng_seed"] = 42
    summary["matlab_ml_rl_v54_identity"]["feature_preprocess"]["result_model_checkpoint_sha256"] = "9" * 64
    result = validate_matlab_ml_rl_v54_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v54_self_consistent_invalid_semantics_are_rejected() -> None:
    summary = deepcopy(_summary_v54())
    replay = summary["matlab_ml_rl_v54_identity"]["prioritized_replay"]
    replay["importance_weight"] = replay["result_importance_weight"] = [1.0] * 4
    preprocess = summary["matlab_ml_rl_v54_identity"]["feature_preprocess"]
    preprocess["data_split_indices"] = preprocess["result_data_split_indices"] = {"train": [0, 1], "validation": [1], "test": [2]}
    result = validate_matlab_ml_rl_v54_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v54_malformed_values_reject_without_raising() -> None:
    summary = deepcopy(_summary_v54())
    summary["matlab_ml_rl_v54_identity"]["prioritized_replay"]["importance_beta"] = [0.4]
    summary["matlab_ml_rl_v54_identity"]["feature_preprocess"]["feature_std"] = [[2.0], 0.5]
    result = validate_matlab_ml_rl_v54_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v54_numeric_sha256_values_are_rejected() -> None:
    summary = _summary_v54()
    numeric_digest = int("9" * 64)
    for row in summary["matlab_ml_rl_v54_identity"].values():
        for name in tuple(row):
            if name.endswith("_sha256"):
                row[name] = numeric_digest
    result = validate_matlab_ml_rl_v54_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v55_positive_identity_is_accepted() -> None:
    result = validate_matlab_ml_rl_v55_identity(_summary_v55())
    assert result and result["status"] == "ok"


def test_v55_frozen_public_counterfactuals_are_rejected() -> None:
    summary = deepcopy(_summary_v55())
    nstep = summary["matlab_ml_rl_v54_identity"]["nstep_return"]
    nstep.update({"result_gamma": 0.5, "result_rewards": [99.0], "result_terminal": True, "result_bootstrap_value": 0.0, "result_n_step_return": -1.0, "result_trajectory_id": "trajectory:stale", "result_policy_owner": "policy:stale"})
    cross_validation = summary["matlab_ml_rl_v54_identity"]["cross_validation"]
    cross_validation.update({"result_fold_id_per_sample": [0] * 6, "result_class_labels": ["class:A"] * 6, "result_preprocess_fit_rows": {"0": [0, 1]}, "result_rng_seed": 42, "result_fold_metrics": [1.0], "result_aggregate_metric": 1.0, "result_model_owner": "model:stale"})
    result = validate_matlab_ml_rl_v55_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v55_self_consistent_invalid_return_and_terminal_semantics_are_rejected() -> None:
    summary = deepcopy(_summary_v55())
    nstep = summary["matlab_ml_rl_v54_identity"]["nstep_return"]
    nstep["terminal"] = nstep["result_terminal"] = True
    result = validate_matlab_ml_rl_v55_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v55_self_consistent_cv_leakage_and_unstratified_folds_are_rejected() -> None:
    summary = deepcopy(_summary_v55())
    cross_validation = summary["matlab_ml_rl_v54_identity"]["cross_validation"]
    labels = ["class:A", "class:A", "class:B", "class:B", "class:A", "class:B"]
    fit_rows = {"0": [0, 1, 2, 3, 4, 5], "1": [0, 1, 4, 5], "2": [0, 1, 2, 3]}
    cross_validation["class_labels"] = cross_validation["result_class_labels"] = labels
    cross_validation["preprocess_fit_rows"] = cross_validation["result_preprocess_fit_rows"] = fit_rows
    result = validate_matlab_ml_rl_v55_identity(summary)
    assert result and result["status"] == "needs_attention"


def test_v55_numeric_digests_are_rejected() -> None:
    summary = deepcopy(_summary_v55())
    numeric_digest = int("1" * 64)
    identity = summary["matlab_ml_rl_v54_identity"]
    for contract_name in ("nstep_return", "cross_validation"):
        identity[contract_name]["result_sha256"] = numeric_digest
        identity[contract_name]["accepted_result_sha256"] = numeric_digest
    result = validate_matlab_ml_rl_v55_identity(summary)
    assert result and result["status"] == "needs_attention"
