"""Payload builders for test_matlab_generalization.py (not collected)."""

from __future__ import annotations

import math
from copy import deepcopy
from math import log10, sqrt

from test_regularized_trace_inverse_gate import _summary, _with_v22_identity


def _summary_v23():
    summary = _with_v22_identity(_summary())
    summary["parallel_pool_worker_path_device_rng_code_generation_identity"] = {
        "pool_generation": "parallel-pool-51",
        "worker_path_pool_generation": "parallel-pool-51",
        "device_pool_generation": "parallel-pool-51",
        "rng_pool_generation": "parallel-pool-51",
        "code_pool_generation": "parallel-pool-51",
        "result_pool_generation": "parallel-pool-51",
        "worker_ids": [1, 2, 3, 4],
        "result_worker_ids": [1, 2, 3, 4],
        "worker_code_paths": ["toolbox/a"] * 4,
        "result_worker_code_paths": ["toolbox/a"] * 4,
        "device_assignments": ["cpu:0", "cpu:1", "cpu:2", "cpu:3"],
        "result_device_assignments": ["cpu:0", "cpu:1", "cpu:2", "cpu:3"],
        "random_stream_seeds": [101, 202, 303, 404],
        "result_random_stream_seeds": [101, 202, 303, 404],
        "worker_code_sha256": "1" * 64,
        "result_worker_code_sha256": "1" * 64,
        "parallel_result_sha256": "2" * 64,
        "assembled_parallel_result_sha256": "2" * 64,
    }
    summary["autodiff_tape_variable_order_mesh_objective_generation_identity"] = {
        "tape_generation": "autodiff-tape-51",
        "variable_order_tape_generation": "autodiff-tape-51",
        "mesh_tape_generation": "autodiff-tape-51",
        "objective_scaling_tape_generation": "autodiff-tape-51",
        "primal_solve_tape_generation": "autodiff-tape-51",
        "gradient_result_tape_generation": "autodiff-tape-51",
        "variable_ids": ["radius", "thickness", "impedance"],
        "gradient_variable_ids": ["radius", "thickness", "impedance"],
        "mesh_sha256": "3" * 64,
        "gradient_mesh_sha256": "3" * 64,
        "objective_id": "radiated_power",
        "gradient_objective_id": "radiated_power",
        "objective_scale": 0.001,
        "gradient_objective_scale": 0.001,
        "primal_state_sha256": "4" * 64,
        "gradient_primal_state_sha256": "4" * 64,
        "gradient_table_sha256": "5" * 64,
        "reported_gradient_table_sha256": "5" * 64,
    }
    return summary


def _summary_v24():
    summary = _summary_v23()
    summary["fembem_trace_normal_interface_node_order_unit_generation_identity"] = {
        "coupling_generation": "fembem-101",
        "trace_coupling_generation": "fembem-101",
        "normal_coupling_generation": "fembem-101",
        "node_order_coupling_generation": "fembem-101",
        "unit_coupling_generation": "fembem-101",
        "operator_coupling_generation": "fembem-101",
        "result_coupling_generation": "fembem-101",
        "trace_orientation": "volume_to_boundary",
        "result_trace_orientation": "volume_to_boundary",
        "outward_normal_convention": "exterior_from_volume",
        "result_outward_normal_convention": "exterior_from_volume",
        "interface_node_ids": [1, 2, 3, 4],
        "result_interface_node_ids": [1, 2, 3, 4],
        "boundary_triangles": [[1, 2, 3], [1, 4, 2]],
        "result_boundary_triangles": [[1, 2, 3], [1, 4, 2]],
        "outward_normals": [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
        "result_outward_normals": [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
        "physical_units": {"pressure": "Pa", "normal_velocity": "m/s"},
        "result_physical_units": {"pressure": "Pa", "normal_velocity": "m/s"},
        "interface_mesh_sha256": "1" * 64,
        "result_interface_mesh_sha256": "1" * 64,
        "coupled_operator_sha256": "2" * 64,
        "result_coupled_operator_sha256": "2" * 64,
    }
    summary[
        "cq_contour_weight_startup_causality_window_result_generation_identity"
    ] = {
        "cq_generation": "cq-time-101",
        "contour_cq_generation": "cq-time-101",
        "weight_cq_generation": "cq-time-101",
        "startup_cq_generation": "cq-time-101",
        "causality_window_cq_generation": "cq-time-101",
        "time_grid_cq_generation": "cq-time-101",
        "result_cq_generation": "cq-time-101",
        "method": "BDF2",
        "result_method": "BDF2",
        "contour_points_ri": [[0.8, 0.0], [0.0, 0.8], [-0.8, 0.0], [0.0, -0.8]],
        "result_contour_points_ri": [[0.8, 0.0], [0.0, 0.8], [-0.8, 0.0], [0.0, -0.8]],
        "cq_weights_ri": [[1.5, 0.0], [-2.0, 0.0], [0.5, 0.0], [0.0, 0.0]],
        "result_cq_weights_ri": [[1.5, 0.0], [-2.0, 0.0], [0.5, 0.0], [0.0, 0.0]],
        "startup_weights_ri": [[1.0, 0.0], [-1.0, 0.0]],
        "result_startup_weights_ri": [[1.0, 0.0], [-1.0, 0.0]],
        "time_samples_s": [0.0, 0.001, 0.002, 0.003],
        "result_time_samples_s": [0.0, 0.001, 0.002, 0.003],
        "causality_window_s": [0.0, 0.003],
        "result_causality_window_s": [0.0, 0.003],
        "prehistory_norm": 0.0,
        "result_prehistory_norm": 0.0,
        "cq_result_sha256": "3" * 64,
        "reported_cq_result_sha256": "3" * 64,
    }
    return summary


def _summary_v25():
    summary = _summary_v24()
    summary[
        "hmatrix_block_tree_admissibility_permutation_tolerance_kernel_mesh_generation_identity"
    ] = {
        "hmatrix_generation": "hmatrix-201",
        "block_tree_hmatrix_generation": "hmatrix-201",
        "admissibility_hmatrix_generation": "hmatrix-201",
        "permutation_hmatrix_generation": "hmatrix-201",
        "tolerance_hmatrix_generation": "hmatrix-201",
        "kernel_hmatrix_generation": "hmatrix-201",
        "mesh_hmatrix_generation": "hmatrix-201",
        "result_hmatrix_generation": "hmatrix-201",
        "matrix_shape": [4, 4],
        "result_matrix_shape": [4, 4],
        "block_tree_sha256": "1" * 64,
        "result_block_tree_sha256": "1" * 64,
        "admissibility_rule": "diameter_le_eta_distance",
        "result_admissibility_rule": "diameter_le_eta_distance",
        "admissibility_eta": 1.5,
        "result_admissibility_eta": 1.5,
        "row_permutation": [2, 0, 3, 1],
        "result_row_permutation": [2, 0, 3, 1],
        "column_permutation": [1, 3, 0, 2],
        "result_column_permutation": [1, 3, 0, 2],
        "relative_tolerance": 1.0e-6,
        "result_relative_tolerance": 1.0e-6,
        "kernel_id": "helmholtz_single_layer_p1",
        "result_kernel_id": "helmholtz_single_layer_p1",
        "boundary_mesh_sha256": "2" * 64,
        "result_boundary_mesh_sha256": "2" * 64,
        "hmatrix_result_sha256": "3" * 64,
        "reported_hmatrix_result_sha256": "3" * 64,
    }
    summary[
        "ad_parameter_tape_material_operator_mesh_objective_primal_gradient_generation_identity"
    ] = {
        "ad_generation": "ad-gradient-201",
        "parameter_tape_ad_generation": "ad-gradient-201",
        "material_law_ad_generation": "ad-gradient-201",
        "operator_ad_generation": "ad-gradient-201",
        "mesh_ad_generation": "ad-gradient-201",
        "objective_ad_generation": "ad-gradient-201",
        "primal_ad_generation": "ad-gradient-201",
        "gradient_ad_generation": "ad-gradient-201",
        "result_ad_generation": "ad-gradient-201",
        "parameter_names": ["density", "bulk_modulus"],
        "result_parameter_names": ["density", "bulk_modulus"],
        "parameter_values": [1.2, 142000.0],
        "result_parameter_values": [1.2, 142000.0],
        "material_law_id": "linear_acoustic_fluid",
        "result_material_law_id": "linear_acoustic_fluid",
        "objective_id": "receiver_pressure_l2",
        "result_objective_id": "receiver_pressure_l2",
        "parameter_tape_sha256": "4" * 64,
        "result_parameter_tape_sha256": "4" * 64,
        "assembled_operator_sha256": "5" * 64,
        "result_assembled_operator_sha256": "5" * 64,
        "mesh_sha256": "6" * 64,
        "result_mesh_sha256": "6" * 64,
        "primal_solution_sha256": "7" * 64,
        "result_primal_solution_sha256": "7" * 64,
        "ad_gradient": [0.25, -0.004],
        "finite_difference_gradient": [0.25000001, -0.0040000001],
        "maximum_gradient_relative_error": 4.0e-8,
        "gradient_relative_tolerance": 1.0e-5,
        "gradient_result_sha256": "8" * 64,
        "reported_gradient_result_sha256": "8" * 64,
    }
    return summary


def _summary_v26():
    summary = _summary_v25()
    generation = "cq-transfer-301"
    summary[
        "cq_contour_radius_timestep_laplace_branch_transfer_operator_inverse_transform_generation_identity"
    ] = {
        "cq_generation": generation,
        "contour_cq_generation": generation,
        "timestep_cq_generation": generation,
        "laplace_branch_cq_generation": generation,
        "transfer_operator_cq_generation": generation,
        "inverse_transform_cq_generation": generation,
        "result_cq_generation": generation,
        "contour_radius": 0.95,
        "result_contour_radius": 0.95,
        "time_step_s": 1.0e-4,
        "result_time_step_s": 1.0e-4,
        "laplace_branch": "principal_sqrt_outgoing",
        "result_laplace_branch": "principal_sqrt_outgoing",
        "laplace_points_ri": [[100.0, 0.0], [80.0, 20.0], [60.0, 35.0], [80.0, -20.0]],
        "result_laplace_points_ri": [[100.0, 0.0], [80.0, 20.0], [60.0, 35.0], [80.0, -20.0]],
        "transfer_operator_id": "helmholtz_calderon_p1",
        "result_transfer_operator_id": "helmholtz_calderon_p1",
        "transfer_operator_sha256": "1" * 64,
        "result_transfer_operator_sha256": "1" * 64,
        "inverse_transform": "fft_conjugate_symmetric",
        "result_inverse_transform": "fft_conjugate_symmetric",
        "time_history_sha256": "2" * 64,
        "reported_time_history_sha256": "2" * 64,
    }
    generation = "fembem-coupling-301"
    summary[
        "fembem_trace_map_normal_material_wavenumber_coupling_matrix_mesh_generation_identity"
    ] = {
        "coupling_generation": generation,
        "trace_map_coupling_generation": generation,
        "normal_coupling_generation": generation,
        "material_coupling_generation": generation,
        "wavenumber_coupling_generation": generation,
        "matrix_coupling_generation": generation,
        "mesh_coupling_generation": generation,
        "result_coupling_generation": generation,
        "trace_map_sha256": "3" * 64,
        "result_trace_map_sha256": "3" * 64,
        "normal_orientation": "volume_outward",
        "result_normal_orientation": "volume_outward",
        "normal_field_sha256": "4" * 64,
        "result_normal_field_sha256": "4" * 64,
        "fluid_density_kg_m3": 1.2,
        "result_fluid_density_kg_m3": 1.2,
        "sound_speed_m_s": 343.0,
        "result_sound_speed_m_s": 343.0,
        "wavenumber_ri_m_inv": [18.318324511, 0.02],
        "result_wavenumber_ri_m_inv": [18.318324511, 0.02],
        "coupling_matrix_sha256": "5" * 64,
        "result_coupling_matrix_sha256": "5" * 64,
        "volume_mesh_sha256": "6" * 64,
        "result_volume_mesh_sha256": "6" * 64,
        "boundary_mesh_sha256": "7" * 64,
        "result_boundary_mesh_sha256": "7" * 64,
        "coupled_result_sha256": "8" * 64,
        "reported_coupled_result_sha256": "8" * 64,
    }
    return summary


def _summary_v27():
    summary = _summary_v26()
    generation = "adaptive-cq-311"
    summary[
        "cq_adaptive_contour_quadrature_order_startup_correction_error_estimator_restart_generation_identity"
    ] = {
        "cq_generation": generation,
        "contour_cq_generation": generation,
        "quadrature_order_cq_generation": generation,
        "startup_correction_cq_generation": generation,
        "error_estimator_cq_generation": generation,
        "restart_cq_generation": generation,
        "result_cq_generation": generation,
        "contour_family": "lubich_bdf2_circle",
        "result_contour_family": "lubich_bdf2_circle",
        "contour_radii": [0.82, 0.9, 0.95],
        "result_contour_radii": [0.82, 0.9, 0.95],
        "quadrature_orders": [16, 32, 64],
        "result_quadrature_orders": [16, 32, 64],
        "startup_correction": "bdf2_consistent_two_step",
        "result_startup_correction": "bdf2_consistent_two_step",
        "error_estimator": "successive_contour_l2_relative",
        "result_error_estimator": "successive_contour_l2_relative",
        "relative_tolerance": 1.0e-5,
        "result_relative_tolerance": 1.0e-5,
        "estimated_relative_errors": [2.0e-3, 1.5e-4, 8.0e-6],
        "result_estimated_relative_errors": [2.0e-3, 1.5e-4, 8.0e-6],
        "restart_step": 40,
        "result_restart_step": 40,
        "restart_state_sha256": "1" * 64,
        "loaded_restart_state_sha256": "1" * 64,
        "time_history_sha256": "2" * 64,
        "accepted_time_history_sha256": "2" * 64,
    }
    generation = "p1-fembem-311"
    summary[
        "p1_fembem_boundary_orientation_quadrature_singular_treatment_trace_matrix_mesh_generation_identity"
    ] = {
        "coupling_generation": generation,
        "boundary_orientation_coupling_generation": generation,
        "quadrature_coupling_generation": generation,
        "singular_treatment_coupling_generation": generation,
        "trace_coupling_generation": generation,
        "matrix_coupling_generation": generation,
        "mesh_coupling_generation": generation,
        "result_coupling_generation": generation,
        "fem_basis_order": 1,
        "bem_basis_order": 1,
        "volume_element": "tet",
        "boundary_element": "tri",
        "boundary_orientation": "volume_outward",
        "result_boundary_orientation": "volume_outward",
        "regular_quadrature": "triangle_degree_4",
        "result_regular_quadrature": "triangle_degree_4",
        "singular_treatment": "duffy_p1_galerkin",
        "result_singular_treatment": "duffy_p1_galerkin",
        "trace_shape": [48, 120],
        "result_trace_shape": [48, 120],
        "trace_matrix_sha256": "3" * 64,
        "result_trace_matrix_sha256": "3" * 64,
        "fem_matrix_sha256": "4" * 64,
        "result_fem_matrix_sha256": "4" * 64,
        "bem_matrix_sha256": "5" * 64,
        "result_bem_matrix_sha256": "5" * 64,
        "volume_mesh_sha256": "6" * 64,
        "result_volume_mesh_sha256": "6" * 64,
        "boundary_mesh_sha256": "7" * 64,
        "result_boundary_mesh_sha256": "7" * 64,
        "coupled_result_sha256": "8" * 64,
        "accepted_coupled_result_sha256": "8" * 64,
    }
    return summary


def _summary_v28():
    summary = _summary_v27()
    generation = "hmatrix-aca-321"
    summary[
        "hmatrix_aca_cluster_permutation_admissibility_rank_tolerance_kernel_mesh_result_generation_identity"
    ] = {
        "hmatrix_generation": generation,
        "cluster_hmatrix_generation": generation,
        "permutation_hmatrix_generation": generation,
        "admissibility_hmatrix_generation": generation,
        "rank_hmatrix_generation": generation,
        "tolerance_hmatrix_generation": generation,
        "kernel_hmatrix_generation": generation,
        "mesh_hmatrix_generation": generation,
        "result_hmatrix_generation": generation,
        "cluster_permutation": [3, 1, 4, 2],
        "result_cluster_permutation": [3, 1, 4, 2],
        "admissibility_rule": "eta-weak",
        "result_admissibility_rule": "eta-weak",
        "admissibility_eta": 2.0,
        "result_admissibility_eta": 2.0,
        "aca_rank": 8,
        "result_aca_rank": 8,
        "relative_tolerance": 1.0e-6,
        "result_relative_tolerance": 1.0e-6,
        "kernel": "helmholtz-single-layer-p1",
        "result_kernel": "helmholtz-single-layer-p1",
        "cluster_tree_sha256": "1" * 64,
        "loaded_cluster_tree_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "result_sha256": "3" * 64,
        "accepted_result_sha256": "3" * 64,
    }
    generation = "calderon-cq-321"
    summary[
        "calderon_cq_operator_v_k_trace_normal_frequency_grid_inverse_transform_mesh_result_generation_identity"
    ] = {
        "calderon_generation": generation,
        "v_calderon_generation": generation,
        "k_calderon_generation": generation,
        "trace_calderon_generation": generation,
        "normal_calderon_generation": generation,
        "frequency_calderon_generation": generation,
        "inverse_calderon_generation": generation,
        "mesh_calderon_generation": generation,
        "result_calderon_generation": generation,
        "v_operator_sha256": "4" * 64,
        "result_v_operator_sha256": "4" * 64,
        "k_operator_sha256": "5" * 64,
        "result_k_operator_sha256": "5" * 64,
        "trace_basis": "p1-nodal-boundary-trace",
        "result_trace_basis": "p1-nodal-boundary-trace",
        "trace_shape": [48, 120],
        "result_trace_shape": [48, 120],
        "boundary_normal": "volume-outward",
        "result_boundary_normal": "volume-outward",
        "laplace_frequency_ri": [[10.0, 0.0], [10.0, 20.0], [10.0, 40.0]],
        "result_laplace_frequency_ri": [[10.0, 0.0], [10.0, 20.0], [10.0, 40.0]],
        "inverse_transform": "bdf2-cq-ifft-real",
        "result_inverse_transform": "bdf2-cq-ifft-real",
        "boundary_mesh_sha256": "6" * 64,
        "result_boundary_mesh_sha256": "6" * 64,
        "result_sha256": "7" * 64,
        "accepted_result_sha256": "7" * 64,
    }
    return summary


def _summary_v29():
    summary = _summary_v28()
    generation = "near-singular-331"
    summary[
        "bem_near_singular_quadrature_distance_element_size_adaptive_order_reference_result_generation_identity"
    ] = {
        "quadrature_generation": generation,
        "target_quadrature_generation": generation,
        "geometry_quadrature_generation": generation,
        "order_quadrature_generation": generation,
        "map_quadrature_generation": generation,
        "kernel_quadrature_generation": generation,
        "reference_quadrature_generation": generation,
        "mesh_quadrature_generation": generation,
        "result_quadrature_generation": generation,
        "target_distance_m": 2.0e-4,
        "result_target_distance_m": 2.0e-4,
        "element_size_m": 1.0e-3,
        "result_element_size_m": 1.0e-3,
        "distance_size_ratio": 0.2,
        "result_distance_size_ratio": 0.2,
        "adaptive_order": 16,
        "result_adaptive_order": 16,
        "quadrature_rule": "adaptive-duffy-p1",
        "result_quadrature_rule": "adaptive-duffy-p1",
        "coordinate_map": "target-aligned-barycentric",
        "result_coordinate_map": "target-aligned-barycentric",
        "kernel": "helmholtz-single-layer-p1",
        "result_kernel": "helmholtz-single-layer-p1",
        "reference_integral_ri": [0.125, -0.03125],
        "computed_integral_ri": [0.1250000001, -0.0312499999],
        "relative_error": 1.1e-9,
        "relative_tolerance": 1.0e-7,
        "element_mesh_sha256": "1" * 64,
        "result_element_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "fembem-energy-331"
    summary[
        "fembem_energy_flux_reciprocity_interface_trace_orientation_frequency_incident_result_generation_identity"
    ] = {
        "coupling_generation": generation,
        "trace_coupling_generation": generation,
        "normal_coupling_generation": generation,
        "frequency_coupling_generation": generation,
        "incident_coupling_generation": generation,
        "reciprocity_coupling_generation": generation,
        "energy_coupling_generation": generation,
        "result_coupling_generation": generation,
        "interface_trace_basis": "p1-nodal-boundary-trace",
        "result_interface_trace_basis": "p1-nodal-boundary-trace",
        "interface_trace_shape": [48, 120],
        "result_interface_trace_shape": [48, 120],
        "normal_orientation": "volume-outward",
        "result_normal_orientation": "volume-outward",
        "frequency_hz": 800.0,
        "result_frequency_hz": 800.0,
        "incident_field_sha256": "3" * 64,
        "result_incident_field_sha256": "3" * 64,
        "reciprocity_pair_ids": ["source-a/receiver-b", "source-b/receiver-a"],
        "result_reciprocity_pair_ids": ["source-a/receiver-b", "source-b/receiver-a"],
        "reciprocity_values_ri": [[0.5, -0.1], [0.5000000002, -0.1000000001]],
        "result_reciprocity_values_ri": [[0.5, -0.1], [0.5000000002, -0.1000000001]],
        "reciprocity_relative_error": 4.4e-10,
        "reciprocity_relative_tolerance": 1.0e-7,
        "fem_outward_power_w": 1.25,
        "bem_radiated_power_w": 1.249999999,
        "energy_flux_relative_error": 8.0e-10,
        "energy_flux_relative_tolerance": 1.0e-7,
        "coupled_result_sha256": "4" * 64,
        "accepted_coupled_result_sha256": "4" * 64,
    }
    return summary


def _summary_v30():
    summary = _summary_v29()
    generation = "hmatrix-recompress-341"
    summary[
        "hmatrix_recompression_svd_tolerance_norm_rank_permutation_operator_mesh_result_identity"
    ] = {
        "hmatrix_generation": generation,
        "svd_hmatrix_generation": generation,
        "tolerance_hmatrix_generation": generation,
        "rank_hmatrix_generation": generation,
        "permutation_hmatrix_generation": generation,
        "operator_hmatrix_generation": generation,
        "mesh_hmatrix_generation": generation,
        "result_hmatrix_generation": generation,
        "svd_basis": "euclidean-orthonormal",
        "result_svd_basis": "euclidean-orthonormal",
        "tolerance": 1.0e-6,
        "result_tolerance": 1.0e-6,
        "tolerance_norm": "spectral-relative",
        "result_tolerance_norm": "spectral-relative",
        "block_ranks_before": [12, 10, 8],
        "result_block_ranks_before": [12, 10, 8],
        "block_ranks_after": [6, 5, 4],
        "result_block_ranks_after": [6, 5, 4],
        "row_permutation": [2, 0, 1],
        "result_row_permutation": [2, 0, 1],
        "column_permutation": [1, 2, 0],
        "result_column_permutation": [1, 2, 0],
        "operator_relative_error": 5.0e-7,
        "result_operator_relative_error": 5.0e-7,
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "cq-restart-341"
    summary[
        "cq_restart_block_history_startup_weight_time_index_sample_contour_operator_result_identity"
    ] = {
        "cq_generation": generation,
        "block_cq_generation": generation,
        "history_cq_generation": generation,
        "startup_cq_generation": generation,
        "time_cq_generation": generation,
        "sample_cq_generation": generation,
        "owner_cq_generation": generation,
        "result_cq_generation": generation,
        "block_size": 16,
        "result_block_size": 16,
        "completed_block_ids": [0, 1, 2],
        "result_completed_block_ids": [0, 1, 2],
        "history_sample_count": 48,
        "result_history_sample_count": 48,
        "startup_weights_ri": [[1.0, 0.0], [0.5, -0.1]],
        "result_startup_weights_ri": [[1.0, 0.0], [0.5, -0.1]],
        "restart_time_index": 48,
        "result_restart_time_index": 48,
        "total_sample_count": 128,
        "result_total_sample_count": 128,
        "contour_owner_sha256": "3" * 64,
        "result_contour_owner_sha256": "3" * 64,
        "operator_owner_sha256": "4" * 64,
        "result_operator_owner_sha256": "4" * 64,
        "history_sha256": "5" * 64,
        "loaded_history_sha256": "5" * 64,
        "result_sha256": "6" * 64,
        "accepted_result_sha256": "6" * 64,
    }
    return summary


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


def _summary_v32():
    payload = deepcopy(_summary_v31())
    generation = "adaptive-cq-361"
    payload[
        "adaptive_cq_timestep_contour_history_interpolation_error_restart_operator_mesh_result_identity"
    ] = {
        "cq_generation": generation,
        **{
            key: generation
            for key in (
                "timestep_cq_generation",
                "contour_cq_generation",
                "history_cq_generation",
                "error_cq_generation",
                "restart_cq_generation",
                "operator_cq_generation",
                "mesh_cq_generation",
                "result_cq_generation",
            )
        },
        "timestep_schedule_s": [1.0e-5, 1.0e-5, 5.0e-6, 5.0e-6],
        "result_timestep_schedule_s": [1.0e-5, 1.0e-5, 5.0e-6, 5.0e-6],
        "contour_rebuild_indices": [0, 2],
        "result_contour_rebuild_indices": [0, 2],
        "history_interpolation": "barycentric_causal",
        "result_history_interpolation": "barycentric_causal",
        "local_error_estimates": [1.0e-4, 8.0e-5, 2.0e-5, 1.0e-5],
        "result_local_error_estimates": [1.0e-4, 8.0e-5, 2.0e-5, 1.0e-5],
        "local_error_tolerance": 1.0e-3,
        "result_local_error_tolerance": 1.0e-3,
        "restart_index": 2,
        "result_restart_index": 2,
        "operator_owner_sha256": "1" * 64,
        "result_operator_owner_sha256": "1" * 64,
        "history_sha256": "2" * 64,
        "loaded_history_sha256": "2" * 64,
        "mesh_sha256": "3" * 64,
        "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    generation = "modal-fembem-361"
    payload[
        "fembem_modal_transient_mass_damping_initial_projection_truncation_energy_mesh_history_result_identity"
    ] = {
        "modal_generation": generation,
        **{
            key: generation
            for key in (
                "mass_modal_generation",
                "damping_modal_generation",
                "initial_modal_generation",
                "truncation_modal_generation",
                "energy_modal_generation",
                "mesh_modal_generation",
                "history_modal_generation",
                "result_modal_generation",
            )
        },
        "mass_normalization": "M_orthonormal",
        "result_mass_normalization": "M_orthonormal",
        "damping_model": "rayleigh",
        "result_damping_model": "rayleigh",
        "rayleigh_coefficients": [0.01, 1.0e-5],
        "result_rayleigh_coefficients": [0.01, 1.0e-5],
        "initial_displacement_projection": [1.0, 0.2, 0.0],
        "result_initial_displacement_projection": [1.0, 0.2, 0.0],
        "initial_velocity_projection": [0.0, 0.0, 0.0],
        "result_initial_velocity_projection": [0.0, 0.0, 0.0],
        "modal_count": 3,
        "result_modal_count": 3,
        "truncation_frequency_hz": 1200.0,
        "result_truncation_frequency_hz": 1200.0,
        "initial_energy_j": 0.5,
        "result_initial_energy_j": 0.5,
        "radiated_energy_j": 0.2,
        "result_radiated_energy_j": 0.2,
        "dissipated_energy_j": 0.1,
        "result_dissipated_energy_j": 0.1,
        "final_energy_j": 0.2,
        "result_final_energy_j": 0.2,
        "mesh_sha256": "5" * 64,
        "result_mesh_sha256": "5" * 64,
        "time_history_owner": "fembem/case-361/modal-transient",
        "accepted_time_history_owner": "fembem/case-361/modal-transient",
        "result_sha256": "6" * 64,
        "accepted_result_sha256": "6" * 64,
    }
    return payload


def _summary_v33():
    payload = deepcopy(_summary_v32())
    generation = "calderon-p1-381"
    payload[
        "calderon_projector_p1_v_k_kt_w_mass_duality_normal_quadrature_mesh_owner_result_identity"
    ] = {
        "calderon_generation": generation,
        **{
            key: generation
            for key in (
                "space_calderon_generation",
                "operator_calderon_generation",
                "mass_calderon_generation",
                "normal_calderon_generation",
                "quadrature_calderon_generation",
                "mesh_calderon_generation",
                "projector_calderon_generation",
                "owner_calderon_generation",
                "result_calderon_generation",
            )
        },
        "trial_space": "P1",
        "result_trial_space": "P1",
        "test_space": "P1",
        "result_test_space": "P1",
        "projector_convention": "interior_calderon_outward",
        "result_projector_convention": "interior_calderon_outward",
        "block_order": ["dirichlet", "neumann"],
        "result_block_order": ["dirichlet", "neumann"],
        "v_sign": -1,
        "result_v_sign": -1,
        "k_sign": 1,
        "result_k_sign": 1,
        "kt_sign": -1,
        "result_kt_sign": -1,
        "w_sign": -1,
        "result_w_sign": -1,
        "mass_duality_residual": 2.0e-11,
        "result_mass_duality_residual": 2.0e-11,
        "mass_duality_tolerance": 1.0e-8,
        "result_mass_duality_tolerance": 1.0e-8,
        "normal_orientation": "outward",
        "result_normal_orientation": "outward",
        "singular_quadrature": "duffy_principal_value_p1",
        "result_singular_quadrature": "duffy_principal_value_p1",
        "projector_residual": 5.0e-10,
        "result_projector_residual": 5.0e-10,
        "projector_tolerance": 1.0e-8,
        "result_projector_tolerance": 1.0e-8,
        "operator_sha256": "1" * 64,
        "result_operator_sha256": "1" * 64,
        "mass_sha256": "2" * 64,
        "result_mass_sha256": "2" * 64,
        "boundary_mesh_sha256": "3" * 64,
        "result_boundary_mesh_sha256": "3" * 64,
        "result_owner": "fembem/calderon/case-381",
        "accepted_result_owner": "fembem/calderon/case-381",
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    generation = "cq-physical-381"
    transfer_samples = [
        [1.0, 0.0],
        [0.8, 0.2],
        [0.5, 0.1],
        [0.5, -0.1],
        [0.8, -0.2],
    ]
    payload[
        "cq_symbol_contour_transfer_conjugate_causal_ifft_parseval_passivity_timestep_operator_result_identity"
    ] = {
        "cq_generation": generation,
        **{
            key: generation
            for key in (
                "symbol_cq_generation",
                "contour_cq_generation",
                "transfer_cq_generation",
                "symmetry_cq_generation",
                "causality_cq_generation",
                "parseval_cq_generation",
                "passivity_cq_generation",
                "timestep_cq_generation",
                "operator_cq_generation",
                "result_cq_generation",
            )
        },
        "multistep_symbol": "BDF2",
        "result_multistep_symbol": "BDF2",
        "symbol_coefficients": [1.5, -2.0, 0.5],
        "result_symbol_coefficients": [1.5, -2.0, 0.5],
        "contour_radius": 0.92,
        "result_contour_radius": 0.92,
        "transfer_samples_ri": transfer_samples,
        "result_transfer_samples_ri": [list(row) for row in transfer_samples],
        "time_response": [0.0, 0.2, 0.1, 0.05],
        "result_time_response": [0.0, 0.2, 0.1, 0.05],
        "negative_time_energy": 0.0,
        "result_negative_time_energy": 0.0,
        "time_domain_work": 0.25,
        "result_time_domain_work": 0.25,
        "frequency_domain_work": 0.25,
        "result_frequency_domain_work": 0.25,
        "parseval_tolerance": 1.0e-10,
        "result_parseval_tolerance": 1.0e-10,
        "minimum_real_transfer": 0.5,
        "result_minimum_real_transfer": 0.5,
        "passivity_sign": "nonnegative_real_transfer",
        "result_passivity_sign": "nonnegative_real_transfer",
        "timestep_s": 1.0e-5,
        "result_timestep_s": 1.0e-5,
        "operator_family": "p1_calderon_bem",
        "result_operator_family": "p1_calderon_bem",
        "operator_sha256": "5" * 64,
        "result_operator_sha256": "5" * 64,
        "result_owner": "fembem/cq/case-381",
        "accepted_result_owner": "fembem/cq/case-381",
        "result_sha256": "6" * 64,
        "accepted_result_sha256": "6" * 64,
    }
    return payload


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


def _summary_v37():
    payload = deepcopy(_summary_v36())
    generation = "hmatrix-benchmark-513"
    mirrored = {
        "boundary_unknown_count": [200, 400, 800],
        "dense_reference_relative_error": [8.0e-4, 6.0e-4, 4.0e-4],
        "relative_tolerance": 1.0e-3,
        "maximum_block_rank": [8, 12, 16],
        "dense_memory_bytes": [640000, 2560000, 10240000],
        "hmatrix_memory_bytes": [80000, 190000, 450000],
        "memory_complexity_exponent": 1.25,
        "rank_complexity_exponent": 0.5,
        "boundary_mesh_sha256": "1" * 64,
    }
    payload["hmatrix_dense_reference_error_tolerance_rank_memory_complexity_mesh_operator_benchmark_result_identity"] = {
        "hmatrix_generation": generation,
        **{key: generation for key in (
            "dense_generation", "tolerance_generation", "rank_generation", "memory_generation",
            "complexity_generation", "mesh_generation", "operator_generation", "benchmark_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "operator_owner": "hmatrix/operator-513", "accepted_operator_owner": "hmatrix/operator-513",
        "benchmark_owner": "hmatrix/benchmark-513", "accepted_benchmark_owner": "hmatrix/benchmark-513",
        "hmatrix_result_sha256": "2" * 64, "accepted_hmatrix_result_sha256": "2" * 64,
    }
    generation = "multifrequency-adjoint-513"
    mirrored = {
        "frequency_hz": [100.0, 200.0, 400.0], "frequency_weights": [0.2, 0.3, 0.5],
        "objective_complex": [[1.0, 0.1], [2.0, 0.2], [3.0, 0.3]],
        "weighted_objective_complex": [2.3, 0.23], "quadrature_order": [4, 4, 6],
        "trace_node_map": [1, 2, 3], "frequency_gradient": [2.0, 4.0, 6.0],
        "accumulated_gradient": 4.6, "finite_difference_gradient": 4.600001,
        "gradient_relative_tolerance": 1.0e-5, "fembem_mesh_sha256": "3" * 64,
    }
    payload["multifrequency_fembem_adjoint_weight_objective_quadrature_trace_gradient_fd_mesh_owner_result_identity"] = {
        "adjoint_generation": generation,
        **{key: generation for key in (
            "frequency_generation", "weight_generation", "objective_generation", "quadrature_generation",
            "trace_generation", "gradient_generation", "fd_generation", "mesh_generation", "owner_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "adjoint_owner": "fembem/adjoint-513", "accepted_adjoint_owner": "fembem/adjoint-513",
        "adjoint_result_sha256": "4" * 64, "accepted_adjoint_result_sha256": "4" * 64,
    }
    return payload


def _summary_v38():
    payload = deepcopy(_summary_v37())
    generation = "simp-topology-614"
    density = [0.4, 0.5, 0.6, 0.5]
    filter_matrix = [
        [0.75, 0.25, 0.0, 0.0],
        [0.25, 0.5, 0.25, 0.0],
        [0.0, 0.25, 0.5, 0.25],
        [0.0, 0.0, 0.25, 0.75],
    ]
    filtered = [
        sum(weight * value for weight, value in zip(row, density))
        for row in filter_matrix
    ]
    beta = [1.0, 2.0, 4.0, 8.0]
    eta = 0.5
    denominator = math.tanh(beta[-1] * eta) + math.tanh(beta[-1] * (1.0 - eta))
    projected = [
        (math.tanh(beta[-1] * eta) + math.tanh(beta[-1] * (value - eta)))
        / denominator
        for value in filtered
    ]
    mirrored = {
        "design_density": density,
        "density_filter_matrix": filter_matrix,
        "filtered_density": filtered,
        "projection_beta_continuation": beta,
        "projection_eta": eta,
        "projected_density": projected,
        "volume_fraction": sum(projected) / len(projected),
        "volume_fraction_limit": 0.5,
        "compliance": 12.0,
        "adjoint_compliance_gradient": [-0.25] * 4,
        "finite_difference_compliance_gradient": [-0.2500001] * 4,
        "volume_gradient": [0.25] * 4,
        "volume_lagrange_multiplier": 1.0,
        "kkt_stationarity_residual": 0.0,
        "gradient_relative_tolerance": 1.0e-5,
        "topology_mesh_sha256": "1" * 64,
    }
    payload["simp_topology_density_filter_projection_volume_compliance_adjoint_fd_kkt_mesh_owner_result_identity"] = {
        "topology_generation": generation,
        **{key: generation for key in (
            "density_generation", "filter_generation", "projection_generation",
            "volume_generation", "compliance_generation", "adjoint_generation",
            "fd_generation", "kkt_generation", "mesh_generation",
            "owner_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "topology_owner": "optimization/topology-614",
        "accepted_topology_owner": "optimization/topology-614",
        "topology_result_sha256": "2" * 64,
        "accepted_topology_result_sha256": "2" * 64,
    }

    generation = "fembem-reduction-614"
    mirrored = {
        "full_order": 4,
        "reduced_order": 2,
        "trial_projection_basis": [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0]],
        "test_projection_basis": [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0]],
        "biorthogonality_gram": [[1.0, 0.0], [0.0, 1.0]],
        "full_model_poles": [[-1.0, 0.0], [-2.0, 0.0], [-3.0, 0.0], [-4.0, 0.0]],
        "reduced_model_poles": [[-1.1, 0.0], [-2.1, 0.0]],
        "minimum_passivity_eigenvalue": 0.1,
        "matched_moments_full": [[1.0, 0.0], [0.5, 0.0]],
        "matched_moments_reduced": [[1.0, 0.0], [0.5, 0.0]],
        "frequency_hz": [100.0, 200.0, 400.0],
        "frequency_response_relative_error": [0.01, 0.005, 0.002],
        "maximum_frequency_response_relative_error": 0.02,
        "reduction_mesh_sha256": "3" * 64,
    }
    payload["fembem_model_reduction_projection_order_stability_passivity_moment_frequency_error_full_mesh_owner_result_identity"] = {
        "reduction_generation": generation,
        **{key: generation for key in (
            "projection_generation", "order_generation", "stability_generation",
            "passivity_generation", "moment_generation", "frequency_generation",
            "error_generation", "full_model_generation", "mesh_generation",
            "owner_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "full_model_owner": "fembem/full-614",
        "accepted_full_model_owner": "fembem/full-614",
        "reduction_owner": "fembem/reduced-614",
        "accepted_reduction_owner": "fembem/reduced-614",
        "reduction_result_sha256": "4" * 64,
        "accepted_reduction_result_sha256": "4" * 64,
    }
    return payload


_NONLINEAR_KEY = (
    "nonlinear_fem_newton_residual_consistent_tangent_linesearch_step_energy_"
    "mesh_owner_result_identity"
)


_CQ_KEY_V39 = (
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
    payload[_CQ_KEY_V39] = {
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


_CQ_KEY_V41 = (
    "cq_acoustic_causality_passivity_timestep_ztransform_energy_history_mesh_"
    "owner_result_identity"
)


_HMATRIX_KEY = (
    "hmatrix_admissibility_cluster_rank_tolerance_matvec_error_memory_mesh_"
    "owner_result_identity"
)


def _summary_v41() -> dict:
    payload = deepcopy(_summary_v40())
    generation = "cq-acoustic-731"
    mirrored = {
        "multistep_method": "bdf2",
        "time_step_s": 2.5e-4,
        "z_transform_radius": 0.94,
        "multistep_symbol_samples": [[0.0, 0.0], [0.5, 0.2], [1.5, 0.0]],
        "laplace_frequency_samples_rad_s": [[20.0, 0.0], [35.0, 80.0], [60.0, 0.0]],
        "excitation_history": [0.0, 1.0, 0.4, 0.1, 0.0],
        "pressure_history": [0.0, 0.2, 0.3, 0.15, 0.04],
        "causal_prefix_length": 1,
        "minimum_passivity_real_part": 0.015,
        "boundary_work_j": 0.012,
        "radiated_energy_j": 0.010,
        "dissipated_energy_j": 0.002,
        "energy_balance_residual_j": 0.0,
        "energy_balance_tolerance_j": 1.0e-8,
        "boundary_mesh_sha256": "1" * 64,
    }
    payload[_CQ_KEY_V41] = {
        "cq_generation": generation,
        **{key: generation for key in (
            "multistep_generation", "timestep_generation", "ztransform_generation",
            "frequency_generation", "history_generation", "passivity_generation",
            "energy_generation", "mesh_generation", "owner_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "cq_owner": "acoustic/cq-731",
        "accepted_cq_owner": "acoustic/cq-731",
        "cq_result_sha256": "2" * 64,
        "accepted_cq_result_sha256": "2" * 64,
    }

    generation = "hmatrix-731"
    mirrored = {
        "cluster_leaf_size": 32,
        "cluster_permutation": [3, 1, 4, 2],
        "admissibility_eta": 2.0,
        "block_partition": [["low_rank", "dense"], ["dense", "low_rank"]],
        "numerical_ranks": [[4, 0], [0, 3]],
        "compression_relative_tolerance": 1.0e-5,
        "measured_matvec_relative_error": 4.0e-6,
        "dense_memory_bytes": 131072,
        "compressed_memory_bytes": 32768,
        "boundary_mesh_sha256": "3" * 64,
    }
    payload[_HMATRIX_KEY] = {
        "hmatrix_generation": generation,
        **{key: generation for key in (
            "cluster_generation", "admissibility_generation", "partition_generation",
            "rank_generation", "tolerance_generation", "matvec_generation",
            "memory_generation", "mesh_generation", "owner_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "hmatrix_owner": "acoustic/hmatrix-731",
        "accepted_hmatrix_owner": "acoustic/hmatrix-731",
        "hmatrix_result_sha256": "4" * 64,
        "accepted_hmatrix_result_sha256": "4" * 64,
    }
    return payload


_LOW_KEY = (
    "lowfrequency_bem_stabilized_kernel_staticlimit_condition_charge_energy_"
    "boundarymesh_result_identity"
)


_DUCT_KEY = (
    "duct_scattering_mode_pressure_transmission_reflection_loss_power_mesh_"
    "result_identity"
)


def _summary_v42() -> dict:
    payload = deepcopy(_summary_v41())
    generation = "low-frequency-bem-842"
    values = {
        "frequency_hz": [1.0e-3, 1.0e-2, 1.0e-1],
        "stabilized_kernel": "static_dynamic_split_p1",
        "dynamic_kernel_correction": [1.0e-6, 1.0e-5, 1.0e-4],
        "static_limit_residual": [1.0e-8, 1.0e-7, 1.0e-6],
        "condition_estimate": [100.0, 80.0, 60.0],
        "boundary_charge_c": [1.0e-9, -1.0e-9],
        "boundary_potential_v": [2.0, -2.0],
        "potential_energy_j": 2.0e-9,
    }
    payload[_LOW_KEY] = {
        "lowfrequency_bem_generation": generation,
        **{key: generation for key in (
            "frequency_generation", "kernel_generation", "staticlimit_generation",
            "condition_generation", "charge_generation", "energy_generation",
            "mesh_generation", "result_generation",
        )},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "boundary_mesh_owner": "boundary-mesh:low-frequency-bem-842",
        "accepted_boundary_mesh_owner": "boundary-mesh:low-frequency-bem-842",
        "bem_result_sha256": "1" * 64,
        "accepted_bem_result_sha256": "1" * 64,
    }

    generation = "duct-scattering-842"
    values = {
        "frequency_hz": 1000.0,
        "incident_mode_pressure": [1.0, 0.0],
        "reflected_mode_pressure": [0.2, 0.0],
        "transmitted_mode_pressure": [sqrt(0.96), 0.0],
        "incident_modal_power_w": 1.0,
        "reflected_modal_power_w": 0.04,
        "transmitted_modal_power_w": 0.96,
        "dissipated_power_w": 0.0,
        "transmission_loss_db": -10.0 * log10(0.96),
        "modal_power_balance_residual_w": 0.0,
    }
    payload[_DUCT_KEY] = {
        "duct_generation": generation,
        **{key: generation for key in (
            "frequency_generation", "mode_generation", "pressure_generation",
            "reflection_generation", "transmission_generation", "loss_generation",
            "power_generation", "mesh_generation", "result_generation",
        )},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "mesh_owner": "mesh:duct-scattering-842",
        "accepted_mesh_owner": "mesh:duct-scattering-842",
        "duct_result_sha256": "2" * 64,
        "accepted_duct_result_sha256": "2" * 64,
    }
    return payload


def _identity_v44() -> dict[str, object]:
    return {
        "supervised": {
            "schema": "cae-ai-lab.matlab-ml-rl-result.v2",
            "task": "regression",
            "matlab_release": "R2026a",
            "session_owner": "matlab:shared-session-v44",
            "random_seed": 844,
            "split_id": "holdout:ml-rl-v44",
            "training_ids": ["train-01", "train-02"],
            "evaluation_ids": ["eval-01"],
            "training_evaluation_disjoint": True,
            "validation_metric": {"name": "rmse", "value": 0.1, "units": "normalized"},
            "result_sha256": "5" * 64,
            "accepted_result_sha256": "5" * 64,
            "timing_breakdown_s": {"train": 0.2, "evaluate": 0.1},
        },
        "reinforcement_learning": {
            "schema": "cae-ai-lab.matlab-ml-rl-result.v2",
            "task": "reinforcement_learning",
            "matlab_release": "R2026a",
            "session_owner": "matlab:shared-session-v44",
            "random_seed": 845,
            "environment_id": "cae:rl-design-control-v44",
            "evaluation_environment_id": "cae:rl-design-control-v44",
            "training_episodes": 100,
            "evaluation_episodes": 20,
            "evaluation_seed": 1845,
            "exploration_during_evaluation": False,
            "training_evaluation_disjoint": True,
            "evaluation_mean_return": 1.5,
            "evaluation_std_return": 0.1,
            "result_sha256": "6" * 64,
            "accepted_result_sha256": "6" * 64,
            "evaluation_result_sha256": "7" * 64,
            "accepted_evaluation_result_sha256": "7" * 64,
            "timing_breakdown_s": {"train": 0.3, "evaluate": 0.2},
        },
    }


def _summary_v44() -> dict[str, object]:
    value = deepcopy(_summary_v42())
    value["matlab_ml_rl_v44_identity"] = _identity_v44()
    return value


def _identity_v45() -> dict[str, object]:
    return {
        "matlab_ml_rl_v45_identity": {
            "supervised": {"schema": "cae-ai-lab.matlab-ml-rl-result.v3", "datastore_id": "ds:845", "split_generation": "split:845", "result_split_generation": "split:845", "normalization_fit_scope": "training_only", "result_normalization_fit_scope": "training_only", "hyperparameter_selection_source": "training_cross_validation", "result_hyperparameter_selection_source": "training_cross_validation", "holdout_id": "holdout:845", "training_ids": ["train:1"], "model_card_release": "R2026a", "release_id": "R2026a", "result_release_id": "R2026a", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64},
            "reinforcement_learning": {"replay_buffer_generation": "replay:845", "result_replay_buffer_generation": "replay:845", "termination_semantics": "environment_defined", "result_termination_semantics": "environment_defined", "discount_factor": 0.99, "result_discount_factor": 0.99, "policy_seed": 845, "evaluation_mode": "greedy_no_exploration", "result_evaluation_mode": "greedy_no_exploration", "release_id": "R2026a", "result_release_id": "R2026a", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64},
            "agentic_toolkit": {"capability_route": "matlab_tool", "result_capability_route": "matlab_tool", "consent_recorded": True, "result_consent_recorded": True, "tool_arguments_sha256": "c" * 64, "result_tool_arguments_sha256": "c" * 64, "session_detection": "existing_shared_matlab", "result_session_detection": "existing_shared_matlab", "release_id": "R2026a", "result_release_id": "R2026a", "owner": "matlab:shared", "result_owner": "matlab:shared"},
            "mlrl_checkpoint": {"checkpoint_generation": "checkpoint:845", "result_checkpoint_generation": "checkpoint:845", "datastore_state": "replayed", "result_datastore_state": "replayed", "optimizer": "adam", "result_optimizer": "adam", "discount_factor": 0.99, "result_discount_factor": 0.99, "evaluation_id": "eval:845", "release_id": "R2026a", "result_release_id": "R2026a", "owner": "matlab:shared", "result_owner": "matlab:shared", "result_sha256": "d" * 64, "accepted_result_sha256": "d" * 64},
        }
    }


def _summary_v46():
    return {"matlab_ml_rl_v46_identity": {
        "supervised": {"nonfinite_policy": "drop_with_count", "result_nonfinite_policy": "drop_with_count", "nonfinite_input_count": 0, "result_nonfinite_input_count": 0, "split_generation": "split:test", "result_split_generation": "split:test", "normalization_fit_scope": "training_only", "result_normalization_fit_scope": "training_only", "worker_seed": 846, "result_worker_seed": 846, "restart_state": "fresh_training", "result_restart_state": "fresh_training", "release_id": "R2026a", "result_release_id": "R2026a", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64},
        "reinforcement_learning": {"episode_timeout_steps": 200, "result_episode_timeout_steps": 200, "termination_semantics": "environment_defined", "result_termination_semantics": "environment_defined", "exploration_mode": "training_only", "result_exploration_mode": "training_only", "checkpoint_generation": "checkpoint:test", "result_checkpoint_generation": "checkpoint:test", "evaluation_mode": "greedy_no_exploration", "result_evaluation_mode": "greedy_no_exploration", "release_id": "R2026a", "result_release_id": "R2026a", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64},
        "agentic_toolkit": {"argument_schema": "json_object", "result_argument_schema": "json_object", "argument_shape_valid": True, "result_argument_shape_valid": True, "session_detection": "existing_shared_matlab", "result_session_detection": "existing_shared_matlab", "timeout_s": 90.0, "result_timeout_s": 90.0, "error_class": "none", "result_error_class": "none", "tool_arguments_sha256": "c" * 64, "result_tool_arguments_sha256": "c" * 64, "release_id": "R2026a", "result_release_id": "R2026a", "owner": "matlab:shared", "result_owner": "matlab:shared"},
        "mlrl_checkpoint": {"worker_seed": 846, "result_worker_seed": 846, "checkpoint_order": [0, 1, 2], "result_checkpoint_order": [0, 1, 2], "optimizer_state": "adam_ready", "result_optimizer_state": "adam_ready", "checkpoint_generation": "checkpoint:test", "result_checkpoint_generation": "checkpoint:test", "release_id": "R2026a", "result_release_id": "R2026a", "owner": "matlab:shared", "result_owner": "matlab:shared", "result_sha256": "d" * 64, "accepted_result_sha256": "d" * 64},
    }}


def _summary_v47() -> dict[str, object]:
    ml_generation = "ml-v47"
    rl_generation = "rl-v47"
    agentic_generation = "agentic-v47"
    trial_generation = "trial-v47"
    partition = {"train": ["sample-1", "sample-2"], "validation": ["sample-3"]}
    return {
        "matlab_ml_rl_v47_identity": {
            "ml_datastore": {
                "generation": ml_generation,
                **{key: ml_generation for key in ("datastore_generation", "label_generation", "preprocess_generation", "partition_generation", "owner_generation", "result_generation")},
                "datastore_order": ["sample-1", "sample-2", "sample-3"], "result_datastore_order": ["sample-1", "sample-2", "sample-3"],
                "labels": ["cat", "dog", "cat"], "result_labels": ["cat", "dog", "cat"],
                "preprocess_sha256": "1" * 64, "result_preprocess_sha256": "1" * 64,
                "partition": partition, "result_partition": partition,
                "partition_owner": "partition:test", "result_partition_owner": "partition:test",
                "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
            },
            "reinforcement_learning": {
                "generation": rl_generation,
                **{key: rl_generation for key in ("observation_generation", "action_generation", "environment_generation", "reset_generation", "policy_generation", "result_generation")},
                "observation_spec_sha256": "3" * 64, "result_observation_spec_sha256": "3" * 64,
                "action_spec_sha256": "4" * 64, "result_action_spec_sha256": "4" * 64,
                "environment_id": "environment:test", "result_environment_id": "environment:test",
                "reset_policy": "deterministic_seeded", "result_reset_policy": "deterministic_seeded",
                "policy_owner": "policy:test", "result_policy_owner": "policy:test",
                "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
            },
            "agentic_toolkit": {
                "generation": agentic_generation,
                **{key: agentic_generation for key in ("call_generation", "correlation_generation", "session_generation", "workspace_generation", "release_generation", "result_generation")},
                "tool_call_id": "call:test", "result_tool_call_id": "call:test",
                "correlation_id": "correlation:test", "result_correlation_id": "correlation:test",
                "session_identity": "matlab:shared", "result_session_identity": "matlab:shared",
                "workspace_sha256": "6" * 64, "result_workspace_sha256": "6" * 64,
                "release_id": "R2026a", "result_release_id": "R2026a",
                "result_sha256": "7" * 64, "accepted_result_sha256": "7" * 64,
            },
            "experiment_trials": {
                "generation": trial_generation,
                **{key: trial_generation for key in ("experiment_generation", "trial_generation", "cv_generation", "result_index_generation", "result_generation")},
                "trial_row_keys": ["trial=0", "trial=1", "trial=2"], "result_trial_row_keys": ["trial=0", "trial=1", "trial=2"],
                "cv_partition_ids": ["fold=0", "fold=1", "fold=2"], "result_cv_partition_ids": ["fold=0", "fold=1", "fold=2"],
                "result_indices": [0, 1, 2], "replayed_result_indices": [0, 1, 2],
                "experiment_owner": "experiment:test", "result_experiment_owner": "experiment:test",
                "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
            },
        }
    }


def _summary_v48() -> dict[str, object]:
    autodiff = "autodiff-v48"
    sequence = "sequence-v48"
    agentic = "agentic-v48"
    experiment = "experiment-v48"
    return {"matlab_ml_rl_v48_identity": {
        "autodiff": {
            "generation": autodiff, **{key: autodiff for key in ("parameter_generation", "tape_generation", "objective_generation", "gradient_generation", "fd_generation", "checkpoint_generation", "result_generation")},
            "parameter_order": ["w1", "w2"], "result_parameter_order": ["w1", "w2"],
            "gradient_tape_sha256": "1" * 64, "result_gradient_tape_sha256": "1" * 64,
            "objective_id": "objective:validation-loss", "result_objective_id": "objective:validation-loss",
            "gradient": [1.5, -0.25], "result_gradient": [1.5, -0.25],
            "fd_spotcheck_indices": [0, 1], "result_fd_spotcheck_indices": [0, 1],
            "fd_gradient": [1.5001, -0.2499], "result_fd_gradient": [1.5001, -0.2499],
            "checkpoint_owner": "checkpoint:autodiff", "result_checkpoint_owner": "checkpoint:autodiff",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
        "sequence_model": {
            "generation": sequence, **{key: sequence for key in ("padding_generation", "mask_generation", "length_generation", "shuffle_generation", "minibatch_generation", "checkpoint_generation", "result_generation")},
            "padding_policy": "right_zero", "result_padding_policy": "right_zero",
            "padding_mask": [[1, 1, 1], [1, 1, 0]], "result_padding_mask": [[1, 1, 1], [1, 1, 0]],
            "sequence_lengths": [3, 2], "result_sequence_lengths": [3, 2],
            "shuffle_order": [1, 0], "result_shuffle_order": [1, 0],
            "minibatch_row_keys": ["sequence:1", "sequence:0"], "result_minibatch_row_keys": ["sequence:1", "sequence:0"],
            "checkpoint_owner": "checkpoint:sequence", "result_checkpoint_owner": "checkpoint:sequence",
            "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
        },
        "agentic_workspace": {
            "generation": agentic, **{key: agentic for key in ("mutation_generation", "diff_generation", "approval_generation", "rollback_generation", "tool_call_generation", "result_generation")},
            "workspace_mutations": ["models/train.m", "tests/test_train.m"], "result_workspace_mutations": ["models/train.m", "tests/test_train.m"],
            "workspace_diff_sha256": "4" * 64, "result_workspace_diff_sha256": "4" * 64,
            "approval_scope": "workspace_write:approved_paths", "result_approval_scope": "workspace_write:approved_paths",
            "rollback_state": "not_required", "result_rollback_state": "not_required",
            "tool_call_owner": "tool-call:test", "result_tool_call_owner": "tool-call:test",
            "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
        },
        "experiment_selection": {
            "generation": experiment, **{key: experiment for key in ("metric_generation", "direction_generation", "encoding_generation", "trial_generation", "selection_generation", "result_generation")},
            "metric_name": "validation_loss", "result_metric_name": "validation_loss",
            "metric_direction": "minimize", "result_metric_direction": "minimize",
            "categorical_encoding": {"optimizer": ["adam", "sgdm"]}, "result_categorical_encoding": {"optimizer": ["adam", "sgdm"]},
            "trial_row_keys": ["trial=0", "trial=1", "trial=2"], "result_trial_row_keys": ["trial=0", "trial=1", "trial=2"],
            "metric_values": [0.3, 0.2, 0.25], "result_metric_values": [0.3, 0.2, 0.25],
            "best_trial_index": 1, "result_best_trial_index": 1,
            "best_result_row": "trial=1", "result_best_result_row": "trial=1",
            "experiment_owner": "experiment:test", "result_experiment_owner": "experiment:test",
            "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64,
        },
    }}


def _summary_v49() -> dict[str, object]:
    rl = "rl-v49"; ml = "ml-v49"; agentic = "agentic-v49"; parallel = "parallel-v49"
    rows = ["transition:0", "transition:1", "transition:2"]
    observations = [[0.0, 1.0], [0.5, 0.8], [1.0, 0.2]]
    splits = {"train": ["sample:0", "sample:1", "sample:2", "sample:3"], "validation": ["sample:4", "sample:5"], "test": ["sample:6", "sample:7"]}
    workers = [1, 2, 3]; seeds = [7101, 7102, 7103]; streams = ["Threefry:1", "Threefry:2", "Threefry:3"]
    return {"matlab_ml_rl_v49_identity": {
        "rl_replay": {
            "generation": rl, **{key: rl for key in ("buffer_generation", "observation_generation", "action_generation", "reward_generation", "terminal_generation", "seed_generation", "checkpoint_generation", "policy_generation", "result_generation")},
            "replay_row_keys": rows, "result_replay_row_keys": rows, "observations": observations, "result_observations": observations,
            "actions": [0, 1, 0], "result_actions": [0, 1, 0], "rewards": [0.1, 1.0, -0.2], "result_rewards": [0.1, 1.0, -0.2],
            "terminals": [False, False, True], "result_terminals": [False, False, True], "episode_seeds": [4101, 4102, 4103], "result_episode_seeds": [4101, 4102, 4103],
            "checkpoint_owner": "checkpoint:rl", "result_checkpoint_owner": "checkpoint:rl", "policy_owner": "policy:rl", "result_policy_owner": "policy:rl",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        "ml_normalization": {
            "generation": ml, **{key: ml for key in ("normalization_generation", "class_generation", "split_generation", "fold_generation", "metric_generation", "model_generation", "result_generation")},
            "normalization_fit_scope": "training_partition_only", "result_normalization_fit_scope": "training_partition_only",
            "normalization_fit_rows": splits["train"], "result_normalization_fit_rows": splits["train"], "class_encoding": {"cold": 0, "hot": 1}, "result_class_encoding": {"cold": 0, "hot": 1},
            "split_row_keys": splits, "result_split_row_keys": splits, "training_fold_ids": [0, 1, 0, 1], "result_training_fold_ids": [0, 1, 0, 1],
            "metric_row_keys": ["fold:0", "fold:1"], "result_metric_row_keys": ["fold:0", "fold:1"], "metric_values": [0.8, 0.85], "result_metric_values": [0.8, 0.85],
            "model_owner": "model:classifier", "result_model_owner": "model:classifier", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
        "agentic_cancel": {
            "generation": agentic, **{key: agentic for key in ("timeout_generation", "cancel_generation", "partial_output_generation", "cleanup_generation", "tool_call_generation", "result_generation")},
            "timeout_s": 120.0, "result_timeout_s": 120.0, "timed_out": True, "result_timed_out": True,
            "cancel_requested": True, "result_cancel_requested": True, "cancel_completed": True, "result_cancel_completed": True,
            "partial_output_policy": "discard", "result_partial_output_policy": "discard", "session_cleanup": "released", "result_session_cleanup": "released",
            "tool_call_owner": "tool-call:cancel", "result_tool_call_owner": "tool-call:cancel", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
        },
        "parallel_resume": {
            "generation": parallel, **{key: parallel for key in ("worker_generation", "seed_generation", "stream_generation", "trial_generation", "resume_generation", "checkpoint_generation", "experiment_generation", "result_generation")},
            "worker_ids": workers, "result_worker_ids": workers, "worker_seeds": seeds, "result_worker_seeds": seeds, "random_streams": streams, "result_random_streams": streams,
            "trial_state": "completed:12", "result_trial_state": "completed:12", "resume_state": "resumed_from:8", "result_resume_state": "resumed_from:8",
            "checkpoint_owner": "checkpoint:parallel", "result_checkpoint_owner": "checkpoint:parallel", "experiment_owner": "experiment:parallel", "result_experiment_owner": "experiment:parallel",
            "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        },
    }}


def _summary_v50() -> dict[str, object]:
    bayes = "bayes-v50"
    mixed = "mixed-v50"
    sandbox = "sandbox-v50"
    tall = "tall-v50"
    constraints = [{"name": "temperature_c", "sense": "<=", "limit": 120.0}]
    acquisition = {
        "name": "expected-improvement-plus",
        "iteration": 18,
        "incumbent": "trial:12",
    }
    optimizer = {
        "name": "adam",
        "iteration": 240,
        "learn_rate": 0.001,
        "moment1_sha256": "1" * 64,
        "moment2_sha256": "2" * 64,
    }
    roots = ["workspace:project", "workspace:artifacts"]
    partitions = ["partition:0", "partition:1", "partition:2", "partition:3"]
    workers = [1, 2, 1, 2]
    return {"matlab_ml_rl_v50_identity": {
        "bayesopt": {
            "generation": bayes,
            **{key: bayes for key in ("objective_generation", "noise_generation", "constraint_generation", "seed_generation", "acquisition_generation", "model_generation", "result_generation")},
            "objective_name": "validation_loss", "result_objective_name": "validation_loss",
            "objective_noise_sigma": 0.02, "result_objective_noise_sigma": 0.02,
            "constraints": constraints, "result_constraints": constraints,
            "rng_seed": 50123, "result_rng_seed": 50123,
            "acquisition_state": acquisition, "result_acquisition_state": acquisition,
            "model_owner": "model:bayes-v50", "result_model_owner": "model:bayes-v50",
            "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
        },
        "mixed_precision": {
            "generation": mixed,
            **{key: mixed for key in ("precision_generation", "loss_scale_generation", "optimizer_generation", "checkpoint_generation", "network_generation", "result_generation")},
            "precision_policy": "mixed-fp16", "result_precision_policy": "mixed-fp16",
            "loss_scale": 1024.0, "result_loss_scale": 1024.0,
            "optimizer_state": optimizer, "result_optimizer_state": optimizer,
            "checkpoint_sha256": "4" * 64, "result_checkpoint_sha256": "4" * 64,
            "network_owner": "network:mixed-v50", "result_network_owner": "network:mixed-v50",
            "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
        },
        "agentic_sandbox": {
            "generation": sandbox,
            **{key: sandbox for key in ("root_generation", "resolution_generation", "traversal_generation", "approval_generation", "tool_call_generation", "result_generation")},
            "allowed_roots": roots, "result_allowed_roots": roots,
            "requested_path": "workspace:project/src/train.m", "result_requested_path": "workspace:project/src/train.m",
            "resolved_path": "workspace:project/src/train.m", "result_resolved_path": "workspace:project/src/train.m",
            "symlink_resolution": "inside_allowed_root", "result_symlink_resolution": "inside_allowed_root",
            "traversal_detected": False, "result_traversal_detected": False,
            "approval_id": "approval:v50", "result_approval_id": "approval:v50",
            "tool_call_owner": "tool-call:v50", "result_tool_call_owner": "tool-call:v50",
            "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64,
        },
        "tall_datastore": {
            "generation": tall,
            **{key: tall for key in ("partition_generation", "worker_generation", "checkpoint_generation", "resume_generation", "datastore_generation", "result_generation")},
            "partition_ids": partitions, "result_partition_ids": partitions,
            "worker_order": workers, "result_worker_order": workers,
            "checkpoint_sha256": "7" * 64, "result_checkpoint_sha256": "7" * 64,
            "resume_cursor": "partition:2/row:4096", "result_resume_cursor": "partition:2/row:4096",
            "datastore_owner": "datastore:v50", "result_datastore_owner": "datastore:v50",
            "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
        },
    }}


def _summary_v51() -> dict[str, object]:
    generation = "matlab-public-v51"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    transitions = ["transition:1", "transition:2", "transition:3", "transition:4"]
    classes = ["normal", "warning", "fault"]
    identity = {
        "rl_replay": {
            "generation": generation, **{name: generation for name in ("transition_generation", "terminal_generation", "discount_generation", "priority_generation", "target_generation", "owner_generation", "result_generation")},
            "transition_ids": transitions, "result_transition_ids": transitions, "terminal_flags": [False, False, False, True],
            "result_terminal_flags": [False, False, False, True], "truncation_flags": [False] * 4,
            "result_truncation_flags": [False] * 4, "discounts": [0.99, 0.99, 0.99, 0.0],
            "result_discounts": [0.99, 0.99, 0.99, 0.0], "priorities": [0.2, 0.5, 0.4, 1.0],
            "result_priorities": [0.2, 0.5, 0.4, 1.0], "target_network_sha256": "1" * 64,
            "result_target_network_sha256": "1" * 64, "buffer_owner": "buffer:prioritized-v51",
            "result_buffer_owner": "buffer:prioritized-v51", **result,
        },
        "classification": {
            "generation": generation, **{name: generation for name in ("class_generation", "label_generation", "split_generation", "confusion_generation", "owner_generation", "result_generation")},
            "class_order": classes, "result_class_order": classes, "label_encoding": {"normal": 0, "warning": 1, "fault": 2},
            "result_label_encoding": {"normal": 0, "warning": 1, "fault": 2},
            "data_split_sha256": {"train": "2" * 64, "validation": "3" * 64, "test": "4" * 64},
            "result_data_split_sha256": {"train": "2" * 64, "validation": "3" * 64, "test": "4" * 64},
            "confusion_matrix": [[18, 1, 0], [2, 15, 1], [0, 1, 12]],
            "result_confusion_matrix": [[18, 1, 0], [2, 15, 1], [0, 1, 12]],
            "confusion_matrix_axes": {"rows": classes, "columns": classes},
            "result_confusion_matrix_axes": {"rows": classes, "columns": classes},
            "model_owner": "model:classifier-v51", "result_model_owner": "model:classifier-v51", **result,
        },
        "agentic_edit": {
            "generation": generation, **{name: generation for name in ("precondition_generation", "patch_generation", "rollback_generation", "tool_call_generation", "owner_generation", "result_generation")},
            "target_path": "workspace:project/src/agent.m", "result_target_path": "workspace:project/src/agent.m",
            "file_precondition_sha256": "5" * 64, "result_file_precondition_sha256": "5" * 64,
            "patch_sha256": "6" * 64, "result_patch_sha256": "6" * 64, "rollback_state": "available",
            "result_rollback_state": "available", "rollback_sha256": "5" * 64, "result_rollback_sha256": "5" * 64,
            "tool_call_owner": "tool-call:edit-v51", "result_tool_call_owner": "tool-call:edit-v51", **result,
        },
        "parallel_rng": {
            "generation": generation, **{name: generation for name in ("rng_generation", "substream_generation", "worker_generation", "reduction_generation", "checkpoint_generation", "owner_generation", "result_generation")},
            "rng_algorithm": "Threefry", "result_rng_algorithm": "Threefry", "rng_seed": 510901, "result_rng_seed": 510901,
            "substream_ids": [1, 2, 3, 4], "result_substream_ids": [1, 2, 3, 4], "worker_map": [1, 2, 1, 2],
            "result_worker_map": [1, 2, 1, 2], "reduction_order": [1, 3, 2, 4], "result_reduction_order": [1, 3, 2, 4],
            "checkpoint_sha256": "7" * 64, "result_checkpoint_sha256": "7" * 64,
            "pool_owner": "pool:threads-v51", "result_pool_owner": "pool:threads-v51", **result,
        },
    }
    return {"matlab_ml_rl_v51_identity": identity}


def _summary_v52() -> dict[str, object]:
    generation = "matlab-v52-test"
    generations = lambda names: {name: generation for name in names}
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    hidden = ["1" * 64, "2" * 64, "3" * 64, "4" * 64]
    hierarchy = ["plant/velocity", "controller/command", "sensor/current"]
    identity = {
        "ml_training": {
            "generation": generation,
            **generations(("gradient_generation", "optimizer_generation", "schedule_generation", "checkpoint_generation", "owner_generation", "result_generation")),
            "gradient_clip_threshold": 1.0, "result_gradient_clip_threshold": 1.0,
            "gradient_clip_method": "global-l2norm", "result_gradient_clip_method": "global-l2norm",
            "optimizer_name": "adam", "result_optimizer_name": "adam",
            "optimizer_step": 240, "result_optimizer_step": 240,
            "optimizer_state_sha256": "4" * 64, "result_optimizer_state_sha256": "4" * 64,
            "learning_rate_schedule": [0.001, 0.0005, 0.0001], "result_learning_rate_schedule": [0.001, 0.0005, 0.0001],
            "checkpoint_sha256": "5" * 64, "result_checkpoint_sha256": "5" * 64,
            "model_owner": "model:training-v52", "result_model_owner": "model:training-v52", **result,
        },
        "recurrent_rl": {
            "generation": generation,
            **generations(("hidden_generation", "reset_generation", "mask_generation", "policy_generation", "owner_generation", "result_generation")),
            "hidden_state_sha256": hidden, "result_hidden_state_sha256": hidden,
            "episode_reset": [True, False, False, False], "result_episode_reset": [True, False, False, False],
            "sequence_mask": [True, True, True, False], "result_sequence_mask": [True, True, True, False],
            "policy_sha256": "6" * 64, "result_policy_sha256": "6" * 64,
            "policy_owner": "policy:recurrent-v52", "result_policy_owner": "policy:recurrent-v52", **result,
        },
        "codegen": {
            "generation": generation,
            **generations(("target_generation", "numeric_generation", "config_generation", "build_generation", "owner_generation", "result_generation")),
            "target_hardware": "ARM-Cortex-M7", "result_target_hardware": "ARM-Cortex-M7",
            "numeric_type": "single", "result_numeric_type": "single",
            "configuration_sha256": "7" * 64, "result_configuration_sha256": "7" * 64,
            "build_sha256": "8" * 64, "result_build_sha256": "8" * 64,
            "build_owner": "build:embedded-v52", "result_build_owner": "build:embedded-v52", **result,
        },
        "simulink_data_inspector": {
            "generation": generation,
            **generations(("run_generation", "signal_generation", "unit_generation", "interpolation_generation", "owner_generation", "result_generation")),
            "run_id": "run:sdi-v52", "result_run_id": "run:sdi-v52",
            "signal_hierarchy": hierarchy, "result_signal_hierarchy": hierarchy,
            "signal_units": ["m/s", "N", "A"], "result_signal_units": ["m/s", "N", "A"],
            "interpolation": ["linear", "zoh", "zoh"], "result_interpolation": ["linear", "zoh", "zoh"],
            "session_owner": "session:sdi-v52", "result_session_owner": "session:sdi-v52", **result,
        },
    }
    return {"matlab_ml_rl_v52_identity": identity}


def _softmax(row: list[float], temperature: float) -> list[float]:
    scaled = [item / temperature for item in row]
    maximum = max(scaled)
    values = [math.exp(item - maximum) for item in scaled]
    total = sum(values)
    return [item / total for item in values]


def _summary_v53() -> dict[str, object]:
    offline_generation = "offline-rl-v53-test"
    calibration_generation = "calibration-v53-test"
    behavior = [0.40, 0.25, 0.20, 0.15]
    target = [0.35, 0.30, 0.20, 0.15]
    weights = [target_item / behavior_item for target_item, behavior_item in zip(target, behavior)]
    logits = [[2.0, 0.5, -1.0], [0.1, 1.2, -0.2]]
    temperature = 1.4
    probabilities = [_softmax(row, temperature) for row in logits]
    identity = {
        "offline_rl": {
            "generation": offline_generation,
            **{name: offline_generation for name in ("behavior_generation", "weight_generation", "support_generation", "dataset_generation", "owner_generation", "result_generation")},
            "behavior_policy_sha256": "1" * 64, "result_behavior_policy_sha256": "1" * 64,
            "target_policy_sha256": "2" * 64, "result_target_policy_sha256": "2" * 64,
            "behavior_action_probability": behavior, "result_behavior_action_probability": behavior,
            "target_action_probability": target, "result_target_action_probability": target,
            "importance_weight": weights, "result_importance_weight": weights,
            "support_coverage_fraction": 1.0, "result_support_coverage_fraction": 1.0,
            "dataset_owner": "dataset:offline-v53", "result_dataset_owner": "dataset:offline-v53",
            "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
        },
        "probability_calibration": {
            "generation": calibration_generation,
            **{name: calibration_generation for name in ("logit_generation", "probability_generation", "class_generation", "temperature_generation", "owner_generation", "result_generation")},
            "logits": logits, "result_logits": logits,
            "calibrated_probability": probabilities, "result_calibrated_probability": probabilities,
            "class_order": ["class:normal", "class:warning", "class:fault"],
            "result_class_order": ["class:normal", "class:warning", "class:fault"],
            "temperature": temperature, "result_temperature": temperature,
            "model_owner": "model:calibration-v53", "result_model_owner": "model:calibration-v53",
            "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        },
    }
    return {"matlab_ml_rl_v53_identity": identity}


def _weights(priorities: list[float], alpha: float, beta: float) -> tuple[list[float], list[float]]:
    scaled = [value**alpha for value in priorities]
    probabilities = [value / sum(scaled) for value in scaled]
    raw = [(len(priorities) * probability) ** (-beta) for probability in probabilities]
    return probabilities, [value / max(raw) for value in raw]


def _summary_v54() -> dict[str, object]:
    replay_generation = "replay-v54-test"
    preprocess_generation = "preprocess-v54-test"
    priorities = [1.0, 2.0, 4.0, 8.0]
    probabilities, weights = _weights(priorities, 0.6, 0.4)
    categories = {"material": ["category:steel", "category:copper"], "state": ["category:normal", "category:fault"]}
    split = {"train": [0, 1, 2, 3, 4, 5], "validation": [6, 7], "test": [8, 9]}
    identity = {
        "prioritized_replay": {
            "generation": replay_generation,
            **{name: replay_generation for name in ("priority_generation", "beta_generation", "seed_generation", "policy_generation", "owner_generation", "result_generation")},
            "priorities": priorities, "result_priorities": priorities,
            "priority_alpha": 0.6, "result_priority_alpha": 0.6,
            "sampling_probability": probabilities, "result_sampling_probability": probabilities,
            "importance_beta": 0.4, "result_importance_beta": 0.4,
            "importance_weight": weights, "result_importance_weight": weights,
            "rng_seed": 8675309, "result_rng_seed": 8675309,
            "policy_checkpoint_sha256": "1" * 64, "result_policy_checkpoint_sha256": "1" * 64,
            "buffer_owner": "buffer:v54", "result_buffer_owner": "buffer:v54",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
        "feature_preprocess": {
            "generation": preprocess_generation,
            **{name: preprocess_generation for name in ("normalization_generation", "category_generation", "split_generation", "model_generation", "owner_generation", "result_generation")},
            "feature_mean": [10.0, 2.0], "result_feature_mean": [10.0, 2.0],
            "feature_std": [2.0, 0.5], "result_feature_std": [2.0, 0.5],
            "category_order": categories, "result_category_order": categories,
            "data_split_indices": split, "result_data_split_indices": split,
            "model_checkpoint_sha256": "3" * 64, "result_model_checkpoint_sha256": "3" * 64,
            "preprocessing_owner": "preprocess:v54", "result_preprocessing_owner": "preprocess:v54",
            "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        },
    }
    return {"matlab_ml_rl_v54_identity": identity}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _summary_v55() -> dict[str, object]:
    gamma = 0.9
    rewards = [1.0, 0.5, -0.2]
    nstep_return = sum(gamma**index * reward for index, reward in enumerate(rewards)) + gamma**len(rewards) * 2.0
    folds = [0, 0, 1, 1, 2, 2]
    labels = ["class:A", "class:B"] * 3
    fit_rows = {"0": [2, 3, 4, 5], "1": [0, 1, 4, 5], "2": [0, 1, 2, 3]}
    metrics = [0.8, 0.7, 0.9]
    return {
        "matlab_ml_rl_v54_identity": {
            "nstep_return": {
                **_generations("nstep-v55-test", ("reward_generation", "gamma_generation", "terminal_generation", "bootstrap_generation", "return_generation", "trajectory_generation", "owner_generation", "result_generation")),
                "gamma": gamma, "result_gamma": gamma,
                "rewards": rewards, "result_rewards": rewards,
                "terminal": False, "result_terminal": False,
                "bootstrap_value": 2.0, "result_bootstrap_value": 2.0,
                "n_step_return": nstep_return, "result_n_step_return": nstep_return,
                "trajectory_id": "trajectory:v55", "result_trajectory_id": "trajectory:v55",
                "policy_owner": "policy:nstep-v55", "result_policy_owner": "policy:nstep-v55",
                "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
            },
            "cross_validation": {
                **_generations("crossval-v55-test", ("fold_generation", "stratification_generation", "preprocess_generation", "seed_generation", "metric_generation", "owner_generation", "result_generation")),
                "fold_id_per_sample": folds, "result_fold_id_per_sample": folds,
                "class_labels": labels, "result_class_labels": labels,
                "preprocess_fit_rows": fit_rows, "result_preprocess_fit_rows": fit_rows,
                "rng_seed": 1729, "result_rng_seed": 1729,
                "fold_metrics": metrics, "result_fold_metrics": metrics,
                "aggregate_metric": sum(metrics) / len(metrics), "result_aggregate_metric": sum(metrics) / len(metrics),
                "model_owner": "model:crossval-v55", "result_model_owner": "model:crossval-v55",
                "result_sha256": "9" * 64, "accepted_result_sha256": "9" * 64,
            },
        }
    }
