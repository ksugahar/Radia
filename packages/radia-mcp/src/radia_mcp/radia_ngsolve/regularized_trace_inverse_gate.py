"""Solver-neutral gate for a regularized FEM/BEM trace inverse path."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from .slot_gates import lcurve_corner_choice, morozov_discrepancy_choice


def regularized_trace_inverse_path_gate(
    summary: Mapping[str, Any],
    *,
    max_solution_relative_error: float = 1.0e-10,
    max_trace_relative_error: float = 1.0e-10,
    max_regularized_objective_relative_error: float = 1.0e-10,
    max_zero_alpha_objective_absolute_error: float = 1.0e-20,
    max_normal_equation_residual: float = 1.0e-9,
    max_gradient_check_absolute_error: float = 2.0e-7,
    max_replay_relative_error: float = 1.0e-12,
    monotonicity_tolerance: float = 1.0e-12,
) -> dict[str, Any]:
    """Validate P1 trace regularization, parameter choices, and replay.

    The input is deliberately solver-neutral.  It records a first-order
    tetrahedron/triangle trace path, two regularization choices, independent
    linear-solver checks, and deterministic replay.  The gate recomputes both
    L-curve and Morozov choices instead of trusting reported indices.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    tolerances = {
        "max_solution_relative_error": _nonnegative(
            max_solution_relative_error, "max_solution_relative_error"
        ),
        "max_trace_relative_error": _nonnegative(
            max_trace_relative_error, "max_trace_relative_error"
        ),
        "max_regularized_objective_relative_error": _nonnegative(
            max_regularized_objective_relative_error,
            "max_regularized_objective_relative_error",
        ),
        "max_zero_alpha_objective_absolute_error": _nonnegative(
            max_zero_alpha_objective_absolute_error,
            "max_zero_alpha_objective_absolute_error",
        ),
        "max_normal_equation_residual": _nonnegative(
            max_normal_equation_residual, "max_normal_equation_residual"
        ),
        "max_gradient_check_absolute_error": _nonnegative(
            max_gradient_check_absolute_error, "max_gradient_check_absolute_error"
        ),
        "max_replay_relative_error": _nonnegative(
            max_replay_relative_error, "max_replay_relative_error"
        ),
        "monotonicity_tolerance": _nonnegative(
            monotonicity_tolerance, "monotonicity_tolerance"
        ),
    }

    mesh = _mapping(summary, "mesh")
    path = _mapping(summary, "path")
    problem = _mapping(summary, "problem")
    reported_lcurve = _mapping(summary, "lcurve")
    reported_morozov = _mapping(summary, "morozov")
    crosscheck = _mapping(summary, "crosscheck")
    replay = _mapping(summary, "replay")

    alphas = _float_list(path, "alphas")
    solution_norms = _float_list(path, "solution_norms")
    trace_residuals = _float_list(path, "trace_residual_norms")
    weighted_residuals = _float_list(path, "weighted_trace_residuals")
    normal_residuals = _float_list(path, "normal_equation_residuals")
    gradient_errors = _float_list(path, "gradient_check_max_abs_errors")
    lengths = {
        len(alphas),
        len(solution_norms),
        len(trace_residuals),
        len(weighted_residuals),
        len(normal_residuals),
        len(gradient_errors),
    }
    if len(lengths) != 1 or not alphas:
        raise ValueError("all regularization path arrays must have one nonzero common length")
    if len(alphas) < 5:
        raise ValueError("at least five regularization path rows are required")
    if any(value < 0.0 for value in alphas + solution_norms + trace_residuals + weighted_residuals):
        raise ValueError("path weights and norms must be non-negative")

    lcurve = lcurve_corner_choice(
        alphas[1:],
        trace_residuals[1:],
        solution_norms[1:],
        tol=tolerances["monotonicity_tolerance"],
    )
    morozov = morozov_discrepancy_choice(
        alphas,
        weighted_residuals,
        _finite_float(problem, "noise_norm"),
        tol=tolerances["monotonicity_tolerance"],
    )

    alpha_path_ok = alphas[0] == 0.0 and all(
        right > left for left, right in zip(alphas[1:], alphas[2:])
    )
    solution_monotone = all(
        right <= left + tolerances["monotonicity_tolerance"]
        for left, right in zip(solution_norms, solution_norms[1:])
    )
    trace_monotone = all(
        right + tolerances["monotonicity_tolerance"] >= left
        for left, right in zip(trace_residuals, trace_residuals[1:])
    )
    weighted_monotone = all(
        right + tolerances["monotonicity_tolerance"] >= left
        for left, right in zip(weighted_residuals, weighted_residuals[1:])
    )
    row_identity_ok = _optional_row_identity_is_aligned(path, len(alphas))
    gradient_resolution_ok = _optional_gradient_check_is_resolved(
        path, len(alphas)
    )
    gradient_generation_ok, solution_run_generation_ok = (
        _optional_generation_ids_are_aligned(path, len(alphas))
    )
    quadrature_generation_ok = _optional_quadrature_generation_ids_are_aligned(
        path, len(alphas)
    )
    curvature_parameterization_ok = (
        _optional_curvature_parameterization_is_aligned(reported_lcurve)
    )
    design_variable_identity_ok = _optional_design_variable_identity_is_aligned(
        summary
    )
    convolution_quadrature_identity_ok = (
        _optional_convolution_quadrature_identity_is_aligned(summary)
    )
    cq_contour_branch_identity_ok = (
        _optional_cq_contour_branch_identity_is_aligned(summary)
    )
    trace_orientation_identity_ok = (
        _optional_trace_orientation_identity_is_aligned(summary)
    )
    cq_starting_weights_identity_ok = (
        _optional_cq_starting_weights_identity_is_aligned(summary)
    )
    hmatrix_permutation_identity_ok = (
        _optional_hmatrix_cluster_permutation_identity_is_aligned(summary)
    )
    trace_surface_normal_identity_ok = (
        _optional_fembem_trace_surface_normal_orientation_is_aligned(summary)
    )
    cq_fft_normalization_identity_ok = (
        _optional_cq_inverse_z_transform_fft_normalization_is_aligned(summary)
    )
    sesquilinear_identity_ok = (
        _optional_fembem_sesquilinear_inner_product_is_aligned(summary)
    )
    hmatrix_tolerance_norm_identity_ok = (
        _optional_hmatrix_low_rank_tolerance_norm_basis_is_aligned(summary)
    )
    complex_symmetry_residual_identity_ok = (
        _optional_complex_operator_symmetry_residual_norm_is_aligned(summary)
    )
    hmatrix_diameter_metric_identity_ok = (
        _optional_hmatrix_admissibility_diameter_metric_is_aligned(summary)
    )
    cq_weight_branch_order_identity_ok = (
        _optional_cq_convolution_weight_branch_order_is_aligned(summary)
    )
    trace_boundary_node_generation_identity_ok = (
        _optional_fembem_trace_boundary_node_generation_is_aligned(summary)
    )
    cq_frequency_contour_identity_ok = (
        _optional_cq_frequency_contour_radius_damping_is_aligned(summary)
    )
    hmatrix_aca_pivot_permutation_identity_ok = (
        _optional_hmatrix_aca_pivot_cluster_permutation_is_aligned(summary)
    )
    cq_conjugate_bin_ownership_identity_ok = (
        _optional_cq_conjugate_frequency_bin_ownership_is_aligned(summary)
    )
    hmatrix_bbox_mesh_scale_identity_ok = (
        _optional_hmatrix_bounding_box_mesh_scale_is_aligned(summary)
    )
    cq_contour_fft_timestep_identity_ok = (
        _optional_cq_contour_fft_timestep_generation_is_aligned(summary)
    )
    hmatrix_aca_tolerance_rank_identity_ok = (
        _optional_hmatrix_aca_tolerance_norm_rank_is_aligned(summary)
    )
    fembem_trace_map_identity_ok = (
        _optional_fembem_trace_orientation_normal_node_map_is_aligned(summary)
    )
    cq_causality_pair_identity_ok = (
        _optional_cq_causality_conjugate_contour_pair_is_aligned(summary)
    )
    p1_boundary_row_identity_ok = (
        _optional_p1_boundary_mass_trace_row_node_identity_is_aligned(summary)
    )
    cq_restart_identity_ok = (
        _optional_cq_restart_history_weight_segment_identity_is_aligned(summary)
    )
    simscape_file_solid_identity_ok = (
        _optional_simscape_file_solid_identity_is_aligned(summary)
    )
    multibody_xml_identity_ok = (
        _optional_multibody_xml_import_identity_is_aligned(summary)
    )
    parallel_pool_identity_ok = (
        _optional_parallel_pool_identity_is_aligned(summary)
    )
    autodiff_tape_identity_ok = (
        _optional_autodiff_tape_identity_is_aligned(summary)
    )
    fembem_interface_identity_ok = (
        _optional_fembem_interface_identity_is_aligned(summary)
    )
    cq_time_history_identity_ok = (
        _optional_cq_time_history_identity_is_aligned(summary)
    )
    hmatrix_block_tree_identity_ok = (
        _optional_hmatrix_block_tree_identity_is_aligned(summary)
    )
    ad_gradient_identity_ok = _optional_ad_gradient_identity_is_aligned(summary)
    cq_transfer_identity_ok = _optional_cq_transfer_identity_is_aligned(summary)
    fembem_coupling_identity_ok = (
        _optional_fembem_coupling_identity_is_aligned(summary)
    )
    adaptive_cq_identity_ok = _optional_adaptive_cq_identity_is_aligned(summary)
    p1_fembem_discretization_identity_ok = (
        _optional_p1_fembem_discretization_identity_is_aligned(summary)
    )
    hmatrix_aca_cluster_identity_ok = (
        _optional_hmatrix_aca_cluster_identity_is_aligned(summary)
    )
    calderon_cq_identity_ok = _optional_calderon_cq_identity_is_aligned(summary)
    near_singular_quadrature_identity_ok = (
        _optional_bem_near_singular_quadrature_identity_is_aligned(summary)
    )
    fembem_energy_reciprocity_identity_ok = (
        _optional_fembem_energy_reciprocity_identity_is_aligned(summary)
    )
    hmatrix_recompression_identity_ok = (
        _optional_hmatrix_recompression_identity_is_aligned(summary)
    )
    cq_block_restart_identity_ok = (
        _optional_cq_block_restart_identity_is_aligned(summary)
    )
    complex_ad_identity_ok = _optional_complex_ad_identity_is_aligned(summary)
    pde_quadratic_vol_identity_ok = (
        _optional_pde_quadratic_vol_identity_is_aligned(summary)
    )
    adaptive_cq_restart_identity_ok = (
        _optional_adaptive_cq_restart_identity_is_aligned(summary)
    )
    modal_fembem_transient_identity_ok = (
        _optional_modal_fembem_transient_identity_is_aligned(summary)
    )
    calderon_projector_identity_ok = (
        _optional_calderon_projector_identity_is_aligned(summary)
    )
    cq_physical_closure_identity_ok = (
        _optional_cq_physical_closure_identity_is_aligned(summary)
    )
    fembem_reciprocity_identity_ok = (
        _optional_fembem_reciprocity_identity_is_aligned(summary)
    )
    nonlinear_eigen_contour_identity_ok = (
        _optional_nonlinear_eigen_contour_identity_is_aligned(summary)
    )
    cq_acoustic_identity_ok = _optional_cq_acoustic_identity_is_aligned(summary)
    fembem_autodiff_identity_ok = (
        _optional_fembem_autodiff_identity_is_aligned(summary)
    )
    adaptive_cq_timestep_identity_ok = (
        _optional_cq_adaptive_timestep_identity_is_aligned(summary)
    )
    fembem_shape_derivative_identity_ok = (
        _optional_fembem_shape_derivative_identity_is_aligned(summary)
    )
    hmatrix_benchmark_identity_ok = (
        _optional_hmatrix_benchmark_identity_is_aligned(summary)
    )
    multifrequency_adjoint_identity_ok = (
        _optional_multifrequency_adjoint_identity_is_aligned(summary)
    )
    simp_topology_identity_ok = _optional_simp_topology_identity_is_aligned(summary)
    fembem_model_reduction_identity_ok = (
        _optional_fembem_model_reduction_identity_is_aligned(summary)
    )
    nonlinear_fem_newton_identity_ok = (
        _optional_nonlinear_fem_newton_identity_is_aligned(summary)
    )
    cq_contour_reconstruction_identity_ok = (
        _optional_cq_contour_reconstruction_identity_is_aligned(summary)
    )
    johnson_nedelec_identity_ok = (
        _optional_johnson_nedelec_identity_is_aligned(summary)
    )
    adjoint_hessian_identity_ok = (
        _optional_adjoint_hessian_identity_is_aligned(summary)
    )
    cq_acoustic_history_identity_ok = (
        _optional_cq_acoustic_history_identity_is_aligned(summary)
    )
    hmatrix_compression_identity_ok = (
        _optional_hmatrix_compression_identity_is_aligned(summary)
    )

    checks = {
        "schema_is_regularized_trace_inverse_v1": (
            str(summary.get("schema", "")) == "regularized_trace_inverse_path/v1"
        ),
        "mesh_is_first_order_tri_tet_trace": (
            str(mesh.get("volume_element", "")) == "tetrahedron"
            and str(mesh.get("boundary_element", "")) == "triangle"
            and _integer(mesh, "polynomial_order") == 1
            and _positive_integer(mesh, "tetrahedra") > 0
            and _positive_integer(mesh, "triangles") > 0
            and _positive_integer(mesh, "surface_nodes")
            < _positive_integer(mesh, "volume_nodes")
            and _positive_integer(mesh, "trace_rows")
            == _positive_integer(mesh, "surface_nodes")
            and _positive_integer(mesh, "fem_unknowns")
            == _positive_integer(mesh, "volume_nodes")
            and _positive_integer(mesh, "trace_nnz")
            == _positive_integer(mesh, "surface_nodes")
        ),
        "zero_then_strictly_increasing_alpha_path": alpha_path_ok,
        "solution_norm_decreases_along_path": solution_monotone,
        "trace_residual_increases_along_path": trace_monotone,
        "weighted_residual_increases_along_path": weighted_monotone,
        "normal_equations_close": (
            max(normal_residuals) <= tolerances["max_normal_equation_residual"]
        ),
        "finite_difference_gradients_close": (
            max(gradient_errors) <= tolerances["max_gradient_check_absolute_error"]
        ),
        "alpha_path_row_identity_is_aligned": row_identity_ok,
        "finite_difference_steps_are_numerically_resolved": gradient_resolution_ok,
        "gradient_uses_current_parameter_generation": gradient_generation_ok,
        "solution_rows_share_regularization_run_generation": (
            solution_run_generation_ok
        ),
        "gradient_uses_current_boundary_quadrature_generation": (
            quadrature_generation_ok
        ),
        "lcurve_curvature_uses_recorded_path_parameterization": (
            curvature_parameterization_ok
        ),
        "adjoint_and_finite_difference_share_design_variable_order": (
            design_variable_identity_ok
        ),
        "cq_weights_match_current_time_grid_and_method": (
            convolution_quadrature_identity_ok
        ),
        "cq_transfer_samples_share_inverse_laplace_contour_branch": (
            cq_contour_branch_identity_ok
        ),
        "fembem_trace_orientation_matches_current_volume_mesh": (
            trace_orientation_identity_ok
        ),
        "cq_starting_weights_match_multistep_startup_scheme": (
            cq_starting_weights_identity_ok
        ),
        "hmatrix_cluster_permutation_matches_boundary_triangle_order": (
            hmatrix_permutation_identity_ok
        ),
        "fembem_trace_and_surface_share_outward_normal_orientation": (
            trace_surface_normal_identity_ok
        ),
        "cq_inverse_z_transform_uses_matching_fft_normalization": (
            cq_fft_normalization_identity_ok
        ),
        "fembem_operators_share_sesquilinear_conjugation_convention": (
            sesquilinear_identity_ok
        ),
        "hmatrix_low_rank_acceptance_uses_calibrated_norm_basis": (
            hmatrix_tolerance_norm_identity_ok
        ),
        "complex_operator_residual_uses_transpose_symmetry_class": (
            complex_symmetry_residual_identity_ok
        ),
        "hmatrix_admissibility_uses_one_cluster_diameter_metric": (
            hmatrix_diameter_metric_identity_ok
        ),
        "cq_convolution_weights_share_ztransform_branch_and_multistep_order": (
            cq_weight_branch_order_identity_ok
        ),
        "fembem_trace_matrix_uses_current_boundary_node_generation": (
            trace_boundary_node_generation_identity_ok
        ),
        "cq_inverse_transform_uses_current_frequency_contour_metadata": (
            cq_frequency_contour_identity_ok
        ),
        "hmatrix_aca_pivots_use_current_cluster_permutation": (
            hmatrix_aca_pivot_permutation_identity_ok
        ),
        "cq_inverse_transform_uses_current_conjugate_frequency_bins": (
            cq_conjugate_bin_ownership_identity_ok
        ),
        "hmatrix_block_clusters_use_current_mesh_scale_bounding_boxes": (
            hmatrix_bbox_mesh_scale_identity_ok
        ),
        "cq_inverse_transform_uses_current_contour_fft_scaling_and_timestep": (
            cq_contour_fft_timestep_identity_ok
        ),
        "hmatrix_aca_uses_current_block_tolerance_norm_and_rank": (
            hmatrix_aca_tolerance_rank_identity_ok
        ),
        "fembem_trace_uses_current_orientation_normals_and_volume_node_map": (
            fembem_trace_map_identity_ok
        ),
        "cq_inverse_transform_uses_current_causal_conjugate_contour_pairs": (
            cq_causality_pair_identity_ok
        ),
        "p1_boundary_mass_and_trace_rows_match_current_node_and_mesh_generation": (
            p1_boundary_row_identity_ok
        ),
        "cq_restart_reuses_current_weight_history_segment_and_time_grid": (
            cq_restart_identity_ok
        ),
        "simscape_file_solid_uses_current_geometry_density_inertia_and_frame": (
            simscape_file_solid_identity_ok
        ),
        "multibody_xml_uses_current_joint_axes_transforms_units_and_geometry": (
            multibody_xml_identity_ok
        ),
        "parallel_results_use_current_worker_paths_devices_rng_and_code": (
            parallel_pool_identity_ok
        ),
        "autodiff_gradients_use_current_tape_variables_mesh_objective_and_primal": (
            autodiff_tape_identity_ok
        ),
        "fembem_trace_uses_current_normals_nodes_units_mesh_and_operator": (
            fembem_interface_identity_ok
        ),
        "cq_time_history_uses_current_contour_weights_startup_and_causality_window": (
            cq_time_history_identity_ok
        ),
        "hmatrix_uses_current_block_tree_admissibility_permutations_tolerance_kernel_and_mesh": (
            hmatrix_block_tree_identity_ok
        ),
        "ad_gradient_uses_current_tape_material_operator_mesh_objective_and_primal": (
            ad_gradient_identity_ok
        ),
        "cq_history_uses_current_contour_timestep_branch_transfer_and_inverse_transform": (
            cq_transfer_identity_ok
        ),
        "fembem_coupling_uses_current_trace_normals_material_wavenumber_matrices_and_mesh": (
            fembem_coupling_identity_ok
        ),
        "adaptive_cq_uses_current_contour_order_startup_error_estimator_restart_and_history": (
            adaptive_cq_identity_ok
        ),
        "p1_fembem_uses_current_boundary_orientation_quadrature_singular_trace_matrices_and_mesh": (
            p1_fembem_discretization_identity_ok
        ),
        "hmatrix_aca_uses_current_clusters_permutation_admissibility_rank_tolerance_kernel_mesh_and_result": (
            hmatrix_aca_cluster_identity_ok
        ),
        "calderon_cq_uses_current_v_k_trace_normals_frequency_inverse_mesh_and_result": (
            calderon_cq_identity_ok
        ),
        "bem_near_singular_quadrature_uses_current_distance_size_order_map_kernel_reference_mesh_and_result": (
            near_singular_quadrature_identity_ok
        ),
        "fembem_energy_flux_and_reciprocity_use_current_trace_normal_frequency_incident_field_and_result": (
            fembem_energy_reciprocity_identity_ok
        ),
        "hmatrix_recompression_uses_current_svd_tolerance_norm_ranks_permutations_operator_mesh_and_result": (
            hmatrix_recompression_identity_ok
        ),
        "cq_block_restart_uses_current_blocks_history_startup_weights_time_samples_owners_and_result": (
            cq_block_restart_identity_ok
        ),
        "complex_ad_uses_current_wirtinger_conjugation_branch_scaling_fd_mesh_and_result": (
            complex_ad_identity_ok
        ),
        "pde_quadratic_vol_uses_current_midnodes_tets_boundary_orientation_regions_order_and_mesh": (
            pde_quadratic_vol_identity_ok
        ),
        "adaptive_cq_uses_current_timesteps_contour_rebuild_history_interpolation_error_restart_operator_mesh_and_result": (
            adaptive_cq_restart_identity_ok
        ),
        "modal_fembem_transient_uses_current_mass_damping_initial_projection_truncation_energy_mesh_history_and_result": (
            modal_fembem_transient_identity_ok
        ),
        "calderon_projector_uses_current_p1_spaces_v_k_kt_w_mass_duality_normals_quadrature_mesh_owner_and_result": (
            calderon_projector_identity_ok
        ),
        "cq_uses_current_symbol_contour_conjugate_transfer_causal_ifft_parseval_passivity_timestep_operator_and_result": (
            cq_physical_closure_identity_ok
        ),
        "fembem_uses_current_reciprocal_transfer_radiation_power_interior_energy_trace_map_frequency_mesh_and_solution": (
            fembem_reciprocity_identity_ok
        ),
        "nonlinear_eigenpairs_use_current_contour_orientation_quadrature_moments_rank_count_residual_biorthogonality_poles_and_result": (
            nonlinear_eigen_contour_identity_ok
        ),
        "cq_acoustics_use_current_bdf2_laplace_contour_weights_passivity_trace_timestep_history_mesh_and_result": (
            cq_acoustic_identity_ok
        ),
        "fembem_autodiff_uses_current_wirtinger_objective_shape_fd_trace_mesh_and_gradient": (
            fembem_autodiff_identity_ok
        ),
        "adaptive_cq_uses_current_timesteps_contour_restart_interpolation_causality_energy_operator_and_result": (
            adaptive_cq_timestep_identity_ok
        ),
        "fembem_shape_derivative_uses_current_morph_normal_velocity_trace_jacobian_fd_mesh_and_result": (
            fembem_shape_derivative_identity_ok
        ),
        "hmatrix_benchmarks_use_current_dense_error_tolerance_rank_memory_complexity_mesh_owners_and_result": (
            hmatrix_benchmark_identity_ok
        ),
        "multifrequency_fembem_adjoints_use_current_weights_objective_quadrature_trace_gradient_fd_mesh_owner_and_result": (
            multifrequency_adjoint_identity_ok
        ),
        "simp_topology_uses_current_density_filter_projection_volume_compliance_adjoint_fd_kkt_mesh_owner_and_result": (
            simp_topology_identity_ok
        ),
        "fembem_model_reduction_uses_current_projection_order_stability_passivity_moments_frequency_error_full_model_mesh_owners_and_result": (
            fembem_model_reduction_identity_ok
        ),
        "nonlinear_fem_uses_current_residual_tangent_newton_linesearch_energy_mesh_owner_and_result": (
            nonlinear_fem_newton_identity_ok
        ),
        "cq_contour_uses_current_nodes_interpolation_aliasing_passivity_reconstruction_time_operator_and_result": (
            cq_contour_reconstruction_identity_ok
        ),
        "johnson_nedelec_uses_current_trace_normals_operators_sign_residual_energy_mesh_owner_and_result": (
            johnson_nedelec_identity_ok
        ),
        "adjoint_hessian_uses_current_design_gradients_constraints_hvp_kkt_fd_model_owner_and_result": (
            adjoint_hessian_identity_ok
        ),
        "cq_acoustic_history_uses_current_causality_passivity_timestep_ztransform_energy_mesh_owner_and_result": (
            cq_acoustic_history_identity_ok
        ),
        "hmatrix_compression_uses_current_clusters_admissibility_ranks_tolerance_matvec_memory_mesh_owner_and_result": (
            hmatrix_compression_identity_ok
        ),
        "lcurve_recomputation_passes": lcurve["status"] == "ok",
        "reported_lcurve_choice_matches": (
            _integer(reported_lcurve, "selected_index") == lcurve["selected_index"]
            and _close(
                _finite_float(reported_lcurve, "selected_alpha"),
                lcurve["selected_alpha"],
            )
        ),
        "morozov_recomputation_passes": morozov["status"] == "ok",
        "reported_morozov_choice_matches": (
            _integer(reported_morozov, "selected_index") == morozov["selected_index"]
            and _close(
                _finite_float(reported_morozov, "selected_alpha"),
                morozov["selected_alpha"],
            )
        ),
        "two_independent_linear_references_close": (
            _positive_integer(crosscheck, "reference_solver_count") >= 2
            and _finite_float(crosscheck, "max_solution_relative_error")
            <= tolerances["max_solution_relative_error"]
            and _finite_float(crosscheck, "max_trace_relative_error")
            <= tolerances["max_trace_relative_error"]
            and _finite_float(
                crosscheck, "max_regularized_objective_relative_error"
            )
            <= tolerances["max_regularized_objective_relative_error"]
            and _finite_float(crosscheck, "zero_alpha_objective_absolute_error")
            <= tolerances["max_zero_alpha_objective_absolute_error"]
        ),
        "deterministic_replay_closes": (
            _positive_integer(replay, "count") >= 2
            and bool(replay.get("selectors_identical", False))
            and _finite_float(replay, "max_relative_error")
            <= tolerances["max_replay_relative_error"]
        ),
    }
    issues = [name for name, passed in checks.items() if not passed]
    return {
        "policy": "regularized_trace_inverse_path_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "lcurve": lcurve,
        "morozov": morozov,
        "metrics": {
            "path_count": len(alphas),
            "max_normal_equation_residual": max(normal_residuals),
            "max_gradient_check_absolute_error": max(gradient_errors),
            "max_solution_relative_error": _finite_float(
                crosscheck, "max_solution_relative_error"
            ),
            "max_replay_relative_error": _finite_float(
                replay, "max_relative_error"
            ),
        },
        "tolerances": tolerances,
        "lesson": (
            "Promote a trace inverse path only after the zero-weight minimum-norm "
            "limit, monotone Tikhonov trade-off, independently recomputed L-curve "
            "and Morozov choices, two linear references, and replay all close."
        ),
    }


def _mapping(parent: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = parent.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"{key} must be a mapping")
    return value


def _float_list(parent: Mapping[str, Any], key: str) -> list[float]:
    value = parent.get(key)
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{key} must be an array")
    out = [float(item) for item in value]
    if not all(math.isfinite(item) for item in out):
        raise ValueError(f"{key} must contain only finite values")
    return out


def _finite_float(parent: Mapping[str, Any], key: str) -> float:
    value = float(parent[key])
    if not math.isfinite(value):
        raise ValueError(f"{key} must be finite")
    return value


def _optional_row_identity_is_aligned(
    path: Mapping[str, Any], expected_length: int
) -> bool:
    names = ("alpha_row_ids", "solution_row_ids", "residual_row_ids")
    if not any(name in path for name in names):
        return True
    rows = []
    for name in names:
        value = path.get(name)
        if (
            not isinstance(value, Sequence)
            or isinstance(value, (str, bytes))
            or len(value) != expected_length
        ):
            return False
        row = [str(item).strip() for item in value]
        if not all(row) or len(set(row)) != expected_length:
            return False
        rows.append(row)
    return rows[0] == rows[1] == rows[2]


def _optional_gradient_check_is_resolved(
    path: Mapping[str, Any], expected_length: int
) -> bool:
    names = (
        "gradient_check_step_sizes",
        "gradient_check_parameter_scales",
        "gradient_check_objective_pair_deltas",
    )
    if not any(name in path for name in names):
        return True
    try:
        steps, scales, objective_deltas = (
            _float_list(path, name) for name in names
        )
    except (KeyError, TypeError, ValueError):
        return False
    if any(
        len(row) != expected_length
        for row in (steps, scales, objective_deltas)
    ):
        return False
    resolution = math.sqrt(math.ulp(1.0))
    return all(
        scale > 0.0
        and step >= resolution * max(scale, 1.0)
        and abs(objective_delta) > 0.0
        for step, scale, objective_delta in zip(steps, scales, objective_deltas)
    )


def _optional_generation_ids_are_aligned(
    path: Mapping[str, Any], expected_length: int
) -> tuple[bool, bool]:
    pairs = (
        ("parameter_generation_ids", "gradient_parameter_generation_ids"),
        ("path_run_generation_ids", "solution_run_generation_ids"),
    )
    results = []
    for left_name, right_name in pairs:
        if left_name not in path and right_name not in path:
            results.append(True)
            continue
        rows = []
        for name in (left_name, right_name):
            value = path.get(name)
            if (
                not isinstance(value, Sequence)
                or isinstance(value, (str, bytes))
                or len(value) != expected_length
            ):
                rows = []
                break
            row = [str(item).strip() for item in value]
            if not all(row):
                rows = []
                break
            rows.append(row)
        results.append(len(rows) == 2 and rows[0] == rows[1])
    return results[0], results[1]


def _optional_quadrature_generation_ids_are_aligned(
    path: Mapping[str, Any], expected_length: int
) -> bool:
    names = (
        "objective_quadrature_generation_ids",
        "gradient_quadrature_generation_ids",
    )
    if not any(name in path for name in names):
        return True
    rows = []
    for name in names:
        value = path.get(name)
        if (
            not isinstance(value, Sequence)
            or isinstance(value, (str, bytes))
            or len(value) != expected_length
        ):
            return False
        row = [str(item).strip() for item in value]
        if not all(row):
            return False
        rows.append(row)
    return rows[0] == rows[1]


def _optional_curvature_parameterization_is_aligned(
    reported_lcurve: Mapping[str, Any],
) -> bool:
    value = reported_lcurve.get("curvature_parameterization")
    if value is None:
        return True
    return (
        isinstance(value, Mapping)
        and value.get("path_coordinate") == "log10_alpha"
        and value.get("curvature_coordinate") == "log10_alpha"
        and value.get("coordinate_transform_recorded") is True
    )


def _optional_design_variable_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("design_variable_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    names = (
        "design_variable_ids",
        "adjoint_gradient_design_variable_ids",
        "finite_difference_design_variable_ids",
    )
    rows = []
    for name in names:
        row = value.get(name)
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)):
            return False
        normalized = [str(item).strip() for item in row]
        if not normalized or not all(normalized) or len(set(normalized)) != len(normalized):
            return False
        rows.append(normalized)
    generations = [
        str(value.get(name, "")).strip()
        for name in (
            "design_generation",
            "adjoint_gradient_design_generation",
            "finite_difference_design_generation",
        )
    ]
    return rows[0] == rows[1] == rows[2] and len(set(generations)) == 1 and bool(
        generations[0]
    )


def _optional_convolution_quadrature_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("convolution_quadrature_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        time_grid_step = _finite_float(value, "time_grid_step_s")
        weight_step = _finite_float(value, "weight_generation_step_s")
    except (KeyError, TypeError, ValueError):
        return False
    time_grid_method = str(value.get("time_grid_method", "")).strip()
    weight_method = str(value.get("weight_generation_method", "")).strip()
    time_grid_generation = str(value.get("time_grid_generation", "")).strip()
    weight_generation = str(value.get("weight_time_grid_generation", "")).strip()
    return (
        time_grid_step > 0.0
        and _close(time_grid_step, weight_step)
        and bool(time_grid_method)
        and time_grid_method == weight_method
        and bool(time_grid_generation)
        and time_grid_generation == weight_generation
    )


def _optional_cq_contour_branch_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("cq_inverse_laplace_contour_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    laplace_ids = value.get("laplace_sample_ids")
    transfer_ids = value.get("transfer_sample_ids")
    branches = value.get("sqrt_branch_conventions")
    if not all(
        isinstance(row, Sequence) and not isinstance(row, (str, bytes))
        for row in (laplace_ids, transfer_ids, branches)
    ):
        return False
    laplace = [str(item).strip() for item in laplace_ids]
    transfer = [str(item).strip() for item in transfer_ids]
    branch_values = [str(item).strip() for item in branches]
    expected_branch = str(
        value.get("inverse_laplace_branch_convention", "")
    ).strip()
    contour_generation = str(value.get("contour_generation", "")).strip()
    return (
        bool(laplace)
        and all(laplace)
        and len(set(laplace)) == len(laplace)
        and transfer == laplace
        and len(branch_values) == len(laplace)
        and bool(expected_branch)
        and all(branch == expected_branch for branch in branch_values)
        and bool(contour_generation)
        and value.get("transfer_sample_contour_generation")
        == contour_generation
    )


def _optional_trace_orientation_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("fembem_trace_orientation_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    volume_digest = str(value.get("volume_mesh_sha256", "")).lower()
    trace_digest = str(value.get("trace_mesh_sha256", "")).lower()
    volume_generation = str(value.get("volume_mesh_generation", "")).strip()
    triangle_digest = str(value.get("trace_boundary_triangle_digest", "")).strip()
    return (
        len(volume_digest) == 64
        and all(character in "0123456789abcdef" for character in volume_digest)
        and trace_digest == volume_digest
        and bool(volume_generation)
        and value.get("trace_mesh_generation") == volume_generation
        and value.get("boundary_orientation_mesh_generation")
        == volume_generation
        and bool(triangle_digest)
        and value.get("oriented_boundary_triangle_digest") == triangle_digest
        and value.get("outward_orientation_verified") is True
    )


def _optional_cq_starting_weights_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("cq_starting_weights_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    method = str(value.get("multistep_method", "")).strip()
    weight_method = str(value.get("convolution_weight_method", "")).strip()
    startup_method = str(value.get("startup_correction_method", "")).strip()
    try:
        order = _positive_integer(value, "multistep_order")
        startup_order = _positive_integer(value, "startup_correction_order")
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("weight_generation", "")).strip()
    return (
        bool(method)
        and weight_method == method
        and startup_method == f"{method}_starting_correction"
        and startup_order == order
        and bool(generation)
        and value.get("startup_correction_weight_generation") == generation
    )


def _optional_hmatrix_cluster_permutation_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("hmatrix_cluster_permutation_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    mesh_digest = str(value.get("boundary_mesh_sha256", "")).lower()
    triangle_digest = str(
        value.get("boundary_triangle_ordering_sha256", "")
    ).lower()
    mesh_generation = str(value.get("boundary_mesh_generation", "")).strip()
    permutation_generation = str(
        value.get("cluster_permutation_generation", "")
    ).strip()
    return (
        len(mesh_digest) == 64
        and all(character in "0123456789abcdef" for character in mesh_digest)
        and len(triangle_digest) == 64
        and all(character in "0123456789abcdef" for character in triangle_digest)
        and value.get("cluster_permutation_triangle_ordering_sha256")
        == triangle_digest
        and bool(mesh_generation)
        and value.get("cluster_permutation_mesh_generation") == mesh_generation
        and bool(permutation_generation)
        and value.get("hmatrix_assembly_permutation_generation")
        == permutation_generation
    )


def _optional_fembem_trace_surface_normal_orientation_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("fembem_trace_surface_normal_orientation_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    mesh_generation = str(value.get("boundary_mesh_generation", "")).strip()
    normal_generation = str(
        value.get("fem_outward_normal_generation", "")
    ).strip()
    fem_triangles = str(value.get("fem_boundary_triangle_sha256", "")).lower()
    bem_triangles = str(value.get("bem_boundary_triangle_sha256", "")).lower()
    try:
        trace_normal_sign = _integer(value, "trace_normal_sign")
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(mesh_generation)
        and value.get("fem_trace_boundary_mesh_generation") == mesh_generation
        and value.get("bem_surface_mesh_generation") == mesh_generation
        and bool(normal_generation)
        and value.get("bem_normal_generation") == normal_generation
        and value.get("trace_operator_normal_generation") == normal_generation
        and value.get("fem_normal_orientation") == "outward"
        and value.get("bem_normal_orientation") == "outward"
        and trace_normal_sign == 1
        and _is_sha256(fem_triangles)
        and bem_triangles == fem_triangles
    )


def _optional_cq_inverse_z_transform_fft_normalization_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("cq_inverse_z_transform_fft_normalization_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        frequency_count = _positive_integer(value, "frequency_sample_count")
        time_count = _positive_integer(value, "time_sample_count")
        inverse_scale = _finite_float(value, "inverse_fft_scale")
    except (KeyError, TypeError, ValueError):
        return False
    frequency_generation = str(
        value.get("cq_frequency_generation", "")
    ).strip()
    convention = str(value.get("transform_convention", "")).strip()
    return (
        frequency_count == time_count
        and value.get("forward_fft_normalization") == "unscaled"
        and value.get("inverse_fft_normalization") == "one_over_n"
        and _close(inverse_scale, 1.0 / frequency_count)
        and bool(frequency_generation)
        and value.get("inverse_transform_frequency_generation")
        == frequency_generation
        and convention
        == "forward_negative_exponent_inverse_positive_exponent"
        and value.get("reconstruction_transform_convention") == convention
    )


def _optional_fembem_sesquilinear_inner_product_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("fembem_sesquilinear_inner_product_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    generation = str(value.get("coupling_generation", "")).strip()
    convention = str(value.get("inner_product_convention", "")).strip()
    operator_digest = str(value.get("operator_metadata_sha256", "")).lower()
    assembled_digest = str(
        value.get("assembled_operator_metadata_sha256", "")
    ).lower()
    return (
        bool(generation)
        and value.get("fem_operator_generation") == generation
        and value.get("bem_operator_generation") == generation
        and value.get("trace_operator_generation") == generation
        and convention == "conjugate_test_linear_trial"
        and value.get("fem_inner_product_convention") == convention
        and value.get("bem_inner_product_convention") == convention
        and value.get("trace_inner_product_convention") == convention
        and value.get("test_argument_conjugated") is True
        and value.get("trial_argument_conjugated") is False
        and _is_sha256(operator_digest)
        and assembled_digest == operator_digest
    )


def _optional_hmatrix_low_rank_tolerance_norm_basis_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("hmatrix_low_rank_tolerance_norm_basis_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        tolerance = _finite_float(value, "compression_tolerance")
        measured_error = _finite_float(value, "measured_relative_error")
    except (KeyError, TypeError, ValueError):
        return False
    assembly_generation = str(
        value.get("hmatrix_assembly_generation", "")
    ).strip()
    tolerance_generation = str(
        value.get("tolerance_calibration_generation", "")
    ).strip()
    norm_basis = str(value.get("tolerance_norm_basis", "")).strip()
    dense_digest = str(value.get("dense_block_sha256", "")).lower()
    source_digest = str(value.get("accepted_block_source_sha256", "")).lower()
    return (
        tolerance >= 0.0
        and measured_error >= 0.0
        and measured_error <= tolerance
        and bool(assembly_generation)
        and value.get("low_rank_factor_generation") == assembly_generation
        and norm_basis
        in {"relative_frobenius", "relative_spectral", "relative_max"}
        and value.get("measured_error_norm_basis") == norm_basis
        and value.get("acceptance_norm_basis") == norm_basis
        and bool(tolerance_generation)
        and value.get("acceptance_tolerance_generation")
        == tolerance_generation
        and _is_sha256(dense_digest)
        and source_digest == dense_digest
    )


def _optional_complex_operator_symmetry_residual_norm_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("complex_operator_symmetry_residual_norm_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        residual_norm = _finite_float(value, "residual_norm")
        reported_norm = _finite_float(value, "reported_residual_norm")
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("operator_generation", "")).strip()
    operator_digest = str(value.get("operator_metadata_sha256", "")).lower()
    residual_digest = str(
        value.get("residual_operator_metadata_sha256", "")
    ).lower()
    return (
        bool(generation)
        and value.get("assembled_operator_generation") == generation
        and value.get("residual_operator_generation") == generation
        and value.get("operator_symmetry_class") == "complex_symmetric"
        and value.get("assembly_symmetry_check") == "transpose"
        and value.get("residual_symmetry_check") == "transpose"
        and value.get("hermitian_symmetry_check_applied") is False
        and value.get("residual_norm_basis") == "complex_euclidean"
        and value.get("reported_residual_norm_basis") == "complex_euclidean"
        and residual_norm >= 0.0
        and _close(reported_norm, residual_norm)
        and _is_sha256(operator_digest)
        and residual_digest == operator_digest
    )


def _optional_hmatrix_admissibility_diameter_metric_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("hmatrix_admissibility_cluster_diameter_metric_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        source_diameter = _finite_float(value, "source_cluster_diameter")
        target_diameter = _finite_float(value, "target_cluster_diameter")
        separation = _finite_float(value, "cluster_separation")
        eta = _finite_float(value, "admissibility_eta")
    except (KeyError, TypeError, ValueError):
        return False
    tree_generation = str(value.get("cluster_tree_generation", "")).strip()
    threshold_generation = str(
        value.get("threshold_calibration_generation", "")
    ).strip()
    diameter_metric = str(value.get("diameter_metric", "")).strip()
    geometry_digest = str(value.get("cluster_geometry_sha256", "")).lower()
    result_digest = str(
        value.get("admissibility_geometry_sha256", "")
    ).lower()
    expected_admissible = max(source_diameter, target_diameter) <= eta * separation
    return (
        bool(tree_generation)
        and value.get("source_cluster_generation") == tree_generation
        and value.get("target_cluster_generation") == tree_generation
        and diameter_metric in {"euclidean_l2", "infinity_linf"}
        and value.get("admissibility_diameter_metric") == diameter_metric
        and value.get("threshold_calibration_diameter_metric") == diameter_metric
        and source_diameter >= 0.0
        and target_diameter >= 0.0
        and separation > 0.0
        and eta > 0.0
        and value.get("admissible") is expected_admissible
        and bool(threshold_generation)
        and value.get("acceptance_threshold_generation") == threshold_generation
        and _is_sha256(geometry_digest)
        and result_digest == geometry_digest
    )


def _optional_cq_convolution_weight_branch_order_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("cq_convolution_weight_ztransform_branch_order_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        order = _positive_integer(value, "multistep_order")
        weight_order = _positive_integer(
            value, "convolution_weight_multistep_order"
        )
        sample_order = _positive_integer(
            value, "transfer_sample_multistep_order"
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("cq_scheme_generation", "")).strip()
    branch = str(value.get("z_transform_branch", "")).strip()
    method = str(value.get("multistep_method", "")).strip().lower()
    scheme_digest = str(value.get("scheme_metadata_sha256", "")).lower()
    weight_digest = str(
        value.get("convolution_weight_scheme_metadata_sha256", "")
    ).lower()
    sample_digest = str(
        value.get("transfer_sample_scheme_metadata_sha256", "")
    ).lower()
    return (
        bool(generation)
        and value.get("convolution_weight_generation") == generation
        and value.get("transfer_sample_generation") == generation
        and branch in {"principal", "continuous_outgoing"}
        and value.get("convolution_weight_z_transform_branch") == branch
        and value.get("transfer_sample_z_transform_branch") == branch
        and method in {"bdf1", "bdf2"}
        and str(value.get("convolution_weight_multistep_method", "")).lower()
        == method
        and str(value.get("transfer_sample_multistep_method", "")).lower()
        == method
        and order in {1, 2}
        and weight_order == order
        and sample_order == order
        and ((method == "bdf1" and order == 1) or (method == "bdf2" and order == 2))
        and _is_sha256(scheme_digest)
        and weight_digest == scheme_digest
        and sample_digest == scheme_digest
    )


def _optional_fembem_trace_boundary_node_generation_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("fembem_trace_matrix_boundary_node_generation_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        boundary_nodes = _positive_integer(value, "boundary_node_count")
        trace_rows = _positive_integer(value, "trace_matrix_row_count")
        trace_nnz = _positive_integer(value, "trace_matrix_nonzero_count")
    except (KeyError, TypeError, ValueError):
        return False
    mesh_generation = str(value.get("volume_mesh_generation", "")).strip()
    node_generation = str(value.get("boundary_node_generation", "")).strip()
    node_digest = str(value.get("boundary_node_ids_sha256", "")).lower()
    fem_digest = str(value.get("fem_boundary_node_ids_sha256", "")).lower()
    bem_digest = str(value.get("bem_boundary_node_ids_sha256", "")).lower()
    trace_digest = str(
        value.get("trace_matrix_boundary_node_ids_sha256", "")
    ).lower()
    return (
        bool(mesh_generation)
        and bool(node_generation)
        and value.get("fem_operator_boundary_node_generation") == node_generation
        and value.get("bem_operator_boundary_node_generation") == node_generation
        and value.get("trace_matrix_boundary_node_generation") == node_generation
        and trace_rows == boundary_nodes
        and trace_nnz == boundary_nodes
        and _is_sha256(node_digest)
        and fem_digest == node_digest
        and bem_digest == node_digest
        and trace_digest == node_digest
    )


def _optional_cq_frequency_contour_radius_damping_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("cq_frequency_contour_radius_damping_generation_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        time_step = _finite_float(value, "time_step_s")
        sample_time_step = _finite_float(value, "transfer_sample_time_step_s")
        sample_count = _positive_integer(value, "sample_count")
        transfer_count = _positive_integer(value, "transfer_sample_count")
        contour_radius = _finite_float(value, "contour_radius")
        inverse_radius = _finite_float(value, "inverse_transform_contour_radius")
        damping = _finite_float(value, "damping_factor")
        inverse_damping = _finite_float(value, "inverse_transform_damping_factor")
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("frequency_sampling_generation", "")).strip()
    contour_digest = str(value.get("contour_metadata_sha256", "")).lower()
    inverse_digest = str(
        value.get("inverse_transform_contour_metadata_sha256", "")
    ).lower()
    return (
        bool(generation)
        and value.get("transfer_sample_frequency_generation") == generation
        and value.get("contour_metadata_frequency_generation") == generation
        and value.get("inverse_transform_frequency_generation") == generation
        and time_step > 0.0
        and _close(sample_time_step, time_step)
        and transfer_count == sample_count
        and 0.0 < contour_radius <= 1.0
        and _close(inverse_radius, contour_radius)
        and 0.0 < damping <= 1.0
        and _close(inverse_damping, damping)
        and _is_sha256(contour_digest)
        and inverse_digest == contour_digest
    )


def _optional_hmatrix_aca_pivot_cluster_permutation_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get("hmatrix_aca_pivot_cluster_permutation_generation_identity")
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        row_pivots = _positive_integer_sequence(value, "row_pivot_global_ids")
        aca_row_pivots = _positive_integer_sequence(
            value, "aca_row_pivot_global_ids"
        )
        column_pivots = _positive_integer_sequence(
            value, "column_pivot_global_ids"
        )
        aca_column_pivots = _positive_integer_sequence(
            value, "aca_column_pivot_global_ids"
        )
    except (KeyError, TypeError, ValueError):
        return False
    tree_generation = str(value.get("cluster_tree_generation", "")).strip()
    source_generation = str(
        value.get("source_cluster_permutation_generation", "")
    ).strip()
    target_generation = str(
        value.get("target_cluster_permutation_generation", "")
    ).strip()
    permutation_digest = str(value.get("cluster_permutation_sha256", "")).lower()
    aca_digest = str(value.get("aca_cluster_permutation_sha256", "")).lower()
    return (
        bool(tree_generation)
        and bool(source_generation)
        and target_generation == source_generation
        and value.get("aca_source_cluster_permutation_generation")
        == source_generation
        and value.get("aca_target_cluster_permutation_generation")
        == target_generation
        and len(row_pivots) == len(aca_row_pivots) > 0
        and len(column_pivots) == len(aca_column_pivots) > 0
        and len(set(row_pivots)) == len(row_pivots)
        and len(set(column_pivots)) == len(column_pivots)
        and aca_row_pivots == row_pivots
        and aca_column_pivots == column_pivots
        and _is_sha256(permutation_digest)
        and aca_digest == permutation_digest
    )


def _optional_cq_conjugate_frequency_bin_ownership_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_inverse_transform_conjugate_frequency_bin_ownership_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        sample_count = _positive_integer(value, "sample_count")
        inverse_count = _positive_integer(value, "inverse_transform_sample_count")
        positive_bins = _positive_integer_sequence(
            value, "positive_frequency_bin_indices"
        )
        conjugate_bins = _positive_integer_sequence(
            value, "conjugate_frequency_bin_indices"
        )
        inverse_conjugate_bins = _positive_integer_sequence(
            value, "inverse_transform_conjugate_frequency_bin_indices"
        )
        dc_bin = _integer(value, "dc_bin_index")
        nyquist_bin = _integer(value, "nyquist_bin_index")
    except (KeyError, TypeError, ValueError):
        return False
    transfer_generation = str(
        value.get("transfer_sample_generation", "")
    ).strip()
    bin_digest = str(value.get("frequency_bin_map_sha256", "")).lower()
    return (
        bool(transfer_generation)
        and value.get("frequency_bin_transfer_sample_generation")
        == transfer_generation
        and value.get("inverse_transform_transfer_sample_generation")
        == transfer_generation
        and sample_count == inverse_count
        and sample_count % 2 == 0
        and len(positive_bins) == len(conjugate_bins) > 0
        and len(set(positive_bins)) == len(positive_bins)
        and len(set(conjugate_bins)) == len(conjugate_bins)
        and all(0 < index < sample_count // 2 for index in positive_bins)
        and conjugate_bins
        == tuple(sample_count - index for index in positive_bins)
        and inverse_conjugate_bins == conjugate_bins
        and dc_bin == 0
        and nyquist_bin == sample_count // 2
        and value.get("real_response_conjugate_symmetry") is True
        and value.get("inverse_transform_conjugate_symmetry") is True
        and _is_sha256(bin_digest)
        and str(
            value.get("inverse_transform_frequency_bin_map_sha256", "")
        ).lower()
        == bin_digest
    )


def _optional_hmatrix_bounding_box_mesh_scale_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "hmatrix_block_cluster_bounding_box_mesh_scale_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        mesh_scale = _finite_float(value, "mesh_length_scale_to_m")
        bbox_scale = _finite_float(
            value, "cluster_bounding_box_length_scale_to_m"
        )
        source_bbox = tuple(
            float(item) for item in value["source_cluster_bounding_box_m"]
        )
        block_bbox = tuple(
            float(item) for item in value["block_source_cluster_bounding_box_m"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    mesh_generation = str(value.get("mesh_generation", "")).strip()
    scale_generation = str(value.get("mesh_scale_generation", "")).strip()
    bbox_digest = str(value.get("cluster_bounding_box_sha256", "")).lower()
    return (
        bool(mesh_generation)
        and value.get("cluster_tree_mesh_generation") == mesh_generation
        and value.get("block_cluster_mesh_generation") == mesh_generation
        and bool(scale_generation)
        and value.get("cluster_bounding_box_mesh_scale_generation")
        == scale_generation
        and value.get("block_admissibility_mesh_scale_generation")
        == scale_generation
        and mesh_scale > 0.0
        and _close(bbox_scale, mesh_scale)
        and len(source_bbox) == 6
        and all(math.isfinite(item) for item in source_bbox)
        and all(source_bbox[index] <= source_bbox[index + 3] for index in range(3))
        and block_bbox == source_bbox
        and _is_sha256(bbox_digest)
        and str(
            value.get("block_admissibility_bounding_box_sha256", "")
        ).lower()
        == bbox_digest
    )


def _optional_cq_contour_fft_timestep_generation_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_laplace_contour_fft_scaling_timestep_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        sample_count = _positive_integer(value, "sample_count")
        fft_count = _positive_integer(value, "fft_sample_count")
        time_step = _finite_float(value, "time_step_s")
        inverse_time_step = _finite_float(
            value, "inverse_transform_time_step_s"
        )
        fft_scaling = _finite_float(value, "fft_scaling")
        inverse_fft_scaling = _finite_float(
            value, "inverse_transform_fft_scaling"
        )
        contour_nodes = tuple(
            tuple(float(component) for component in node)
            for node in value["laplace_contour_nodes_re_im"]
        )
        inverse_nodes = tuple(
            tuple(float(component) for component in node)
            for node in value[
                "inverse_transform_laplace_contour_nodes_re_im"
            ]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("cq_generation", "")).strip()
    digest = str(value.get("contour_fft_timestep_sha256", "")).lower()
    return (
        bool(generation)
        and value.get("laplace_contour_cq_generation") == generation
        and value.get("fft_scaling_cq_generation") == generation
        and value.get("timestep_cq_generation") == generation
        and value.get("inverse_transform_cq_generation") == generation
        and sample_count == fft_count
        and time_step > 0.0
        and _close(inverse_time_step, time_step)
        and fft_scaling > 0.0
        and _close(fft_scaling, 1.0 / sample_count)
        and _close(inverse_fft_scaling, fft_scaling)
        and bool(contour_nodes)
        and all(
            len(node) == 2 and all(math.isfinite(component) for component in node)
            for node in contour_nodes
        )
        and inverse_nodes == contour_nodes
        and _is_sha256(digest)
        and str(
            value.get("inverse_transform_contour_fft_timestep_sha256", "")
        ).lower()
        == digest
    )


def _optional_hmatrix_aca_tolerance_norm_rank_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "hmatrix_aca_tolerance_norm_block_rank_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        block_ids = _positive_integer_sequence(value, "block_ids")
        aca_block_ids = _positive_integer_sequence(value, "aca_block_ids")
        tolerances = tuple(float(item) for item in value["aca_tolerances"])
        applied_tolerances = tuple(
            float(item) for item in value["applied_aca_tolerances"]
        )
        norm_estimates = tuple(
            float(item) for item in value["block_norm_estimates"]
        )
        aca_norm_estimates = tuple(
            float(item) for item in value["aca_block_norm_estimates"]
        )
        ranks = _positive_integer_sequence(value, "retained_ranks")
        assembled_ranks = _positive_integer_sequence(
            value, "assembled_block_ranks"
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("block_generation", "")).strip()
    table_digest = str(value.get("aca_block_table_sha256", "")).lower()
    return (
        bool(generation)
        and value.get("aca_tolerance_block_generation") == generation
        and value.get("norm_estimate_block_generation") == generation
        and value.get("retained_rank_block_generation") == generation
        and value.get("assembled_block_generation") == generation
        and bool(block_ids)
        and len(set(block_ids)) == len(block_ids)
        and aca_block_ids == block_ids
        and len(tolerances) == len(norm_estimates) == len(ranks) == len(block_ids)
        and all(math.isfinite(item) and 0.0 < item < 1.0 for item in tolerances)
        and applied_tolerances == tolerances
        and all(math.isfinite(item) and item > 0.0 for item in norm_estimates)
        and aca_norm_estimates == norm_estimates
        and assembled_ranks == ranks
        and _is_sha256(table_digest)
        and str(value.get("assembled_aca_block_table_sha256", "")).lower()
        == table_digest
    )


def _optional_fembem_trace_orientation_normal_node_map_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "fembem_trace_orientation_normal_node_map_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        triangle_ids = _positive_integer_sequence(value, "boundary_triangle_ids")
        trace_ids = _positive_integer_sequence(value, "trace_triangle_ids")
        orientations = tuple(int(item) for item in value["triangle_orientations"])
        applied_orientations = tuple(
            int(item) for item in value["applied_triangle_orientations"]
        )
        normal_digests = tuple(str(item).lower() for item in value["normal_sha256"])
        applied_normal_digests = tuple(
            str(item).lower() for item in value["applied_normal_sha256"]
        )
        volume_nodes = _positive_integer_sequence(value, "volume_boundary_node_ids")
        trace_nodes = _positive_integer_sequence(value, "trace_volume_node_ids")
    except (KeyError, TypeError, ValueError):
        return False
    volume_generation = str(value.get("volume_mesh_generation", "")).strip()
    boundary_generation = str(value.get("boundary_mesh_generation", "")).strip()
    map_digest = str(value.get("trace_map_sha256", "")).lower()
    return (
        bool(volume_generation)
        and value.get("boundary_trace_volume_mesh_generation") == volume_generation
        and bool(boundary_generation)
        and all(
            value.get(key) == boundary_generation
            for key in (
                "triangle_orientation_boundary_mesh_generation",
                "normal_boundary_mesh_generation",
                "node_map_boundary_mesh_generation",
                "fem_trace_boundary_mesh_generation",
                "bem_trace_boundary_mesh_generation",
            )
        )
        and bool(triangle_ids)
        and len(set(triangle_ids)) == len(triangle_ids)
        and trace_ids == triangle_ids
        and len(orientations) == len(triangle_ids)
        and all(item in {-1, 1} for item in orientations)
        and applied_orientations == orientations
        and len(normal_digests) == len(triangle_ids)
        and all(_is_sha256(item) for item in normal_digests)
        and applied_normal_digests == normal_digests
        and bool(volume_nodes)
        and len(set(volume_nodes)) == len(volume_nodes)
        and trace_nodes == volume_nodes
        and _is_sha256(map_digest)
        and str(value.get("assembled_trace_map_sha256", "")).lower() == map_digest
    )


def _optional_cq_causality_conjugate_contour_pair_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_causality_conjugate_symmetry_contour_pair_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        positive_ids = _positive_integer_sequence(value, "positive_frequency_ids")
        conjugate_ids = _positive_integer_sequence(value, "conjugate_frequency_ids")
        inverse_ids = _positive_integer_sequence(
            value, "inverse_transform_conjugate_frequency_ids"
        )
        contour_indices = _positive_integer_sequence(value, "contour_indices")
        inverse_contour = _positive_integer_sequence(
            value, "inverse_transform_contour_indices"
        )
        causality_window = tuple(int(item) for item in value["causality_window_samples"])
        inverse_window = tuple(
            int(item) for item in value["inverse_transform_causality_window_samples"]
        )
        precausal_max = _finite_float(value, "precausal_max_abs")
        precausal_tolerance = _finite_float(value, "precausal_tolerance")
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("cq_generation", "")).strip()
    pair_digest = str(value.get("cq_pair_table_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "frequency_pair_cq_generation",
                "contour_index_cq_generation",
                "causality_window_cq_generation",
                "inverse_transform_cq_generation",
            )
        )
        and bool(positive_ids)
        and len(positive_ids) == len(conjugate_ids)
        and inverse_ids == conjugate_ids
        and contour_indices == inverse_contour
        and set(positive_ids).issubset(contour_indices)
        and set(conjugate_ids).issubset(contour_indices)
        and len(causality_window) == 2
        and 0 <= causality_window[0] < causality_window[1]
        and inverse_window == causality_window
        and value.get("real_time_response") is True
        and value.get("inverse_transform_conjugate_symmetry") is True
        and 0.0 <= precausal_max <= precausal_tolerance
        and precausal_tolerance > 0.0
        and _is_sha256(pair_digest)
        and str(value.get("inverse_transform_pair_table_sha256", "")).lower()
        == pair_digest
    )


def _optional_p1_boundary_mass_trace_row_node_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "p1_boundary_mass_trace_row_node_mesh_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        boundary_nodes = _positive_integer_sequence(value, "boundary_node_ids")
        mass_nodes = _positive_integer_sequence(value, "mass_row_node_ids")
        trace_nodes = _positive_integer_sequence(value, "trace_row_node_ids")
        mass_rows = _positive_integer_sequence(value, "boundary_mass_row_ids")
        assembled_mass_rows = _positive_integer_sequence(
            value, "assembled_boundary_mass_row_ids"
        )
        trace_rows = _positive_integer_sequence(value, "trace_row_ids")
        assembled_trace_rows = _positive_integer_sequence(
            value, "assembled_trace_row_ids"
        )
    except (KeyError, TypeError, ValueError):
        return False
    volume_generation = str(value.get("volume_mesh_generation", "")).strip()
    boundary_generation = str(value.get("boundary_mesh_generation", "")).strip()
    mass_digest = str(value.get("boundary_mass_matrix_sha256", "")).lower()
    trace_digest = str(value.get("p1_trace_matrix_sha256", "")).lower()
    return (
        bool(volume_generation)
        and value.get("boundary_volume_mesh_generation") == volume_generation
        and bool(boundary_generation)
        and all(
            value.get(key) == boundary_generation
            for key in (
                "mass_boundary_mesh_generation",
                "trace_boundary_mesh_generation",
                "node_map_boundary_mesh_generation",
            )
        )
        and bool(boundary_nodes)
        and len(set(boundary_nodes)) == len(boundary_nodes)
        and mass_nodes == boundary_nodes
        and trace_nodes == boundary_nodes
        and len(mass_rows) == len(boundary_nodes)
        and len(set(mass_rows)) == len(mass_rows)
        and assembled_mass_rows == mass_rows
        and trace_rows == mass_rows
        and assembled_trace_rows == trace_rows
        and _is_sha256(mass_digest)
        and str(value.get("assembled_boundary_mass_matrix_sha256", "")).lower()
        == mass_digest
        and _is_sha256(trace_digest)
        and str(value.get("assembled_p1_trace_matrix_sha256", "")).lower()
        == trace_digest
    )


def _optional_cq_restart_history_weight_segment_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_restart_history_weight_segment_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        weights = tuple(
            tuple(float(component) for component in row)
            for row in value["convolution_weights_re_im"]
        )
        restart_weights = tuple(
            tuple(float(component) for component in row)
            for row in value["restart_convolution_weights_re_im"]
        )
        history = tuple(
            str(item).lower() for item in value["history_vector_sha256"]
        )
        restart_history = tuple(
            str(item).lower()
            for item in value["restart_history_vector_sha256"]
        )
        segment_offset = _integer(value, "segment_offset")
        restart_segment_offset = _integer(value, "restart_segment_offset")
        time_grid = tuple(float(item) for item in value["time_grid_s"])
        restart_time_grid = tuple(
            float(item) for item in value["restart_time_grid_s"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    cq_generation = str(value.get("cq_generation", "")).strip()
    checkpoint_generation = str(value.get("checkpoint_generation", "")).strip()
    table_digest = str(value.get("cq_restart_table_sha256", "")).lower()
    return (
        bool(cq_generation)
        and value.get("restart_cq_generation") == cq_generation
        and bool(checkpoint_generation)
        and all(
            value.get(key) == checkpoint_generation
            for key in (
                "weight_checkpoint_generation",
                "history_checkpoint_generation",
                "segment_checkpoint_generation",
                "time_grid_checkpoint_generation",
            )
        )
        and bool(weights)
        and all(
            len(row) == 2 and all(math.isfinite(component) for component in row)
            for row in weights
        )
        and restart_weights == weights
        and len(history) == len(weights)
        and all(_is_sha256(item) for item in history)
        and restart_history == history
        and segment_offset >= 0
        and restart_segment_offset == segment_offset
        and len(time_grid) == len(weights) + 1
        and all(math.isfinite(item) and item >= 0.0 for item in time_grid)
        and all(right > left for left, right in zip(time_grid, time_grid[1:]))
        and restart_time_grid == time_grid
        and _is_sha256(table_digest)
        and str(value.get("assembled_cq_restart_table_sha256", "")).lower()
        == table_digest
    )


def _optional_simscape_file_solid_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "simscape_file_solid_geometry_density_inertia_frame_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        density = _finite_float(value, "density_kg_per_m3")
        result_density = _finite_float(value, "result_density_kg_per_m3")
        center = tuple(float(item) for item in value["center_of_mass_m"])
        result_center = tuple(float(item) for item in value["result_center_of_mass_m"])
        inertia = tuple(
            tuple(float(component) for component in row)
            for row in value["inertia_tensor_kg_m2"]
        )
        result_inertia = tuple(
            tuple(float(component) for component in row)
            for row in value["result_inertia_tensor_kg_m2"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("assembly_generation", "")).strip()
    geometry_digest = str(value.get("geometry_sha256", "")).lower()
    table_digest = str(value.get("file_solid_table_sha256", "")).lower()
    inertia_is_symmetric = (
        len(inertia) == 3
        and all(len(row) == 3 for row in inertia)
        and all(
            math.isfinite(inertia[row][column])
            and math.isclose(
                inertia[row][column], inertia[column][row], rel_tol=0.0, abs_tol=1e-15
            )
            for row in range(3)
            for column in range(3)
        )
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "geometry_assembly_generation",
                "density_assembly_generation",
                "inertia_assembly_generation",
                "frame_assembly_generation",
                "result_assembly_generation",
            )
        )
        and bool(str(value.get("geometry_file", "")).strip())
        and value.get("result_geometry_file") == value.get("geometry_file")
        and _is_sha256(geometry_digest)
        and str(value.get("result_geometry_sha256", "")).lower() == geometry_digest
        and density > 0.0
        and result_density == density
        and bool(str(value.get("center_of_mass_frame", "")).strip())
        and value.get("result_center_of_mass_frame") == value.get("center_of_mass_frame")
        and len(center) == 3
        and all(math.isfinite(component) for component in center)
        and result_center == center
        and inertia_is_symmetric
        and result_inertia == inertia
        and _is_sha256(table_digest)
        and str(value.get("result_file_solid_table_sha256", "")).lower()
        == table_digest
    )


def _optional_multibody_xml_import_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "multibody_xml_joint_axis_transform_unit_geometry_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        joint_ids = tuple(str(item) for item in value["joint_ids"])
        result_joint_ids = tuple(str(item) for item in value["result_joint_ids"])
        axes = tuple(tuple(float(component) for component in row) for row in value["joint_axes"])
        result_axes = tuple(
            tuple(float(component) for component in row)
            for row in value["result_joint_axes"]
        )
        transform = tuple(
            tuple(float(component) for component in row)
            for row in value["rigid_transforms"]
        )
        result_transform = tuple(
            tuple(float(component) for component in row)
            for row in value["result_rigid_transforms"]
        )
        scale = _finite_float(value, "length_scale_to_m")
        result_scale = _finite_float(value, "result_length_scale_to_m")
        geometry_files = tuple(str(item) for item in value["geometry_files"])
        result_geometry_files = tuple(str(item) for item in value["result_geometry_files"])
        geometry_digests = tuple(str(item).lower() for item in value["geometry_digests"])
        result_geometry_digests = tuple(
            str(item).lower() for item in value["result_geometry_digests"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("assembly_generation", "")).strip()
    table_digest = str(value.get("xml_import_table_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "joint_axis_assembly_generation",
                "rigid_transform_assembly_generation",
                "length_unit_assembly_generation",
                "geometry_file_assembly_generation",
                "result_assembly_generation",
            )
        )
        and bool(joint_ids)
        and all(joint_ids)
        and len(set(joint_ids)) == len(joint_ids)
        and result_joint_ids == joint_ids
        and len(axes) == len(joint_ids)
        and all(
            len(axis) == 3
            and all(math.isfinite(component) for component in axis)
            and any(component != 0.0 for component in axis)
            for axis in axes
        )
        and result_axes == axes
        and len(transform) == 4
        and all(len(row) == 4 for row in transform)
        and all(math.isfinite(component) for row in transform for component in row)
        and result_transform == transform
        and bool(str(value.get("length_unit", "")).strip())
        and value.get("result_length_unit") == value.get("length_unit")
        and scale > 0.0
        and result_scale == scale
        and bool(geometry_files)
        and all(geometry_files)
        and result_geometry_files == geometry_files
        and len(geometry_digests) == len(geometry_files)
        and all(_is_sha256(digest) for digest in geometry_digests)
        and result_geometry_digests == geometry_digests
        and _is_sha256(table_digest)
        and str(value.get("result_xml_import_table_sha256", "")).lower()
        == table_digest
    )


def _optional_parallel_pool_identity_is_aligned(summary: Mapping[str, Any]) -> bool:
    value = summary.get(
        "parallel_pool_worker_path_device_rng_code_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        worker_ids = tuple(int(item) for item in value["worker_ids"])
        result_worker_ids = tuple(int(item) for item in value["result_worker_ids"])
        paths = tuple(str(item) for item in value["worker_code_paths"])
        result_paths = tuple(str(item) for item in value["result_worker_code_paths"])
        devices = tuple(str(item) for item in value["device_assignments"])
        result_devices = tuple(str(item) for item in value["result_device_assignments"])
        seeds = tuple(int(item) for item in value["random_stream_seeds"])
        result_seeds = tuple(int(item) for item in value["result_random_stream_seeds"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("pool_generation", "")).strip()
    code_digest = str(value.get("worker_code_sha256", "")).lower()
    result_digest = str(value.get("parallel_result_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "worker_path_pool_generation",
                "device_pool_generation",
                "rng_pool_generation",
                "code_pool_generation",
                "result_pool_generation",
            )
        )
        and bool(worker_ids)
        and all(item > 0 for item in worker_ids)
        and len(set(worker_ids)) == len(worker_ids)
        and result_worker_ids == worker_ids
        and len(paths) == len(worker_ids)
        and all(paths)
        and result_paths == paths
        and len(devices) == len(worker_ids)
        and all(devices)
        and result_devices == devices
        and len(seeds) == len(worker_ids)
        and all(item >= 0 for item in seeds)
        and len(set(seeds)) == len(seeds)
        and result_seeds == seeds
        and _is_sha256(code_digest)
        and str(value.get("result_worker_code_sha256", "")).lower() == code_digest
        and _is_sha256(result_digest)
        and str(value.get("assembled_parallel_result_sha256", "")).lower()
        == result_digest
    )


def _optional_autodiff_tape_identity_is_aligned(summary: Mapping[str, Any]) -> bool:
    value = summary.get(
        "autodiff_tape_variable_order_mesh_objective_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        variable_ids = tuple(str(item) for item in value["variable_ids"])
        gradient_ids = tuple(str(item) for item in value["gradient_variable_ids"])
        scale = float(value["objective_scale"])
        gradient_scale = float(value["gradient_objective_scale"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("tape_generation", "")).strip()
    mesh_digest = str(value.get("mesh_sha256", "")).lower()
    primal_digest = str(value.get("primal_state_sha256", "")).lower()
    gradient_digest = str(value.get("gradient_table_sha256", "")).lower()
    objective_id = str(value.get("objective_id", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "variable_order_tape_generation",
                "mesh_tape_generation",
                "objective_scaling_tape_generation",
                "primal_solve_tape_generation",
                "gradient_result_tape_generation",
            )
        )
        and bool(variable_ids)
        and all(variable_ids)
        and len(set(variable_ids)) == len(variable_ids)
        and gradient_ids == variable_ids
        and _is_sha256(mesh_digest)
        and str(value.get("gradient_mesh_sha256", "")).lower() == mesh_digest
        and bool(objective_id)
        and value.get("gradient_objective_id") == objective_id
        and math.isfinite(scale)
        and scale > 0.0
        and gradient_scale == scale
        and _is_sha256(primal_digest)
        and str(value.get("gradient_primal_state_sha256", "")).lower()
        == primal_digest
        and _is_sha256(gradient_digest)
        and str(value.get("reported_gradient_table_sha256", "")).lower()
        == gradient_digest
    )


def _optional_fembem_interface_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "fembem_trace_normal_interface_node_order_unit_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        node_ids = tuple(int(item) for item in value["interface_node_ids"])
        result_node_ids = tuple(
            int(item) for item in value["result_interface_node_ids"]
        )
        triangles = tuple(
            tuple(int(item) for item in row) for row in value["boundary_triangles"]
        )
        result_triangles = tuple(
            tuple(int(item) for item in row)
            for row in value["result_boundary_triangles"]
        )
        normals = tuple(
            tuple(float(item) for item in row) for row in value["outward_normals"]
        )
        result_normals = tuple(
            tuple(float(item) for item in row)
            for row in value["result_outward_normals"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("coupling_generation", "")).strip()
    units = value.get("physical_units")
    mesh_digest = str(value.get("interface_mesh_sha256", "")).lower()
    operator_digest = str(value.get("coupled_operator_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "trace_coupling_generation",
                "normal_coupling_generation",
                "node_order_coupling_generation",
                "unit_coupling_generation",
                "operator_coupling_generation",
                "result_coupling_generation",
            )
        )
        and value.get("trace_orientation") == "volume_to_boundary"
        and value.get("result_trace_orientation") == value.get("trace_orientation")
        and value.get("outward_normal_convention") == "exterior_from_volume"
        and value.get("result_outward_normal_convention")
        == value.get("outward_normal_convention")
        and bool(node_ids)
        and node_ids[0] >= 1
        and all(left < right for left, right in zip(node_ids, node_ids[1:]))
        and result_node_ids == node_ids
        and bool(triangles)
        and all(
            len(row) == 3 and len(set(row)) == 3 and all(item in node_ids for item in row)
            for row in triangles
        )
        and result_triangles == triangles
        and len(normals) == len(triangles)
        and all(
            len(row) == 3
            and all(math.isfinite(item) for item in row)
            and math.isclose(sum(item * item for item in row), 1.0, abs_tol=1.0e-12)
            for row in normals
        )
        and result_normals == normals
        and units == {"pressure": "Pa", "normal_velocity": "m/s"}
        and value.get("result_physical_units") == units
        and _is_sha256(mesh_digest)
        and str(value.get("result_interface_mesh_sha256", "")).lower()
        == mesh_digest
        and _is_sha256(operator_digest)
        and str(value.get("result_coupled_operator_sha256", "")).lower()
        == operator_digest
    )


def _optional_cq_time_history_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_contour_weight_startup_causality_window_result_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        contour = tuple(
            tuple(float(item) for item in row) for row in value["contour_points_ri"]
        )
        result_contour = tuple(
            tuple(float(item) for item in row)
            for row in value["result_contour_points_ri"]
        )
        weights = tuple(
            tuple(float(item) for item in row) for row in value["cq_weights_ri"]
        )
        result_weights = tuple(
            tuple(float(item) for item in row)
            for row in value["result_cq_weights_ri"]
        )
        startup = tuple(
            tuple(float(item) for item in row) for row in value["startup_weights_ri"]
        )
        result_startup = tuple(
            tuple(float(item) for item in row)
            for row in value["result_startup_weights_ri"]
        )
        times = tuple(float(item) for item in value["time_samples_s"])
        result_times = tuple(float(item) for item in value["result_time_samples_s"])
        window = tuple(float(item) for item in value["causality_window_s"])
        result_window = tuple(
            float(item) for item in value["result_causality_window_s"]
        )
        prehistory = float(value["prehistory_norm"])
        result_prehistory = float(value["result_prehistory_norm"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("cq_generation", "")).strip()
    digest = str(value.get("cq_result_sha256", "")).lower()
    complex_rows = contour + weights + startup
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "contour_cq_generation",
                "weight_cq_generation",
                "startup_cq_generation",
                "causality_window_cq_generation",
                "time_grid_cq_generation",
                "result_cq_generation",
            )
        )
        and value.get("method") == "BDF2"
        and value.get("result_method") == value.get("method")
        and len(contour) >= 4
        and len(weights) == len(contour)
        and len(startup) == 2
        and all(
            len(row) == 2 and all(math.isfinite(item) for item in row)
            for row in complex_rows
        )
        and result_contour == contour
        and result_weights == weights
        and result_startup == startup
        and len(times) == len(contour)
        and times[0] == 0.0
        and all(math.isfinite(item) for item in times)
        and all(left < right for left, right in zip(times, times[1:]))
        and result_times == times
        and len(window) == 2
        and window == (times[0], times[-1])
        and result_window == window
        and math.isfinite(prehistory)
        and prehistory == 0.0
        and result_prehistory == prehistory
        and _is_sha256(digest)
        and str(value.get("reported_cq_result_sha256", "")).lower() == digest
    )


def _optional_hmatrix_block_tree_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "hmatrix_block_tree_admissibility_permutation_tolerance_kernel_mesh_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        shape = tuple(int(item) for item in value["matrix_shape"])
        result_shape = tuple(int(item) for item in value["result_matrix_shape"])
        row_permutation = tuple(int(item) for item in value["row_permutation"])
        result_row_permutation = tuple(
            int(item) for item in value["result_row_permutation"]
        )
        column_permutation = tuple(
            int(item) for item in value["column_permutation"]
        )
        result_column_permutation = tuple(
            int(item) for item in value["result_column_permutation"]
        )
        eta = float(value["admissibility_eta"])
        result_eta = float(value["result_admissibility_eta"])
        tolerance = float(value["relative_tolerance"])
        result_tolerance = float(value["result_relative_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("hmatrix_generation", "")).strip()
    block_digest = str(value.get("block_tree_sha256", "")).lower()
    mesh_digest = str(value.get("boundary_mesh_sha256", "")).lower()
    result_digest = str(value.get("hmatrix_result_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "block_tree_hmatrix_generation",
                "admissibility_hmatrix_generation",
                "permutation_hmatrix_generation",
                "tolerance_hmatrix_generation",
                "kernel_hmatrix_generation",
                "mesh_hmatrix_generation",
                "result_hmatrix_generation",
            )
        )
        and len(shape) == 2
        and all(item > 0 for item in shape)
        and result_shape == shape
        and _is_sha256(block_digest)
        and str(value.get("result_block_tree_sha256", "")).lower() == block_digest
        and value.get("admissibility_rule") == "diameter_le_eta_distance"
        and value.get("result_admissibility_rule")
        == value.get("admissibility_rule")
        and math.isfinite(eta)
        and eta > 0.0
        and result_eta == eta
        and sorted(row_permutation) == list(range(shape[0]))
        and result_row_permutation == row_permutation
        and sorted(column_permutation) == list(range(shape[1]))
        and result_column_permutation == column_permutation
        and math.isfinite(tolerance)
        and 0.0 < tolerance < 1.0
        and result_tolerance == tolerance
        and bool(str(value.get("kernel_id", "")).strip())
        and value.get("result_kernel_id") == value.get("kernel_id")
        and _is_sha256(mesh_digest)
        and str(value.get("result_boundary_mesh_sha256", "")).lower()
        == mesh_digest
        and _is_sha256(result_digest)
        and str(value.get("reported_hmatrix_result_sha256", "")).lower()
        == result_digest
    )


def _optional_ad_gradient_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "ad_parameter_tape_material_operator_mesh_objective_primal_gradient_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        names = tuple(str(item) for item in value["parameter_names"])
        result_names = tuple(str(item) for item in value["result_parameter_names"])
        parameters = tuple(float(item) for item in value["parameter_values"])
        result_parameters = tuple(
            float(item) for item in value["result_parameter_values"]
        )
        gradient = tuple(float(item) for item in value["ad_gradient"])
        reference = tuple(
            float(item) for item in value["finite_difference_gradient"]
        )
        reported_error = float(value["maximum_gradient_relative_error"])
        tolerance = float(value["gradient_relative_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("ad_generation", "")).strip()
    relative_errors = [
        abs(left - right) / max(abs(left), abs(right), 1.0e-300)
        for left, right in zip(gradient, reference)
    ]
    digests = (
        ("parameter_tape_sha256", "result_parameter_tape_sha256"),
        ("assembled_operator_sha256", "result_assembled_operator_sha256"),
        ("mesh_sha256", "result_mesh_sha256"),
        ("primal_solution_sha256", "result_primal_solution_sha256"),
        ("gradient_result_sha256", "reported_gradient_result_sha256"),
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "parameter_tape_ad_generation",
                "material_law_ad_generation",
                "operator_ad_generation",
                "mesh_ad_generation",
                "objective_ad_generation",
                "primal_ad_generation",
                "gradient_ad_generation",
                "result_ad_generation",
            )
        )
        and bool(names)
        and all(names)
        and len(set(names)) == len(names)
        and result_names == names
        and len(parameters) == len(names)
        and all(math.isfinite(item) for item in parameters)
        and result_parameters == parameters
        and bool(str(value.get("material_law_id", "")).strip())
        and value.get("result_material_law_id") == value.get("material_law_id")
        and bool(str(value.get("objective_id", "")).strip())
        and value.get("result_objective_id") == value.get("objective_id")
        and all(
            _is_sha256(str(value.get(source, "")).lower())
            and str(value.get(target, "")).lower()
            == str(value.get(source, "")).lower()
            for source, target in digests
        )
        and len(gradient) == len(reference) == len(names)
        and all(math.isfinite(item) for item in gradient + reference)
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and math.isfinite(reported_error)
        and relative_errors
        and math.isclose(reported_error, max(relative_errors), rel_tol=1.0e-6)
        and reported_error <= tolerance
    )


def _optional_cq_transfer_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_contour_radius_timestep_laplace_branch_transfer_operator_inverse_transform_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        radius = float(value["contour_radius"])
        result_radius = float(value["result_contour_radius"])
        timestep = float(value["time_step_s"])
        result_timestep = float(value["result_time_step_s"])
        points = tuple(
            tuple(float(item) for item in pair)
            for pair in value["laplace_points_ri"]
        )
        result_points = tuple(
            tuple(float(item) for item in pair)
            for pair in value["result_laplace_points_ri"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("cq_generation", "")).strip()
    operator_digest = str(value.get("transfer_operator_sha256", "")).lower()
    result_digest = str(value.get("time_history_sha256", "")).lower()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "contour_cq_generation",
                "timestep_cq_generation",
                "laplace_branch_cq_generation",
                "transfer_operator_cq_generation",
                "inverse_transform_cq_generation",
                "result_cq_generation",
            )
        )
        and math.isfinite(radius)
        and 0.0 < radius < 1.0
        and result_radius == radius
        and math.isfinite(timestep)
        and timestep > 0.0
        and result_timestep == timestep
        and value.get("laplace_branch") == "principal_sqrt_outgoing"
        and value.get("result_laplace_branch") == value.get("laplace_branch")
        and len(points) >= 4
        and all(
            len(pair) == 2 and all(math.isfinite(item) for item in pair)
            for pair in points
        )
        and all(pair[0] > 0.0 for pair in points)
        and result_points == points
        and bool(str(value.get("transfer_operator_id", "")).strip())
        and value.get("result_transfer_operator_id")
        == value.get("transfer_operator_id")
        and _is_sha256(operator_digest)
        and str(value.get("result_transfer_operator_sha256", "")).lower()
        == operator_digest
        and value.get("inverse_transform") == "fft_conjugate_symmetric"
        and value.get("result_inverse_transform")
        == value.get("inverse_transform")
        and _is_sha256(result_digest)
        and str(value.get("reported_time_history_sha256", "")).lower()
        == result_digest
    )


def _optional_fembem_coupling_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "fembem_trace_map_normal_material_wavenumber_coupling_matrix_mesh_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        density = float(value["fluid_density_kg_m3"])
        result_density = float(value["result_fluid_density_kg_m3"])
        sound_speed = float(value["sound_speed_m_s"])
        result_sound_speed = float(value["result_sound_speed_m_s"])
        wavenumber = tuple(float(item) for item in value["wavenumber_ri_m_inv"])
        result_wavenumber = tuple(
            float(item) for item in value["result_wavenumber_ri_m_inv"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("coupling_generation", "")).strip()
    digests = (
        ("trace_map_sha256", "result_trace_map_sha256"),
        ("normal_field_sha256", "result_normal_field_sha256"),
        ("coupling_matrix_sha256", "result_coupling_matrix_sha256"),
        ("volume_mesh_sha256", "result_volume_mesh_sha256"),
        ("boundary_mesh_sha256", "result_boundary_mesh_sha256"),
        ("coupled_result_sha256", "reported_coupled_result_sha256"),
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "trace_map_coupling_generation",
                "normal_coupling_generation",
                "material_coupling_generation",
                "wavenumber_coupling_generation",
                "matrix_coupling_generation",
                "mesh_coupling_generation",
                "result_coupling_generation",
            )
        )
        and value.get("normal_orientation") == "volume_outward"
        and value.get("result_normal_orientation")
        == value.get("normal_orientation")
        and math.isfinite(density)
        and density > 0.0
        and result_density == density
        and math.isfinite(sound_speed)
        and sound_speed > 0.0
        and result_sound_speed == sound_speed
        and len(wavenumber) == 2
        and all(math.isfinite(item) for item in wavenumber)
        and wavenumber[0] > 0.0
        and wavenumber[1] >= 0.0
        and result_wavenumber == wavenumber
        and all(
            _is_sha256(str(value.get(source, "")).lower())
            and str(value.get(target, "")).lower()
            == str(value.get(source, "")).lower()
            for source, target in digests
        )
    )


def _optional_adaptive_cq_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_adaptive_contour_quadrature_order_startup_correction_error_estimator_restart_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        radii = tuple(float(item) for item in value["contour_radii"])
        result_radii = tuple(float(item) for item in value["result_contour_radii"])
        orders = tuple(_integer({"value": item}, "value") for item in value["quadrature_orders"])
        result_orders = tuple(
            _integer({"value": item}, "value")
            for item in value["result_quadrature_orders"]
        )
        tolerance = float(value["relative_tolerance"])
        result_tolerance = float(value["result_relative_tolerance"])
        errors = tuple(float(item) for item in value["estimated_relative_errors"])
        result_errors = tuple(
            float(item) for item in value["result_estimated_relative_errors"]
        )
        restart_step = _integer(value, "restart_step")
        result_restart_step = _integer(value, "result_restart_step")
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("cq_generation", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "contour_cq_generation",
                "quadrature_order_cq_generation",
                "startup_correction_cq_generation",
                "error_estimator_cq_generation",
                "restart_cq_generation",
                "result_cq_generation",
            )
        )
        and value.get("contour_family") == "lubich_bdf2_circle"
        and value.get("result_contour_family") == value.get("contour_family")
        and len(radii) >= 2
        and all(math.isfinite(item) and 0.0 < item < 1.0 for item in radii)
        and result_radii == radii
        and len(orders) == len(radii)
        and all(item > 0 for item in orders)
        and all(left < right for left, right in zip(orders, orders[1:]))
        and result_orders == orders
        and value.get("startup_correction") == "bdf2_consistent_two_step"
        and value.get("result_startup_correction") == value.get("startup_correction")
        and value.get("error_estimator") == "successive_contour_l2_relative"
        and value.get("result_error_estimator") == value.get("error_estimator")
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and result_tolerance == tolerance
        and len(errors) == len(orders)
        and all(math.isfinite(item) and item >= 0.0 for item in errors)
        and result_errors == errors
        and errors[-1] <= tolerance
        and restart_step >= 0
        and result_restart_step == restart_step
        and _is_sha256(str(value.get("restart_state_sha256", "")).lower())
        and str(value.get("loaded_restart_state_sha256", "")).lower()
        == str(value.get("restart_state_sha256", "")).lower()
        and _is_sha256(str(value.get("time_history_sha256", "")).lower())
        and str(value.get("accepted_time_history_sha256", "")).lower()
        == str(value.get("time_history_sha256", "")).lower()
    )


def _optional_p1_fembem_discretization_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "p1_fembem_boundary_orientation_quadrature_singular_treatment_trace_matrix_mesh_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        fem_order = _integer(value, "fem_basis_order")
        bem_order = _integer(value, "bem_basis_order")
        trace_shape = tuple(
            _integer({"value": item}, "value") for item in value["trace_shape"]
        )
        result_trace_shape = tuple(
            _integer({"value": item}, "value")
            for item in value["result_trace_shape"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("coupling_generation", "")).strip()
    digests = (
        ("trace_matrix_sha256", "result_trace_matrix_sha256"),
        ("fem_matrix_sha256", "result_fem_matrix_sha256"),
        ("bem_matrix_sha256", "result_bem_matrix_sha256"),
        ("volume_mesh_sha256", "result_volume_mesh_sha256"),
        ("boundary_mesh_sha256", "result_boundary_mesh_sha256"),
        ("coupled_result_sha256", "accepted_coupled_result_sha256"),
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "boundary_orientation_coupling_generation",
                "quadrature_coupling_generation",
                "singular_treatment_coupling_generation",
                "trace_coupling_generation",
                "matrix_coupling_generation",
                "mesh_coupling_generation",
                "result_coupling_generation",
            )
        )
        and fem_order == bem_order == 1
        and value.get("volume_element") == "tet"
        and value.get("boundary_element") == "tri"
        and value.get("boundary_orientation") == "volume_outward"
        and value.get("result_boundary_orientation") == value.get("boundary_orientation")
        and value.get("regular_quadrature") == "triangle_degree_4"
        and value.get("result_regular_quadrature") == value.get("regular_quadrature")
        and value.get("singular_treatment") == "duffy_p1_galerkin"
        and value.get("result_singular_treatment") == value.get("singular_treatment")
        and len(trace_shape) == 2
        and all(item > 0 for item in trace_shape)
        and trace_shape[0] <= trace_shape[1]
        and result_trace_shape == trace_shape
        and all(
            _is_sha256(str(value.get(source, "")).lower())
            and str(value.get(target, "")).lower()
            == str(value.get(source, "")).lower()
            for source, target in digests
        )
    )


def _optional_hmatrix_aca_cluster_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "hmatrix_aca_cluster_permutation_admissibility_rank_tolerance_kernel_mesh_result_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        permutation = tuple(
            _integer({"value": item}, "value")
            for item in value["cluster_permutation"]
        )
        result_permutation = tuple(
            _integer({"value": item}, "value")
            for item in value["result_cluster_permutation"]
        )
        eta = float(value["admissibility_eta"])
        result_eta = float(value["result_admissibility_eta"])
        rank = _integer(value, "aca_rank")
        result_rank = _integer(value, "result_aca_rank")
        tolerance = float(value["relative_tolerance"])
        result_tolerance = float(value["result_relative_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("hmatrix_generation", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "cluster_hmatrix_generation",
                "permutation_hmatrix_generation",
                "admissibility_hmatrix_generation",
                "rank_hmatrix_generation",
                "tolerance_hmatrix_generation",
                "kernel_hmatrix_generation",
                "mesh_hmatrix_generation",
                "result_hmatrix_generation",
            )
        )
        and len(permutation) >= 2
        and sorted(permutation) == list(range(1, len(permutation) + 1))
        and result_permutation == permutation
        and value.get("admissibility_rule") == "eta-weak"
        and value.get("result_admissibility_rule") == "eta-weak"
        and math.isfinite(eta)
        and eta > 0.0
        and result_eta == eta
        and rank > 0
        and result_rank == rank
        and math.isfinite(tolerance)
        and 0.0 < tolerance < 1.0
        and result_tolerance == tolerance
        and value.get("kernel") == "helmholtz-single-layer-p1"
        and value.get("result_kernel") == "helmholtz-single-layer-p1"
        and _is_sha256(str(value.get("cluster_tree_sha256", "")).lower())
        and value.get("loaded_cluster_tree_sha256")
        == value.get("cluster_tree_sha256")
        and _is_sha256(str(value.get("mesh_sha256", "")).lower())
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_calderon_cq_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "calderon_cq_operator_v_k_trace_normal_frequency_grid_inverse_transform_mesh_result_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        trace_shape = tuple(
            _integer({"value": item}, "value") for item in value["trace_shape"]
        )
        result_trace_shape = tuple(
            _integer({"value": item}, "value")
            for item in value["result_trace_shape"]
        )
        frequencies = tuple(
            tuple(float(component) for component in row)
            for row in value["laplace_frequency_ri"]
        )
        result_frequencies = tuple(
            tuple(float(component) for component in row)
            for row in value["result_laplace_frequency_ri"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("calderon_generation", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "v_calderon_generation",
                "k_calderon_generation",
                "trace_calderon_generation",
                "normal_calderon_generation",
                "frequency_calderon_generation",
                "inverse_calderon_generation",
                "mesh_calderon_generation",
                "result_calderon_generation",
            )
        )
        and _is_sha256(str(value.get("v_operator_sha256", "")).lower())
        and value.get("result_v_operator_sha256") == value.get("v_operator_sha256")
        and _is_sha256(str(value.get("k_operator_sha256", "")).lower())
        and value.get("result_k_operator_sha256") == value.get("k_operator_sha256")
        and value.get("trace_basis") == "p1-nodal-boundary-trace"
        and value.get("result_trace_basis") == "p1-nodal-boundary-trace"
        and len(trace_shape) == 2
        and all(item > 0 for item in trace_shape)
        and trace_shape[0] <= trace_shape[1]
        and result_trace_shape == trace_shape
        and value.get("boundary_normal") == "volume-outward"
        and value.get("result_boundary_normal") == "volume-outward"
        and len(frequencies) >= 3
        and all(
            len(pair) == 2
            and all(math.isfinite(item) for item in pair)
            and pair[0] > 0.0
            for pair in frequencies
        )
        and result_frequencies == frequencies
        and value.get("inverse_transform") == "bdf2-cq-ifft-real"
        and value.get("result_inverse_transform") == "bdf2-cq-ifft-real"
        and _is_sha256(str(value.get("boundary_mesh_sha256", "")).lower())
        and value.get("result_boundary_mesh_sha256")
        == value.get("boundary_mesh_sha256")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_bem_near_singular_quadrature_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "bem_near_singular_quadrature_distance_element_size_adaptive_order_reference_result_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        distance = float(value["target_distance_m"])
        result_distance = float(value["result_target_distance_m"])
        element_size = float(value["element_size_m"])
        result_element_size = float(value["result_element_size_m"])
        ratio = float(value["distance_size_ratio"])
        result_ratio = float(value["result_distance_size_ratio"])
        order = _integer(value, "adaptive_order")
        result_order = _integer(value, "result_adaptive_order")
        reference = tuple(float(item) for item in value["reference_integral_ri"])
        computed = tuple(float(item) for item in value["computed_integral_ri"])
        reported_error = float(value["relative_error"])
        tolerance = float(value["relative_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("quadrature_generation", "")).strip()
    actual_error = math.hypot(
        computed[0] - reference[0], computed[1] - reference[1]
    ) / max(math.hypot(*reference), 1.0e-300) if len(reference) == len(computed) == 2 else math.inf
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "target_quadrature_generation",
                "geometry_quadrature_generation",
                "order_quadrature_generation",
                "map_quadrature_generation",
                "kernel_quadrature_generation",
                "reference_quadrature_generation",
                "mesh_quadrature_generation",
                "result_quadrature_generation",
            )
        )
        and all(math.isfinite(item) and item > 0.0 for item in (distance, element_size))
        and result_distance == distance
        and result_element_size == element_size
        and math.isfinite(ratio)
        and ratio > 0.0
        and math.isclose(ratio, distance / element_size, rel_tol=1.0e-12)
        and result_ratio == ratio
        and order >= 8
        and result_order == order
        and value.get("quadrature_rule") == "adaptive-duffy-p1"
        and value.get("result_quadrature_rule") == value.get("quadrature_rule")
        and value.get("coordinate_map") == "target-aligned-barycentric"
        and value.get("result_coordinate_map") == value.get("coordinate_map")
        and value.get("kernel") == "helmholtz-single-layer-p1"
        and value.get("result_kernel") == value.get("kernel")
        and len(reference) == len(computed) == 2
        and all(math.isfinite(item) for item in reference + computed)
        and math.isfinite(reported_error)
        and math.isclose(reported_error, actual_error, rel_tol=0.25, abs_tol=1.0e-15)
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and reported_error <= tolerance
        and _is_sha256(str(value.get("element_mesh_sha256", "")).lower())
        and value.get("result_element_mesh_sha256") == value.get("element_mesh_sha256")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_fembem_energy_reciprocity_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "fembem_energy_flux_reciprocity_interface_trace_orientation_frequency_incident_result_generation_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        trace_shape = tuple(
            _integer({"value": item}, "value")
            for item in value["interface_trace_shape"]
        )
        result_trace_shape = tuple(
            _integer({"value": item}, "value")
            for item in value["result_interface_trace_shape"]
        )
        frequency = float(value["frequency_hz"])
        result_frequency = float(value["result_frequency_hz"])
        pair_ids = tuple(str(item) for item in value["reciprocity_pair_ids"])
        result_pair_ids = tuple(
            str(item) for item in value["result_reciprocity_pair_ids"]
        )
        reciprocity = tuple(
            tuple(float(component) for component in row)
            for row in value["reciprocity_values_ri"]
        )
        result_reciprocity = tuple(
            tuple(float(component) for component in row)
            for row in value["result_reciprocity_values_ri"]
        )
        reciprocity_error = float(value["reciprocity_relative_error"])
        reciprocity_tolerance = float(value["reciprocity_relative_tolerance"])
        fem_power = float(value["fem_outward_power_w"])
        bem_power = float(value["bem_radiated_power_w"])
        energy_error = float(value["energy_flux_relative_error"])
        energy_tolerance = float(value["energy_flux_relative_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("coupling_generation", "")).strip()
    actual_energy_error = abs(fem_power - bem_power) / max(
        abs(fem_power), 1.0e-300
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "trace_coupling_generation",
                "normal_coupling_generation",
                "frequency_coupling_generation",
                "incident_coupling_generation",
                "reciprocity_coupling_generation",
                "energy_coupling_generation",
                "result_coupling_generation",
            )
        )
        and value.get("interface_trace_basis") == "p1-nodal-boundary-trace"
        and value.get("result_interface_trace_basis")
        == value.get("interface_trace_basis")
        and len(trace_shape) == 2
        and all(item > 0 for item in trace_shape)
        and trace_shape[0] <= trace_shape[1]
        and result_trace_shape == trace_shape
        and value.get("normal_orientation") == "volume-outward"
        and value.get("result_normal_orientation")
        == value.get("normal_orientation")
        and math.isfinite(frequency)
        and frequency > 0.0
        and result_frequency == frequency
        and _is_sha256(str(value.get("incident_field_sha256", "")).lower())
        and value.get("result_incident_field_sha256")
        == value.get("incident_field_sha256")
        and len(pair_ids)
        == len(result_pair_ids)
        == len(reciprocity)
        == len(result_reciprocity)
        == 2
        and len(set(pair_ids)) == 2
        and result_pair_ids == pair_ids
        and all(
            len(row) == 2 and all(math.isfinite(item) for item in row)
            for row in reciprocity
        )
        and result_reciprocity == reciprocity
        and math.isfinite(reciprocity_error)
        and 0.0 <= reciprocity_error <= reciprocity_tolerance
        and math.isfinite(reciprocity_tolerance)
        and reciprocity_tolerance > 0.0
        and math.isfinite(fem_power)
        and fem_power > 0.0
        and math.isfinite(bem_power)
        and bem_power >= 0.0
        and math.isfinite(energy_error)
        and math.isclose(
            energy_error, actual_energy_error, rel_tol=1.0e-6, abs_tol=1.0e-15
        )
        and math.isfinite(energy_tolerance)
        and energy_tolerance > 0.0
        and energy_error <= energy_tolerance
        and _is_sha256(str(value.get("coupled_result_sha256", "")).lower())
        and value.get("accepted_coupled_result_sha256")
        == value.get("coupled_result_sha256")
    )


def _optional_hmatrix_recompression_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "hmatrix_recompression_svd_tolerance_norm_rank_permutation_operator_mesh_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        tolerance = float(value["tolerance"])
        result_tolerance = float(value["result_tolerance"])
        ranks_before = tuple(
            _integer({"value": item}, "value")
            for item in value["block_ranks_before"]
        )
        result_ranks_before = tuple(
            _integer({"value": item}, "value")
            for item in value["result_block_ranks_before"]
        )
        ranks_after = tuple(
            _integer({"value": item}, "value")
            for item in value["block_ranks_after"]
        )
        result_ranks_after = tuple(
            _integer({"value": item}, "value")
            for item in value["result_block_ranks_after"]
        )
        row_permutation = tuple(
            _integer({"value": item}, "value")
            for item in value["row_permutation"]
        )
        result_row_permutation = tuple(
            _integer({"value": item}, "value")
            for item in value["result_row_permutation"]
        )
        column_permutation = tuple(
            _integer({"value": item}, "value")
            for item in value["column_permutation"]
        )
        result_column_permutation = tuple(
            _integer({"value": item}, "value")
            for item in value["result_column_permutation"]
        )
        operator_error = float(value["operator_relative_error"])
        result_operator_error = float(value["result_operator_relative_error"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("hmatrix_generation", "")).strip()
    valid_row_permutation = sorted(row_permutation) == list(
        range(len(row_permutation))
    )
    valid_column_permutation = sorted(column_permutation) == list(
        range(len(column_permutation))
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "svd_hmatrix_generation",
                "tolerance_hmatrix_generation",
                "rank_hmatrix_generation",
                "permutation_hmatrix_generation",
                "operator_hmatrix_generation",
                "mesh_hmatrix_generation",
                "result_hmatrix_generation",
            )
        )
        and value.get("svd_basis") == "euclidean-orthonormal"
        and value.get("result_svd_basis") == value.get("svd_basis")
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and result_tolerance == tolerance
        and value.get("tolerance_norm") == "spectral-relative"
        and value.get("result_tolerance_norm") == value.get("tolerance_norm")
        and bool(ranks_before)
        and len(ranks_before) == len(ranks_after)
        and all(before > 0 and 0 < after <= before for before, after in zip(ranks_before, ranks_after))
        and result_ranks_before == ranks_before
        and result_ranks_after == ranks_after
        and bool(row_permutation)
        and valid_row_permutation
        and result_row_permutation == row_permutation
        and bool(column_permutation)
        and valid_column_permutation
        and result_column_permutation == column_permutation
        and math.isfinite(operator_error)
        and 0.0 <= operator_error <= tolerance
        and result_operator_error == operator_error
        and _is_sha256(str(value.get("mesh_sha256", "")).lower())
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_cq_block_restart_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_restart_block_history_startup_weight_time_index_sample_contour_operator_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        block_size = _positive_integer(value, "block_size")
        result_block_size = _positive_integer(value, "result_block_size")
        block_ids = tuple(
            _integer({"value": item}, "value")
            for item in value["completed_block_ids"]
        )
        result_block_ids = tuple(
            _integer({"value": item}, "value")
            for item in value["result_completed_block_ids"]
        )
        history_count = _positive_integer(value, "history_sample_count")
        result_history_count = _positive_integer(
            value, "result_history_sample_count"
        )
        startup_weights = tuple(
            tuple(float(component) for component in row)
            for row in value["startup_weights_ri"]
        )
        result_startup_weights = tuple(
            tuple(float(component) for component in row)
            for row in value["result_startup_weights_ri"]
        )
        restart_index = _integer(value, "restart_time_index")
        result_restart_index = _integer(value, "result_restart_time_index")
        total_count = _positive_integer(value, "total_sample_count")
        result_total_count = _positive_integer(value, "result_total_sample_count")
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("cq_generation", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "block_cq_generation",
                "history_cq_generation",
                "startup_cq_generation",
                "time_cq_generation",
                "sample_cq_generation",
                "owner_cq_generation",
                "result_cq_generation",
            )
        )
        and result_block_size == block_size
        and bool(block_ids)
        and block_ids == tuple(range(len(block_ids)))
        and result_block_ids == block_ids
        and history_count == block_size * len(block_ids)
        and result_history_count == history_count
        and bool(startup_weights)
        and all(
            len(row) == 2 and all(math.isfinite(component) for component in row)
            for row in startup_weights
        )
        and result_startup_weights == startup_weights
        and restart_index == history_count
        and result_restart_index == restart_index
        and total_count > restart_index
        and result_total_count == total_count
        and _is_sha256(str(value.get("contour_owner_sha256", "")).lower())
        and value.get("result_contour_owner_sha256")
        == value.get("contour_owner_sha256")
        and _is_sha256(str(value.get("operator_owner_sha256", "")).lower())
        and value.get("result_operator_owner_sha256")
        == value.get("operator_owner_sha256")
        and _is_sha256(str(value.get("history_sha256", "")).lower())
        and value.get("loaded_history_sha256") == value.get("history_sha256")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_complex_ad_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "complex_ad_wirtinger_conjugation_branch_scaling_fd_mesh_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        scaling = tuple(float(item) for item in value["design_variable_scaling"])
        result_scaling = tuple(
            float(item) for item in value["result_design_variable_scaling"]
        )
        gradient = tuple(
            tuple(float(component) for component in row)
            for row in value["gradient_ri"]
        )
        fd_gradient = tuple(
            tuple(float(component) for component in row)
            for row in value["finite_difference_gradient_ri"]
        )
        fd_error = float(value["finite_difference_relative_error"])
        fd_tolerance = float(value["finite_difference_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("ad_generation", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "wirtinger_ad_generation",
                "conjugation_ad_generation",
                "branch_ad_generation",
                "scaling_ad_generation",
                "finite_difference_ad_generation",
                "mesh_ad_generation",
                "result_ad_generation",
            )
        )
        and value.get("wirtinger_convention") == "dJ_dconj_z"
        and value.get("result_wirtinger_convention")
        == value.get("wirtinger_convention")
        and value.get("adjoint_conjugation") == "conjugate_transpose"
        and value.get("result_adjoint_conjugation")
        == value.get("adjoint_conjugation")
        and value.get("objective_branch") == "real_objective"
        and value.get("result_objective_branch") == value.get("objective_branch")
        and bool(scaling)
        and all(math.isfinite(item) and item > 0.0 for item in scaling)
        and result_scaling == scaling
        and len(gradient) == len(fd_gradient) == len(scaling)
        and all(
            len(row) == 2 and all(math.isfinite(component) for component in row)
            for row in gradient + fd_gradient
        )
        and math.isfinite(fd_error)
        and fd_error >= 0.0
        and math.isfinite(fd_tolerance)
        and fd_tolerance > 0.0
        and fd_error <= fd_tolerance
        and _is_sha256(str(value.get("mesh_sha256", "")).lower())
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_pde_quadratic_vol_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "pde_quadratic_curved_vol_midnode_tet_boundary_region_order_mesh_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        order = _integer(value, "geometry_order")
        result_order = _integer(value, "result_geometry_order")
        tets = tuple(
            tuple(_positive_integer({"value": node}, "value") for node in row)
            for row in value["tet_connectivity"]
        )
        result_tets = tuple(
            tuple(_positive_integer({"value": node}, "value") for node in row)
            for row in value["result_tet_connectivity"]
        )
        triangles = tuple(
            tuple(_positive_integer({"value": node}, "value") for node in row)
            for row in value["boundary_tri_connectivity"]
        )
        result_triangles = tuple(
            tuple(_positive_integer({"value": node}, "value") for node in row)
            for row in value["result_boundary_tri_connectivity"]
        )
        orientations = tuple(
            _integer({"value": item}, "value")
            for item in value["boundary_orientation"]
        )
        result_orientations = tuple(
            _integer({"value": item}, "value")
            for item in value["result_boundary_orientation"]
        )
        tet_regions = tuple(
            _positive_integer({"value": item}, "value")
            for item in value["tet_region_labels"]
        )
        result_tet_regions = tuple(
            _positive_integer({"value": item}, "value")
            for item in value["result_tet_region_labels"]
        )
        boundary_regions = tuple(
            _positive_integer({"value": item}, "value")
            for item in value["boundary_region_labels"]
        )
        result_boundary_regions = tuple(
            _positive_integer({"value": item}, "value")
            for item in value["result_boundary_region_labels"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("mesh_generation", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "midnode_mesh_generation",
                "tet_mesh_generation",
                "boundary_mesh_generation",
                "region_mesh_generation",
                "order_mesh_generation",
                "result_mesh_generation",
            )
        )
        and order == 2
        and result_order == order
        and bool(tets)
        and all(len(row) == 10 and len(set(row)) == 10 for row in tets)
        and result_tets == tets
        and _is_sha256(str(value.get("curved_midnode_sha256", "")).lower())
        and value.get("result_curved_midnode_sha256")
        == value.get("curved_midnode_sha256")
        and bool(triangles)
        and all(len(row) == 6 and len(set(row)) == 6 for row in triangles)
        and result_triangles == triangles
        and len(orientations) == len(triangles)
        and all(item in {-1, 1} for item in orientations)
        and result_orientations == orientations
        and len(tet_regions) == len(tets)
        and result_tet_regions == tet_regions
        and len(boundary_regions) == len(triangles)
        and result_boundary_regions == boundary_regions
        and _is_sha256(str(value.get("mesh_sha256", "")).lower())
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
    )


def _optional_adaptive_cq_restart_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "adaptive_cq_timestep_contour_history_interpolation_error_restart_operator_mesh_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        timesteps = tuple(float(item) for item in value["timestep_schedule_s"])
        result_timesteps = tuple(
            float(item) for item in value["result_timestep_schedule_s"]
        )
        rebuild_indices = tuple(
            _integer({"value": item}, "value")
            for item in value["contour_rebuild_indices"]
        )
        result_rebuild_indices = tuple(
            _integer({"value": item}, "value")
            for item in value["result_contour_rebuild_indices"]
        )
        errors = tuple(float(item) for item in value["local_error_estimates"])
        result_errors = tuple(
            float(item) for item in value["result_local_error_estimates"]
        )
        tolerance = float(value["local_error_tolerance"])
        result_tolerance = float(value["result_local_error_tolerance"])
        restart_index = _integer(value, "restart_index")
        result_restart_index = _integer(value, "result_restart_index")
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("cq_generation", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
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
        )
        and len(timesteps) >= 3
        and all(math.isfinite(item) and item > 0.0 for item in timesteps)
        and result_timesteps == timesteps
        and rebuild_indices
        and rebuild_indices[0] == 0
        and tuple(sorted(set(rebuild_indices))) == rebuild_indices
        and rebuild_indices[-1] < len(timesteps)
        and result_rebuild_indices == rebuild_indices
        and all(
            timesteps[index] == timesteps[index - 1]
            or index in rebuild_indices
            for index in range(1, len(timesteps))
        )
        and value.get("history_interpolation") == "barycentric_causal"
        and value.get("result_history_interpolation")
        == value.get("history_interpolation")
        and len(errors) == len(timesteps)
        and all(math.isfinite(item) and item >= 0.0 for item in errors)
        and result_errors == errors
        and math.isfinite(tolerance)
        and tolerance > 0.0
        and max(errors) <= tolerance
        and result_tolerance == tolerance
        and restart_index in rebuild_indices
        and result_restart_index == restart_index
        and _is_sha256(str(value.get("operator_owner_sha256", "")).lower())
        and value.get("result_operator_owner_sha256")
        == value.get("operator_owner_sha256")
        and _is_sha256(str(value.get("history_sha256", "")).lower())
        and value.get("loaded_history_sha256") == value.get("history_sha256")
        and _is_sha256(str(value.get("mesh_sha256", "")).lower())
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_modal_fembem_transient_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "fembem_modal_transient_mass_damping_initial_projection_truncation_energy_mesh_history_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        rayleigh = tuple(float(item) for item in value["rayleigh_coefficients"])
        result_rayleigh = tuple(
            float(item) for item in value["result_rayleigh_coefficients"]
        )
        displacement = tuple(
            float(item) for item in value["initial_displacement_projection"]
        )
        result_displacement = tuple(
            float(item) for item in value["result_initial_displacement_projection"]
        )
        velocity = tuple(float(item) for item in value["initial_velocity_projection"])
        result_velocity = tuple(
            float(item) for item in value["result_initial_velocity_projection"]
        )
        modal_count = _positive_integer(value, "modal_count")
        result_modal_count = _positive_integer(value, "result_modal_count")
        truncation_frequency = float(value["truncation_frequency_hz"])
        result_truncation_frequency = float(value["result_truncation_frequency_hz"])
        initial_energy = float(value["initial_energy_j"])
        result_initial_energy = float(value["result_initial_energy_j"])
        radiated_energy = float(value["radiated_energy_j"])
        result_radiated_energy = float(value["result_radiated_energy_j"])
        dissipated_energy = float(value["dissipated_energy_j"])
        result_dissipated_energy = float(value["result_dissipated_energy_j"])
        final_energy = float(value["final_energy_j"])
        result_final_energy = float(value["result_final_energy_j"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("modal_generation", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
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
        )
        and value.get("mass_normalization") == "M_orthonormal"
        and value.get("result_mass_normalization") == "M_orthonormal"
        and value.get("damping_model") == "rayleigh"
        and value.get("result_damping_model") == "rayleigh"
        and len(rayleigh) == 2
        and all(math.isfinite(item) and item >= 0.0 for item in rayleigh)
        and result_rayleigh == rayleigh
        and len(displacement) == len(velocity) == modal_count
        and all(math.isfinite(item) for item in displacement + velocity)
        and result_displacement == displacement
        and result_velocity == velocity
        and result_modal_count == modal_count
        and math.isfinite(truncation_frequency)
        and truncation_frequency > 0.0
        and result_truncation_frequency == truncation_frequency
        and all(
            math.isfinite(item) and item >= 0.0
            for item in (initial_energy, radiated_energy, dissipated_energy, final_energy)
        )
        and math.isclose(
            radiated_energy + dissipated_energy + final_energy,
            initial_energy,
            rel_tol=1.0e-10,
            abs_tol=1.0e-12,
        )
        and result_initial_energy == initial_energy
        and result_radiated_energy == radiated_energy
        and result_dissipated_energy == dissipated_energy
        and result_final_energy == final_energy
        and _is_sha256(str(value.get("mesh_sha256", "")).lower())
        and value.get("result_mesh_sha256") == value.get("mesh_sha256")
        and bool(str(value.get("time_history_owner", "")).strip())
        and value.get("accepted_time_history_owner") == value.get("time_history_owner")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_calderon_projector_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "calderon_projector_p1_v_k_kt_w_mass_duality_normal_quadrature_mesh_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        signs = tuple(_integer(value, key) for key in ("v_sign", "k_sign", "kt_sign", "w_sign"))
        result_signs = tuple(
            _integer(value, key)
            for key in ("result_v_sign", "result_k_sign", "result_kt_sign", "result_w_sign")
        )
        mass_residual = float(value["mass_duality_residual"])
        result_mass_residual = float(value["result_mass_duality_residual"])
        mass_tolerance = float(value["mass_duality_tolerance"])
        result_mass_tolerance = float(value["result_mass_duality_tolerance"])
        projector_residual = float(value["projector_residual"])
        result_projector_residual = float(value["result_projector_residual"])
        projector_tolerance = float(value["projector_tolerance"])
        result_projector_tolerance = float(value["result_projector_tolerance"])
        block_order = tuple(str(item) for item in value["block_order"])
        result_block_order = tuple(str(item) for item in value["result_block_order"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("calderon_generation", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
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
        )
        and value.get("trial_space") == "P1"
        and value.get("result_trial_space") == "P1"
        and value.get("test_space") == "P1"
        and value.get("result_test_space") == "P1"
        and value.get("projector_convention") == "interior_calderon_outward"
        and value.get("result_projector_convention")
        == value.get("projector_convention")
        and block_order == ("dirichlet", "neumann")
        and result_block_order == block_order
        and signs == (-1, 1, -1, -1)
        and result_signs == signs
        and all(
            math.isfinite(item) and item >= 0.0
            for item in (
                mass_residual,
                result_mass_residual,
                projector_residual,
                result_projector_residual,
            )
        )
        and all(
            math.isfinite(item) and item > 0.0
            for item in (
                mass_tolerance,
                result_mass_tolerance,
                projector_tolerance,
                result_projector_tolerance,
            )
        )
        and mass_residual <= mass_tolerance
        and result_mass_residual == mass_residual
        and result_mass_tolerance == mass_tolerance
        and value.get("normal_orientation") == "outward"
        and value.get("result_normal_orientation") == "outward"
        and value.get("singular_quadrature") == "duffy_principal_value_p1"
        and value.get("result_singular_quadrature")
        == value.get("singular_quadrature")
        and projector_residual <= projector_tolerance
        and result_projector_residual == projector_residual
        and result_projector_tolerance == projector_tolerance
        and _is_sha256(str(value.get("operator_sha256", "")).lower())
        and value.get("result_operator_sha256") == value.get("operator_sha256")
        and _is_sha256(str(value.get("mass_sha256", "")).lower())
        and value.get("result_mass_sha256") == value.get("mass_sha256")
        and _is_sha256(str(value.get("boundary_mesh_sha256", "")).lower())
        and value.get("result_boundary_mesh_sha256")
        == value.get("boundary_mesh_sha256")
        and bool(str(value.get("result_owner", "")).strip())
        and value.get("accepted_result_owner") == value.get("result_owner")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_cq_physical_closure_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_symbol_contour_transfer_conjugate_causal_ifft_parseval_passivity_timestep_operator_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        coefficients = tuple(float(item) for item in value["symbol_coefficients"])
        result_coefficients = tuple(
            float(item) for item in value["result_symbol_coefficients"]
        )
        radius = float(value["contour_radius"])
        result_radius = float(value["result_contour_radius"])
        samples = tuple(
            (float(row[0]), float(row[1])) for row in value["transfer_samples_ri"]
        )
        result_samples = tuple(
            (float(row[0]), float(row[1]))
            for row in value["result_transfer_samples_ri"]
        )
        response = tuple(float(item) for item in value["time_response"])
        result_response = tuple(float(item) for item in value["result_time_response"])
        negative_energy = float(value["negative_time_energy"])
        result_negative_energy = float(value["result_negative_time_energy"])
        time_work = float(value["time_domain_work"])
        result_time_work = float(value["result_time_domain_work"])
        frequency_work = float(value["frequency_domain_work"])
        result_frequency_work = float(value["result_frequency_domain_work"])
        parseval_tolerance = float(value["parseval_tolerance"])
        result_parseval_tolerance = float(value["result_parseval_tolerance"])
        minimum_real = float(value["minimum_real_transfer"])
        result_minimum_real = float(value["result_minimum_real_transfer"])
        timestep = float(value["timestep_s"])
        result_timestep = float(value["result_timestep_s"])
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    generation = str(value.get("cq_generation", "")).strip()
    conjugate_pairs_close = bool(samples) and math.isclose(
        samples[0][1], 0.0, rel_tol=0.0, abs_tol=1.0e-12
    )
    for index in range(1, (len(samples) + 1) // 2):
        conjugate_pairs_close = conjugate_pairs_close and math.isclose(
            samples[index][0], samples[-index][0], rel_tol=1.0e-12, abs_tol=1.0e-12
        ) and math.isclose(
            samples[index][1], -samples[-index][1], rel_tol=1.0e-12, abs_tol=1.0e-12
        )
    expected_minimum_real = min((item[0] for item in samples), default=math.nan)
    return (
        bool(generation)
        and all(
            value.get(key) == generation
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
        )
        and value.get("multistep_symbol") == "BDF2"
        and value.get("result_multistep_symbol") == "BDF2"
        and coefficients == (1.5, -2.0, 0.5)
        and result_coefficients == coefficients
        and math.isfinite(radius)
        and 0.0 < radius < 1.0
        and result_radius == radius
        and len(samples) >= 5
        and all(math.isfinite(part) for sample in samples for part in sample)
        and result_samples == samples
        and conjugate_pairs_close
        and bool(response)
        and all(math.isfinite(item) for item in response)
        and result_response == response
        and math.isfinite(negative_energy)
        and negative_energy >= 0.0
        and negative_energy <= 1.0e-12
        and result_negative_energy == negative_energy
        and all(
            math.isfinite(item) and item >= 0.0
            for item in (time_work, result_time_work, frequency_work, result_frequency_work)
        )
        and math.isfinite(parseval_tolerance)
        and parseval_tolerance > 0.0
        and result_parseval_tolerance == parseval_tolerance
        and math.isclose(
            time_work,
            frequency_work,
            rel_tol=parseval_tolerance,
            abs_tol=parseval_tolerance,
        )
        and result_time_work == time_work
        and result_frequency_work == frequency_work
        and math.isfinite(minimum_real)
        and minimum_real >= 0.0
        and math.isclose(
            minimum_real,
            expected_minimum_real,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        and result_minimum_real == minimum_real
        and value.get("passivity_sign") == "nonnegative_real_transfer"
        and value.get("result_passivity_sign") == value.get("passivity_sign")
        and math.isfinite(timestep)
        and timestep > 0.0
        and result_timestep == timestep
        and value.get("operator_family") == "p1_calderon_bem"
        and value.get("result_operator_family") == value.get("operator_family")
        and _is_sha256(str(value.get("operator_sha256", "")).lower())
        and value.get("result_operator_sha256") == value.get("operator_sha256")
        and bool(str(value.get("result_owner", "")).strip())
        and value.get("accepted_result_owner") == value.get("result_owner")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_fembem_reciprocity_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "fembem_reciprocity_radiation_power_interior_energy_trace_orientation_boundary_volume_map_frequency_mesh_solution_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        frequency = float(value["frequency_hz"])
        result_frequency = float(value["result_frequency_hz"])
        transfer_ab = complex(*[float(item) for item in value["transfer_ab_ri"]])
        result_transfer_ab = complex(
            *[float(item) for item in value["result_transfer_ab_ri"]]
        )
        transfer_ba = complex(*[float(item) for item in value["transfer_ba_ri"]])
        result_transfer_ba = complex(
            *[float(item) for item in value["result_transfer_ba_ri"]]
        )
        reciprocity_tolerance = float(value["reciprocity_tolerance"])
        result_reciprocity_tolerance = float(value["result_reciprocity_tolerance"])
        radiated_power = float(value["radiated_power_w"])
        result_radiated_power = float(value["result_radiated_power_w"])
        flux_power = float(value["boundary_flux_power_w"])
        result_flux_power = float(value["result_boundary_flux_power_w"])
        interior_energy = float(value["interior_energy_j"])
        result_interior_energy = float(value["result_interior_energy_j"])
        node_map = tuple(_integer({"item": item}, "item") for item in value["boundary_volume_node_map"])
        result_node_map = tuple(
            _integer({"item": item}, "item")
            for item in value["result_boundary_volume_node_map"]
        )
        trace_nodes = tuple(_integer({"item": item}, "item") for item in value["trace_node_ids"])
        result_trace_nodes = tuple(
            _integer({"item": item}, "item") for item in value["result_trace_node_ids"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("fembem_generation", "")).strip()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "transfer_fembem_generation",
                "radiation_fembem_generation",
                "interior_fembem_generation",
                "trace_fembem_generation",
                "map_fembem_generation",
                "frequency_fembem_generation",
                "mesh_fembem_generation",
                "solution_fembem_generation",
                "result_fembem_generation",
            )
        )
        and math.isfinite(frequency)
        and frequency > 0.0
        and result_frequency == frequency
        and math.isfinite(reciprocity_tolerance)
        and 0.0 < reciprocity_tolerance <= 1.0e-6
        and result_reciprocity_tolerance == reciprocity_tolerance
        and abs(transfer_ab - transfer_ba) <= reciprocity_tolerance
        and result_transfer_ab == transfer_ab
        and result_transfer_ba == transfer_ba
        and all(
            math.isfinite(item) and item >= 0.0
            for item in (
                radiated_power,
                result_radiated_power,
                flux_power,
                result_flux_power,
                interior_energy,
                result_interior_energy,
            )
        )
        and math.isclose(
            radiated_power, flux_power, rel_tol=1.0e-10, abs_tol=1.0e-12
        )
        and result_radiated_power == radiated_power
        and result_flux_power == flux_power
        and result_interior_energy == interior_energy
        and value.get("trace_orientation") == "outward_volume_to_boundary"
        and value.get("result_trace_orientation") == value.get("trace_orientation")
        and bool(node_map)
        and all(item > 0 for item in node_map)
        and len(set(node_map)) == len(node_map)
        and trace_nodes == node_map
        and result_node_map == node_map
        and result_trace_nodes == trace_nodes
        and bool(str(value.get("mesh_owner", "")).strip())
        and value.get("accepted_mesh_owner") == value.get("mesh_owner")
        and _is_sha256(str(value.get("mesh_sha256", "")).lower())
        and value.get("accepted_mesh_sha256") == value.get("mesh_sha256")
        and bool(str(value.get("solution_owner", "")).strip())
        and value.get("accepted_solution_owner") == value.get("solution_owner")
        and _is_sha256(str(value.get("solution_sha256", "")).lower())
        and value.get("accepted_solution_sha256") == value.get("solution_sha256")
    )


def _optional_nonlinear_eigen_contour_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "nonlinear_eigen_contour_orientation_quadrature_moment_rank_count_residual_biorthogonality_pole_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        contour = tuple(
            complex(float(row[0]), float(row[1])) for row in value["contour_points_ri"]
        )
        result_contour = tuple(
            complex(float(row[0]), float(row[1]))
            for row in value["result_contour_points_ri"]
        )
        moment_ranks = tuple(_integer({"item": item}, "item") for item in value["moment_ranks"])
        result_moment_ranks = tuple(
            _integer({"item": item}, "item") for item in value["result_moment_ranks"]
        )
        rank = _integer(value, "numerical_rank")
        result_rank = _integer(value, "result_numerical_rank")
        count = _integer(value, "enclosed_eigenvalue_count")
        result_count = _integer(value, "result_enclosed_eigenvalue_count")
        eigenvalues = tuple(
            complex(float(row[0]), float(row[1])) for row in value["eigenvalues_ri"]
        )
        result_eigenvalues = tuple(
            complex(float(row[0]), float(row[1]))
            for row in value["result_eigenvalues_ri"]
        )
        residuals = tuple(float(item) for item in value["residual_norms"])
        result_residuals = tuple(float(item) for item in value["result_residual_norms"])
        gram = tuple(
            tuple(complex(float(item[0]), float(item[1])) for item in row)
            for row in value["biorthogonality_gram_ri"]
        )
        result_gram = tuple(
            tuple(complex(float(item[0]), float(item[1])) for item in row)
            for row in value["result_biorthogonality_gram_ri"]
        )
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    generation = str(value.get("nonlinear_eigen_generation", "")).strip()
    signed_area = 0.5 * sum(
        left.real * right.imag - right.real * left.imag
        for left, right in zip(contour, contour[1:] + contour[:1])
    )
    gram_ok = (
        len(gram) == rank
        and all(len(row) == rank for row in gram)
        and all(
            abs(gram[i][j] - (1.0 if i == j else 0.0)) <= 1.0e-8
            for i in range(rank)
            for j in range(rank)
        )
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "contour_eigen_generation",
                "quadrature_eigen_generation",
                "moment_eigen_generation",
                "rank_eigen_generation",
                "count_eigen_generation",
                "residual_eigen_generation",
                "biorthogonality_eigen_generation",
                "pole_eigen_generation",
                "result_eigen_generation",
            )
        )
        and value.get("contour_orientation") == "counterclockwise"
        and value.get("result_contour_orientation") == "counterclockwise"
        and len(contour) >= 4
        and result_contour == contour
        and signed_area > 0.0
        and value.get("quadrature_rule") == "trapezoidal_periodic"
        and value.get("result_quadrature_rule") == value.get("quadrature_rule")
        and bool(moment_ranks)
        and all(item == rank for item in moment_ranks)
        and result_moment_ranks == moment_ranks
        and rank > 0
        and result_rank == rank
        and count == rank
        and result_count == count
        and len(eigenvalues) == count
        and result_eigenvalues == eigenvalues
        and len(residuals) == count
        and all(math.isfinite(item) and 0.0 <= item <= 1.0e-6 for item in residuals)
        and result_residuals == residuals
        and gram_ok
        and result_gram == gram
        and bool(str(value.get("pole_owner", "")).strip())
        and value.get("accepted_pole_owner") == value.get("pole_owner")
        and _is_sha256(str(value.get("pole_sha256", "")).lower())
        and value.get("accepted_pole_sha256") == value.get("pole_sha256")
        and bool(str(value.get("result_owner", "")).strip())
        and value.get("accepted_result_owner") == value.get("result_owner")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_cq_acoustic_identity_is_aligned(summary: Mapping[str, Any]) -> bool:
    value = summary.get(
        "cq_acoustic_laplace_contour_weight_passivity_trace_timestep_history_mesh_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        radius = float(value["contour_radius"])
        result_radius = float(value["result_contour_radius"])
        zeta = tuple(complex(float(row[0]), float(row[1])) for row in value["zeta_points_ri"])
        result_zeta = tuple(
            complex(float(row[0]), float(row[1]))
            for row in value["result_zeta_points_ri"]
        )
        laplace = tuple(
            complex(float(row[0]), float(row[1]))
            for row in value["laplace_points_ri"]
        )
        result_laplace = tuple(
            complex(float(row[0]), float(row[1]))
            for row in value["result_laplace_points_ri"]
        )
        weights = tuple(float(item) for item in value["cq_weights"])
        result_weights = tuple(float(item) for item in value["result_cq_weights"])
        impedances = tuple(
            complex(float(row[0]), float(row[1]))
            for row in value["boundary_impedance_ri"]
        )
        result_impedances = tuple(
            complex(float(row[0]), float(row[1]))
            for row in value["result_boundary_impedance_ri"]
        )
        fem_nodes = tuple(_integer({"item": item}, "item") for item in value["fem_trace_node_ids"])
        result_fem_nodes = tuple(
            _integer({"item": item}, "item") for item in value["result_fem_trace_node_ids"]
        )
        bem_nodes = tuple(_integer({"item": item}, "item") for item in value["bem_trace_node_ids"])
        result_bem_nodes = tuple(
            _integer({"item": item}, "item") for item in value["result_bem_trace_node_ids"]
        )
        trace_sign = _integer(value, "trace_sign")
        result_trace_sign = _integer(value, "result_trace_sign")
        timestep = float(value["time_step_s"])
        result_timestep = float(value["result_time_step_s"])
        history_length = _integer(value, "history_length")
        result_history_length = _integer(value, "result_history_length")
        times = tuple(float(item) for item in value["time_samples_s"])
        result_times = tuple(float(item) for item in value["result_time_samples_s"])
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    generation = str(value.get("cq_generation", "")).strip()
    expected_laplace = tuple(
        (1.5 - 2.0 * point + 0.5 * point * point) / timestep
        for point in zeta
    ) if math.isfinite(timestep) and timestep > 0.0 else ()
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "contour_cq_generation", "weight_cq_generation",
                "passivity_cq_generation", "trace_cq_generation",
                "timestep_cq_generation", "history_cq_generation",
                "mesh_cq_generation", "owner_cq_generation",
                "result_cq_generation",
            )
        )
        and value.get("cq_method") == "bdf2"
        and value.get("result_cq_method") == value.get("cq_method")
        and math.isfinite(radius) and 0.0 < radius < 1.0
        and result_radius == radius
        and len(zeta) >= 2 and result_zeta == zeta
        and all(math.isclose(abs(point), radius, rel_tol=1.0e-12, abs_tol=1.0e-12) for point in zeta)
        and len(laplace) == len(zeta) and result_laplace == laplace
        and all(abs(actual - expected) <= 1.0e-10 * max(abs(expected), 1.0) for actual, expected in zip(laplace, expected_laplace))
        and len(weights) == history_length and result_weights == weights
        and all(math.isfinite(item) and item >= 0.0 for item in weights)
        and math.isclose(sum(weights), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and len(impedances) == len(zeta) and result_impedances == impedances
        and all(math.isfinite(item.real) and math.isfinite(item.imag) and item.real >= 0.0 for item in impedances)
        and value.get("trace_orientation") == "outward_volume_to_boundary"
        and value.get("result_trace_orientation") == value.get("trace_orientation")
        and bool(fem_nodes) and fem_nodes == bem_nodes
        and all(item > 0 for item in fem_nodes) and len(set(fem_nodes)) == len(fem_nodes)
        and result_fem_nodes == fem_nodes and result_bem_nodes == bem_nodes
        and trace_sign == 1 and result_trace_sign == trace_sign
        and math.isfinite(timestep) and timestep > 0.0 and result_timestep == timestep
        and history_length >= 2 and result_history_length == history_length
        and len(times) == history_length and result_times == times
        and all(math.isfinite(item) for item in times)
        and all(math.isclose(item, index * timestep, rel_tol=1.0e-12, abs_tol=1.0e-15) for index, item in enumerate(times))
        and bool(str(value.get("mesh_owner", "")).strip())
        and value.get("accepted_mesh_owner") == value.get("mesh_owner")
        and _is_sha256(str(value.get("mesh_sha256", "")).lower())
        and value.get("accepted_mesh_sha256") == value.get("mesh_sha256")
        and bool(str(value.get("result_owner", "")).strip())
        and value.get("accepted_result_owner") == value.get("result_owner")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_fembem_autodiff_identity_is_aligned(summary: Mapping[str, Any]) -> bool:
    value = summary.get(
        "fembem_autodiff_wirtinger_objective_shape_fd_trace_mesh_owner_gradient_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        design = tuple(float(item) for item in value["complex_design_ri"])
        result_design = tuple(float(item) for item in value["result_complex_design_ri"])
        objective = float(value["objective_value"])
        result_objective = float(value["result_objective_value"])
        wirtinger = tuple(float(item) for item in value["wirtinger_gradient_ri"])
        result_wirtinger = tuple(float(item) for item in value["result_wirtinger_gradient_ri"])
        real_gradient = tuple(float(item) for item in value["real_gradient_ri"])
        result_real_gradient = tuple(float(item) for item in value["result_real_gradient_ri"])
        direction = tuple(float(item) for item in value["shape_direction_ri"])
        result_direction = tuple(float(item) for item in value["result_shape_direction_ri"])
        step = float(value["shape_step"])
        result_step = float(value["result_shape_step"])
        finite_difference = float(value["finite_difference_directional_derivative"])
        result_finite_difference = float(value["result_finite_difference_directional_derivative"])
        trace_map = tuple(_integer({"item": item}, "item") for item in value["trace_node_map"])
        result_trace_map = tuple(
            _integer({"item": item}, "item") for item in value["result_trace_node_map"]
        )
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("autodiff_generation", "")).strip()
    finite_vectors = all(
        math.isfinite(item)
        for vector in (design, wirtinger, real_gradient, direction)
        for item in vector
    )
    return (
        bool(generation)
        and all(
            value.get(key) == generation
            for key in (
                "complex_autodiff_generation", "wirtinger_autodiff_generation",
                "objective_autodiff_generation", "shape_autodiff_generation",
                "finite_difference_autodiff_generation", "trace_autodiff_generation",
                "mesh_autodiff_generation", "owner_autodiff_generation",
                "gradient_autodiff_generation", "result_autodiff_generation",
            )
        )
        and len(design) == 2 and result_design == design and finite_vectors
        and value.get("objective_scaling") == "one_half_l2_squared"
        and value.get("result_objective_scaling") == value.get("objective_scaling")
        and math.isfinite(objective)
        and math.isclose(objective, 0.5 * sum(item * item for item in design), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and result_objective == objective
        and value.get("wirtinger_convention") == "dJ_dconjugate_z"
        and value.get("result_wirtinger_convention") == value.get("wirtinger_convention")
        and len(wirtinger) == len(design)
        and all(math.isclose(actual, 0.5 * expected, rel_tol=1.0e-12, abs_tol=1.0e-12) for actual, expected in zip(wirtinger, design))
        and result_wirtinger == wirtinger
        and real_gradient == design and result_real_gradient == real_gradient
        and len(direction) == len(design) and result_direction == direction
        and math.isclose(sum(item * item for item in direction), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isfinite(step) and step > 0.0 and result_step == step
        and math.isfinite(finite_difference)
        and math.isclose(finite_difference, sum(gradient * tangent for gradient, tangent in zip(real_gradient, direction)), rel_tol=1.0e-10, abs_tol=1.0e-12)
        and result_finite_difference == finite_difference
        and bool(trace_map) and all(item > 0 for item in trace_map)
        and len(set(trace_map)) == len(trace_map) and result_trace_map == trace_map
        and bool(str(value.get("mesh_owner", "")).strip())
        and value.get("accepted_mesh_owner") == value.get("mesh_owner")
        and _is_sha256(str(value.get("mesh_sha256", "")).lower())
        and value.get("accepted_mesh_sha256") == value.get("mesh_sha256")
        and bool(str(value.get("gradient_owner", "")).strip())
        and value.get("accepted_gradient_owner") == value.get("gradient_owner")
        and _is_sha256(str(value.get("gradient_sha256", "")).lower())
        and value.get("accepted_gradient_sha256") == value.get("gradient_sha256")
    )


def _optional_cq_adaptive_timestep_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_adaptive_timestep_contour_restart_interpolation_causality_energy_operator_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["adaptive_cq_generation"]).strip()
        time_steps = tuple(float(item) for item in value["time_step_history_s"])
        result_time_steps = tuple(
            float(item) for item in value["result_time_step_history_s"]
        )
        time_samples = tuple(float(item) for item in value["time_samples_s"])
        result_time_samples = tuple(
            float(item) for item in value["result_time_samples_s"]
        )
        anchors = tuple(float(item) for item in value["laplace_anchor_real_per_s"])
        result_anchors = tuple(
            float(item) for item in value["result_laplace_anchor_real_per_s"]
        )
        radius = float(value["contour_radius"])
        result_radius = float(value["result_contour_radius"])
        restart_step = _integer(value, "restart_step")
        result_restart_step = _integer(value, "result_restart_step")
        restart_state = tuple(float(item) for item in value["restart_history_state"])
        result_restart_state = tuple(
            float(item) for item in value["result_restart_history_state"]
        )
        prehistory = float(value["prehistory_max_abs"])
        result_prehistory = float(value["result_prehistory_max_abs"])
        energy = tuple(float(item) for item in value["discrete_energy_j"])
        result_energy = tuple(float(item) for item in value["result_discrete_energy_j"])
    except (KeyError, TypeError, ValueError):
        return False
    expected_times = [0.0]
    for time_step in time_steps:
        expected_times.append(expected_times[-1] + time_step)
    anchor_products = tuple(anchor * step for anchor, step in zip(anchors, time_steps))
    return (
        bool(generation)
        and all(
            value.get(key) == generation
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
        )
        and value.get("cq_method") == "bdf2"
        and value.get("result_cq_method") == value.get("cq_method")
        and bool(time_steps)
        and all(math.isfinite(item) and item > 0.0 for item in time_steps)
        and result_time_steps == time_steps
        and len(time_samples) == len(time_steps) + 1
        and all(math.isfinite(item) for item in time_samples)
        and all(
            math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-15)
            for actual, expected in zip(time_samples, expected_times)
        )
        and result_time_samples == time_samples
        and math.isfinite(radius)
        and 0.0 < radius < 1.0
        and result_radius == radius
        and len(anchors) == len(time_steps)
        and all(math.isfinite(item) and item > 0.0 for item in anchors)
        and result_anchors == anchors
        and bool(anchor_products)
        and all(
            math.isclose(item, anchor_products[0], rel_tol=1.0e-12, abs_tol=1.0e-15)
            for item in anchor_products
        )
        and 0 < restart_step <= len(time_steps)
        and result_restart_step == restart_step
        and bool(restart_state)
        and all(math.isfinite(item) for item in restart_state)
        and result_restart_state == restart_state
        and value.get("history_interpolation") == "piecewise_linear_causal"
        and value.get("result_history_interpolation")
        == value.get("history_interpolation")
        and math.isfinite(prehistory)
        and prehistory == 0.0
        and result_prehistory == prehistory
        and len(energy) == len(time_samples)
        and all(math.isfinite(item) and item >= 0.0 for item in energy)
        and all(right <= left for left, right in zip(energy, energy[1:]))
        and result_energy == energy
        and bool(str(value.get("operator_owner", "")).strip())
        and value.get("accepted_operator_owner") == value.get("operator_owner")
        and _is_sha256(str(value.get("operator_sha256", "")).lower())
        and value.get("accepted_operator_sha256") == value.get("operator_sha256")
        and _is_sha256(str(value.get("result_sha256", "")).lower())
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def _optional_hmatrix_benchmark_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "hmatrix_dense_reference_error_tolerance_rank_memory_complexity_mesh_operator_benchmark_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["hmatrix_generation"]).strip()
        sizes = tuple(_integer({"item": item}, "item") for item in value["boundary_unknown_count"])
        errors = tuple(float(item) for item in value["dense_reference_relative_error"])
        tolerance = float(value["relative_tolerance"])
        ranks = tuple(_integer({"item": item}, "item") for item in value["maximum_block_rank"])
        dense_memory = tuple(_integer({"item": item}, "item") for item in value["dense_memory_bytes"])
        hmatrix_memory = tuple(_integer({"item": item}, "item") for item in value["hmatrix_memory_bytes"])
        memory_exponent = float(value["memory_complexity_exponent"])
        rank_exponent = float(value["rank_complexity_exponent"])
    except (KeyError, TypeError, ValueError):
        return False
    empirical_memory_exponent = (
        math.log(hmatrix_memory[-1] / hmatrix_memory[0])
        / math.log(sizes[-1] / sizes[0])
        if len(sizes) >= 2 and min(sizes[0], hmatrix_memory[0]) > 0 else math.nan
    )
    empirical_rank_exponent = (
        math.log(ranks[-1] / ranks[0]) / math.log(sizes[-1] / sizes[0])
        if len(sizes) >= 2 and min(sizes[0], ranks[0]) > 0 else math.nan
    )
    mirrored = (
        "boundary_unknown_count", "dense_reference_relative_error",
        "relative_tolerance", "maximum_block_rank", "dense_memory_bytes",
        "hmatrix_memory_bytes", "memory_complexity_exponent",
        "rank_complexity_exponent", "boundary_mesh_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "dense_generation", "tolerance_generation", "rank_generation",
            "memory_generation", "complexity_generation", "mesh_generation",
            "operator_generation", "benchmark_generation", "result_generation",
        ))
        and len(sizes) == len(errors) == len(ranks) == len(dense_memory) == len(hmatrix_memory) >= 3
        and all(item > 0 for item in sizes)
        and all(left < right for left, right in zip(sizes, sizes[1:]))
        and math.isfinite(tolerance) and 0.0 < tolerance < 1.0
        and all(math.isfinite(item) and 0.0 <= item <= tolerance for item in errors)
        and all(right <= left for left, right in zip(errors, errors[1:]))
        and all(0 < rank < size for rank, size in zip(ranks, sizes))
        and all(left <= right for left, right in zip(ranks, ranks[1:]))
        and dense_memory == tuple(16 * size * size for size in sizes)
        and all(0 < compressed < dense for compressed, dense in zip(hmatrix_memory, dense_memory))
        and all(left < right for left, right in zip(hmatrix_memory, hmatrix_memory[1:]))
        and math.isfinite(empirical_memory_exponent) and 1.0 <= empirical_memory_exponent < 2.0
        and math.isclose(memory_exponent, empirical_memory_exponent, rel_tol=0.02, abs_tol=0.02)
        and math.isfinite(empirical_rank_exponent) and 0.0 <= empirical_rank_exponent <= 1.0
        and math.isclose(rank_exponent, empirical_rank_exponent, rel_tol=0.02, abs_tol=0.02)
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _is_sha256(str(value.get("boundary_mesh_sha256", "")).lower())
        and bool(str(value.get("operator_owner", "")).strip())
        and value.get("accepted_operator_owner") == value.get("operator_owner")
        and bool(str(value.get("benchmark_owner", "")).strip())
        and value.get("accepted_benchmark_owner") == value.get("benchmark_owner")
        and _is_sha256(str(value.get("hmatrix_result_sha256", "")).lower())
        and value.get("accepted_hmatrix_result_sha256") == value.get("hmatrix_result_sha256")
    )


def _optional_multifrequency_adjoint_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "multifrequency_fembem_adjoint_weight_objective_quadrature_trace_gradient_fd_mesh_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["adjoint_generation"]).strip()
        frequencies = tuple(float(item) for item in value["frequency_hz"])
        weights = tuple(float(item) for item in value["frequency_weights"])
        objectives = tuple(complex(float(row[0]), float(row[1])) for row in value["objective_complex"])
        weighted_objective = complex(*[float(item) for item in value["weighted_objective_complex"]])
        quadrature = tuple(_integer({"item": item}, "item") for item in value["quadrature_order"])
        trace_map = tuple(_integer({"item": item}, "item") for item in value["trace_node_map"])
        gradients = tuple(float(item) for item in value["frequency_gradient"])
        accumulated = float(value["accumulated_gradient"])
        finite_difference = float(value["finite_difference_gradient"])
        tolerance = float(value["gradient_relative_tolerance"])
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    count = len(frequencies)
    expected_objective = sum(weight * objective for weight, objective in zip(weights, objectives))
    expected_gradient = sum(weight * gradient for weight, gradient in zip(weights, gradients))
    mirrored = (
        "frequency_hz", "frequency_weights", "objective_complex",
        "weighted_objective_complex", "quadrature_order", "trace_node_map",
        "frequency_gradient", "accumulated_gradient", "finite_difference_gradient",
        "gradient_relative_tolerance", "fembem_mesh_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "frequency_generation", "weight_generation", "objective_generation",
            "quadrature_generation", "trace_generation", "gradient_generation",
            "fd_generation", "mesh_generation", "owner_generation", "result_generation",
        ))
        and count == len(weights) == len(objectives) == len(quadrature) == len(gradients) >= 3
        and all(math.isfinite(item) and item > 0.0 for item in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and all(math.isfinite(item) and item >= 0.0 for item in weights)
        and math.isclose(sum(weights), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and all(math.isfinite(item.real) and math.isfinite(item.imag) for item in objectives)
        and abs(weighted_objective - expected_objective) <= 1.0e-12
        and all(item >= 2 and item % 2 == 0 for item in quadrature)
        and bool(trace_map) and all(item > 0 for item in trace_map)
        and len(set(trace_map)) == len(trace_map)
        and all(math.isfinite(item) for item in gradients)
        and math.isclose(accumulated, expected_gradient, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isfinite(tolerance) and 0.0 < tolerance <= 1.0e-2
        and math.isfinite(finite_difference)
        and abs(accumulated - finite_difference) / max(abs(accumulated), abs(finite_difference), 1.0e-300) <= tolerance
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _is_sha256(str(value.get("fembem_mesh_sha256", "")).lower())
        and bool(str(value.get("adjoint_owner", "")).strip())
        and value.get("accepted_adjoint_owner") == value.get("adjoint_owner")
        and _is_sha256(str(value.get("adjoint_result_sha256", "")).lower())
        and value.get("accepted_adjoint_result_sha256") == value.get("adjoint_result_sha256")
    )


def _optional_fembem_shape_derivative_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "fembem_shape_derivative_morph_normal_velocity_trace_jacobian_objective_fd_mesh_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["shape_generation"]).strip()
        step = float(value["shape_step"])
        result_step = float(value["result_shape_step"])
        reference = tuple(tuple(float(item) for item in row) for row in value["reference_nodes_m"])
        result_reference = tuple(
            tuple(float(item) for item in row) for row in value["result_reference_nodes_m"]
        )
        velocity = tuple(float(item) for item in value["normal_velocity_m"])
        result_velocity = tuple(float(item) for item in value["result_normal_velocity_m"])
        morphed = tuple(tuple(float(item) for item in row) for row in value["morphed_nodes_m"])
        result_morphed = tuple(
            tuple(float(item) for item in row) for row in value["result_morphed_nodes_m"]
        )
        trace_map = tuple(_integer({"item": item}, "item") for item in value["trace_node_map"])
        result_trace_map = tuple(
            _integer({"item": item}, "item") for item in value["result_trace_node_map"]
        )
        jacobians = tuple(float(item) for item in value["geometry_jacobian_determinant"])
        result_jacobians = tuple(
            float(item) for item in value["result_geometry_jacobian_determinant"]
        )
        derivative = float(value["objective_directional_derivative"])
        result_derivative = float(value["result_objective_directional_derivative"])
        objective_minus = float(value["objective_minus"])
        result_objective_minus = float(value["result_objective_minus"])
        objective_plus = float(value["objective_plus"])
        result_objective_plus = float(value["result_objective_plus"])
    except (KeyError, TypeError, ValueError):
        return False
    expected_morphed = tuple(
        (node[0], node[1], node[2] + step * normal_velocity)
        for node, normal_velocity in zip(reference, velocity)
    )
    central_difference = (objective_plus - objective_minus) / (2.0 * step) if step else math.nan
    return (
        bool(generation)
        and all(
            value.get(key) == generation
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
        )
        and math.isfinite(step)
        and step > 0.0
        and result_step == step
        and bool(reference)
        and all(len(row) == 3 and all(math.isfinite(item) for item in row) for row in reference)
        and result_reference == reference
        and len(velocity) == len(reference)
        and all(math.isfinite(item) for item in velocity)
        and result_velocity == velocity
        and len(morphed) == len(reference)
        and all(len(row) == 3 and all(math.isfinite(item) for item in row) for row in morphed)
        and all(
            all(math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-15) for actual, expected in zip(actual_row, expected_row))
            for actual_row, expected_row in zip(morphed, expected_morphed)
        )
        and result_morphed == morphed
        and len(trace_map) == len(reference)
        and all(item > 0 for item in trace_map)
        and len(set(trace_map)) == len(trace_map)
        and result_trace_map == trace_map
        and len(jacobians) == len(reference)
        and all(math.isfinite(item) and item > 0.0 for item in jacobians)
        and result_jacobians == jacobians
        and math.isfinite(derivative)
        and result_derivative == derivative
        and all(math.isfinite(item) for item in (objective_minus, objective_plus))
        and result_objective_minus == objective_minus
        and result_objective_plus == objective_plus
        and math.isclose(central_difference, derivative, rel_tol=1.0e-10, abs_tol=1.0e-12)
        and bool(str(value.get("mesh_owner", "")).strip())
        and value.get("accepted_mesh_owner") == value.get("mesh_owner")
        and _is_sha256(str(value.get("mesh_sha256", "")).lower())
        and value.get("accepted_mesh_sha256") == value.get("mesh_sha256")
        and _is_sha256(str(value.get("shape_result_sha256", "")).lower())
        and value.get("accepted_shape_result_sha256") == value.get("shape_result_sha256")
    )


def _optional_simp_topology_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "simp_topology_density_filter_projection_volume_compliance_adjoint_fd_kkt_mesh_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["topology_generation"]).strip()
        density = tuple(float(item) for item in value["design_density"])
        filter_matrix = tuple(
            tuple(float(item) for item in row)
            for row in value["density_filter_matrix"]
        )
        filtered = tuple(float(item) for item in value["filtered_density"])
        beta = tuple(float(item) for item in value["projection_beta_continuation"])
        eta = float(value["projection_eta"])
        projected = tuple(float(item) for item in value["projected_density"])
        volume_fraction = float(value["volume_fraction"])
        volume_limit = float(value["volume_fraction_limit"])
        compliance = float(value["compliance"])
        adjoint = tuple(float(item) for item in value["adjoint_compliance_gradient"])
        finite_difference = tuple(
            float(item) for item in value["finite_difference_compliance_gradient"]
        )
        volume_gradient = tuple(float(item) for item in value["volume_gradient"])
        multiplier = float(value["volume_lagrange_multiplier"])
        kkt_residual = float(value["kkt_stationarity_residual"])
        gradient_tolerance = float(value["gradient_relative_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    count = len(density)
    if count == 0 or len(filter_matrix) != count:
        return False
    if any(len(row) != count for row in filter_matrix):
        return False
    expected_filtered = tuple(
        sum(weight * item for weight, item in zip(row, density))
        for row in filter_matrix
    )
    if not beta:
        return False
    denominator = math.tanh(beta[-1] * eta) + math.tanh(beta[-1] * (1.0 - eta))
    if not math.isfinite(denominator) or denominator <= 0.0:
        return False
    expected_projected = tuple(
        (
            math.tanh(beta[-1] * eta)
            + math.tanh(beta[-1] * (item - eta))
        )
        / denominator
        for item in filtered
    )
    expected_volume = sum(projected) / count
    expected_kkt = max(
        abs(gradient + multiplier * constraint_gradient)
        for gradient, constraint_gradient in zip(adjoint, volume_gradient)
    )
    mirrored = (
        "design_density", "density_filter_matrix", "filtered_density",
        "projection_beta_continuation", "projection_eta", "projected_density",
        "volume_fraction", "volume_fraction_limit", "compliance",
        "adjoint_compliance_gradient", "finite_difference_compliance_gradient",
        "volume_gradient", "volume_lagrange_multiplier",
        "kkt_stationarity_residual", "gradient_relative_tolerance",
        "topology_mesh_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "density_generation", "filter_generation", "projection_generation",
            "volume_generation", "compliance_generation", "adjoint_generation",
            "fd_generation", "kkt_generation", "mesh_generation",
            "owner_generation", "result_generation",
        ))
        and all(math.isfinite(item) and 0.0 <= item <= 1.0 for item in density)
        and all(
            math.isfinite(weight) and weight >= 0.0
            for row in filter_matrix for weight in row
        )
        and all(math.isclose(sum(row), 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12) for row in filter_matrix)
        and len(filtered) == len(projected) == count
        and all(math.isclose(left, right, rel_tol=1.0e-12, abs_tol=1.0e-12) for left, right in zip(filtered, expected_filtered))
        and all(math.isfinite(item) and item > 0.0 for item in beta)
        and all(left < right for left, right in zip(beta, beta[1:]))
        and math.isfinite(eta) and 0.0 < eta < 1.0
        and all(math.isclose(left, right, rel_tol=1.0e-12, abs_tol=1.0e-12) for left, right in zip(projected, expected_projected))
        and all(0.0 <= item <= 1.0 for item in projected)
        and math.isclose(volume_fraction, expected_volume, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isfinite(volume_limit) and 0.0 < volume_limit <= 1.0
        and volume_fraction <= volume_limit + max(1.0e-2, gradient_tolerance)
        and math.isfinite(compliance) and compliance > 0.0
        and len(adjoint) == len(finite_difference) == len(volume_gradient) == count
        and all(math.isfinite(item) for item in (*adjoint, *finite_difference, *volume_gradient))
        and math.isfinite(gradient_tolerance) and 0.0 < gradient_tolerance <= 1.0e-2
        and all(
            abs(left - right) / max(abs(left), abs(right), 1.0e-300) <= gradient_tolerance
            for left, right in zip(adjoint, finite_difference)
        )
        and math.isfinite(multiplier) and multiplier >= 0.0
        and math.isclose(kkt_residual, expected_kkt, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and kkt_residual <= gradient_tolerance
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _is_sha256(str(value.get("topology_mesh_sha256", "")).lower())
        and bool(str(value.get("topology_owner", "")).strip())
        and value.get("accepted_topology_owner") == value.get("topology_owner")
        and _is_sha256(str(value.get("topology_result_sha256", "")).lower())
        and value.get("accepted_topology_result_sha256") == value.get("topology_result_sha256")
    )


def _optional_fembem_model_reduction_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "fembem_model_reduction_projection_order_stability_passivity_moment_frequency_error_full_mesh_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["reduction_generation"]).strip()
        full_order = _integer(value, "full_order")
        reduced_order = _integer(value, "reduced_order")
        trial = tuple(tuple(float(item) for item in row) for row in value["trial_projection_basis"])
        test = tuple(tuple(float(item) for item in row) for row in value["test_projection_basis"])
        gram = tuple(tuple(float(item) for item in row) for row in value["biorthogonality_gram"])
        full_poles = tuple(complex(float(row[0]), float(row[1])) for row in value["full_model_poles"])
        reduced_poles = tuple(complex(float(row[0]), float(row[1])) for row in value["reduced_model_poles"])
        passivity = float(value["minimum_passivity_eigenvalue"])
        full_moments = tuple(complex(float(row[0]), float(row[1])) for row in value["matched_moments_full"])
        reduced_moments = tuple(complex(float(row[0]), float(row[1])) for row in value["matched_moments_reduced"])
        frequencies = tuple(float(item) for item in value["frequency_hz"])
        errors = tuple(float(item) for item in value["frequency_response_relative_error"])
        maximum_error = float(value["maximum_frequency_response_relative_error"])
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    if not (
        0 < reduced_order < full_order
        and len(trial) == len(test) == full_order
        and all(len(row) == reduced_order for row in (*trial, *test))
        and len(gram) == reduced_order
        and all(len(row) == reduced_order for row in gram)
    ):
        return False
    expected_gram = tuple(
        tuple(sum(test[row][left] * trial[row][right] for row in range(full_order)) for right in range(reduced_order))
        for left in range(reduced_order)
    )
    mirrored = (
        "full_order", "reduced_order", "trial_projection_basis",
        "test_projection_basis", "biorthogonality_gram", "full_model_poles",
        "reduced_model_poles", "minimum_passivity_eigenvalue",
        "matched_moments_full", "matched_moments_reduced", "frequency_hz",
        "frequency_response_relative_error",
        "maximum_frequency_response_relative_error", "reduction_mesh_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "projection_generation", "order_generation", "stability_generation",
            "passivity_generation", "moment_generation", "frequency_generation",
            "error_generation", "full_model_generation", "mesh_generation",
            "owner_generation", "result_generation",
        ))
        and all(math.isfinite(item) for row in (*trial, *test, *gram) for item in row)
        and all(math.isclose(gram[i][j], expected_gram[i][j], rel_tol=1.0e-12, abs_tol=1.0e-12) for i in range(reduced_order) for j in range(reduced_order))
        and all(math.isclose(gram[i][j], 1.0 if i == j else 0.0, rel_tol=1.0e-12, abs_tol=1.0e-12) for i in range(reduced_order) for j in range(reduced_order))
        and len(full_poles) == full_order and len(reduced_poles) == reduced_order
        and all(math.isfinite(item.real) and math.isfinite(item.imag) and item.real < 0.0 for item in (*full_poles, *reduced_poles))
        and math.isfinite(passivity) and passivity >= 0.0
        and len(full_moments) == len(reduced_moments) >= reduced_order
        and all(math.isfinite(item.real) and math.isfinite(item.imag) for item in (*full_moments, *reduced_moments))
        and all(abs(left - right) <= 1.0e-12 for left, right in zip(full_moments, reduced_moments))
        and len(frequencies) == len(errors) >= 3
        and all(math.isfinite(item) and item > 0.0 for item in frequencies)
        and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and all(math.isfinite(item) and item >= 0.0 for item in errors)
        and all(right <= left for left, right in zip(errors, errors[1:]))
        and math.isfinite(maximum_error) and 0.0 < maximum_error < 1.0
        and max(errors) <= maximum_error
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _is_sha256(str(value.get("reduction_mesh_sha256", "")).lower())
        and bool(str(value.get("full_model_owner", "")).strip())
        and value.get("accepted_full_model_owner") == value.get("full_model_owner")
        and bool(str(value.get("reduction_owner", "")).strip())
        and value.get("accepted_reduction_owner") == value.get("reduction_owner")
        and _is_sha256(str(value.get("reduction_result_sha256", "")).lower())
        and value.get("accepted_reduction_result_sha256") == value.get("reduction_result_sha256")
    )


def _is_symmetric_positive_definite_matrix(
    matrix: tuple[tuple[float, ...], ...],
) -> bool:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        return False
    if not all(math.isfinite(item) for row in matrix for item in row):
        return False
    scale = max(1.0, max(abs(item) for row in matrix for item in row))
    tolerance = 1.0e-12 * scale
    if any(
        abs(matrix[row][column] - matrix[column][row]) > tolerance
        for row in range(size)
        for column in range(size)
    ):
        return False
    lower = [[0.0] * size for _ in range(size)]
    for row in range(size):
        for column in range(row + 1):
            remainder = matrix[row][column] - sum(
                lower[row][index] * lower[column][index]
                for index in range(column)
            )
            if row == column:
                if remainder <= tolerance:
                    return False
                lower[row][column] = math.sqrt(remainder)
            else:
                lower[row][column] = remainder / lower[column][column]
    return True


def _optional_nonlinear_fem_newton_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "nonlinear_fem_newton_residual_consistent_tangent_linesearch_step_energy_mesh_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["nonlinear_generation"]).strip()
        residuals = tuple(float(item) for item in value["residual_norm_history"])
        tangent = tuple(
            tuple(float(item) for item in row)
            for row in value["consistent_tangent_matrix"]
        )
        tangent_product = tuple(
            float(item) for item in value["directional_tangent_product"]
        )
        finite_difference = tuple(
            float(item) for item in value["finite_difference_directional_derivative"]
        )
        tangent_tolerance = float(value["tangent_relative_tolerance"])
        steps = tuple(
            tuple(float(item) for item in row)
            for row in value["newton_step_history"]
        )
        alphas = tuple(float(item) for item in value["line_search_alpha_history"])
        trial_residuals = tuple(
            float(item) for item in value["line_search_trial_residual_norm"]
        )
        armijo = float(value["line_search_armijo_constant"])
        energies = tuple(float(item) for item in value["strain_energy_history_j"])
        external_work = float(value["external_work_final_j"])
        energy_residual = float(value["energy_balance_residual_j"])
        convergence_tolerance = float(value["convergence_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    size = len(tangent)
    if size == 0 or any(len(row) != size for row in tangent):
        return False
    positive_definite = _is_symmetric_positive_definite_matrix(tangent)
    mirrored = (
        "nonlinear_formulation", "residual_norm_history",
        "consistent_tangent_matrix", "directional_tangent_product",
        "finite_difference_directional_derivative", "tangent_relative_tolerance",
        "newton_step_history", "line_search_alpha_history",
        "line_search_trial_residual_norm", "line_search_armijo_constant",
        "strain_energy_history_j", "external_work_final_j",
        "energy_balance_residual_j", "convergence_tolerance",
        "nonlinear_mesh_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "residual_generation", "tangent_generation", "step_generation",
            "linesearch_generation", "iteration_generation", "energy_generation",
            "mesh_generation", "owner_generation", "result_generation",
        ))
        and value.get("nonlinear_formulation") == "total_lagrangian_hyperelastic"
        and len(residuals) >= 3
        and all(math.isfinite(item) and item >= 0.0 for item in residuals)
        and all(left > right for left, right in zip(residuals, residuals[1:]))
        and math.isfinite(convergence_tolerance) and convergence_tolerance > 0.0
        and residuals[-1] <= convergence_tolerance
        and positive_definite
        and len(tangent_product) == len(finite_difference) == size
        and math.isfinite(tangent_tolerance) and 0.0 < tangent_tolerance <= 1.0e-4
        and all(
            abs(left - right) / max(abs(left), abs(right), 1.0e-300)
            <= tangent_tolerance
            for left, right in zip(tangent_product, finite_difference)
        )
        and len(steps) == len(alphas) == len(trial_residuals) == len(residuals) - 1
        and all(len(row) == size and all(math.isfinite(item) for item in row) for row in steps)
        and all(math.isfinite(item) and 0.0 < item <= 1.0 for item in alphas)
        and math.isfinite(armijo) and 0.0 < armijo < 1.0
        and all(
            math.isfinite(trial) and trial <= previous * (1.0 - armijo * alpha)
            for previous, trial, alpha in zip(residuals, trial_residuals, alphas)
        )
        and trial_residuals == residuals[1:]
        and len(energies) == len(residuals)
        and all(math.isfinite(item) and item >= 0.0 for item in energies)
        and all(left <= right for left, right in zip(energies, energies[1:]))
        and math.isfinite(external_work) and external_work >= 0.0
        and math.isclose(external_work, energies[-1], rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isfinite(energy_residual) and abs(energy_residual) <= 1.0e-12
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _is_sha256(str(value.get("nonlinear_mesh_sha256", "")).lower())
        and str(value.get("nonlinear_owner", "")).startswith("fem/")
        and value.get("accepted_nonlinear_owner") == value.get("nonlinear_owner")
        and _is_sha256(str(value.get("nonlinear_result_sha256", "")).lower())
        and value.get("accepted_nonlinear_result_sha256")
        == value.get("nonlinear_result_sha256")
    )


def _optional_cq_contour_reconstruction_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_contour_frequency_interpolation_aliasing_passivity_reconstruction_time_operator_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["cq_generation"]).strip()
        timestep = float(value["time_step_s"])
        count = _integer(value, "time_step_count")
        radius = float(value["contour_radius"])
        contour = tuple(
            complex(float(row[0]), float(row[1]))
            for row in value["contour_nodes_complex"]
        )
        interpolation_error = float(value["frequency_interpolation_relative_error"])
        interpolation_limit = float(value["maximum_frequency_interpolation_relative_error"])
        aliasing_error = float(value["aliasing_error_bound"])
        aliasing_limit = float(value["maximum_aliasing_error"])
        passivity = float(value["minimum_transfer_passivity_eigenvalue"])
        reconstruction_error = float(value["time_reconstruction_relative_error"])
        reconstruction_limit = float(value["maximum_time_reconstruction_relative_error"])
        history = tuple(float(item) for item in value["time_history"])
        reconstructed = tuple(float(item) for item in value["reconstructed_time_history"])
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    expected_contour = tuple(
        radius * complex(
            math.cos(2.0 * math.pi * index / count),
            math.sin(2.0 * math.pi * index / count),
        )
        for index in range(count)
    ) if count > 0 else ()
    mirrored = (
        "cq_method", "time_step_s", "time_step_count", "contour_radius",
        "contour_nodes_complex", "frequency_interpolation_relative_error",
        "maximum_frequency_interpolation_relative_error", "aliasing_error_bound",
        "maximum_aliasing_error", "minimum_transfer_passivity_eigenvalue",
        "time_reconstruction_relative_error",
        "maximum_time_reconstruction_relative_error", "time_history",
        "reconstructed_time_history", "cq_operator_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "contour_generation", "frequency_generation",
            "interpolation_generation", "aliasing_generation",
            "passivity_generation", "reconstruction_generation",
            "time_generation", "operator_generation", "result_generation",
        ))
        and value.get("cq_method") == "bdf2"
        and math.isfinite(timestep) and timestep > 0.0
        and count >= 4 and 0.0 < radius < 1.0 and len(contour) == count
        and all(
            abs(actual - expected) <= 1.0e-12
            for actual, expected in zip(contour, expected_contour)
        )
        and all(math.isfinite(item) and item >= 0.0 for item in (
            interpolation_error, interpolation_limit, aliasing_error,
            aliasing_limit, passivity, reconstruction_error, reconstruction_limit,
        ))
        and interpolation_error <= interpolation_limit
        and aliasing_error <= aliasing_limit
        and reconstruction_error <= reconstruction_limit
        and len(history) == len(reconstructed) == count
        and all(math.isfinite(item) for item in (*history, *reconstructed))
        and history == reconstructed
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _is_sha256(str(value.get("cq_operator_sha256", "")).lower())
        and str(value.get("operator_owner", "")).startswith("cq/")
        and value.get("accepted_operator_owner") == value.get("operator_owner")
        and _is_sha256(str(value.get("cq_result_sha256", "")).lower())
        and value.get("accepted_cq_result_sha256") == value.get("cq_result_sha256")
    )


def _optional_johnson_nedelec_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "johnson_nedelec_volume_trace_normal_single_double_layer_sign_residual_energy_mesh_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["fembem_generation"]).strip()
        trace = tuple(tuple(float(item) for item in row) for row in value["volume_trace_matrix"])
        normals = tuple(tuple(float(item) for item in row) for row in value["boundary_normals"])
        outward = tuple(tuple(float(item) for item in row) for row in value["outward_reference_vectors"])
        single_layer = tuple(tuple(float(item) for item in row) for row in value["single_layer_matrix"])
        double_layer = tuple(tuple(float(item) for item in row) for row in value["double_layer_matrix"])
        coupling_sign = float(value["coupling_sign"])
        residual = tuple(float(item) for item in value["interface_residual_vector"])
        residual_tolerance = float(value["interface_residual_tolerance"])
        interior_flux = float(value["interior_energy_flux_w"])
        exterior_flux = float(value["exterior_energy_flux_w"])
        energy_residual = float(value["energy_flux_residual_w"])
    except (KeyError, TypeError, ValueError):
        return False
    boundary_size = len(single_layer)
    trace_ok = (
        boundary_size > 0
        and len(trace) == boundary_size
        and all(len(row) > boundary_size and all(math.isfinite(item) for item in row) for row in trace)
        and all(sum(abs(item) > 1.0e-15 for item in row) == 1 for row in trace)
        and all(math.isclose(sum(row), 1.0, rel_tol=0.0, abs_tol=1.0e-12) for row in trace)
    )
    normals_ok = (
        len(normals) == len(outward) == boundary_size
        and all(len(row) == 3 and all(math.isfinite(item) for item in row) for row in normals + outward)
        and all(math.isclose(math.sqrt(sum(item * item for item in row)), 1.0, rel_tol=0.0, abs_tol=1.0e-12) for row in normals)
        and all(sum(n * reference for n, reference in zip(normal, reference_row)) > 0.0 for normal, reference_row in zip(normals, outward))
    )
    operators_ok = (
        _is_symmetric_positive_definite_matrix(single_layer)
        and len(double_layer) == boundary_size
        and all(len(row) == boundary_size and all(math.isfinite(item) for item in row) for row in double_layer)
    )
    mirrored = (
        "volume_trace_matrix", "boundary_normals", "outward_reference_vectors",
        "single_layer_matrix", "double_layer_matrix", "coupling_sign",
        "interface_residual_vector", "interface_residual_tolerance",
        "interior_energy_flux_w", "exterior_energy_flux_w",
        "energy_flux_residual_w", "fembem_mesh_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "trace_generation", "normal_generation", "single_layer_generation",
            "double_layer_generation", "coupling_generation", "residual_generation",
            "energy_generation", "mesh_generation", "owner_generation",
            "result_generation",
        ))
        and trace_ok and normals_ok and operators_ok
        and coupling_sign == -1.0
        and len(residual) == boundary_size
        and all(math.isfinite(item) for item in residual)
        and math.isfinite(residual_tolerance) and 0.0 < residual_tolerance <= 1.0e-6
        and math.sqrt(sum(item * item for item in residual)) <= residual_tolerance
        and all(math.isfinite(item) for item in (interior_flux, exterior_flux, energy_residual))
        and math.isclose(energy_residual, interior_flux + exterior_flux, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and abs(energy_residual) <= 1.0e-12
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _is_sha256(str(value.get("fembem_mesh_sha256", "")).lower())
        and str(value.get("fembem_owner", "")).startswith("acoustic/")
        and value.get("accepted_fembem_owner") == value.get("fembem_owner")
        and _is_sha256(str(value.get("fembem_result_sha256", "")).lower())
        and value.get("accepted_fembem_result_sha256") == value.get("fembem_result_sha256")
    )


def _optional_adjoint_hessian_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "adjoint_hessian_design_objective_constraint_hvp_kkt_fd_model_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["optimization_generation"]).strip()
        design = tuple(float(item) for item in value["design_variables"])
        objective_gradient = tuple(float(item) for item in value["objective_gradient"])
        jacobian = tuple(tuple(float(item) for item in row) for row in value["constraint_jacobian"])
        multipliers = tuple(float(item) for item in value["lagrange_multipliers"])
        constraints = tuple(float(item) for item in value["constraint_values"])
        direction = tuple(float(item) for item in value["hessian_vector_direction"])
        adjoint_hvp = tuple(float(item) for item in value["adjoint_hessian_vector_product"])
        finite_difference_hvp = tuple(float(item) for item in value["finite_difference_hessian_vector_product"])
        hvp_tolerance = float(value["hessian_vector_relative_tolerance"])
        kkt_residual = tuple(float(item) for item in value["kkt_stationarity_residual"])
        kkt_tolerance = float(value["kkt_residual_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    size = len(design)
    expected_stationarity = tuple(
        objective_gradient[column]
        + sum(jacobian[row][column] * multipliers[row] for row in range(len(jacobian)))
        for column in range(size)
    ) if size and len(jacobian) == len(multipliers) else ()
    hvp_scale = max(
        math.sqrt(sum(item * item for item in adjoint_hvp)),
        math.sqrt(sum(item * item for item in finite_difference_hvp)),
        1.0e-300,
    )
    hvp_error = math.sqrt(sum((left - right) ** 2 for left, right in zip(adjoint_hvp, finite_difference_hvp))) / hvp_scale
    mirrored = (
        "design_variables", "objective_gradient", "constraint_jacobian",
        "lagrange_multipliers", "constraint_values", "hessian_vector_direction",
        "adjoint_hessian_vector_product", "finite_difference_hessian_vector_product",
        "hessian_vector_relative_tolerance", "kkt_stationarity_residual",
        "kkt_residual_tolerance", "optimization_model_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "design_generation", "gradient_generation", "constraint_generation",
            "adjoint_generation", "hessian_generation", "kkt_generation",
            "finite_difference_generation", "model_generation", "owner_generation",
            "result_generation",
        ))
        and size >= 2
        and len(objective_gradient) == len(direction) == len(adjoint_hvp) == len(finite_difference_hvp) == len(kkt_residual) == size
        and all(math.isfinite(item) for item in design + objective_gradient + direction + adjoint_hvp + finite_difference_hvp + kkt_residual)
        and len(jacobian) == len(multipliers) == len(constraints) >= 1
        and all(len(row) == size and all(math.isfinite(item) for item in row) for row in jacobian)
        and all(math.isfinite(item) for item in multipliers + constraints)
        and all(abs(item) <= 1.0e-10 for item in constraints)
        and len(expected_stationarity) == size
        and all(math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-12) for actual, expected in zip(kkt_residual, expected_stationarity))
        and math.isfinite(kkt_tolerance) and 0.0 < kkt_tolerance <= 1.0e-6
        and math.sqrt(sum(item * item for item in kkt_residual)) <= kkt_tolerance
        and math.isfinite(hvp_tolerance) and 0.0 < hvp_tolerance <= 1.0e-6
        and len(adjoint_hvp) == len(finite_difference_hvp)
        and hvp_error <= hvp_tolerance
        and sum(left * right for left, right in zip(direction, adjoint_hvp)) > 0.0
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _is_sha256(str(value.get("optimization_model_sha256", "")).lower())
        and str(value.get("model_owner", "")).startswith("optimization/")
        and value.get("accepted_model_owner") == value.get("model_owner")
        and _is_sha256(str(value.get("optimization_result_sha256", "")).lower())
        and value.get("accepted_optimization_result_sha256") == value.get("optimization_result_sha256")
    )


def _optional_cq_acoustic_history_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "cq_acoustic_causality_passivity_timestep_ztransform_energy_history_mesh_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["cq_generation"]).strip()
        timestep = float(value["time_step_s"])
        radius = float(value["z_transform_radius"])
        symbols = tuple(
            complex(float(row[0]), float(row[1]))
            for row in value["multistep_symbol_samples"]
        )
        frequencies = tuple(
            complex(float(row[0]), float(row[1]))
            for row in value["laplace_frequency_samples_rad_s"]
        )
        excitation = tuple(float(item) for item in value["excitation_history"])
        pressure = tuple(float(item) for item in value["pressure_history"])
        causal_prefix = _integer(value, "causal_prefix_length")
        passivity = float(value["minimum_passivity_real_part"])
        work = float(value["boundary_work_j"])
        radiated = float(value["radiated_energy_j"])
        dissipated = float(value["dissipated_energy_j"])
        residual = float(value["energy_balance_residual_j"])
        tolerance = float(value["energy_balance_tolerance_j"])
    except (KeyError, TypeError, ValueError, IndexError):
        return False
    mirrored = (
        "multistep_method", "time_step_s", "z_transform_radius",
        "multistep_symbol_samples", "laplace_frequency_samples_rad_s",
        "excitation_history", "pressure_history", "causal_prefix_length",
        "minimum_passivity_real_part", "boundary_work_j", "radiated_energy_j",
        "dissipated_energy_j", "energy_balance_residual_j",
        "energy_balance_tolerance_j", "boundary_mesh_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "multistep_generation", "timestep_generation", "ztransform_generation",
            "frequency_generation", "history_generation", "passivity_generation",
            "energy_generation", "mesh_generation", "owner_generation",
            "result_generation",
        ))
        and value.get("multistep_method") == "bdf2"
        and math.isfinite(timestep) and timestep > 0.0
        and math.isfinite(radius) and 0.0 < radius < 1.0
        and len(symbols) == len(frequencies) >= 3
        and all(math.isfinite(item.real) and math.isfinite(item.imag) for item in symbols)
        and all(
            math.isfinite(item.real) and math.isfinite(item.imag) and item.real > 0.0
            for item in frequencies
        )
        and len(excitation) == len(pressure) >= 3
        and all(math.isfinite(item) for item in excitation + pressure)
        and 1 <= causal_prefix < len(pressure)
        and all(abs(item) <= 1.0e-12 for item in pressure[:causal_prefix])
        and math.isfinite(passivity) and passivity >= 0.0
        and all(math.isfinite(item) and item >= 0.0 for item in (
            work, radiated, dissipated, tolerance,
        ))
        and tolerance > 0.0
        and math.isfinite(residual)
        and math.isclose(
            residual, work - radiated - dissipated,
            rel_tol=1.0e-10, abs_tol=1.0e-12,
        )
        and abs(residual) <= tolerance
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _is_sha256(str(value.get("boundary_mesh_sha256", "")).lower())
        and str(value.get("cq_owner", "")).startswith("acoustic/")
        and value.get("accepted_cq_owner") == value.get("cq_owner")
        and _is_sha256(str(value.get("cq_result_sha256", "")).lower())
        and value.get("accepted_cq_result_sha256") == value.get("cq_result_sha256")
    )


def _optional_hmatrix_compression_identity_is_aligned(
    summary: Mapping[str, Any],
) -> bool:
    value = summary.get(
        "hmatrix_admissibility_cluster_rank_tolerance_matvec_error_memory_mesh_owner_result_identity"
    )
    if value is None:
        return True
    if not isinstance(value, Mapping):
        return False
    try:
        generation = str(value["hmatrix_generation"]).strip()
        leaf_size = _integer(value, "cluster_leaf_size")
        permutation = tuple(_integer({"item": item}, "item") for item in value["cluster_permutation"])
        eta = float(value["admissibility_eta"])
        partition = tuple(tuple(str(item) for item in row) for row in value["block_partition"])
        ranks = tuple(tuple(_integer({"item": item}, "item") for item in row) for row in value["numerical_ranks"])
        tolerance = float(value["compression_relative_tolerance"])
        matvec_error = float(value["measured_matvec_relative_error"])
        dense_memory = _integer(value, "dense_memory_bytes")
        compressed_memory = _integer(value, "compressed_memory_bytes")
    except (KeyError, TypeError, ValueError):
        return False
    block_count = len(partition)
    partition_ok = (
        block_count > 0
        and all(len(row) == block_count for row in partition)
        and len(ranks) == block_count
        and all(len(row) == block_count for row in ranks)
        and all(item in {"low_rank", "dense"} for row in partition for item in row)
        and all(
            (rank > 0 if kind == "low_rank" else rank == 0)
            and rank <= leaf_size
            for kinds, rank_row in zip(partition, ranks)
            for kind, rank in zip(kinds, rank_row)
        )
    )
    mirrored = (
        "cluster_leaf_size", "cluster_permutation", "admissibility_eta",
        "block_partition", "numerical_ranks", "compression_relative_tolerance",
        "measured_matvec_relative_error", "dense_memory_bytes",
        "compressed_memory_bytes", "boundary_mesh_sha256",
    )
    return (
        bool(generation)
        and all(value.get(key) == generation for key in (
            "cluster_generation", "admissibility_generation", "partition_generation",
            "rank_generation", "tolerance_generation", "matvec_generation",
            "memory_generation", "mesh_generation", "owner_generation",
            "result_generation",
        ))
        and leaf_size > 0
        and len(permutation) >= 2
        and all(item > 0 for item in permutation)
        and len(set(permutation)) == len(permutation)
        and math.isfinite(eta) and eta > 0.0
        and partition_ok
        and math.isfinite(tolerance) and 0.0 < tolerance <= 1.0e-2
        and math.isfinite(matvec_error) and 0.0 <= matvec_error <= tolerance
        and dense_memory > 0 and 0 < compressed_memory < dense_memory
        and all(value.get(f"result_{field}") == value.get(field) for field in mirrored)
        and _is_sha256(str(value.get("boundary_mesh_sha256", "")).lower())
        and str(value.get("hmatrix_owner", "")).startswith("acoustic/")
        and value.get("accepted_hmatrix_owner") == value.get("hmatrix_owner")
        and _is_sha256(str(value.get("hmatrix_result_sha256", "")).lower())
        and value.get("accepted_hmatrix_result_sha256") == value.get("hmatrix_result_sha256")
    )


def _is_sha256(value: str) -> bool:
    return len(value) == 64 and all(
        character in "0123456789abcdef" for character in value
    )


def _integer(parent: Mapping[str, Any], key: str) -> int:
    value = parent[key]
    if isinstance(value, bool) or int(value) != float(value):
        raise ValueError(f"{key} must be an integer")
    return int(value)


def _positive_integer_sequence(parent: Mapping[str, Any], key: str) -> tuple[int, ...]:
    values = parent[key]
    if not isinstance(values, (list, tuple)) or not values:
        raise ValueError(f"{key} must be a nonempty integer sequence")
    result = tuple(_integer({"value": value}, "value") for value in values)
    if any(value <= 0 for value in result):
        raise ValueError(f"{key} must contain positive integers")
    return result


def _positive_integer(parent: Mapping[str, Any], key: str) -> int:
    value = _integer(parent, key)
    if value <= 0:
        raise ValueError(f"{key} must be positive")
    return value


def _nonnegative(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and non-negative")
    return result


def _close(left: float, right: float, tol: float = 1.0e-12) -> bool:
    return abs(left - right) <= tol * max(abs(left), abs(right), 1.0)
