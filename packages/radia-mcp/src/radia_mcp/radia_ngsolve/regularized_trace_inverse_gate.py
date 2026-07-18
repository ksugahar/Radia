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
        trace_shape = tuple(_integer({"value": item}, "value") for item in value["interface_trace_shape"])
        result_trace_shape = tuple(_integer({"value": item}, "value") for item in value["result_interface_trace_shape"])
        frequency = float(value["frequency_hz"])
        result_frequency = float(value["result_frequency_hz"])
        pair_ids = tuple(str(item) for item in value["reciprocity_pair_ids"])
        result_pair_ids = tuple(str(item) for item in value["result_reciprocity_pair_ids"])
        reciprocity = tuple(tuple(float(component) for component in row) for row in value["reciprocity_values_ri"])
        result_reciprocity = tuple(tuple(float(component) for component in row) for row in value["result_reciprocity_values_ri"])
        reciprocity_error = float(value["reciprocity_relative_error"])
        reciprocity_tolerance = float(value["reciprocity_relative_tolerance"])
        fem_power = float(value["fem_outward_power_w"])
        bem_power = float(value["bem_radiated_power_w"])
        energy_error = float(value["energy_flux_relative_error"])
        energy_tolerance = float(value["energy_flux_relative_tolerance"])
    except (KeyError, TypeError, ValueError):
        return False
    generation = str(value.get("coupling_generation", "")).strip()
    actual_energy_error = abs(fem_power - bem_power) / max(abs(fem_power), 1.0e-300)
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
        and value.get("result_interface_trace_basis") == value.get("interface_trace_basis")
        and len(trace_shape) == 2
        and all(item > 0 for item in trace_shape)
        and trace_shape[0] <= trace_shape[1]
        and result_trace_shape == trace_shape
        and value.get("normal_orientation") == "volume-outward"
        and value.get("result_normal_orientation") == value.get("normal_orientation")
        and math.isfinite(frequency)
        and frequency > 0.0
        and result_frequency == frequency
        and _is_sha256(str(value.get("incident_field_sha256", "")).lower())
        and value.get("result_incident_field_sha256") == value.get("incident_field_sha256")
        and len(pair_ids) == len(result_pair_ids) == len(reciprocity) == len(result_reciprocity) == 2
        and len(set(pair_ids)) == 2
        and result_pair_ids == pair_ids
        and all(len(row) == 2 and all(math.isfinite(item) for item in row) for row in reciprocity)
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
        and math.isclose(energy_error, actual_energy_error, rel_tol=1.0e-6, abs_tol=1.0e-15)
        and math.isfinite(energy_tolerance)
        and energy_tolerance > 0.0
        and energy_error <= energy_tolerance
        and _is_sha256(str(value.get("coupled_result_sha256", "")).lower())
        and value.get("accepted_coupled_result_sha256") == value.get("coupled_result_sha256")
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
