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
