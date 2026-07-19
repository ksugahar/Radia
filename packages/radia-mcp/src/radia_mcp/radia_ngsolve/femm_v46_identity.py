from __future__ import annotations

import math
from collections.abc import Mapping


def _sha(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _closed(row: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(row.get("generation", "")).strip()
    return bool(generation) and all(row.get(field) == generation for field in fields)


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    force = identity.get(
        "v46_public_axisymmetric_force_weighted_stress_contour_open_boundary_partial_mismatch"
    )
    if isinstance(force, Mapping):
        checks["femm_v46_force_generation_closure"] = _closed(
            force,
            ("force_generation", "contour_generation", "open_boundary_generation", "solve_generation", "result_generation"),
        )
        checks["femm_v46_force_method_boundary_state"] = (
            force.get("force_method") == force.get("result_force_method") == "weighted_stress_tensor"
            and force.get("open_boundary") == force.get("result_open_boundary") == "outer_kelvin"
            and force.get("solve_state") == force.get("result_solve_state") == "complete"
            and force.get("partial_solve_status") == force.get("result_partial_solve_status") == "none"
        )
        checks["femm_v46_force_values_owner_digest"] = (
            isinstance(force.get("force_n"), list)
            and force.get("force_n") == force.get("result_force_n")
            and all(math.isfinite(float(value)) for value in force["force_n"])
            and math.isclose(float(force.get("axisymmetric_factor")), 2.0 * math.pi, rel_tol=1.0e-12)
            and force.get("result_axisymmetric_factor") == force.get("axisymmetric_factor")
            and str(force.get("contour_owner", "")).startswith("contour:")
            and force.get("result_contour_owner") == force.get("contour_owner")
            and _sha(force.get("result_sha256"))
            and force.get("accepted_result_sha256") == force.get("result_sha256")
        )

    nonlinear = identity.get(
        "v46_public_nonlinear_bh_curve_unit_scale_convergence_status_nan_mismatch"
    )
    if isinstance(nonlinear, Mapping):
        checks["femm_v46_bh_generation_closure"] = _closed(
            nonlinear,
            ("bh_curve_generation", "unit_generation", "convergence_generation", "result_generation"),
        )
        points = nonlinear.get("bh_points")
        checks["femm_v46_bh_units_convergence_finite"] = (
            nonlinear.get("bh_unit") == nonlinear.get("result_bh_unit") == "tesla_ampere_per_meter"
            and nonlinear.get("convergence_status") == nonlinear.get("result_convergence_status") == "converged"
            and nonlinear.get("nonfinite_status") == nonlinear.get("result_nonfinite_status") == "none"
            and isinstance(points, list)
            and points == nonlinear.get("result_bh_points")
            and all(isinstance(pair, list) and len(pair) == 2 and all(math.isfinite(float(value)) for value in pair) for pair in points)
        )
        checks["femm_v46_bh_owner_digest"] = (
            str(nonlinear.get("material_owner", "")).startswith("material:")
            and nonlinear.get("result_material_owner") == nonlinear.get("material_owner")
            and _sha(nonlinear.get("result_sha256"))
            and nonlinear.get("accepted_result_sha256") == nonlinear.get("result_sha256")
        )
    return checks


def validate_source_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    problem = identity.get("v46_source_tool_fem_problem_depth_axis_boundary_material_reload_mismatch")
    if isinstance(problem, Mapping):
        checks["femm_v46_source_problem_generation_closure"] = _closed(
            problem,
            ("problem_generation", "depth_generation", "axis_generation", "boundary_generation", "material_generation", "reload_generation", "restart_generation", "result_generation"),
        )
        checks["femm_v46_source_problem_replay"] = (
            problem.get("problem_type") == problem.get("replayed_problem_type") == "axisymmetric"
            and float(problem.get("depth_m")) > 0.0
            and problem.get("replayed_depth_m") == problem.get("depth_m")
            and problem.get("axis_boundary") == problem.get("replayed_axis_boundary") == "outer_kelvin"
            and problem.get("material_reload_status") == problem.get("replayed_material_reload_status") == "reloaded"
            and problem.get("restart_state") == problem.get("replayed_restart_state") == "cold_start"
        )
        checks["femm_v46_source_problem_owner_digest"] = (
            str(problem.get("document_owner", "")).startswith("document:")
            and problem.get("replayed_document_owner") == problem.get("document_owner")
            and str(problem.get("solution_owner", "")).startswith("solution:")
            and problem.get("replayed_solution_owner") == problem.get("solution_owner")
            and _sha(problem.get("result_sha256"))
            and problem.get("accepted_result_sha256") == problem.get("result_sha256")
        )

    stress = identity.get("v46_source_tool_weighted_stress_contour_sampling_partial_solution_owner_mismatch")
    if isinstance(stress, Mapping):
        checks["femm_v46_source_stress_generation_closure"] = _closed(
            stress,
            ("weighted_stress_generation", "contour_sampling_generation", "partial_solution_generation", "owner_generation", "result_generation"),
        )
        checks["femm_v46_source_stress_sampling_solution"] = (
            stress.get("weighted_stress_tensor") == stress.get("replayed_weighted_stress_tensor") == "mo_blockintegral(19)"
            and int(stress.get("contour_sample_count")) >= 5
            and stress.get("replayed_contour_sample_count") == stress.get("contour_sample_count")
            and stress.get("partial_solution_status") == stress.get("replayed_partial_solution_status") == "complete"
        )
        checks["femm_v46_source_stress_owner_digest"] = (
            str(stress.get("solution_owner", "")).startswith("solution:")
            and stress.get("replayed_solution_owner") == stress.get("solution_owner")
            and _sha(stress.get("result_sha256"))
            and stress.get("accepted_result_sha256") == stress.get("result_sha256")
        )
    return checks
