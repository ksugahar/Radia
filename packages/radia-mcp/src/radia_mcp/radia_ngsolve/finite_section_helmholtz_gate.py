"""Validation gate for a finite-section Helmholtz coil and mixed Omega solve."""
from __future__ import annotations

import math
from functools import lru_cache
from typing import Any, Mapping, Sequence

MU0 = 4.0e-7 * math.pi


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    if positive and number <= 0.0:
        raise ValueError(f"{name} must be positive")
    return number


def _vector(values: object, name: str, count: int) -> list[float]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    result = [_finite(value, f"{name}[{index}]") for index, value in enumerate(values)]
    if len(result) != count:
        raise ValueError(f"{name} must contain {count} values")
    return result


def _field(values: object, name: str, count: int) -> list[list[float]]:
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    result = [_vector(value, f"{name}[{index}]", 3) for index, value in enumerate(values)]
    if len(result) != count:
        raise ValueError(f"{name} must contain {count} vectors")
    return result


def _positive_int(value: object, name: str) -> int:
    number = _finite(value, name, positive=True)
    integer = int(number)
    if float(integer) != number:
        raise ValueError(f"{name} must be an integer")
    return integer


@lru_cache(maxsize=8)
def _gauss_legendre(order: int) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if order < 8 or order > 128:
        raise ValueError("quadrature order must be between 8 and 128")
    nodes = [0.0] * order
    weights = [0.0] * order
    half = (order + 1) // 2
    for index in range(half):
        root = math.cos(math.pi * (index + 0.75) / (order + 0.5))
        derivative = 0.0
        for _ in range(64):
            p0, p1 = 1.0, root
            for degree in range(2, order + 1):
                p0, p1 = p1, ((2 * degree - 1) * root * p1 - (degree - 1) * p0) / degree
            derivative = order * (root * p1 - p0) / (root * root - 1.0)
            update = p1 / derivative
            root -= update
            if abs(update) <= 4.0 * math.ulp(1.0):
                break
        weight = 2.0 / ((1.0 - root * root) * derivative * derivative)
        nodes[index] = -root
        nodes[order - 1 - index] = root
        weights[index] = weight
        weights[order - 1 - index] = weight
    return tuple(nodes), tuple(weights)


def finite_section_helmholtz_axis_bz(
    z_m: float,
    *,
    radial_bounds_m: Sequence[float],
    positive_axial_bounds_m: Sequence[float],
    current_density_A_per_m2: float,
    quadrature_order: int = 48,
) -> float:
    """Return the exact-on-axis integral for a mirrored rectangular coil pair."""
    z = _finite(z_m, "z_m")
    radial = _vector(radial_bounds_m, "radial_bounds_m", 2)
    axial = _vector(positive_axial_bounds_m, "positive_axial_bounds_m", 2)
    current_density = _finite(
        current_density_A_per_m2, "current_density_A_per_m2", positive=True
    )
    if not (0.0 < radial[0] < radial[1]):
        raise ValueError("radial_bounds_m must be strictly increasing and positive")
    if not (0.0 < axial[0] < axial[1]):
        raise ValueError("positive_axial_bounds_m must be strictly increasing and positive")
    nodes, weights = _gauss_legendre(_positive_int(quadrature_order, "quadrature_order"))
    r_mid, r_half = 0.5 * sum(radial), 0.5 * (radial[1] - radial[0])
    z_mid, z_half = 0.5 * sum(axial), 0.5 * (axial[1] - axial[0])
    integral = 0.0
    for r_node, r_weight in zip(nodes, weights):
        radius = r_mid + r_half * r_node
        for z_node, z_weight in zip(nodes, weights):
            source_z = z_mid + z_half * z_node
            positive = radius * radius / (radius * radius + (z - source_z) ** 2) ** 1.5
            negative = radius * radius / (radius * radius + (z + source_z) ** 2) ** 1.5
            integral += r_weight * z_weight * (positive + negative)
    return 0.5 * MU0 * current_density * r_half * z_half * integral


def _relative_errors(values: Sequence[float], reference: Sequence[float]) -> list[float]:
    return [abs(value - exact) / abs(exact) for value, exact in zip(values, reference)]


def divergence_free_source_assembly_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Gate a volume-to-boundary source-load acceleration by solved fields."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    required = ("manufactured_load_relative_error", "solved_field_relative_difference",
                "allowed_solved_field_relative_difference", "volume_rhs_seconds",
                "surface_rhs_seconds")
    values = {name: _finite(summary.get(name), name, positive=True) for name in required}
    checks = {
        "schema_is_supported": summary.get("schema") == "radia.divergence-free-source-assembly.v1",
        "scope_is_affine_p1_constant_mu": summary.get("scope") == "affine_p1_constant_mu",
        "source_is_divergence_free": summary.get("source_identity") == "divergence_free",
        "boundary_orientation_is_element_derived": summary.get("boundary_orientation") == "opposite_vertex",
        "manufactured_load_matches": values["manufactured_load_relative_error"] <= 1.0e-12,
        "solved_field_matches": values["solved_field_relative_difference"] <= values["allowed_solved_field_relative_difference"],
        "surface_rhs_is_faster": values["surface_rhs_seconds"] < values["volume_rhs_seconds"],
    }
    return {
        "policy": "divergence_free_source_assembly_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "speedup": values["volume_rhs_seconds"] / values["surface_rhs_seconds"],
        "claim_limit": "Accepted only for affine P1, constant permeability and a verified divergence-free source; load-vector L2 agreement alone is not an acceptance gate.",
    }


def finite_section_helmholtz_mixed_omega_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Gate source construction, analytic accuracy, open boundary, and p/h convergence."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    units = summary.get("units")
    geometry = summary.get("geometry")
    source = summary.get("source")
    open_boundary = summary.get("open_boundary")
    tolerances = summary.get("gate_tolerances")
    mixed = summary.get("mixed_total_reduced")
    if not all(isinstance(item, Mapping) for item in (units, geometry, source, open_boundary, tolerances, mixed)):
        raise ValueError("units, geometry, source, open_boundary, gate_tolerances, and mixed_total_reduced must be objects")

    raw_axis = summary.get("axis_z_m")
    if not isinstance(raw_axis, Sequence) or isinstance(raw_axis, (str, bytes)):
        raise ValueError("axis_z_m must be an array")
    axis = _vector(raw_axis, "axis_z_m", len(raw_axis))
    if len(axis) < 3 or any(right <= left for left, right in zip(axis, axis[1:])):
        raise ValueError("axis_z_m must contain at least three strictly increasing samples")
    radial = _vector(geometry.get("radial_bounds_m"), "geometry.radial_bounds_m", 2)
    axial = _vector(geometry.get("positive_axial_bounds_m"), "geometry.positive_axial_bounds_m", 2)
    current_density = _finite(
        geometry.get("current_density_A_per_m2"),
        "geometry.current_density_A_per_m2",
        positive=True,
    )
    exact = [
        finite_section_helmholtz_axis_bz(
            z,
            radial_bounds_m=radial,
            positive_axial_bounds_m=axial,
            current_density_A_per_m2=current_density,
        )
        for z in axis
    ]
    supplied_reference = _vector(summary.get("reference_bz_T"), "reference_bz_T", len(axis))
    coil = _field(summary.get("coil_builder_field_T"), "coil_builder_field_T", len(axis))
    final_field = _field(mixed.get("field_T"), "mixed_total_reduced.field_T", len(axis))
    if any(value[2] == 0.0 for value in [*coil, *final_field]):
        raise ValueError("coil_builder_field_T and mixed_total_reduced.field_T require nonzero Bz")
    reference_error = _relative_errors(supplied_reference, exact)
    coil_error = _relative_errors([value[2] for value in coil], exact)
    mixed_error = _relative_errors([value[2] for value in final_field], exact)
    coil_transverse = max(math.hypot(value[0], value[1]) / abs(value[2]) for value in coil)
    mixed_transverse = max(math.hypot(value[0], value[1]) / abs(value[2]) for value in final_field)

    levels = mixed.get("refinement_levels")
    if not isinstance(levels, list) or len(levels) < 3:
        raise ValueError("mixed_total_reduced.refinement_levels must contain at least three levels")
    parsed_levels = []
    for index, level in enumerate(levels):
        if not isinstance(level, Mapping):
            raise ValueError(f"refinement_levels[{index}] must be an object")
        parsed_levels.append(
            {
                "mesh_maxh_m": _finite(level.get("mesh_maxh_m"), f"refinement_levels[{index}].mesh_maxh_m", positive=True),
                "order": _positive_int(level.get("order"), f"refinement_levels[{index}].order"),
                "cell_count": _positive_int(level.get("cell_count"), f"refinement_levels[{index}].cell_count"),
                "ndof": _positive_int(level.get("ndof"), f"refinement_levels[{index}].ndof"),
                "maximum_relative_error": _finite(level.get("maximum_relative_error"), f"refinement_levels[{index}].maximum_relative_error", positive=True),
            }
        )
    reported_errors = [level["maximum_relative_error"] for level in parsed_levels]
    final_reported_consistency = abs(reported_errors[-1] - max(mixed_error)) <= max(1.0e-12, 1.0e-6 * max(mixed_error))

    limits = {
        "reference": _finite(tolerances.get("maximum_reference_relative_error"), "maximum_reference_relative_error", positive=True),
        "coil": _finite(tolerances.get("maximum_coil_builder_relative_error"), "maximum_coil_builder_relative_error", positive=True),
        "mixed": _finite(tolerances.get("maximum_mixed_relative_error"), "maximum_mixed_relative_error", positive=True),
        "coil_transverse": _finite(tolerances.get("maximum_coil_transverse_ratio"), "maximum_coil_transverse_ratio", positive=True),
        "mixed_transverse": _finite(tolerances.get("maximum_mixed_transverse_ratio"), "maximum_mixed_transverse_ratio", positive=True),
    }
    checks = {
        "schema_is_supported": summary.get("schema") == "radia.finite-section-helmholtz-mixed-omega.v1",
        "units_are_explicit_si": units.get("length") == "m" and units.get("field") == "T" and units.get("current_density") == "A/m^2",
        "geometry_is_mirrored_finite_section": geometry.get("mirror_pair") is True and 0.0 < radial[0] < radial[1] and 0.0 < axial[0] < axial[1],
        "source_is_coil_builder_finite_section": source.get("representation") == "coil_builder_finite_cross_section" and source.get("closed_paths") is True,
        "formulation_is_mixed_total_reduced_omega": mixed.get("formulation") == "mixed_total_reduced_omega",
        "open_boundary_uses_kelvin_transform": open_boundary.get("method") == "kelvin_transform" and _finite(open_boundary.get("radius_m"), "open_boundary.radius_m", positive=True) > max(radial[1], axial[1]),
        "stored_reference_matches_analytic_integral": max(reference_error) <= limits["reference"],
        "coil_builder_matches_analytic_integral": max(coil_error) <= limits["coil"],
        "coil_builder_preserves_axis_symmetry": coil_transverse <= limits["coil_transverse"],
        "mixed_solution_matches_analytic_integral": max(mixed_error) <= limits["mixed"],
        "mixed_solution_preserves_axis_symmetry": mixed_transverse <= limits["mixed_transverse"],
        "refinement_effort_increases": all(right["mesh_maxh_m"] <= left["mesh_maxh_m"] and right["order"] >= left["order"] and right["ndof"] > left["ndof"] for left, right in zip(parsed_levels, parsed_levels[1:])),
        "refinement_error_decreases": all(right < left for left, right in zip(reported_errors, reported_errors[1:])),
        "final_refinement_matches_field": final_reported_consistency,
    }
    return {
        "policy": "finite_section_helmholtz_mixed_omega_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "sample_count": len(axis),
            "analytic_bz_T": exact,
            "maximum_reference_relative_error": max(reference_error),
            "maximum_coil_builder_relative_error": max(coil_error),
            "maximum_coil_transverse_ratio": coil_transverse,
            "maximum_mixed_relative_error": max(mixed_error),
            "maximum_mixed_transverse_ratio": mixed_transverse,
            "refinement_maximum_relative_errors": reported_errors,
            "final_order": parsed_levels[-1]["order"],
            "final_ndof": parsed_levels[-1]["ndof"],
        },
        "lesson": (
            "A mixed total/reduced Omega result is accepted only when a finite-cross-section CoilBuilder source, "
            "an independently integrated on-axis field, Kelvin open-boundary treatment, axis symmetry, and "
            "strictly improving p/h evidence agree under one SI geometry identity."
        ),
    }
