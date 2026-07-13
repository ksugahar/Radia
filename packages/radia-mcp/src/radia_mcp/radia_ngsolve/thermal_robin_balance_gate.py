"""Solver-neutral validation for multi-boundary steady Robin heat flow."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


_ROLES = ("plateau_a", "plateau_b", "fine", "fine_repeat")


def _finite(value: object, default: float = math.nan) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def _mapping(value: object) -> Mapping[str, object]:
    return value if isinstance(value, Mapping) else {}


def _rows(value: object) -> list[Mapping[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        return []
    return [row for row in value if isinstance(row, Mapping)]


def _relative_error(left: float, right: float) -> float:
    if not math.isfinite(left) or not math.isfinite(right):
        return math.inf
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def _maximum_finite(mapping: Mapping[str, object]) -> float:
    values = [_finite(value, math.inf) for value in mapping.values()]
    return max(values, default=math.inf)


def thermal_robin_boundary_balance_gate(
    summary: Mapping[str, object],
    *,
    maximum_balance_relative: float = 1.0e-4,
    maximum_internal_cut_relative_error: float = 1.0e-4,
    maximum_refinement_relative_change: float = 5.0e-5,
    maximum_repeat_relative_error: float = 1.0e-12,
    maximum_symmetry_relative_error: float = 1.0e-4,
    maximum_constitutive_relative_error: float = 1.0e-12,
    maximum_reflection_temperature_relative_error: float = 1.0e-8,
    maximum_reflection_flux_relative_error: float = 2.0e-8,
) -> dict[str, object]:
    """Gate a solved 2-D/3-D Robin heat-flow result using independent identities.

    ``boundary_groups`` must contain signed outward heat rates. ``mesh_ladder``
    carries four role-labelled rows: two possible adaptive-mesher plateau rows,
    a refined row, and its exact replay. This explicitly avoids assuming that a
    smaller requested element size always changes the generated mesh.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be an object")
    boundaries = _rows(summary.get("boundary_groups"))
    ladder_rows = _rows(summary.get("mesh_ladder"))
    by_role = {str(row.get("role", "")): row for row in ladder_rows}
    if len(boundaries) < 2:
        raise ValueError("boundary_groups must contain at least two rows")
    if set(by_role) != set(_ROLES) or len(ladder_rows) != len(_ROLES):
        raise ValueError("mesh_ladder must contain each required role exactly once")

    rates = [_finite(row.get("heat_rate_W")) for row in boundaries]
    if not all(math.isfinite(value) for value in rates):
        raise ValueError("every boundary heat_rate_W must be finite")
    throughput = 0.5 * sum(abs(value) for value in rates)
    if throughput <= 0.0:
        raise ValueError("boundary heat-rate throughput must be positive")
    balance_relative = abs(sum(rates)) / throughput
    internal_cut = abs(_finite(summary.get("internal_cut_heat_rate_W")))
    internal_cut_relative_error = _relative_error(internal_cut, throughput)

    plateau_a = by_role["plateau_a"]
    plateau_b = by_role["plateau_b"]
    fine = by_role["fine"]
    repeat = by_role["fine_repeat"]

    def integer(row: Mapping[str, object], key: str) -> int:
        value = _finite(row.get(key), -1.0)
        return int(value) if value >= 0.0 and value.is_integer() else -1

    def observable_errors(
        left: Mapping[str, object], right: Mapping[str, object]
    ) -> dict[str, float]:
        return {
            key: _relative_error(_finite(left.get(key)), _finite(right.get(key)))
            for key in ("average_temperature_K", "robin_throughput_W")
        }

    plateau_errors = observable_errors(plateau_a, plateau_b)
    refinement_errors = observable_errors(plateau_b, fine)
    repeat_errors = observable_errors(fine, repeat)
    reflection = _mapping(summary.get("temperature_reflection"))
    reflection_temperature = _maximum_finite(
        _mapping(reflection.get("temperature_relative_errors"))
    )
    reflection_flux = _maximum_finite(
        _mapping(reflection.get("flux_relative_errors"))
    )
    fine_balance = _finite(fine.get("balance_relative"), math.inf)
    plateau_balance = _finite(plateau_b.get("balance_relative"), math.inf)
    symmetry_error = _finite(summary.get("symmetry_relative_error"), math.inf)
    constitutive_error = _finite(
        summary.get("constitutive_relative_error"), math.inf
    )

    checks = {
        "opposed_boundary_heat_rate_signs_exist": min(rates) < 0.0 < max(rates),
        "signed_robin_energy_balance_closes": balance_relative
        <= float(maximum_balance_relative),
        "internal_cut_matches_boundary_throughput": internal_cut_relative_error
        <= float(maximum_internal_cut_relative_error),
        "adaptive_mesh_plateau_is_explicit": integer(plateau_a, "node_count")
        == integer(plateau_b, "node_count")
        and integer(plateau_a, "element_count")
        == integer(plateau_b, "element_count")
        and max(plateau_errors.values()) <= maximum_repeat_relative_error,
        "refined_mesh_leaves_plateau_and_stabilizes_observables": integer(
            fine, "node_count"
        )
        > integer(plateau_b, "node_count")
        and integer(fine, "element_count") > integer(plateau_b, "element_count")
        and max(refinement_errors.values()) <= maximum_refinement_relative_change
        and fine_balance <= 0.5 * plateau_balance,
        "refined_replay_is_exact": integer(fine, "node_count")
        == integer(repeat, "node_count")
        and integer(fine, "element_count") == integer(repeat, "element_count")
        and max(repeat_errors.values()) <= maximum_repeat_relative_error,
        "exact_boundary_nonfinite_flux_is_rejected": summary.get(
            "exact_boundary_flux_status"
        )
        == "rejected_nonfinite_boundary_flux",
        "symmetry_identity_closes": symmetry_error
        <= float(maximum_symmetry_relative_error),
        "constitutive_flux_gradient_identity_closes": constitutive_error
        <= float(maximum_constitutive_relative_error),
        "temperature_reflection_covariance_closes": reflection_temperature
        <= float(maximum_reflection_temperature_relative_error)
        and reflection_flux <= float(maximum_reflection_flux_relative_error),
    }
    issues = [name for name, passed in checks.items() if not passed]
    return {
        "policy": "thermal_robin_boundary_balance_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "boundary_group_count": len(boundaries),
            "robin_throughput_W": throughput,
            "balance_relative": balance_relative,
            "internal_cut_relative_error": internal_cut_relative_error,
            "maximum_plateau_relative_error": max(plateau_errors.values()),
            "maximum_refinement_relative_change": max(refinement_errors.values()),
            "maximum_repeat_relative_error": max(repeat_errors.values()),
            "maximum_reflection_temperature_relative_error": reflection_temperature,
            "maximum_reflection_flux_relative_error": reflection_flux,
        },
        "notes": [
            "Integrate Robin heat rates from h*(T-T_inf); do not trust a nonfinite exact-boundary F.n sample.",
            "Treat an unchanged adaptive mesh as an explicit plateau, then prove that a finer mesh leaves it.",
            "Use signed balance, an independent internal cut, exact replay, and affine temperature reflection together.",
        ],
    }
