"""Finite-difference width selection for motor virtual-work torque checks."""

from __future__ import annotations

import math
from typing import Mapping


def motor_virtual_work_width_ladder_gate(
    summary: Mapping[str, object],
    *,
    max_selected_relative_error: float = 1.0e-2,
    max_mesh_count_relative_span: float = 2.0e-2,
    current_balance_tolerance: float = 1.0e-10,
) -> dict[str, object]:
    """Select a virtual-work angle width against independent direct torque.

    The smallest perturbation is not automatically best when every displaced
    geometry is remeshed.  Require an interior optimum bracketed by a smaller
    noise-dominated width and a larger truncation-dominated width.
    """

    max_error = float(max_selected_relative_error)
    max_mesh_span = float(max_mesh_count_relative_span)
    current_tolerance = float(current_balance_tolerance)
    if any(
        not math.isfinite(value) or value < 0.0
        for value in (max_error, max_mesh_span, current_tolerance)
    ):
        raise ValueError("tolerances must be finite and nonnegative")

    raw_rows = summary.get("virtual_work") or []
    if not isinstance(raw_rows, list) or len(raw_rows) < 3:
        raise ValueError("virtual_work must contain at least three width rows")
    rows = []
    for raw in raw_rows:
        delta = float(raw.get("delta_deg", math.nan))
        direct = float(raw.get("weighted_stress_torque_Nm", math.nan))
        virtual = float(raw.get("coenergy_derivative_torque_Nm", math.nan))
        error = abs(virtual - direct) / max(abs(direct), 1.0e-30)
        rows.append(
            {
                "delta_deg": delta,
                "weighted_stress_torque_Nm": direct,
                "coenergy_derivative_torque_Nm": virtual,
                "relative_error": error,
                "reported_relative_error": float(raw.get("relative_error", math.nan)),
            }
        )
    rows.sort(key=lambda row: row["delta_deg"])
    best = min(rows, key=lambda row: row["relative_error"])
    widths = [float(row["delta_deg"]) for row in rows]
    direct_values = [float(row["weighted_stress_torque_Nm"]) for row in rows]
    excitation = summary.get("excitation") or {}
    phase_currents = excitation.get("phase_currents_A") or {}
    current_sum = sum(float(value) for value in phase_currents.values())
    current_square_sum = sum(float(value) ** 2 for value in phase_currents.values())
    expected_square_sum = float(excitation.get("expected_square_sum_A2", math.nan))
    selected_delta = float(summary.get("selected_virtual_work_delta_deg", math.nan))
    selected_index = widths.index(best["delta_deg"])
    mesh_span = float(summary.get("mesh_element_count_relative_span", math.inf))
    checks = {
        "widths_are_positive_unique_and_increasing": all(width > 0.0 for width in widths)
        and len(set(widths)) == len(widths),
        "independent_direct_torque_is_nonzero_and_consistent": all(
            math.isfinite(value) and abs(value) > 1.0e-9 for value in direct_values
        )
        and max(direct_values) - min(direct_values)
        <= 1.0e-12 * max(abs(direct_values[0]), 1.0),
        "reported_errors_match_recomputed_errors": all(
            math.isfinite(row["reported_relative_error"])
            and abs(row["reported_relative_error"] - row["relative_error"]) <= 1.0e-12
            for row in rows
        ),
        "balanced_three_phase_current_sum": len(phase_currents) == 3
        and abs(current_sum) <= current_tolerance,
        "balanced_three_phase_square_sum": math.isfinite(expected_square_sum)
        and abs(current_square_sum - expected_square_sum) <= current_tolerance,
        "selected_width_is_error_minimum": math.isfinite(selected_delta)
        and selected_delta == best["delta_deg"],
        "selected_width_is_interior": 0 < selected_index < len(rows) - 1,
        "smaller_width_exposes_remesh_noise": 0 < selected_index
        and rows[0]["relative_error"] > best["relative_error"],
        "larger_width_exposes_truncation": selected_index < len(rows) - 1
        and rows[-1]["relative_error"] > best["relative_error"],
        "selected_virtual_work_matches_direct_torque": best["relative_error"] <= max_error,
        "mesh_count_span_is_bounded": mesh_span <= max_mesh_span,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "motor_virtual_work_width_ladder_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "selected_delta_deg": best["delta_deg"],
        "selected_relative_error": best["relative_error"],
        "mesh_element_count_relative_span": mesh_span,
        "rows": rows,
        "notes": [
            "Choose displacement width by agreement with an independent torque route, not by smallest delta.",
            "A smaller central difference can be worse when each displaced geometry is remeshed independently.",
            "Keep at least one width on either side of the selected optimum to expose noise and truncation.",
        ],
    }
