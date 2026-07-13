"""Rayleigh-monotonicity gate for conductive contact networks."""
from __future__ import annotations

import math


def _relative(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def conductive_network_resistance_monotonicity_gate(
    summary: dict,
    *,
    current_balance_rtol: float = 1.0e-4,
    power_balance_rtol: float = 1.0e-6,
    replay_rtol: float = 1.0e-9,
    minimum_drop_to_error_ratio: float = 5.0,
    maximum_log10_residual: float = -5.0,
) -> dict:
    """Gate resistance decrease as conductive contact paths are added.

    Rayleigh monotonicity applies only when terminal definitions, material
    properties and boundary conditions are fixed.  Each topology must first
    close current conservation and terminal power against Joule loss.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be an object")
    cases = summary.get("cases")
    if not isinstance(cases, list) or len(cases) < 3:
        raise ValueError("cases must contain at least three topology rows")
    tolerances = (
        current_balance_rtol,
        power_balance_rtol,
        replay_rtol,
        minimum_drop_to_error_ratio,
    )
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")
    if not math.isfinite(maximum_log10_residual):
        raise ValueError("maximum_log10_residual must be finite")

    rows = []
    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            raise ValueError(f"case {index} must be an object")
        case_id = str(case.get("case_id") or "").strip()
        if not case_id:
            raise ValueError(f"case {index} needs case_id")
        count = int(case.get("contacting_conductor_count", 0))
        resistance = float(case.get("effective_resistance_ohm", math.nan))
        terminal_power = float(case.get("terminal_power_W", math.nan))
        joule_loss = float(case.get("joule_loss_W", math.nan))
        current_error = float(case.get("current_balance_relative_error", math.nan))
        adaptive_error = float(case.get("adaptive_relative_error", math.nan))
        residual = float(case.get("final_log10_residual", math.nan))
        dof = int(case.get("final_dof", 0))
        solve_s = float(case.get("solve_s", math.nan))
        finite = all(
            math.isfinite(value)
            for value in (
                resistance,
                terminal_power,
                joule_loss,
                current_error,
                adaptive_error,
                residual,
                solve_s,
            )
        )
        power_error = (
            _relative(terminal_power, joule_loss) if finite else math.inf
        )
        checks = {
            "physical_values_are_finite": finite,
            "contact_count_resistance_and_effort_are_positive": count > 0
            and resistance > 0.0
            and dof > 0
            and solve_s > 0.0,
            "two_terminal_current_conserves": current_error <= current_balance_rtol,
            "terminal_power_closes_joule_loss": power_error <= power_balance_rtol,
            "adaptive_error_is_nonnegative": adaptive_error >= 0.0,
            "linear_residual_is_below_limit": residual <= maximum_log10_residual,
        }
        rows.append(
            {
                "case_id": case_id,
                "contacting_conductor_count": count,
                "effective_resistance_ohm": resistance,
                "terminal_power_loss_relative_error": power_error,
                "current_balance_relative_error": current_error,
                "adaptive_relative_error": adaptive_error,
                "final_log10_residual": residual,
                "final_dof": dof,
                "solve_s": solve_s,
                "checks": checks,
                "status": "ok" if all(checks.values()) else "needs_attention",
            }
        )

    contact_counts = [row["contacting_conductor_count"] for row in rows]
    resistances = [row["effective_resistance_ohm"] for row in rows]
    adjacent = []
    for left, right in zip(rows, rows[1:]):
        drop = (
            left["effective_resistance_ohm"] - right["effective_resistance_ohm"]
        ) / left["effective_resistance_ohm"]
        error_scale = max(
            left["adaptive_relative_error"], right["adaptive_relative_error"], 1.0e-300
        )
        adjacent.append(
            {
                "from": left["case_id"],
                "to": right["case_id"],
                "resistance_drop_relative": drop,
                "drop_to_adaptive_error_ratio": drop / error_scale,
            }
        )
    replay = float(summary.get("replay_max_relative_error", math.inf))
    timing = summary.get("timing_breakdown_s")
    timing_ok = False
    if isinstance(timing, dict) and len(timing) == 4:
        try:
            timing_ok = all(
                math.isfinite(float(value)) and float(value) >= 0.0
                for value in timing.values()
            )
        except (TypeError, ValueError):
            timing_ok = False
    checks = {
        "case_ids_are_unique": len({row["case_id"] for row in rows}) == len(rows),
        "all_topologies_are_internally_valid": all(row["status"] == "ok" for row in rows),
        "contacting_conductor_count_strictly_increases": all(
            right > left for left, right in zip(contact_counts, contact_counts[1:])
        ),
        "effective_resistance_strictly_decreases": all(
            right < left for left, right in zip(resistances, resistances[1:])
        ),
        "all_resistance_drops_exceed_adaptive_error": all(
            row["drop_to_adaptive_error_ratio"] >= minimum_drop_to_error_ratio
            for row in adjacent
        ),
        "independent_replay_is_deterministic": math.isfinite(replay)
        and replay <= replay_rtol,
        "exactly_four_timing_stages": timing_ok,
    }
    return {
        "policy": "conductive_network_resistance_monotonicity_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "cases": rows,
        "adjacent_topology_changes": adjacent,
        "metrics": {
            "maximum_current_balance_relative_error": max(
                row["current_balance_relative_error"] for row in rows
            ),
            "maximum_terminal_power_loss_relative_error": max(
                row["terminal_power_loss_relative_error"] for row in rows
            ),
            "minimum_drop_to_adaptive_error_ratio": min(
                row["drop_to_adaptive_error_ratio"] for row in adjacent
            ),
            "maximum_replay_relative_error": replay,
            "dof_growth": rows[-1]["final_dof"] / rows[0]["final_dof"],
            "solve_time_growth": rows[-1]["solve_s"] / rows[0]["solve_s"],
        },
        "lesson": (
            "With fixed terminals, materials and boundaries, adding conductive "
            "contact paths must not increase effective resistance. First close "
            "current and Joule power, then require each resistance drop to exceed "
            "the adaptive-discretization error and verify a fresh replay."
        ),
    }
