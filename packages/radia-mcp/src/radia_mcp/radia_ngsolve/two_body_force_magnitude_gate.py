"""Solver-neutral gate for unsigned two-body force tables and fresh replay."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(result) or (positive and result <= 0.0):
        qualifier = "positive and finite" if positive else "finite"
        raise ValueError(f"{name} must be {qualifier}")
    return result


def _relative(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def two_body_force_magnitude_replay_gate(
    summary: Mapping[str, object],
    *,
    max_body_balance_relative_error: float = 5.0e-4,
    max_replay_relative_error: float = 1.0e-9,
    max_current_relative_error: float = 1.0e-9,
) -> dict[str, Any]:
    """Gate magnitude balance without inventing a missing force-vector sign."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    limits = {
        "body_balance": _finite(
            max_body_balance_relative_error, "max_body_balance_relative_error"
        ),
        "replay": _finite(max_replay_relative_error, "max_replay_relative_error"),
        "current": _finite(max_current_relative_error, "max_current_relative_error"),
    }
    if any(value < 0.0 for value in limits.values()):
        raise ValueError("relative tolerances must be nonnegative")
    commanded_current = _finite(
        summary.get("commanded_current_a"), "commanded_current_a", positive=True
    )
    runs = summary.get("runs")
    if not isinstance(runs, Sequence) or isinstance(runs, (str, bytes)) or len(runs) != 2:
        raise ValueError("runs must contain exactly two replay records")

    parsed: list[dict[str, float | int | bool]] = []
    metadata_ok = True
    for index, raw in enumerate(runs):
        if not isinstance(raw, Mapping):
            raise ValueError(f"runs[{index}] must be an object")
        parsed.append(
            {
                "replay": int(raw.get("replay", -1)),
                "body_a_force_n": _finite(
                    raw.get("body_a_force_magnitude_n"),
                    f"runs[{index}].body_a_force_magnitude_n",
                    positive=True,
                ),
                "body_b_force_n": _finite(
                    raw.get("body_b_force_magnitude_n"),
                    f"runs[{index}].body_b_force_magnitude_n",
                    positive=True,
                ),
                "current_a": _finite(raw.get("current_a"), f"runs[{index}].current_a"),
                "flux_wb": _finite(
                    raw.get("flux_wb"), f"runs[{index}].flux_wb", positive=True
                ),
                "element_count": int(raw.get("element_count", 0)),
                "vertex_count": int(raw.get("vertex_count", 0)),
                "solver_runtime_s": _finite(
                    raw.get("solver_runtime_s"),
                    f"runs[{index}].solver_runtime_s",
                    positive=True,
                ),
                "fresh_result": raw.get("fresh_result") is True,
            }
        )
        metadata_ok = metadata_ok and raw.get("force_unit") == "N"
        metadata_ok = metadata_ok and raw.get("flux_unit") == "Wb"
        metadata_ok = metadata_ok and raw.get("current_unit") == "A"

    balance_errors = [
        _relative(float(row["body_a_force_n"]), float(row["body_b_force_n"]))
        for row in parsed
    ]
    replay_keys = ("body_a_force_n", "body_b_force_n", "current_a", "flux_wb")
    replay_errors = [
        _relative(float(parsed[0][key]), float(parsed[1][key])) for key in replay_keys
    ]
    current_errors = [
        _relative(float(row["current_a"]), commanded_current) for row in parsed
    ]
    checks = {
        "unsigned_force_quantity_recorded": summary.get("force_quantity")
        == "unsigned_magnitude",
        "action_reaction_sign_not_inferred": summary.get("action_reaction_sign_inferred")
        is False,
        "force_flux_current_units_recorded": metadata_ok,
        "two_ordered_replays_recorded": [row["replay"] for row in parsed] == [1, 2],
        "both_force_magnitudes_are_positive": all(
            float(row["body_a_force_n"]) > 0.0 and float(row["body_b_force_n"]) > 0.0
            for row in parsed
        ),
        "two_body_force_magnitudes_balance": max(balance_errors)
        <= limits["body_balance"],
        "current_matches_command": max(current_errors) <= limits["current"],
        "force_current_flux_replay_is_stable": max(replay_errors) <= limits["replay"],
        "same_positive_mesh_inventory_reused": min(
            int(row["element_count"]) for row in parsed
        )
        > 0
        and min(int(row["vertex_count"]) for row in parsed) > 0
        and len(
            {(int(row["element_count"]), int(row["vertex_count"])) for row in parsed}
        )
        == 1,
        "both_results_are_fresh_solver_outputs": all(
            bool(row["fresh_result"]) and float(row["solver_runtime_s"]) > 1.0
            for row in parsed
        ),
    }
    return {
        "policy": "two_body_force_magnitude_replay_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "maximum_body_balance_relative_error": max(balance_errors),
            "maximum_replay_relative_error": max(replay_errors),
            "maximum_current_relative_error": max(current_errors),
            "element_count": int(parsed[0]["element_count"]),
            "vertex_count": int(parsed[0]["vertex_count"]),
            "solver_runtime_s": [float(row["solver_runtime_s"]) for row in parsed],
        },
        "lesson": (
            "A two-body force table that reports magnitudes can gate equal-and-opposite "
            "magnitude balance, but it cannot prove vector sign. Record that limitation, "
            "require two fresh solves on the same mesh, and replay force, current, and flux."
        ),
    }
