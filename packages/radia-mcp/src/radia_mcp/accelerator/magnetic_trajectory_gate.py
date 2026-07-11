"""Solver-independent magnetic trajectory-pair conservation gate."""
from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence


def _number(value: object, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _vector(value: object, name: str) -> tuple[float, float, float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3:
        raise ValueError(f"{name} must contain three components")
    return tuple(_number(component, f"{name} component") for component in value)  # type: ignore[return-value]


def _norm(value: Sequence[float]) -> float:
    return math.sqrt(sum(component * component for component in value))


def _trajectory_map(case: Mapping[str, object], name: str) -> dict[tuple[int, int], dict]:
    rows = case.get("trajectories")
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
        raise ValueError(f"{name}.trajectories must be an array")
    result: dict[tuple[int, int], dict] = {}
    for index, raw in enumerate(rows):
        if not isinstance(raw, Mapping):
            raise ValueError(f"{name}.trajectories[{index}] must be an object")
        try:
            key = (int(raw["emission_id"]), int(raw["particle_id"]))
            point_count = int(raw["point_count"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{name}.trajectories[{index}] needs integer ids and point_count") from exc
        if key in result:
            raise ValueError(f"{name} has duplicate trajectory key {key}")
        if point_count < 2:
            raise ValueError(f"{name}.trajectories[{index}] needs at least two points")
        result[key] = {
            "initial_position": _vector(raw.get("position_initial"), f"{name}.initial_position"),
            "final_position": _vector(raw.get("position_final"), f"{name}.final_position"),
            "final_speed": _number(raw.get("speed_final_m_per_s"), f"{name}.final_speed"),
            "final_energy": _number(raw.get("energy_final_ev"), f"{name}.final_energy"),
        }
    return result


def evaluate_magnetic_trajectory_pair(
    summary: Mapping[str, object],
    min_trajectory_count: int = 3,
    initial_position_absolute_tolerance: float = 1.0e-12,
    minimum_endpoint_deflection: float = 1.0e-8,
    longitudinal_endpoint_absolute_tolerance: float = 1.0e-9,
    maximum_speed_relative_difference: float = 1.0e-6,
    maximum_energy_relative_difference: float = 1.0e-6,
    maximum_current_closure_relative_error: float = 1.0e-12,
    maximum_collision_power_relative_difference: float = 1.0e-6,
) -> dict:
    """Compare otherwise identical trajectories with magnetic field off/on.

    A magnetic field may rotate velocity and displace the endpoint, but the
    Lorentz magnetic term does no work. Therefore speed, kinetic energy,
    transported current, hit count, and collision power should remain closed
    while at least one trajectory changes transversely.
    """
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    if min_trajectory_count < 1:
        raise ValueError("min_trajectory_count must be positive")
    tolerances = (
        initial_position_absolute_tolerance,
        minimum_endpoint_deflection,
        longitudinal_endpoint_absolute_tolerance,
        maximum_speed_relative_difference,
        maximum_energy_relative_difference,
        maximum_current_closure_relative_error,
        maximum_collision_power_relative_difference,
    )
    if any(value < 0.0 for value in tolerances):
        raise ValueError("tolerances must be nonnegative")
    units = {
        "position": str(summary.get("position_unit") or ""),
        "speed": str(summary.get("speed_unit") or ""),
        "energy": str(summary.get("energy_unit") or ""),
    }
    off = summary.get("magnetic_off")
    on = summary.get("magnetic_on")
    if not isinstance(off, Mapping) or not isinstance(on, Mapping):
        raise ValueError("magnetic_off and magnetic_on must be objects")
    off_rows = _trajectory_map(off, "magnetic_off")
    on_rows = _trajectory_map(on, "magnetic_on")
    keys_match = set(off_rows) == set(on_rows)
    common = sorted(set(off_rows) & set(on_rows))
    if not common:
        raise ValueError("the two cases have no common trajectory ids")

    pair_rows = []
    for key in common:
        a = off_rows[key]
        b = on_rows[key]
        initial_delta = tuple(y - x for x, y in zip(a["initial_position"], b["initial_position"]))
        endpoint_delta = tuple(y - x for x, y in zip(a["final_position"], b["final_position"]))
        speed_scale = max(abs(a["final_speed"]), abs(b["final_speed"]), 1.0e-300)
        energy_scale = max(abs(a["final_energy"]), abs(b["final_energy"]), 1.0e-300)
        pair_rows.append({
            "initial_position_difference": _norm(initial_delta),
            "endpoint_difference": _norm(endpoint_delta),
            "endpoint_transverse_difference": math.hypot(endpoint_delta[0], endpoint_delta[1]),
            "endpoint_longitudinal_difference": abs(endpoint_delta[2]),
            "final_speed_relative_difference": abs(a["final_speed"] - b["final_speed"]) / speed_scale,
            "final_energy_relative_difference": abs(a["final_energy"] - b["final_energy"]) / energy_scale,
        })

    accounting = {}
    for name, case in (("magnetic_off", off), ("magnetic_on", on)):
        source_current = _number(case.get("source_current_a"), f"{name}.source_current_a")
        collision_current = _number(case.get("collision_current_a"), f"{name}.collision_current_a")
        collision_power = _number(case.get("collision_power_w"), f"{name}.collision_power_w")
        try:
            hit_count = int(case["boundary_hit_count"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{name}.boundary_hit_count must be an integer") from exc
        accounting[name] = {
            "source_current_a": source_current,
            "collision_current_a": collision_current,
            "collision_power_w": collision_power,
            "boundary_hit_count": hit_count,
            "current_closure_relative_error": abs(collision_current - source_current) / max(abs(source_current), 1.0e-300),
            "power_per_current_v": collision_power / collision_current if abs(collision_current) > 1.0e-300 else math.nan,
        }

    endpoint_max = max(row["endpoint_difference"] for row in pair_rows)
    endpoint_rms = math.sqrt(sum(row["endpoint_difference"] ** 2 for row in pair_rows) / len(pair_rows))
    power_scale = max(abs(accounting["magnetic_off"]["collision_power_w"]), abs(accounting["magnetic_on"]["collision_power_w"]), 1.0e-300)
    power_difference = abs(accounting["magnetic_off"]["collision_power_w"] - accounting["magnetic_on"]["collision_power_w"]) / power_scale
    metrics = {
        "trajectory_count": len(pair_rows),
        "max_initial_position_difference": max(row["initial_position_difference"] for row in pair_rows),
        "endpoint_difference_max": endpoint_max,
        "endpoint_difference_rms": endpoint_rms,
        "endpoint_transverse_difference_max": max(row["endpoint_transverse_difference"] for row in pair_rows),
        "max_endpoint_longitudinal_difference": max(row["endpoint_longitudinal_difference"] for row in pair_rows),
        "max_final_speed_relative_difference": max(row["final_speed_relative_difference"] for row in pair_rows),
        "max_final_energy_relative_difference": max(row["final_energy_relative_difference"] for row in pair_rows),
        "collision_power_relative_difference": power_difference,
        "accounting": accounting,
    }
    checks = {
        "units_explicit": all(units.values()),
        "trajectory_ids_match": keys_match,
        "trajectory_count_sufficient": len(pair_rows) >= min_trajectory_count,
        "initial_positions_match": metrics["max_initial_position_difference"] <= initial_position_absolute_tolerance,
        "magnetic_field_changes_trajectory": metrics["endpoint_transverse_difference_max"] >= minimum_endpoint_deflection,
        "same_longitudinal_endpoint_plane": metrics["max_endpoint_longitudinal_difference"] <= longitudinal_endpoint_absolute_tolerance,
        "magnetic_field_preserves_speed": metrics["max_final_speed_relative_difference"] <= maximum_speed_relative_difference,
        "magnetic_field_preserves_energy": metrics["max_final_energy_relative_difference"] <= maximum_energy_relative_difference,
        "source_collision_current_closes": max(value["current_closure_relative_error"] for value in accounting.values()) <= maximum_current_closure_relative_error,
        "boundary_hit_count_matches_trajectories": all(value["boundary_hit_count"] == len(pair_rows) for value in accounting.values()),
        "collision_power_positive_and_preserved": all(value["collision_power_w"] > 0.0 for value in accounting.values()) and power_difference <= maximum_collision_power_relative_difference,
    }
    return {
        "schema": "radia-accelerator-magnetic-trajectory-pair/v1",
        "policy": "magnetic_deflection_with_no_magnetic_work_and_particle_accounting",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "units": units,
        "metrics": metrics,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "tolerances": {
            "min_trajectory_count": min_trajectory_count,
            "initial_position_absolute": initial_position_absolute_tolerance,
            "minimum_endpoint_deflection": minimum_endpoint_deflection,
            "longitudinal_endpoint_absolute": longitudinal_endpoint_absolute_tolerance,
            "maximum_speed_relative_difference": maximum_speed_relative_difference,
            "maximum_energy_relative_difference": maximum_energy_relative_difference,
            "maximum_current_closure_relative_error": maximum_current_closure_relative_error,
            "maximum_collision_power_relative_difference": maximum_collision_power_relative_difference,
        },
    }


def magnetic_trajectory_pair_gate(summary_json: str) -> str:
    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    return json.dumps(evaluate_magnetic_trajectory_pair(summary), indent=2, sort_keys=True)
