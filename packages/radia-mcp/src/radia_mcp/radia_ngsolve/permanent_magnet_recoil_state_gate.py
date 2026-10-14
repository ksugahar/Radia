"""Solver-neutral replay gate for permanent-magnet demagnetization states."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _relative(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def permanent_magnet_recoil_state_gate(
    summary: Mapping[str, object],
    *,
    maximum_replay_relative_error: float = 1.0e-3,
    minimum_initial_axis_concentration: float = 100.0,
    maximum_open_axis_concentration: float = 2.0,
    minimum_recoil_axis_concentration: float = 10.0,
) -> dict[str, object]:
    """Gate nonlinear, open-circuit, and linear-recoil PM state evidence.

    The three states must share geometry, mesh, and observation points.  The
    gate intentionally uses only field ordering, spatial concentration, and
    replay agreement; it does not assume a particular solver or material name.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    tolerance = float(maximum_replay_relative_error)
    concentration_limits = (
        float(minimum_initial_axis_concentration),
        float(maximum_open_axis_concentration),
        float(minimum_recoil_axis_concentration),
    )
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("maximum_replay_relative_error must be finite and nonnegative")
    if not all(math.isfinite(value) and value > 0.0 for value in concentration_limits):
        raise ValueError("concentration thresholds must be finite and positive")

    contract = summary.get("model_contract")
    states = summary.get("states")
    units = summary.get("units")
    if not isinstance(contract, Mapping) or not isinstance(states, Mapping):
        raise ValueError("model_contract and states must be mappings")
    if not isinstance(units, Mapping):
        raise ValueError("units must be a mapping")

    expected_roles = {
        "initial": "nonlinear_in_circuit",
        "out_of_circuit": "nonlinear_open_circuit",
        "recoil_return": "linear_recoil_in_circuit",
    }
    parsed: dict[str, dict[str, float]] = {}
    for state_name in expected_roles:
        state = states.get(state_name)
        if not isinstance(state, Mapping):
            raise ValueError(f"states.{state_name} must be a mapping")
        parsed[state_name] = {
            "on_axis": float(state.get("on_axis", math.nan)),
            "off_axis": float(state.get("off_axis", math.nan)),
        }
    state_values = [value for state in parsed.values() for value in state.values()]

    stored = summary.get("stored_reference")
    fresh = summary.get("fresh_replay")
    if not isinstance(stored, Sequence) or isinstance(stored, (str, bytes)) or len(stored) != 6:
        raise ValueError("stored_reference must contain exactly six values")
    if not isinstance(fresh, Sequence) or isinstance(fresh, (str, bytes)) or len(fresh) != 6:
        raise ValueError("fresh_replay must contain exactly six values")
    stored_values = [float(value) for value in stored]
    fresh_values = [float(value) for value in fresh]
    replay_errors = [_relative(left, right) for left, right in zip(stored_values, fresh_values)]

    initial = parsed["initial"]
    opened = parsed["out_of_circuit"]
    recoil = parsed["recoil_return"]
    concentrations = {
        name: value["on_axis"] / max(value["off_axis"], 1.0e-300)
        for name, value in parsed.items()
    }
    recovery = {
        "recoil_to_open_on_axis": recoil["on_axis"] / max(opened["on_axis"], 1.0e-300),
        "recoil_to_initial_on_axis": recoil["on_axis"] / max(initial["on_axis"], 1.0e-300),
    }
    checks = {
        "state_roles_are_explicit": all(contract.get(key) == role for key, role in expected_roles.items()),
        "geometry_mesh_and_points_shared": all(
            contract.get(key) is True
            for key in ("same_geometry", "same_mesh", "same_observation_points")
        ),
        "magnetic_flux_density_unit_is_tesla": str(units.get("magnetic_flux_density", "")).strip().lower()
        in {"t", "tesla"},
        "all_state_fields_finite_positive": all(
            math.isfinite(value) and value > 0.0 for value in state_values
        ),
        "stored_and_fresh_values_finite_positive": all(
            math.isfinite(value) and value > 0.0
            for value in stored_values + fresh_values
        ),
        "saved_and_fresh_replay_close": max(replay_errors) <= tolerance,
        "on_axis_state_order_initial_recoil_open": initial["on_axis"] > recoil["on_axis"] > opened["on_axis"],
        "initial_field_is_axis_concentrated": concentrations["initial"] >= concentration_limits[0],
        "open_field_is_spatially_spread": 1.0 / concentration_limits[1]
        <= concentrations["out_of_circuit"]
        <= concentration_limits[1],
        "recoil_field_is_axis_concentrated": concentrations["recoil_return"] >= concentration_limits[2],
        "recoil_is_partial_not_full_recovery": recovery["recoil_to_open_on_axis"] > 1.0
        and 0.0 < recovery["recoil_to_initial_on_axis"] < 1.0,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "permanent_magnet_recoil_state_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "axis_concentration": concentrations,
            "recovery_ratios": recovery,
            "replay_relative_errors": replay_errors,
            "maximum_replay_relative_error": max(replay_errors),
        },
        "tolerances": {
            "maximum_replay_relative_error": tolerance,
            "minimum_initial_axis_concentration": concentration_limits[0],
            "maximum_open_axis_concentration": concentration_limits[1],
            "minimum_recoil_axis_concentration": concentration_limits[2],
        },
        "lesson": (
            "Demagnetization evidence needs distinct nonlinear in-circuit, open-circuit, "
            "and recoil-return states. Full recovery after the open state is not implied: "
            "accept only a reproducible partial return with unchanged probes and mesh."
        ),
    }
