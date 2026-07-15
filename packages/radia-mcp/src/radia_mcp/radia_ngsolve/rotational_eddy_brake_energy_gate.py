"""Solver-neutral energy gate for a freely decelerating eddy-current brake."""

from __future__ import annotations

import math
from typing import Any


_UNITS = {
    "time": "s",
    "angular_velocity": "rad/s",
    "torque": "N*m",
    "power": "W",
    "inertia": "kg*m^2",
    "energy": "J",
    "density": "kg/m^3",
    "length": "m",
}


def _series(row: dict[str, Any], key: str) -> list[float]:
    value = row.get(key)
    if not isinstance(value, list):
        raise ValueError(f"{key} must be a list")
    result = [float(item) for item in value]
    if not result or not all(math.isfinite(item) for item in result):
        raise ValueError(f"{key} must contain finite values")
    return result


def _increasing(values: list[float]) -> bool:
    return len(values) >= 2 and all(
        right > left for left, right in zip(values, values[1:])
    )


def _integral(values: list[float], times: list[float]) -> float:
    return sum(
        0.5 * (left + right) * (t_right - t_left)
        for left, right, t_left, t_right in zip(
            values[:-1], values[1:], times[:-1], times[1:], strict=True
        )
    )


def _cumulative_integral(values: list[float], times: list[float]) -> list[float]:
    result = [0.0]
    for left, right, t_left, t_right in zip(
        values[:-1], values[1:], times[:-1], times[1:], strict=True
    ):
        result.append(result[-1] + 0.5 * (left + right) * (t_right - t_left))
    return result


def _span_error(left: list[float], right: list[float], reference: list[float]) -> float:
    span = max(reference) - min(reference)
    if len(left) != len(right) or span <= 0.0:
        return math.inf
    return max(abs(a - b) for a, b in zip(left, right, strict=True)) / span


def _parse_replay(row: Any) -> dict[str, Any]:
    if not isinstance(row, dict):
        raise ValueError("replay entries must be mappings")
    parsed = {
        key: _series(row, key)
        for key in (
            "time_s",
            "angular_velocity_rad_s",
            "braking_torque_nm",
            "joule_loss_w",
        )
    }
    parsed["solve_seconds"] = float(row.get("solve_seconds", math.nan))
    return parsed


def rotational_eddy_brake_energy_gate(
    summary: dict[str, Any],
    *,
    maximum_inertia_relative_error: float = 1.0e-5,
    maximum_angular_impulse_residual: float = 0.01,
    maximum_total_energy_residual: float = 0.01,
    maximum_replay_error_over_span: float = 1.0e-6,
    maximum_field_energy_adjacent_jump_fraction: float = 0.1,
    minimum_decay_fraction: float = 0.5,
) -> dict[str, Any]:
    """Gate free rotational braking with angular momentum and full energy storage.

    The energy identity includes magnetic field storage.  Torque times speed is
    reported only as a diagnostic because it is not generally equal to Joule
    loss when the field-energy rate is unavailable.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    contract = summary.get("contract")
    units = summary.get("units")
    disc = summary.get("disc")
    replays = summary.get("replays")
    energy_row = summary.get("energy_replay")
    timing = summary.get("timing_breakdown_s")
    if not all(isinstance(value, dict) for value in (contract, units, disc, energy_row)):
        raise ValueError("contract, units, disc, and energy_replay must be mappings")
    if not isinstance(replays, list) or len(replays) < 2:
        raise ValueError("at least two fresh replays are required")

    density = float(disc.get("density_kg_m3", math.nan))
    radius = float(disc.get("radius_m", math.nan))
    thickness = float(disc.get("thickness_m", math.nan))
    reported_inertia = float(summary.get("reported_inertia_kg_m2", math.nan))
    if not all(
        math.isfinite(value) and value > 0.0
        for value in (density, radius, thickness, reported_inertia)
    ):
        raise ValueError("disc dimensions, density, and inertia must be positive")
    analytic_inertia = 0.5 * density * math.pi * radius**4 * thickness
    inertia_error = abs(reported_inertia - analytic_inertia) / analytic_inertia

    parsed = [_parse_replay(row) for row in replays]
    all_cardinalities = True
    all_times_increase = True
    nonnegative_dissipation = True
    monotone_decay = True
    impulse_errors: list[float] = []
    decay_fractions: list[float] = []
    power_diagnostics: list[float] = []
    for replay in parsed:
        lengths = {len(value) for key, value in replay.items() if isinstance(value, list)}
        all_cardinalities &= len(lengths) == 1 and next(iter(lengths), 0) >= 20
        times = replay["time_s"]
        omega = replay["angular_velocity_rad_s"]
        torque = replay["braking_torque_nm"]
        joule = replay["joule_loss_w"]
        all_times_increase &= _increasing(times)
        nonnegative_dissipation &= min(torque) >= 0.0 and min(joule) >= 0.0
        monotone_decay &= omega[0] > 0.0 and all(
            right <= left + 1.0e-10 * max(abs(left), 1.0)
            for left, right in zip(omega, omega[1:])
        )
        angular_scale = reported_inertia * abs(omega[0] - omega[-1])
        impulse = _cumulative_integral(torque, times)
        residual = [
            reported_inertia * (value - omega[0]) + integrated
            for value, integrated in zip(omega, impulse, strict=True)
        ]
        impulse_errors.append(
            max(abs(value) for value in residual) / angular_scale
            if angular_scale > 0.0
            else math.inf
        )
        decay_fractions.append(1.0 - omega[-1] / omega[0])
        power_scale = max(joule)
        power_diagnostics.append(
            max(abs(t * w - q) for t, w, q in zip(torque, omega, joule, strict=True))
            / power_scale
            if power_scale > 0.0
            else math.inf
        )

    reference = parsed[0]
    replay_time_errors: list[float] = []
    replay_field_errors: list[float] = []
    for replay in parsed[1:]:
        replay_time_errors.append(
            max(
                (
                    abs(a - b)
                    for a, b in zip(
                        reference["time_s"], replay["time_s"], strict=True
                    )
                ),
                default=math.inf,
            )
            if len(reference["time_s"]) == len(replay["time_s"])
            else math.inf
        )
        replay_field_errors.extend(
            _span_error(reference[key], replay[key], reference[key])
            for key in (
                "angular_velocity_rad_s",
                "braking_torque_nm",
                "joule_loss_w",
            )
        )

    energy = _parse_replay(energy_row)
    field_time = _series(energy_row, "field_energy_time_s")
    magnetic_energy = _series(energy_row, "magnetic_energy_j")
    field_energy_scale = max(abs(value) for value in magnetic_energy)
    maximum_field_energy_adjacent_jump_fraction_observed = (
        max(
            abs(right - left)
            for left, right in zip(magnetic_energy, magnetic_energy[1:])
        )
        / field_energy_scale
        if field_energy_scale > 0.0
        else math.inf
    )
    energy_times = energy["time_s"]
    energy_omega = energy["angular_velocity_rad_s"]
    energy_joule = energy["joule_loss_w"]
    energy_cardinality = (
        len({len(value) for key, value in energy.items() if isinstance(value, list)}) == 1
        and len(field_time) == len(magnetic_energy) >= 2
    )
    field_time_alignment = (
        _increasing(field_time)
        and abs(field_time[0] - energy_times[0]) <= 1.0e-12
        and abs(field_time[-1] - energy_times[-1]) <= 1.0e-12
    )
    kinetic_drop = 0.5 * reported_inertia * (
        energy_omega[0] ** 2 - energy_omega[-1] ** 2
    )
    magnetic_drop = magnetic_energy[0] - magnetic_energy[-1]
    total_stored_drop = kinetic_drop + magnetic_drop
    joule_energy = _integral(energy_joule, energy_times)
    total_energy_error = (
        abs(total_stored_drop - joule_energy) / abs(total_stored_drop)
        if total_stored_drop > 0.0
        else math.inf
    )
    energy_replay_error = max(
        _span_error(energy[key], reference[key], reference[key])
        for key in (
            "angular_velocity_rad_s",
            "braking_torque_nm",
            "joule_loss_w",
        )
    )

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
        "si_units_explicit": units == _UNITS,
        "free_brake_contract_recorded": contract.get("body")
        == "uniform_solid_conducting_disc"
        and contract.get("inertia_reference") == "analytic_uniform_solid_disc"
        and contract.get("angular_momentum_balance")
        == "inertia_delta_angular_velocity_plus_integrated_braking_torque_equals_zero",
        "instantaneous_power_is_diagnostic_only": contract.get(
            "instantaneous_power_comparison"
        )
        == "diagnostic_only_when_field_energy_rate_is_not_sampled_on_the_probe_grid",
        "full_energy_storage_contract_recorded": contract.get("energy_balance")
        == "initial_kinetic_plus_magnetic_equals_final_kinetic_plus_magnetic_plus_joule",
        "analytic_disc_inertia_matches_reported": inertia_error
        <= float(maximum_inertia_relative_error),
        "replay_series_are_complete": all_cardinalities,
        "time_axes_strictly_increase": all_times_increase,
        "braking_torque_and_joule_loss_nonnegative": nonnegative_dissipation,
        "angular_velocity_monotonically_decays": monotone_decay
        and min(decay_fractions, default=-math.inf) >= float(minimum_decay_fraction),
        "angular_impulse_balance_closes": max(impulse_errors, default=math.inf)
        <= float(maximum_angular_impulse_residual),
        "fresh_replay_time_axes_match": max(replay_time_errors, default=math.inf)
        <= 1.0e-12,
        "fresh_replay_fields_match": max(replay_field_errors, default=math.inf)
        <= float(maximum_replay_error_over_span),
        "field_energy_history_is_present_and_aligned": energy_cardinality
        and field_time_alignment,
        "field_energy_history_is_nonnegative_and_has_no_isolated_jump": min(
            magnetic_energy
        )
        >= 0.0
        and maximum_field_energy_adjacent_jump_fraction_observed
        <= float(maximum_field_energy_adjacent_jump_fraction),
        "field_energy_run_replays_primary_history": energy_replay_error
        <= float(maximum_replay_error_over_span),
        "kinetic_magnetic_joule_energy_closes": total_energy_error
        <= float(maximum_total_energy_residual),
        "exactly_four_timing_stages": timing_ok,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "rotational_eddy_brake_energy_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "analytic_inertia_kg_m2": analytic_inertia,
            "reported_inertia_kg_m2": reported_inertia,
            "inertia_relative_error": inertia_error,
            "maximum_angular_impulse_residual_relative": max(
                impulse_errors, default=math.inf
            ),
            "maximum_replay_error_over_span": max(
                replay_field_errors + [energy_replay_error], default=math.inf
            ),
            "minimum_angular_velocity_decay_fraction": min(
                decay_fractions, default=-math.inf
            ),
            "kinetic_energy_drop_j": kinetic_drop,
            "magnetic_energy_drop_j": magnetic_drop,
            "maximum_field_energy_adjacent_jump_fraction": (
                maximum_field_energy_adjacent_jump_fraction_observed
            ),
            "integrated_joule_loss_j": joule_energy,
            "total_energy_residual_relative": total_energy_error,
            "mechanical_joule_relative_diagnostic": max(
                power_diagnostics, default=math.inf
            ),
        },
        "tolerances": {
            "maximum_inertia_relative_error": float(maximum_inertia_relative_error),
            "maximum_angular_impulse_residual": float(
                maximum_angular_impulse_residual
            ),
            "maximum_total_energy_residual": float(maximum_total_energy_residual),
            "maximum_replay_error_over_span": float(maximum_replay_error_over_span),
            "maximum_field_energy_adjacent_jump_fraction": float(
                maximum_field_energy_adjacent_jump_fraction
            ),
            "minimum_decay_fraction": float(minimum_decay_fraction),
        },
        "lesson": (
            "For a freely decelerating conductor, close angular impulse and the "
            "combined kinetic-plus-magnetic energy balance. Torque times speed is "
            "only diagnostic unless the magnetic-energy rate is represented."
        ),
    }
