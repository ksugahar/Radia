"""Solver-neutral cylindrical conductor skin-effect and Bessel gate."""

from __future__ import annotations

import math
from typing import Any


MU0 = 4.0 * math.pi * 1.0e-7


def _pair(value: object, name: str) -> complex:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{name} must contain [real, imag]")
    parsed = complex(float(value[0]), float(value[1]))
    if not math.isfinite(parsed.real) or not math.isfinite(parsed.imag):
        raise ValueError(f"{name} must be finite")
    return parsed


def _relative(left: complex | float, right: complex | float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-30)


def _bessel_j0(value: complex) -> complex:
    total = 1.0 + 0.0j
    term = total
    factor = -(value * value) / 4.0
    for order in range(1, 256):
        term *= factor / (order * order)
        total += term
        if abs(term) <= 1.0e-15 * max(abs(total), 1.0):
            return total
    raise ValueError("J0 series did not converge")


def _bessel_j1(value: complex) -> complex:
    total = value / 2.0
    term = total
    factor = -(value * value) / 4.0
    for order in range(1, 256):
        term *= factor / (order * (order + 1))
        total += term
        if abs(term) <= 1.0e-15 * max(abs(total), 1.0):
            return total
    raise ValueError("J1 series did not converge")


def cylindrical_conductor_skin_bessel_gate(
    summary: dict[str, Any],
    *,
    identity_rtol: float = 1.0e-8,
    dc_resistance_rtol: float = 1.0e-5,
    ac_resistance_rtol: float = 5.0e-4,
    profile_l2_rtol: float = 0.02,
    replay_rtol: float = 1.0e-10,
) -> dict[str, Any]:
    """Validate harmonic cylindrical skin effect against exact Bessel structure."""

    if not isinstance(summary, dict):
        raise ValueError("summary must be an object")
    model = summary.get("model_contract")
    replays = summary.get("replays")
    if not isinstance(model, dict):
        raise ValueError("model_contract must be an object")
    if not isinstance(replays, list) or len(replays) != 2:
        raise ValueError("replays must contain exactly two entries")
    frequency = float(model.get("frequency_hz", math.nan))
    radius = float(model.get("radius_m", math.nan))
    length = float(model.get("length_m", math.nan))
    conductivity = float(model.get("conductivity_s_per_m", math.nan))
    if not all(math.isfinite(value) and value > 0.0 for value in (frequency, radius, length, conductivity)):
        raise ValueError("frequency, radius, length, and conductivity must be finite and positive")
    omega = 2.0 * math.pi * frequency
    skin_depth = math.sqrt(2.0 / (omega * MU0 * conductivity))
    kappa = (1.0 + 1.0j) / skin_depth
    analytic_dc_resistance = length / (conductivity * math.pi * radius * radius)
    analytic_internal = (
        kappa
        * _bessel_j0(kappa * radius)
        / (2.0 * math.pi * radius * conductivity * _bessel_j1(kappa * radius))
        * length
    )

    metrics_per_replay = []
    normalized_replays = []
    for replay in replays:
        if not isinstance(replay, dict):
            raise ValueError("replay entries must be objects")
        current = _pair(replay.get("current_a"), "current_a")
        flux = _pair(replay.get("flux_linkage_wb"), "flux_linkage_wb")
        impedance = _pair(replay.get("impedance_ohm"), "impedance_ohm")
        power = _pair(replay.get("power_va"), "power_va")
        voltage = _pair(replay.get("voltage_v"), "voltage_v")
        induced = _pair(replay.get("induced_voltage_v"), "induced_voltage_v")
        energy = float(replay.get("energy_j", math.nan))
        loss = float(replay.get("loss_w", math.nan))
        residual = float(replay.get("final_log10_relative_residual", math.nan))
        radii = replay.get("profile_radii_m")
        density_rows = replay.get("current_density_a_per_m2")
        if not isinstance(radii, list) or not isinstance(density_rows, list) or len(radii) != len(density_rows) or len(radii) < 10:
            raise ValueError("profile radii and current density must have equal length >= 10")
        parsed_radii = [float(value) for value in radii]
        density = [_pair(value, "current_density_a_per_m2") for value in density_rows]
        if not all(math.isfinite(value) and value >= 0.0 for value in parsed_radii):
            raise ValueError("profile radii must be finite and nonnegative")
        if not all(right > left for left, right in zip(parsed_radii, parsed_radii[1:])) or parsed_radii[-1] > radius:
            raise ValueError("profile radii must be increasing and inside the conductor")
        shape = [_bessel_j0(kappa * value) for value in parsed_radii]
        denominator = sum(abs(value) ** 2 for value in shape)
        scale = sum(expected.conjugate() * actual for expected, actual in zip(shape, density, strict=True)) / denominator
        profile_l2 = math.sqrt(
            sum(abs(scale * expected - actual) ** 2 for expected, actual in zip(shape, density, strict=True))
            / sum(abs(actual) ** 2 for actual in density)
        )
        extracted_dc = ((voltage - induced) / current).real
        row_metrics = {
            "voltage_impedance_current_relative_error": _relative(voltage, impedance * current),
            "complex_power_relative_error": _relative(power, 0.5 * voltage * current.conjugate()),
            "real_power_loss_relative_error": _relative(power.real, loss),
            "reactive_power_energy_relative_error": _relative(power.imag, 2.0 * omega * energy),
            "faraday_flux_voltage_relative_error": _relative(induced, 1.0j * omega * flux),
            "dc_resistance_relative_error": _relative(extracted_dc, analytic_dc_resistance),
            "ac_resistance_relative_error": _relative(impedance.real, analytic_internal.real),
            "bessel_profile_relative_l2_error": profile_l2,
            "surface_to_center_current_density_ratio": abs(density[-1]) / max(abs(density[0]), 1.0e-30),
            "final_log10_relative_residual": residual,
        }
        metrics_per_replay.append(row_metrics)
        normalized_replays.append(
            {
                "scalars": [current, flux, impedance, power, voltage, induced, complex(energy), complex(loss)],
                "radii": parsed_radii,
                "density": density,
            }
        )

    replay_errors = []
    left, right = normalized_replays
    for left_value, right_value in zip(left["scalars"], right["scalars"], strict=True):
        replay_errors.append(_relative(left_value, right_value))
    for left_value, right_value in zip(left["radii"], right["radii"], strict=True):
        replay_errors.append(_relative(left_value, right_value))
    for left_value, right_value in zip(left["density"], right["density"], strict=True):
        replay_errors.append(_relative(left_value, right_value))

    timing = summary.get("timing_breakdown_s")
    timing_ok = False
    if isinstance(timing, dict) and len(timing) == 4:
        try:
            timing_ok = all(math.isfinite(float(value)) and float(value) >= 0.0 for value in timing.values())
        except (TypeError, ValueError):
            timing_ok = False

    maximums = {
        key: max(row[key] for row in metrics_per_replay)
        for key in (
            "voltage_impedance_current_relative_error",
            "complex_power_relative_error",
            "real_power_loss_relative_error",
            "reactive_power_energy_relative_error",
            "faraday_flux_voltage_relative_error",
            "dc_resistance_relative_error",
            "ac_resistance_relative_error",
            "bessel_profile_relative_l2_error",
        )
    }
    checks = {
        "port_and_energy_identities_close": max(
            maximums["voltage_impedance_current_relative_error"],
            maximums["complex_power_relative_error"],
            maximums["real_power_loss_relative_error"],
            maximums["reactive_power_energy_relative_error"],
            maximums["faraday_flux_voltage_relative_error"],
        )
        <= identity_rtol,
        "dc_resistance_matches_cylinder": maximums["dc_resistance_relative_error"] <= dc_resistance_rtol,
        "ac_resistance_matches_bessel_internal_impedance": maximums["ac_resistance_relative_error"] <= ac_resistance_rtol,
        "current_density_matches_bessel_profile": maximums["bessel_profile_relative_l2_error"] <= profile_l2_rtol,
        "skin_crowding_is_resolved": min(row["surface_to_center_current_density_ratio"] for row in metrics_per_replay) >= 10.0,
        "independent_replays_are_deterministic": max(replay_errors, default=math.inf) <= replay_rtol,
        "both_linear_residuals_are_below_1e_10": all(row["final_log10_relative_residual"] <= -10.0 for row in metrics_per_replay),
        "exactly_four_timing_stages": timing_ok,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cylindrical_conductor_skin_bessel_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            **maximums,
            "maximum_replay_relative_error": max(replay_errors),
            "minimum_surface_to_center_current_density_ratio": min(row["surface_to_center_current_density_ratio"] for row in metrics_per_replay),
            "skin_depth_m": skin_depth,
            "analytic_dc_resistance_ohm": analytic_dc_resistance,
            "analytic_ac_resistance_ohm": analytic_internal.real,
        },
        "lesson": (
            "For a round conductor, close port and energy identities, then compare "
            "Rdc to length/(sigma*pi*a^2), Rac to the exact complex Bessel internal "
            "impedance, and the full Jz(r) profile to J0((1+j)r/delta)."
        ),
    }
