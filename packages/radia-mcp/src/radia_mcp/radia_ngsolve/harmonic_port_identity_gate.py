"""Solver-neutral identities for a one-port harmonic MQS result."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _complex(value: object, name: str) -> complex:
    if isinstance(value, Mapping):
        result = complex(float(value.get("real", 0.0)), float(value.get("imag", 0.0)))
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) != 2:
            raise ValueError(f"{name} must contain real and imaginary components")
        result = complex(float(value[0]), float(value[1]))
    else:
        result = complex(value)
    if not math.isfinite(result.real) or not math.isfinite(result.imag):
        raise ValueError(f"{name} must be finite")
    return result


def _relative_error(left: complex | float, right: complex | float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def _profile(value: object, name: str) -> list[tuple[float, float]]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{name} must be a nonempty list")
    result = []
    for index, row in enumerate(value):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or len(row) != 2:
            raise ValueError(f"{name}[{index}] must contain coordinate and magnitude")
        pair = (float(row[0]), float(row[1]))
        if not all(math.isfinite(component) for component in pair):
            raise ValueError(f"{name}[{index}] must be finite")
        result.append(pair)
    return result


def harmonic_current_port_power_energy_identity_gate(
    summary: Mapping[str, object],
    *,
    maximum_identity_relative_error: float = 1.0e-9,
    maximum_cross_run_relative_error: float = 1.0e-9,
) -> dict[str, object]:
    """Gate peak-phasor port, loss, energy, flux, and profile identities."""
    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    identity_tolerance = float(maximum_identity_relative_error)
    cross_tolerance = float(maximum_cross_run_relative_error)
    if not all(
        math.isfinite(value) and value >= 0.0
        for value in (identity_tolerance, cross_tolerance)
    ):
        raise ValueError("tolerances must be finite and nonnegative")

    frequency = float(summary.get("frequency_hz", 0.0))
    if not math.isfinite(frequency) or frequency <= 0.0:
        raise ValueError("frequency_hz must be finite and positive")
    omega = 2.0 * math.pi * frequency
    records = summary.get("runs")
    if not isinstance(records, list) or len(records) != 2:
        raise ValueError("runs must contain exactly two records")

    parsed = []
    identity_errors = {}
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise ValueError(f"runs[{index}] must be a mapping")
        label = str(record.get("label") or f"run_{index}")
        current = _complex(record.get("current_a"), f"runs[{index}].current_a")
        if abs(current) == 0.0:
            raise ValueError(f"runs[{index}].current_a must be nonzero")
        voltage = _complex(record.get("voltage_v"), f"runs[{index}].voltage_v")
        impedance = _complex(record.get("impedance_ohm"), f"runs[{index}].impedance_ohm")
        power = _complex(record.get("complex_power_w"), f"runs[{index}].complex_power_w")
        flux = _complex(record.get("flux_linkage_wb"), f"runs[{index}].flux_linkage_wb")
        energy = float(record.get("magnetic_energy_j", math.nan))
        loss = float(record.get("loss_w", math.nan))
        conductor_loss = float(record.get("conductor_loss_w", math.nan))
        if not all(math.isfinite(value) and value >= 0.0 for value in (energy, loss, conductor_loss)):
            raise ValueError(f"runs[{index}] energy and losses must be finite and nonnegative")
        profile = _profile(record.get("current_density_profile"), f"runs[{index}].current_density_profile")

        errors = {
            "v_equals_z_i": _relative_error(voltage, impedance * current),
            "peak_complex_power": _relative_error(power, 0.5 * voltage * current.conjugate()),
            "real_power_equals_loss": _relative_error(power.real, loss),
            "loss_is_conductor_owned": _relative_error(loss, conductor_loss),
            "reactive_power_energy": _relative_error(energy, power.imag / (2.0 * omega)),
            "reactance_flux_linkage": _relative_error(
                impedance.imag, omega * (flux / current).real
            ),
        }
        identity_errors[label] = errors
        parsed.append(
            {
                "label": label,
                "current": current,
                "voltage": voltage,
                "impedance": impedance,
                "power": power,
                "flux": flux,
                "energy": energy,
                "loss": loss,
                "profile": profile,
            }
        )

    left, right = parsed
    profile_grid_matches = len(left["profile"]) == len(right["profile"]) and all(
        _relative_error(a[0], b[0]) <= cross_tolerance
        for a, b in zip(left["profile"], right["profile"])
    )
    profile_error = math.inf
    if profile_grid_matches:
        profile_error = max(
            _relative_error(a[1] / abs(left["current"]), b[1] / abs(right["current"]))
            for a, b in zip(left["profile"], right["profile"])
        )
    cross_errors = {
        "impedance": _relative_error(left["impedance"], right["impedance"]),
        "flux_per_ampere": _relative_error(
            left["flux"] / left["current"], right["flux"] / right["current"]
        ),
        "loss_per_ampere_squared": _relative_error(
            left["loss"] / abs(left["current"]) ** 2,
            right["loss"] / abs(right["current"]) ** 2,
        ),
        "energy_per_ampere_squared": _relative_error(
            left["energy"] / abs(left["current"]) ** 2,
            right["energy"] / abs(right["current"]) ** 2,
        ),
        "current_density_per_ampere": profile_error,
    }
    checks = {
        "peak_phasor_convention_explicit": summary.get("amplitude_convention") == "peak_phasor",
        "all_port_identities_close": max(
            error for errors in identity_errors.values() for error in errors.values()
        )
        <= identity_tolerance,
        "profile_grids_match": profile_grid_matches,
        "normalized_cross_run_observables_close": max(cross_errors.values())
        <= cross_tolerance,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "harmonic_current_port_power_energy_identity_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "identity_relative_errors": identity_errors,
            "cross_run_relative_errors": cross_errors,
            "profile_sample_count": len(left["profile"]),
        },
        "tolerances": {
            "maximum_identity_relative_error": identity_tolerance,
            "maximum_cross_run_relative_error": cross_tolerance,
        },
        "lesson": (
            "For peak phasors, S=0.5*V*conj(I), loss=Re(S), and stored magnetic "
            "energy=Im(S)/(2*omega). Excitation semantics and normalization must come "
            "from the result tree, never from a project filename."
        ),
    }
