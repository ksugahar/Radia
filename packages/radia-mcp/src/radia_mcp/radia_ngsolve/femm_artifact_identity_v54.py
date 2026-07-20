"""Harmonic-power and axisymmetric-force artifact identity checks for v54."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .femm_artifact_identity_v55 import validate_public_identity as validate_public_v55_identity


POWER = "harmonic_complexpower_active_reactive_loss_frequency_circuit_owner_identity"
FORCE = "axisymmetric_weightedstress_force_radius_measure_selection_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _power_ok(row: Mapping[str, object]) -> bool:
    complex_power = row.get("complex_power_va")
    active = row.get("active_power_w")
    reactive = row.get("reactive_power_var")
    losses = row.get("loss_components_w")
    complex_ok = isinstance(complex_power, Mapping) and set(complex_power) == {"real", "imag"} and all(_finite(value) for value in complex_power.values())
    losses_ok = isinstance(losses, Mapping) and bool(losses) and all(isinstance(name, str) and name and _finite(value) and float(value) >= 0.0 for name, value in losses.items())
    return (
        _generations(row, "power_generation", "loss_generation", "frequency_generation", "circuit_generation", "owner_generation", "result_generation")
        and complex_ok
        and _finite(active)
        and float(active) >= 0.0
        and math.isclose(float(complex_power["real"]), float(active), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and _finite(reactive)
        and math.isclose(float(complex_power["imag"]), float(reactive), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and losses_ok
        and math.isclose(sum(float(value) for value in losses.values()), float(active), rel_tol=1.0e-10, abs_tol=1.0e-12)
        and row.get("result_complex_power_va") == complex_power
        and row.get("result_active_power_w") == active
        and row.get("result_reactive_power_var") == reactive
        and row.get("result_loss_components_w") == losses
        and _finite(row.get("frequency_hz"))
        and float(row["frequency_hz"]) > 0.0
        and row.get("result_frequency_hz") == row.get("frequency_hz")
        and str(row.get("circuit_id") or "").startswith("circuit:")
        and row.get("result_circuit_id") == row.get("circuit_id")
        and str(row.get("circuit_owner") or "").startswith("circuit-owner:")
        and row.get("result_circuit_owner") == row.get("circuit_owner")
        and _result(row)
    )


def _force_ok(row: Mapping[str, object]) -> bool:
    selection = row.get("integration_selection")
    direction = row.get("force_direction_rz")
    selection_ok = (
        isinstance(selection, Sequence)
        and not isinstance(selection, (str, bytes))
        and bool(selection)
        and all(isinstance(item, str) and item.startswith(("block:", "interface:")) for item in selection)
        and len(selection) == len(set(selection))
    )
    direction_ok = (
        isinstance(direction, Sequence)
        and not isinstance(direction, (str, bytes))
        and len(direction) == 2
        and all(_finite(value) for value in direction)
        and math.isclose(sum(float(value) ** 2 for value in direction), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
    )
    return (
        _generations(row, "radius_generation", "stress_generation", "selection_generation", "direction_generation", "owner_generation", "result_generation")
        and row.get("radius_weighting") == "2*pi*r"
        and row.get("result_radius_weighting") == row.get("radius_weighting")
        and row.get("stress_measure") == "weighted_stress_tensor"
        and row.get("result_stress_measure") == row.get("stress_measure")
        and selection_ok
        and row.get("result_integration_selection") == selection
        and direction_ok
        and row.get("result_force_direction_rz") == direction
        and _finite(row.get("force_n"))
        and row.get("result_force_n") == row.get("force_n")
        and str(row.get("mesh_owner") or "").startswith("mesh:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = validate_public_v55_identity(identity)
    power = identity.get(POWER)
    force = identity.get(FORCE)
    if power is not None:
        checks["v54_harmonic_complex_power_loss_frequency_circuit_owner"] = isinstance(power, Mapping) and _power_ok(power)
    if force is not None:
        checks["v54_axisymmetric_weighted_stress_selection_direction_owner"] = isinstance(force, Mapping) and _force_ok(force)
    return checks
