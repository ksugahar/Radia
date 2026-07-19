"""Solver-neutral far-field and EM-to-thermal semantic identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


FARFIELD = "farfield_polarization_basis_normalization_radiated_power_angular_grid_owner_identity"
THERMAL = "em_thermal_loss_mapping_mesh_interpolation_time_average_frequency_owner_identity"


def _sha(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _digest(row: Mapping[str, object]) -> bool:
    return _sha(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_sequence(value: object, *, minimum: int = 1) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) >= minimum
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _strictly_increasing(value: object) -> bool:
    return _finite_sequence(value) and all(float(left) < float(right) for left, right in zip(value, value[1:]))


def _field_samples(value: object, count: int) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) == count
        and all(_finite_sequence(sample, minimum=2) and len(sample) == 2 for sample in value)
    )


def _farfield_ok(row: Mapping[str, object]) -> bool:
    theta = row.get("theta_deg")
    phi = row.get("phi_deg")
    samples = row.get("field_samples")
    power = row.get("radiated_power_w")
    return (
        _generation(row, ("polarization_generation", "normalization_generation", "power_generation", "angular_grid_generation", "monitor_generation", "result_generation"))
        and row.get("polarization_basis") == row.get("result_polarization_basis") == "spherical_theta_phi"
        and row.get("normalization") == row.get("result_normalization") == "accepted_radiated_power"
        and isinstance(power, (int, float))
        and math.isfinite(float(power))
        and float(power) > 0.0
        and row.get("result_radiated_power_w") == power
        and _strictly_increasing(theta)
        and all(0.0 <= float(value) <= 180.0 for value in theta)
        and row.get("result_theta_deg") == theta
        and _strictly_increasing(phi)
        and all(0.0 <= float(value) < 360.0 for value in phi)
        and row.get("result_phi_deg") == phi
        and _field_samples(samples, len(theta))
        and row.get("result_field_samples") == samples
        and str(row.get("monitor_owner") or "").startswith("monitor:")
        and row.get("result_monitor_owner") == row.get("monitor_owner")
        and _digest(row)
    )


def _thermal_ok(row: Mapping[str, object]) -> bool:
    components = row.get("loss_component_ids")
    losses = row.get("loss_w")
    frequency = row.get("frequency_hz")
    return (
        _generation(row, ("loss_mapping_generation", "source_mesh_generation", "target_mesh_generation", "interpolation_generation", "time_average_generation", "frequency_generation", "task_generation", "result_generation"))
        and isinstance(components, list)
        and bool(components)
        and len(components) == len(set(components))
        and all(isinstance(component, str) and component for component in components)
        and row.get("mapped_loss_component_ids") == components
        and _finite_sequence(losses, minimum=len(components))
        and len(losses) == len(components)
        and all(float(loss) >= 0.0 for loss in losses)
        and row.get("mapped_loss_w") == losses
        and _sha(row.get("source_mesh_sha256"))
        and row.get("mapped_source_mesh_sha256") == row.get("source_mesh_sha256")
        and _sha(row.get("target_mesh_sha256"))
        and row.get("mapped_target_mesh_sha256") == row.get("target_mesh_sha256")
        and row.get("interpolation_method") == row.get("mapped_interpolation_method") == "conservative_nodal"
        and row.get("time_average") == row.get("mapped_time_average") == "cycle_average"
        and isinstance(frequency, (int, float))
        and math.isfinite(float(frequency))
        and float(frequency) > 0.0
        and row.get("mapped_frequency_hz") == frequency
        and str(row.get("task_owner") or "").startswith("task:")
        and row.get("mapped_task_owner") == row.get("task_owner")
        and _digest(row)
    )


def validate_public_v48_identity(payload: object) -> dict[str, bool]:
    if not isinstance(payload, Mapping):
        return {}
    rows = [row for row in (payload.get("runs") or []) if isinstance(row, Mapping)]
    checks: dict[str, bool] = {}
    farfields = [row[FARFIELD] for row in rows if FARFIELD in row]
    thermal = [row[THERMAL] for row in rows if THERMAL in row]
    if farfields:
        checks["network_v48_farfield_polarization_power_grid_owner"] = len(farfields) == len(rows) and all(isinstance(row, Mapping) and _farfield_ok(row) for row in farfields)
    if thermal:
        checks["network_v48_em_thermal_mapping_mesh_average_owner"] = len(thermal) == len(rows) and all(isinstance(row, Mapping) and _thermal_ok(row) for row in thermal)
    return checks
