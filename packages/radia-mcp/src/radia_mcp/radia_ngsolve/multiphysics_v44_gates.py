"""Public-safe COMSOL-derived identity gates for generalization v44."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _sha(value: object) -> bool:
    text = str(value or "")
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text.lower())


def _equal_numbers(identity: Mapping[str, object], fields: Sequence[str]) -> bool:
    try:
        return all(
            math.isclose(float(identity[field]), float(identity[f"result_{field}"]), rel_tol=1e-12, abs_tol=1e-15)
            for field in fields
        )
    except (KeyError, TypeError, ValueError):
        return False


def microwave_boundaryport_v44_ok(summary: Mapping[str, object]) -> bool:
    identity = summary.get(
        "microwave_boundaryport_sparameter_power_normalization_temperature_coupling_restart_owner_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("microwave_port_generation") or "")
    fields = (
        "frequency_hz", "port_normalization_ohm", "incident_power_w",
        "reflected_power_w", "transmitted_power_w", "absorbed_power_w",
        "s11_power_fraction", "s21_power_fraction", "thermal_coupling_power_w",
        "temperature_rise_k",
    )
    try:
        values = {field: float(identity[field]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(identity.get(field) == generation for field in (
            "port_generation", "power_generation", "thermal_generation",
            "restart_generation", "owner_generation", "result_generation",
        ))
        and all(math.isfinite(value) and value > 0.0 for value in values.values())
        and values["reflected_power_w"] < values["incident_power_w"]
        and values["transmitted_power_w"] < values["incident_power_w"]
        and math.isclose(values["absorbed_power_w"], values["incident_power_w"] - values["reflected_power_w"] - values["transmitted_power_w"], rel_tol=1e-12)
        and math.isclose(values["s11_power_fraction"], values["reflected_power_w"] / values["incident_power_w"], rel_tol=1e-12)
        and math.isclose(values["s21_power_fraction"], values["transmitted_power_w"] / values["incident_power_w"], rel_tol=1e-12)
        and math.isclose(values["thermal_coupling_power_w"], values["absorbed_power_w"], rel_tol=1e-12)
        and _equal_numbers(identity, fields)
        and identity.get("result_restart_checkpoint_id") == identity.get("restart_checkpoint_id")
        and bool(str(identity.get("mesh_owner") or ""))
        and identity.get("accepted_mesh_owner") == identity.get("mesh_owner")
        and _sha(identity.get("mesh_sha256"))
        and identity.get("accepted_mesh_sha256") == identity.get("mesh_sha256")
        and _sha(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def acoustic_poroelastic_v44_ok(summary: Mapping[str, object]) -> bool:
    identity = summary.get(
        "acoustics_poroelastic_impedance_phase_flux_energy_timewindow_dataset_owner_result_identity"
    )
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    generation = str(identity.get("acoustic_generation") or "")
    fields = (
        "frequency_hz", "impedance_magnitude_pa_s_per_m", "impedance_phase_deg",
        "pressure_flux_w", "mechanical_energy_j", "fluid_energy_j",
        "dissipated_power_w", "time_window_s", "energy_balance_residual_w",
    )
    try:
        values = {field: float(identity[field]) for field in fields}
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(generation)
        and all(identity.get(field) == generation for field in (
            "impedance_generation", "phase_generation", "flux_generation",
            "energy_generation", "dataset_generation", "owner_generation",
            "result_generation",
        ))
        and all(math.isfinite(value) for value in values.values())
        and values["frequency_hz"] > 0.0
        and values["impedance_magnitude_pa_s_per_m"] > 0.0
        and -180.0 <= values["impedance_phase_deg"] <= 180.0
        and values["pressure_flux_w"] >= 0.0
        and values["mechanical_energy_j"] >= 0.0
        and values["fluid_energy_j"] >= 0.0
        and values["dissipated_power_w"] >= 0.0
        and values["time_window_s"] > 0.0
        and abs(values["energy_balance_residual_w"]) <= 1e-12
        and math.isclose(values["pressure_flux_w"], values["dissipated_power_w"], rel_tol=1e-12)
        and _equal_numbers(identity, fields)
        and identity.get("result_dataset_tag") == identity.get("dataset_tag")
        and bool(str(identity.get("boundary_owner") or ""))
        and identity.get("accepted_boundary_owner") == identity.get("boundary_owner")
        and _sha(identity.get("boundary_sha256"))
        and identity.get("accepted_boundary_sha256") == identity.get("boundary_sha256")
        and _sha(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )
