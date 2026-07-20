"""Wave-energy balance and eigenmode-Q identity checks for v52."""

from __future__ import annotations

import math
from collections.abc import Mapping

from .wave_port_identity_v53 import validate_public_v53_identity


ENERGY_BALANCE = "energy_balance_incident_reflected_transmitted_absorbed_dissipated_run_owner_identity"
EIGENMODE_Q = "eigenmode_q_stored_energy_boundary_loss_normalization_mode_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _energy_balance_ok(row: Mapping[str, object]) -> bool:
    incident = row.get("incident_power_w")
    reflected = row.get("reflected_power_w")
    transmitted = row.get("transmitted_power_w")
    absorbed = row.get("absorbed_power_w")
    dissipated = row.get("dissipated_power_w")
    powers_ok = all(_finite(value) and float(value) >= 0.0 for value in (incident, reflected, transmitted, absorbed))
    dissipated_ok = (
        isinstance(dissipated, Mapping)
        and bool(dissipated)
        and all(isinstance(name, str) and bool(name) and _finite(value) and float(value) >= 0.0 for name, value in dissipated.items())
    )
    return (
        _generation(row, ("incident_generation", "scattered_generation", "absorbed_generation", "dissipated_generation", "owner_generation", "result_generation"))
        and powers_ok
        and float(incident) > 0.0
        and math.isclose(float(incident), float(reflected) + float(transmitted) + float(absorbed), rel_tol=1.0e-10, abs_tol=1.0e-12)
        and dissipated_ok
        and math.isclose(float(absorbed), sum(float(value) for value in dissipated.values()), rel_tol=1.0e-10, abs_tol=1.0e-12)
        and row.get("result_incident_power_w") == incident
        and row.get("result_reflected_power_w") == reflected
        and row.get("result_transmitted_power_w") == transmitted
        and row.get("result_absorbed_power_w") == absorbed
        and row.get("result_dissipated_power_w") == dissipated
        and str(row.get("run_owner") or "").startswith("run:")
        and row.get("result_run_owner") == row.get("run_owner")
        and _result(row)
    )


def _eigenmode_q_ok(row: Mapping[str, object]) -> bool:
    frequency = row.get("frequency_hz")
    stored_energy = row.get("stored_energy_j")
    boundary_loss = row.get("boundary_loss_w")
    volume_loss = row.get("volume_loss_w")
    q_factor = row.get("q_factor")
    numeric_ok = all(_finite(value) for value in (frequency, stored_energy, boundary_loss, volume_loss, q_factor))
    total_loss = float(boundary_loss) + float(volume_loss) if numeric_ok else 0.0
    expected_q = 2.0 * math.pi * float(frequency) * float(stored_energy) / total_loss if total_loss > 0.0 else None
    return (
        _generation(row, ("frequency_generation", "energy_generation", "loss_generation", "normalization_generation", "owner_generation", "result_generation"))
        and numeric_ok
        and float(frequency) > 0.0
        and float(stored_energy) > 0.0
        and float(boundary_loss) >= 0.0
        and float(volume_loss) >= 0.0
        and expected_q is not None
        and math.isclose(float(q_factor), expected_q, rel_tol=1.0e-10, abs_tol=1.0e-10)
        and row.get("result_frequency_hz") == frequency
        and row.get("result_stored_energy_j") == stored_energy
        and row.get("result_boundary_loss_w") == boundary_loss
        and row.get("result_volume_loss_w") == volume_loss
        and row.get("result_q_factor") == q_factor
        and row.get("normalization") == "physical_stored_energy"
        and row.get("result_normalization") == row.get("normalization")
        and str(row.get("mode_owner") or "").startswith("mode:")
        and row.get("result_mode_owner") == row.get("mode_owner")
        and _result(row)
    )


def validate_public_v52_identity(payload: object) -> dict[str, bool]:
    if not isinstance(payload, Mapping):
        return {}
    rows = [row for row in (payload.get("runs") or []) if isinstance(row, Mapping)]
    checks = validate_public_v53_identity(payload)
    balances = [row[ENERGY_BALANCE] for row in rows if ENERGY_BALANCE in row]
    q_factors = [row[EIGENMODE_Q] for row in rows if EIGENMODE_Q in row]
    if balances:
        checks["wave_v52_power_conservation_dissipation_run_owner"] = len(balances) == len(rows) and all(isinstance(row, Mapping) and _energy_balance_ok(row) for row in balances)
    if q_factors:
        checks["wave_v52_eigenmode_q_energy_loss_normalization_owner"] = len(q_factors) == len(rows) and all(isinstance(row, Mapping) and _eigenmode_q_ok(row) for row in q_factors)
    return checks
