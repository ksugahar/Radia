"""Acoustic-modal and induction-heating artifact identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


MODAL = "acoustic_modalparticipation_response_normalization_energy_owner_identity"
INDUCTION = "inductionheating_input_joule_thermal_time_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _number(value: object, *, positive: bool = False, nonnegative: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    number = float(value)
    return math.isfinite(number) and (not positive or number > 0.0) and (not nonnegative or number >= 0.0)


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-12)


def _generation(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result_identity(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _modal_ok(row: Mapping[str, object]) -> bool:
    modes = row.get("modal_rows")
    response = row.get("frequency_response")
    if not isinstance(modes, Sequence) or isinstance(modes, (str, bytes)) or len(modes) < 2:
        return False
    mode_names: list[str] = []
    frequencies: list[float] = []
    participation: list[float] = []
    energies: list[float] = []
    for mode in modes:
        if not isinstance(mode, Mapping) or set(mode) != {"mode", "frequency_hz", "participation", "normalized_energy"}:
            return False
        name = str(mode.get("mode") or "")
        if not name.startswith("mode:") or not _number(mode.get("frequency_hz"), positive=True) or not _number(mode.get("participation"), nonnegative=True) or not _number(mode.get("normalized_energy"), nonnegative=True):
            return False
        mode_names.append(name)
        frequencies.append(float(mode["frequency_hz"]))
        participation.append(float(mode["participation"]))
        energies.append(float(mode["normalized_energy"]))
    response_ok = isinstance(response, Sequence) and not isinstance(response, (str, bytes)) and len(response) == len(modes)
    response_frequencies: list[float] = []
    if response_ok:
        for record in response:
            if not isinstance(record, Mapping) or set(record) != {"frequency_hz", "pressure_pa"} or not _number(record.get("frequency_hz"), positive=True) or not _number(record.get("pressure_pa")):
                response_ok = False
                break
            response_frequencies.append(float(record["frequency_hz"]))
    total_energy = row.get("total_normalized_energy")
    return (
        _generation(row, "participation_generation", "response_generation", "normalization_generation", "energy_generation", "owner_generation", "result_generation")
        and len(mode_names) == len(set(mode_names)) and all(left < right for left, right in zip(frequencies, frequencies[1:]))
        and _close(sum(participation), 1.0) and _close(sum(energies), 1.0)
        and response_ok and response_frequencies == frequencies
        and row.get("normalization") == "unit_total_modal_energy"
        and _number(total_energy, positive=True) and _close(total_energy, sum(energies))
        and all(row.get("result_" + field) == row.get(field) for field in ("modal_rows", "frequency_response", "normalization", "total_normalized_energy"))
        and str(row.get("solution_owner") or "").startswith("solution:") and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity(row)
    )


def _induction_ok(row: Mapping[str, object]) -> bool:
    input_energy = row.get("coil_input_energy_j")
    joule = row.get("joule_heat_energy_j")
    stored = row.get("stored_thermal_energy_j")
    heat_loss = row.get("boundary_heat_loss_j")
    times = row.get("time_s")
    temperatures = row.get("average_temperature_k")
    history_ok = (
        isinstance(times, Sequence) and not isinstance(times, (str, bytes)) and len(times) >= 2
        and isinstance(temperatures, Sequence) and not isinstance(temperatures, (str, bytes)) and len(temperatures) == len(times)
        and all(_number(value, nonnegative=True) for value in times)
        and all(_number(value, positive=True) for value in temperatures)
    )
    if history_ok:
        time_values = [float(value) for value in times]
        temperature_values = [float(value) for value in temperatures]
        history_ok = time_values[0] == 0.0 and all(left < right for left, right in zip(time_values, time_values[1:])) and all(left <= right for left, right in zip(temperature_values, temperature_values[1:])) and temperature_values[-1] > temperature_values[0]
    energy_ok = all(_number(value, nonnegative=True) for value in (joule, stored, heat_loss)) and _number(input_energy, positive=True)
    return (
        _generation(row, "input_generation", "joule_generation", "thermal_generation", "time_generation", "owner_generation", "result_generation")
        and energy_ok and _close(input_energy, float(joule) + float(stored) + float(heat_loss))
        and history_ok
        and all(row.get("result_" + field) == row.get(field) for field in ("coil_input_energy_j", "joule_heat_energy_j", "stored_thermal_energy_j", "boundary_heat_loss_j", "time_s", "average_temperature_k"))
        and str(row.get("solution_owner") or "").startswith("solution:") and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity(row)
    )


def validate_public_v56_identity(payload: object) -> dict[str, object]:
    """Validate optional v56 acoustic-modal and induction-heating identities."""
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    modal = payload.get(MODAL)
    induction = payload.get(INDUCTION)
    if modal is not None:
        checks["v56_acoustic_modal_participation_response_energy_owner"] = isinstance(modal, Mapping) and _modal_ok(modal)
    if induction is not None:
        checks["v56_induction_input_joule_thermal_time_owner"] = isinstance(induction, Mapping) and _induction_ok(induction)
    if not checks:
        return {}
    return {
        "policy": "multiphysics_identity_v56",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }
