"""Resonator-Q and antenna-power artifact checks for v55."""

from __future__ import annotations

import math
from collections.abc import Mapping


RESONATOR = "resonator_loaded_unloaded_q_coupling_linewidth_energy_owner_identity"
ANTENNA = "antenna_efficiency_accepted_radiated_loss_gain_directivity_owner_identity"


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _close(value: object, expected: float) -> bool:
    return _finite(value) and math.isclose(float(value), expected, rel_tol=1.0e-10, abs_tol=1.0e-12)


def _resonator_ok(row: Mapping[str, object]) -> bool:
    frequency = row.get("resonance_frequency_hz"); loaded = row.get("loaded_q"); unloaded = row.get("unloaded_q"); external = row.get("external_q"); power = row.get("dissipated_power_w")
    primary_ok = all(_finite(value) and float(value) > 0.0 for value in (frequency, loaded, unloaded, external, power))
    if not primary_ok:
        return False
    loaded_expected = 1.0 / (1.0 / float(unloaded) + 1.0 / float(external))
    beta = float(unloaded) / float(external); linewidth = float(frequency) / float(loaded); energy = float(loaded) * float(power) / (2.0 * math.pi * float(frequency))
    return (_generations(row, "q_generation", "coupling_generation", "linewidth_generation", "energy_generation", "owner_generation", "result_generation") and _close(loaded, loaded_expected) and _close(row.get("coupling_beta"), beta) and _close(row.get("linewidth_hz"), linewidth) and _close(row.get("stored_energy_j"), energy) and row.get("result_resonance_frequency_hz") == frequency and row.get("result_loaded_q") == loaded and row.get("result_unloaded_q") == unloaded and row.get("result_external_q") == external and row.get("result_coupling_beta") == row.get("coupling_beta") and row.get("result_linewidth_hz") == row.get("linewidth_hz") and row.get("result_dissipated_power_w") == power and row.get("result_stored_energy_j") == row.get("stored_energy_j") and str(row.get("monitor_owner") or "").startswith("monitor:") and row.get("result_monitor_owner") == row.get("monitor_owner") and _result(row))


def _antenna_ok(row: Mapping[str, object]) -> bool:
    accepted = row.get("accepted_power_w"); radiated = row.get("radiated_power_w"); loss = row.get("loss_power_w"); efficiency = row.get("radiation_efficiency"); directivity = row.get("directivity_dbi"); gain = row.get("gain_dbi")
    powers_ok = all(_finite(value) for value in (accepted, radiated, loss)) and float(accepted) > 0.0 and float(radiated) >= 0.0 and float(loss) >= 0.0 and math.isclose(float(accepted), float(radiated) + float(loss), rel_tol=1.0e-10, abs_tol=1.0e-12)
    if not powers_ok:
        return False
    eta = float(radiated) / float(accepted)
    return (_generations(row, "power_generation", "efficiency_generation", "gain_generation", "directivity_generation", "owner_generation", "result_generation") and _close(efficiency, eta) and 0.0 < eta <= 1.0 and _finite(directivity) and _close(gain, float(directivity) + 10.0 * math.log10(eta)) and row.get("result_accepted_power_w") == accepted and row.get("result_radiated_power_w") == radiated and row.get("result_loss_power_w") == loss and row.get("result_radiation_efficiency") == efficiency and row.get("result_directivity_dbi") == directivity and row.get("result_gain_dbi") == gain and str(row.get("farfield_owner") or "").startswith("farfield:") and row.get("result_farfield_owner") == row.get("farfield_owner") and _result(row))


def validate_public_v55_identity(payload: object) -> dict[str, bool]:
    if not isinstance(payload, Mapping):
        return {}
    rows = [row for row in (payload.get("runs") or []) if isinstance(row, Mapping)]; checks: dict[str, bool] = {}
    resonators = [row[RESONATOR] for row in rows if RESONATOR in row]; antennas = [row[ANTENNA] for row in rows if ANTENNA in row]
    if resonators:
        checks["wave_v55_resonator_q_coupling_linewidth_energy_owner"] = len(resonators) == len(rows) and all(isinstance(item, Mapping) and _resonator_ok(item) for item in resonators)
    if antennas:
        checks["wave_v55_antenna_power_efficiency_gain_directivity_owner"] = len(antennas) == len(rows) and all(isinstance(item, Mapping) and _antenna_ok(item) for item in antennas)
    return checks
