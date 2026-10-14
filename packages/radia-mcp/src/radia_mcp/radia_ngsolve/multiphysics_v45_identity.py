"""Optional identity gates for coupled wave, port, and exported-result cards."""

from __future__ import annotations

import math
from collections.abc import Mapping


def _same(identity: Mapping[str, object], *names: str) -> bool:
    return all(identity.get(f"result_{name}") == identity.get(name) for name in names)


def _sha(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text.lower())


def microwave_v45_ok(summary: Mapping[str, object]) -> bool:
    identity = summary.get("microwave_sparameter_port_reference_plane_deembed_complex_power_mesh_result_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        values = [
            float(identity["frequency_hz"]), float(identity["reference_plane_m"]),
            float(identity["deembed_length_m"]), float(identity["complex_power_w"]),
            float(identity["s11_complex"]["real"]), float(identity["s11_complex"]["imag"]),
            float(identity["s21_complex"]["real"]), float(identity["s21_complex"]["imag"]),
        ]
    except (KeyError, TypeError, ValueError):
        return False
    generations = ("generation", "reference_plane_m_generation", "deembed_length_m_generation", "complex_power_w_generation", "mesh_generation_generation", "port_mode_generation", "s11_complex_generation", "s21_complex_generation")
    return (
        all(bool(str(identity.get(name) or "")) for name in generations)
        and all(math.isfinite(value) for value in values)
        and values[0] > 0.0 and values[1] >= 0.0 and values[2] >= 0.0 and values[3] >= 0.0
        and values[4] ** 2 + values[5] ** 2 <= 1.0
        and values[6] ** 2 + values[7] ** 2 <= 1.0
        and _same(identity, "reference_plane_m", "deembed_length_m", "complex_power_w", "port_mode", "s11_complex", "s21_complex")
        and bool(str(identity.get("owner") or ""))
        and identity.get("accepted_owner") == identity.get("owner")
        and _sha(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def acoustic_v45_ok(summary: Mapping[str, object]) -> bool:
    identity = summary.get("acoustic_impedance_absorption_phase_energy_flux_farfield_window_dataset_result_identity")
    if identity is None:
        return True
    if not isinstance(identity, Mapping):
        return False
    try:
        values = [
            float(identity["frequency_hz"]), float(identity["impedance_magnitude_pa_s_per_m"]),
            float(identity["impedance_phase_deg"]), float(identity["absorption_coefficient"]),
            float(identity["normal_energy_flux_w"]), float(identity["farfield_radius_m"]),
            float(identity["time_window_s"]),
        ]
    except (KeyError, TypeError, ValueError):
        return False
    generations = ("generation", "frequency_hz_generation", "impedance_magnitude_pa_s_per_m_generation", "impedance_phase_deg_generation", "absorption_coefficient_generation", "normal_energy_flux_w_generation", "farfield_radius_m_generation", "time_window_s_generation", "dataset_tag_generation")
    return (
        all(bool(str(identity.get(name) or "")) for name in generations)
        and all(math.isfinite(value) for value in values)
        and values[0] > 0.0 and values[1] > 0.0 and -180.0 <= values[2] <= 180.0
        and 0.0 <= values[3] <= 1.0 and values[4] >= 0.0 and values[5] > 0.0 and values[6] > 0.0
        and _same(identity, "frequency_hz", "impedance_magnitude_pa_s_per_m", "impedance_phase_deg", "absorption_coefficient", "normal_energy_flux_w", "farfield_radius_m", "time_window_s", "dataset_tag")
        and bool(str(identity.get("owner") or ""))
        and identity.get("accepted_owner") == identity.get("owner")
        and _sha(identity.get("result_sha256"))
        and identity.get("accepted_result_sha256") == identity.get("result_sha256")
    )


def multiphysics_v45_ok(summary: Mapping[str, object]) -> bool:
    return microwave_v45_ok(summary) and acoustic_v45_ok(summary)
