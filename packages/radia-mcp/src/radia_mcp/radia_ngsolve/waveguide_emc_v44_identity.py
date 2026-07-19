"""Neutral identity checks for waveguide and EMC result rows."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


_WAVEGUIDE = "waveguide_modal_cutoff_impedance_groupdelay_power_orthogonality_mesh_result_identity"
_EMC = "emc_probe_interpolation_timewindow_fft_parseval_coordinate_monitor_owner_identity"


def _same(row: Mapping[str, object], left: str, right: str) -> bool:
    return row.get(left) == row.get(right)


def _seq(value: object) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and bool(value)


def _finite(value: object) -> bool:
    return _seq(value) and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _waveguide_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("waveguide_generation", "")).strip()
    links = ("frequency_generation", "cutoff_generation", "impedance_generation", "groupdelay_generation", "power_generation", "orthogonality_generation", "mesh_generation", "result_generation")
    arrays = (("frequency_hz", "result_frequency_hz"), ("modal_impedance_ohm", "result_modal_impedance_ohm"), ("group_delay_s", "result_group_delay_s"), ("power_normalization_w", "result_power_normalization_w"))
    return (
        bool(generation)
        and all(row.get(key) == generation for key in links)
        and all(_finite(row.get(left)) and _same(row, left, right) for left, right in arrays)
        and isinstance(row.get("cutoff_frequency_hz"), (int, float))
        and float(row["cutoff_frequency_hz"]) > 0.0
        and row.get("result_cutoff_frequency_hz") == row.get("cutoff_frequency_hz")
        and all(float(value) > 0.0 for value in row.get("modal_impedance_ohm", []))
        and all(float(value) >= 0.0 for value in row.get("power_normalization_w", []))
        and isinstance(row.get("orthogonality_matrix"), list)
        and row.get("result_orthogonality_matrix") == row.get("orthogonality_matrix")
        and str(row.get("mesh_owner", "")).startswith("mesh:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _digest(row.get("waveguide_result_sha256"))
        and row.get("accepted_waveguide_result_sha256") == row.get("waveguide_result_sha256")
    )


def _emc_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("emc_probe_generation", "")).strip()
    links = ("coordinate_generation", "interpolation_generation", "timewindow_generation", "fft_generation", "parseval_generation", "monitor_generation", "result_generation")
    return (
        bool(generation)
        and all(row.get(key) == generation for key in links)
        and row.get("coordinate_system") in {"global_cartesian", "global_cylindrical"}
        and row.get("result_coordinate_system") == row.get("coordinate_system")
        and row.get("interpolation_order") == row.get("result_interpolation_order")
        and isinstance(row.get("interpolation_order"), int)
        and row.get("interpolation_order") >= 1
        and row.get("time_window_s") == row.get("result_time_window_s")
        and isinstance(row.get("time_window_s"), list)
        and len(row["time_window_s"]) == 2
        and row["time_window_s"][0] <= row["time_window_s"][1]
        and row.get("fft_normalization") == row.get("result_fft_normalization") == "parseval_unitary"
        and row.get("parseval_time_energy_j") == row.get("result_parseval_time_energy_j")
        and row.get("parseval_frequency_energy_j") == row.get("result_parseval_frequency_energy_j")
        and math.isclose(float(row["parseval_time_energy_j"]), float(row["parseval_frequency_energy_j"]), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and str(row.get("monitor_owner", "")).startswith("monitor:")
        and row.get("result_monitor_owner") == row.get("monitor_owner")
        and _digest(row.get("emc_probe_result_sha256"))
        and row.get("accepted_emc_probe_result_sha256") == row.get("emc_probe_result_sha256")
    )


def validate_public_identity(payload: object) -> dict[str, bool]:
    if not isinstance(payload, Mapping):
        return {}
    runs = payload.get("runs")
    rows = [row for row in (runs or []) if isinstance(row, Mapping)]
    checks: dict[str, bool] = {}
    waveguide = [row[_WAVEGUIDE] for row in rows if _WAVEGUIDE in row]
    emc = [row[_EMC] for row in rows if _EMC in row]
    if waveguide:
        checks["waveguide_v44_modal_identity"] = len(waveguide) == len(rows) and all(isinstance(row, Mapping) and _waveguide_ok(row) for row in waveguide)
    if emc:
        checks["waveguide_v44_emc_probe_identity"] = len(emc) == len(rows) and all(isinstance(row, Mapping) and _emc_ok(row) for row in emc)
    return checks
