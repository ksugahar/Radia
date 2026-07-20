"""Neutral CST v45 waveguide/EMC identity checks.

The gate is intentionally solver-neutral: it checks lineage, units, owners,
and replayed values without importing CST or exposing private result files.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


_WAVEGUIDE = "waveguide_cutoff_impedance_group_delay_power_orthogonality_mesh_owner_identity"
_EMC = "emc_probe_coordinate_interpolation_window_fft_parseval_monitor_result_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _finite_sequence(value: object) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _waveguide_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("waveguide_generation", "")).strip()
    links = (
        "frequency_generation", "cutoff_generation", "impedance_generation",
        "groupdelay_generation", "power_generation", "orthogonality_generation",
        "mesh_generation", "result_generation",
    )
    arrays = (
        ("frequency_hz", "result_frequency_hz"),
        ("modal_impedance_ohm", "result_modal_impedance_ohm"),
        ("group_delay_s", "result_group_delay_s"),
        ("power_normalization_w", "result_power_normalization_w"),
    )
    matrix = row.get("orthogonality_matrix")
    return (
        bool(generation)
        and all(row.get(key) == generation for key in links)
        and all(_finite_sequence(row.get(left)) and row.get(left) == row.get(right) for left, right in arrays)
        and isinstance(row.get("cutoff_frequency_hz"), (int, float))
        and math.isfinite(float(row["cutoff_frequency_hz"]))
        and float(row["cutoff_frequency_hz"]) > 0.0
        and row.get("result_cutoff_frequency_hz") == row.get("cutoff_frequency_hz")
        and all(float(value) > 0.0 for value in row.get("modal_impedance_ohm", []))
        and all(float(value) >= 0.0 for value in row.get("power_normalization_w", []))
        and isinstance(matrix, list)
        and len(matrix) == 2
        and all(isinstance(line, list) and len(line) == 2 for line in matrix)
        and row.get("result_orthogonality_matrix") == matrix
        and str(row.get("mesh_owner", "")).startswith("mesh:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and row.get("release_id") == row.get("result_release_id") == "cst-v45"
        and _digest(row.get("waveguide_result_sha256"))
        and row.get("accepted_waveguide_result_sha256") == row.get("waveguide_result_sha256")
    )


def _emc_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("emc_probe_generation", "")).strip()
    links = (
        "coordinate_generation", "interpolation_generation", "timewindow_generation",
        "fft_generation", "parseval_generation", "monitor_generation", "result_generation",
    )
    window = row.get("time_window_s")
    return (
        bool(generation)
        and all(row.get(key) == generation for key in links)
        and row.get("coordinate_system") in {"global_cartesian", "global_cylindrical"}
        and row.get("result_coordinate_system") == row.get("coordinate_system")
        and isinstance(row.get("interpolation_order"), int)
        and row.get("interpolation_order", 0) >= 1
        and row.get("result_interpolation_order") == row.get("interpolation_order")
        and isinstance(window, list)
        and len(window) == 2
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in window)
        and window[0] <= window[1]
        and row.get("result_time_window_s") == window
        and row.get("fft_normalization") == row.get("result_fft_normalization") == "parseval_unitary"
        and isinstance(row.get("parseval_time_energy_j"), (int, float))
        and isinstance(row.get("parseval_frequency_energy_j"), (int, float))
        and math.isclose(float(row["parseval_time_energy_j"]), float(row["parseval_frequency_energy_j"]), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and row.get("parseval_time_energy_j") == row.get("result_parseval_time_energy_j")
        and row.get("parseval_frequency_energy_j") == row.get("result_parseval_frequency_energy_j")
        and str(row.get("monitor_owner", "")).startswith("monitor:")
        and row.get("result_monitor_owner") == row.get("monitor_owner")
        and row.get("release_id") == row.get("result_release_id") == "cst-v45"
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
        checks["cst_v45_waveguide_identity"] = len(waveguide) == len(rows) and all(_waveguide_ok(row) for row in waveguide)
    if emc:
        checks["cst_v45_emc_probe_identity"] = len(emc) == len(rows) and all(_emc_ok(row) for row in emc)
    return checks
