"""Waveguide power and SAR-support identity checks for v54."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


CUTOFF = "waveguide_cutoff_mode_normalization_power_impedance_port_owner_identity"
SAR = "sar_average_mass_density_voxel_frequency_field_owner_identity"
_C0 = 299792458.0
_ETA0 = 376.730313668


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _cutoff_ok(row: Mapping[str, object]) -> bool:
    width = row.get("waveguide_width_m")
    cutoff = row.get("cutoff_frequency_hz")
    frequency = row.get("sample_frequency_hz")
    power = row.get("modal_power_w")
    impedance = row.get("mode_impedance_ohm")
    dimensions_ok = _finite(width) and float(width) > 0.0
    cutoff_ok = dimensions_ok and _finite(cutoff) and math.isclose(
        float(cutoff),
        _C0 / (2.0 * float(width)),
        rel_tol=1.0e-12,
        abs_tol=1.0e-6,
    )
    frequency_ok = cutoff_ok and _finite(frequency) and float(frequency) > float(cutoff)
    expected_impedance = (
        _ETA0 / math.sqrt(1.0 - (float(cutoff) / float(frequency)) ** 2)
        if frequency_ok
        else math.nan
    )
    return (
        _generations(
            row,
            "cutoff_generation",
            "mode_generation",
            "normalization_generation",
            "power_generation",
            "impedance_generation",
            "owner_generation",
            "result_generation",
        )
        and row.get("mode_id") == "TE10"
        and row.get("result_mode_id") == row.get("mode_id")
        and dimensions_ok
        and row.get("result_waveguide_width_m") == width
        and cutoff_ok
        and row.get("result_cutoff_frequency_hz") == cutoff
        and frequency_ok
        and row.get("result_sample_frequency_hz") == frequency
        and row.get("field_normalization") == "unit_modal_power_w"
        and row.get("result_field_normalization") == row.get("field_normalization")
        and _finite(power)
        and math.isclose(float(power), 1.0, rel_tol=0.0, abs_tol=1.0e-12)
        and row.get("result_modal_power_w") == power
        and _finite(impedance)
        and math.isclose(float(impedance), expected_impedance, rel_tol=1.0e-12, abs_tol=1.0e-9)
        and row.get("result_mode_impedance_ohm") == impedance
        and str(row.get("port_owner") or "").startswith("port:")
        and row.get("result_port_owner") == row.get("port_owner")
        and _result(row)
    )


def _sar_ok(row: Mapping[str, object]) -> bool:
    mass = row.get("averaging_mass_kg")
    density = row.get("tissue_density_kg_m3")
    voxels = row.get("voxel_support_m3")
    fields = row.get("electric_field_rms_v_m")
    mass_ok = _finite(mass) and any(
        math.isclose(float(mass), expected, rel_tol=0.0, abs_tol=1.0e-12)
        for expected in (0.001, 0.01)
    )
    density_ok = _finite(density) and float(density) > 0.0
    voxels_ok = (
        isinstance(voxels, Sequence)
        and not isinstance(voxels, (str, bytes))
        and bool(voxels)
        and all(_finite(value) and float(value) > 0.0 for value in voxels)
    )
    support_ok = mass_ok and density_ok and voxels_ok and math.isclose(
        sum(float(value) for value in voxels) * float(density),
        float(mass),
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    )
    fields_ok = (
        isinstance(fields, Sequence)
        and not isinstance(fields, (str, bytes))
        and voxels_ok
        and len(fields) == len(voxels)
        and all(_finite(value) and float(value) >= 0.0 for value in fields)
    )
    return (
        _generations(
            row,
            "mass_generation",
            "density_generation",
            "voxel_generation",
            "frequency_generation",
            "field_generation",
            "owner_generation",
            "result_generation",
        )
        and support_ok
        and row.get("result_averaging_mass_kg") == mass
        and row.get("result_tissue_density_kg_m3") == density
        and row.get("result_voxel_support_m3") == voxels
        and fields_ok
        and row.get("result_electric_field_rms_v_m") == fields
        and _finite(row.get("frequency_hz"))
        and float(row["frequency_hz"]) > 0.0
        and row.get("result_frequency_hz") == row.get("frequency_hz")
        and _digest(row.get("field_solution_sha256"))
        and row.get("result_field_solution_sha256") == row.get("field_solution_sha256")
        and str(row.get("monitor_owner") or "").startswith("monitor:")
        and row.get("result_monitor_owner") == row.get("monitor_owner")
        and _result(row)
    )


def validate_public_v54_identity(payload: object) -> dict[str, bool]:
    if not isinstance(payload, Mapping):
        return {}
    rows = [row for row in (payload.get("runs") or []) if isinstance(row, Mapping)]
    checks: dict[str, bool] = {}
    cutoffs = [row[CUTOFF] for row in rows if CUTOFF in row]
    sar_rows = [row[SAR] for row in rows if SAR in row]
    if cutoffs:
        checks["wave_v54_cutoff_mode_power_impedance_port_owner"] = len(cutoffs) == len(rows) and all(
            isinstance(row, Mapping) and _cutoff_ok(row) for row in cutoffs
        )
    if sar_rows:
        checks["wave_v54_sar_mass_density_voxel_frequency_field_owner"] = len(sar_rows) == len(rows) and all(
            isinstance(row, Mapping) and _sar_ok(row) for row in sar_rows
        )
    return checks
