"""Waveguide-port and far-field identity checks for v53."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .wave_sar_identity_v54 import validate_public_v54_identity


WAVEGUIDE = "waveguide_mode_cutoff_normalization_referenceplane_port_owner_identity"
FARFIELD = "farfield_realizedgain_polarization_basis_angulargrid_monitor_owner_identity"
_C0 = 299792458.0


def _digest(value: object) -> bool:
    text = str(value or "").lower(); return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    generation = str(row.get("generation") or ""); return bool(generation) and all(row.get(name) == generation for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _waveguide_ok(row: Mapping[str, object]) -> bool:
    width = row.get("waveguide_width_m"); cutoff = row.get("cutoff_frequency_hz"); sample = row.get("sample_frequency_hz")
    expected = _C0 / (2.0 * float(width)) if _finite(width) and float(width) > 0.0 else math.nan
    return (
        _generation(row, ("mode_generation", "cutoff_generation", "normalization_generation", "plane_generation", "owner_generation", "result_generation"))
        and row.get("mode_id") == "TE10" and row.get("result_mode_id") == row.get("mode_id")
        and _finite(width) and float(width) > 0.0 and row.get("result_waveguide_width_m") == width
        and _finite(cutoff) and math.isclose(float(cutoff), expected, rel_tol=1.0e-12, abs_tol=1.0e-6) and row.get("result_cutoff_frequency_hz") == cutoff
        and _finite(sample) and float(sample) > float(cutoff) and row.get("result_sample_frequency_hz") == sample
        and row.get("normalization") == "unit_incident_power_w" and row.get("result_normalization") == row.get("normalization")
        and _finite(row.get("reference_plane_offset_m")) and row.get("result_reference_plane_offset_m") == row.get("reference_plane_offset_m")
        and str(row.get("port_owner") or "").startswith("port:") and row.get("result_port_owner") == row.get("port_owner") and _result(row)
    )


def _grid(value: object, lower: float, upper: float, include_upper: bool) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) >= 2 and all(_finite(item) and lower <= float(item) <= upper if include_upper else _finite(item) and lower <= float(item) < upper for item in value) and all(float(left) < float(right) for left, right in zip(value, value[1:]))


def _farfield_ok(row: Mapping[str, object]) -> bool:
    theta = row.get("theta_deg"); phi = row.get("phi_deg"); gain = row.get("realized_gain_dbi")
    theta_ok = _grid(theta, 0.0, 180.0, True) and math.isclose(float(theta[0]), 0.0, abs_tol=1.0e-12) and math.isclose(float(theta[-1]), 180.0, abs_tol=1.0e-12)
    phi_ok = _grid(phi, 0.0, 360.0, False)
    gain_ok = isinstance(gain, Sequence) and not isinstance(gain, (str, bytes)) and theta_ok and phi_ok and len(gain) == len(theta) and all(isinstance(values, Sequence) and not isinstance(values, (str, bytes)) and len(values) == len(phi) and all(_finite(value) for value in values) for values in gain)
    return (
        _generation(row, ("gain_generation", "polarization_generation", "grid_generation", "owner_generation", "result_generation"))
        and theta_ok and row.get("result_theta_deg") == theta and phi_ok and row.get("result_phi_deg") == phi
        and gain_ok and row.get("result_realized_gain_dbi") == gain
        and row.get("polarization_basis") == "ludwig3" and row.get("result_polarization_basis") == row.get("polarization_basis")
        and str(row.get("monitor_owner") or "").startswith("monitor:") and row.get("result_monitor_owner") == row.get("monitor_owner") and _result(row)
    )


def validate_public_v53_identity(payload: object) -> dict[str, bool]:
    if not isinstance(payload, Mapping): return {}
    rows = [row for row in (payload.get("runs") or []) if isinstance(row, Mapping)]; checks = validate_public_v54_identity(payload)
    waveguides = [row[WAVEGUIDE] for row in rows if WAVEGUIDE in row]; farfields = [row[FARFIELD] for row in rows if FARFIELD in row]
    if waveguides: checks["wave_v53_waveguide_cutoff_normalization_plane_port_owner"] = len(waveguides) == len(rows) and all(isinstance(row, Mapping) and _waveguide_ok(row) for row in waveguides)
    if farfields: checks["wave_v53_farfield_gain_polarization_grid_monitor_owner"] = len(farfields) == len(rows) and all(isinstance(row, Mapping) and _farfield_ok(row) for row in farfields)
    return checks
