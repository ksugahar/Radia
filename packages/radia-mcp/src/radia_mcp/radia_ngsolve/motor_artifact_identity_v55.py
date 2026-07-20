"""Dq operating-point and iron-loss artifact checks for v55."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


DQ = "dq_fluxlinkage_current_angle_saliency_torque_owner_identity"
IRON = "ironloss_hysteresis_eddy_excess_frequency_fluxdensity_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _dq_pair(value: object, *, positive: bool = False) -> bool:
    return (
        isinstance(value, Mapping)
        and set(value) == {"d", "q"}
        and all(_finite(component) for component in value.values())
        and (not positive or all(float(component) > 0.0 for component in value.values()))
    )


def _close(left: object, right: float) -> bool:
    return _finite(left) and math.isclose(float(left), right, rel_tol=1.0e-10, abs_tol=1.0e-12)


def _dq_ok(row: Mapping[str, object]) -> bool:
    pole_pairs = row.get("pole_pairs")
    current = row.get("current_dq_a")
    inductance = row.get("inductance_dq_h")
    flux = row.get("flux_linkage_dq_wb")
    pm_flux = row.get("pm_flux_linkage_wb")
    values_ok = (
        isinstance(pole_pairs, int)
        and not isinstance(pole_pairs, bool)
        and pole_pairs > 0
        and _dq_pair(current)
        and _dq_pair(inductance, positive=True)
        and _dq_pair(flux)
        and _finite(pm_flux)
        and float(pm_flux) >= 0.0
    )
    if not values_ok:
        return False
    expected_flux_d = float(pm_flux) + float(inductance["d"]) * float(current["d"])
    expected_flux_q = float(inductance["q"]) * float(current["q"])
    expected_angle = math.degrees(math.atan2(float(current["q"]), float(current["d"]))) % 360.0
    expected_saliency = 1.5 * pole_pairs * (float(inductance["d"]) - float(inductance["q"])) * float(current["d"]) * float(current["q"])
    expected_torque = 1.5 * pole_pairs * (float(flux["d"]) * float(current["q"]) - float(flux["q"]) * float(current["d"]))
    return (
        _generations(row, "flux_generation", "current_generation", "angle_generation", "saliency_generation", "torque_generation", "owner_generation", "result_generation")
        and _close(flux["d"], expected_flux_d)
        and _close(flux["q"], expected_flux_q)
        and _close(row.get("current_electrical_angle_deg"), expected_angle)
        and _close(row.get("saliency_torque_nm"), expected_saliency)
        and _close(row.get("torque_nm"), expected_torque)
        and row.get("result_pole_pairs") == pole_pairs
        and row.get("result_current_dq_a") == current
        and row.get("result_inductance_dq_h") == inductance
        and row.get("result_pm_flux_linkage_wb") == pm_flux
        and row.get("result_flux_linkage_dq_wb") == flux
        and row.get("result_current_electrical_angle_deg") == row.get("current_electrical_angle_deg")
        and row.get("result_saliency_torque_nm") == row.get("saliency_torque_nm")
        and row.get("result_torque_nm") == row.get("torque_nm")
        and str(row.get("result_owner") or "").startswith("result:")
        and row.get("accepted_result_owner") == row.get("result_owner")
        and _result(row)
    )


def _iron_ok(row: Mapping[str, object]) -> bool:
    components = row.get("loss_components")
    waveform = row.get("flux_density_waveform_t")
    components_ok = (
        isinstance(components, Mapping)
        and set(components) == {"hysteresis_w", "eddy_w", "excess_w"}
        and all(_finite(value) and float(value) >= 0.0 for value in components.values())
    )
    waveform_ok = (
        isinstance(waveform, Sequence)
        and not isinstance(waveform, (str, bytes))
        and len(waveform) >= 3
        and all(_finite(value) for value in waveform)
        and math.isclose(float(waveform[0]), float(waveform[-1]), rel_tol=0.0, abs_tol=1.0e-12)
        and min(float(value) for value in waveform) < 0.0 < max(float(value) for value in waveform)
    )
    if not components_ok or not waveform_ok:
        return False
    total = sum(float(value) for value in components.values())
    peak = max(abs(float(value)) for value in waveform)
    return (
        _generations(row, "component_generation", "frequency_generation", "waveform_generation", "material_generation", "owner_generation", "result_generation")
        and _close(row.get("total_iron_loss_w"), total)
        and _finite(row.get("frequency_hz"))
        and float(row["frequency_hz"]) > 0.0
        and _close(row.get("peak_flux_density_t"), peak)
        and peak > 0.0
        and row.get("result_loss_components") == components
        and row.get("result_total_iron_loss_w") == row.get("total_iron_loss_w")
        and row.get("result_frequency_hz") == row.get("frequency_hz")
        and row.get("result_flux_density_waveform_t") == waveform
        and row.get("result_peak_flux_density_t") == row.get("peak_flux_density_t")
        and bool(str(row.get("material_revision") or ""))
        and row.get("result_material_revision") == row.get("material_revision")
        and str(row.get("material_owner") or "").startswith("material:")
        and row.get("result_material_owner") == row.get("material_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    dq = identity.get(DQ)
    iron = identity.get(IRON)
    if dq is not None:
        checks["motor_v55_dq_flux_current_angle_saliency_torque_owner"] = isinstance(dq, Mapping) and _dq_ok(dq)
    if iron is not None:
        checks["motor_v55_ironloss_components_frequency_waveform_material_owner"] = isinstance(iron, Mapping) and _iron_ok(iron)
    return checks
