"""Motor dq and electrothermal artifact identity checks for v50."""

from __future__ import annotations

import math
from collections.abc import Mapping


DQ = "dq_park_angle_current_phase_saliency_operating_point_owner_identity"
THERMAL = "thermal_loss_map_convection_temperature_material_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _numeric_map(value: object, keys: set[str], *, positive: bool = False) -> bool:
    return (
        isinstance(value, Mapping)
        and set(value) == keys
        and all(_finite(item) and (not positive or float(item) > 0.0) for item in value.values())
    )


def _dq_ok(row: Mapping[str, object]) -> bool:
    angle = row.get("park_angle_electrical_deg")
    currents = row.get("dq_currents")
    saliency = row.get("saliency_parameters")
    operating_point = str(row.get("operating_point_id") or "")
    owner = str(row.get("result_owner") or "")
    currents_ok = (
        isinstance(currents, Mapping)
        and set(currents) == {"id_a", "iq_a", "phase_order"}
        and _finite(currents.get("id_a"))
        and _finite(currents.get("iq_a"))
        and currents.get("phase_order") in {"uvw", "vwu", "wuv"}
    )
    return (
        _generations(row, "angle_generation", "current_generation", "saliency_generation", "operating_point_generation", "result_generation")
        and _finite(angle)
        and -360.0 <= float(angle) <= 360.0
        and row.get("result_park_angle_electrical_deg") == angle
        and currents_ok
        and row.get("result_dq_currents") == currents
        and _numeric_map(saliency, {"ld_h", "lq_h", "psi_pm_wb"}, positive=True)
        and row.get("result_saliency_parameters") == saliency
        and operating_point.startswith("operating-point:")
        and row.get("result_operating_point_id") == operating_point
        and owner.startswith("dq-result:")
        and row.get("accepted_result_owner") == owner
        and _result(row)
    )


def _loss_map(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(str(name).endswith("_w") for name in value)
        and all(_finite(loss) and float(loss) >= 0.0 for loss in value.values())
    )


def _convection(value: object) -> bool:
    if not isinstance(value, list) or not value or not all(isinstance(row, Mapping) for row in value):
        return False
    boundaries = [str(row.get("boundary") or "") for row in value]
    return (
        all(boundaries)
        and len(set(boundaries)) == len(boundaries)
        and all(_finite(row.get("h_w_m2k")) and float(row["h_w_m2k"]) >= 0.0 for row in value)
        and all(_finite(row.get("ambient_c")) and float(row["ambient_c"]) > -273.15 for row in value)
    )


def _temperatures(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(_finite(item) and float(item) > -273.15 for item in value.values())
    )


def _materials(value: object) -> bool:
    return isinstance(value, Mapping) and bool(value) and all(bool(str(name)) and ":" in str(revision) for name, revision in value.items())


def _thermal_ok(row: Mapping[str, object]) -> bool:
    loss = row.get("loss_map")
    convection = row.get("convection_boundaries")
    temperatures = row.get("temperature_map")
    materials = row.get("material_revisions")
    owner = str(row.get("thermal_owner") or "")
    return (
        _generations(row, "loss_generation", "boundary_generation", "temperature_generation", "material_generation", "result_generation")
        and _loss_map(loss)
        and row.get("replayed_loss_map") == loss
        and _convection(convection)
        and row.get("replayed_convection_boundaries") == convection
        and _temperatures(temperatures)
        and row.get("replayed_temperature_map") == temperatures
        and _materials(materials)
        and row.get("replayed_material_revisions") == materials
        and owner.startswith("thermal-result:")
        and row.get("replayed_thermal_owner") == owner
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    """Validate optional v50 dq and electrothermal identity records."""
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    dq = identity.get(DQ)
    thermal = identity.get(THERMAL)
    if dq is not None:
        checks["motor_v50_dq_park_current_saliency_operating_point_owner"] = isinstance(dq, Mapping) and _dq_ok(dq)
    if thermal is not None:
        checks["motor_v50_thermal_loss_convection_temperature_material_owner"] = isinstance(thermal, Mapping) and _thermal_ok(thermal)
    return checks
