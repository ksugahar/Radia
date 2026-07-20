"""MOS operating-point and transmission-line replay identity gates."""

from __future__ import annotations

import math
from collections.abc import Mapping


MOS = "mosfet_operatingpoint_region_gm_gds_capacitance_temperature_owner_identity"
TLINE = "transmissionline_delay_impedance_termination_reflection_event_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(contract: Mapping[str, object], *names: str) -> bool:
    generation = str(contract.get("generation_id") or "")
    return bool(generation) and all(contract.get(name) == generation for name in names)


def _number(value: object, *, positive: bool = False, nonnegative: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(float(value)):
        return False
    return (not positive or float(value) > 0.0) and (not nonnegative or float(value) >= 0.0)


def _same(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-15)


def _mos_ok(contract: Mapping[str, object]) -> bool:
    vgs, vds, vth, drain_current = (contract.get(name) for name in ("vgs_v", "vds_v", "vth_v", "id_a"))
    gm, gds, temperature = (contract.get(name) for name in ("gm_s", "gds_s", "model_temperature_c"))
    capacitances = contract.get("capacitances")
    if not all(_number(value) for value in (vgs, vds, vth, drain_current, gm, gds, temperature)):
        return False
    expected_region = "cutoff" if float(vgs) <= float(vth) else ("linear" if float(vds) < float(vgs) - float(vth) else "saturation")
    capacitance_ok = isinstance(capacitances, Mapping) and set(capacitances) == {"cgs_f", "cgd_f", "cds_f"} and all(_number(value, nonnegative=True) for value in capacitances.values())
    operating_point_ok = (expected_region == "cutoff" and math.isclose(float(drain_current), 0.0, abs_tol=1.0e-15)) or (expected_region != "cutoff" and float(drain_current) > 0.0 and float(gm) > 0.0 and float(gds) >= 0.0)
    return (
        _generation(contract, "bias_generation_id", "region_generation_id", "small_signal_generation_id", "capacitance_generation_id", "temperature_generation_id", "owner_generation_id", "result_generation_id")
        and contract.get("region") == expected_region
        and operating_point_ok
        and (expected_region != "saturation" or float(gm) > float(gds))
        and capacitance_ok
        and float(temperature) > -273.15
        and all(contract.get("result_" + name) == contract.get(name) for name in ("vgs_v", "vds_v", "vth_v", "id_a", "region", "gm_s", "gds_s", "capacitances", "model_temperature_c"))
        and str(contract.get("device_owner") or "").startswith("device:")
        and contract.get("result_device_owner") == contract.get("device_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def _tline_ok(contract: Mapping[str, object]) -> bool:
    delay = contract.get("propagation_delay_s")
    impedance = contract.get("characteristic_impedance_ohm")
    termination = contract.get("termination_ohm")
    coefficient = contract.get("reflection_coefficient")
    event = contract.get("reflection_event")
    if not (_number(delay, positive=True) and _number(impedance, positive=True) and _number(termination, nonnegative=True) and _number(coefficient) and isinstance(event, Mapping)):
        return False
    expected = (float(termination) - float(impedance)) / (float(termination) + float(impedance))
    launch = event.get("launch_time_s")
    reflection = event.get("reflection_time_s")
    incident = event.get("incident_v")
    reflected = event.get("reflected_v")
    event_ok = all(_number(value) for value in (launch, reflection, incident, reflected)) and _same(reflection, float(launch) + 2.0 * float(delay)) and _same(reflected, float(incident) * expected)
    return (
        _generation(contract, "delay_generation_id", "impedance_generation_id", "termination_generation_id", "reflection_generation_id", "event_generation_id", "owner_generation_id", "result_generation_id")
        and _same(coefficient, expected)
        and event_ok
        and all(contract.get("result_" + name) == contract.get(name) for name in ("propagation_delay_s", "characteristic_impedance_ohm", "termination_ohm", "reflection_coefficient", "reflection_event"))
        and str(contract.get("waveform_owner") or "").startswith("waveform:")
        and contract.get("result_waveform_owner") == contract.get("waveform_owner")
        and _digest(contract.get("result_sha256"))
        and contract.get("accepted_result_sha256") == contract.get("result_sha256")
    )


def validate_ltspice_v52_identity(positive: Mapping[str, object]) -> bool:
    if not isinstance(positive, Mapping):
        return False
    mos = positive.get(MOS)
    tline = positive.get(TLINE)
    if mos is None and tline is None:
        return True
    return isinstance(mos, Mapping) and isinstance(tline, Mapping) and _mos_ok(mos) and _tline_ok(tline)
