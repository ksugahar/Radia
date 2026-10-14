from __future__ import annotations

import json
import math
from copy import deepcopy

from radia_mcp.radia_ngsolve.server import harmonic_current_port_power_energy_identity_gate


def _summary() -> dict:
    frequency = 1000.0
    omega = 2.0 * math.pi * frequency
    current = {"real": 1.0, "imag": 0.0}
    impedance = {"real": 3.2e-7, "imag": 4.8e-7}
    voltage = dict(impedance)
    power = {"real": 1.6e-7, "imag": 2.4e-7}
    flux = {"real": 4.8e-7 / omega, "imag": -1.0e-11}
    record = {
        "current_a": current,
        "voltage_v": voltage,
        "impedance_ohm": impedance,
        "complex_power_w": power,
        "flux_linkage_wb": flux,
        "magnetic_energy_j": 2.4e-7 / (2.0 * omega),
        "loss_w": 1.6e-7,
        "conductor_loss_w": 1.6e-7,
        "current_density_profile": [[0.0, 2.0], [1.0, 3.0]],
    }
    return {
        "frequency_hz": frequency,
        "amplitude_convention": "peak_phasor",
        "runs": [
            dict(deepcopy(record), label="single_point"),
            dict(deepcopy(record), label="sweep_point"),
        ],
    }


def _call(summary: dict) -> dict:
    return json.loads(harmonic_current_port_power_energy_identity_gate(json.dumps(summary)))


def test_harmonic_port_identity_gate_accepts_closed_peak_phasor_pair():
    result = _call(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["profile_sample_count"] == 2


def test_harmonic_port_identity_gate_rejects_rms_label_for_peak_values():
    summary = _summary()
    summary["amplitude_convention"] = "rms_phasor"
    result = _call(summary)
    assert result["status"] == "needs_attention"
    assert "peak_phasor_convention_explicit" in result["issues"]


def test_harmonic_port_identity_gate_rejects_unclosed_loss():
    summary = _summary()
    summary["runs"][1]["loss_w"] *= 1.1
    result = _call(summary)
    assert result["status"] == "needs_attention"
    assert "all_port_identities_close" in result["issues"]


def test_harmonic_port_identity_gate_rejects_changed_profile():
    summary = _summary()
    summary["runs"][1]["current_density_profile"][1][1] *= 1.2
    result = _call(summary)
    assert result["status"] == "needs_attention"
    assert "normalized_cross_run_observables_close" in result["issues"]
