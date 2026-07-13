from __future__ import annotations

import copy
import json
import math

from radia_mcp.radia_ngsolve.cylindrical_conductor_skin_gate import (
    MU0,
    _bessel_j0,
    _bessel_j1,
    cylindrical_conductor_skin_bessel_gate,
)
from radia_mcp.radia_ngsolve.server import cylindrical_conductor_skin_bessel_gate as mcp_gate


def _summary() -> dict:
    frequency, radius, length, sigma = 10000.0, 5.0e-3, 1.0e-3, 5.8e7
    omega = 2.0 * math.pi * frequency
    delta = math.sqrt(2.0 / (omega * MU0 * sigma))
    kappa = (1.0 + 1.0j) / delta
    z_internal = kappa * _bessel_j0(kappa * radius) / (2.0 * math.pi * radius * sigma * _bessel_j1(kappa * radius)) * length
    rdc = length / (sigma * math.pi * radius * radius)
    impedance = complex(z_internal.real, abs(z_internal.imag) + 1.0e-4)
    induced = impedance - rdc
    flux = induced / (1.0j * omega)
    power = 0.5 * impedance
    radii = [radius * index / 20.0 for index in range(20)]
    density = [10.0 * _bessel_j0(kappa * value) for value in radii]
    replay = {
        "current_a": [1.0, 0.0], "flux_linkage_wb": [flux.real, flux.imag],
        "impedance_ohm": [impedance.real, impedance.imag], "power_va": [power.real, power.imag],
        "voltage_v": [impedance.real, impedance.imag], "induced_voltage_v": [induced.real, induced.imag],
        "energy_j": power.imag / (2.0 * omega), "loss_w": power.real,
        "final_log10_relative_residual": -12.0, "profile_radii_m": radii,
        "current_density_a_per_m2": [[value.real, value.imag] for value in density],
    }
    return {
        "model_contract": {"frequency_hz": frequency, "radius_m": radius, "length_m": length, "conductivity_s_per_m": sigma},
        "replays": [copy.deepcopy(replay), copy.deepcopy(replay)],
        "timing_breakdown_s": {"stage": 1.0, "solve": 2.0, "read": 1.0, "verify": 1.0},
    }


def test_accepts_bessel_skin_solution_and_dispatches():
    result = cylindrical_conductor_skin_bessel_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["bessel_profile_relative_l2_error"] < 1.0e-12
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_resistance_profile_and_replay_drift():
    row = copy.deepcopy(_summary())
    for key in ("impedance_ohm", "voltage_v", "power_va"):
        row["replays"][0][key][0] *= 1.2
    row["replays"][0]["loss_w"] *= 1.2
    row["replays"][1]["current_density_a_per_m2"][-1][0] *= 0.1
    row["replays"][1]["current_density_a_per_m2"][-1][1] *= 0.1
    result = cylindrical_conductor_skin_bessel_gate(row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["ac_resistance_matches_bessel_internal_impedance"] is False
    assert result["checks"]["current_density_matches_bessel_profile"] is False
    assert result["checks"]["independent_replays_are_deterministic"] is False


def test_rejects_residual_and_timing_failures():
    row = _summary()
    row["replays"][1]["final_log10_relative_residual"] = -5.0
    row["timing_breakdown_s"].pop("verify")
    result = cylindrical_conductor_skin_bessel_gate(row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["both_linear_residuals_are_below_1e_10"] is False
    assert result["checks"]["exactly_four_timing_stages"] is False
