from __future__ import annotations

import copy
import json
import math

from radia_mcp.radia_ngsolve.rotational_eddy_brake_energy_gate import (
    rotational_eddy_brake_energy_gate as gate,
)
from radia_mcp.radia_ngsolve.server import rotational_eddy_brake_energy_gate as mcp_gate


def _summary() -> dict:
    density, radius, thickness = 7800.0, 0.1, 0.01
    inertia = 0.5 * density * math.pi * radius**4 * thickness
    times = [0.05 * index for index in range(201)]
    rate = 0.2
    omega = [100.0 * math.exp(-rate * value) for value in times]
    torque = [inertia * rate * value for value in omega]
    joule = [force * speed for force, speed in zip(torque, omega, strict=True)]

    def replay(label: str) -> dict:
        return {
            "label": label,
            "solve_seconds": 2.0,
            "time_s": times,
            "angular_velocity_rad_s": omega,
            "braking_torque_nm": torque,
            "joule_loss_w": joule,
        }

    return {
        "contract": {
            "body": "uniform_solid_conducting_disc",
            "inertia_reference": "analytic_uniform_solid_disc",
            "angular_momentum_balance": "inertia_delta_angular_velocity_plus_integrated_braking_torque_equals_zero",
            "instantaneous_power_comparison": "diagnostic_only_when_field_energy_rate_is_not_sampled_on_the_probe_grid",
            "energy_balance": "initial_kinetic_plus_magnetic_equals_final_kinetic_plus_magnetic_plus_joule",
        },
        "units": {
            "time": "s",
            "angular_velocity": "rad/s",
            "torque": "N*m",
            "power": "W",
            "inertia": "kg*m^2",
            "energy": "J",
            "density": "kg/m^3",
            "length": "m",
        },
        "disc": {
            "density_kg_m3": density,
            "radius_m": radius,
            "thickness_m": thickness,
        },
        "reported_inertia_kg_m2": inertia,
        "replays": [replay("one"), replay("two")],
        "energy_replay": {
            **replay("energy"),
            "field_energy_time_s": times,
            "magnetic_energy_j": [0.3 for _ in times],
        },
        "timing_breakdown_s": {
            "attach": 0.1,
            "solve": 6.0,
            "extract": 0.1,
            "verify": 0.1,
        },
    }


def test_accepts_free_brake_with_field_energy_and_mcp_dispatch() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["kinetic_magnetic_joule_energy_closes"] is True
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_wrong_inertia_and_angular_impulse() -> None:
    summary = copy.deepcopy(_summary())
    summary["reported_inertia_kg_m2"] *= 1.25
    summary["contract"]["inertia_reference"] = "unchecked_reported_value"
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["analytic_disc_inertia_matches_reported"] is False
    assert result["checks"]["angular_impulse_balance_closes"] is False


def test_rejects_missing_field_storage_and_false_power_claim() -> None:
    summary = copy.deepcopy(_summary())
    summary["contract"]["instantaneous_power_comparison"] = "always_equal_to_joule"
    summary["energy_replay"]["magnetic_energy_j"] = []
    try:
        result = gate(summary)
    except ValueError as exc:
        assert "magnetic_energy_j" in str(exc)
    else:
        assert result["status"] == "needs_attention"


def test_rejects_nonreplaying_torque_history() -> None:
    summary = copy.deepcopy(_summary())
    summary["replays"][1]["braking_torque_nm"] = list(
        summary["replays"][1]["braking_torque_nm"]
    )
    summary["replays"][1]["braking_torque_nm"][50] *= 1.2
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["fresh_replay_fields_match"] is False
