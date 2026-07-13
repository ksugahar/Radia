from __future__ import annotations

import copy
import json
import math

from radia_mcp.radia_ngsolve.motion_coupled_levitation_gate import (
    motion_coupled_eddy_levitation_transient_gate as gate,
)
from radia_mcp.radia_ngsolve.server import (
    motion_coupled_eddy_levitation_transient_gate as mcp_gate,
)


def _summary() -> dict:
    output_time = [0.01 * index for index in range(101)]
    probe_time = [0.001 * index for index in range(1001)]
    mass = 0.1
    gravity = 9.8

    def replay(offset: float) -> dict:
        probe_force = [
            0.4 + 0.3 * math.cos(2.0 * math.pi * 100.0 * time_s) + offset
            for time_s in probe_time
        ]
        probe_displacement = [
            0.01 + 0.002 * math.sin(2.0 * math.pi * time_s) + offset * 1.0e-4
            for time_s in probe_time
        ]
        output_force = [0.7 + offset for _ in output_time]
        acceleration = [(force - mass * gravity) / mass for force in output_force]
        return {
            "label": "replay",
            "solve_seconds": 2.0,
            "output_time_s": output_time,
            "output_displacement_m": [
                0.01 + 0.002 * math.sin(2.0 * math.pi * time_s)
                for time_s in output_time
            ],
            "output_velocity_m_s": [
                0.004 * math.pi * math.cos(2.0 * math.pi * time_s)
                for time_s in output_time
            ],
            "output_acceleration_m_s2": acceleration,
            "output_lift_force_n": output_force,
            "output_gravity_force_n": [mass * gravity for _ in output_time],
            "probe_time_s": probe_time,
            "probe_displacement_m": probe_displacement,
            "probe_lift_force_n": probe_force,
            "adaptive_probe_row_count": 2201,
            "adaptive_probe_median_samples_per_force_period": 22.0,
        }

    first = replay(0.0)
    second = replay(1.0e-4)
    experiment_time = [0.01 * index for index in range(101)]
    return {
        "contract": {
            "drive_frequency_hz": 50.0,
            "expected_force_harmonic_hz": 100.0,
            "force_observation": "adaptive_internal_steps_interpolated_to_fixed_grid",
            "motion_equation": "mass_times_acceleration_equals_lift_minus_gravity",
            "experimental_comparison": "displacement_time_history",
        },
        "units": {
            "time": "s",
            "displacement": "m",
            "velocity": "m/s",
            "acceleration": "m/s^2",
            "force": "N",
            "mass": "kg",
        },
        "mass_kg": mass,
        "gravity_m_s2": gravity,
        "replays": [first, second],
        "experiment": {
            "time_s": experiment_time,
            "displacement_m": [
                0.01 + 0.002 * math.sin(2.0 * math.pi * time_s)
                for time_s in experiment_time
            ],
        },
        "timing_breakdown_s": {
            "stage": 0.1,
            "solve": 4.0,
            "extract": 0.2,
            "verify": 0.1,
        },
    }


def test_accepts_alias_aware_motion_coupled_replay_and_mcp_dispatch() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["requested_output_alias_risk"] is True
    assert result["metrics"]["maximum_requested_output_to_probe_force_span_ratio"] < 0.01
    assert json.loads(mcp_gate(json.dumps(_summary())))["status"] == "ok"


def test_rejects_requested_output_substituted_for_adaptive_probe() -> None:
    summary = _summary()
    summary["contract"]["force_observation"] = "requested_output_only"
    for replay in summary["replays"]:
        replay["probe_time_s"] = replay["output_time_s"]
        replay["probe_displacement_m"] = replay["output_displacement_m"]
        replay["probe_lift_force_n"] = replay["output_lift_force_n"]
        replay["adaptive_probe_row_count"] = len(replay["output_time_s"])
        replay["adaptive_probe_median_samples_per_force_period"] = 1.0
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["force_sampling_strategy_resolves_alias"] is False
    assert result["checks"]["aliased_output_force_not_substituted_for_probe"] is False


def test_rejects_motion_equation_replay_and_experiment_failures() -> None:
    summary = copy.deepcopy(_summary())
    summary["replays"][1]["output_acceleration_m_s2"][20] += 2.0
    summary["replays"][1]["probe_lift_force_n"][30] += 0.2
    summary["experiment"]["displacement_m"] = [
        0.02 - value for value in summary["experiment"]["displacement_m"]
    ]
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert result["checks"]["motion_equation_closes"] is False
    assert result["checks"]["force_history_replays"] is False
    assert result["checks"]["experimental_shape_is_correlated"] is False
