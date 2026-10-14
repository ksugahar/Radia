import copy
import json

from radia_mcp.motor.server import motor_transient_no_load_load_cycle_gate
from radia_mcp.motor.transient_no_load_load_gate import (
    motor_transient_no_load_load_cycle_gate as build_gate,
)


def _payload():
    def cycle(loaded):
        return {
            "sample_count": 193,
            "cycle_duration_s": 1.0 / 60.0,
            "torque_nm": {
                "mean": 1.9217177651 if loaded else 2.0039728815e-7,
                "endpoint_relative_error": 2.8594529618e-5 if loaded else 1.0356795395e-5,
            },
            "phase_current_a": {
                "rms": [2.8284271383, 2.8284271092, 2.8284271325]
                if loaded
                else [4.317110181e-7, 3.367132768e-7, 4.408253246e-7],
                "rms_spread_relative": 1.0297148355e-8 if loaded else 0.2582892219,
                "peak_abs": 4.0000001479 if loaded else 3.4379620418e-6,
                "sum_max_abs": 1.7906432206e-6 if loaded else 2.9498194846e-6,
                "endpoint_relative_error_max": 7.1146067182e-8 if loaded else 0.0536869876,
            },
            "phase_voltage_v": {
                "rms": [52.5123398091, 52.4346800279, 52.5085666381]
                if loaded
                else [44.7733669687, 44.7382790755, 44.7349306127],
                "rms_spread_relative": 0.0014796512 if loaded else 0.0008589349,
            },
            "phase_flux_wb": {
                "rms": [0.1323286977, 0.1323267524, 0.1323272155]
                if loaded
                else [0.1148825173, 0.1148820431, 0.1148821876],
                "rms_spread_relative": 1.4700325238e-5 if loaded else 4.1282031462e-6,
                "endpoint_relative_error_max": 6.4304364638e-6 if loaded else 1.3733090477e-5,
            },
            "kinematics": {
                "right_endpoint_speed_angle_max_error_deg": 0.0,
                "trapezoid_speed_angle_max_error_deg": 0.46875,
                "first_step_angle_deg": 0.9375,
            },
            "power_w": {
                "electrical_cycle_mean": 385.4947362674 if loaded else -5.5e-7,
                "copper_loss_cycle_mean": 19.5360000157 if loaded else 2.3e-13,
                "mechanical_cycle_mean": 361.2919437842 if loaded else 3.8e-5,
                "balance_relative_error": 0.0121059823 if loaded else 69.31,
            },
        }

    return {
        "schema": "motor.transient-no-load-load-cycle.v1",
        "electrical_frequency_hz": 60.0,
        "mechanical_angle_span_deg": 180.0,
        "final_speed_rpm": 1800.0,
        "no_load": cycle(False),
        "loaded": cycle(True),
    }


def test_no_load_load_cycle_accepts_balanced_periodic_power_closed_pair():
    result = build_gate(_payload())
    assert result["status"] == "ok"
    assert result["metrics"]["inferred_pole_pairs"] == 2.0
    assert result["metrics"]["loaded_power_balance_relative_error"] < 0.013
    assert json.loads(motor_transient_no_load_load_cycle_gate(json.dumps(_payload())))["status"] == "ok"


def test_no_load_load_cycle_rejects_false_trapezoid_timing_and_power_fit():
    payload = copy.deepcopy(_payload())
    payload["loaded"]["kinematics"]["right_endpoint_speed_angle_max_error_deg"] = 0.46875
    payload["loaded"]["power_w"]["balance_relative_error"] = 0.12
    result = build_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["right_endpoint_speed_angle_convention"] is False
    assert result["checks"]["loaded_cycle_power_balance_within_3pct"] is False


def test_no_load_load_cycle_rejects_unbalanced_three_phase_current():
    payload = copy.deepcopy(_payload())
    payload["loaded"]["phase_current_a"]["sum_max_abs"] = 0.2
    payload["loaded"]["phase_current_a"]["rms"][2] *= 0.9
    result = build_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["loaded_four_amp_three_phase_balance"] is False
