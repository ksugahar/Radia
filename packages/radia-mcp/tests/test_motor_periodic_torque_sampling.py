import json
import pytest

from radia_mcp.motor.periodic_torque_sampling_gate import periodic_torque_sampling_gate
from radia_mcp.motor.server import motor_periodic_torque_sampling_gate


def test_periodic_torque_sampling_gate_accepts_plot_endpoint_and_fft_exclusion():
    result = periodic_torque_sampling_gate(
        period_deg=180.0,
        sample_count=361,
        endpoint_included=True,
        spectrum_excludes_duplicate_endpoint=True,
        torque_min_Nm=1.678,
        torque_max_Nm=2.544,
        speed_rps=17.5,
        expected_step_deg=0.5,
    )
    assert result["status"] == "ok"
    assert result["interval_count"] == 360
    assert result["unique_spectrum_sample_count"] == 360
    assert result["step_deg"] == 0.5
    assert result["torque_ripple_peak_to_peak_Nm"] == pytest.approx(0.866)


def test_periodic_torque_sampling_gate_rejects_duplicate_fft_endpoint():
    result = periodic_torque_sampling_gate(
        period_deg=180.0,
        sample_count=361,
        endpoint_included=True,
        spectrum_excludes_duplicate_endpoint=False,
        torque_min_Nm=1.678,
        torque_max_Nm=2.544,
        speed_rps=17.5,
        expected_step_deg=0.5,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["duplicate_endpoint_excluded_from_spectrum"] is False


def test_motor_periodic_torque_sampling_mcp_tool_dispatches_json():
    result = json.loads(motor_periodic_torque_sampling_gate(
        180.0, 361, True, True, 1.678, 2.544, 17.5, 0.5
    ))
    assert result["status"] == "ok"
    assert result["unique_spectrum_sample_count"] == 360
