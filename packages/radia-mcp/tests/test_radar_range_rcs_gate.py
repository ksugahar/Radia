import json

from radia_mcp.radia_ngsolve.radar_range_rcs_gate import radar_range_rcs_profile_gate
from radia_mcp.radia_ngsolve.server import radar_range_rcs_profile_gate as mcp_gate


FREQUENCY_HZ = [75.0e9 + 1.0e6 * index for index in range(1001)]
GOOD = {
    "target_range_m": 5.0,
    "radar_peak_range_m": 5.006897966438633,
    "radar_peak_rcs_m2": 4.807180419774478,
    "generalized_peak_range_m": 5.006897966438633,
    "generalized_peak_rcs_m2": 4.807116452413412,
    "analytic_peak_rcs_m2": 4.981296338340539,
    "profile_relative_l2": 1.2036785780859605e-5,
}


def test_range_rcs_gate_accepts_resolved_consistent_profiles():
    result = radar_range_rcs_profile_gate(FREQUENCY_HZ, **GOOD)
    assert result["status"] == "ok"
    assert 0.149 < result["metrics"]["physical_range_resolution_m"] < 0.151
    assert result["metrics"]["unambiguous_range_m"] > 149.0


def test_range_rcs_gate_rejects_stale_target_and_analytic_peak():
    result = radar_range_rcs_profile_gate(
        FREQUENCY_HZ,
        **{**GOOD, "target_range_m": 5.5, "analytic_peak_rcs_m2": 9.0},
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["radar_peak_localizes_target"] is False
    assert result["checks"]["radar_peak_matches_analytic_reference"] is False


def test_range_rcs_gate_rejects_frequency_jitter_and_method_drift():
    frequency = FREQUENCY_HZ.copy()
    frequency[500] += 0.25e6
    result = radar_range_rcs_profile_gate(
        frequency,
        **{**GOOD, "generalized_peak_rcs_m2": 3.0, "profile_relative_l2": 0.2},
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["frequency_grid_equispaced"] is False
    assert result["checks"]["two_reconstruction_profiles_agree"] is False
    assert result["checks"]["two_reconstruction_peaks_agree"] is False


def test_range_rcs_mcp_dispatches_physical_resolution_gate():
    result = json.loads(mcp_gate(FREQUENCY_HZ, **GOOD))
    assert result["status"] == "ok"
    assert result["policy"] == "radar_range_rcs_profile_gate_v1"
