import copy

from radia.ltspice.bridge_rectifier_gate import bridge_rectifier_gate
from radia.ltspice.mcp_server import bridge_rectifier_gate as mcp_bridge_rectifier_gate


def good_summary():
    return {
        "topology": "single_phase_full_wave_bridge",
        "units": {"voltage": "V", "current": "A", "frequency": "Hz", "time": "s"},
        "analysis_window_s": [80.0e-6, 100.0e-6],
        "input_frequency_hz": 100000.0,
        "ripple_frequency_hz": 200000.0,
        "vout_average_v": 4.8114,
        "vout_min_v": 0.00889,
        "load_average_a": 0.048114,
        "capacitor_average_a": -1.63e-8,
        "diode_average_sum_a": 0.096228,
        "diagonal_pair_a_waveform_relative_error": 1.1e-11,
        "diagonal_pair_b_waveform_relative_error": 5.6e-12,
        "alternate_pair_overlap_fraction": 0.0,
        "kcl_max_relative_error": 9.2e-8,
        "tolerances": {
            "ripple_frequency_ratio_relative_error": 1.0e-6,
            "diode_sum_to_twice_load_relative_error": 1.0e-5,
            "diagonal_pair_waveform_relative_error": 1.0e-8,
            "alternate_pair_overlap_fraction": 1.0e-4,
            "kcl_max_relative_error": 1.0e-5,
            "capacitor_average_to_load_relative_error": 1.0e-5,
        },
    }


def test_accepts_frequency_pair_and_current_closure():
    result = bridge_rectifier_gate(good_summary())
    assert result["status"] == "ok"
    assert result["checks"]["ripple_frequency_is_twice_input"] is True
    assert mcp_bridge_rectifier_gate(good_summary())["status"] == "ok"


def test_rejects_half_wave_frequency_and_bad_current_sum():
    bad = copy.deepcopy(good_summary())
    bad["ripple_frequency_hz"] = bad["input_frequency_hz"]
    bad["diode_average_sum_a"] = bad["load_average_a"]
    result = bridge_rectifier_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["ripple_frequency_is_twice_input"] is False
    assert result["checks"]["four_diode_average_sum_is_twice_load"] is False


def test_rejects_overlapping_diagonal_pairs():
    bad = good_summary()
    bad["alternate_pair_overlap_fraction"] = 0.2
    result = bridge_rectifier_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["opposite_pairs_do_not_overlap"] is False
