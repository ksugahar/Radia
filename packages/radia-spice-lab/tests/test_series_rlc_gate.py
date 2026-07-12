from __future__ import annotations

from copy import deepcopy

from ltspice_converter.mcp_server import series_rlc_complex_impedance_gate as mcp_gate
from ltspice_converter.series_rlc_gate import series_rlc_complex_impedance_gate


def good_summary() -> dict:
    return {
        "analysis": "ac",
        "source_orientation": "current_leaves_driven_node",
        "resistance_ohm": 100.0,
        "inductance_h": 0.01,
        "capacitance_f": 1.0e-6,
        "points": 401,
        "complex_relative_l2": {
            "input": 2.49e-7,
            "intermediate": 2.49e-7,
            "capacitor": 1.5e-12,
        },
        "maximum_input_pointwise_relative_error": 1.0e-5,
        "analytic_resonance_frequency_hz": 1591.5494309189535,
        "resonance_bracket_hz": [1548.816618912481, 1621.81009735893],
        "minimum_impedance_ohm": 100.00451283332686,
        "converted_and_reference_raw_are_equivalent": True,
        "converted_netlist_semantics_verified": True,
        "known_broken_schematic_is_rejected": True,
    }


def test_accepts_full_complex_series_rlc_recovery_evidence():
    result = series_rlc_complex_impedance_gate(good_summary())
    assert result["status"] == "ok"
    assert result["checks"]["three_full_complex_traces_match"] is True
    assert mcp_gate(good_summary())["status"] == "ok"


def test_rejects_magnitude_only_or_unverified_conversion_semantics():
    bad = deepcopy(good_summary())
    bad["complex_relative_l2"] = {"input_magnitude": 1.0e-12}
    bad["converted_netlist_semantics_verified"] = False
    try:
        series_rlc_complex_impedance_gate(bad)
    except ValueError as exc:
        assert "three traces" in str(exc)
    else:
        raise AssertionError("magnitude-only evidence must be rejected")

    bad = deepcopy(good_summary())
    bad["source_orientation"] = "unknown"
    bad["known_broken_schematic_is_rejected"] = False
    result = series_rlc_complex_impedance_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "source_orientation_preserved",
        "known_broken_schematic_is_rejected",
    }
