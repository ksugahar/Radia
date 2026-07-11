import copy
import json

from radia_mcp.radia_ngsolve.conductor_frequency_gate import (
    homogenized_bundle_impedance_comparison_gate,
)
from radia_mcp.radia_ngsolve.server import (
    homogenized_bundle_impedance_comparison_gate as mcp_gate,
)


ROWS = [
    {
        "model_role": "homogenized",
        "frequency_hz": 50000.0,
        "current_a_complex": [1.0, 0.0],
        "voltage_v_complex": [8.36452174059618, 46.64697852963744],
        "resistance_ohm": 8.36452174059618,
        "inductance_h": 0.0001484819442658662,
        "element_count": 12953,
        "solve_time_s": 0.40457660000538453,
    },
    {
        "model_role": "explicit_reference",
        "frequency_hz": 50000.0,
        "current_a_complex": [1.0, 0.0],
        "voltage_v_complex": [8.579232608733186, 46.79852460556663],
        "resistance_ohm": 8.579232608733186,
        "inductance_h": 0.0001489643304076724,
        "element_count": 90604,
        "solve_time_s": 3.9139793000067584,
    },
]


def test_homogenized_bundle_accepts_live_shape_accuracy_and_speedup():
    result = homogenized_bundle_impedance_comparison_gate(ROWS)
    assert result["status"] == "ok"
    assert result["metrics"]["resistance_relative_error"] < 0.026
    assert result["metrics"]["solve_time_speedup"] > 9.0
    assert json.loads(mcp_gate(ROWS))["status"] == "ok"


def test_homogenized_bundle_rejects_active_or_inaccurate_fast_surrogate():
    bad = copy.deepcopy(ROWS)
    bad[0]["voltage_v_complex"][0] = -1.0
    bad[0]["resistance_ohm"] = -1.0
    result = homogenized_bundle_impedance_comparison_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["passive_positive_resistance"] is False
    assert result["checks"]["homogenized_resistance_accurate"] is False
