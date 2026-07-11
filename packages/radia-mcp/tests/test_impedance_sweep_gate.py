import json

from radia_mcp.radia_ngsolve.impedance_sweep_gate import multiport_impedance_sweep_gate
from radia_mcp.radia_ngsolve.server import multiport_impedance_sweep_gate as mcp_gate


FREQUENCIES = [1.0e5, 1.0e6, 1.0e7, 1.0e8, 1.0e9]


def test_multiport_impedance_sweep_accepts_common_passive_profiles():
    result = multiport_impedance_sweep_gate(
        [FREQUENCIES, FREQUENCIES],
        [[10.0, 5.0, 2.0, 1.0, 0.5], [8.0, 4.0, 2.0, 1.0, 0.25]],
        [[-100.0, -20.0, 0.0, 10.0, 50.0], [-80.0, -10.0, 1.0, 8.0, 20.0]],
        port_ids=["supply", "load"],
    )
    assert result["status"] == "ok"
    assert result["metrics"]["port_count"] == 2


def test_multiport_impedance_sweep_rejects_grid_drift_and_active_real_part():
    shifted = FREQUENCIES.copy()
    shifted[2] *= 1.1
    result = multiport_impedance_sweep_gate(
        [FREQUENCIES, shifted],
        [[1.0] * 5, [1.0, 1.0, -0.5, 1.0, 1.0]],
        [[1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 2.0, 3.0, 4.0, 5.0]],
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["common_frequency_grid"] is False
    assert result["checks"]["impedance_positive_real_with_tolerance"] is False


def test_multiport_impedance_sweep_mcp_dispatches_nested_rows():
    result = json.loads(mcp_gate(
        [FREQUENCIES],
        [[1.0, 2.0, 3.0, 4.0, 5.0]],
        [[5.0, 4.0, 3.0, 2.0, 1.0]],
        ["input"],
    ))
    assert result["status"] == "ok"
