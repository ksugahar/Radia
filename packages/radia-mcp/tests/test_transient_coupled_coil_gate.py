import json

from radia_mcp.radia_ngsolve.server import transient_coupled_coil_response_gate as mcp_gate
from radia_mcp.radia_ngsolve.transient_coupled_coil_gate import transient_coupled_coil_response_gate


def _synthetic_history():
    times = [index * 0.001 for index in range(16)]
    primary = [0.0, 0.2, 0.7, 1.0, 0.6, 0.0, -0.6, -1.0, -0.7, 0.0, 0.7, 1.0, 0.5, 0.0, -0.5, -0.9]
    secondary = [0.0]
    for old, before, after in zip(secondary, primary, primary[1:]):
        secondary.append(0.82 * old - 0.31 * (after - before))
    return times, primary, secondary


def test_transient_coupled_coil_gate_accepts_passive_first_order_history_and_dispatches():
    times, primary, secondary = _synthetic_history()
    result = transient_coupled_coil_response_gate(
        times,
        primary,
        secondary,
        secondary_resistance_ohm=0.4,
        secondary_turns=80,
    )
    assert result["status"] == "ok"
    assert result["metrics"]["coupling_polarity"] == -1
    assert result["metrics"]["maximum_relative_residual"] < 1.0e-12
    assert json.loads(mcp_gate(times, primary, secondary, 0.4, 80))["status"] == "ok"


def test_transient_coupled_coil_gate_rejects_active_memory_and_waveform_corruption():
    times, primary, secondary = _synthetic_history()
    corrupted = list(secondary)
    for index in range(1, len(corrupted)):
        corrupted[index] = 1.03 * corrupted[index - 1] - 0.31 * (primary[index] - primary[index - 1])
    corrupted[9] += 0.2
    result = transient_coupled_coil_response_gate(
        times,
        primary,
        corrupted,
        secondary_resistance_ohm=0.4,
        secondary_turns=80,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["passive_first_order_memory"] is False
    assert result["checks"]["first_order_response_matches"] is False
