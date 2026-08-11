from radia.ltspice.hysteresis_gate import hysteretic_inductor_cycle_gate
from radia.ltspice.mcp_server import hysteretic_inductor_cycle_gate as mcp_gate


ROWS = [
    {"cycle_index": 2, "current_peak_a": 1.19995, "total_energy_j": 0.00731825,
     "copper_energy_j": 0.00719528, "hysteresis_energy_j": 0.000122973,
     "flux_loop_energy_j": 0.000123010, "flux_closure_relative": 0.000422},
    {"cycle_index": 3, "current_peak_a": 1.19982, "total_energy_j": 0.00731842,
     "copper_energy_j": 0.00719528, "hysteresis_energy_j": 0.000123135,
     "flux_loop_energy_j": 0.000122998, "flux_closure_relative": 0.000305},
]


def test_hysteretic_cycle_gate_accepts_closed_repeatable_energy_loops():
    result = hysteretic_inductor_cycle_gate(
        ROWS,
        expected_current_peak_a=1.2,
        expected_copper_energy_j=0.0072,
        voltage_thd=0.04348,
    )
    assert result["status"] == "ok"
    assert result["metrics"]["steady_hysteresis_span_relative"] < 0.002


def test_hysteretic_cycle_gate_rejects_open_loop_and_missing_nonlinearity():
    bad = [dict(row) for row in ROWS]
    bad[1]["flux_closure_relative"] = 0.2
    bad[1]["flux_loop_energy_j"] *= -1
    result = hysteretic_inductor_cycle_gate(
        bad,
        expected_current_peak_a=1.2,
        expected_copper_energy_j=0.0072,
        voltage_thd=0.0,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["hysteresis_energy_is_positive"] is False
    assert result["checks"]["flux_loop_closes_each_steady_cycle"] is False
    assert result["checks"]["voltage_contains_nonlinear_harmonics"] is False


def test_hysteretic_cycle_mcp_dispatches_nested_rows():
    result = mcp_gate(ROWS, 1.2, 0.0072, 0.04348)
    assert result["status"] == "ok"
