import copy
import json

from radia_mcp.radia_ngsolve.leakage_inductance_closure_gate import (
    leakage_inductance_closure_gate,
)
from radia_mcp.radia_ngsolve.server import leakage_inductance_closure_gate as mcp_gate


def _summary():
    return {
        "turns": [2.0, 1.0],
        "compensated_currents_A": [1.0, -2.0],
        "compensated_flux_linkage_Wb": [2.0, -0.11],
        "matrix_H": [[5.0, 1.5], [1.49, 0.8]],
        "replay_max_abs_Wb": 0.0,
    }


def test_leakage_closure_accepts_energy_and_matrix_routes():
    result = leakage_inductance_closure_gate(_summary())
    assert result["status"] == "ok"
    assert result["observables"]["leakage_inductance_direct_H"] == 2.22
    assert abs(result["observables"]["direct_matrix_relative_error"]) < 1e-14
    assert 0.0 < result["observables"]["short_circuit_inductance_H"] < 5.0


def test_leakage_closure_rejects_uncompensated_and_inconsistent_direct_route():
    bad = copy.deepcopy(_summary())
    bad["compensated_currents_A"][1] = -1.5
    bad["compensated_flux_linkage_Wb"] = [9.0, 9.0]
    result = leakage_inductance_closure_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["compensated_ampere_turns_close"] is False
    assert result["checks"]["direct_and_matrix_leakage_close"] is False


def test_leakage_closure_mcp_dispatches_and_handles_bad_shape():
    result = json.loads(mcp_gate(json.dumps(_summary())))
    assert result["status"] == "ok"
    bad = json.loads(mcp_gate(json.dumps({"matrix_H": [[1.0]]})))
    assert bad["status"] == "invalid_input"
