import copy
import json
import math

from radia_mcp.radia_ngsolve.fem_bem_capstone_gate import fem_bem_capstone_suite_gate
from radia_mcp.radia_ngsolve.server import fem_bem_capstone_suite_gate as mcp_gate


def _payload():
    details = {
        91: {"meshVertices": 4, "meshElements": 1},
        92: {"error": 0.0},
        93: {"error": 0.0},
        94: {"hasLaplaceSL": True},
        95: {"hasHelmholtzSL": True},
        96: {"lowFrequencyValue": {"real": 1.0 / (4.0 * math.pi), "imag": 1.0e-10}},
        97: {"meanPotential": 0.6},
        98: {"meanAmplitude": 0.5},
        99: {"traceRows": 4},
        100: {"h1Error": 1.0e-16, "hcurlError": 8.0e-16},
    }
    return {
        "cases": [
            {"id": f"GYP-{index:03d}", "passed": True, "failures": [], "details": details[index]}
            for index in range(91, 101)
        ],
        "capabilities": {
            "ok": True,
            "mesh_vertices": 4,
            "mesh_elements": 1,
            "h1_dofs": 4,
            "hcurl_dofs": 6,
            "has_bem": True,
            "has_laplace_sl": True,
            "has_helmholtz_sl": True,
        },
    }


def test_capstone_gate_accepts_complete_first_order_suite_and_dispatches():
    payload = _payload()
    assert fem_bem_capstone_suite_gate(payload)["status"] == "ok"
    assert json.loads(mcp_gate(payload))["status"] == "ok"


def test_capstone_gate_rejects_trace_and_low_frequency_drift():
    payload = copy.deepcopy(_payload())
    payload["cases"][5]["details"]["lowFrequencyValue"]["real"] *= 1.1
    payload["cases"][8]["details"]["traceRows"] = 3
    result = fem_bem_capstone_suite_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["low_frequency_kernel_reaches_laplace_limit"] is False
    assert result["checks"]["p1_trace_has_four_boundary_rows"] is False
