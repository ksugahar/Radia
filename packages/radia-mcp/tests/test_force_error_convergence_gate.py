import json

from radia_mcp.radia_ngsolve.force_error_convergence_gate import (
    dual_formulation_force_error_convergence_gate,
)
from radia_mcp.radia_ngsolve.server import (
    dual_formulation_force_error_convergence_gate as mcp_gate,
)


ROWS = [
    {
        "id": "boundary",
        "refinement_levels": [0.5, 1.0, 2.0],
        "observable_ids": ["body", "probe"],
        "error_profiles": [
            [0.2971919, 0.0260253, 0.0050251],
            [0.0029537, 0.0020911, 0.0019800],
        ],
    },
    {
        "id": "volume",
        "refinement_levels": [0.5, 1, 2, 3, 4, 5],
        "observable_ids": ["body", "probe"],
        "error_profiles": [
            [0.0355826, 0.0606084, 0.0314470, 0.0068225, 0.0082182, 0.0082078],
            [0.0028072, 0.0026989, 0.0024433, 0.0026169, 0.0024112, 0.0024115],
        ],
    },
]


def test_dual_force_error_gate_accepts_nonmonotone_convergence_envelopes():
    result = dual_formulation_force_error_convergence_gate(ROWS, reference_force=1.0)
    assert result["status"] == "ok"
    assert result["metrics"]["nonmonotone_profile_count"] == 2
    assert result["metrics"]["maximum_final_relative_error"] < 0.01


def test_dual_force_error_gate_rejects_bad_final_and_lost_plateau():
    bad = json.loads(json.dumps(ROWS))
    bad[1]["error_profiles"][0][-1] = 0.2
    result = dual_formulation_force_error_convergence_gate(bad, reference_force=1.0)
    assert result["status"] == "needs_attention"
    assert result["checks"]["final_force_errors_meet_reference_band"] is False
    assert result["checks"]["final_errors_remain_near_best"] is False
    assert result["checks"]["long_sweeps_reach_tail_plateau"] is False


def test_dual_force_error_mcp_dispatches_nested_formulations():
    result = json.loads(mcp_gate(ROWS, 1.0))
    assert result["status"] == "ok"
