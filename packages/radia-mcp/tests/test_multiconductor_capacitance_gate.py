from __future__ import annotations

import copy
import json

import numpy as np

from radia_mcp.radia_ngsolve.multiconductor_capacitance_gate import (
    multiconductor_capacitance_cross_formulation_gate,
)
from radia_mcp.radia_ngsolve.server import (
    multiconductor_capacitance_cross_formulation_gate as mcp_gate,
)


def good_summary() -> dict:
    first = np.array(
        [[4.0, -1.0, -0.5], [-1.0, 3.0, -0.3], [-0.5, -0.3, 2.0]]
    )
    last = 0.99 * first
    return {
        "capacitance_unit": "pF",
        "positions": [0.02, 0.28],
        "formulations": {
            "volume_fem": [first.tolist(), last.tolist()],
            "boundary_integral": [(0.95 * first).tolist(), (0.95 * last).tolist()],
        },
    }


def test_accepts_reciprocal_passive_multiconductor_matrix_families():
    result = multiconductor_capacitance_cross_formulation_gate(good_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["conductor_count"] == 3
    assert json.loads(mcp_gate(json.dumps(good_summary())))["status"] == "ok"


def test_rejects_nonreciprocity_and_opposite_position_sensitivity():
    bad = copy.deepcopy(good_summary())
    bad["formulations"]["boundary_integral"][0][1][0] *= 0.7
    bad["formulations"]["boundary_integral"][1] = (
        1.02 * np.asarray(bad["formulations"]["boundary_integral"][0])
    ).tolist()
    result = multiconductor_capacitance_cross_formulation_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_matrices_reciprocal"] is False
    assert result["checks"]["position_sensitivity_direction_agrees"] is False


def test_rejects_scalar_or_two_conductor_evidence():
    bad = good_summary()
    bad["formulations"]["volume_fem"] = [[[1.0, -0.2], [-0.2, 1.0]]] * 2
    bad["formulations"]["boundary_integral"] = [[[1.0, -0.2], [-0.2, 1.0]]] * 2
    try:
        multiconductor_capacitance_cross_formulation_gate(bad)
    except ValueError as exc:
        assert "at least three conductors" in str(exc)
    else:
        raise AssertionError("multi-conductor gate must reject a 2x2 substitute")
