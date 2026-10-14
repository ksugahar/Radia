from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.eddy_loss_formulation_gate import (
    alternate_eddy_loss_formulation_gate,
)
from radia_mcp.radia_ngsolve.server import alternate_eddy_loss_formulation_gate as mcp_gate


def _summary() -> dict:
    return {
        "frequency_hz": 60.0,
        "combine_requested": False,
        "volume": {
            "dataset_id": "resolved",
            "solution_id": "volume_solution",
            "selection_dimension": 3,
            "native_loss_w": 8.0,
            "builtin_integral_w": 8.0,
            "jdot_e_integral_w": 8.0,
        },
        "surface": {
            "dataset_id": "reduced",
            "solution_id": "surface_solution",
            "selection_dimension": 2,
            "native_loss_w": 35.0,
            "builtin_integral_w": 35.0,
            "rerun_relative_change": 2.0e-8,
        },
    }


def test_accepts_independently_closed_alternate_formulations():
    result = alternate_eddy_loss_formulation_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["surface_to_volume_loss_ratio"] == 35.0 / 8.0


def test_rejects_cross_dataset_addition_and_wrong_selection_dimension():
    bad = copy.deepcopy(_summary())
    bad["combine_requested"] = True
    bad["surface"]["selection_dimension"] = 3
    result = alternate_eddy_loss_formulation_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "alternate_formulations_not_added",
        "surface_selection_is_boundary",
    }


def test_mcp_wrapper_reports_invalid_json_without_raising():
    result = json.loads(mcp_gate("{"))
    assert result["status"] == "invalid_input"
