from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.manual_auto_mixed_mesh_gate import (
    manual_auto_mixed_mesh_preservation_gate as gate,
)
from radia_mcp.radia_ngsolve.server import manual_auto_mixed_mesh_preservation_gate


def _row(total, vertices, manual, automatic, manual_size, keep, partial):
    return {
        "total_elements": total,
        "vertices": vertices,
        "region_elements": {
            "manual_region": manual,
            "automatic_region": automatic,
        },
        "auto_mesh_size_m": {
            "manual_region": manual_size,
            "automatic_region": 0.00091,
        },
        "keep_existing_mesh": keep,
        "partial_mesh_before": partial,
        "physics_result_after": False,
        "source_preserved": True,
        "temporary_work_copy": True,
        "pass_marker": True,
        "owned_processes_after": 0,
    }


def _summary():
    auto = _row(958, 517, 34, 341, 0.00091, False, False)
    mixed = _row(832, 462, 16, 235, 0.0, True, True)
    return {
        "model_contract": {
            "dimension": 2,
            "same_two_part_geometry": True,
            "manual_region_element_family": "quadrilateral",
            "automatic_region_element_family": "triangle",
            "manual_automatic_interface": "conformal",
        },
        "archived_reference": {
            "automatic_only": _row(938, 500, 34, 323, 0.00091, False, False),
            "manual_plus_automatic": _row(844, 470, 16, 247, 0.0, True, True),
        },
        "fresh_replays": {
            "automatic_only": [copy.deepcopy(auto), copy.deepcopy(auto)],
            "manual_plus_automatic": [copy.deepcopy(mixed), copy.deepcopy(mixed)],
        },
        "classification": "manual_region_preserved_auto_region_version_drift_recorded",
        "solver_ready": False,
    }


def test_accepts_exact_manual_preservation_with_bounded_auto_drift():
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["manual_region_live_elements"] == 16
    assert result["metrics"]["mixed_automatic_region_relative_drift"] < 0.05


def test_rejects_manual_region_replacement_and_missing_keep_flag():
    summary = _summary()
    summary["fresh_replays"]["manual_plus_automatic"][1]["region_elements"][
        "manual_region"
    ] = 17
    summary["fresh_replays"]["manual_plus_automatic"][0][
        "keep_existing_mesh"
    ] = False
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "two_deterministic_mixed_replays" in result["issues"]
    assert "manual_region_is_preserved_exactly" in result["issues"]
    assert "keep_existing_mesh_and_partial_mesh_are_explicit" in result["issues"]


def test_rejects_unbounded_automatic_mesher_drift():
    summary = _summary()
    for row in summary["fresh_replays"]["manual_plus_automatic"]:
        row["region_elements"]["automatic_region"] = 100
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "automatic_region_version_drift_is_bounded" in result["issues"]


def test_mcp_wrapper_returns_structured_invalid_input():
    result = json.loads(manual_auto_mixed_mesh_preservation_gate("{}"))
    assert result["status"] == "invalid_input"
