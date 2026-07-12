from __future__ import annotations

import copy
import json

from radia_mcp.cubit.server import (
    cubit_symmetric_swept_mixed_mesh_gate,
    cubit_symmetric_swept_source_replay_gate,
)
from radia_mcp.cubit.symmetric_mixed_gate import (
    cubit_symmetric_swept_mixed_mesh_gate as public_gate,
    cubit_symmetric_swept_source_replay_gate as source_gate,
)


def _public() -> dict:
    return {
        "geometry": {"volume_count": 5, "cad_total_volume": 2718.47, "left_right_aggregate_relative_error": 1.0e-16},
        "element_counts": {"hex": 390, "pyramid": 39, "tet": 3518, "wedge": 0},
        "connectivity_sizes": {"hex": [8], "pyramid": [5], "tet": [4]},
        "quality": {"hex_scaled_jacobian": {"min": 0.469}, "tet_scaled_jacobian": {"min": 0.246}},
        "gmsh": {
            "mesh_format": "4.1",
            "binary": False,
            "element_counts": {"hex": 390, "pyramid": 39, "tet": 3518, "other_3d": 0},
            "volume_by_type": {"hex": 1359.0, "pyramid": 10.0, "tet": 1349.0},
            "cad_volume_relative_error": 9.9e-5,
        },
    }


def _source() -> dict:
    return {
        "source_native_journal": True,
        "promotion": "mirrored_quad_sections_to_swept_hex_pyramid_tet",
        "source_sha256": "a" * 64,
        "execution_mode": "headless_combined_journal_then_python_inventory",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "process": {"acceptable": True, "result_artifact_fresh": True, "unexpected_error_lines": [], "owned_processes_remaining": 0, "known_headless_diagnostics_only": True, "exit_code": 3},
        "timing_breakdown_s": {"a": 1.0, "b": 2.0, "c": 3.0, "d": 4.0},
        "public_gate": public_gate(_public()),
        "public_negative_status": "needs_attention",
    }


def test_public_gate_accepts_symmetric_mixed_transition() -> None:
    assert public_gate(_public())["status"] == "ok"


def test_public_gate_rejects_missing_pyramids_and_symmetry_drift() -> None:
    summary = copy.deepcopy(_public())
    summary["element_counts"]["pyramid"] = 0
    summary["geometry"]["left_right_aggregate_relative_error"] = 0.1
    assert public_gate(summary)["status"] == "needs_attention"


def test_source_gate_accepts_classified_headless_replay() -> None:
    assert source_gate(_source())["status"] == "ok"


def test_mcp_wrappers_serialize_results() -> None:
    assert json.loads(cubit_symmetric_swept_mixed_mesh_gate(_public()))["status"] == "ok"
    assert json.loads(cubit_symmetric_swept_source_replay_gate(_source()))["status"] == "ok"
