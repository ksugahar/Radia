from __future__ import annotations

import copy
import json

from radia_mcp.cubit.pyramid_export_replay_gate import (
    cubit_pyramid_mixed_export_gate as public_gate,
    cubit_pyramid_source_plugin_replay_gate as source_gate,
)
from radia_mcp.cubit.server import (
    cubit_pyramid_mixed_export_gate,
    cubit_pyramid_source_plugin_replay_gate,
)


def _public() -> dict:
    return {
        "element_counts": {"hex": 8, "pyramid": 4, "tet": 115, "wedge": 0},
        "connectivity_sizes": {
            "hex": [8],
            "pyramid": [5],
            "tet": [4],
            "wedge": [],
        },
        "quality": {
            "hex": {"scaled_jacobian": {"available": True, "count": 8, "min": 1.0}},
            "tet": {"scaled_jacobian": {"available": True, "count": 115, "min": 0.36}},
            "pyramid": {"scaled_jacobian": {"available": False, "count": 0}},
        },
        "gmsh_header": {"status": "ok", "mesh_format": "4.1", "binary": False},
        "gmsh_independent": {
            "element_counts": {"hex": 8, "pyramid": 4, "tet": 115, "other_3d": 0},
            "volume_by_type": {"hex": 1000.0, "pyramid": 117.8, "tet": 882.2},
        },
        "nastran_cards": {"CHEXA": 8, "CPYRAM": 4, "CTETRA": 115},
        "gmsh_companions": {"geo": True, "geo_opt": True, "msh_opt": True},
        "gmsh_cad_volume_relative_error": 0.0,
    }


def _source() -> dict:
    return {
        "source_native_example": True,
        "source_extension": ".py",
        "source_sha256": "a" * 64,
        "execution_mode": "combined_journal_headless",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "ordinary_python_plugin_command_rejected": True,
        "combined_journal_plugin_command_succeeded": True,
        "compatibility_transforms": {
            "cubit_runtime": "2024.3 Python init -> 2025.12 combined journal startup",
            "legacy_nastran_function": "cubit_mesh_export.export_3D_Nastran -> export jmag_nastran",
        },
        "version": "2025.12",
        "process": {
            "exit_code": 3,
            "error_lines": ["ERROR: known startup diagnostic"],
            "unexpected_error_lines": [],
            "known_headless_diagnostics_only": True,
            "acceptable": True,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
        },
        "source_commands_present": True,
        "public_gate": {"policy": "cubit_pyramid_mixed_export_gate_v1", "status": "ok"},
        "public_negative_status": "needs_attention",
        "timing_breakdown_s": {"a": 1.0, "b": 2.0, "c": 3.0, "d": 4.0},
    }


def test_public_gate_accepts_explicit_pyramid_exports() -> None:
    result = public_gate(_public())
    assert result["status"] == "ok"
    assert result["metrics"]["nastran_card_counts"]["CPYRAM"] == 4
    assert json.loads(cubit_pyramid_mixed_export_gate(_public()))["status"] == "ok"


def test_public_gate_rejects_dropped_pyramids_and_volume_drift() -> None:
    summary = copy.deepcopy(_public())
    summary["nastran_cards"]["CPYRAM"] = 0
    summary["gmsh_cad_volume_relative_error"] = 0.1
    result = public_gate(summary)
    assert result["status"] == "needs_attention"
    assert "nastran_preserves_explicit_cpyram" in result["issues"]
    assert "gmsh_volume_matches_cad" in result["issues"]


def test_source_gate_accepts_combined_journal_plugin_startup() -> None:
    result = source_gate(_source())
    assert result["status"] == "ok"
    assert json.loads(cubit_pyramid_source_plugin_replay_gate(_source()))["status"] == "ok"


def test_source_gate_rejects_ordinary_python_and_unexpected_exit() -> None:
    summary = copy.deepcopy(_source())
    summary["combined_journal_plugin_command_succeeded"] = False
    summary["process"]["unexpected_error_lines"] = ["ERROR: export failed"]
    summary["process"]["acceptable"] = False
    result = source_gate(summary)
    assert result["status"] == "needs_attention"
    assert "plugin_registration_path_was_demonstrated" in result["issues"]
    assert "known_nonzero_headless_exit_classified" in result["issues"]
