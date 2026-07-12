from __future__ import annotations

import copy
import json

from radia_mcp.build123d.server import (
    build123d_repeated_cavity_dual_api_gate,
    build123d_repeated_cavity_source_replay_gate,
)


def public_summary() -> dict:
    topology = {
        "solid_count": 1,
        "shell_count": 1,
        "face_count": 40,
        "edge_count": 90,
        "vertex_count": 60,
        "bbox_extent": [48.0, 16.0, 11.4],
    }
    base = {
        "volume": 3200.0,
        "area": 5600.0,
        "bbox_extent": topology["bbox_extent"],
        **{key: topology[key] for key in ("solid_count", "shell_count", "face_count", "edge_count", "vertex_count")},
    }
    native = {}
    for mode in ("builder", "algebra"):
        native[mode] = {**base, "self_roundtrip": copy.deepcopy(base)}
    external = []
    for source_mode in ("builder", "algebra"):
        for import_mode in ("noheal", "heal"):
            external.append({
                "source_mode": source_mode,
                "import_mode": import_mode,
                "volume": 3200.0,
                "area": 5600.0,
                "volume_count": 1,
                "surface_count": 40,
                "curve_count": 90,
                "vertex_count": 60,
                "volume_relative_error": 1.0e-10,
                "area_relative_error": 1.0e-10,
                "bbox_extent_relative_error": 0.0,
            })
    return {
        "official_expected_volume": 3200.0,
        "feature_contract": {
            "repeated_top_feature_count": 12,
            "internal_support_count": 5,
            "internal_cavity_present": True,
        },
        "topology_contract": topology,
        "native": native,
        "external_imports": external,
    }


def call_public(summary: dict) -> dict:
    return json.loads(build123d_repeated_cavity_dual_api_gate(json.dumps(summary)))


def source_summary(public_gate: dict) -> dict:
    return {
        "source_kind": "upstream_native_examples",
        "source_commit": "a" * 40,
        "source_files_preserved": True,
        "sources": [
            {"mode": "builder", "sha256": "b" * 64},
            {"mode": "algebra", "sha256": "c" * 64},
        ],
        "viewer_stub_only": True,
        "cubit_batch_entry_mode": "exec_compile_wrapper",
        "pip_count_x": 6,
        "pip_count_y": 2,
        "repeated_top_feature_count": 12,
        "internal_cavity_present": True,
        "official_volume_assertion_reproduced": True,
        "step_artifacts": [
            {"mode": "builder", "sha256": "d" * 64, "bytes": 1000},
            {"mode": "algebra", "sha256": "e" * 64, "bytes": 1000},
        ],
        "external_process": {
            "execution_mode": "python_api_headless",
            "headless_flags": ["-nographics", "-batch"],
            "gui_daemon_enabled": False,
            "acceptable": True,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
        },
        "timing_breakdown_s": {"builder": 0.1, "algebra": 0.1, "external": 0.1, "review": 0.1},
        "public_gate": public_gate,
    }


def test_accepts_dual_api_cavity_and_four_imports():
    result = call_public(public_summary())
    assert result["status"] == "ok"


def test_rejects_cavity_topology_loss_even_when_volume_matches():
    summary = public_summary()
    summary["external_imports"][0]["surface_count"] -= 2
    result = call_public(summary)
    assert result["status"] == "needs_attention"
    assert "external_imports_preserve_topology" in result["issues"]


def test_source_gate_accepts_immutable_dual_replay():
    result = json.loads(
        build123d_repeated_cavity_source_replay_gate(
            json.dumps(source_summary(call_public(public_summary())))
        )
    )
    assert result["status"] == "ok"


def test_source_gate_rejects_geometry_edit_and_stale_external_result():
    summary = source_summary(call_public(public_summary()))
    summary["viewer_stub_only"] = False
    summary["cubit_batch_entry_mode"] = "direct_multiline_python"
    summary["repeated_top_feature_count"] = 10
    summary["external_process"]["result_artifact_fresh"] = False
    result = json.loads(build123d_repeated_cavity_source_replay_gate(json.dumps(summary)))
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "viewer_suppression_is_stub_only",
        "cubit_batch_python_uses_exec_compile_entry",
        "source_parameter_and_feature_contract_bound",
        "fresh_result_and_no_owned_process_leak",
    }
