from __future__ import annotations

import copy
import json

from radia_mcp.build123d.faceted_edit_gate import (
    build123d_faceted_edit_portability_gate as public_gate,
    build123d_faceted_source_replay_gate as source_gate,
)
from radia_mcp.build123d.server import (
    build123d_faceted_edit_portability_gate,
    build123d_faceted_source_replay_gate,
)


def _shape(volume=15163.2, area=9256.3, faces=918):
    return {
        "volume_mm3": volume,
        "area_mm2": area,
        "solid_count": 1,
        "face_count": faces,
        "edge_count": 1399,
        "vertex_count": 476,
        "is_valid": True,
    }


def _public() -> dict:
    native_runs = []
    for index in (1, 2):
        native_runs.append(
            {
                "run_index": index,
                "native": _shape(),
                "imported_stl_before_edit": _shape(15161.3, 9238.7, 980),
                "self_roundtrip": _shape(),
            }
        )
    external = [
        {
            "run_index": run,
            "import_mode": mode,
            "volume_count": 1,
            "positive_volume_count": 1,
            "surface_count": 918,
            "curve_count": 1399,
            "vertex_count": 476,
            "native_volume_relative_error": 2.0e-6,
            "native_area_relative_error": 3.0e-6,
        }
        for run in (1, 2)
        for mode in ("heal", "noheal")
    ]
    return {
        "native_runs": native_runs,
        "external_imports": external,
        "mesh_evidence": {
            "connectivity_orders": [4],
            "cad_volume_relative_error": 2.6e-4,
            "min_scaled_jacobian": 8.6e-5,
            "tet_count": 33582,
            "gmsh_companions": {"geo": True, "geo_opt": True, "msh_opt": True},
        },
        "solver_ready_claimed": False,
        "mesh_quality_classification": "diagnostic_only",
    }


def _source(public_result: dict) -> dict:
    return {
        "source_kind": "upstream_native_example_tag_v0.10.0",
        "upstream_tag": "v0.10.0",
        "upstream_commit": "a" * 40,
        "source_sha256": "b" * 64,
        "dependent_asset_bound": True,
        "stl_sha256": "c" * 64,
        "source_preserved": True,
        "stl_preserved": True,
        "execution_mode": "exact_source_with_viewer_stub_only",
        "viewer_stub_only": True,
        "build123d_version": "0.10.0",
        "source_run_count": 2,
        "mesh_imported_solid_edited": True,
        "face_count_delta": -62,
        "external_cad_version": "2025.12",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "external_process": {
            "acceptable": True,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
            "unexpected_error_lines": [],
        },
        "public_gate": public_result,
        "public_negative_status": "needs_attention",
        "timing_breakdown_s": {"a": 1.0, "b": 2.0, "c": 3.0, "d": 4.0},
    }


def test_public_gate_accepts_portable_but_not_solver_ready_case() -> None:
    result = public_gate(_public())
    assert result["status"] == "ok"
    assert result["solver_ready"] is False
    assert json.loads(build123d_faceted_edit_portability_gate(json.dumps(_public())))[
        "status"
    ] == "ok"


def test_public_gate_rejects_solver_ready_overclaim_and_volume_drift() -> None:
    summary = copy.deepcopy(_public())
    summary["solver_ready_claimed"] = True
    summary["mesh_quality_classification"] = "solver_ready"
    summary["external_imports"][0]["native_volume_relative_error"] = 0.01
    result = public_gate(summary)
    assert result["status"] == "needs_attention"
    assert "solver_ready_claim_matches_quality" in result["issues"]
    assert "external_mass_properties_close" in result["issues"]


def test_source_gate_accepts_bound_stl_and_headless_replay() -> None:
    public_result = public_gate(_public())
    result = source_gate(_source(public_result))
    assert result["status"] == "ok"
    assert json.loads(
        build123d_faceted_source_replay_gate(json.dumps(_source(public_result)))
    )["status"] == "ok"


def test_source_gate_rejects_missing_asset_and_viewer_execution() -> None:
    summary = _source(public_gate(_public()))
    summary["dependent_asset_bound"] = False
    summary["stl_sha256"] = ""
    summary["viewer_stub_only"] = False
    result = source_gate(summary)
    assert result["status"] == "needs_attention"
    assert "dependent_stl_identity_recorded" in result["issues"]
    assert "exact_source_viewer_stub_only" in result["issues"]
