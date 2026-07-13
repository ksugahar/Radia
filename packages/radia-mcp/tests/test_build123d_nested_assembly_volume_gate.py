import copy
import json

from radia_mcp.build123d.server import (
    build123d_nested_assembly_volume_gate,
    build123d_stud_wall_source_replay_gate,
)


SOURCE_SHA = "f1f707eaa299d9d858c30ccdf0f01f6e9374ac47e5195b2b935d3d1901fcdc00"
SOURCE_COMMIT = "fa8e93687c2e6069d0eae0e4b0b8ae128e33de1f"


def _shape(*, direct_children=2, parent_volume=0.0, leaf_sum=195_425_934.0):
    return {
        "type": "Compound",
        "direct_child_count": direct_children,
        "solid_count": 23,
        "face_count": 230,
        "edge_count": 552,
        "vertex_count": 368,
        "top_level_volume_mm3": parent_volume,
        "leaf_positive_volume_count": 23 if leaf_sum > 0.0 else 0,
        "leaf_volume_sum_mm3": leaf_sum,
        "is_valid": True,
    }


def _external():
    rows = []
    for index in (1, 2):
        rows.append(
            {
                "index": index,
                "import_command": {"returned": True, "exception": None},
                "imported": {
                    "volume_count": 23,
                    "positive_volume_count": 23,
                    "volume_sum_mm3": 195_426_000.0,
                    "surface_count": 230,
                    "curve_count": 552,
                    "vertex_count": 368,
                },
                "unite_command": {"returned": True, "exception": None},
                "united": {
                    "volume_count": 1,
                    "positive_volume_count": 1,
                    "volume_sum_mm3": 195_426_000.0,
                },
            }
        )
    return {
        "execution_mode": "headless_python_api_synchronous_commands",
        "headless_flags": ["-nographics", "-batch"],
        "gui_daemon_enabled": False,
        "replays": rows,
        "cad_handoff_ready": True,
        "mesh_attempted": False,
        "solver_ready": False,
        "process": {
            "acceptable": True,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
        },
    }


def generic_summary():
    native = _shape()
    return {
        "native": native,
        "same_kernel_roundtrips": {
            "step": _shape(direct_children=23),
            "brep": _shape(direct_children=0),
        },
        "external": _external(),
    }


def source_summary():
    assembly = _shape()
    x_wall = {"solid_count": 13}
    y_wall = {"solid_count": 10}
    replay = {
        "x_wall": x_wall,
        "y_wall": y_wall,
        "assembly": assembly,
        "x_wall_joint_names": ["end0", "inside0"],
        "y_wall_joint_names": ["end0", "inside0"],
        "wall_child_types": {
            "x_wall": ["Stud"] * 13,
            "y_wall": ["Stud"] * 10,
        },
    }
    return {
        "build": {
            "source": {
                "repository": "gumyr/build123d",
                "tag": "v0.10.0",
                "commit": SOURCE_COMMIT,
                "path": "examples/stud_wall.py",
                "sha256": SOURCE_SHA,
                "copy_sha256": SOURCE_SHA,
                "source_preserved": True,
                "display_stubbed_only": True,
            },
            "replays": [replay, copy.deepcopy(replay)],
            "files": {
                "step": {"sha256": "a" * 64},
                "brep": {"sha256": "b" * 64},
            },
        },
        "external": _external(),
        "nested_gate_status": "ok",
        "nested_gate_diagnosis": "valid_nested_compound_zero_parent_scalar",
        "timing_breakdown_s": {
            "source_replays": 1.0,
            "neutral_cad_roundtrips": 0.3,
            "external_cad_replays": 1.5,
            "mcp_evidence": 0.1,
        },
    }


def test_nested_assembly_volume_gate_accepts_zero_parent_with_positive_leaves():
    result = json.loads(build123d_nested_assembly_volume_gate(json.dumps(generic_summary())))
    assert result["status"] == "ok"
    assert result["diagnosis"] == "valid_nested_compound_zero_parent_scalar"
    assert result["parent_volume_zero_is_not_empty"] is True
    assert result["cad_handoff_ready"] is True
    assert result["solver_ready"] is False


def test_nested_assembly_volume_gate_rejects_empty_leaf_inventory():
    payload = generic_summary()
    payload["native"]["leaf_volume_sum_mm3"] = 0.0
    payload["native"]["leaf_positive_volume_count"] = 0
    result = json.loads(build123d_nested_assembly_volume_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["parent_volume_zero_is_not_empty"] is False
    assert result["cad_handoff_ready"] is False


def test_nested_assembly_volume_gate_rejects_external_volume_loss():
    payload = generic_summary()
    payload["external"]["replays"][1]["imported"]["volume_sum_mm3"] *= 0.9
    result = json.loads(build123d_nested_assembly_volume_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["external_leaf_volume_sum_matches"] is False
    assert result["checks"]["external_replays_are_deterministic"] is False


def test_stud_wall_source_replay_gate_accepts_exact_headless_replay():
    result = json.loads(build123d_stud_wall_source_replay_gate(json.dumps(source_summary())))
    assert result["status"] == "ok"
    assert result["cad_handoff_ready"] is True
    assert result["solver_ready"] is False


def test_stud_wall_source_replay_gate_rejects_source_drift():
    payload = source_summary()
    payload["build"]["source"]["sha256"] = "0" * 64
    result = json.loads(build123d_stud_wall_source_replay_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["tagged_upstream_source_identity_is_exact"] is False


def test_stud_wall_source_replay_gate_rejects_solver_ready_overclaim():
    payload = source_summary()
    payload["external"]["mesh_attempted"] = True
    payload["external"]["solver_ready"] = True
    result = json.loads(build123d_stud_wall_source_replay_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["cad_handoff_not_solver_readiness"] is False
