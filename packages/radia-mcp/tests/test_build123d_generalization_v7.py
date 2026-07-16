import copy
import json

import pytest

from radia_mcp.build123d.server import (
    build123d_jointed_assembly_source_replay_gate,
    build123d_mass_property_crosscheck,
)


def _box(name, size):
    x, y, z = size
    return {
        "name": name,
        "type": "Solid",
        "is_valid": True,
        "volume": x * y * z,
        "area": 2.0 * (x * y + y * z + x * z),
        "faces": 6,
        "edges": 12,
        "vertices": 8,
        "solids": 1,
        "bounding_box": {
            "min": [0.0, 0.0, 0.0],
            "max": [x, y, z],
            "center": [x / 2.0, y / 2.0, z / 2.0],
            "size": [x, y, z],
            "diagonal": (x * x + y * y + z * z) ** 0.5,
        },
    }


def _public():
    reference = [_box("frame", (2.0, 3.0, 4.0)), _box("insert", (1.0, 2.0, 5.0))]
    return reference, {"external_cad": copy.deepcopy(reference)}


def _source():
    return {
        "source_kind": "upstream_source_native_example_with_display_stub_only",
        "source_sha256": "a" * 64,
        "source_url": "https://example.invalid/project/blob/v0.10.0/examples/model.py",
        "source_preserved": True,
        "display_stubbed_only": True,
        "components": [
            {"name": "frame", "joint_names": ["frame_joint"]},
            {"name": "insert", "joint_names": ["insert_joint"]},
        ],
        "joint_connections": [
            {"from": "frame_joint", "to": "insert_joint", "kind": "rigid"}
        ],
        "external_execution": {
            "mode": "python_api_headless_synchronous_commands",
            "headless_flags": ["-nographics", "-batch"],
            "gui_daemon_enabled": False,
            "result_artifact_fresh": True,
            "owned_processes_remaining": 0,
        },
        "diagnosis_gate_status": "ok",
        "diagnosis": "component_solid_closure_loss",
        "solver_ready": False,
        "timing_breakdown_s": {
            "source_replay": 0.2,
            "neutral_cad_export": 0.1,
            "external_replay": 0.3,
            "identity_validation": 0.05,
        },
        "replay_identity": {
            "source_commit": "b" * 40,
            "replayed_source_commit": "b" * 40,
            "cad_artifacts": [
                {
                    "name": "assembly.step",
                    "sha256": "c" * 64,
                    "fresh": True,
                    "source_commit": "b" * 40,
                }
            ],
            "external_kernel": {
                "name": "OCCT",
                "claimed_version": "7.8.1",
                "replay_versions": ["7.8.1", "7.8.1"],
            },
        },
    }


def _public_result(reference, measured):
    return json.loads(
        build123d_mass_property_crosscheck(
            json.dumps(reference),
            json.dumps(measured),
            rtol=1.0e-10,
            bbox_atol=1.0e-10,
        )
    )


@pytest.mark.parametrize(
    "case_id",
    [
        "v7_public_consistent_unit_scale_error",
        "v7_public_missing_body_mass_compensation",
    ],
)
def test_generalization_v7_public(case_id):
    reference, measured = _public()
    rows = measured["external_cad"]
    if case_id == "v7_public_consistent_unit_scale_error":
        for row in rows:
            row["volume"] *= 1.0e9
            row["area"] *= 1.0e6
            for key in ("min", "max", "center", "size"):
                row["bounding_box"][key] = [value * 1.0e3 for value in row["bounding_box"][key]]
            row["bounding_box"]["diagonal"] *= 1.0e3
    else:
        removed = rows.pop()
        rows[0]["volume"] += removed["volume"]
        rows[0]["area"] += removed["area"]
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"


def test_source_identity_contract_accepts_consistent_replays():
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(_source())))
    assert result["status"] == "ok"
    assert result["warnings"] == []


@pytest.mark.parametrize("identity", [[], {"cad_artifacts": ["not-a-row"]}])
def test_source_identity_contract_rejects_malformed_evidence(identity):
    row = _source()
    row["replay_identity"] = identity
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"


@pytest.mark.parametrize(
    ("case_id", "failed_check"),
    [
        (
            "v7_source_stale_commit_fresh_cad_digest",
            "neutral_cad_artifacts_bind_current_source_commit",
        ),
        (
            "v7_source_external_kernel_version_drift",
            "external_kernel_versions_are_replay_invariant",
        ),
    ],
)
def test_generalization_v7_source(case_id, failed_check):
    row = _source()
    if case_id == "v7_source_stale_commit_fresh_cad_digest":
        row["replay_identity"]["cad_artifacts"][0]["source_commit"] = "d" * 40
    else:
        row["replay_identity"]["external_kernel"]["replay_versions"][1] = "7.9.0"
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][failed_check] is False
