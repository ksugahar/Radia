from __future__ import annotations

import copy
import json

import pytest

from radia_mcp.build123d.server import (
    build123d_jointed_assembly_source_replay_gate,
    build123d_mass_property_crosscheck,
)


def _box(name: str, size: tuple[float, float, float]) -> dict:
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


def _public() -> tuple[list[dict], dict[str, list[dict]]]:
    reference = [
        _box("frame", (2.0, 3.0, 4.0)),
        _box("insert", (1.0, 2.0, 5.0)),
    ]
    measured = copy.deepcopy(reference)
    child_revisions = {"frame": "child-frame-5", "insert": "child-insert-3"}
    for rows in (reference, measured):
        for row in rows:
            name = row["name"]
            revision = f"brep-{name}-42"
            digest = ("d" if name == "frame" else "e") * 64
            row["brep_identity"] = {"revision": revision, "sha256": digest}
            row["mass_property_identity"] = {
                "brep_revision": revision,
                "brep_sha256": digest,
            }
            row["assembly_identity"] = {
                "generation": "assembly-generation-42",
                "child_revisions": copy.deepcopy(child_revisions),
            }
    return reference, {"external_cad": measured}


def _public_result(reference: list[dict], measured: dict[str, list[dict]]) -> dict:
    return json.loads(
        build123d_mass_property_crosscheck(
            json.dumps(reference),
            json.dumps(measured),
            rtol=1.0e-10,
            bbox_atol=1.0e-10,
        )
    )


def _source() -> dict:
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
            "source_replay_started_utc": "2026-07-16T02:00:00Z",
            "cad_artifacts": [
                {
                    "name": "assembly.step",
                    "sha256": "c" * 64,
                    "fresh": True,
                    "source_commit": "b" * 40,
                    "export_completed_utc": "2026-07-16T02:00:01Z",
                }
            ],
            "external_kernel": {
                "name": "OCCT",
                "claimed_version": "7.8.1",
                "replay_versions": ["7.8.1", "7.8.1"],
                "claimed_session_generation": "occt-session-42",
                "replay_sessions": [
                    {
                        "session_generation": "occt-session-42",
                        "process_start_utc": "2026-07-16T01:59:00Z",
                    },
                    {
                        "session_generation": "occt-session-42",
                        "process_start_utc": "2026-07-16T01:59:00Z",
                    },
                ],
            },
        },
    }


def test_v8_positive_revision_and_session_contracts() -> None:
    reference, measured = _public()
    assert _public_result(reference, measured)["status"] == "ok"
    source = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source()))
    )
    assert source["status"] == "ok"


@pytest.mark.parametrize(
    "case_id",
    [
        "v8_public_mass_properties_older_than_brep",
        "v8_public_assembly_child_revision_mix",
    ],
)
def test_generalization_v8_public(case_id: str) -> None:
    reference, measured = _public()
    rows = measured["external_cad"]
    if case_id == "v8_public_mass_properties_older_than_brep":
        rows[0]["mass_property_identity"].update(
            {"brep_revision": "brep-frame-41", "brep_sha256": "f" * 64}
        )
        expected = "mass_properties_bind_current_brep_revision"
    else:
        rows[1]["assembly_identity"]["child_revisions"][
            "insert"
        ] = "child-insert-2"
        expected = "assembly_children_match_reference_revision_map"
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v8_source_export_precedes_source_replay",
        "v8_source_kernel_session_restarts_between_replays",
    ],
)
def test_generalization_v8_source(case_id: str) -> None:
    row = _source()
    identity = row["replay_identity"]
    if case_id == "v8_source_export_precedes_source_replay":
        identity["cad_artifacts"][0][
            "export_completed_utc"
        ] = "2026-07-16T01:59:59Z"
        expected = "neutral_cad_export_follows_source_replay"
    else:
        identity["external_kernel"]["replay_sessions"][1].update(
            {
                "session_generation": "occt-session-43",
                "process_start_utc": "2026-07-16T02:00:30Z",
            }
        )
        expected = "external_kernel_session_generation_is_continuous"
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(row))
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False
