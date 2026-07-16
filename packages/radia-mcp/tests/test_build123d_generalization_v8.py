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
            row["mass_property_frame_identity"] = {
                "frame_id": "assembly-global-frame-42",
                "transform_generation": "assembly-transform-42",
            }
            row["topology_identity"] = {
                "brep_revision": revision,
                "brep_sha256": digest,
                "face_adjacency_sha256": ("1" if name == "frame" else "2") * 64,
            }
            row["compound_volume_identity"] = {
                "topology_kind": "physical_union_solid",
                "reported_volume_basis": "physical_union",
                "overlap_volume": 0.0,
                "topology_generation": "compound-generation-42",
                "volume_generation": "compound-generation-42",
            }
            row["placement_transform_identity"] = {
                "center_of_mass_frame": "assembly-global-frame-42",
                "center_of_mass_transform_generation": "assembly-transform-42",
                "final_placement_transform_generation": "assembly-transform-42",
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
            "topology_replay_identity": {
                "source_topology_sha256": "3" * 64,
                "imports": [
                    {"mode": "noheal", "topology_sha256": "3" * 64},
                    {"mode": "heal", "topology_sha256": "3" * 64},
                ],
            },
            "unit_conversion_identity": {
                "source_geometry_unit": "mm",
                "target_geometry_unit": "m",
                "length_scale_to_target": 0.001,
                "external_measurement_stage": "after_unit_conversion",
                "external_volume_unit": "m^3",
                "declared_volume_scale_to_target": 1.0,
            },
            "boolean_clean_identity": {
                "boolean_result_sha256": "4" * 64,
                "shape_clean_input_sha256": "4" * 64,
                "cleaned_topology_sha256": "5" * 64,
                "export_topology_sha256": "5" * 64,
                "shape_generation": "shape-generation-42",
                "export_shape_generation": "shape-generation-42",
            },
            "tessellation_identity": {
                "shape_generation": "shape-generation-42",
                "tolerance_shape_generation": "shape-generation-42",
                "linear_deflection": 0.01,
                "angular_tolerance_rad": 0.1,
                "tessellation_generation": "tessellation-generation-42",
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


@pytest.mark.parametrize(
    "case_id",
    [
        "v9_public_center_of_mass_reference_frame_mismatch",
        "v9_public_face_adjacency_before_fillet",
    ],
)
def test_generalization_v9_public(case_id: str) -> None:
    reference, measured = _public()
    rows = measured["external_cad"]
    if case_id == "v9_public_center_of_mass_reference_frame_mismatch":
        rows[0]["mass_property_frame_identity"].update(
            {
                "frame_id": "component-local-frame",
                "transform_generation": "component-transform-9",
            }
        )
        expected = "mass_property_centers_share_reference_frames"
    else:
        rows[0]["topology_identity"].update(
            {
                "brep_revision": "brep-frame-41",
                "face_adjacency_sha256": "4" * 64,
            }
        )
        expected = "face_adjacency_matches_current_brep_revision"
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v9_source_import_heal_topology_digest_omitted",
        "v9_source_external_volume_before_unit_conversion",
    ],
)
def test_generalization_v9_source(case_id: str) -> None:
    row = _source()
    identity = row["replay_identity"]
    if case_id == "v9_source_import_heal_topology_digest_omitted":
        identity["topology_replay_identity"]["imports"][1].pop("topology_sha256")
        expected = "heal_and_noheal_imports_record_topology_identity"
    else:
        identity["unit_conversion_identity"].update(
            {
                "external_measurement_stage": "before_unit_conversion",
                "external_volume_unit": "mm^3",
                "declared_volume_scale_to_target": 1.0e-9,
            }
        )
        expected = "external_volume_is_measured_after_unit_conversion"
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(row))
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v10_public_compound_overlap_double_count_volume",
        "v10_public_center_of_mass_frame_volume_frame_mismatch",
    ],
)
def test_generalization_v10_public(case_id: str) -> None:
    reference, measured = _public()
    row = measured["external_cad"][0]
    if case_id == "v10_public_compound_overlap_double_count_volume":
        row["compound_volume_identity"].update(
            {"reported_volume_basis": "child_volume_sum", "overlap_volume": 0.125}
        )
        expected = "compound_volume_uses_physical_union_not_child_sum"
    else:
        row["placement_transform_identity"][
            "center_of_mass_transform_generation"
        ] = "assembly-transform-41"
        expected = "center_of_mass_uses_final_placement_transform"
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    "case_id",
    [
        "v10_source_boolean_result_digest_before_shape_clean",
        "v10_source_tessellation_tolerance_previous_shape",
    ],
)
def test_generalization_v10_source(case_id: str) -> None:
    row = _source()
    identity = row["replay_identity"]
    if case_id == "v10_source_boolean_result_digest_before_shape_clean":
        identity["boolean_clean_identity"].update(
            {
                "shape_clean_input_sha256": "3" * 64,
                "export_shape_generation": "shape-generation-43",
            }
        )
        expected = "boolean_export_follows_shape_clean_identity"
    else:
        identity["tessellation_identity"].update(
            {
                "tolerance_shape_generation": "shape-generation-41",
                "tessellation_generation": "tessellation-generation-41",
            }
        )
        expected = "tessellation_tolerances_belong_to_current_shape"
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(row))
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False
