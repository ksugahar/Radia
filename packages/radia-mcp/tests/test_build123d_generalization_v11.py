from __future__ import annotations

import json

import pytest

from radia_mcp.build123d.server import (
    build123d_jointed_assembly_source_replay_gate,
)
from test_build123d_generalization_v8 import _public, _public_result, _source


def _public_v11():
    reference, measured = _public()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            final_digest = ("6" if index == 0 else "7") * 64
            row["shape_healing_identity"] = {
                "pre_heal_brep_sha256": ("4" if index == 0 else "5") * 64,
                "healed_brep_sha256": final_digest,
                "final_brep_sha256": final_digest,
                "mass_property_brep_sha256": final_digest,
                "final_shape_generation": "shape-generation-43",
                "mass_property_shape_generation": "shape-generation-43",
            }
            frame = row["mass_property_frame_identity"]["frame_id"]
            transform = row["placement_transform_identity"][
                "final_placement_transform_generation"
            ]
            row["inertia_tensor_identity"] = {
                "tensor_frame_id": frame,
                "center_of_mass_frame_id": frame,
                "tensor_transform_generation": transform,
                "final_placement_transform_generation": transform,
                "mirror_transform_applied": True,
                "mirror_transform_determinant": -1.0,
                "tensor_basis_handedness": "right_handed",
            }
    return reference, measured


def _source_v11():
    row = _source()
    identity = row["replay_identity"]
    identity["step_export_tolerance_identity"] = {
        "kernel_session_generation": "occt-session-42",
        "tolerance_kernel_session_generation": "occt-session-42",
        "shape_generation": "shape-generation-42",
        "tolerance_shape_generation": "shape-generation-42",
        "sewing_tolerance": 1.0e-7,
        "brep_tolerance": 1.0e-7,
        "export_artifact_sha256": "8" * 64,
    }
    identity["assembly_replacement_identity"] = {
        "assembly_generation": "assembly-generation-43",
        "components": [
            {
                "slot_id": "insert-slot",
                "replacement_generation": "replacement-generation-43",
                "removed_instance_uuid": "instance-uuid-42",
                "current_instance_uuid": "instance-uuid-43",
                "removed_shape_sha256": "9" * 64,
                "current_shape_sha256": "a" * 64,
                "placement_shape_sha256": "a" * 64,
                "placement_assembly_generation": "assembly-generation-43",
            }
        ],
    }
    return row


def test_v11_positive_healing_frame_and_session_contracts() -> None:
    reference, measured = _public_v11()
    assert _public_result(reference, measured)["status"] == "ok"
    source = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v11()))
    )
    assert source["status"] == "ok"


@pytest.mark.parametrize(
    ("case_id", "expected"),
    [
        (
            "v11_public_shape_healed_after_mass_properties_stale",
            "mass_properties_follow_final_healed_brep",
        ),
        (
            "v11_public_mirrored_part_inertia_tensor_frame_mismatch",
            "mirrored_inertia_tensor_uses_final_global_frame",
        ),
    ],
)
def test_generalization_v11_public(case_id: str, expected: str) -> None:
    reference, measured = _public_v11()
    row = measured["external_cad"][0]
    if case_id == "v11_public_shape_healed_after_mass_properties_stale":
        row["shape_healing_identity"].update(
            {
                "mass_property_brep_sha256": "4" * 64,
                "mass_property_shape_generation": "shape-generation-42",
            }
        )
    else:
        row["inertia_tensor_identity"].update(
            {
                "tensor_frame_id": "source-part-frame",
                "tensor_transform_generation": "source-transform-42",
            }
        )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False


@pytest.mark.parametrize(
    ("case_id", "expected"),
    [
        (
            "v11_source_step_export_tolerance_previous_kernel_session",
            "step_export_tolerances_belong_to_current_kernel_and_shape",
        ),
        (
            "v11_source_assembly_instance_uuid_reused_after_replace",
            "replacement_rotates_instance_uuid_and_rebinds_placement",
        ),
    ],
)
def test_generalization_v11_source(case_id: str, expected: str) -> None:
    row = _source_v11()
    identity = row["replay_identity"]
    if case_id == "v11_source_step_export_tolerance_previous_kernel_session":
        identity["step_export_tolerance_identity"].update(
            {
                "tolerance_kernel_session_generation": "occt-session-41",
                "tolerance_shape_generation": "shape-generation-41",
            }
        )
    else:
        component = identity["assembly_replacement_identity"]["components"][0]
        component.update(
            {
                "current_instance_uuid": component["removed_instance_uuid"],
                "placement_shape_sha256": component["removed_shape_sha256"],
            }
        )
    result = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(row))
    )
    assert result["status"] == "needs_attention"
    assert result["checks"][expected] is False
