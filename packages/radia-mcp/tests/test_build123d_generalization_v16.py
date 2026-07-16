from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v15 import _public_v15, _source_v15


def _public_v16():
    reference, measured = _public_v15()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["mass_inertia_reference_frame_placement_identity"] = {
                "shape_generation": "shape-50",
                "mass_property_shape_generation": "shape-50",
                "placement_generation": "placement-50",
                "mass_property_placement_generation": "placement-50",
                "mass_reference_frame": "world",
                "inertia_reference_frame": "world",
                "center_of_mass_reference_frame": "world",
                "placed_shape_sha256": digest,
                "mass_property_shape_sha256": digest,
            }
            row["loft_wire_correspondence_seam_identity"] = {
                "loft_generation": "loft-50",
                "section_wire_loft_generation": "loft-50",
                "seam_normalization_generation": "seam-50",
                "wire_correspondence_seam_generation": "seam-50",
                "section_wire_ids": [11, 12, 13],
                "loft_section_wire_ids": [11, 12, 13],
                "wire_correspondence_sha256": digest,
                "loft_wire_correspondence_sha256": digest,
            }
    return reference, measured


def _source_v16():
    row = _source_v15()
    identity = row["replay_identity"]
    identity["step_assembly_child_parent_unit_identity"] = {
        "step_import_generation": "step-import-50",
        "child_placement_import_generation": "step-import-50",
        "parent_placement_import_generation": "step-import-50",
        "assembly_length_unit": "mm",
        "child_placement_length_unit": "mm",
        "parent_placement_length_unit": "mm",
        "assembly_scale_to_m": 1.0e-3,
        "child_placement_scale_to_m": 1.0e-3,
        "parent_placement_scale_to_m": 1.0e-3,
        "assembly_placement_sha256": "3" * 64,
        "resolved_assembly_placement_sha256": "3" * 64,
    }
    identity["selector_normal_world_frame_identity"] = {
        "shape_generation": "selector-shape-50",
        "selector_shape_generation": "selector-shape-50",
        "placement_generation": "selector-placement-50",
        "selector_placement_generation": "selector-placement-50",
        "normal_predicate_frame": "world",
        "evaluated_normal_frame": "world",
        "selected_face_ids": [21, 22],
        "resolved_face_ids": [21, 22],
        "normal_table_sha256": "4" * 64,
        "evaluated_normal_table_sha256": "4" * 64,
    }
    return row


def test_v16_positive_public_and_source_lineage():
    reference, measured = _public_v16()
    assert _public_result(reference, measured)["status"] == "ok"
    source = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v16()))
    )
    assert source["status"] == "ok"


def test_v16_public_mass_inertia_tensor_reference_frame_placement_generation_mismatch():
    reference, measured = _public_v16()
    measured["external_cad"][0][
        "mass_inertia_reference_frame_placement_identity"
    ].update(
        {
            "mass_property_placement_generation": "placement-49",
            "inertia_reference_frame": "local",
            "mass_property_shape_sha256": "5" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["mass_inertia_uses_final_world_placement_frame"] is False


def test_v16_public_loft_wire_correspondence_seam_normalization_generation_mismatch():
    reference, measured = _public_v16()
    measured["external_cad"][0]["loft_wire_correspondence_seam_identity"].update(
        {
            "wire_correspondence_seam_generation": "seam-49",
            "loft_section_wire_ids": [11, 13, 12],
            "loft_wire_correspondence_sha256": "5" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "loft_sections_use_current_seam_normalized_correspondence"
        ]
        is False
    )


def test_v16_source_step_assembly_child_parent_length_unit_metadata_mismatch():
    row = _source_v16()
    row["replay_identity"]["step_assembly_child_parent_unit_identity"].update(
        {
            "child_placement_import_generation": "step-import-49",
            "child_placement_length_unit": "m",
            "child_placement_scale_to_m": 1.0,
            "resolved_assembly_placement_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "step_assembly_child_parent_placements_share_length_unit"
        ]
        is False
    )


def test_v16_source_selector_normal_predicate_preplacement_frame_mismatch():
    row = _source_v16()
    row["replay_identity"]["selector_normal_world_frame_identity"].update(
        {
            "selector_placement_generation": "selector-placement-49",
            "evaluated_normal_frame": "local",
            "resolved_face_ids": [20, 21],
            "evaluated_normal_table_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["selector_normals_use_final_world_placement_frame"] is False
