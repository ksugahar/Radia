from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v17 import _public_v17, _source_v17


def _public_v18():
    reference, measured = _public_v17()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("6" if index == 0 else "7") * 64
            row["nested_assembly_location_transform_composition_identity"] = {
                "assembly_generation": "assembly-52",
                "location_tree_assembly_generation": "assembly-52",
                "composition_assembly_generation": "assembly-52",
                "transform_order_generation": "transform-order-52",
                "location_transform_order_generation": "transform-order-52",
                "child_paths": ["root/frame", "root/frame/shaft"],
                "location_child_paths": ["root/frame", "root/frame/shaft"],
                "local_transform_sha256": ["1" * 64, "2" * 64],
                "composed_transform_sha256": ["3" * 64, "4" * 64],
                "resolved_composed_transform_sha256": ["3" * 64, "4" * 64],
                "composition_order": "parent_then_child",
                "resolved_composition_order": "parent_then_child",
                "location_tree_sha256": digest,
                "resolved_location_tree_sha256": digest,
            }
            row["boolean_retained_face_name_history_refine_identity"] = {
                "boolean_generation": "boolean-52",
                "refine_generation": "refine-52",
                "retained_name_boolean_generation": "boolean-52",
                "retained_name_refine_generation": "refine-52",
                "topology_history_refine_generation": "refine-52",
                "retained_face_names": ["inlet", "outlet"],
                "resolved_face_names": ["inlet", "outlet"],
                "retained_face_ids": [41, 42],
                "resolved_face_ids": [41, 42],
                "topology_history_sha256": digest,
                "resolved_topology_history_sha256": digest,
            }
    return reference, measured


def _source_v18():
    row = _source_v17()
    identity = row["replay_identity"]
    identity["step_occurrence_color_material_inheritance_identity"] = {
        "step_import_generation": "step-import-52",
        "assembly_generation": "assembly-52",
        "occurrence_metadata_import_generation": "step-import-52",
        "occurrence_metadata_assembly_generation": "assembly-52",
        "color_inheritance_assembly_generation": "assembly-52",
        "material_inheritance_assembly_generation": "assembly-52",
        "occurrence_ids": ["root/frame", "root/frame/shaft"],
        "metadata_occurrence_ids": ["root/frame", "root/frame/shaft"],
        "parent_occurrence_ids": ["root", "root/frame"],
        "metadata_parent_occurrence_ids": ["root", "root/frame"],
        "inherited_colors_rgb": [[0.7, 0.7, 0.7], [0.2, 0.2, 0.2]],
        "imported_colors_rgb": [[0.7, 0.7, 0.7], [0.2, 0.2, 0.2]],
        "inherited_material_names": ["aluminum", "steel"],
        "imported_material_names": ["aluminum", "steel"],
        "occurrence_metadata_sha256": "8" * 64,
        "imported_occurrence_metadata_sha256": "8" * 64,
    }
    identity["stl_tolerance_model_length_unit_generation_identity"] = {
        "model_length_unit_generation": "model-unit-52",
        "tessellation_model_unit_generation": "model-unit-52",
        "tolerance_conversion_model_unit_generation": "model-unit-52",
        "model_length_unit": "mm",
        "chordal_tolerance_length_unit": "mm",
        "angular_tolerance_unit": "deg",
        "chordal_tolerance_value": 0.01,
        "chordal_tolerance_si_m": 1.0e-5,
        "tessellator_chordal_tolerance_si_m": 1.0e-5,
        "angular_tolerance_value": 5.0,
        "angular_tolerance_rad": 0.08726646259971647,
        "tessellator_angular_tolerance_rad": 0.08726646259971647,
        "tolerance_contract_sha256": "9" * 64,
        "tessellator_tolerance_contract_sha256": "9" * 64,
    }
    return row


def test_v18_positive_public_and_source_lineage():
    reference, measured = _public_v18()
    public = _public_result(reference, measured)
    assert public["status"] == "ok"
    assert public["checks"][
        "nested_assembly_locations_use_current_transform_composition_order"
    ]
    assert public["checks"][
        "boolean_retained_face_names_follow_post_refine_topology_history"
    ]
    source = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v18()))
    )
    assert source["status"] == "ok"
    assert source["checks"][
        "step_occurrence_inheritance_uses_current_assembly_generation"
    ]
    assert source["checks"][
        "stl_tolerances_use_current_model_length_unit_generation"
    ]


def test_v18_public_nested_assembly_location_transform_composition_order_mismatch():
    reference, measured = _public_v18()
    measured["external_cad"][0][
        "nested_assembly_location_transform_composition_identity"
    ].update(
        {
            "location_transform_order_generation": "transform-order-51",
            "resolved_composed_transform_sha256": ["4" * 64, "3" * 64],
            "resolved_composition_order": "child_then_parent",
            "resolved_location_tree_sha256": "5" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "nested_assembly_locations_use_current_transform_composition_order"
    ] is False


def test_v18_public_boolean_retained_face_name_history_refine_generation_mismatch():
    reference, measured = _public_v18()
    measured["external_cad"][0][
        "boolean_retained_face_name_history_refine_identity"
    ].update(
        {
            "retained_name_refine_generation": "refine-51",
            "topology_history_refine_generation": "refine-51",
            "resolved_face_names": ["outlet", "inlet"],
            "resolved_face_ids": [42, 41],
            "resolved_topology_history_sha256": "5" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "boolean_retained_face_names_follow_post_refine_topology_history"
    ] is False


def test_v18_source_step_occurrence_color_material_inheritance_generation_mismatch():
    row = _source_v18()
    row["replay_identity"][
        "step_occurrence_color_material_inheritance_identity"
    ].update(
        {
            "occurrence_metadata_assembly_generation": "assembly-51",
            "color_inheritance_assembly_generation": "assembly-51",
            "material_inheritance_assembly_generation": "assembly-51",
            "metadata_parent_occurrence_ids": ["root/old", "root/old/frame"],
            "imported_colors_rgb": [[0.2, 0.2, 0.2], [0.7, 0.7, 0.7]],
            "imported_material_names": ["steel", "aluminum"],
            "imported_occurrence_metadata_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "step_occurrence_inheritance_uses_current_assembly_generation"
    ] is False


def test_v18_source_stl_chordal_angular_tolerance_model_unit_generation_mismatch():
    row = _source_v18()
    row["replay_identity"][
        "stl_tolerance_model_length_unit_generation_identity"
    ].update(
        {
            "tolerance_conversion_model_unit_generation": "model-unit-51",
            "chordal_tolerance_length_unit": "m",
            "tessellator_chordal_tolerance_si_m": 0.01,
            "angular_tolerance_unit": "rad",
            "tessellator_angular_tolerance_rad": 5.0,
            "tessellator_tolerance_contract_sha256": "5" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"][
        "stl_tolerances_use_current_model_length_unit_generation"
    ] is False
