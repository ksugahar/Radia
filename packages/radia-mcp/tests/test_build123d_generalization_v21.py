from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v20 import _public_v20, _source_v20


def _public_v21():
    reference, measured = _public_v20()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["boolean_result_solid_orientation_location_label_generation_identity"] = {
                "boolean_generation": "boolean-71",
                "result_boolean_generation": "boolean-71",
                "orientation_boolean_generation": "boolean-71",
                "location_boolean_generation": "boolean-71",
                "label_boolean_generation": "boolean-71",
                "operand_generations": ["operand-a-71", "operand-b-71"],
                "result_operand_generations": ["operand-a-71", "operand-b-71"],
                "solid_orientation": "forward",
                "resolved_solid_orientation": "forward",
                "semantic_labels": ["body", "interface"],
                "resolved_semantic_labels": ["body", "interface"],
                "result_location_sha256": "3" * 64,
                "resolved_result_location_sha256": "3" * 64,
                "boolean_result_sha256": digest,
                "resolved_boolean_result_sha256": digest,
            }
            row["tessellation_chord_angle_unit_location_generation_identity"] = {
                "shape_generation": "shape-71",
                "tessellation_shape_generation": "shape-71",
                "tessellation_generation": "tessellation-71",
                "metric_tessellation_generation": "tessellation-71",
                "location_tessellation_generation": "tessellation-71",
                "chord_tolerance": 1.0e-4,
                "evaluated_chord_tolerance": 1.0e-4,
                "angular_tolerance_deg": 12.0,
                "evaluated_angular_tolerance_deg": 12.0,
                "length_unit": "m",
                "evaluated_length_unit": "m",
                "object_location_sha256": "4" * 64,
                "evaluated_object_location_sha256": "4" * 64,
                "tessellation_sha256": digest,
                "evaluated_tessellation_sha256": digest,
            }
    return reference, measured


def _source_v21():
    row = _source_v20()
    identity = row["replay_identity"]
    identity["step_assembly_product_id_color_location_generation_identity"] = {
        "step_export_generation": "step-export-71",
        "decoder_step_export_generation": "step-export-71",
        "assembly_generation": "assembly-71",
        "product_id_assembly_generation": "assembly-71",
        "color_assembly_generation": "assembly-71",
        "hierarchy_assembly_generation": "assembly-71",
        "location_assembly_generation": "assembly-71",
        "product_ids": ["root", "rotor", "housing"],
        "decoded_product_ids": ["root", "rotor", "housing"],
        "parent_product_ids": ["", "root", "root"],
        "decoded_parent_product_ids": ["", "root", "root"],
        "colors_rgb": [[180, 180, 180], [220, 40, 40], [80, 100, 140]],
        "decoded_colors_rgb": [[180, 180, 180], [220, 40, 40], [80, 100, 140]],
        "component_location_sha256": ["5" * 64, "6" * 64, "7" * 64],
        "decoded_component_location_sha256": ["5" * 64, "6" * 64, "7" * 64],
        "assembly_metadata_sha256": "8" * 64,
        "decoded_assembly_metadata_sha256": "8" * 64,
    }
    identity["sketch_constraint_entity_id_solver_order_generation_identity"] = {
        "sketch_generation": "sketch-71",
        "entity_table_sketch_generation": "sketch-71",
        "constraint_table_sketch_generation": "sketch-71",
        "solver_order_sketch_generation": "sketch-71",
        "entity_ids": [101, 102, 103],
        "replay_entity_ids": [101, 102, 103],
        "constraint_ids": [201, 202],
        "solver_constraint_order": [201, 202],
        "replay_solver_constraint_order": [201, 202],
        "constraint_entity_ids": [[101, 102], [102, 103]],
        "replay_constraint_entity_ids": [[101, 102], [102, 103]],
        "entity_table_sha256": "9" * 64,
        "replay_entity_table_sha256": "9" * 64,
        "constraint_table_sha256": "a" * 64,
        "replay_constraint_table_sha256": "a" * 64,
    }
    return row


def test_v21_positive_public_and_source_identity():
    reference, measured = _public_v21()
    assert _public_result(reference, measured)["status"] == "ok"
    source = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v21()))
    )
    assert source["status"] == "ok"


def test_v21_public_boolean_result_solid_orientation_location_label_generation_mismatch():
    reference, measured = _public_v21()
    measured["external_cad"][0][
        "boolean_result_solid_orientation_location_label_generation_identity"
    ].update(
        {
            "orientation_boolean_generation": "boolean-70",
            "location_boolean_generation": "boolean-69",
            "label_boolean_generation": "boolean-68",
            "result_operand_generations": ["operand-a-70", "operand-b-70"],
            "resolved_solid_orientation": "reversed",
            "resolved_semantic_labels": ["old_body", "old_interface"],
            "resolved_result_location_sha256": "b" * 64,
            "resolved_boolean_result_sha256": "c" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "boolean_result_uses_current_orientation_location_and_labels"
    ]


def test_v21_public_tessellation_chord_angle_unit_location_generation_mismatch():
    reference, measured = _public_v21()
    measured["external_cad"][0][
        "tessellation_chord_angle_unit_location_generation_identity"
    ].update(
        {
            "tessellation_shape_generation": "shape-70",
            "metric_tessellation_generation": "tessellation-70",
            "location_tessellation_generation": "tessellation-69",
            "evaluated_chord_tolerance": 0.1,
            "evaluated_angular_tolerance_deg": 25.0,
            "evaluated_length_unit": "mm",
            "evaluated_object_location_sha256": "d" * 64,
            "evaluated_tessellation_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "tessellation_uses_current_tolerances_units_and_object_location"
    ]


def test_v21_source_step_assembly_product_id_color_location_generation_mismatch():
    row = _source_v21()
    row["replay_identity"][
        "step_assembly_product_id_color_location_generation_identity"
    ].update(
        {
            "decoder_step_export_generation": "step-export-70",
            "color_assembly_generation": "assembly-70",
            "hierarchy_assembly_generation": "assembly-69",
            "location_assembly_generation": "assembly-68",
            "decoded_product_ids": ["root", "housing", "rotor"],
            "decoded_parent_product_ids": ["", "root", "housing"],
            "decoded_component_location_sha256": ["5" * 64, "7" * 64, "6" * 64],
            "decoded_assembly_metadata_sha256": "f" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_assembly_uses_current_product_colors_hierarchy_and_locations"
    ]


def test_v21_source_sketch_constraint_entity_id_solver_order_generation_mismatch():
    row = _source_v21()
    row["replay_identity"][
        "sketch_constraint_entity_id_solver_order_generation_identity"
    ].update(
        {
            "entity_table_sketch_generation": "sketch-70",
            "solver_order_sketch_generation": "sketch-69",
            "replay_entity_ids": [101, 104, 103],
            "replay_solver_constraint_order": [202, 201],
            "replay_constraint_entity_ids": [[101, 104], [104, 103]],
            "replay_entity_table_sha256": "b" * 64,
            "replay_constraint_table_sha256": "c" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "sketch_constraints_use_current_entity_ids_and_solver_order"
    ]
