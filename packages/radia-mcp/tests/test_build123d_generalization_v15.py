from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v14 import _public_v14, _source_v14


def _public_v15():
    reference, measured = _public_v14()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("b" if index == 0 else "c") * 64
            row["boolean_tolerance_length_unit_identity"] = {
                "boolean_generation": "boolean-47",
                "result_geometry_generation": "boolean-47",
                "model_length_unit": "mm",
                "tolerance_value": 1.0e-6,
                "tolerance_unit": "mm",
                "tolerance_scale_to_m": 1.0e-3,
                "kernel_tolerance_m": 1.0e-9,
                "input_shape_sha256": digest,
                "boolean_input_shape_sha256": digest,
            }
            row["nested_assembly_placement_order_identity"] = {
                "assembly_generation": "assembly-47",
                "parent_placement_generation": "assembly-47",
                "child_placement_generation": "assembly-47",
                "world_placement_generation": "assembly-47",
                "multiplication_order": "parent_then_child",
                "applied_multiplication_order": "parent_then_child",
                "placement_chain_sha256": digest,
                "world_transform_sha256": digest,
            }
    return reference, measured


def _source_v15():
    row = _source_v14()
    identity = row["replay_identity"]
    identity["step_color_label_topology_identity"] = {
        "step_import_generation": "step-import-47",
        "topology_generation": "step-topology-47",
        "attribute_map_topology_generation": "step-topology-47",
        "face_ids": [31, 32, 33],
        "attribute_face_ids": [31, 32, 33],
        "labels": ["housing", "shaft", "terminal"],
        "colors_rgb": [[0.8, 0.8, 0.8], [0.2, 0.2, 0.2], [0.8, 0.1, 0.1]],
        "attribute_map_sha256": "d" * 64,
        "resolved_attribute_map_sha256": "d" * 64,
    }
    identity["brep_surface_parameter_orientation_identity"] = {
        "serialization_generation": "brep-serialization-47",
        "surface_parameter_generation": "brep-serialization-47",
        "surface_ids": [41, 42],
        "exported_surface_ids": [41, 42],
        "parameter_orientation": "u_cross_v_outward",
        "exported_parameter_orientation": "u_cross_v_outward",
        "parameter_range_sha256": "e" * 64,
        "exported_parameter_range_sha256": "e" * 64,
    }
    return row


def test_v15_positive_public_and_source_lineage():
    reference, measured = _public_v15()
    assert _public_result(reference, measured)["status"] == "ok"
    source = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v15()))
    )
    assert source["status"] == "ok"


def test_v15_public_boolean_tolerance_model_length_unit_mismatch():
    reference, measured = _public_v15()
    measured["external_cad"][0]["boolean_tolerance_length_unit_identity"].update(
        {
            "tolerance_unit": "m",
            "tolerance_scale_to_m": 1.0,
            "kernel_tolerance_m": 1.0e-6,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["boolean_tolerance_uses_one_physical_model_length_basis"] is False


def test_v15_public_nested_assembly_placement_multiplication_order_mismatch():
    reference, measured = _public_v15()
    measured["external_cad"][0]["nested_assembly_placement_order_identity"].update(
        {
            "applied_multiplication_order": "child_then_parent",
            "world_placement_generation": "assembly-46",
            "world_transform_sha256": "f" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert result["checks"]["nested_assembly_placements_use_parent_then_child_order"] is False


def test_v15_source_step_color_label_map_previous_topology_generation():
    row = _source_v15()
    row["replay_identity"]["step_color_label_topology_identity"].update(
        {
            "attribute_map_topology_generation": "step-topology-46",
            "attribute_face_ids": [30, 32, 33],
            "resolved_attribute_map_sha256": "f" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["step_color_labels_follow_current_topology_generation"] is False


def test_v15_source_brep_surface_parameter_orientation_range_mismatch():
    row = _source_v15()
    row["replay_identity"]["brep_surface_parameter_orientation_identity"].update(
        {
            "surface_parameter_generation": "brep-serialization-46",
            "exported_parameter_orientation": "v_cross_u_outward",
            "exported_parameter_range_sha256": "f" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert (
        result["checks"][
            "brep_surface_parameter_ranges_preserve_outward_orientation"
        ]
        is False
    )
