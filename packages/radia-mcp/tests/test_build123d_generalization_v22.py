from __future__ import annotations

import json

from radia_mcp.build123d.server import build123d_jointed_assembly_source_replay_gate
from test_build123d_generalization_v11 import _public_result
from test_build123d_generalization_v21 import _public_v21, _source_v21


def _public_v22():
    reference, measured = _public_v21()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["joint_connector_frame_labeled_face_subshape_generation_identity"] = {
                "shape_generation": "shape-81",
                "label_table_shape_generation": "shape-81",
                "connector_shape_generation": "shape-81",
                "location_shape_generation": "shape-81",
                "labeled_face_subshape_id": "face:mount_a",
                "resolved_labeled_face_subshape_id": "face:mount_a",
                "labeled_face_geometry_sha256": digest,
                "resolved_labeled_face_geometry_sha256": digest,
                "connector_origin": [0.0, 0.0, 0.0],
                "evaluated_connector_origin": [0.0, 0.0, 0.0],
                "connector_axis": [0.0, 0.0, 1.0],
                "evaluated_connector_axis": [0.0, 0.0, 1.0],
                "parent_location_sha256": "3" * 64,
                "evaluated_parent_location_sha256": "3" * 64,
                "connector_frame_sha256": "4" * 64,
                "evaluated_connector_frame_sha256": "4" * 64,
            }
            row[
                "inertia_tensor_principal_axes_density_unit_location_generation_identity"
            ] = {
                "shape_generation": "shape-81",
                "density_shape_generation": "shape-81",
                "mass_property_shape_generation": "shape-81",
                "location_shape_generation": "shape-81",
                "principal_axis_shape_generation": "shape-81",
                "density_value": 7800.0,
                "evaluated_density_value": 7800.0,
                "density_unit": "kg/m^3",
                "evaluated_density_unit": "kg/m^3",
                "shape_location_sha256": "5" * 64,
                "evaluated_shape_location_sha256": "5" * 64,
                "center_of_mass": [0.01, 0.02, 0.03],
                "evaluated_center_of_mass": [0.01, 0.02, 0.03],
                "inertia_tensor_sha256": "6" * 64,
                "evaluated_inertia_tensor_sha256": "6" * 64,
                "principal_axes_sha256": "7" * 64,
                "evaluated_principal_axes_sha256": "7" * 64,
                "mass_property_sha256": digest,
                "evaluated_mass_property_sha256": digest,
            }
    return reference, measured


def _source_v22():
    row = _source_v21()
    identity = row["replay_identity"]
    identity["step_ap242_component_transform_name_product_generation_identity"] = {
        "step_export_generation": "step-export-81",
        "decoder_step_export_generation": "step-export-81",
        "assembly_generation": "assembly-81",
        "product_assembly_generation": "assembly-81",
        "name_assembly_generation": "assembly-81",
        "transform_assembly_generation": "assembly-81",
        "component_product_ids": ["root", "rotor", "housing"],
        "decoded_component_product_ids": ["root", "rotor", "housing"],
        "component_names": ["machine", "rotor", "housing"],
        "decoded_component_names": ["machine", "rotor", "housing"],
        "nested_transform_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "decoded_nested_transform_sha256": ["1" * 64, "2" * 64, "3" * 64],
        "ap242_product_map_sha256": "4" * 64,
        "decoded_ap242_product_map_sha256": "4" * 64,
    }
    identity["curved_mesh_export_edge_chord_surface_label_generation_identity"] = {
        "shape_generation": "shape-81",
        "mesh_shape_generation": "shape-81",
        "edge_curve_shape_generation": "shape-81",
        "surface_label_shape_generation": "shape-81",
        "mesh_export_generation": "mesh-export-81",
        "metric_mesh_export_generation": "mesh-export-81",
        "edge_curve_sha256": ["5" * 64, "6" * 64],
        "exported_edge_curve_sha256": ["5" * 64, "6" * 64],
        "chordal_tolerance": 1.0e-4,
        "evaluated_chordal_tolerance": 1.0e-4,
        "length_unit": "m",
        "evaluated_length_unit": "m",
        "boundary_surface_labels": ["outer", "interface"],
        "exported_boundary_surface_labels": ["outer", "interface"],
        "curved_mesh_sha256": "7" * 64,
        "exported_curved_mesh_sha256": "7" * 64,
    }
    return row


def test_v22_positive_public_and_source_identity():
    reference, measured = _public_v22()
    assert _public_result(reference, measured)["status"] == "ok"
    source = json.loads(
        build123d_jointed_assembly_source_replay_gate(json.dumps(_source_v22()))
    )
    assert source["status"] == "ok"


def test_v22_public_joint_connector_frame_labeled_face_subshape_generation_mismatch():
    reference, measured = _public_v22()
    measured["external_cad"][0][
        "joint_connector_frame_labeled_face_subshape_generation_identity"
    ].update(
        {
            "label_table_shape_generation": "shape-80",
            "connector_shape_generation": "shape-79",
            "location_shape_generation": "shape-78",
            "resolved_labeled_face_subshape_id": "face:stale_mount",
            "resolved_labeled_face_geometry_sha256": "8" * 64,
            "evaluated_connector_axis": [0.0, 1.0, 0.0],
            "evaluated_parent_location_sha256": "9" * 64,
            "evaluated_connector_frame_sha256": "a" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "joint_connectors_use_current_labeled_face_frame_and_parent_location"
    ]


def test_v22_public_inertia_tensor_principal_axes_density_unit_location_generation_mismatch():
    reference, measured = _public_v22()
    measured["external_cad"][0][
        "inertia_tensor_principal_axes_density_unit_location_generation_identity"
    ].update(
        {
            "density_shape_generation": "shape-80",
            "mass_property_shape_generation": "shape-79",
            "location_shape_generation": "shape-78",
            "principal_axis_shape_generation": "shape-77",
            "evaluated_density_value": 7.8,
            "evaluated_density_unit": "g/cm^3",
            "evaluated_shape_location_sha256": "b" * 64,
            "evaluated_center_of_mass": [10.0, 20.0, 30.0],
            "evaluated_inertia_tensor_sha256": "c" * 64,
            "evaluated_principal_axes_sha256": "d" * 64,
            "evaluated_mass_property_sha256": "e" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "inertia_uses_current_density_unit_location_and_principal_axes"
    ]


def test_v22_source_step_ap242_component_transform_name_product_generation_mismatch():
    row = _source_v22()
    row["replay_identity"][
        "step_ap242_component_transform_name_product_generation_identity"
    ].update(
        {
            "decoder_step_export_generation": "step-export-80",
            "product_assembly_generation": "assembly-80",
            "name_assembly_generation": "assembly-79",
            "transform_assembly_generation": "assembly-78",
            "decoded_component_product_ids": ["root", "housing", "rotor"],
            "decoded_component_names": ["machine", "old_housing", "old_rotor"],
            "decoded_nested_transform_sha256": ["1" * 64, "3" * 64, "2" * 64],
            "decoded_ap242_product_map_sha256": "f" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_ap242_components_use_current_products_names_and_nested_transforms"
    ]


def test_v22_source_curved_mesh_export_edge_chord_surface_label_generation_mismatch():
    row = _source_v22()
    row["replay_identity"][
        "curved_mesh_export_edge_chord_surface_label_generation_identity"
    ].update(
        {
            "mesh_shape_generation": "shape-80",
            "edge_curve_shape_generation": "shape-79",
            "surface_label_shape_generation": "shape-78",
            "metric_mesh_export_generation": "mesh-export-80",
            "exported_edge_curve_sha256": ["8" * 64, "9" * 64],
            "evaluated_chordal_tolerance": 0.1,
            "evaluated_length_unit": "mm",
            "exported_boundary_surface_labels": ["stale_outer", "stale_interface"],
            "exported_curved_mesh_sha256": "a" * 64,
        }
    )
    result = json.loads(build123d_jointed_assembly_source_replay_gate(json.dumps(row)))
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "curved_mesh_export_uses_current_edges_chord_and_surface_labels"
    ]
