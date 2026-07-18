from __future__ import annotations

import math

from test_build123d_generalization_v28 import (
    _public_result,
    _public_v28,
    _source_result,
    _source_v28,
)


_PROMOTED_CASE_IDS = (
    "v29_public_sheet_metal_bend_allowance_kfactor_neutral_axis_flat_pattern_area_mismatch",
    "v29_public_joint_kinematic_loop_dof_limit_frame_closure_swept_volume_mismatch",
    "v29_source_dxf_arc_spline_layer_plane_unit_closed_wire_face_digest_mismatch",
    "v29_source_3mf_component_transform_triangle_winding_material_watertight_volume_digest_mismatch",
)


def _public_v29():
    reference, measured = _public_v28()
    neutral_radius = 0.003 + 0.42 * 0.0015
    bend_allowance = math.radians(90.0) * neutral_radius
    flat_area = 0.05 * (0.1 + 0.08 + bend_allowance)
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            sheet_digest = ("1" if index == 0 else "2") * 64
            joint_digest = ("3" if index == 0 else "4") * 64
            generation = "sheet-metal-161"
            row[
                "sheet_metal_bend_allowance_kfactor_neutral_axis_relief_thickness_flat_pattern_area_generation_identity"
            ] = {
                "sheet_generation": generation,
                "bend_sheet_generation": generation,
                "neutral_axis_sheet_generation": generation,
                "relief_sheet_generation": generation,
                "thickness_sheet_generation": generation,
                "pattern_sheet_generation": generation,
                "area_sheet_generation": generation,
                "result_sheet_generation": generation,
                "bend_radius_m": 0.003,
                "result_bend_radius_m": 0.003,
                "bend_angle_deg": 90.0,
                "result_bend_angle_deg": 90.0,
                "k_factor": 0.42,
                "result_k_factor": 0.42,
                "neutral_axis_radius_m": neutral_radius,
                "result_neutral_axis_radius_m": neutral_radius,
                "bend_allowance_m": bend_allowance,
                "result_bend_allowance_m": bend_allowance,
                "relief_type": "rectangular",
                "result_relief_type": "rectangular",
                "relief_width_m": 0.002,
                "result_relief_width_m": 0.002,
                "thickness_m": 0.0015,
                "result_thickness_m": 0.0015,
                "strip_width_m": 0.05,
                "result_strip_width_m": 0.05,
                "straight_lengths_m": [0.1, 0.08],
                "result_straight_lengths_m": [0.1, 0.08],
                "flat_pattern_area_m2": flat_area,
                "result_flat_pattern_area_m2": flat_area,
                "flat_pattern_wire_closed": True,
                "result_flat_pattern_wire_closed": True,
                "flat_pattern_sha256": sheet_digest,
                "result_flat_pattern_sha256": sheet_digest,
            }
            generation = "joint-loop-161"
            row[
                "joint_kinematic_loop_graph_dof_limit_connector_frame_closure_configuration_swept_volume_generation_identity"
            ] = {
                "joint_generation": generation,
                "graph_joint_generation": generation,
                "dof_joint_generation": generation,
                "limit_joint_generation": generation,
                "frame_joint_generation": generation,
                "closure_joint_generation": generation,
                "configuration_joint_generation": generation,
                "swept_joint_generation": generation,
                "result_joint_generation": generation,
                "joint_graph_edges": [
                    ["ground", "revolute-a", "link-a"],
                    ["link-a", "prismatic-b", "slider-b"],
                    ["slider-b", "fixed-c", "ground"],
                ],
                "result_joint_graph_edges": [
                    ["ground", "revolute-a", "link-a"],
                    ["link-a", "prismatic-b", "slider-b"],
                    ["slider-b", "fixed-c", "ground"],
                ],
                "dof_names": ["theta_a_deg", "travel_b_m"],
                "result_dof_names": ["theta_a_deg", "travel_b_m"],
                "dof_types": ["revolute", "prismatic"],
                "result_dof_types": ["revolute", "prismatic"],
                "lower_limits": [-30.0, 0.0],
                "result_lower_limits": [-30.0, 0.0],
                "upper_limits": [30.0, 0.02],
                "result_upper_limits": [30.0, 0.02],
                "configuration_values": [10.0, 0.012],
                "result_configuration_values": [10.0, 0.012],
                "connector_frame_sha256": "5" * 64,
                "result_connector_frame_sha256": "5" * 64,
                "loop_closure_error_m": 2.0e-10,
                "result_loop_closure_error_m": 2.0e-10,
                "loop_closure_tolerance_m": 1.0e-8,
                "result_loop_closure_tolerance_m": 1.0e-8,
                "configuration_id": "pose-10deg-12mm",
                "result_configuration_id": "pose-10deg-12mm",
                "swept_volume_m3": 0.00073,
                "result_swept_volume_m3": 0.00073,
                "swept_shape_sha256": joint_digest,
                "result_swept_shape_sha256": joint_digest,
            }
    return reference, measured


def _source_v29():
    row = _source_v28()
    identity = row["replay_identity"]
    generation = "dxf-face-161"
    identity[
        "dxf_arc_spline_layer_plane_unit_closed_wire_orientation_face_digest_generation_identity"
    ] = {
        "dxf_generation": generation,
        "arc_dxf_generation": generation,
        "spline_dxf_generation": generation,
        "layer_dxf_generation": generation,
        "plane_dxf_generation": generation,
        "unit_dxf_generation": generation,
        "wire_dxf_generation": generation,
        "face_dxf_generation": generation,
        "result_dxf_generation": generation,
        "arc_parameters": [["arc-1", 0.0, 0.0, 10.0, 0.0, 180.0]],
        "decoded_arc_parameters": [["arc-1", 0.0, 0.0, 10.0, 0.0, 180.0]],
        "spline_parameters": [["spline-1", 3, 4, "6" * 64]],
        "decoded_spline_parameters": [["spline-1", 3, 4, "6" * 64]],
        "entity_layer_map": [["arc-1", "profile"], ["spline-1", "profile"]],
        "decoded_entity_layer_map": [["arc-1", "profile"], ["spline-1", "profile"]],
        "workplane_matrix": [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        "decoded_workplane_matrix": [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        "length_unit": "mm",
        "decoded_length_unit": "mm",
        "wire_closed": True,
        "decoded_wire_closed": True,
        "wire_orientation": "counterclockwise",
        "decoded_wire_orientation": "counterclockwise",
        "face_area_mm2": 314.1592653589793,
        "decoded_face_area_mm2": 314.1592653589793,
        "dxf_sha256": "7" * 64,
        "decoded_dxf_sha256": "7" * 64,
        "face_sha256": "8" * 64,
        "decoded_face_sha256": "8" * 64,
    }
    generation = "3mf-import-161"
    identity[
        "three_mf_component_transform_triangle_winding_material_watertight_volume_unit_digest_generation_identity"
    ] = {
        "three_mf_generation": generation,
        "component_three_mf_generation": generation,
        "transform_three_mf_generation": generation,
        "winding_three_mf_generation": generation,
        "material_three_mf_generation": generation,
        "watertight_three_mf_generation": generation,
        "volume_three_mf_generation": generation,
        "file_three_mf_generation": generation,
        "result_three_mf_generation": generation,
        "component_names": ["housing", "insert"],
        "decoded_component_names": ["housing", "insert"],
        "component_transform_sha256": [["housing", "9" * 64], ["insert", "a" * 64]],
        "decoded_component_transform_sha256": [["housing", "9" * 64], ["insert", "a" * 64]],
        "triangle_winding": "outward_counterclockwise",
        "decoded_triangle_winding": "outward_counterclockwise",
        "material_id_map": [["housing", 1], ["insert", 2]],
        "decoded_material_id_map": [["housing", 1], ["insert", 2]],
        "watertight_component_names": ["housing", "insert"],
        "decoded_watertight_component_names": ["housing", "insert"],
        "component_volumes_mm3": [["housing", 1200.0], ["insert", 300.0]],
        "decoded_component_volumes_mm3": [["housing", 1200.0], ["insert", 300.0]],
        "total_volume_mm3": 1500.0,
        "decoded_total_volume_mm3": 1500.0,
        "length_unit": "mm",
        "decoded_length_unit": "mm",
        "three_mf_sha256": "b" * 64,
        "decoded_three_mf_sha256": "b" * 64,
        "triangle_mesh_sha256": "c" * 64,
        "decoded_triangle_mesh_sha256": "c" * 64,
    }
    return row


def test_v29_positive_public_and_source_identity():
    reference, measured = _public_v29()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v29())["status"] == "ok"


def test_v29_public_sheet_metal_bend_allowance_kfactor_neutral_axis_flat_pattern_area_mismatch():
    reference, measured = _public_v29()
    measured["external_cad"][0][
        "sheet_metal_bend_allowance_kfactor_neutral_axis_relief_thickness_flat_pattern_area_generation_identity"
    ].update(
        {
            "bend_sheet_generation": "sheet-metal-160",
            "pattern_sheet_generation": "sheet-metal-159",
            "result_sheet_generation": "sheet-metal-158",
            "result_bend_radius_m": 0.004,
            "result_bend_angle_deg": 80.0,
            "result_k_factor": 0.5,
            "result_neutral_axis_radius_m": 0.00475,
            "result_bend_allowance_m": 0.0066,
            "result_relief_type": "none",
            "result_relief_width_m": 0.0,
            "result_thickness_m": 0.002,
            "result_strip_width_m": 0.06,
            "result_straight_lengths_m": [0.1, 0.07],
            "result_flat_pattern_area_m2": 0.0101,
            "result_flat_pattern_wire_closed": False,
            "result_flat_pattern_sha256": "d" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "sheet_metal_flat_patterns_use_current_bends_neutral_axis_relief_thickness_and_area"
    ]


def test_v29_public_joint_kinematic_loop_dof_limit_frame_closure_swept_volume_mismatch():
    reference, measured = _public_v29()
    measured["external_cad"][0][
        "joint_kinematic_loop_graph_dof_limit_connector_frame_closure_configuration_swept_volume_generation_identity"
    ].update(
        {
            "graph_joint_generation": "joint-loop-160",
            "closure_joint_generation": "joint-loop-159",
            "result_joint_generation": "joint-loop-158",
            "result_joint_graph_edges": [["ground", "revolute-a", "slider-b"]],
            "result_dof_names": ["travel_b_m", "theta_a_deg"],
            "result_dof_types": ["prismatic", "fixed"],
            "result_lower_limits": [0.02, 30.0],
            "result_upper_limits": [0.0, -30.0],
            "result_configuration_values": [0.03, 45.0],
            "result_connector_frame_sha256": "e" * 64,
            "result_loop_closure_error_m": 1.0e-3,
            "result_loop_closure_tolerance_m": 1.0e-6,
            "result_configuration_id": "stale-pose",
            "result_swept_volume_m3": 0.00051,
            "result_swept_shape_sha256": "f" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "joint_loops_use_current_graph_dofs_limits_frames_closure_configuration_and_swept_volume"
    ]


def test_v29_source_dxf_arc_spline_layer_plane_unit_closed_wire_face_digest_mismatch():
    row = _source_v29()
    row["replay_identity"][
        "dxf_arc_spline_layer_plane_unit_closed_wire_orientation_face_digest_generation_identity"
    ].update(
        {
            "arc_dxf_generation": "dxf-face-160",
            "face_dxf_generation": "dxf-face-159",
            "result_dxf_generation": "dxf-face-158",
            "decoded_arc_parameters": [["arc-1", 0.0, 0.0, 9.0, 0.0, 170.0]],
            "decoded_spline_parameters": [["spline-1", 2, 3, "0" * 64]],
            "decoded_entity_layer_map": [["arc-1", "construction"]],
            "decoded_workplane_matrix": [
                [1.0, 0.0, 0.0, 10.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
            "decoded_length_unit": "in",
            "decoded_wire_closed": False,
            "decoded_wire_orientation": "clockwise",
            "decoded_face_area_mm2": 280.0,
            "decoded_dxf_sha256": "1" * 64,
            "decoded_face_sha256": "2" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "dxf_faces_use_current_arcs_splines_layers_plane_units_closure_orientation_and_digests"
    ]


def test_v29_source_3mf_component_transform_triangle_winding_material_watertight_volume_digest_mismatch():
    row = _source_v29()
    row["replay_identity"][
        "three_mf_component_transform_triangle_winding_material_watertight_volume_unit_digest_generation_identity"
    ].update(
        {
            "component_three_mf_generation": "3mf-import-160",
            "winding_three_mf_generation": "3mf-import-159",
            "result_three_mf_generation": "3mf-import-158",
            "decoded_component_names": ["insert", "housing"],
            "decoded_component_transform_sha256": [["housing", "3" * 64]],
            "decoded_triangle_winding": "inward_clockwise",
            "decoded_material_id_map": [["housing", 2], ["insert", 1]],
            "decoded_watertight_component_names": ["housing"],
            "decoded_component_volumes_mm3": [["housing", 1000.0], ["insert", -20.0]],
            "decoded_total_volume_mm3": 980.0,
            "decoded_length_unit": "m",
            "decoded_three_mf_sha256": "4" * 64,
            "decoded_triangle_mesh_sha256": "5" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "three_mf_imports_use_current_components_transforms_winding_materials_watertight_volumes_units_and_digests"
    ]
