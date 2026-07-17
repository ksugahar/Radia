from __future__ import annotations

from test_build123d_generalization_v27 import (
    _public_result,
    _public_v27,
    _source_result,
    _source_v27,
)


def _public_v28():
    reference, measured = _public_v27()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["shell_offset_face_normal_thickness_join_self_intersection_mass_generation_identity"] = {
                "shell_generation": "shell-offset-151", "face_shell_generation": "shell-offset-151",
                "normal_shell_generation": "shell-offset-151", "thickness_shell_generation": "shell-offset-151",
                "join_shell_generation": "shell-offset-151", "intersection_shell_generation": "shell-offset-151",
                "shape_shell_generation": "shell-offset-151", "mass_shell_generation": "shell-offset-151",
                "result_shell_generation": "shell-offset-151",
                "selected_face_names": ["top-face", "bottom-face"],
                "result_selected_face_names": ["top-face", "bottom-face"],
                "normal_direction": "outward", "result_normal_direction": "outward",
                "thickness_m": 0.002, "result_thickness_m": 0.002,
                "join_mode": "arc", "result_join_mode": "arc",
                "self_intersection": False, "result_self_intersection": False,
                "shape_generation": "shell-shape-151", "result_shape_generation": "shell-shape-151",
                "is_valid_solid": True, "result_is_valid_solid": True,
                "volume_m3": 0.00042, "result_volume_m3": 0.00042,
                "mass_property_sha256": digest, "result_mass_property_sha256": digest,
            }
            row["path_sweep_frame_transition_profile_orientation_solid_volume_generation_identity"] = {
                "sweep_generation": "path-sweep-151", "path_sweep_generation": "path-sweep-151",
                "frame_sweep_generation": "path-sweep-151", "transition_sweep_generation": "path-sweep-151",
                "profile_sweep_generation": "path-sweep-151", "orientation_sweep_generation": "path-sweep-151",
                "solid_sweep_generation": "path-sweep-151", "volume_sweep_generation": "path-sweep-151",
                "result_sweep_generation": "path-sweep-151",
                "path_edge_names": ["path-0", "path-1", "path-2"],
                "result_path_edge_names": ["path-0", "path-1", "path-2"],
                "moving_frame": "parallel_transport", "result_moving_frame": "parallel_transport",
                "transition_mode": "round", "result_transition_mode": "round",
                "profile_wire_names": ["profile-outer"], "result_profile_wire_names": ["profile-outer"],
                "profile_orientation_deg": [0.0, 12.0, 24.0],
                "result_profile_orientation_deg": [0.0, 12.0, 24.0],
                "is_valid_solid": True, "result_is_valid_solid": True,
                "volume_m3": 0.00031, "result_volume_m3": 0.00031,
                "sweep_shape_sha256": digest, "result_sweep_shape_sha256": digest,
            }
    return reference, measured


def _source_v28():
    row = _source_v27()
    identity = row["replay_identity"]
    identity["step_ap242_context_product_uuid_unit_color_placement_shape_file_generation_identity"] = {
        "step_generation": "ap242-import-151", "context_step_generation": "ap242-import-151",
        "product_step_generation": "ap242-import-151", "unit_step_generation": "ap242-import-151",
        "color_step_generation": "ap242-import-151", "placement_step_generation": "ap242-import-151",
        "shape_step_generation": "ap242-import-151", "file_step_generation": "ap242-import-151",
        "result_step_generation": "ap242-import-151", "representation_context": "mechanical_design_3d",
        "decoded_representation_context": "mechanical_design_3d",
        "product_uuid_map": [["base", "11111111-1111-4111-8111-111111111111"], ["cover", "22222222-2222-4222-8222-222222222222"]],
        "decoded_product_uuid_map": [["base", "11111111-1111-4111-8111-111111111111"], ["cover", "22222222-2222-4222-8222-222222222222"]],
        "length_unit": "mm", "decoded_length_unit": "mm",
        "colors_rgb": [["base", 0.8, 0.1, 0.1], ["cover", 0.1, 0.1, 0.8]],
        "decoded_colors_rgb": [["base", 0.8, 0.1, 0.1], ["cover", 0.1, 0.1, 0.8]],
        "placements": [["base", 0.0, 0.0, 0.0], ["cover", 0.0, 0.0, 10.0]],
        "decoded_placements": [["base", 0.0, 0.0, 0.0], ["cover", 0.0, 0.0, 10.0]],
        "shape_owner_map": [["base-shape", "base"], ["cover-shape", "cover"]],
        "decoded_shape_owner_map": [["base-shape", "base"], ["cover-shape", "cover"]],
        "step_sha256": "3" * 64, "decoded_step_sha256": "3" * 64,
        "shape_map_sha256": "4" * 64, "decoded_shape_map_sha256": "4" * 64,
    }
    identity["occt_kernel_shape_location_tolerance_triangulation_serialization_cache_generation_identity"] = {
        "cache_generation": "occt-cache-151", "kernel_cache_generation": "occt-cache-151",
        "shape_cache_generation": "occt-cache-151", "location_cache_generation": "occt-cache-151",
        "tolerance_cache_generation": "occt-cache-151", "triangulation_cache_generation": "occt-cache-151",
        "serialization_cache_generation": "occt-cache-151", "result_cache_generation": "occt-cache-151",
        "kernel_version": "OCCT-7.8.1", "decoded_kernel_version": "OCCT-7.8.1",
        "shape_sha256": "5" * 64, "decoded_shape_sha256": "5" * 64,
        "location_transform": [[1.0, 0.0, 0.0, 0.01], [0.0, 1.0, 0.0, 0.02], [0.0, 0.0, 1.0, 0.03]],
        "decoded_location_transform": [[1.0, 0.0, 0.0, 0.01], [0.0, 1.0, 0.0, 0.02], [0.0, 0.0, 1.0, 0.03]],
        "modeling_tolerance_m": 1.0e-7, "decoded_modeling_tolerance_m": 1.0e-7,
        "triangulation_parameters": [0.001, 0.2, 1.0],
        "decoded_triangulation_parameters": [0.001, 0.2, 1.0],
        "triangulation_sha256": "6" * 64, "decoded_triangulation_sha256": "6" * 64,
        "serialization_format": "brep-binary-v1", "decoded_serialization_format": "brep-binary-v1",
        "serialization_sha256": "7" * 64, "decoded_serialization_sha256": "7" * 64,
        "cache_sha256": "8" * 64, "decoded_cache_sha256": "8" * 64,
    }
    return row


def test_v28_positive_public_and_source_identity():
    reference, measured = _public_v28()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v28())["status"] == "ok"


def test_v28_public_shell_offset_face_normal_thickness_join_self_intersection_mass_property_mismatch():
    reference, measured = _public_v28()
    measured["external_cad"][0]["shell_offset_face_normal_thickness_join_self_intersection_mass_generation_identity"].update(
        {"face_shell_generation": "shell-offset-150", "result_selected_face_names": ["side-face"],
         "result_normal_direction": "inward", "result_thickness_m": 0.003,
         "result_join_mode": "intersection", "result_self_intersection": True,
         "result_shape_generation": "shell-shape-150", "result_is_valid_solid": False,
         "result_volume_m3": 0.00037, "result_mass_property_sha256": "9" * 64}
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["shell_offsets_use_current_faces_normals_thickness_join_intersection_shape_and_mass"]


def test_v28_public_path_sweep_frame_transition_profile_orientation_solid_validity_volume_mismatch():
    reference, measured = _public_v28()
    measured["external_cad"][0]["path_sweep_frame_transition_profile_orientation_solid_volume_generation_identity"].update(
        {"path_sweep_generation": "path-sweep-150", "result_path_edge_names": ["path-0", "path-2", "path-1"],
         "result_moving_frame": "frenet", "result_transition_mode": "right",
         "result_profile_wire_names": ["profile-inner"], "result_profile_orientation_deg": [0.0, -12.0, -24.0],
         "result_is_valid_solid": False, "result_volume_m3": 0.00025,
         "result_sweep_shape_sha256": "a" * 64}
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["path_sweeps_use_current_path_frame_transition_profile_orientation_solid_and_volume"]


def test_v28_source_step_ap242_representation_context_product_uuid_unit_color_placement_mismatch():
    row = _source_v28()
    row["replay_identity"]["step_ap242_context_product_uuid_unit_color_placement_shape_file_generation_identity"].update(
        {"context_step_generation": "ap242-import-150", "decoded_representation_context": "geometric_curve_set",
         "decoded_product_uuid_map": [["base", "33333333-3333-4333-8333-333333333333"]],
         "decoded_length_unit": "m", "decoded_shape_owner_map": [["base-shape", "cover"]],
         "decoded_step_sha256": "b" * 64}
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["step_ap242_import_uses_current_context_products_units_colors_placements_shapes_and_file"]


def test_v28_source_occt_kernel_shape_hash_location_tolerance_triangulation_serialization_cache_mismatch():
    row = _source_v28()
    row["replay_identity"]["occt_kernel_shape_location_tolerance_triangulation_serialization_cache_generation_identity"].update(
        {"kernel_cache_generation": "occt-cache-150", "decoded_kernel_version": "OCCT-7.7.0",
         "decoded_shape_sha256": "d" * 64, "decoded_modeling_tolerance_m": 1.0e-4,
         "decoded_triangulation_parameters": [0.01, 0.5, 0.0],
         "decoded_serialization_format": "step", "decoded_cache_sha256": "0" * 64}
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["occt_shape_cache_uses_current_kernel_shape_location_tolerance_triangulation_and_serialization"]
