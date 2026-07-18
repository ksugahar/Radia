from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v34 import _public_v34, _source_v34


_PROMOTED_CASE_IDS = (
    "v35_public_mirrored_pattern_instance_handedness_transform_mass_inertia_mismatch",
    "v35_public_offset_thicken_curvature_selfintersection_wall_volume_topology_mismatch",
    "v35_source_step_ap242_pmi_name_color_unit_occurrence_transform_roundtrip_mismatch",
    "v35_source_dxf_profile_arc_bulge_winding_plane_unit_extrusion_roundtrip_mismatch",
)


def _public_v35():
    reference, measured = _public_v34()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            suffix = "1" if index == 0 else "2"
            generation = "mirrored-pattern-221"
            row["mirrored_pattern_occurrence_handedness_transform_suppression_volume_mass_center_inertia_owner_shape_result_generation_identity"] = {
                "pattern_generation": generation,
                **{key: generation for key in ("occurrence_generation", "handedness_generation", "transform_generation", "suppression_generation", "volume_generation", "mass_generation", "inertia_generation", "owner_generation", "shape_generation", "result_generation")},
                "occurrence_ids": ["part:0", "part:mirror_x"], "result_occurrence_ids": ["part:0", "part:mirror_x"],
                "transform_determinants": [1.0, -1.0], "result_transform_determinants": [1.0, -1.0],
                "handedness": ["right", "left"], "result_handedness": ["right", "left"],
                "suppressed": [False, False], "result_suppressed": [False, False],
                "occurrence_volumes_m3": [1.0, 1.0], "result_occurrence_volumes_m3": [1.0, 1.0],
                "occurrence_masses_kg": [2.0, 2.0], "result_occurrence_masses_kg": [2.0, 2.0],
                "occurrence_centers_m": [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], "result_occurrence_centers_m": [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
                "occurrence_inertia_principal_kg_m2": [[0.1, 0.2, 0.3], [0.1, 0.2, 0.3]], "result_occurrence_inertia_principal_kg_m2": [[0.1, 0.2, 0.3], [0.1, 0.2, 0.3]],
                "assembly_volume_m3": 2.0, "result_assembly_volume_m3": 2.0,
                "assembly_mass_kg": 4.0, "result_assembly_mass_kg": 4.0,
                "assembly_center_m": [0.0, 0.0, 0.0], "result_assembly_center_m": [0.0, 0.0, 0.0],
                "assembly_owner": "pattern/root", "result_assembly_owner": "pattern/root",
                "mirrored_shape_sha256": suffix * 64, "accepted_mirrored_shape_sha256": suffix * 64,
            }
            generation = "offset-thicken-221"
            digest = ("3" if index == 0 else "4") * 64
            row["offset_thicken_curvature_sign_selfintersection_repair_thickness_volume_topology_convergence_brep_result_generation_identity"] = {
                "thicken_generation": generation,
                **{key: generation for key in ("curvature_generation", "offset_generation", "intersection_generation", "repair_generation", "thickness_generation", "volume_generation", "topology_generation", "convergence_generation", "brep_generation", "result_generation")},
                "minimum_curvature_radius_m": 0.05, "result_minimum_curvature_radius_m": 0.05,
                "offset_m": -0.01, "result_offset_m": -0.01,
                "self_intersection_count": 0, "result_self_intersection_count": 0,
                "repair_mode": "none_required", "result_repair_mode": "none_required",
                "wall_thickness_samples_m": [0.01, 0.010001, 0.009999], "result_wall_thickness_samples_m": [0.01, 0.010001, 0.009999],
                "original_volume_m3": 1.0, "result_original_volume_m3": 1.0,
                "removed_volume_m3": 0.7, "result_removed_volume_m3": 0.7,
                "thickened_volume_m3": 0.3, "result_thickened_volume_m3": 0.3,
                "solid_count": 1, "result_solid_count": 1, "shell_count": 1, "result_shell_count": 1,
                "convergence_tolerances_m": [1.0e-4, 1.0e-5, 1.0e-6], "result_convergence_tolerances_m": [1.0e-4, 1.0e-5, 1.0e-6],
                "convergence_volumes_m3": [0.3001, 0.300001, 0.3], "result_convergence_volumes_m3": [0.3001, 0.300001, 0.3],
                "thicken_brep_sha256": digest, "accepted_thicken_brep_sha256": digest,
            }
    return reference, measured


def _source_v35():
    row = _source_v34()
    identity = row["replay_identity"]
    generation = "step-ap242-pmi-221"
    identity["step_ap242_pmi_unit_name_color_occurrence_transform_validity_owner_file_result_generation_identity"] = {
        "ap242_generation": generation,
        **{key: generation for key in ("pmi_generation", "unit_generation", "name_generation", "color_generation", "occurrence_generation", "transform_generation", "validity_generation", "owner_generation", "file_generation", "result_generation")},
        "schema": "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING", "decoded_schema": "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING",
        "length_unit": "m", "decoded_length_unit": "m",
        "pmi_annotations": [[1, "linear_dimension", "mm", 25.0], [2, "diameter", "mm", 10.0]], "decoded_pmi_annotations": [[1, "linear_dimension", "mm", 25.0], [2, "diameter", "mm", 10.0]],
        "product_names": [[1, "base"], [2, "arm"]], "decoded_product_names": [[1, "base"], [2, "arm"]],
        "colors_rgb": [[1, [0.2, 0.3, 0.4]], [2, [0.8, 0.1, 0.1]]], "decoded_colors_rgb": [[1, [0.2, 0.3, 0.4]], [2, [0.8, 0.1, 0.1]]],
        "occurrence_paths": [[1, "root/base"], [2, "root/arm"]], "decoded_occurrence_paths": [[1, "root/base"], [2, "root/arm"]],
        "transform_sha256": [[1, "5" * 64], [2, "6" * 64]], "decoded_transform_sha256": [[1, "5" * 64], [2, "6" * 64]],
        "solid_validity": [[1, True], [2, True]], "decoded_solid_validity": [[1, True], [2, True]],
        "source_owner": "assembly/root", "decoded_source_owner": "assembly/root",
        "ap242_file_sha256": "7" * 64, "decoded_ap242_file_sha256": "7" * 64,
    }
    generation = "dxf-profile-221"
    identity["dxf_profile_unit_plane_layer_arc_bulge_loop_winding_extrusion_topology_owner_file_result_generation_identity"] = {
        "dxf_generation": generation,
        **{key: generation for key in ("unit_generation", "plane_generation", "layer_generation", "arc_generation", "loop_generation", "winding_generation", "extrusion_generation", "topology_generation", "owner_generation", "file_generation", "result_generation")},
        "length_unit": "mm", "decoded_length_unit": "mm", "sketch_plane": "XY", "decoded_sketch_plane": "XY",
        "layer_names": ["OUTER", "HOLES"], "decoded_layer_names": ["OUTER", "HOLES"],
        "arc_bulges": [[11, 0.41421356237], [12, -0.41421356237]], "decoded_arc_bulges": [[11, 0.41421356237], [12, -0.41421356237]],
        "loop_winding_signs": [["outer", 1], ["hole", -1]], "decoded_loop_winding_signs": [["outer", 1], ["hole", -1]],
        "closed_loop_count": 2, "decoded_closed_loop_count": 2,
        "extrusion_height_m": 0.01, "decoded_extrusion_height_m": 0.01,
        "profile_area_m2": 0.1, "decoded_profile_area_m2": 0.1,
        "solid_count": 1, "decoded_solid_count": 1,
        "extruded_volume_m3": 0.001, "decoded_extruded_volume_m3": 0.001,
        "profile_owner": "dxf/profile1", "decoded_profile_owner": "dxf/profile1",
        "dxf_file_sha256": "8" * 64, "decoded_dxf_file_sha256": "8" * 64,
    }
    return row


def test_v35_positive_all_four_contracts():
    reference, measured = _public_v35()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v35())["status"] == "ok"


def test_v35_public_mirrored_pattern_instance_handedness_transform_mass_inertia_mismatch():
    reference, measured = _public_v35()
    identity = measured["external_cad"][0]["mirrored_pattern_occurrence_handedness_transform_suppression_volume_mass_center_inertia_owner_shape_result_generation_identity"]
    identity.update({"handedness_generation": "mirrored-pattern-220", "mass_generation": "mirrored-pattern-219", "result_generation": "mirrored-pattern-218", "result_occurrence_ids": ["part:0"], "result_transform_determinants": [1.0, 1.0], "result_handedness": ["right", "right"], "result_suppressed": [False, True], "result_occurrence_volumes_m3": [1.0, 0.5], "result_occurrence_masses_kg": [2.0, 1.0], "result_occurrence_centers_m": [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], "result_occurrence_inertia_principal_kg_m2": [[0.1, 0.2, -0.3]], "result_assembly_volume_m3": 1.5, "result_assembly_mass_kg": 3.0, "result_assembly_center_m": [1.5, 0.0, 0.0], "result_assembly_owner": "stale/root", "accepted_mirrored_shape_sha256": "9" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["mirrored_patterns_use_current_occurrences_handedness_transforms_suppression_mass_center_inertia_owner_and_shape"]


def test_v35_public_offset_thicken_curvature_selfintersection_wall_volume_topology_mismatch():
    reference, measured = _public_v35()
    identity = measured["external_cad"][0]["offset_thicken_curvature_sign_selfintersection_repair_thickness_volume_topology_convergence_brep_result_generation_identity"]
    identity.update({"curvature_generation": "offset-thicken-220", "topology_generation": "offset-thicken-219", "result_generation": "offset-thicken-218", "result_minimum_curvature_radius_m": 0.005, "result_offset_m": 0.02, "result_self_intersection_count": 3, "result_repair_mode": "discard_faces", "result_wall_thickness_samples_m": [0.001, 0.02, 0.04], "result_removed_volume_m3": 0.6, "result_thickened_volume_m3": 0.5, "result_solid_count": 2, "result_shell_count": 3, "result_convergence_tolerances_m": [1.0e-6, 1.0e-5, 1.0e-4], "result_convergence_volumes_m3": [0.3, 0.4, 0.6], "accepted_thicken_brep_sha256": "a" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["offset_thickens_use_current_curvature_sign_intersections_repair_thickness_volume_topology_convergence_and_brep"]


def test_v35_source_step_ap242_pmi_name_color_unit_occurrence_transform_roundtrip_mismatch():
    row = _source_v35(); identity = row["replay_identity"]["step_ap242_pmi_unit_name_color_occurrence_transform_validity_owner_file_result_generation_identity"]
    identity.update({"pmi_generation": "step-ap242-pmi-220", "occurrence_generation": "step-ap242-pmi-219", "result_generation": "step-ap242-pmi-218", "decoded_schema": "AP203", "decoded_length_unit": "mm", "decoded_pmi_annotations": [[1, "linear_dimension", "in", 25.0]], "decoded_product_names": [[1, "old"]], "decoded_colors_rgb": [[1, [0.0, 0.0, 0.0]]], "decoded_occurrence_paths": [[1, "old/root"]], "decoded_transform_sha256": [[1, "b" * 64]], "decoded_solid_validity": [[1, False]], "decoded_source_owner": "stale/root", "decoded_ap242_file_sha256": "c" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"
    assert not result["checks"]["step_ap242_roundtrips_use_current_pmi_units_names_colors_occurrences_transforms_validity_owner_and_file"]


def test_v35_source_dxf_profile_arc_bulge_winding_plane_unit_extrusion_roundtrip_mismatch():
    row = _source_v35(); identity = row["replay_identity"]["dxf_profile_unit_plane_layer_arc_bulge_loop_winding_extrusion_topology_owner_file_result_generation_identity"]
    identity.update({"arc_generation": "dxf-profile-220", "extrusion_generation": "dxf-profile-219", "result_generation": "dxf-profile-218", "decoded_length_unit": "in", "decoded_sketch_plane": "XZ", "decoded_layer_names": ["0"], "decoded_arc_bulges": [[11, -0.41421356237]], "decoded_loop_winding_signs": [["outer", -1], ["hole", 1]], "decoded_closed_loop_count": 1, "decoded_extrusion_height_m": 0.1, "decoded_profile_area_m2": 0.01, "decoded_solid_count": 2, "decoded_extruded_volume_m3": 0.1, "decoded_profile_owner": "stale/profile", "decoded_dxf_file_sha256": "d" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"
    assert not result["checks"]["dxf_profiles_use_current_units_plane_layers_arcs_bulges_winding_extrusion_topology_owner_and_file"]


def test_v35_rejects_self_consistent_mirror_handedness_error():
    reference, measured = _public_v35(); identity = measured["external_cad"][0]["mirrored_pattern_occurrence_handedness_transform_suppression_volume_mass_center_inertia_owner_shape_result_generation_identity"]
    identity["handedness"] = ["right", "right"]; identity["result_handedness"] = ["right", "right"]
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_offset_volume_error():
    reference, measured = _public_v35(); identity = measured["external_cad"][0]["offset_thicken_curvature_sign_selfintersection_repair_thickness_volume_topology_convergence_brep_result_generation_identity"]
    identity["thickened_volume_m3"] = 0.4; identity["result_thickened_volume_m3"] = 0.4; identity["convergence_volumes_m3"][-1] = 0.4; identity["result_convergence_volumes_m3"][-1] = 0.4
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_ap242_invalid_solid():
    row = _source_v35(); identity = row["replay_identity"]["step_ap242_pmi_unit_name_color_occurrence_transform_validity_owner_file_result_generation_identity"]
    identity["solid_validity"] = [[1, True], [2, False]]; identity["decoded_solid_validity"] = [[1, True], [2, False]]
    assert _source_result(row)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_dxf_volume_error():
    row = _source_v35(); identity = row["replay_identity"]["dxf_profile_unit_plane_layer_arc_bulge_loop_winding_extrusion_topology_owner_file_result_generation_identity"]
    identity["extruded_volume_m3"] = 0.01; identity["decoded_extruded_volume_m3"] = 0.01
    assert _source_result(row)["status"] == "needs_attention"
