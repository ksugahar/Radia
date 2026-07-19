from __future__ import annotations

import math

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v39 import _generation
from test_build123d_generalization_v41 import _public_v41, _source_v41


_HELIX = (
    "helicalsweep_pitch_turns_profile_frame_frenet_torsion_validity_volume_"
    "centroid_owner_brep_generation_identity"
)
_BOOLEAN = (
    "boolean_fuzzy_tolerance_sliver_topology_validity_volume_surfacearea_"
    "owner_brep_generation_identity"
)
_SELECTOR = (
    "selector_cache_topology_renumber_geometry_predicate_feature_parent_"
    "owner_result_generation_identity"
)
_STEP = (
    "step_assembly_unit_color_name_transform_shape_hierarchy_export_owner_"
    "file_generation_identity"
)
_PROMOTED_CASE_IDS = (
    "v42_public_helicalsweep_pitch_turns_profile_frame_torsion_volume_centroid_brep_mismatch",
    "v42_public_boolean_fuzzy_tolerance_sliver_face_topology_volume_surfacearea_brep_mismatch",
    "v42_source_selector_cache_topology_renumber_geometry_generation_feature_owner_mismatch",
    "v42_source_step_assembly_units_colors_names_transforms_partshape_export_owner_mismatch",
)


def _public_v42():
    reference, measured = _public_v41()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "helical-sweep-725"
            radius, pitch, turns, area = 0.01, 0.02, 5.0, 1.0e-4
            rise = pitch * turns
            path = turns * math.sqrt((2.0 * math.pi * radius) ** 2 + pitch**2)
            b = pitch / (2.0 * math.pi)
            torsion = b / (radius**2 + b**2)
            row[_HELIX] = {
                "helical_generation": generation,
                **_generation(generation, "pitch_generation", "turn_generation", "profile_generation", "frame_generation", "torsion_generation", "validity_generation", "mass_generation", "owner_generation", "brep_generation", "result_generation"),
                "helix_radius_m": radius, "result_helix_radius_m": radius,
                "pitch_m": pitch, "result_pitch_m": pitch,
                "turns": turns, "result_turns": turns,
                "axial_rise_m": rise, "result_axial_rise_m": rise,
                "profile_area_m2": area, "result_profile_area_m2": area,
                "profile_frame": "frenet", "result_profile_frame": "frenet",
                "frenet_transport": True, "result_frenet_transport": True,
                "path_length_m": path, "result_path_length_m": path,
                "torsion_per_m": torsion, "result_torsion_per_m": torsion,
                "solid_valid": True, "result_solid_valid": True,
                "volume_m3": area * path, "result_volume_m3": area * path,
                "centroid_m": [0.0, 0.0, rise / 2.0], "result_centroid_m": [0.0, 0.0, rise / 2.0],
                "shape_owner": "part:helical-sweep-725", "result_shape_owner": "part:helical-sweep-725",
                "helical_brep_sha256": suffix * 64, "accepted_helical_brep_sha256": suffix * 64,
            }
            generation = "fuzzy-boolean-725"
            topology = {"solid": 1, "shell": 1, "face": 10, "edge": 24, "vertex": 16}
            row[_BOOLEAN] = {
                "boolean_generation": generation,
                **_generation(generation, "tolerance_generation", "sliver_generation", "topology_generation", "validity_generation", "mass_generation", "owner_generation", "brep_generation", "result_generation"),
                "operation": "cut", "result_operation": "cut",
                "fuzzy_tolerance_m": 1.0e-7, "result_fuzzy_tolerance_m": 1.0e-7,
                "minimum_feature_size_m": 1.0e-4, "result_minimum_feature_size_m": 1.0e-4,
                "sliver_face_count_before": 2, "result_sliver_face_count_before": 2,
                "sliver_face_count_after": 0, "result_sliver_face_count_after": 0,
                "topology_signature": topology, "result_topology_signature": topology,
                "solid_valid": True, "result_solid_valid": True,
                "volume_m3": 9.9e-4, "result_volume_m3": 9.9e-4,
                "surface_area_m2": 6.1e-2, "result_surface_area_m2": 6.1e-2,
                "shape_owner": "part:fuzzy-boolean-725", "result_shape_owner": "part:fuzzy-boolean-725",
                "boolean_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_boolean_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v42():
    row = _source_v41()
    identity = row["replay_identity"]
    generation = "selector-cache-725"
    renumber = {"face:11": "face:21", "face:12": "face:22"}
    identity[_SELECTOR] = {
        "selector_generation": generation,
        **_generation(generation, "topology_generation", "renumber_generation", "geometry_generation", "predicate_generation", "feature_generation", "owner_generation", "result_generation"),
        "topology_renumber_map": renumber, "replayed_topology_renumber_map": renumber,
        "geometry_generation_id": 725, "replayed_geometry_generation_id": 725,
        "selector_predicate": "Axis.Z and Area>1e-4", "replayed_selector_predicate": "Axis.Z and Area>1e-4",
        "selected_feature_ids": ["face:21", "face:22"], "replayed_selected_feature_ids": ["face:21", "face:22"],
        "parent_shape_owner": "headless:selector-parent-725", "replayed_parent_shape_owner": "headless:selector-parent-725",
        "selector_result_sha256": "5" * 64, "accepted_selector_result_sha256": "5" * 64,
    }
    generation = "step-assembly-725"
    transforms = {
        "part:1": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
        "part:2": [[1, 0, 0, 100], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
    }
    names = {"part:1": "base", "part:2": "slider"}
    colors = {"part:1": [0.2, 0.4, 0.8], "part:2": [0.8, 0.4, 0.2]}
    shape_ids = {"part:1": "shape:101", "part:2": "shape:102"}
    hierarchy = {"assembly:1": ["part:1", "part:2"]}
    identity[_STEP] = {
        "assembly_generation": generation,
        **_generation(generation, "unit_generation", "color_generation", "name_generation", "transform_generation", "shape_generation", "hierarchy_generation", "owner_generation", "file_generation", "result_generation"),
        "length_unit": "mm", "replayed_length_unit": "mm",
        "unit_scale_to_m": 1.0e-3, "replayed_unit_scale_to_m": 1.0e-3,
        "part_names": names, "replayed_part_names": names,
        "part_colors_rgb": colors, "replayed_part_colors_rgb": colors,
        "part_transforms_in_source_units": transforms, "replayed_part_transforms_in_source_units": transforms,
        "part_shape_ids": shape_ids, "replayed_part_shape_ids": shape_ids,
        "assembly_hierarchy": hierarchy, "replayed_assembly_hierarchy": hierarchy,
        "export_owner": "headless:step-assembly-725", "replayed_export_owner": "headless:step-assembly-725",
        "step_file_sha256": "6" * 64, "replayed_step_file_sha256": "6" * 64,
        "assembly_result_sha256": "7" * 64, "accepted_assembly_result_sha256": "7" * 64,
    }
    return row


def test_v42_positive_contracts():
    reference, measured = _public_v42()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v42())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 4


def test_v42_public_helical_sweep_mismatch():
    reference, measured = _public_v42()
    measured["external_cad"][0][_HELIX].update({"frame_generation": "helical-sweep-724", "result_pitch_m": -1.0, "result_turns": 0.0, "result_axial_rise_m": -1.0, "result_profile_frame": "fixed", "result_frenet_transport": False, "result_torsion_per_m": -1.0, "result_solid_valid": False, "result_volume_m3": -1.0, "result_centroid_m": [9.0], "result_shape_owner": "stale:helix", "accepted_helical_brep_sha256": "a" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["helical_sweeps_use_current_pitch_turns_frame_transport_torsion_validity_mass_owner_and_brep"]


def test_v42_public_fuzzy_boolean_mismatch():
    reference, measured = _public_v42()
    measured["external_cad"][0][_BOOLEAN].update({"sliver_generation": "fuzzy-boolean-724", "result_fuzzy_tolerance_m": 1.0e-2, "result_sliver_face_count_after": 9, "result_topology_signature": {"solid": 0}, "result_solid_valid": False, "result_volume_m3": -1.0, "result_surface_area_m2": -1.0, "result_shape_owner": "stale:boolean", "accepted_boolean_brep_sha256": "b" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["fuzzy_booleans_use_current_tolerance_slivers_topology_validity_mass_owner_and_brep"]


def test_v42_source_selector_cache_mismatch():
    row = _source_v42()
    row["replay_identity"][_SELECTOR].update({"renumber_generation": "selector-cache-724", "replayed_topology_renumber_map": {}, "replayed_geometry_generation_id": 724, "replayed_selector_predicate": "Area<0", "replayed_selected_feature_ids": ["face:11"], "replayed_parent_shape_owner": "stale:parent", "accepted_selector_result_sha256": "c" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["selector_caches_use_current_renumbering_geometry_predicate_features_parent_owner_and_result"]


def test_v42_source_step_assembly_mismatch():
    row = _source_v42()
    row["replay_identity"][_STEP].update({"unit_generation": "step-assembly-724", "replayed_length_unit": "m", "replayed_unit_scale_to_m": 1.0, "replayed_part_names": {"part:1": "old"}, "replayed_part_colors_rgb": {"part:1": [2.0, -1.0, 0.0]}, "replayed_part_transforms_in_source_units": {}, "replayed_part_shape_ids": {}, "replayed_assembly_hierarchy": {}, "replayed_export_owner": "stale:export", "replayed_step_file_sha256": "d" * 64, "accepted_assembly_result_sha256": "e" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["step_assembly_exports_use_current_units_colors_names_transforms_shapes_hierarchy_owner_and_file"]


def test_v42_rejects_self_consistent_wrong_helical_path_length():
    reference, measured = _public_v42()
    for rows in [reference, *measured.values()]:
        for row in rows:
            value = row[_HELIX]
            wrong = 2.0 * value["path_length_m"]
            value["path_length_m"] = value["result_path_length_m"] = wrong
            value["volume_m3"] = value["result_volume_m3"] = value["profile_area_m2"] * wrong
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v42_rejects_self_consistent_oversized_fuzzy_tolerance():
    reference, measured = _public_v42()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_BOOLEAN]["fuzzy_tolerance_m"] = row[_BOOLEAN]["result_fuzzy_tolerance_m"] = 5.0e-5
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v42_rejects_self_consistent_selector_using_pre_renumber_features():
    row = _source_v42()
    value = row["replay_identity"][_SELECTOR]
    value["selected_feature_ids"] = value["replayed_selected_feature_ids"] = ["face:11", "face:12"]
    assert _source_result(row)["status"] == "needs_attention"


def test_v42_rejects_self_consistent_nonrigid_step_transform():
    row = _source_v42()
    value = row["replay_identity"][_STEP]
    transform = [[2, 0, 0, 100], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    value["part_transforms_in_source_units"]["part:2"] = transform
    value["replayed_part_transforms_in_source_units"]["part:2"] = transform
    assert _source_result(row)["status"] == "needs_attention"
