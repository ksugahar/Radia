from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v38 import _public_v38, _source_v38


_LOFT = "loft_multiwire_hole_correspondence_section_orientation_selfintersection_volume_centroid_euler_owner_brep_generation_identity"
_OFFSET = "offset_shell_signed_thickness_join_repair_area_volume_mass_centroid_inertia_owner_brep_generation_identity"
_HEAL = "occt_heal_tolerance_sew_orientation_degenerate_stablename_roundtrip_owner_digest_generation_identity"
_ASSEMBLY = "assembly_hierarchy_location_joint_axis_collision_quantity_bom_owner_result_generation_identity"
_PROMOTED_CASE_IDS = (
    "v39_public_loft_multiwire_hole_correspondence_section_orientation_volume_topology_mismatch",
    "v39_public_offset_shell_thickness_join_selfintersection_area_mass_owner_mismatch",
    "v39_source_occt_heal_tolerance_sew_orientation_stable_name_roundtrip_mismatch",
    "v39_source_assembly_hierarchy_location_joint_axis_collision_bom_owner_mismatch",
)


def _generation(generation: str, *keys: str) -> dict[str, str]:
    return {key: generation for key in keys}


def _public_v39():
    reference, measured = _public_v38()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "multiwire-loft-271"
            row[_LOFT] = {
                "loft_generation": generation,
                **_generation(generation, "wire_generation", "correspondence_generation", "orientation_generation", "hole_generation", "intersection_generation", "mass_generation", "topology_generation", "owner_generation", "brep_generation", "result_generation"),
                "outer_wire_ids": [11, 21, 31],
                "result_outer_wire_ids": [11, 21, 31],
                "inner_wire_ids": [12, 22, 32],
                "result_inner_wire_ids": [12, 22, 32],
                "section_order": [0, 1, 2],
                "result_section_order": [0, 1, 2],
                "outer_orientation": ["ccw", "ccw", "ccw"],
                "result_outer_orientation": ["ccw", "ccw", "ccw"],
                "inner_orientation": ["cw", "cw", "cw"],
                "result_inner_orientation": ["cw", "cw", "cw"],
                "hole_continuity": [[12, 22], [22, 32]],
                "result_hole_continuity": [[12, 22], [22, 32]],
                "self_intersection_free": True,
                "result_self_intersection_free": True,
                "volume_m3": 2.0e-3,
                "result_volume_m3": 2.0e-3,
                "centroid_m": [0.0, 0.0, 5.0e-2],
                "result_centroid_m": [0.0, 0.0, 5.0e-2],
                "solid_count": 1,
                "result_solid_count": 1,
                "boundary_euler_characteristic": 0,
                "result_boundary_euler_characteristic": 0,
                "shape_owner": "part:multiwire-loft-271",
                "result_shape_owner": "part:multiwire-loft-271",
                "loft_brep_sha256": suffix * 64,
                "accepted_loft_brep_sha256": suffix * 64,
            }
            generation = "offset-shell-271"
            row[_OFFSET] = {
                "offset_generation": generation,
                **_generation(generation, "thickness_generation", "join_generation", "repair_generation", "area_generation", "volume_generation", "mass_generation", "owner_generation", "brep_generation", "result_generation"),
                "signed_thickness_m": -2.0e-3,
                "result_signed_thickness_m": -2.0e-3,
                "offset_direction": "inward",
                "result_offset_direction": "inward",
                "join_mode": "arc",
                "result_join_mode": "arc",
                "self_intersection_detected": True,
                "result_self_intersection_detected": True,
                "self_intersection_repaired": True,
                "result_self_intersection_repaired": True,
                "outer_area_m2": 1.0,
                "result_outer_area_m2": 1.0,
                "inner_area_m2": 0.8,
                "result_inner_area_m2": 0.8,
                "enclosed_volume_m3": 1.6e-3,
                "result_enclosed_volume_m3": 1.6e-3,
                "density_kg_per_m3": 1000.0,
                "result_density_kg_per_m3": 1000.0,
                "mass_kg": 1.6,
                "result_mass_kg": 1.6,
                "centroid_m": [0.0, 0.0, 0.0],
                "result_centroid_m": [0.0, 0.0, 0.0],
                "principal_inertia_kg_m2": [0.01, 0.01, 0.02],
                "result_principal_inertia_kg_m2": [0.01, 0.01, 0.02],
                "shape_owner": "part:offset-shell-271",
                "result_shape_owner": "part:offset-shell-271",
                "offset_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_offset_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v39():
    row = _source_v38()
    identity = row["replay_identity"]
    generation = "occt-heal-271"
    identity[_HEAL] = {
        "heal_generation": generation,
        **_generation(generation, "tolerance_generation", "sew_generation", "orientation_generation", "degenerate_generation", "name_generation", "roundtrip_generation", "owner_generation", "result_generation"),
        "healing_tolerance_m": 1.0e-7,
        "replayed_healing_tolerance_m": 1.0e-7,
        "sewn_shell_count": 1,
        "replayed_sewn_shell_count": 1,
        "solid_count": 1,
        "replayed_solid_count": 1,
        "face_orientations": [1, 1, 1, 1, 1, 1],
        "replayed_face_orientations": [1, 1, 1, 1, 1, 1],
        "removed_degenerate_edge_ids": [91, 92],
        "replayed_removed_degenerate_edge_ids": [91, 92],
        "stable_subshape_names": {"face:top": 5, "face:bottom": 6},
        "replayed_stable_subshape_names": {"face:top": 5, "face:bottom": 6},
        "roundtrip_subshape_counts": {"solid": 1, "shell": 1, "face": 6, "edge": 12},
        "replayed_roundtrip_subshape_counts": {"solid": 1, "shell": 1, "face": 6, "edge": 12},
        "shape_owner": "headless:occt-heal-271",
        "replayed_shape_owner": "headless:occt-heal-271",
        "healed_brep_sha256": "5" * 64,
        "replayed_healed_brep_sha256": "5" * 64,
        "heal_result_sha256": "6" * 64,
        "accepted_heal_result_sha256": "6" * 64,
    }
    generation = "assembly-hierarchy-271"
    identity[_ASSEMBLY] = {
        "assembly_generation": generation,
        **_generation(generation, "hierarchy_generation", "location_generation", "joint_generation", "collision_generation", "quantity_generation", "bom_generation", "owner_generation", "result_generation"),
        "hierarchy": {"root": ["frame", "shaft"], "shaft": ["rotor"]},
        "replayed_hierarchy": {"root": ["frame", "shaft"], "shaft": ["rotor"]},
        "local_locations_m": {"frame": [0.0, 0.0, 0.0], "shaft": [0.0, 0.0, 0.1], "rotor": [0.0, 0.0, 0.2]},
        "replayed_local_locations_m": {"frame": [0.0, 0.0, 0.0], "shaft": [0.0, 0.0, 0.1], "rotor": [0.0, 0.0, 0.2]},
        "global_locations_m": {"frame": [0.0, 0.0, 0.0], "shaft": [0.0, 0.0, 0.1], "rotor": [0.0, 0.0, 0.3]},
        "replayed_global_locations_m": {"frame": [0.0, 0.0, 0.0], "shaft": [0.0, 0.0, 0.1], "rotor": [0.0, 0.0, 0.3]},
        "joint_axis": [0.0, 0.0, 1.0],
        "replayed_joint_axis": [0.0, 0.0, 1.0],
        "collision_pairs": [["frame", "rotor"]],
        "replayed_collision_pairs": [["frame", "rotor"]],
        "component_quantities": {"frame": 1, "shaft": 1, "rotor": 1},
        "replayed_component_quantities": {"frame": 1, "shaft": 1, "rotor": 1},
        "bom_identity": {"frame": "FRAME-001", "shaft": "SHAFT-001", "rotor": "ROTOR-001"},
        "replayed_bom_identity": {"frame": "FRAME-001", "shaft": "SHAFT-001", "rotor": "ROTOR-001"},
        "assembly_owner": "headless:assembly-271",
        "replayed_assembly_owner": "headless:assembly-271",
        "assembly_result_sha256": "7" * 64,
        "accepted_assembly_result_sha256": "7" * 64,
    }
    return row


def test_v39_positive_contracts():
    reference, measured = _public_v39()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v39())["status"] == "ok"


def test_v39_public_loft_multiwire_hole_correspondence_section_orientation_volume_topology_mismatch():
    reference, measured = _public_v39()
    value = measured["external_cad"][0][_LOFT]
    value.update({"correspondence_generation": "multiwire-loft-270", "topology_generation": "multiwire-loft-269", "result_generation": "multiwire-loft-268", "result_outer_wire_ids": [31, 21, 11], "result_inner_wire_ids": [12, 32], "result_section_order": [0, 2, 1], "result_outer_orientation": ["cw", "ccw", "ccw"], "result_inner_orientation": ["ccw", "cw", "cw"], "result_hole_continuity": [[12, 32]], "result_self_intersection_free": False, "result_volume_m3": -2.0e-3, "result_centroid_m": [1.0, 0.0, 0.0], "result_solid_count": 2, "result_boundary_euler_characteristic": 2, "result_shape_owner": "stale:loft", "accepted_loft_brep_sha256": "9" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["multiwire_lofts_use_current_wire_correspondence_orientation_holes_intersection_mass_topology_owner_and_brep"]


def test_v39_public_offset_shell_thickness_join_selfintersection_area_mass_owner_mismatch():
    reference, measured = _public_v39()
    value = measured["external_cad"][0][_OFFSET]
    value.update({"thickness_generation": "offset-shell-270", "mass_generation": "offset-shell-269", "result_generation": "offset-shell-268", "result_signed_thickness_m": 2.0e-3, "result_offset_direction": "outward", "result_join_mode": "intersection", "result_self_intersection_repaired": False, "result_outer_area_m2": 0.8, "result_inner_area_m2": 1.0, "result_enclosed_volume_m3": -1.6e-3, "result_density_kg_per_m3": 1.0, "result_mass_kg": -1.6, "result_centroid_m": [1.0, 0.0, 0.0], "result_principal_inertia_kg_m2": [-0.01, 0.01, 0.02], "result_shape_owner": "stale:offset", "accepted_offset_brep_sha256": "a" * 64})
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["offset_shells_use_current_signed_thickness_join_repair_areas_mass_inertia_owner_and_brep"]


def test_v39_source_occt_heal_tolerance_sew_orientation_stable_name_roundtrip_mismatch():
    row = _source_v39()
    value = row["replay_identity"][_HEAL]
    value.update({"tolerance_generation": "occt-heal-270", "name_generation": "occt-heal-269", "result_generation": "occt-heal-268", "replayed_healing_tolerance_m": 1.0e-2, "replayed_sewn_shell_count": 2, "replayed_solid_count": 0, "replayed_face_orientations": [1, -1], "replayed_removed_degenerate_edge_ids": [92, 99], "replayed_stable_subshape_names": {"face:top": 6}, "replayed_roundtrip_subshape_counts": {"solid": 0, "shell": 2, "face": 5, "edge": 10}, "replayed_shape_owner": "stale:heal", "replayed_healed_brep_sha256": "b" * 64, "accepted_heal_result_sha256": "c" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["occt_heals_use_current_tolerance_sewing_orientation_degenerate_edges_names_roundtrip_owner_and_digests"]


def test_v39_source_assembly_hierarchy_location_joint_axis_collision_bom_owner_mismatch():
    row = _source_v39()
    value = row["replay_identity"][_ASSEMBLY]
    value.update({"hierarchy_generation": "assembly-hierarchy-270", "location_generation": "assembly-hierarchy-269", "result_generation": "assembly-hierarchy-268", "replayed_hierarchy": {"root": ["rotor"]}, "replayed_local_locations_m": {"rotor": [1.0, 0.0, 0.0]}, "replayed_global_locations_m": {"rotor": [0.0, 0.0, 0.2]}, "replayed_joint_axis": [1.0, 1.0, 0.0], "replayed_collision_pairs": [["shaft", "rotor"]], "replayed_component_quantities": {"rotor": 2}, "replayed_bom_identity": {"rotor": "OLD-001"}, "replayed_assembly_owner": "stale:assembly", "accepted_assembly_result_sha256": "d" * 64})
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["assemblies_use_current_hierarchy_locations_joint_axis_collisions_quantities_bom_owner_and_result"]


def test_v39_rejects_self_consistent_loft_hole_discontinuity():
    reference, measured = _public_v39()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_LOFT]["hole_continuity"] = [[12, 32]]
            row[_LOFT]["result_hole_continuity"] = [[12, 32]]
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v39_rejects_self_consistent_offset_mass_error():
    reference, measured = _public_v39()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_OFFSET]["mass_kg"] = 2.0
            row[_OFFSET]["result_mass_kg"] = 2.0
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v39_rejects_self_consistent_heal_orientation_flip():
    row = _source_v39()
    row["replay_identity"][_HEAL]["face_orientations"] = [1, 1, -1, 1]
    row["replay_identity"][_HEAL]["replayed_face_orientations"] = [1, 1, -1, 1]
    assert _source_result(row)["status"] == "needs_attention"


def test_v39_rejects_self_consistent_assembly_location_error():
    row = _source_v39()
    row["replay_identity"][_ASSEMBLY]["global_locations_m"]["rotor"] = [0.0, 0.0, 0.2]
    row["replay_identity"][_ASSEMBLY]["replayed_global_locations_m"]["rotor"] = [0.0, 0.0, 0.2]
    assert _source_result(row)["status"] == "needs_attention"
