from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v39 import _generation, _public_v39, _source_v39


_DRAFT = (
    "draft_neutral_plane_pull_direction_signed_angle_face_tangent_topology_"
    "volume_centroid_owner_brep_generation_identity"
)
_THREAD = (
    "thread_pitch_handedness_flank_profile_diameter_runout_selfintersection_"
    "volume_owner_brep_generation_identity"
)
_OBB = (
    "oriented_bounding_box_principal_inertia_axis_com_frame_transform_unit_"
    "owner_result_generation_identity"
)
_TESSELLATION = (
    "tessellation_linear_angular_deflection_triangle_vertex_orientation_"
    "watertight_stl_brep_owner_result_generation_identity"
)
_PROMOTED_CASE_IDS = (
    "v40_public_draft_neutral_plane_pull_direction_angle_face_topology_volume_mismatch",
    "v40_public_thread_pitch_handedness_flank_profile_runout_selfintersection_volume_mismatch",
    "v40_source_oriented_bounding_box_principal_inertia_axis_transform_unit_owner_mismatch",
    "v40_source_tessellation_linear_angular_deflection_triangle_orientation_stl_owner_mismatch",
)


def _public_v40():
    reference, measured = _public_v39()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "draft-solid-311"
            row[_DRAFT] = {
                "draft_generation": generation,
                **_generation(
                    generation,
                    "plane_generation",
                    "pull_generation",
                    "angle_generation",
                    "face_generation",
                    "tangent_generation",
                    "topology_generation",
                    "mass_generation",
                    "owner_generation",
                    "brep_generation",
                    "result_generation",
                ),
                "neutral_plane_origin_m": [0.0, 0.0, 0.0],
                "result_neutral_plane_origin_m": [0.0, 0.0, 0.0],
                "neutral_plane_normal": [0.0, 0.0, 1.0],
                "result_neutral_plane_normal": [0.0, 0.0, 1.0],
                "pull_direction": [0.0, 0.0, 1.0],
                "result_pull_direction": [0.0, 0.0, 1.0],
                "signed_draft_angle_rad": 0.05235987755982988,
                "result_signed_draft_angle_rad": 0.05235987755982988,
                "selected_face_ids": [1, 2, 3, 4],
                "result_selected_face_ids": [1, 2, 3, 4],
                "tangent_continuity": True,
                "result_tangent_continuity": True,
                "topology_signature": {
                    "solid": 1,
                    "shell": 1,
                    "face": 10,
                    "edge": 24,
                    "vertex": 16,
                },
                "result_topology_signature": {
                    "solid": 1,
                    "shell": 1,
                    "face": 10,
                    "edge": 24,
                    "vertex": 16,
                },
                "volume_m3": 1.9e-3,
                "result_volume_m3": 1.9e-3,
                "centroid_m": [0.0, 0.0, 5.0e-2],
                "result_centroid_m": [0.0, 0.0, 5.0e-2],
                "shape_owner": "part:draft-solid-311",
                "result_shape_owner": "part:draft-solid-311",
                "draft_brep_sha256": suffix * 64,
                "accepted_draft_brep_sha256": suffix * 64,
            }

            generation = "modeled-thread-311"
            row[_THREAD] = {
                "thread_generation": generation,
                **_generation(
                    generation,
                    "pitch_generation",
                    "handedness_generation",
                    "profile_generation",
                    "diameter_generation",
                    "runout_generation",
                    "intersection_generation",
                    "mass_generation",
                    "owner_generation",
                    "brep_generation",
                    "result_generation",
                ),
                "pitch_m": 2.0e-3,
                "result_pitch_m": 2.0e-3,
                "handedness": "right",
                "result_handedness": "right",
                "flank_angle_rad": 1.0471975511965976,
                "result_flank_angle_rad": 1.0471975511965976,
                "profile_type": "iso_v",
                "result_profile_type": "iso_v",
                "major_diameter_m": 2.0e-2,
                "result_major_diameter_m": 2.0e-2,
                "minor_diameter_m": 1.70e-2,
                "result_minor_diameter_m": 1.70e-2,
                "thread_length_m": 2.0e-2,
                "result_thread_length_m": 2.0e-2,
                "turn_count": 10.0,
                "result_turn_count": 10.0,
                "runout_length_m": 4.0e-3,
                "result_runout_length_m": 4.0e-3,
                "self_intersection_free": True,
                "result_self_intersection_free": True,
                "volume_m3": 5.0e-6,
                "result_volume_m3": 5.0e-6,
                "shape_owner": "part:modeled-thread-311",
                "result_shape_owner": "part:modeled-thread-311",
                "thread_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_thread_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v40():
    row = _source_v39()
    identity = row["replay_identity"]
    generation = "obb-inertia-311"
    identity[_OBB] = {
        "obb_generation": generation,
        **_generation(
            generation,
            "box_generation",
            "inertia_generation",
            "axis_generation",
            "com_generation",
            "transform_generation",
            "unit_generation",
            "owner_generation",
            "result_generation",
        ),
        "obb_center_m": [1.0, 2.0, 3.0],
        "replayed_obb_center_m": [1.0, 2.0, 3.0],
        "obb_half_extents_m": [0.5, 0.25, 0.125],
        "replayed_obb_half_extents_m": [0.5, 0.25, 0.125],
        "obb_axes": [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        "replayed_obb_axes": [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        "principal_moments_kg_m2": [0.01, 0.02, 0.025],
        "replayed_principal_moments_kg_m2": [0.01, 0.02, 0.025],
        "principal_axes": [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        "replayed_principal_axes": [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        "center_of_mass_local_m": [0.0, 0.0, 0.0],
        "replayed_center_of_mass_local_m": [0.0, 0.0, 0.0],
        "center_of_mass_world_m": [1.0, 2.0, 3.0],
        "replayed_center_of_mass_world_m": [1.0, 2.0, 3.0],
        "local_to_world_transform": [
            [0.0, -1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 2.0],
            [0.0, 0.0, 1.0, 3.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        "replayed_local_to_world_transform": [
            [0.0, -1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 2.0],
            [0.0, 0.0, 1.0, 3.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        "length_unit": "m",
        "replayed_length_unit": "m",
        "inertia_unit": "kg*m^2",
        "replayed_inertia_unit": "kg*m^2",
        "shape_owner": "headless:obb-inertia-311",
        "replayed_shape_owner": "headless:obb-inertia-311",
        "obb_result_sha256": "5" * 64,
        "accepted_obb_result_sha256": "5" * 64,
    }

    generation = "tessellation-311"
    vertices = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ]
    identity[_TESSELLATION] = {
        "tessellation_generation": generation,
        **_generation(
            generation,
            "linear_generation",
            "angular_generation",
            "triangle_generation",
            "vertex_generation",
            "orientation_generation",
            "watertight_generation",
            "stl_generation",
            "brep_generation",
            "owner_generation",
            "result_generation",
        ),
        "linear_deflection_m": 1.0e-4,
        "replayed_linear_deflection_m": 1.0e-4,
        "angular_deflection_rad": 0.2,
        "replayed_angular_deflection_rad": 0.2,
        "triangle_count": 12,
        "replayed_triangle_count": 12,
        "vertex_coordinates_m": vertices,
        "replayed_vertex_coordinates_m": vertices,
        "outward_orientation": True,
        "replayed_outward_orientation": True,
        "watertight": True,
        "replayed_watertight": True,
        "nonmanifold_edge_count": 0,
        "replayed_nonmanifold_edge_count": 0,
        "stl_owner": "headless:tessellation-311",
        "replayed_stl_owner": "headless:tessellation-311",
        "source_brep_sha256": "6" * 64,
        "replayed_source_brep_sha256": "6" * 64,
        "stl_sha256": "7" * 64,
        "replayed_stl_sha256": "7" * 64,
        "tessellation_result_sha256": "8" * 64,
        "accepted_tessellation_result_sha256": "8" * 64,
    }
    return row


def test_v40_positive_contracts():
    reference, measured = _public_v40()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v40())["status"] == "ok"


def test_v40_public_draft_mismatch():
    reference, measured = _public_v40()
    value = measured["external_cad"][0][_DRAFT]
    value.update(
        {
            "plane_generation": "draft-solid-310",
            "topology_generation": "draft-solid-309",
            "result_generation": "draft-solid-308",
            "result_neutral_plane_origin_m": [0.0, 0.0, 1.0],
            "result_neutral_plane_normal": [0.0, 1.0, 0.0],
            "result_pull_direction": [0.0, 0.0, -1.0],
            "result_signed_draft_angle_rad": -0.05235987755982988,
            "result_selected_face_ids": [1, 3],
            "result_tangent_continuity": False,
            "result_topology_signature": {"solid": 2, "face": 8},
            "result_volume_m3": -1.9e-3,
            "result_centroid_m": [1.0, 0.0, 0.0],
            "result_shape_owner": "stale:draft",
            "accepted_draft_brep_sha256": "9" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "drafts_use_current_neutral_plane_pull_angle_faces_tangency_topology_mass_owner_and_brep"
    ]


def test_v40_public_thread_mismatch():
    reference, measured = _public_v40()
    value = measured["external_cad"][0][_THREAD]
    value.update(
        {
            "pitch_generation": "modeled-thread-310",
            "profile_generation": "modeled-thread-309",
            "result_generation": "modeled-thread-308",
            "result_pitch_m": 1.0e-3,
            "result_handedness": "left",
            "result_flank_angle_rad": 0.5,
            "result_profile_type": "square",
            "result_major_diameter_m": 1.6e-2,
            "result_minor_diameter_m": 2.1e-2,
            "result_thread_length_m": 1.5e-2,
            "result_turn_count": 4.0,
            "result_runout_length_m": 3.0e-2,
            "result_self_intersection_free": False,
            "result_volume_m3": -5.0e-6,
            "result_shape_owner": "stale:thread",
            "accepted_thread_brep_sha256": "a" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "threads_use_current_pitch_handedness_profile_diameters_runout_intersection_mass_owner_and_brep"
    ]


def test_v40_source_obb_inertia_mismatch():
    row = _source_v40()
    value = row["replay_identity"][_OBB]
    value.update(
        {
            "box_generation": "obb-inertia-310",
            "transform_generation": "obb-inertia-309",
            "result_generation": "obb-inertia-308",
            "replayed_obb_center_m": [3.0, 2.0, 1.0],
            "replayed_obb_half_extents_m": [0.5, -0.25, 0.125],
            "replayed_obb_axes": [[1.0, 0.0, 0.0]] * 3,
            "replayed_principal_moments_kg_m2": [0.025, -0.02, 0.01],
            "replayed_principal_axes": [[1.0, 0.0, 0.0]] * 3,
            "replayed_center_of_mass_world_m": [0.0, 0.0, 0.0],
            "replayed_local_to_world_transform": [[1.0, 0.0], [0.0, 1.0]],
            "replayed_length_unit": "mm",
            "replayed_inertia_unit": "g*mm^2",
            "replayed_shape_owner": "stale:obb",
            "accepted_obb_result_sha256": "b" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "mass_property_replays_use_current_obb_principal_inertia_axes_com_transform_units_owner_and_result"
    ]


def test_v40_source_tessellation_mismatch():
    row = _source_v40()
    value = row["replay_identity"][_TESSELLATION]
    value.update(
        {
            "linear_generation": "tessellation-310",
            "orientation_generation": "tessellation-309",
            "result_generation": "tessellation-308",
            "replayed_linear_deflection_m": 1.0e-2,
            "replayed_angular_deflection_rad": 1.2,
            "replayed_triangle_count": 10,
            "replayed_vertex_coordinates_m": [[0.0, 0.0, 0.0]],
            "replayed_outward_orientation": False,
            "replayed_watertight": False,
            "replayed_nonmanifold_edge_count": 4,
            "replayed_stl_owner": "stale:stl",
            "replayed_source_brep_sha256": "c" * 64,
            "replayed_stl_sha256": "d" * 64,
            "accepted_tessellation_result_sha256": "e" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "tessellation_replays_use_current_deflections_triangles_vertices_orientation_watertight_stl_brep_owner_and_result"
    ]


def test_v40_rejects_self_consistent_draft_frame_misalignment():
    reference, measured = _public_v40()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_DRAFT]["pull_direction"] = [1.0, 0.0, 0.0]
            row[_DRAFT]["result_pull_direction"] = [1.0, 0.0, 0.0]
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v40_rejects_self_consistent_thread_length_pitch_error():
    reference, measured = _public_v40()
    for rows in [reference, *measured.values()]:
        for row in rows:
            row[_THREAD]["turn_count"] = 9.0
            row[_THREAD]["result_turn_count"] = 9.0
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v40_rejects_self_consistent_left_handed_obb_axes():
    row = _source_v40()
    axes = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    row["replay_identity"][_OBB]["obb_axes"] = axes
    row["replay_identity"][_OBB]["replayed_obb_axes"] = axes
    assert _source_result(row)["status"] == "needs_attention"


def test_v40_rejects_self_consistent_duplicate_tessellation_vertex():
    row = _source_v40()
    vertices = row["replay_identity"][_TESSELLATION]["vertex_coordinates_m"]
    vertices[-1] = vertices[0]
    row["replay_identity"][_TESSELLATION]["replayed_vertex_coordinates_m"] = vertices
    assert _source_result(row)["status"] == "needs_attention"
