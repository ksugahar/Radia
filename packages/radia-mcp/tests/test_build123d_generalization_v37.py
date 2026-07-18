from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v36 import _public_v36, _source_v36


_PROMOTED_CASE_IDS = (
    "v37_public_mass_properties_centroid_inertia_principal_axes_placement_density_owner_mismatch",
    "v37_public_shell_offset_thickness_normal_side_open_faces_topology_volume_owner_mismatch",
    "v37_source_sketch_constraint_dof_solver_status_reference_geometry_owner_digest_mismatch",
    "v37_source_topological_naming_edge_face_history_ocp_version_shape_owner_mismatch",
)


def _public_v37():
    reference, measured = _public_v36()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1); generation = "mass-placement-contract-241"
            row["mass_properties_centroid_inertia_principal_axes_placement_density_shape_brep_generation_identity"] = {
                "mass_generation": generation, **{key: generation for key in ("density_generation", "centroid_generation", "inertia_generation", "principal_generation", "placement_generation", "owner_generation", "brep_generation", "result_generation")},
                "density_kg_m3": 7800.0, "result_density_kg_m3": 7800.0, "volume_m3": 1.0e-3, "result_volume_m3": 1.0e-3, "mass_kg": 7.8, "result_mass_kg": 7.8,
                "centroid_world_m": [1.0, 2.0, 3.0], "result_centroid_world_m": [1.0, 2.0, 3.0],
                "inertia_world_kg_m2": [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]], "result_inertia_world_kg_m2": [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
                "principal_moments_kg_m2": [1.0, 2.0, 3.0], "result_principal_moments_kg_m2": [1.0, 2.0, 3.0],
                "principal_axes_world": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], "result_principal_axes_world": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "placement_transform": [[1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 2.0], [0.0, 0.0, 1.0, 3.0], [0.0, 0.0, 0.0, 1.0]], "result_placement_transform": [[1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 2.0], [0.0, 0.0, 1.0, 3.0], [0.0, 0.0, 0.0, 1.0]],
                "shape_owner": "part:placed-solid-241", "result_shape_owner": "part:placed-solid-241", "shape_brep_sha256": suffix * 64, "accepted_shape_brep_sha256": suffix * 64,
            }
            generation = "shell-offset-contract-241"
            row["shell_offset_thickness_normal_side_removed_face_topology_volume_input_owner_brep_generation_identity"] = {
                "shell_generation": generation, **{key: generation for key in ("thickness_generation", "normal_generation", "side_generation", "removed_generation", "topology_generation", "volume_generation", "owner_generation", "brep_generation", "result_generation")},
                "thickness_m": 0.002, "result_thickness_m": 0.002, "offset_side": "inside", "result_offset_side": "inside",
                "face_normals": [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]], "result_face_normals": [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
                "removed_face_ids": [6], "result_removed_face_ids": [6], "wall_topology_v_e_f": [16, 24, 10], "result_wall_topology_v_e_f": [16, 24, 10],
                "analytical_shell_volume_m3": 0.012, "result_analytical_shell_volume_m3": 0.012, "input_owner": "part:shell-input-241", "result_input_owner": "part:shell-input-241",
                "shell_brep_sha256": ("3" if index == 0 else "4") * 64, "accepted_shell_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v37():
    row = _source_v36(); identity = row["replay_identity"]; generation = "sketch-solve-contract-241"
    identity["sketch_constraint_dof_solver_reference_unit_owner_source_result_generation_identity"] = {
        "sketch_generation": generation, **{key: generation for key in ("constraint_generation", "dof_generation", "solver_generation", "reference_generation", "unit_generation", "owner_generation", "source_generation", "result_generation")},
        "constraint_ids": ["horizontal:1", "vertical:2", "distance:3"], "solved_constraint_ids": ["horizontal:1", "vertical:2", "distance:3"], "remaining_dof": 0, "solved_remaining_dof": 0,
        "solver_status": "fully_constrained", "solved_solver_status": "fully_constrained", "reference_geometry_ids": ["axis:x", "origin"], "solved_reference_geometry_ids": ["axis:x", "origin"],
        "length_unit": "m", "solved_length_unit": "m", "sketch_owner": "sketch:profile-241", "solved_sketch_owner": "sketch:profile-241", "sketch_source_sha256": "5" * 64, "solved_sketch_source_sha256": "5" * 64, "sketch_result_sha256": "6" * 64, "accepted_sketch_result_sha256": "6" * 64,
    }
    generation = "toponame-contract-241"
    identity["topological_naming_edge_face_history_ocp_selector_shape_feature_source_brep_generation_identity"] = {
        "toponame_generation": generation, **{key: generation for key in ("edge_generation", "face_generation", "history_generation", "ocp_generation", "selector_generation", "shape_generation", "owner_generation", "source_generation", "brep_generation", "result_generation")},
        "edge_names": ["edge:fillet:0", "edge:fillet:1"], "replayed_edge_names": ["edge:fillet:0", "edge:fillet:1"], "face_names": ["face:top", "face:side"], "replayed_face_names": ["face:top", "face:side"],
        "operation_history": ["box", "fillet", "select:face:top"], "replayed_operation_history": ["box", "fillet", "select:face:top"], "ocp_version": "7.8.1", "replayed_ocp_version": "7.8.1", "selector_result": ["face:top"], "replayed_selector_result": ["face:top"],
        "shape_generation_id": 41, "replayed_shape_generation_id": 41, "feature_owner": "feature:fillet-241", "replayed_feature_owner": "feature:fillet-241", "feature_source_sha256": "7" * 64, "replayed_feature_source_sha256": "7" * 64, "feature_brep_sha256": "8" * 64, "accepted_feature_brep_sha256": "8" * 64,
    }
    return row


def test_v37_positive_contracts():
    reference, measured = _public_v37(); assert _public_result(reference, measured)["status"] == "ok"; assert _source_result(_source_v37())["status"] == "ok"


def test_v37_public_mass_properties_centroid_inertia_principal_axes_placement_density_owner_mismatch():
    reference, measured = _public_v37(); identity = measured["external_cad"][0]["mass_properties_centroid_inertia_principal_axes_placement_density_shape_brep_generation_identity"]
    identity.update({"density_generation": "mass-placement-contract-240", "result_density_kg_m3": -7800.0, "result_mass_kg": -7.8, "result_centroid_world_m": [3.0, 2.0, 1.0], "result_inertia_world_kg_m2": [[-1.0, 2.0]], "result_principal_moments_kg_m2": [3.0, 2.0, 1.0], "result_principal_axes_world": [[-1.0, 0.0, 0.0]], "result_placement_transform": [[1.0]], "result_shape_owner": "stale:shape", "accepted_shape_brep_sha256": "9" * 64})
    result = _public_result(reference, measured); assert result["status"] == "needs_attention"; assert not result["checks"]["placed_mass_properties_use_current_density_mass_centroid_inertia_principal_axes_placement_owner_and_brep"]


def test_v37_public_shell_offset_thickness_normal_side_open_faces_topology_volume_owner_mismatch():
    reference, measured = _public_v37(); identity = measured["external_cad"][0]["shell_offset_thickness_normal_side_removed_face_topology_volume_input_owner_brep_generation_identity"]
    identity.update({"thickness_generation": "shell-offset-contract-240", "result_thickness_m": -0.002, "result_offset_side": "outside", "result_face_normals": [[2.0, 0.0, 0.0]], "result_removed_face_ids": [6, 6], "result_wall_topology_v_e_f": [0, 0, 0], "result_analytical_shell_volume_m3": -0.012, "result_input_owner": "stale:input", "accepted_shell_brep_sha256": "a" * 64})
    result = _public_result(reference, measured); assert result["status"] == "needs_attention"; assert not result["checks"]["shell_offsets_use_current_thickness_side_normals_removed_faces_topology_volume_owner_and_brep"]


def test_v37_source_sketch_constraint_dof_solver_status_reference_geometry_owner_digest_mismatch():
    row = _source_v37(); identity = row["replay_identity"]["sketch_constraint_dof_solver_reference_unit_owner_source_result_generation_identity"]
    identity.update({"constraint_generation": "sketch-solve-contract-240", "solved_constraint_ids": ["distance:old"], "solved_remaining_dof": 3, "solved_solver_status": "under_constrained", "solved_reference_geometry_ids": ["stale"], "solved_length_unit": "mm", "solved_sketch_owner": "stale:sketch", "solved_sketch_source_sha256": "b" * 64, "accepted_sketch_result_sha256": "c" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"; assert not result["checks"]["sketch_solves_use_current_constraints_dof_status_references_units_owner_source_and_result"]


def test_v37_source_topological_naming_edge_face_history_ocp_version_shape_owner_mismatch():
    row = _source_v37(); identity = row["replay_identity"]["topological_naming_edge_face_history_ocp_selector_shape_feature_source_brep_generation_identity"]
    identity.update({"edge_generation": "toponame-contract-240", "replayed_edge_names": ["edge:unknown"], "replayed_face_names": ["face:old"], "replayed_operation_history": ["box", "cut"], "replayed_ocp_version": "7.7.0", "replayed_selector_result": ["face:side"], "replayed_shape_generation_id": 40, "replayed_feature_owner": "stale:feature", "replayed_feature_source_sha256": "d" * 64, "accepted_feature_brep_sha256": "e" * 64})
    result = _source_result(row); assert result["status"] == "needs_attention"; assert not result["checks"]["topological_names_use_current_edges_faces_history_ocp_selector_shape_owner_source_and_brep"]


def test_v37_rejects_self_consistent_mass_density_closure_error():
    reference, measured = _public_v37(); identity = measured["external_cad"][0]["mass_properties_centroid_inertia_principal_axes_placement_density_shape_brep_generation_identity"]; identity["mass_kg"] = identity["result_mass_kg"] = 8.0; assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v37_rejects_self_consistent_duplicate_topological_name():
    row = _source_v37(); identity = row["replay_identity"]["topological_naming_edge_face_history_ocp_selector_shape_feature_source_brep_generation_identity"]; identity["edge_names"] = identity["replayed_edge_names"] = ["edge:fillet:0", "edge:fillet:0"]; assert _source_result(row)["status"] == "needs_attention"
