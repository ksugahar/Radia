from __future__ import annotations

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v33 import _public_v33, _source_v33


_PROMOTED_CASE_IDS = (
    "v34_public_assembly_mate_transform_cycle_frame_closure_mass_inertia_mismatch",
    "v34_public_shell_fillet_topology_euler_thickness_volume_area_inertia_convergence_mismatch",
    "v34_source_step_unit_entity_color_assembly_transform_shape_owner_roundtrip_mismatch",
    "v34_source_brep_periodic_seam_face_orientation_edge_pcurve_manifold_digest_mismatch",
)


def _public_v34():
    reference, measured = _public_v33()
    identity4 = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    rotations = [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    ]
    inertia_local = [
        [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]],
        [[0.4, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.6]],
        [[0.7, 0.0, 0.0], [0.0, 0.8, 0.0], [0.0, 0.0, 0.9]],
    ]
    inertia_global = [
        [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]],
        [[0.5, 0.0, 0.0], [0.0, 0.4, 0.0], [0.0, 0.0, 0.6]],
        [[0.7, 0.0, 0.0], [0.0, 0.8, 0.0], [0.0, 0.0, 0.9]],
    ]
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            suffix = "1" if index == 0 else "2"
            generation = "assembly-mate-211"
            row[
                "assembly_mate_transform_cycle_frame_handedness_mass_center_inertia_owner_shape_result_generation_identity"
            ] = {
                "assembly_generation": generation,
                **{
                    key: generation
                    for key in (
                        "mate_generation",
                        "cycle_generation",
                        "frame_generation",
                        "mass_generation",
                        "center_generation",
                        "inertia_generation",
                        "owner_generation",
                        "shape_generation",
                        "result_generation",
                    )
                },
                "part_ids": ["base", "arm", "tool"],
                "result_part_ids": ["base", "arm", "tool"],
                "mate_edges": [["base", "arm"], ["arm", "tool"], ["tool", "base"]],
                "result_mate_edges": [["base", "arm"], ["arm", "tool"], ["tool", "base"]],
                "mate_cycle_transform": identity4,
                "result_mate_cycle_transform": identity4,
                "frame_determinants": [1.0, 1.0, 1.0],
                "result_frame_determinants": [1.0, 1.0, 1.0],
                "part_masses_kg": [2.0, 3.0, 1.0],
                "result_part_masses_kg": [2.0, 3.0, 1.0],
                "part_centers_m": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
                "result_part_centers_m": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
                "assembly_mass_kg": 6.0,
                "result_assembly_mass_kg": 6.0,
                "assembly_center_of_mass_m": [2.0 / 3.0, 1.0 / 6.0, 0.0],
                "result_assembly_center_of_mass_m": [2.0 / 3.0, 1.0 / 6.0, 0.0],
                "part_rotation_matrices": rotations,
                "result_part_rotation_matrices": rotations,
                "part_inertia_local_kg_m2": inertia_local,
                "result_part_inertia_local_kg_m2": inertia_local,
                "part_inertia_global_kg_m2": inertia_global,
                "result_part_inertia_global_kg_m2": inertia_global,
                "assembly_owner": "assembly/root",
                "result_assembly_owner": "assembly/root",
                "assembly_shape_sha256": suffix * 64,
                "accepted_assembly_shape_sha256": suffix * 64,
            }
            generation = "shell-fillet-211"
            shell_suffix = "3" if index == 0 else "4"
            row[
                "shell_fillet_topology_euler_manifold_thickness_volume_area_inertia_convergence_brep_result_generation_identity"
            ] = {
                "shell_fillet_generation": generation,
                **{
                    key: generation
                    for key in (
                        "topology_generation",
                        "thickness_generation",
                        "volume_generation",
                        "area_generation",
                        "inertia_generation",
                        "convergence_generation",
                        "brep_generation",
                        "result_generation",
                    )
                },
                "vertex_count": 16,
                "result_vertex_count": 16,
                "edge_count": 24,
                "result_edge_count": 24,
                "face_count": 10,
                "result_face_count": 10,
                "euler_characteristic": 2,
                "result_euler_characteristic": 2,
                "edge_face_incidence_counts": [2] * 24,
                "result_edge_face_incidence_counts": [2] * 24,
                "nominal_wall_thickness_m": 0.002,
                "result_nominal_wall_thickness_m": 0.002,
                "wall_thickness_samples_m": [0.002, 0.0020001, 0.0019999],
                "result_wall_thickness_samples_m": [0.002, 0.0020001, 0.0019999],
                "original_volume_m3": 1.0,
                "result_original_volume_m3": 1.0,
                "removed_volume_m3": 0.8,
                "result_removed_volume_m3": 0.8,
                "shell_volume_m3": 0.2,
                "result_shell_volume_m3": 0.2,
                "surface_area_m2": 5.0,
                "result_surface_area_m2": 5.0,
                "inertia_tensor_kg_m2": [[0.2, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.4]],
                "result_inertia_tensor_kg_m2": [[0.2, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.4]],
                "convergence_tolerances_m": [1.0e-4, 1.0e-5, 1.0e-6],
                "result_convergence_tolerances_m": [1.0e-4, 1.0e-5, 1.0e-6],
                "convergence_volumes_m3": [0.2001, 0.200001, 0.2],
                "result_convergence_volumes_m3": [0.2001, 0.200001, 0.2],
                "shell_brep_sha256": shell_suffix * 64,
                "accepted_shell_brep_sha256": shell_suffix * 64,
            }
    return reference, measured


def _source_v34():
    row = _source_v33()
    identity = row["replay_identity"]
    generation = "step-semantic-roundtrip-211"
    identity[
        "step_unit_product_entity_color_assembly_transform_shape_validity_owner_file_result_generation_identity"
    ] = {
        "step_generation": generation,
        **{
            key: generation
            for key in (
                "unit_generation",
                "entity_generation",
                "color_generation",
                "transform_generation",
                "shape_generation",
                "validity_generation",
                "owner_generation",
                "file_generation",
                "result_generation",
            )
        },
        "length_unit": "m",
        "decoded_length_unit": "m",
        "product_entities": [[1, "base"], [2, "arm"]],
        "decoded_product_entities": [[1, "base"], [2, "arm"]],
        "entity_colors_rgb": [[1, [0.2, 0.3, 0.4]], [2, [0.8, 0.1, 0.1]]],
        "decoded_entity_colors_rgb": [[1, [0.2, 0.3, 0.4]], [2, [0.8, 0.1, 0.1]]],
        "assembly_transform_sha256": [[1, "5" * 64], [2, "6" * 64]],
        "decoded_assembly_transform_sha256": [[1, "5" * 64], [2, "6" * 64]],
        "shape_count": 2,
        "decoded_shape_count": 2,
        "solid_validity": [[1, True], [2, True]],
        "decoded_solid_validity": [[1, True], [2, True]],
        "source_owner": "assembly/root",
        "decoded_source_owner": "assembly/root",
        "step_file_sha256": "7" * 64,
        "decoded_step_file_sha256": "7" * 64,
    }
    generation = "brep-periodic-seam-211"
    identity[
        "brep_periodic_face_seam_edge_orientation_pcurve_tolerance_manifold_serializer_shape_result_generation_identity"
    ] = {
        "periodic_brep_generation": generation,
        **{
            key: generation
            for key in (
                "face_generation",
                "seam_generation",
                "orientation_generation",
                "pcurve_generation",
                "tolerance_generation",
                "manifold_generation",
                "serializer_generation",
                "shape_generation",
                "result_generation",
            )
        },
        "periodic_face_ids": [1, 2],
        "decoded_periodic_face_ids": [1, 2],
        "seam_edge_multiplicity": [[1, 11, 2], [2, 12, 2]],
        "decoded_seam_edge_multiplicity": [[1, 11, 2], [2, 12, 2]],
        "face_orientation_signs": [[1, 1], [2, -1]],
        "decoded_face_orientation_signs": [[1, 1], [2, -1]],
        "edge_pcurve_max_deviation_m": [[11, 1.0e-8], [12, 2.0e-8]],
        "decoded_edge_pcurve_max_deviation_m": [[11, 1.0e-8], [12, 2.0e-8]],
        "vertex_tolerances_m": [[101, 1.0e-7], [102, 2.0e-7]],
        "decoded_vertex_tolerances_m": [[101, 1.0e-7], [102, 2.0e-7]],
        "edge_face_incidence_counts": [[11, 2], [12, 2]],
        "decoded_edge_face_incidence_counts": [[11, 2], [12, 2]],
        "serializer_version": "occt-brep-v3",
        "decoded_serializer_version": "occt-brep-v3",
        "periodic_brep_shape_sha256": "8" * 64,
        "decoded_periodic_brep_shape_sha256": "8" * 64,
    }
    return row


def test_v34_positive_public_and_source_identity():
    reference, measured = _public_v34()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v34())["status"] == "ok"


def test_v34_public_assembly_mate_transform_cycle_frame_closure_mass_inertia_mismatch():
    reference, measured = _public_v34()
    identity = measured["external_cad"][0][
        "assembly_mate_transform_cycle_frame_handedness_mass_center_inertia_owner_shape_result_generation_identity"
    ]
    identity.update(
        {
            "mate_generation": "assembly-mate-210",
            "mass_generation": "assembly-mate-209",
            "result_generation": "assembly-mate-208",
            "result_part_ids": ["base", "tool"],
            "result_mate_edges": [["base", "arm"], ["tool", "base"]],
            "result_mate_cycle_transform": [[1.0, 0.0, 0.0, 0.1], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, -1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
            "result_frame_determinants": [1.0, -1.0, 1.0],
            "result_part_masses_kg": [2.0, 3.0, 2.0],
            "result_part_centers_m": [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 1.0, 0.0]],
            "result_assembly_mass_kg": 5.0,
            "result_assembly_center_of_mass_m": [0.0, 0.0, 0.0],
            "result_part_inertia_global_kg_m2": identity["part_inertia_local_kg_m2"],
            "result_assembly_owner": "stale/root",
            "accepted_assembly_shape_sha256": "9" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "assembly_mates_close_current_cycles_frames_mass_centers_rotated_inertia_owner_and_shape"
    ]


def test_v34_public_shell_fillet_topology_euler_thickness_volume_area_inertia_convergence_mismatch():
    reference, measured = _public_v34()
    identity = measured["external_cad"][0][
        "shell_fillet_topology_euler_manifold_thickness_volume_area_inertia_convergence_brep_result_generation_identity"
    ]
    identity.update(
        {
            "topology_generation": "shell-fillet-210",
            "volume_generation": "shell-fillet-209",
            "result_generation": "shell-fillet-208",
            "result_vertex_count": 15,
            "result_edge_count": 22,
            "result_face_count": 8,
            "result_euler_characteristic": 1,
            "result_edge_face_incidence_counts": [2] * 20 + [1, 3],
            "result_nominal_wall_thickness_m": 0.004,
            "result_wall_thickness_samples_m": [0.001, 0.004, 0.006],
            "result_removed_volume_m3": 0.7,
            "result_shell_volume_m3": 0.4,
            "result_surface_area_m2": 4.0,
            "result_inertia_tensor_kg_m2": [[0.2, 0.1, 0.0], [0.0, -0.3, 0.0], [0.0, 0.0, 0.4]],
            "result_convergence_tolerances_m": [1.0e-4, 1.0e-3, 1.0e-2],
            "result_convergence_volumes_m3": [0.2, 0.3, 0.5],
            "accepted_shell_brep_sha256": "a" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "shell_fillets_use_current_euler_manifold_thickness_volume_area_inertia_convergence_and_brep"
    ]


def test_v34_source_step_unit_entity_color_assembly_transform_shape_owner_roundtrip_mismatch():
    row = _source_v34()
    row["replay_identity"][
        "step_unit_product_entity_color_assembly_transform_shape_validity_owner_file_result_generation_identity"
    ].update(
        {
            "unit_generation": "step-semantic-roundtrip-210",
            "entity_generation": "step-semantic-roundtrip-209",
            "result_generation": "step-semantic-roundtrip-208",
            "decoded_length_unit": "mm",
            "decoded_product_entities": [[1, "old-base"]],
            "decoded_entity_colors_rgb": [[1, [0.0, 0.0, 0.0]]],
            "decoded_assembly_transform_sha256": [[1, "b" * 64]],
            "decoded_shape_count": 1,
            "decoded_solid_validity": [[1, False]],
            "decoded_source_owner": "stale/root",
            "decoded_step_file_sha256": "c" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "step_roundtrips_use_current_units_products_colors_transforms_shapes_validity_owner_and_file"
    ]


def test_v34_source_brep_periodic_seam_face_orientation_edge_pcurve_manifold_digest_mismatch():
    row = _source_v34()
    row["replay_identity"][
        "brep_periodic_face_seam_edge_orientation_pcurve_tolerance_manifold_serializer_shape_result_generation_identity"
    ].update(
        {
            "seam_generation": "brep-periodic-seam-210",
            "serializer_generation": "brep-periodic-seam-209",
            "result_generation": "brep-periodic-seam-208",
            "decoded_periodic_face_ids": [1],
            "decoded_seam_edge_multiplicity": [[1, 11, 1], [2, 12, 3]],
            "decoded_face_orientation_signs": [[1, -1], [2, 1]],
            "decoded_edge_pcurve_max_deviation_m": [[11, 1.0e-2]],
            "decoded_vertex_tolerances_m": [[101, 1.0e-2]],
            "decoded_edge_face_incidence_counts": [[11, 1], [12, 3]],
            "decoded_serializer_version": "occt-brep-v2",
            "decoded_periodic_brep_shape_sha256": "d" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "periodic_brep_roundtrips_use_current_seams_orientations_pcurves_tolerances_manifold_serializer_and_shape"
    ]


def test_v34_rejects_self_consistent_nonclosing_mate_cycle():
    reference, measured = _public_v34()
    identity = measured["external_cad"][0][
        "assembly_mate_transform_cycle_frame_handedness_mass_center_inertia_owner_shape_result_generation_identity"
    ]
    cycle = [[1.0, 0.0, 0.0, 0.1], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    identity["mate_cycle_transform"] = cycle
    identity["result_mate_cycle_transform"] = cycle
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_shell_volume_without_subtraction_closure():
    reference, measured = _public_v34()
    identity = measured["external_cad"][0][
        "shell_fillet_topology_euler_manifold_thickness_volume_area_inertia_convergence_brep_result_generation_identity"
    ]
    identity["shell_volume_m3"] = 0.3
    identity["result_shell_volume_m3"] = 0.3
    identity["convergence_volumes_m3"][-1] = 0.3
    identity["result_convergence_volumes_m3"][-1] = 0.3
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_step_invalid_solid():
    row = _source_v34()
    identity = row["replay_identity"][
        "step_unit_product_entity_color_assembly_transform_shape_validity_owner_file_result_generation_identity"
    ]
    identity["solid_validity"] = [[1, True], [2, False]]
    identity["decoded_solid_validity"] = [[1, True], [2, False]]
    assert _source_result(row)["status"] == "needs_attention"


def test_v34_rejects_self_consistent_nonmanifold_periodic_seam():
    row = _source_v34()
    identity = row["replay_identity"][
        "brep_periodic_face_seam_edge_orientation_pcurve_tolerance_manifold_serializer_shape_result_generation_identity"
    ]
    identity["edge_face_incidence_counts"] = [[11, 1], [12, 3]]
    identity["decoded_edge_face_incidence_counts"] = [[11, 1], [12, 3]]
    assert _source_result(row)["status"] == "needs_attention"
