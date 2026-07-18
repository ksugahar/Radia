from __future__ import annotations

import math

from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v37 import _public_v37, _source_v37

_PROMOTED_CASE_IDS = (
    "v38_public_revolve_axis_angle_profile_crossing_orientation_volume_centroid_topology_mismatch",
    "v38_public_involute_gear_module_teeth_pressure_angle_backlash_pitch_volume_owner_mismatch",
    "v38_source_brep_roundtrip_tolerance_ocp_version_subshape_counts_bounds_volume_digest_mismatch",
    "v38_source_svg_path_fillrule_curve_transform_unit_wire_face_extrusion_owner_mismatch",
)


def _public_v38():
    reference, measured = _public_v37()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "revolve-contract-258"
            pappus_volume = 1.0e-2 * 2.0 * math.pi * 5.0e-2
            row[
                "revolve_axis_angle_profile_crossing_orientation_volume_centroid_topology_shape_brep_generation_identity"
            ] = {
                "revolve_generation": generation,
                **{
                    key: generation
                    for key in (
                        "axis_generation",
                        "angle_generation",
                        "profile_generation",
                        "orientation_generation",
                        "volume_generation",
                        "centroid_generation",
                        "topology_generation",
                        "owner_generation",
                        "brep_generation",
                        "result_generation",
                    )
                },
                "axis_origin_m": [0.0, 0.0, 0.0],
                "result_axis_origin_m": [0.0, 0.0, 0.0],
                "axis_direction": [0.0, 0.0, 1.0],
                "result_axis_direction": [0.0, 0.0, 1.0],
                "sweep_angle_deg": 360.0,
                "result_sweep_angle_deg": 360.0,
                "profile_axis_crossing": False,
                "result_profile_axis_crossing": False,
                "profile_loop_orientation": "counterclockwise_outward",
                "result_profile_loop_orientation": "counterclockwise_outward",
                "profile_area_m2": 1.0e-2,
                "result_profile_area_m2": 1.0e-2,
                "profile_centroid_radius_m": 5.0e-2,
                "result_profile_centroid_radius_m": 5.0e-2,
                "analytic_volume_m3": pappus_volume,
                "result_analytic_volume_m3": pappus_volume,
                "volume_tolerance_m3": 1.0e-12,
                "result_volume_tolerance_m3": 1.0e-12,
                "centroid_world_m": [0.0, 0.0, 0.1],
                "result_centroid_world_m": [0.0, 0.0, 0.1],
                "solid_count": 1,
                "result_solid_count": 1,
                "shell_count": 1,
                "result_shell_count": 1,
                "boundary_genus": 1,
                "result_boundary_genus": 1,
                "boundary_euler_characteristic": 0,
                "result_boundary_euler_characteristic": 0,
                "shape_owner": "part:revolve-258",
                "result_shape_owner": "part:revolve-258",
                "revolve_brep_sha256": suffix * 64,
                "accepted_revolve_brep_sha256": suffix * 64,
            }

            generation = "involute-gear-contract-258"
            module_m = 2.0e-3
            teeth = 20
            pressure_angle_deg = 20.0
            pitch_diameter_m = module_m * teeth
            base_diameter_m = pitch_diameter_m * math.cos(math.radians(pressure_angle_deg))
            addendum_diameter_m = pitch_diameter_m + 2.0 * module_m
            row[
                "involute_gear_module_teeth_pressure_backlash_pitch_base_addendum_periodicity_volume_shape_brep_generation_identity"
            ] = {
                "gear_generation": generation,
                **{
                    key: generation
                    for key in (
                        "module_generation",
                        "tooth_generation",
                        "pressure_generation",
                        "backlash_generation",
                        "diameter_generation",
                        "periodicity_generation",
                        "volume_generation",
                        "owner_generation",
                        "brep_generation",
                        "result_generation",
                    )
                },
                "module_m": module_m,
                "result_module_m": module_m,
                "tooth_count": teeth,
                "result_tooth_count": teeth,
                "pressure_angle_deg": pressure_angle_deg,
                "result_pressure_angle_deg": pressure_angle_deg,
                "backlash_m": 1.0e-4,
                "result_backlash_m": 1.0e-4,
                "pitch_diameter_m": pitch_diameter_m,
                "result_pitch_diameter_m": pitch_diameter_m,
                "base_diameter_m": base_diameter_m,
                "result_base_diameter_m": base_diameter_m,
                "addendum_diameter_m": addendum_diameter_m,
                "result_addendum_diameter_m": addendum_diameter_m,
                "tooth_period_angle_deg": 360.0 / teeth,
                "result_tooth_period_angle_deg": 360.0 / teeth,
                "gear_volume_m3": 8.0e-6,
                "result_gear_volume_m3": 8.0e-6,
                "shape_owner": "part:involute-gear-258",
                "result_shape_owner": "part:involute-gear-258",
                "gear_brep_sha256": ("3" if index == 0 else "4") * 64,
                "accepted_gear_brep_sha256": ("3" if index == 0 else "4") * 64,
            }
    return reference, measured


def _source_v38():
    row = _source_v37()
    identity = row["replay_identity"]
    generation = "brep-roundtrip-contract-258"
    identity[
        "brep_roundtrip_tolerance_ocp_version_subshape_count_bounds_volume_shape_owner_source_restored_generation_identity"
    ] = {
        "brep_generation": generation,
        **{
            key: generation
            for key in (
                "tolerance_generation",
                "ocp_generation",
                "topology_generation",
                "bounds_generation",
                "volume_generation",
                "owner_generation",
                "source_generation",
                "restored_generation",
                "result_generation",
            )
        },
        "serialization_tolerance_m": 1.0e-7,
        "restored_serialization_tolerance_m": 1.0e-7,
        "ocp_version": "7.8.1",
        "restored_ocp_version": "7.8.1",
        "subshape_counts": {"solid": 1, "shell": 1, "face": 6, "edge": 12, "vertex": 8},
        "restored_subshape_counts": {"solid": 1, "shell": 1, "face": 6, "edge": 12, "vertex": 8},
        "bounding_box_min_m": [-0.5, -0.25, 0.0],
        "restored_bounding_box_min_m": [-0.5, -0.25, 0.0],
        "bounding_box_max_m": [0.5, 0.25, 0.1],
        "restored_bounding_box_max_m": [0.5, 0.25, 0.1],
        "volume_m3": 5.0e-2,
        "restored_volume_m3": 5.0e-2,
        "shape_owner": "part:brep-roundtrip-258",
        "restored_shape_owner": "part:brep-roundtrip-258",
        "source_brep_sha256": "5" * 64,
        "restored_source_brep_sha256": "5" * 64,
        "restored_brep_sha256": "6" * 64,
        "accepted_restored_brep_sha256": "6" * 64,
    }

    generation = "svg-extrusion-contract-258"
    identity[
        "svg_path_fillrule_curve_transform_unit_wire_face_extrusion_source_digest_generation_identity"
    ] = {
        "svg_generation": generation,
        **{
            key: generation
            for key in (
                "path_generation",
                "fill_generation",
                "transform_generation",
                "unit_generation",
                "wire_generation",
                "face_generation",
                "volume_generation",
                "owner_generation",
                "digest_generation",
                "result_generation",
            )
        },
        "path_commands": ["M", "C", "C", "Z"],
        "replayed_path_commands": ["M", "C", "C", "Z"],
        "fill_rule": "nonzero",
        "replayed_fill_rule": "nonzero",
        "curve_transform": [[1.0, 0.0, 10.0], [0.0, 1.0, 20.0], [0.0, 0.0, 1.0]],
        "replayed_curve_transform": [[1.0, 0.0, 10.0], [0.0, 1.0, 20.0], [0.0, 0.0, 1.0]],
        "document_unit": "mm",
        "replayed_document_unit": "mm",
        "document_to_meter_scale": 1.0e-3,
        "replayed_document_to_meter_scale": 1.0e-3,
        "wire_closed": True,
        "replayed_wire_closed": True,
        "face_orientation": "counterclockwise_positive",
        "replayed_face_orientation": "counterclockwise_positive",
        "profile_area_m2": 2.0e-3,
        "replayed_profile_area_m2": 2.0e-3,
        "extrusion_distance_m": 1.0e-2,
        "replayed_extrusion_distance_m": 1.0e-2,
        "extrusion_volume_m3": 2.0e-5,
        "replayed_extrusion_volume_m3": 2.0e-5,
        "source_owner": "svg:path-logo-258",
        "replayed_source_owner": "svg:path-logo-258",
        "svg_source_sha256": "7" * 64,
        "replayed_svg_source_sha256": "7" * 64,
        "svg_result_sha256": "8" * 64,
        "accepted_svg_result_sha256": "8" * 64,
    }
    return row


def test_v38_positive_contracts():
    reference, measured = _public_v38()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v38())["status"] == "ok"


def test_v38_public_revolve_axis_angle_profile_crossing_orientation_volume_centroid_topology_mismatch():
    reference, measured = _public_v38()
    identity = measured["external_cad"][0][
        "revolve_axis_angle_profile_crossing_orientation_volume_centroid_topology_shape_brep_generation_identity"
    ]
    identity.update(
        {
            "axis_generation": "revolve-contract-257",
            "topology_generation": "revolve-contract-256",
            "result_generation": "revolve-contract-255",
            "result_axis_direction": [1.0, 0.0, 0.0],
            "result_sweep_angle_deg": 180.0,
            "result_profile_axis_crossing": True,
            "result_profile_loop_orientation": "clockwise_inward",
            "result_analytic_volume_m3": -1.0,
            "result_centroid_world_m": [1.0, 0.0, 0.0],
            "result_solid_count": 2,
            "result_boundary_genus": 0,
            "result_boundary_euler_characteristic": 2,
            "result_shape_owner": "stale:revolve",
            "accepted_revolve_brep_sha256": "9" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "revolved_solids_use_current_axis_angle_profile_orientation_pappus_volume_centroid_topology_owner_and_brep"
    ]


def test_v38_public_involute_gear_module_teeth_pressure_angle_backlash_pitch_volume_owner_mismatch():
    reference, measured = _public_v38()
    identity = measured["external_cad"][0][
        "involute_gear_module_teeth_pressure_backlash_pitch_base_addendum_periodicity_volume_shape_brep_generation_identity"
    ]
    identity.update(
        {
            "module_generation": "involute-gear-contract-257",
            "diameter_generation": "involute-gear-contract-256",
            "result_generation": "involute-gear-contract-255",
            "result_module_m": -2.0e-3,
            "result_tooth_count": 19,
            "result_pressure_angle_deg": 45.0,
            "result_backlash_m": -1.0e-4,
            "result_pitch_diameter_m": 0.1,
            "result_base_diameter_m": 0.2,
            "result_addendum_diameter_m": 0.03,
            "result_tooth_period_angle_deg": 20.0,
            "result_gear_volume_m3": -8.0e-6,
            "result_shape_owner": "stale:gear",
            "accepted_gear_brep_sha256": "a" * 64,
        }
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "involute_gears_use_current_module_teeth_pressure_backlash_diameters_periodicity_volume_owner_and_brep"
    ]


def test_v38_source_brep_roundtrip_tolerance_ocp_version_subshape_counts_bounds_volume_digest_mismatch():
    row = _source_v38()
    identity = row["replay_identity"][
        "brep_roundtrip_tolerance_ocp_version_subshape_count_bounds_volume_shape_owner_source_restored_generation_identity"
    ]
    identity.update(
        {
            "tolerance_generation": "brep-roundtrip-contract-257",
            "topology_generation": "brep-roundtrip-contract-256",
            "result_generation": "brep-roundtrip-contract-255",
            "restored_serialization_tolerance_m": 1.0e-2,
            "restored_ocp_version": "7.7.0",
            "restored_subshape_counts": {
                "solid": 0,
                "shell": 1,
                "face": 5,
                "edge": 12,
                "vertex": 8,
            },
            "restored_bounding_box_min_m": [0.5, 0.25, 0.1],
            "restored_bounding_box_max_m": [-0.5, -0.25, 0.0],
            "restored_volume_m3": -1.0,
            "restored_shape_owner": "stale:brep",
            "restored_source_brep_sha256": "b" * 64,
            "accepted_restored_brep_sha256": "c" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "brep_roundtrips_use_current_tolerance_ocp_topology_bounds_volume_owner_source_and_restored_shape"
    ]


def test_v38_source_svg_path_fillrule_curve_transform_unit_wire_face_extrusion_owner_mismatch():
    row = _source_v38()
    identity = row["replay_identity"][
        "svg_path_fillrule_curve_transform_unit_wire_face_extrusion_source_digest_generation_identity"
    ]
    identity.update(
        {
            "path_generation": "svg-extrusion-contract-257",
            "unit_generation": "svg-extrusion-contract-256",
            "result_generation": "svg-extrusion-contract-255",
            "replayed_path_commands": ["M", "L"],
            "replayed_fill_rule": "evenodd",
            "replayed_curve_transform": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            "replayed_document_unit": "px",
            "replayed_document_to_meter_scale": 1.0,
            "replayed_wire_closed": False,
            "replayed_face_orientation": "clockwise_negative",
            "replayed_extrusion_volume_m3": -2.0e-5,
            "replayed_source_owner": "stale:svg",
            "accepted_svg_result_sha256": "d" * 64,
        }
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "svg_extrusions_use_current_paths_fill_transform_units_wire_face_volume_owner_and_digests"
    ]


def test_v38_rejects_self_consistent_non_pappus_volume():
    reference, measured = _public_v38()
    identity = measured["external_cad"][0][
        "revolve_axis_angle_profile_crossing_orientation_volume_centroid_topology_shape_brep_generation_identity"
    ]
    identity["analytic_volume_m3"] = 1.0
    identity["result_analytic_volume_m3"] = 1.0
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v38_rejects_self_consistent_involute_pitch_error():
    reference, measured = _public_v38()
    identity = measured["external_cad"][0][
        "involute_gear_module_teeth_pressure_backlash_pitch_base_addendum_periodicity_volume_shape_brep_generation_identity"
    ]
    identity["pitch_diameter_m"] = 0.1
    identity["result_pitch_diameter_m"] = 0.1
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v38_rejects_self_consistent_brep_euler_error():
    row = _source_v38()
    identity = row["replay_identity"][
        "brep_roundtrip_tolerance_ocp_version_subshape_count_bounds_volume_shape_owner_source_restored_generation_identity"
    ]
    counts = {"solid": 1, "shell": 1, "face": 5, "edge": 12, "vertex": 8}
    identity["subshape_counts"] = counts
    identity["restored_subshape_counts"] = counts
    assert _source_result(row)["status"] == "needs_attention"


def test_v38_rejects_self_consistent_svg_reflection():
    row = _source_v38()
    identity = row["replay_identity"][
        "svg_path_fillrule_curve_transform_unit_wire_face_extrusion_source_digest_generation_identity"
    ]
    reflection = [[-1.0, 0.0, 10.0], [0.0, 1.0, 20.0], [0.0, 0.0, 1.0]]
    identity["curve_transform"] = reflection
    identity["replayed_curve_transform"] = reflection
    assert _source_result(row)["status"] == "needs_attention"
