from __future__ import annotations

from test_build123d_generalization_v42 import _public_v42, _source_v42
from test_build123d_generalization_v29 import _public_result, _source_result
from test_build123d_generalization_v39 import _generation


_GEAR = "involute_gear_module_teeth_pressureangle_backlash_volume_inertia_brep_generation_identity"
_SHEET = "sheetmetal_bend_radius_kfactor_thickness_neutralaxis_volume_brep_generation_identity"
_LOCATION = "location_composition_local_global_rotation_order_subshape_owner_result_generation_identity"
_REBUILD = "parametric_rebuild_dependency_cache_property_invalidation_export_owner_result_generation_identity"
_PROMOTED_CASE_IDS = (
    "v43_public_involutegear_module_teeth_pressureangle_backlash_volume_inertia_brep_mismatch",
    "v43_public_sheetmetal_bend_radius_kfactor_thickness_neutralaxis_volume_brep_mismatch",
    "v43_source_location_composition_local_global_rotation_order_subshape_owner_mismatch",
    "v43_source_parametric_rebuild_dependency_cache_property_invalidation_export_owner_mismatch",
)


def _public_v43():
    reference, measured = _public_v42()
    for rows in [reference, *measured.values()]:
        for index, row in enumerate(rows):
            suffix = str(index + 1)
            generation = "involute-gear-726"
            row[_GEAR] = {
                "gear_generation": generation,
                **_generation(generation, "module_generation", "teeth_generation", "pressure_angle_generation", "backlash_generation", "volume_generation", "inertia_generation", "owner_generation", "brep_generation", "result_generation"),
                "module_m": 0.002, "result_module_m": 0.002,
                "tooth_count": 24, "result_tooth_count": 24,
                "pressure_angle_deg": 20.0, "result_pressure_angle_deg": 20.0,
                "backlash_m": 2.0e-5, "result_backlash_m": 2.0e-5,
                "pitch_diameter_m": 0.048, "result_pitch_diameter_m": 0.048,
                "base_diameter_m": 0.048 * __import__("math").cos(__import__("math").radians(20.0)),
                "result_base_diameter_m": 0.048 * __import__("math").cos(__import__("math").radians(20.0)),
                "gear_volume_m3": 2.4e-5, "result_gear_volume_m3": 2.4e-5,
                "inertia_tensor_kg_m2": [[1.0e-8, 0.0, 0.0], [0.0, 2.0e-8, 0.0], [0.0, 0.0, 3.0e-8]],
                "result_inertia_tensor_kg_m2": [[1.0e-8, 0.0, 0.0], [0.0, 2.0e-8, 0.0], [0.0, 0.0, 3.0e-8]],
                "shape_owner": "part:involute-gear-726", "result_shape_owner": "part:involute-gear-726",
                "gear_brep_sha256": suffix * 64, "accepted_gear_brep_sha256": suffix * 64,
            }
            generation = "sheetmetal-bend-726"
            row[_SHEET] = {
                "sheet_generation": generation,
                **_generation(generation, "bend_radius_generation", "kfactor_generation", "thickness_generation", "neutral_axis_generation", "volume_generation", "owner_generation", "brep_generation", "result_generation"),
                "bend_radius_m": 0.003, "result_bend_radius_m": 0.003,
                "bend_angle_deg": 90.0, "result_bend_angle_deg": 90.0,
                "k_factor": 0.42, "result_k_factor": 0.42,
                "thickness_m": 0.001, "result_thickness_m": 0.001,
                "neutral_axis_radius_m": 0.00342, "result_neutral_axis_radius_m": 0.00342,
                "volume_m3": 1.25e-5, "result_volume_m3": 1.25e-5,
                "surface_area_m2": 0.025, "result_surface_area_m2": 0.025,
                "shape_owner": "part:sheetmetal-726", "result_shape_owner": "part:sheetmetal-726",
                "sheet_brep_sha256": suffix * 64, "accepted_sheet_brep_sha256": suffix * 64,
            }
    return reference, measured


def _source_v43():
    row = _source_v42()
    identity = row["replay_identity"]
    generation = "location-compose-726"
    identity[_LOCATION] = {
        "location_generation": generation,
        **_generation(generation, "local_generation", "global_generation", "rotation_generation", "order_generation", "subshape_generation", "owner_generation", "result_generation"),
        "local_location": {"translation_m": [0.01, 0.0, 0.0], "rotation_deg": [0.0, 0.0, 30.0]},
        "global_location": {"translation_m": [0.01, 0.02, 0.0], "rotation_deg": [0.0, 0.0, 30.0]},
        "replayed_local_location": {"translation_m": [0.01, 0.0, 0.0], "rotation_deg": [0.0, 0.0, 30.0]},
        "replayed_global_location": {"translation_m": [0.01, 0.02, 0.0], "rotation_deg": [0.0, 0.0, 30.0]},
        "rotation_order": ["Z", "Y", "X"], "replayed_rotation_order": ["Z", "Y", "X"],
        "translation_frame": "parent", "replayed_translation_frame": "parent",
        "subshape_id": "face:31", "replayed_subshape_id": "face:31",
        "shape_owner": "headless:location-compose-726", "replayed_shape_owner": "headless:location-compose-726",
        "location_result_sha256": "8" * 64, "accepted_location_result_sha256": "8" * 64,
    }
    generation = "parametric-rebuild-726"
    identity[_REBUILD] = {
        "rebuild_generation": generation,
        **_generation(generation, "dependency_generation", "cache_generation", "property_generation", "invalidation_generation", "topology_generation", "export_generation", "owner_generation", "result_generation"),
        "dependency_values": {"length_m": 0.1, "radius_m": 0.02},
        "replayed_dependency_values": {"length_m": 0.1, "radius_m": 0.02},
        "dependency_identity_sha256": "04d499f658a26ae1536b0756e01599404c8e2694d14fe7149437857a10ad60d0",
        "replayed_dependency_identity_sha256": "04d499f658a26ae1536b0756e01599404c8e2694d14fe7149437857a10ad60d0",
        "cache_key": "length=0.1;radius=0.02", "replayed_cache_key": "length=0.1;radius=0.02",
        "invalidated_properties": ["volume", "surface_area", "center_of_mass"],
        "replayed_invalidated_properties": ["volume", "surface_area", "center_of_mass"],
        "topology_signature": {"solid": 1, "face": 8}, "replayed_topology_signature": {"solid": 1, "face": 8},
        "export_owner": "headless:parametric-rebuild-726", "replayed_export_owner": "headless:parametric-rebuild-726",
        "rebuild_result_sha256": "9" * 64, "accepted_rebuild_result_sha256": "9" * 64,
    }
    return row


def test_v43_positive_contracts():
    reference, measured = _public_v43()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v43())["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 4


def test_v43_public_gear_mismatch():
    reference, measured = _public_v43()
    measured["external_cad"][0][_GEAR]["result_pressure_angle_deg"] = 25.0
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["involute_gears_use_current_module_teeth_pressure_angle_backlash_volume_inertia_owner_and_brep"]


def test_v43_public_sheet_mismatch():
    reference, measured = _public_v43()
    measured["external_cad"][0][_SHEET]["result_k_factor"] = 1.2
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["sheet_metal_bends_use_current_radius_kfactor_thickness_neutral_axis_volume_area_owner_and_brep"]


def test_v43_source_location_mismatch():
    row = _source_v43()
    row["replay_identity"][_LOCATION]["replayed_translation_frame"] = "world"
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["locations_use_current_local_global_composition_rotation_order_subshape_owner_and_result"]


def test_v43_source_rebuild_mismatch():
    row = _source_v43()
    row["replay_identity"][_REBUILD]["replayed_invalidated_properties"] = []
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["parametric_rebuilds_use_current_dependencies_cache_invalidation_topology_owner_and_result"]


def test_v43_rejects_self_consistent_wrong_gear_inertia():
    reference, measured = _public_v43()
    for rows in [reference, *measured.values()]:
        for row in rows:
            wrong = [[1.0e-8, 0.0, 0.0], [0.0, 2.0e-8, 0.0], [0.0, 0.0, -3.0e-8]]
            row[_GEAR]["inertia_tensor_kg_m2"] = wrong
            row[_GEAR]["result_inertia_tensor_kg_m2"] = wrong
    assert _public_result(reference, measured)["status"] == "needs_attention"


def test_v43_rejects_self_consistent_stale_parametric_cache():
    row = _source_v43()
    value = row["replay_identity"][_REBUILD]
    value["cache_key"] = value["replayed_cache_key"] = "length=0.2;radius=0.02"
    value["dependency_values"] = value["replayed_dependency_values"] = {"length_m": 0.2, "radius_m": 0.02}
    value["dependency_identity_sha256"] = value["replayed_dependency_identity_sha256"] = "a" * 64
    assert _source_result(row)["status"] == "needs_attention"
