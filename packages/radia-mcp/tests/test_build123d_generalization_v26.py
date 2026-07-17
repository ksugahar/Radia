from __future__ import annotations

from test_build123d_generalization_v25 import (
    _public_result,
    _public_v25,
    _source_result,
    _source_v25,
)


def _public_v26():
    reference, measured = _public_v25()
    for rows in (reference, *measured.values()):
        for index, row in enumerate(rows):
            digest = ("1" if index == 0 else "2") * 64
            row["fillet_chamfer_edge_selector_topology_naming_tolerance_shape_generation_identity"] = {
                "feature_generation": "fillet-131", "selector_feature_generation": "fillet-131",
                "topology_feature_generation": "fillet-131", "tolerance_feature_generation": "fillet-131",
                "shape_feature_generation": "fillet-131", "result_feature_generation": "fillet-131",
                "feature_type": "fillet", "result_feature_type": "fillet",
                "edge_selector_names": ["outer-top-1", "outer-top-2"],
                "result_edge_selector_names": ["outer-top-1", "outer-top-2"],
                "persistent_edge_ids": [101, 102], "result_persistent_edge_ids": [101, 102],
                "feature_radius_m": 0.002, "result_feature_radius_m": 0.002,
                "linear_tolerance_m": 1.0e-7, "result_linear_tolerance_m": 1.0e-7,
                "topology_signature": {"solids": 1, "faces": 18, "edges": 36},
                "result_topology_signature": {"solids": 1, "faces": 18, "edges": 36},
                "input_shape_sha256": "3" * 64, "result_input_shape_sha256": "3" * 64,
                "feature_shape_sha256": digest, "result_feature_shape_sha256": digest,
            }
            row["mass_density_center_inertia_reference_frame_assembly_generation_identity"] = {
                "mass_generation": "mass-131", "density_mass_generation": "mass-131",
                "center_mass_generation": "mass-131", "inertia_mass_generation": "mass-131",
                "frame_mass_generation": "mass-131", "assembly_mass_generation": "mass-131",
                "result_mass_generation": "mass-131", "density_kg_m3": 7800.0,
                "result_density_kg_m3": 7800.0, "volume_m3": 0.001, "result_volume_m3": 0.001,
                "mass_kg": 7.8, "result_mass_kg": 7.8,
                "center_of_mass_m": [0.01, 0.02, 0.03], "result_center_of_mass_m": [0.01, 0.02, 0.03],
                "inertia_tensor_kg_m2": [[0.2, 0.01, 0.0], [0.01, 0.3, 0.02], [0.0, 0.02, 0.4]],
                "result_inertia_tensor_kg_m2": [[0.2, 0.01, 0.0], [0.01, 0.3, 0.02], [0.0, 0.02, 0.4]],
                "reference_frame": "assembly-root", "result_reference_frame": "assembly-root",
                "assembly_transform_sha256": "4" * 64, "result_assembly_transform_sha256": "4" * 64,
                "mass_property_sha256": digest, "result_mass_property_sha256": digest,
            }
    return reference, measured


def _source_v26():
    row = _source_v25()
    identity = row["replay_identity"]
    identity["stl_chord_tolerance_triangle_normal_orientation_component_digest_generation_identity"] = {
        "stl_generation": "stl-131", "tolerance_stl_generation": "stl-131",
        "triangle_stl_generation": "stl-131", "normal_stl_generation": "stl-131",
        "component_stl_generation": "stl-131", "result_stl_generation": "stl-131",
        "linear_deflection_m": 0.0001, "decoded_linear_deflection_m": 0.0001,
        "angular_deflection_rad": 0.1, "decoded_angular_deflection_rad": 0.1,
        "triangle_count": 240, "decoded_triangle_count": 240,
        "normal_orientation": "outward", "decoded_normal_orientation": "outward",
        "component_ids": ["base", "rotor", "cover"], "decoded_component_ids": ["base", "rotor", "cover"],
        "source_shape_sha256": "5" * 64, "tessellated_source_shape_sha256": "5" * 64,
        "stl_sha256": "6" * 64, "decoded_stl_sha256": "6" * 64,
    }
    identity["builder_context_workplane_local_frame_part_identity_cache_generation_identity"] = {
        "context_generation": "builder-131", "stack_context_generation": "builder-131",
        "workplane_context_generation": "builder-131", "frame_context_generation": "builder-131",
        "part_context_generation": "builder-131", "cache_context_generation": "builder-131",
        "result_context_generation": "builder-131", "context_stack": ["BuildPart", "BuildSketch", "Locations"],
        "result_context_stack": ["BuildPart", "BuildSketch", "Locations"],
        "workplane_origin_m": [0.0, 0.0, 0.01], "result_workplane_origin_m": [0.0, 0.0, 0.01],
        "workplane_normal": [0.0, 0.0, 1.0], "result_workplane_normal": [0.0, 0.0, 1.0],
        "local_frame_transform": [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.01]],
        "result_local_frame_transform": [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.01]],
        "part_ids": ["base", "rotor", "cover"], "result_part_ids": ["base", "rotor", "cover"],
        "builder_cache_sha256": "7" * 64, "result_builder_cache_sha256": "7" * 64,
        "builder_result_sha256": "8" * 64, "result_builder_result_sha256": "8" * 64,
    }
    return row


def test_v26_positive_public_and_source_identity():
    reference, measured = _public_v26()
    assert _public_result(reference, measured)["status"] == "ok"
    assert _source_result(_source_v26())["status"] == "ok"


def test_v26_public_fillet_chamfer_edge_selector_topology_naming_tolerance_shape_generation_mismatch():
    reference, measured = _public_v26()
    measured["external_cad"][0]["fillet_chamfer_edge_selector_topology_naming_tolerance_shape_generation_identity"].update(
        {"selector_feature_generation": "fillet-130", "result_feature_type": "chamfer",
         "result_edge_selector_names": ["outer-top-2", "stale-edge"], "result_persistent_edge_ids": [102, 999],
         "result_linear_tolerance_m": 1.0e-3, "result_feature_shape_sha256": "b" * 64}
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["fillet_chamfer_features_use_current_selectors_topology_names_tolerance_and_shape"]


def test_v26_public_mass_density_center_inertia_reference_frame_assembly_generation_mismatch():
    reference, measured = _public_v26()
    measured["external_cad"][0]["mass_density_center_inertia_reference_frame_assembly_generation_identity"].update(
        {"density_mass_generation": "mass-130", "result_density_kg_m3": 2700.0,
         "result_mass_kg": 5.4, "result_center_of_mass_m": [0.03, 0.02, 0.01],
         "result_reference_frame": "part-local", "result_mass_property_sha256": "d" * 64}
    )
    result = _public_result(reference, measured)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["mass_properties_use_current_density_center_inertia_frame_and_assembly"]


def test_v26_source_stl_chord_tolerance_triangle_normal_orientation_component_digest_mismatch():
    row = _source_v26()
    row["replay_identity"]["stl_chord_tolerance_triangle_normal_orientation_component_digest_generation_identity"].update(
        {"tolerance_stl_generation": "stl-130", "decoded_linear_deflection_m": 0.001,
         "decoded_triangle_count": 120, "decoded_normal_orientation": "mixed",
         "decoded_component_ids": ["base", "cover", "rotor-old"], "decoded_stl_sha256": "f" * 64}
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["stl_handoff_uses_current_tolerances_triangles_normals_components_and_digests"]


def test_v26_source_builder_context_workplane_local_frame_part_identity_cache_generation_mismatch():
    row = _source_v26()
    row["replay_identity"]["builder_context_workplane_local_frame_part_identity_cache_generation_identity"].update(
        {"stack_context_generation": "builder-130", "result_context_stack": ["BuildSketch", "BuildPart", "Locations"],
         "result_workplane_normal": [0.0, 1.0, 0.0], "result_part_ids": ["base", "cover", "rotor-old"],
         "result_builder_cache_sha256": "0" * 64}
    )
    result = _source_result(row)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["builder_replay_uses_current_context_workplane_frame_parts_and_cache"]
