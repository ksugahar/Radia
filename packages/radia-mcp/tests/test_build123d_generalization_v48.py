from __future__ import annotations

from copy import deepcopy

from radia_mcp.build123d.semantic_cad_identity_v48 import (
    IMPORT,
    LOFT,
    SELECTOR,
    SUPPRESSION,
    validate_public_identity,
    validate_source_identity,
)


PROMOTED_CASE_IDS = {
    "v48_public_feature_suppression_configuration_mass_cache_owner_mismatch",
    "v48_public_loft_profile_orientation_wire_seam_correspondence_self_intersection_owner_mismatch",
    "v48_source_tool_import_unit_inference_layer_color_subshape_mapping_owner_mismatch",
    "v48_source_tool_topology_selector_query_cardinality_witness_feature_history_mismatch",
}


def _payloads() -> tuple[dict[str, object], dict[str, object]]:
    suppression_generation = "feature-config-mass-v48-901"
    loft_generation = "loft-correspondence-v48-901"
    configuration = {"length_mm": 40.0, "hole_enabled": False, "fillet_enabled": True}
    profiles = ["wire:base", "wire:mid", "wire:top"]
    orientation = {profile: 1 for profile in profiles}
    seams = {"wire:base": "vertex:1", "wire:mid": "vertex:5", "wire:top": "vertex:9"}
    correspondence = [["vertex:1", "vertex:5", "vertex:9"], ["vertex:2", "vertex:6", "vertex:10"]]
    row = {
        SUPPRESSION: {
            "generation": suppression_generation,
            "suppression_generation": suppression_generation,
            "configuration_generation": suppression_generation,
            "shape_generation": suppression_generation,
            "cache_generation": suppression_generation,
            "result_generation": suppression_generation,
            "suppressed_features": ["hole:1"],
            "result_suppressed_features": ["hole:1"],
            "configuration": configuration,
            "result_configuration": configuration,
            "shape_sha256": "6" * 64,
            "cache_shape_sha256": "6" * 64,
            "mass_kg": 1.25,
            "cached_mass_kg": 1.25,
            "result_mass_kg": 1.25,
            "owner": "build:feature-config-v48-901",
            "cache_owner": "build:feature-config-v48-901",
            "accepted_owner": "build:feature-config-v48-901",
            "result_sha256": "7" * 64,
            "accepted_result_sha256": "7" * 64,
        },
        LOFT: {
            "generation": loft_generation,
            "profile_generation": loft_generation,
            "orientation_generation": loft_generation,
            "seam_generation": loft_generation,
            "correspondence_generation": loft_generation,
            "diagnostic_generation": loft_generation,
            "result_generation": loft_generation,
            "profile_order": profiles,
            "result_profile_order": profiles,
            "profile_orientation": orientation,
            "result_profile_orientation": orientation,
            "wire_seams": seams,
            "result_wire_seams": seams,
            "profile_correspondence": correspondence,
            "result_profile_correspondence": correspondence,
            "self_intersection_count": 0,
            "result_self_intersection_count": 0,
            "owner": "build:loft-v48-901",
            "accepted_owner": "build:loft-v48-901",
            "result_sha256": "8" * 64,
            "accepted_result_sha256": "8" * 64,
        },
    }
    import_generation = "import-metadata-v48-901"
    selector_generation = "selector-witness-v48-901"
    layers = {"body:1": "mechanical", "body:2": "insulation"}
    colors = {"body:1": [0.8, 0.2, 0.1], "body:2": [0.2, 0.4, 0.9]}
    mapping = {"source:body-1": "body:1", "source:body-2": "body:2"}
    entities = ["face:11", "face:12"]
    witnesses = {"face:11": [0.0, 0.0, 1.0], "face:12": [0.0, 0.0, -1.0]}
    source = {
        "replay_identity": {
            IMPORT: {
                "generation": import_generation,
                "unit_generation": import_generation,
                "metadata_generation": import_generation,
                "subshape_generation": import_generation,
                "result_generation": import_generation,
                "source_revision": "source-rev-v48-901",
                "result_source_revision": "source-rev-v48-901",
                "inferred_unit": "mm",
                "result_inferred_unit": "mm",
                "unit_scale_to_m": 0.001,
                "result_unit_scale_to_m": 0.001,
                "layer_map": layers,
                "result_layer_map": layers,
                "color_map": colors,
                "result_color_map": colors,
                "persistent_subshape_map": mapping,
                "result_persistent_subshape_map": mapping,
                "owner": "import:source-rev-v48-901",
                "accepted_owner": "import:source-rev-v48-901",
                "result_sha256": "9" * 64,
                "accepted_result_sha256": "9" * 64,
            },
            SELECTOR: {
                "generation": selector_generation,
                "query_generation": selector_generation,
                "witness_generation": selector_generation,
                "history_generation": selector_generation,
                "result_generation": selector_generation,
                "query": "faces().filter_by(Axis.Z)",
                "result_query": "faces().filter_by(Axis.Z)",
                "expected_cardinality": 2,
                "result_cardinality": 2,
                "selected_entities": entities,
                "result_selected_entities": entities,
                "witness_points": witnesses,
                "result_witness_points": witnesses,
                "feature_history_owner": "history:feature-v48-901",
                "result_feature_history_owner": "history:feature-v48-901",
                "result_sha256": "a" * 64,
                "accepted_result_sha256": "a" * 64,
            },
        }
    }
    return {"reference": [row], "measured": {"cad": [deepcopy(row)]}}, source


def test_v48_positive_public_and_source_replays_are_accepted() -> None:
    public, source = _payloads()
    assert validate_public_identity(public)["status"] == "ok"
    assert validate_source_identity(source)["status"] == "ok"


def test_v48_configuration_cache_and_loft_mutations_are_rejected() -> None:
    public, _ = _payloads()
    public["reference"][0][SUPPRESSION]["cache_owner"] = "build:feature-config-v48-old"
    public["reference"][0][LOFT]["result_self_intersection_count"] = 1
    result = validate_public_identity(public)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) == {
        "v48_feature_configuration_mass_cache_owner",
        "v48_loft_profile_seam_correspondence_owner",
    }


def test_v48_import_metadata_and_selector_mutations_are_rejected() -> None:
    _, source = _payloads()
    source["replay_identity"][IMPORT]["result_source_revision"] = "source-rev-v48-old"
    source["replay_identity"][SELECTOR]["result_cardinality"] = 1
    result = validate_source_identity(source)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) == {
        "v48_import_unit_metadata_subshape_owner",
        "v48_selector_cardinality_witness_history",
    }
