from copy import deepcopy

from radia_mcp.build123d.assembly_replay_identity_v49 import ASSEMBLY, BOOLEAN, SKETCH, STEP, validate_public_identity, validate_source_identity


PROMOTED_CASE_IDS = {
    "v49_public_assembly_mass_density_material_occurrence_transform_suppression_owner_mismatch",
    "v49_public_sketch_constraint_dof_plane_unit_profile_wire_owner_mismatch",
    "v49_source_tool_step_schema_assembly_color_layer_unit_tolerance_owner_mismatch",
    "v49_source_tool_boolean_deleted_subshape_selector_adjacency_mass_cache_owner_mismatch",
}


def _payloads() -> tuple[dict[str, object], dict[str, object]]:
    generation = "build123d-v49"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    row = {
        ASSEMBLY: {
            "generation": generation, "density_generation": generation, "material_generation": generation,
            "occurrence_generation": generation, "transform_generation": generation, "suppression_generation": generation,
            "mass_generation": generation, "result_generation": generation,
            "material_by_occurrence": {"occurrence:a": "steel"}, "result_material_by_occurrence": {"occurrence:a": "steel"},
            "density_kg_m3_by_occurrence": {"occurrence:a": 7850.0}, "result_density_kg_m3_by_occurrence": {"occurrence:a": 7850.0},
            "occurrence_transforms": {"occurrence:a": [0.0, 0.0, 0.0, 1.0]}, "result_occurrence_transforms": {"occurrence:a": [0.0, 0.0, 0.0, 1.0]},
            "suppressed_occurrences": [], "result_suppressed_occurrences": [], "assembly_mass_kg": 12.5, "result_assembly_mass_kg": 12.5,
            "assembly_owner": "assembly:v49", "result_assembly_owner": "assembly:v49", **result,
        },
        SKETCH: {
            "generation": generation, "constraint_generation": generation, "dof_generation": generation, "plane_generation": generation,
            "unit_generation": generation, "wire_generation": generation, "result_generation": generation,
            "constraints": ["horizontal:e1", "distance:e1=40mm"], "result_constraints": ["horizontal:e1", "distance:e1=40mm"],
            "remaining_dof": 0, "result_remaining_dof": 0, "work_plane": "Plane.XY", "result_work_plane": "Plane.XY",
            "length_unit": "mm", "result_length_unit": "mm", "profile_wires": ["wire:outer"], "result_profile_wires": ["wire:outer"],
            "sketch_owner": "sketch:v49", "result_sketch_owner": "sketch:v49", **result,
        },
    }
    step = {
        "generation": generation, "schema_generation": generation, "assembly_generation": generation, "metadata_generation": generation,
        "unit_generation": generation, "tolerance_generation": generation, "result_generation": generation,
        "ap_schema": "AP242", "result_ap_schema": "AP242", "assembly_structure": {"assembly:root": ["part:a"]},
        "result_assembly_structure": {"assembly:root": ["part:a"]}, "color_map": {"part:a": [0.8, 0.2, 0.1]},
        "result_color_map": {"part:a": [0.8, 0.2, 0.1]}, "layer_map": {"part:a": "structure"}, "result_layer_map": {"part:a": "structure"},
        "length_unit": "mm", "result_length_unit": "mm", "model_tolerance_m": 1e-7, "result_model_tolerance_m": 1e-7,
        "import_owner": "import:v49", "result_import_owner": "import:v49", **result,
    }
    boolean = {
        "generation": generation, "deletion_generation": generation, "selector_generation": generation, "adjacency_generation": generation,
        "mass_generation": generation, "cache_generation": generation, "result_generation": generation,
        "deleted_subshapes": ["face:old"], "result_deleted_subshapes": ["face:old"],
        "selector_results": {"selector:top": ["face:new"]}, "result_selector_results": {"selector:top": ["face:new"]},
        "adjacency_map": {"face:new": ["edge:new"]}, "result_adjacency_map": {"face:new": ["edge:new"]},
        "mass_kg": 3.75, "cached_mass_kg": 3.75, "result_mass_kg": 3.75,
        "cache_shape_sha256": "e" * 64, "result_cache_shape_sha256": "e" * 64,
        "history_owner": "history:v49", "result_history_owner": "history:v49", **result,
    }
    return {"reference": [row], "measured": {"cad": [deepcopy(row)]}}, {"replay_identity": {STEP: step, BOOLEAN: boolean}}


def test_v49_positive_public_and_source_replays_are_accepted() -> None:
    public, source = _payloads()
    assert validate_public_identity(public)["status"] == "ok"
    assert validate_source_identity(source)["status"] == "ok"


def test_v49_public_assembly_and_sketch_mutations_are_rejected() -> None:
    public, _ = _payloads()
    public["reference"][0][ASSEMBLY]["result_assembly_mass_kg"] = 13.0
    public["reference"][0][SKETCH]["result_remaining_dof"] = 2
    result = validate_public_identity(public)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2


def test_v49_source_step_and_boolean_mutations_are_rejected() -> None:
    _, source = _payloads()
    source["replay_identity"][STEP]["result_ap_schema"] = "AP203"
    source["replay_identity"][BOOLEAN]["result_cache_shape_sha256"] = "a" * 64
    result = validate_source_identity(source)
    assert result["status"] == "needs_attention"
    assert len(result["issues"]) == 2

