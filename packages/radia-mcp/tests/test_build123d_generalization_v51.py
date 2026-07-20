from copy import deepcopy

from radia_mcp.build123d.mass_sweep_exchange_identity_v51 import validate_public_identity, validate_source_identity


PROMOTED_CASE_IDS = {
    "v51_public_mass_properties_frame_inertia_tensor_parallel_axis_density_shape_owner_mismatch",
    "v51_public_sweep_profile_path_trihedron_transition_selfintersection_history_owner_mismatch",
    "v51_source_tool_step_external_reference_occurrence_name_schema_unit_color_owner_mismatch",
    "v51_source_tool_brep_occt_version_location_precision_triangulation_cache_owner_mismatch",
}


def _records() -> tuple[dict[str, object], dict[str, object]]:
    mass = "mass-v51"
    sweep = "sweep-v51"
    step = "step-v51"
    brep = "brep-v51"
    inertia = [[0.10, 0.0, 0.0], [0.0, 0.20, 0.0], [0.0, 0.0, 0.25]]
    shifted = [[0.10, 0.0, 0.0], [0.0, 0.22, 0.0], [0.0, 0.0, 0.27]]
    history = {"profile:1": ["face:11"], "path:1": ["edge:21"], "result": ["solid:31"]}
    public_row = {
        "mass_properties_frame_inertia_parallel_axis_density_shape_owner_identity": {"generation": mass, **{name: mass for name in ("frame_generation", "inertia_generation", "shift_generation", "density_generation", "owner_generation", "result_generation")}, "coordinate_frame": "global_cartesian", "result_coordinate_frame": "global_cartesian", "mass_kg": 2.0, "result_mass_kg": 2.0, "centroidal_inertia_kg_m2": inertia, "result_centroidal_inertia_kg_m2": inertia, "parallel_axis_shift_m": [0.1, 0.0, 0.0], "result_parallel_axis_shift_m": [0.1, 0.0, 0.0], "shifted_inertia_kg_m2": shifted, "result_shifted_inertia_kg_m2": shifted, "density_kg_m3": 7800.0, "result_density_kg_m3": 7800.0, "shape_owner": "shape:mass-v51", "result_shape_owner": "shape:mass-v51", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64},
        "sweep_profile_path_trihedron_transition_selfintersection_history_owner_identity": {"generation": sweep, **{name: sweep for name in ("profile_generation", "path_generation", "trihedron_generation", "transition_generation", "intersection_generation", "history_generation", "owner_generation", "result_generation")}, "profile_id": "profile:1", "result_profile_id": "profile:1", "path_id": "path:1", "result_path_id": "path:1", "trihedron": "corrected_frenet", "result_trihedron": "corrected_frenet", "transition": "transformed", "result_transition": "transformed", "self_intersections": [], "result_self_intersections": [], "topology_history": history, "result_topology_history": history, "shape_owner": "shape:sweep-v51", "result_shape_owner": "shape:sweep-v51", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64},
    }
    references = {"document:main": ["part:shaft", "part:housing"]}
    names = {"occurrence:1": "shaft", "occurrence:2": "housing"}
    colors = {"occurrence:1": [0.8, 0.8, 0.8], "occurrence:2": [0.2, 0.4, 0.8]}
    location = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    source = {"replay_identity": {
        "step_external_reference_occurrence_name_schema_unit_color_owner_identity": {"generation": step, **{name: step for name in ("reference_generation", "name_generation", "schema_generation", "unit_generation", "color_generation", "owner_generation", "result_generation")}, "external_references": references, "result_external_references": references, "occurrence_names": names, "result_occurrence_names": names, "step_schema": "AP242", "result_step_schema": "AP242", "length_unit": "m", "result_length_unit": "m", "occurrence_colors": colors, "result_occurrence_colors": colors, "document_owner": "document:step-v51", "result_document_owner": "document:step-v51", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64},
        "brep_occt_version_location_precision_triangulation_cache_owner_identity": {"generation": brep, **{name: brep for name in ("version_generation", "location_generation", "precision_generation", "triangulation_generation", "owner_generation", "result_generation")}, "occt_version": "7.9.0", "result_occt_version": "7.9.0", "shape_location": location, "result_shape_location": location, "model_precision_m": 1e-7, "result_model_precision_m": 1e-7, "triangulation_cache_sha256": "4" * 64, "result_triangulation_cache_sha256": "4" * 64, "shape_owner": "shape:brep-v51", "result_shape_owner": "shape:brep-v51", "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64},
    }}
    return {"reference": [public_row], "measured": {"cad": [deepcopy(public_row)]}}, source


def test_v51_positive_public_and_source_replays_are_accepted() -> None:
    public, source = _records()
    assert validate_public_identity(public)["status"] == "ok"
    assert validate_source_identity(source)["status"] == "ok"


def test_v51_public_mutations_are_rejected() -> None:
    public, _ = _records()
    public["reference"][0]["mass_properties_frame_inertia_parallel_axis_density_shape_owner_identity"]["result_density_kg_m3"] = 2700.0
    public["reference"][0]["sweep_profile_path_trihedron_transition_selfintersection_history_owner_identity"]["result_self_intersections"] = ["edge:99"]
    assert validate_public_identity(public)["status"] == "needs_attention"


def test_v51_source_mutations_are_rejected() -> None:
    _, source = _records()
    source["replay_identity"]["step_external_reference_occurrence_name_schema_unit_color_owner_identity"]["result_step_schema"] = "AP203"
    source["replay_identity"]["brep_occt_version_location_precision_triangulation_cache_owner_identity"]["result_occt_version"] = "7.8.0"
    assert validate_source_identity(source)["status"] == "needs_attention"


def test_v51_invalid_canonical_records_are_rejected() -> None:
    public, source = _records()
    public["reference"][0]["mass_properties_frame_inertia_parallel_axis_density_shape_owner_identity"]["shifted_inertia_kg_m2"] = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.25]]
    source["replay_identity"]["step_external_reference_occurrence_name_schema_unit_color_owner_identity"]["occurrence_colors"] = {"occurrence:1": [2.0, 0.0, 0.0]}
    assert validate_public_identity(public)["status"] == "needs_attention"
    assert validate_source_identity(source)["status"] == "needs_attention"
