from copy import deepcopy

from radia_mcp.build123d.assembly_exchange_identity_v53 import HISTORY, MASS, MATE, STEP, validate_public_identity, validate_source_identity


CASE_IDS = {
    "v53_public_assembly_mate_frame_handedness_axis_offset_component_owner_mismatch",
    "v53_public_step_import_unit_tolerance_color_component_hierarchy_owner_mismatch",
    "v53_source_tool_boolean_history_face_ancestry_fillet_chamfer_shape_owner_mismatch",
    "v53_source_tool_massproperty_inertia_tensor_frame_centroid_density_owner_mismatch",
}


def _generations(generation: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{name: generation for name in names}}


def _payloads():
    frame = {"origin": [0.0, 0.0, 0.0], "x_dir": [1.0, 0.0, 0.0], "z_dir": [0.0, 0.0, 1.0]}
    pair = ["component:a", "component:b"]
    mate = {**_generations("mate-v53", ("frame_generation", "handedness_generation", "axis_generation", "offset_generation", "component_generation", "owner_generation", "result_generation")), "mate_frame": frame, "result_mate_frame": frame, "handedness": "right", "result_handedness": "right", "mate_axis": [0.0, 0.0, 1.0], "result_mate_axis": [0.0, 0.0, 1.0], "offset_m": -0.001, "result_offset_m": -0.001, "component_pair": pair, "result_component_pair": pair, "assembly_owner": "assembly:v53", "result_assembly_owner": "assembly:v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    colors = {"component:a": [0.1, 0.2, 0.3], "component:b": [0.7, 0.8, 0.9]}; hierarchy = {"assembly:root": pair}
    step = {**_generations("step-v53", ("unit_generation", "tolerance_generation", "color_generation", "hierarchy_generation", "owner_generation", "result_generation")), "length_unit": "mm", "result_length_unit": "mm", "linear_tolerance_m": 1.0e-6, "result_linear_tolerance_m": 1.0e-6, "component_colors": colors, "result_component_colors": colors, "component_hierarchy": hierarchy, "result_component_hierarchy": hierarchy, "document_owner": "document:v53", "result_document_owner": "document:v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    ancestry = {"face:r1": ["face:a1", "face:b1"]}
    history = {**_generations("history-v53", ("boolean_generation", "ancestry_generation", "fillet_generation", "chamfer_generation", "owner_generation", "result_generation")), "boolean_operation": "cut", "replayed_boolean_operation": "cut", "face_ancestry": ancestry, "replayed_face_ancestry": ancestry, "fillet_edges": ["edge:1"], "replayed_fillet_edges": ["edge:1"], "chamfer_edges": ["edge:2"], "replayed_chamfer_edges": ["edge:2"], "shape_owner": "shape:v53", "replayed_shape_owner": "shape:v53", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64}
    tensor = [[0.02, 0.001, 0.0], [0.001, 0.03, 0.0], [0.0, 0.0, 0.04]]
    mass = {**_generations("mass-v53", ("tensor_generation", "frame_generation", "centroid_generation", "density_generation", "owner_generation", "result_generation")), "inertia_tensor_kg_m2": tensor, "replayed_inertia_tensor_kg_m2": tensor, "reference_frame": frame, "replayed_reference_frame": frame, "centroid_m": [0.1, 0.0, 0.0], "replayed_centroid_m": [0.1, 0.0, 0.0], "density_kg_m3": 7850.0, "replayed_density_kg_m3": 7850.0, "solid_owner": "solid:v53", "replayed_solid_owner": "solid:v53", "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64}
    return {"reference": [{MATE: mate, STEP: step}], "measured": {}}, {"replay_identity": {HISTORY: history, MASS: mass}}


def test_v53_positive_public_and_source_replays_are_accepted():
    public, source = _payloads()
    assert validate_public_identity(public)["status"] == "ok"
    assert validate_source_identity(source)["status"] == "ok"


def test_v53_frozen_mutations_are_rejected():
    public, source = _payloads(); public = deepcopy(public); source = deepcopy(source)
    public["reference"][0][MATE]["result_handedness"] = "left"
    public["reference"][0][STEP]["result_length_unit"] = "inch"
    source["replay_identity"][HISTORY]["replayed_face_ancestry"] = {"face:r1": ["face:stale"]}
    source["replay_identity"][MASS]["replayed_density_kg_m3"] = 2700.0
    assert validate_public_identity(public)["status"] == "needs_attention"
    assert validate_source_identity(source)["status"] == "needs_attention"


def test_v53_self_consistent_nonphysical_records_are_rejected():
    public, source = _payloads(); public = deepcopy(public); source = deepcopy(source)
    bad_colors = {"component:a": [1.2, 0.0, 0.0], "component:b": [0.0, 0.0, 0.0]}
    public["reference"][0][STEP]["component_colors"] = public["reference"][0][STEP]["result_component_colors"] = bad_colors
    bad_tensor = [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.05]]
    source["replay_identity"][MASS]["inertia_tensor_kg_m2"] = source["replay_identity"][MASS]["replayed_inertia_tensor_kg_m2"] = bad_tensor
    assert validate_public_identity(public)["status"] == "needs_attention"
    assert validate_source_identity(source)["status"] == "needs_attention"
