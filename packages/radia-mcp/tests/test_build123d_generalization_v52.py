from copy import deepcopy

from radia_mcp.build123d.selector_exchange_identity_v52 import BREP, GLTF, SELECTOR, WORKPLANE, validate_public_identity, validate_source_identity


CASE_IDS = {
    "v52_public_selector_query_order_stability_topology_change_label_owner_mismatch",
    "v52_public_workplane_localcoordinate_pendingedge_wireclosure_owner_mismatch",
    "v52_source_tool_brep_occversion_location_tshape_serialization_owner_mismatch",
    "v52_source_tool_gltf_axis_scale_material_instance_scene_owner_mismatch",
}


def _generation(prefix: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": prefix, **{name: prefix for name in names}}


def _payloads():
    selected = ["face:1", "face:2"]; labels = {"face:1": "a", "face:2": "b"}
    selector = {**_generation("sel-v52", ("query_generation", "order_generation", "topology_generation", "label_generation", "owner_generation", "result_generation")), "selector_query": "faces().sort_by(Axis.X)", "result_selector_query": "faces().sort_by(Axis.X)", "selected_topology_order": selected, "result_selected_topology_order": selected, "topology_revision": "topology:v52", "result_topology_revision": "topology:v52", "topology_labels": labels, "result_topology_labels": labels, "shape_owner": "shape:sel-v52", "result_shape_owner": "shape:sel-v52", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64}
    frame = {"origin": [0.0, 0.0, 0.0], "x_dir": [1.0, 0.0, 0.0], "z_dir": [0.0, 0.0, 1.0]}; edges = ["edge:1", "edge:2", "edge:3"]
    workplane = {**_generation("wp-v52", ("coordinate_generation", "edge_generation", "closure_generation", "owner_generation", "result_generation")), "local_frame": frame, "result_local_frame": frame, "pending_edge_order": edges, "result_pending_edge_order": edges, "wire_closed": True, "result_wire_closed": True, "builder_owner": "builder:wp-v52", "result_builder_owner": "builder:wp-v52", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64}
    public = {"reference": [{SELECTOR: selector, WORKPLANE: workplane}], "measured": {}}
    location = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    brep = {**_generation("brep-v52", ("version_generation", "location_generation", "tshape_generation", "serialization_generation", "owner_generation", "result_generation")), "occt_version": "7.9.0", "replayed_occt_version": "7.9.0", "shape_location": location, "replayed_shape_location": location, "tshape_sha256": "c" * 64, "replayed_tshape_sha256": "c" * 64, "serialization_sha256": "d" * 64, "replayed_serialization_sha256": "d" * 64, "shape_owner": "shape:brep-v52", "replayed_shape_owner": "shape:brep-v52", "result_sha256": "e" * 64, "accepted_result_sha256": "e" * 64}
    materials = {"part:a": "material:steel"}; instances = {"instance:a": location}
    gltf = {**_generation("gltf-v52", ("axis_generation", "scale_generation", "material_generation", "instance_generation", "owner_generation", "result_generation")), "axis_convention": "Y_up_right_handed", "replayed_axis_convention": "Y_up_right_handed", "length_scale_to_m": 1.0, "replayed_length_scale_to_m": 1.0, "part_materials": materials, "replayed_part_materials": materials, "instance_transforms": instances, "replayed_instance_transforms": instances, "scene_owner": "scene:gltf-v52", "replayed_scene_owner": "scene:gltf-v52", "result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    return public, {"replay_identity": {BREP: brep, GLTF: gltf}}


def test_v52_positive_public_and_source_replays_are_accepted():
    public, source = _payloads()
    assert validate_public_identity(public)["status"] == "ok"
    assert validate_source_identity(source)["status"] == "ok"


def test_v52_public_mutations_are_rejected():
    public, _ = _payloads(); value = deepcopy(public); value["reference"][0][SELECTOR]["result_topology_revision"] = "topology:stale"
    assert validate_public_identity(value)["status"] == "needs_attention"


def test_v52_source_mutations_are_rejected():
    _, source = _payloads(); value = deepcopy(source); value["replay_identity"][GLTF]["replayed_axis_convention"] = "Z_up_left_handed"
    assert validate_source_identity(value)["status"] == "needs_attention"


def test_v52_invalid_canonical_records_are_rejected():
    public, _ = _payloads(); value = deepcopy(public); value["reference"][0][WORKPLANE]["wire_closed"] = False; value["reference"][0][WORKPLANE]["result_wire_closed"] = False
    assert validate_public_identity(value)["status"] == "needs_attention"
