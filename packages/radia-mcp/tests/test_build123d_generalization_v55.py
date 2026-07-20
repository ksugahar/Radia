from copy import deepcopy
import math

from radia_mcp.build123d.thread_sheet_identity_v55 import (
    BREP,
    PMI,
    SHEET,
    THREAD,
    validate_public_identity,
    validate_source_identity,
)


CASE_IDS = {
    "v55_public_screwthread_pitch_handedness_start_topology_volume_owner_mismatch",
    "v55_public_sheetmetal_thickness_bendallowance_neutralaxis_flatpattern_owner_mismatch",
    "v55_source_tool_step_pmi_tolerance_datum_unit_product_owner_mismatch",
    "v55_source_tool_brep_repair_tolerance_sewnshell_orientation_volume_owner_mismatch",
}


def _payloads():
    generation = "build123d-v55-test"
    generations = lambda names: {name: generation for name in names}
    topology = {"solids": 1, "shells": 1, "faces": 12, "edges": 26, "vertices": 16}
    thread = {"generation": generation, **generations(("pitch_generation", "handedness_generation", "start_generation", "topology_generation", "volume_generation", "owner_generation", "result_generation")), "pitch_m": 2.0e-3, "result_pitch_m": 2.0e-3, "handedness": "right", "result_handedness": "right", "start_count": 1, "result_start_count": 1, "thread_topology": topology, "result_thread_topology": topology, "volume_m3": 1.25e-5, "result_volume_m3": 1.25e-5, "shape_owner": "shape:thread-v55", "result_shape_owner": "shape:thread-v55", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    thickness = 1.0e-3; radius = 5.0e-3; angle = math.pi / 2.0; k = 0.4; allowance = angle * (radius + k * thickness)
    sheet = {"generation": generation, **generations(("thickness_generation", "bend_generation", "neutralaxis_generation", "pattern_generation", "area_generation", "owner_generation", "result_generation")), "thickness_m": thickness, "result_thickness_m": thickness, "inside_bend_radius_m": radius, "result_inside_bend_radius_m": radius, "bend_angle_rad": angle, "result_bend_angle_rad": angle, "neutral_axis_factor": k, "result_neutral_axis_factor": k, "bend_allowance_m": allowance, "result_bend_allowance_m": allowance, "flat_pattern_area_m2": 0.020, "result_flat_pattern_area_m2": 0.020, "folded_surface_area_m2": 0.020, "result_folded_surface_area_m2": 0.020, "shape_owner": "shape:sheet-v55", "result_shape_owner": "shape:sheet-v55", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    tolerances = [{"feature": "feature:hole-1", "kind": "position", "value": 0.05, "datum_refs": ["A", "B", "C"]}]
    pmi = {"generation": generation, **generations(("tolerance_generation", "datum_generation", "unit_generation", "product_generation", "revision_generation", "owner_generation", "result_generation")), "geometric_tolerances": tolerances, "replayed_geometric_tolerances": tolerances, "datum_frame": {"A": "face:base", "B": "face:side", "C": "axis:hole"}, "replayed_datum_frame": {"A": "face:base", "B": "face:side", "C": "axis:hole"}, "length_unit": "mm", "replayed_length_unit": "mm", "product_association": {"feature:hole-1": "part:housing"}, "replayed_product_association": {"feature:hole-1": "part:housing"}, "document_revision": "document:v55-r3", "replayed_document_revision": "document:v55-r3", "document_owner": "document:pmi-v55", "replayed_document_owner": "document:pmi-v55", "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64}
    shells = [{"shell": "shell:1", "faces": ["face:1", "face:2", "face:3", "face:4"], "closed": True}]; orientations = {"face:1": 1, "face:2": 1, "face:3": 1, "face:4": 1}
    brep = {"generation": generation, **generations(("tolerance_generation", "shell_generation", "orientation_generation", "volume_generation", "owner_generation", "result_generation")), "healing_tolerance_m": 1.0e-6, "replayed_healing_tolerance_m": 1.0e-6, "sewn_shells": shells, "replayed_sewn_shells": shells, "face_orientations": orientations, "replayed_face_orientations": orientations, "closed_volume_count": 1, "replayed_closed_volume_count": 1, "volume_m3": 8.0e-5, "replayed_volume_m3": 8.0e-5, "shape_owner": "shape:brep-v55", "replayed_shape_owner": "shape:brep-v55", "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64}
    return {"reference": [{THREAD: thread, SHEET: sheet}], "measured": {}}, {"replay_identity": {PMI: pmi, BREP: brep}}


def test_v55_positive_public_and_source_identities_are_accepted():
    public, source = _payloads()
    assert validate_public_identity(public)["status"] == "ok"
    assert validate_source_identity(source)["status"] == "ok"


def test_v55_frozen_mutations_are_rejected():
    public, source = _payloads(); public = deepcopy(public); source = deepcopy(source)
    public["reference"][0][THREAD]["result_pitch_m"] = 3.0e-3
    public["reference"][0][SHEET]["result_neutral_axis_factor"] = 1.5
    source["replay_identity"][PMI]["replayed_length_unit"] = "inch"
    source["replay_identity"][BREP]["replayed_closed_volume_count"] = 0
    assert validate_public_identity(public)["status"] == "needs_attention"
    assert validate_source_identity(source)["status"] == "needs_attention"


def test_v55_self_consistent_bad_topology_or_bend_allowance_is_rejected():
    public, _ = _payloads(); public = deepcopy(public)
    public["reference"][0][THREAD]["thread_topology"]["edges"] = 23
    public["reference"][0][THREAD]["result_thread_topology"] = public["reference"][0][THREAD]["thread_topology"]
    public["reference"][0][SHEET]["bend_allowance_m"] = public["reference"][0][SHEET]["result_bend_allowance_m"] = 0.1
    assert validate_public_identity(public)["status"] == "needs_attention"


def test_v55_self_consistent_dangling_pmi_or_open_brep_is_rejected():
    _, source = _payloads(); source = deepcopy(source)
    source["replay_identity"][PMI]["product_association"] = source["replay_identity"][PMI]["replayed_product_association"] = {"feature:other": "part:housing"}
    source["replay_identity"][BREP]["sewn_shells"][0]["closed"] = False
    source["replay_identity"][BREP]["replayed_sewn_shells"] = source["replay_identity"][BREP]["sewn_shells"]
    assert validate_source_identity(source)["status"] == "needs_attention"
