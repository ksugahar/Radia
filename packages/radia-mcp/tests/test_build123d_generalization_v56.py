from copy import deepcopy
import math

from radia_mcp.build123d.gear_pipe_exchange_identity_v56 import (
    GEAR,
    MESH,
    PIPE,
    STEP,
    validate_public_identity,
    validate_source_identity,
)


CASE_IDS = {
    "v56_public_involutegear_module_toothcount_pressureangle_pitchdiameter_volume_owner_mismatch",
    "v56_public_pipesweep_pathlength_frame_twist_selfintersection_volume_owner_mismatch",
    "v56_source_tool_step_occurrence_transform_matrix_unit_product_owner_mismatch",
    "v56_source_tool_meshformat_watertight_manifold_unit_signedvolume_owner_mismatch",
}


def _payloads() -> tuple[dict[str, object], dict[str, object]]:
    generation = "build123d-v56-test"
    generations = lambda names: {name: generation for name in names}
    module = 2.0e-3
    teeth = 24
    pressure = math.radians(20.0)
    gear = {
        "generation": generation,
        **generations(("module_generation", "tooth_generation", "pressure_generation", "pitch_generation", "volume_generation", "owner_generation", "result_generation")),
        "module_m": module, "result_module_m": module,
        "tooth_count": teeth, "result_tooth_count": teeth,
        "pressure_angle_rad": pressure, "result_pressure_angle_rad": pressure,
        "pitch_diameter_m": module * teeth, "result_pitch_diameter_m": module * teeth,
        "volume_m3": 1.62e-5, "result_volume_m3": 1.62e-5,
        "shape_owner": "shape:gear-v56", "result_shape_owner": "shape:gear-v56",
        "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64,
    }
    frame = {"method": "parallel_transport", "samples": 65, "closed": False}
    pipe = {
        "generation": generation,
        **generations(("path_generation", "frame_generation", "twist_generation", "intersection_generation", "volume_generation", "owner_generation", "result_generation")),
        "path_length_m": 0.125, "result_path_length_m": 0.125,
        "transport_frame": frame, "result_transport_frame": frame,
        "net_twist_rad": 0.0, "result_net_twist_rad": 0.0,
        "self_intersection": False, "result_self_intersection": False,
        "volume_m3": 2.45e-5, "result_volume_m3": 2.45e-5,
        "shape_owner": "shape:pipe-v56", "result_shape_owner": "shape:pipe-v56",
        "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64,
    }
    transform = [[1.0, 0.0, 0.0, 0.01], [0.0, 1.0, 0.0, 0.02], [0.0, 0.0, 1.0, 0.03], [0.0, 0.0, 0.0, 1.0]]
    occurrence = {"id": "occurrence:gear-1", "product": "product:gear"}
    assembly_frame = {"parent": "assembly:root", "child": "occurrence:gear-1"}
    step = {
        "generation": generation,
        **generations(("transform_generation", "unit_generation", "product_generation", "frame_generation", "owner_generation", "result_generation")),
        "occurrence_transform_4x4": transform, "replayed_occurrence_transform_4x4": transform,
        "length_unit": "m", "replayed_length_unit": "m",
        "product_occurrence": occurrence, "replayed_product_occurrence": occurrence,
        "assembly_frame": assembly_frame, "replayed_assembly_frame": assembly_frame,
        "document_owner": "document:step-v56", "replayed_document_owner": "document:step-v56",
        "result_sha256": "7" * 64, "accepted_result_sha256": "7" * 64,
    }
    mesh = {
        "generation": generation,
        **generations(("format_generation", "watertight_generation", "manifold_generation", "unit_generation", "volume_generation", "owner_generation", "result_generation")),
        "mesh_format": "3mf", "replayed_mesh_format": "3mf",
        "watertight": True, "replayed_watertight": True,
        "manifold": True, "replayed_manifold": True,
        "length_unit": "mm", "replayed_length_unit": "mm",
        "unit_scale_to_m": 1.0e-3, "replayed_unit_scale_to_m": 1.0e-3,
        "signed_volume_m3": 4.8e-5, "replayed_signed_volume_m3": 4.8e-5,
        "mesh_owner": "mesh:3mf-v56", "replayed_mesh_owner": "mesh:3mf-v56",
        "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
    }
    return {"reference": [{GEAR: gear, PIPE: pipe}], "measured": {}}, {"replay_identity": {STEP: step, MESH: mesh}}


def test_v56_positive_public_and_source_identities_are_accepted() -> None:
    public, source = _payloads()
    assert validate_public_identity(public)["status"] == "ok"
    assert validate_source_identity(source)["status"] == "ok"


def test_v56_frozen_mutations_are_rejected() -> None:
    public, source = _payloads()
    public = deepcopy(public)
    source = deepcopy(source)
    public["reference"][0][GEAR]["result_pitch_diameter_m"] = 0.060
    public["reference"][0][PIPE]["result_self_intersection"] = True
    source["replay_identity"][STEP]["replayed_length_unit"] = "inch"
    source["replay_identity"][MESH]["replayed_signed_volume_m3"] = -4.8e-5
    assert validate_public_identity(public)["status"] == "needs_attention"
    assert validate_source_identity(source)["status"] == "needs_attention"


def test_v56_self_consistent_geometry_contradictions_are_rejected() -> None:
    public, _ = _payloads()
    public = deepcopy(public)
    gear = public["reference"][0][GEAR]
    gear["pitch_diameter_m"] = gear["result_pitch_diameter_m"] = 0.060
    pipe = public["reference"][0][PIPE]
    pipe["self_intersection"] = pipe["result_self_intersection"] = True
    assert validate_public_identity(public)["status"] == "needs_attention"


def test_v56_self_consistent_exchange_contradictions_are_rejected() -> None:
    _, source = _payloads()
    source = deepcopy(source)
    step = source["replay_identity"][STEP]
    step["occurrence_transform_4x4"] = step["replayed_occurrence_transform_4x4"] = [[1.0]]
    mesh = source["replay_identity"][MESH]
    mesh["unit_scale_to_m"] = mesh["replayed_unit_scale_to_m"] = 0.0254
    assert validate_source_identity(source)["status"] == "needs_attention"
