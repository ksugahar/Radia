from copy import deepcopy

from radia_mcp.build123d.assembly_tessellation_identity_v54 import (
    ASSEMBLY,
    LOFT,
    STEP,
    TESSELLATION,
    validate_public_identity,
    validate_source_identity,
)


CASE_IDS = {
    "v54_public_assembly_massproperty_density_location_center_inertia_owner_mismatch",
    "v54_public_loft_section_orientation_parameter_seam_topology_owner_mismatch",
    "v54_source_tool_step_ap242_unit_productstructure_color_layer_owner_mismatch",
    "v54_source_tool_tessellation_deflection_angle_orientation_index_owner_mismatch",
}


def _generations(generation: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{name: generation for name in names}}


def _payloads():
    locations = {"solid:a": {"translation_m": [0.0, 0.0, 0.0], "quaternion_wxyz": [1.0, 0.0, 0.0, 0.0]}}
    densities = {"solid:a": 7850.0}
    tensor = [[0.02, 0.001, 0.0], [0.001, 0.03, 0.0], [0.0, 0.0, 0.04]]
    assembly = {**_generations("assembly-v54", ("density_generation", "location_generation", "center_generation", "inertia_generation", "owner_generation", "result_generation")), "solid_densities_kg_m3": densities, "result_solid_densities_kg_m3": densities, "located_solids": locations, "result_located_solids": locations, "center_of_mass_m": [0.0, 0.0, 0.0], "result_center_of_mass_m": [0.0, 0.0, 0.0], "inertia_tensor_kg_m2": tensor, "result_inertia_tensor_kg_m2": tensor, "assembly_owner": "assembly:v54", "result_assembly_owner": "assembly:v54", "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64}
    sections = [{"section": "wire:1", "orientation": 1, "parameter_start": 0.0, "seam": "vertex:1"}, {"section": "wire:2", "orientation": 1, "parameter_start": 0.0, "seam": "vertex:2"}]
    topology = {"solids": 1, "shells": 1, "faces": 6, "edges": 12, "vertices": 8}
    loft = {**_generations("loft-v54", ("section_generation", "parameter_generation", "seam_generation", "topology_generation", "owner_generation", "result_generation")), "section_correspondence": sections, "result_section_correspondence": sections, "resulting_topology": topology, "result_resulting_topology": topology, "shape_owner": "shape:loft-v54", "result_shape_owner": "shape:loft-v54", "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64}
    structure = {"product:root": ["part:a"], "part:a": ["solid:a"]}; colors = {"part:a": [0.2, 0.3, 0.4]}; layers = {"part:a": "layer:main"}
    step = {**_generations("step-v54", ("unit_generation", "structure_generation", "color_generation", "layer_generation", "owner_generation", "result_generation")), "schema": "AP242", "replayed_schema": "AP242", "length_unit": "mm", "replayed_length_unit": "mm", "product_structure": structure, "replayed_product_structure": structure, "component_colors": colors, "replayed_component_colors": colors, "component_layers": layers, "replayed_component_layers": layers, "document_owner": "document:v54", "replayed_document_owner": "document:v54", "result_sha256": "7" * 64, "accepted_result_sha256": "7" * 64}
    vertices = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]; triangles = [[0, 1, 2]]; orientations = [1]
    tessellation = {**_generations("tess-v54", ("deflection_generation", "angle_generation", "orientation_generation", "index_generation", "owner_generation", "result_generation")), "linear_deflection_m": 1.0e-4, "replayed_linear_deflection_m": 1.0e-4, "angular_deflection_rad": 0.2, "replayed_angular_deflection_rad": 0.2, "vertices_m": vertices, "replayed_vertices_m": vertices, "triangle_indices": triangles, "replayed_triangle_indices": triangles, "triangle_orientations": orientations, "replayed_triangle_orientations": orientations, "shape_owner": "shape:tess-v54", "replayed_shape_owner": "shape:tess-v54", "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64}
    return {"reference": [{ASSEMBLY: assembly, LOFT: loft}], "measured": {}}, {"replay_identity": {STEP: step, TESSELLATION: tessellation}}


def test_v54_positive_public_and_source_identities_are_accepted():
    public, source = _payloads()
    assert validate_public_identity(public)["status"] == "ok"
    assert validate_source_identity(source)["status"] == "ok"


def test_v54_frozen_mutations_are_rejected():
    public, source = _payloads(); public = deepcopy(public); source = deepcopy(source)
    public["reference"][0][ASSEMBLY]["result_center_of_mass_m"] = [1.0, 0.0, 0.0]
    public["reference"][0][LOFT]["result_section_correspondence"] = []
    source["replay_identity"][STEP]["replayed_schema"] = "AP203"
    source["replay_identity"][TESSELLATION]["replayed_triangle_indices"] = [[0, 1, 9]]
    assert validate_public_identity(public)["status"] == "needs_attention"
    assert validate_source_identity(source)["status"] == "needs_attention"


def test_v54_self_consistent_nonphysical_records_are_rejected():
    public, source = _payloads(); public = deepcopy(public); source = deepcopy(source)
    public["reference"][0][ASSEMBLY]["located_solids"]["solid:a"]["quaternion_wxyz"] = [2.0, 0.0, 0.0, 0.0]
    public["reference"][0][ASSEMBLY]["result_located_solids"] = public["reference"][0][ASSEMBLY]["located_solids"]
    public["reference"][0][LOFT]["resulting_topology"]["edges"] = 11
    public["reference"][0][LOFT]["result_resulting_topology"] = public["reference"][0][LOFT]["resulting_topology"]
    source["replay_identity"][STEP]["component_colors"]["part:a"] = [1.2, 0.0, 0.0]
    source["replay_identity"][STEP]["replayed_component_colors"] = source["replay_identity"][STEP]["component_colors"]
    source["replay_identity"][TESSELLATION]["triangle_indices"] = source["replay_identity"][TESSELLATION]["replayed_triangle_indices"] = [[0, 1, 9]]
    assert validate_public_identity(public)["status"] == "needs_attention"
    assert validate_source_identity(source)["status"] == "needs_attention"


def test_v54_malformed_nested_values_reject_without_raising():
    public, source = _payloads(); public = deepcopy(public); source = deepcopy(source)
    public["reference"][0][ASSEMBLY]["solid_densities_kg_m3"] = {"solid:a": [7850.0]}
    source["replay_identity"][TESSELLATION]["triangle_indices"] = [[[0], 1, 2]]
    source["replay_identity"][TESSELLATION]["triangle_orientations"] = [[1]]
    assert validate_public_identity(public)["status"] == "needs_attention"
    assert validate_source_identity(source)["status"] == "needs_attention"
