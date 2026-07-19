from __future__ import annotations

from copy import deepcopy

from radia_mcp.build123d.cross_artifact_cad_lineage_v47 import ASSEMBLY, COMPOUND, EXTERNAL, ROUNDTRIP, validate_public_identity, validate_source_identity


PROMOTED_CASE_IDS = {
    "v47_public_compound_child_permutation_mass_property_aggregate_mismatch",
    "v47_public_assembly_mate_hierarchy_transform_owner_chain_mismatch",
    "v47_source_tool_step_brep_label_uuid_roundtrip_duplicate_owner_mismatch",
    "v47_source_tool_sketch_external_reference_dependency_cycle_revision_mismatch",
}


def _payloads() -> tuple[dict[str, object], dict[str, object]]:
    compound_generation = "compound-v47"
    assembly_generation = "assembly-v47"
    children = ["body:a", "body:b", "body:c"]
    masses = [1.0, 2.0, 3.0]
    hierarchy = ["root", "root/arm", "root/arm/tool"]
    owners = {name: "assembly:a1" for name in hierarchy}
    row = {
        COMPOUND: {"generation":compound_generation,"child_generation":compound_generation,"mass_property_generation":compound_generation,"result_generation":compound_generation,"child_ids":children,"result_child_ids":children,"child_masses_kg":masses,"result_child_masses_kg":masses,"aggregate_mass_kg":6.0,"result_aggregate_mass_kg":6.0,"owner":"compound:a1","accepted_owner":"compound:a1","result_sha256":"1"*64,"accepted_result_sha256":"1"*64},
        ASSEMBLY: {"generation":assembly_generation,"mate_generation":assembly_generation,"hierarchy_generation":assembly_generation,"transform_generation":assembly_generation,"result_generation":assembly_generation,"mate_hierarchy":hierarchy,"result_mate_hierarchy":hierarchy,"transform_owner_map":owners,"result_transform_owner_map":owners,"local_to_world_transform_sha256":"2"*64,"result_local_to_world_transform_sha256":"2"*64,"owner":"assembly:a1","accepted_owner":"assembly:a1","result_sha256":"3"*64,"accepted_result_sha256":"3"*64},
    }
    public = {"reference":[row],"measured":{"cad":[deepcopy(row)]}}
    roundtrip_generation = "roundtrip-v47"
    external_generation = "external-v47"
    entities = [{"label":"body","uuid":"uuid-body","owner":"body:1"},{"label":"face-top","uuid":"uuid-top","owner":"body:1"}]
    references = ["edge:10@rev7", "vertex:2@rev7"]
    edges = [["sketch:1","edge:10"],["edge:10","geometry:rev7"]]
    source = {"replay_identity":{
        ROUNDTRIP:{"generation":roundtrip_generation,"roundtrip_generation":roundtrip_generation,"result_generation":roundtrip_generation,"entity_identities":entities,"result_entity_identities":entities,"duplicate_label_count":0,"result_duplicate_label_count":0,"duplicate_uuid_count":0,"result_duplicate_uuid_count":0,"result_sha256":"4"*64,"accepted_result_sha256":"4"*64},
        EXTERNAL:{"generation":external_generation,"reference_generation":external_generation,"dependency_generation":external_generation,"result_generation":external_generation,"geometry_revision":"rev7","result_geometry_revision":"rev7","external_references":references,"result_external_references":references,"dependency_edges":edges,"result_dependency_edges":edges,"dependency_cycle_count":0,"result_dependency_cycle_count":0,"result_sha256":"5"*64,"accepted_result_sha256":"5"*64},
    }}
    return public, source


def test_v47_positive_replays_are_accepted() -> None:
    public, source = _payloads()
    assert validate_public_identity(public)["status"] == "ok"
    assert validate_source_identity(source)["status"] == "ok"


def test_v47_compound_and_assembly_mutations_are_rejected() -> None:
    public, _ = _payloads()
    public["reference"][0][COMPOUND]["result_child_ids"] = ["body:c", "body:a", "body:b"]
    public["reference"][0][ASSEMBLY]["result_local_to_world_transform_sha256"] = "a" * 64
    assert validate_public_identity(public)["status"] == "needs_attention"


def test_v47_roundtrip_and_external_reference_mutations_are_rejected() -> None:
    _, source = _payloads()
    source["replay_identity"][ROUNDTRIP]["result_duplicate_uuid_count"] = 1
    source["replay_identity"][EXTERNAL]["result_dependency_cycle_count"] = 1
    assert validate_source_identity(source)["status"] == "needs_attention"
