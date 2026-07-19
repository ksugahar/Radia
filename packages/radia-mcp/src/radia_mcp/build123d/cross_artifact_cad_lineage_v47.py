"""Cross-artifact CAD ownership and dependency checks for v47."""

from __future__ import annotations

import math
from collections.abc import Mapping


COMPOUND = "compound_child_permutation_mass_property_aggregate_identity"
ASSEMBLY = "assembly_mate_hierarchy_transform_owner_chain_identity"
ROUNDTRIP = "step_brep_label_uuid_roundtrip_duplicate_owner_identity"
EXTERNAL = "sketch_external_reference_dependency_cycle_revision_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _result_ok(row: Mapping[str, object], *, owner: bool = False) -> bool:
    return (
        (not owner or bool(str(row.get("owner") or "")) and row.get("accepted_owner") == row.get("owner"))
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _compound_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation") or "")
    children = row.get("child_ids")
    masses = row.get("child_masses_kg")
    aggregate = row.get("aggregate_mass_kg")
    return (
        bool(generation)
        and row.get("child_generation") == generation
        and row.get("mass_property_generation") == generation
        and row.get("result_generation") == generation
        and isinstance(children, list)
        and bool(children)
        and len(set(children)) == len(children)
        and row.get("result_child_ids") == children
        and isinstance(masses, list)
        and len(masses) == len(children)
        and all(isinstance(value, (int, float)) and math.isfinite(float(value)) and float(value) >= 0.0 for value in masses)
        and row.get("result_child_masses_kg") == masses
        and isinstance(aggregate, (int, float))
        and math.isclose(float(aggregate), sum(float(value) for value in masses), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and row.get("result_aggregate_mass_kg") == aggregate
        and _result_ok(row, owner=True)
    )


def _assembly_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation") or "")
    hierarchy = row.get("mate_hierarchy")
    owners = row.get("transform_owner_map")
    digest = row.get("local_to_world_transform_sha256")
    return (
        bool(generation)
        and row.get("mate_generation") == generation
        and row.get("hierarchy_generation") == generation
        and row.get("transform_generation") == generation
        and row.get("result_generation") == generation
        and isinstance(hierarchy, list)
        and bool(hierarchy)
        and len(set(hierarchy)) == len(hierarchy)
        and row.get("result_mate_hierarchy") == hierarchy
        and isinstance(owners, Mapping)
        and set(owners) == set(hierarchy)
        and all(value == row.get("owner") for value in owners.values())
        and row.get("result_transform_owner_map") == owners
        and _digest(digest)
        and row.get("result_local_to_world_transform_sha256") == digest
        and _result_ok(row, owner=True)
    )


def _roundtrip_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation") or "")
    entities = row.get("entity_identities")
    valid_entities = (
        isinstance(entities, list)
        and bool(entities)
        and all(isinstance(item, Mapping) and all(isinstance(item.get(key), str) and bool(item.get(key).strip()) for key in ("label", "uuid", "owner")) for item in entities)
    )
    labels = [item["label"] for item in entities] if valid_entities else []
    uuids = [item["uuid"] for item in entities] if valid_entities else []
    return (
        bool(generation)
        and row.get("roundtrip_generation") == generation
        and row.get("result_generation") == generation
        and valid_entities
        and len(set(labels)) == len(labels)
        and len(set(uuids)) == len(uuids)
        and row.get("result_entity_identities") == entities
        and row.get("duplicate_label_count") == row.get("result_duplicate_label_count") == 0
        and row.get("duplicate_uuid_count") == row.get("result_duplicate_uuid_count") == 0
        and _result_ok(row)
    )


def _acyclic(edges: object) -> bool:
    if not isinstance(edges, list):
        return False
    graph: dict[str, list[str]] = {}
    for edge in edges:
        if not isinstance(edge, list) or len(edge) != 2 or not all(isinstance(node, str) and node for node in edge):
            return False
        graph.setdefault(edge[0], []).append(edge[1])
    active: set[str] = set()
    done: set[str] = set()
    def visit(node: str) -> bool:
        if node in active:
            return False
        if node in done:
            return True
        active.add(node)
        if not all(visit(child) for child in graph.get(node, [])):
            return False
        active.remove(node)
        done.add(node)
        return True
    return all(visit(node) for node in tuple(graph))


def _external_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation") or "")
    revision = str(row.get("geometry_revision") or "")
    references = row.get("external_references")
    edges = row.get("dependency_edges")
    return (
        bool(generation)
        and row.get("reference_generation") == generation
        and row.get("dependency_generation") == generation
        and row.get("result_generation") == generation
        and bool(revision)
        and row.get("result_geometry_revision") == revision
        and isinstance(references, list)
        and bool(references)
        and len(set(references)) == len(references)
        and all(isinstance(reference, str) and reference.endswith(f"@{revision}") for reference in references)
        and row.get("result_external_references") == references
        and _acyclic(edges)
        and row.get("result_dependency_edges") == edges
        and row.get("dependency_cycle_count") == row.get("result_dependency_cycle_count") == 0
        and _result_ok(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, ok in checks.items() if not ok]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    rows: list[Mapping[str, object]] = []
    if isinstance(payload.get("reference"), list):
        rows.extend(row for row in payload["reference"] if isinstance(row, Mapping))
    measured = payload.get("measured")
    if isinstance(measured, Mapping):
        for values in measured.values():
            if isinstance(values, list):
                rows.extend(row for row in values if isinstance(row, Mapping))
    checks: dict[str, bool] = {}
    compounds = [row.get(COMPOUND) for row in rows if COMPOUND in row]
    assemblies = [row.get(ASSEMBLY) for row in rows if ASSEMBLY in row]
    if compounds:
        checks["v47_compound_child_mass_aggregation"] = len(compounds) == len(rows) and all(isinstance(item, Mapping) and _compound_ok(item) for item in compounds)
    if assemblies:
        checks["v47_assembly_mate_transform_owner_chain"] = len(assemblies) == len(rows) and all(isinstance(item, Mapping) and _assembly_ok(item) for item in assemblies)
    return _report("build123d_v47_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping) or not isinstance(payload.get("replay_identity"), Mapping):
        return {}
    identity = payload["replay_identity"]
    checks: dict[str, bool] = {}
    roundtrip = identity.get(ROUNDTRIP)
    external = identity.get(EXTERNAL)
    if roundtrip is not None:
        checks["v47_step_brep_label_uuid_roundtrip"] = isinstance(roundtrip, Mapping) and _roundtrip_ok(roundtrip)
    if external is not None:
        checks["v47_sketch_external_reference_dependency"] = isinstance(external, Mapping) and _external_ok(external)
    return _report("build123d_v47_source_identity_v1", checks) if checks else {}
