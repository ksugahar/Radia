"""Headless mesh-topology and transaction replay identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .periodic_boundary_identity_v55 import (
    validate_public_identity as validate_public_v55_identity,
    validate_source_identity as validate_source_v55_identity,
)


HEX = "hex_sweep_edgecorrespondence_highorder_jacobian_block_owner_identity"
TRANSITION = "tethex_pyramid_transition_facediagonal_orientation_quality_owner_identity"
JOURNAL = "journal_undo_transaction_entityallocator_checkpoint_owner_identity"
EXODUS = "exodus_sideset_distributionfactor_topology_qa_owner_identity"
_SIDE_COUNTS = {"HEX8": 6, "TET4": 4, "PYRAMID5": 5}


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _positive_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) and float(value) > 0.0


def _hex_ok(row: Mapping[str, object]) -> bool:
    correspondence = row.get("edge_correspondence")
    jacobians = row.get("high_order_jacobian_samples")
    pairs_ok = isinstance(correspondence, Sequence) and not isinstance(correspondence, (str, bytes)) and bool(correspondence)
    if pairs_ok:
        pairs = [tuple(pair) for pair in correspondence if isinstance(pair, Sequence) and not isinstance(pair, (str, bytes))]
        pairs_ok = (
            len(pairs) == len(correspondence)
            and all(len(pair) == 2 and all(isinstance(node, int) and not isinstance(node, bool) and node > 0 for node in pair) for pair in pairs)
            and len({pair[0] for pair in pairs}) == len(pairs)
            and len({pair[1] for pair in pairs}) == len(pairs)
        )
    return (
        _generation(row, "edge_generation", "jacobian_generation", "block_generation", "owner_generation", "result_generation")
        and pairs_ok
        and isinstance(jacobians, Sequence) and not isinstance(jacobians, (str, bytes)) and len(jacobians) == len(correspondence)
        and all(_positive_number(value) for value in jacobians)
        and row.get("result_edge_correspondence") == correspondence
        and row.get("result_high_order_jacobian_samples") == jacobians
        and str(row.get("element_block") or "").startswith("block:")
        and row.get("result_element_block") == row.get("element_block")
        and str(row.get("mesh_owner") or "").startswith("headless:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _result(row)
    )


def _transition_ok(row: Mapping[str, object]) -> bool:
    transitions = row.get("pyramid_transitions")
    transitions_ok = isinstance(transitions, Sequence) and not isinstance(transitions, (str, bytes)) and bool(transitions)
    if transitions_ok:
        seen: set[int] = set()
        for item in transitions:
            if not isinstance(item, Mapping) or set(item) != {"element", "base_nodes", "face_diagonal", "orientation", "scaled_jacobian"}:
                transitions_ok = False; break
            element, base, diagonal = item["element"], item["base_nodes"], item["face_diagonal"]
            if not (
                isinstance(element, int)
                and not isinstance(element, bool)
                and element > 0
                and element not in seen
                and isinstance(base, Sequence)
                and not isinstance(base, (str, bytes))
                and len(base) == 4
                and all(isinstance(node, int) and not isinstance(node, bool) and node > 0 for node in base)
                and len(set(base)) == 4
                and isinstance(diagonal, Sequence)
                and not isinstance(diagonal, (str, bytes))
                and len(diagonal) == 2
                and all(isinstance(node, int) and not isinstance(node, bool) and node > 0 for node in diagonal)
                and set(diagonal) in ({base[0], base[2]}, {base[1], base[3]})
                and item["orientation"] == 1
                and _positive_number(item["scaled_jacobian"])
            ):
                transitions_ok = False; break
            seen.add(element)
    return (
        _generation(row, "diagonal_generation", "orientation_generation", "quality_generation", "topology_generation", "owner_generation", "result_generation")
        and transitions_ok and row.get("interface_topology") == "tet-pyramid-hex"
        and row.get("result_pyramid_transitions") == transitions
        and row.get("result_interface_topology") == row.get("interface_topology")
        and str(row.get("mesh_owner") or "").startswith("headless:")
        and row.get("result_mesh_owner") == row.get("mesh_owner") and _result(row)
    )


def _journal_ok(row: Mapping[str, object]) -> bool:
    before = row.get("entity_allocator_before")
    after = row.get("entity_allocator_after")
    depth = row.get("undo_depth")
    checkpoint = row.get("command_checkpoint")
    return (
        _generation(row, "transaction_generation", "allocator_generation", "checkpoint_generation", "revision_generation", "owner_generation", "result_generation")
        and str(row.get("transaction_id") or "").startswith("transaction:")
        and isinstance(depth, int) and not isinstance(depth, bool) and depth >= 1
        and isinstance(before, int) and not isinstance(before, bool) and before >= 1 and after == before
        and isinstance(checkpoint, int) and not isinstance(checkpoint, bool) and checkpoint >= 1
        and str(row.get("database_revision") or "").startswith("database:")
        and all(row.get("replayed_" + name) == row.get(name) for name in ("transaction_id", "undo_depth", "entity_allocator_before", "entity_allocator_after", "command_checkpoint", "database_revision"))
        and str(row.get("journal_owner") or "").startswith("headless:")
        and row.get("replayed_journal_owner") == row.get("journal_owner") and _result(row)
    )


def _exodus_ok(row: Mapping[str, object]) -> bool:
    topology = row.get("sideset_topology")
    factors = row.get("distribution_factors")
    qa = row.get("qa_record")
    topology_ok = isinstance(topology, Sequence) and not isinstance(topology, (str, bytes)) and bool(topology)
    if topology_ok:
        seen: set[tuple[int, int]] = set()
        for item in topology:
            if not isinstance(item, Mapping) or set(item) != {"element", "side", "topology"}:
                topology_ok = False; break
            topology_name = item["topology"]
            ids_ok = all(
                isinstance(item[name], int)
                and not isinstance(item[name], bool)
                and item[name] >= 1
                for name in ("element", "side")
            )
            if not (
                isinstance(topology_name, str)
                and topology_name in _SIDE_COUNTS
                and ids_ok
                and item["side"] <= _SIDE_COUNTS[topology_name]
            ):
                topology_ok = False; break
            key = (item["element"], item["side"])
            if key in seen:
                topology_ok = False; break
            seen.add(key)
    qa_ok = isinstance(qa, Mapping) and set(qa) == {"application", "version", "date", "time"} and all(isinstance(value, str) and value for value in qa.values())
    return (
        _generation(row, "sideset_generation", "factor_generation", "topology_generation", "qa_generation", "owner_generation", "result_generation")
        and topology_ok
        and isinstance(factors, Sequence) and not isinstance(factors, (str, bytes)) and len(factors) == len(topology) and all(_positive_number(value) for value in factors)
        and qa_ok and str(row.get("mesh_revision") or "").startswith("mesh:")
        and all(row.get("replayed_" + name) == row.get(name) for name in ("sideset_topology", "distribution_factors", "qa_record", "mesh_revision"))
        and str(row.get("export_owner") or "").startswith("headless:")
        and row.get("replayed_export_owner") == row.get("export_owner") and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, accepted in checks.items() if not accepted]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping): return {}
    checks: dict[str, bool] = {}
    v55 = validate_public_v55_identity(payload)
    if v55: checks.update(v55["checks"])
    if payload.get(HEX) is not None: checks["v54_hex_edge_jacobian_block_owner"] = isinstance(payload[HEX], Mapping) and _hex_ok(payload[HEX])
    if payload.get(TRANSITION) is not None: checks["v54_pyramid_diagonal_orientation_quality_owner"] = isinstance(payload[TRANSITION], Mapping) and _transition_ok(payload[TRANSITION])
    return _report("cubit_v54_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping): return {}
    checks: dict[str, bool] = {}
    v55 = validate_source_v55_identity(payload)
    if v55: checks.update(v55["checks"])
    if payload.get(JOURNAL) is not None: checks["v54_journal_undo_allocator_checkpoint_owner"] = isinstance(payload[JOURNAL], Mapping) and _journal_ok(payload[JOURNAL])
    if payload.get(EXODUS) is not None: checks["v54_exodus_factor_topology_qa_owner"] = isinstance(payload[EXODUS], Mapping) and _exodus_ok(payload[EXODUS])
    return _report("cubit_v54_source_identity_v1", checks) if checks else {}
