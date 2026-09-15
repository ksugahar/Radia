"""Topology and source-replay identity checks for Coreform v49 summaries."""

from __future__ import annotations

import math
from collections.abc import Mapping


_HEX = "hex_sweep_source_target_topology_layer_interval_bias_periodic_owner_identity"
_PYRAMID = "pyramid_transition_apex_orientation_jacobian_boundary_block_owner_identity"
_JOURNAL = "journal_undo_checkpoint_entity_allocator_replay_cursor_revision_owner_identity"
_NETGEN = "netgen_export_order_curved_node_boundary_name_index_digest_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _pairs(value: object) -> bool:
    return (
        isinstance(value, list)
        and bool(value)
        and all(isinstance(pair, list) and len(pair) == 2 and all(isinstance(node, int) and not isinstance(node, bool) and node > 0 for node in pair) for pair in value)
        and len({pair[0] for pair in value}) == len(value)
        and len({pair[1] for pair in value}) == len(value)
    )


def _hex_ok(row: Mapping[str, object]) -> bool:
    topology = row.get("source_target_topology")
    layers = row.get("layer_count")
    bias = row.get("interval_bias")
    pairs = row.get("periodic_node_pairs")
    return (
        _generations(row, "topology_generation", "layer_generation", "interval_generation", "periodic_generation", "result_generation")
        and isinstance(topology, Mapping)
        and set(topology) == {"source", "target", "side_count"}
        and topology.get("source") != topology.get("target")
        and int(topology.get("side_count", 0)) > 0
        and row.get("result_source_target_topology") == topology
        and isinstance(layers, int)
        and not isinstance(layers, bool)
        and layers > 0
        and row.get("result_layer_count") == layers
        and isinstance(bias, (int, float))
        and math.isfinite(float(bias))
        and float(bias) > 0.0
        and row.get("result_interval_bias") == bias
        and _pairs(pairs)
        and row.get("result_periodic_node_pairs") == pairs
        and str(row.get("mesh_owner") or "").startswith("headless:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _result(row)
    )


def _pyramid_ok(row: Mapping[str, object]) -> bool:
    orientations = row.get("face_orientations")
    jacobian = row.get("minimum_scaled_jacobian")
    blocks = row.get("boundary_blocks")
    return (
        _generations(row, "apex_generation", "orientation_generation", "jacobian_generation", "boundary_generation", "result_generation")
        and isinstance(row.get("apex_node_id"), int)
        and int(row["apex_node_id"]) > 0
        and row.get("result_apex_node_id") == row.get("apex_node_id")
        and isinstance(orientations, list)
        and len(orientations) == 4
        and all(value in {-1, 1} for value in orientations)
        and row.get("result_face_orientations") == orientations
        and isinstance(jacobian, (int, float))
        and math.isfinite(float(jacobian))
        and float(jacobian) > 0.0
        and row.get("result_minimum_scaled_jacobian") == jacobian
        and isinstance(blocks, Mapping)
        and bool(blocks)
        and all(
            str(name)
            and isinstance(value, int)
            and not isinstance(value, bool)
            and value > 0
            for name, value in blocks.items()
        )
        and row.get("result_boundary_blocks") == blocks
        and str(row.get("mesh_owner") or "").startswith("headless:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _result(row)
    )


def _journal_ok(row: Mapping[str, object]) -> bool:
    allocator = row.get("entity_id_allocator")
    cursor = row.get("replay_cursor")
    return (
        _generations(row, "checkpoint_generation", "allocator_generation", "cursor_generation", "revision_generation", "result_generation")
        and str(row.get("undo_checkpoint") or "").startswith("checkpoint:")
        and row.get("result_undo_checkpoint") == row.get("undo_checkpoint")
        and isinstance(allocator, Mapping)
        and bool(allocator)
        and all(
            str(name)
            and isinstance(value, int)
            and not isinstance(value, bool)
            and value > 0
            for name, value in allocator.items()
        )
        and row.get("result_entity_id_allocator") == allocator
        and isinstance(cursor, int)
        and not isinstance(cursor, bool)
        and cursor >= 0
        and row.get("result_replay_cursor") == cursor
        and bool(str(row.get("model_revision") or ""))
        and row.get("result_model_revision") == row.get("model_revision")
        and str(row.get("journal_owner") or "").startswith("headless:")
        and row.get("result_journal_owner") == row.get("journal_owner")
        and _result(row)
    )


def _node_map(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(str(entity) and isinstance(nodes, list) and bool(nodes) and all(isinstance(node, int) and not isinstance(node, bool) and node > 0 for node in nodes) for entity, nodes in value.items())
    )


def _netgen_ok(row: Mapping[str, object]) -> bool:
    order = row.get("element_order")
    nodes = row.get("curved_node_map")
    boundaries = row.get("boundary_names")
    return (
        _generations(row, "order_generation", "curved_node_generation", "boundary_generation", "index_generation", "file_generation", "result_generation")
        and isinstance(order, int)
        and not isinstance(order, bool)
        and order >= 2
        and row.get("result_element_order") == order
        and _node_map(nodes)
        and row.get("result_curved_node_map") == nodes
        and isinstance(boundaries, Mapping)
        and bool(boundaries)
        and all(str(identifier) and isinstance(name, str) and name for identifier, name in boundaries.items())
        and len(set(boundaries.values())) == len(boundaries)
        and row.get("result_boundary_names") == boundaries
        and row.get("index_base") == 1
        and row.get("result_index_base") == 1
        and _digest(row.get("export_sha256"))
        and row.get("result_export_sha256") == row.get("export_sha256")
        and str(row.get("export_owner") or "").startswith("headless:")
        and row.get("result_export_owner") == row.get("export_owner")
        and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, ok in checks.items() if not ok]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    hex_sweep = payload.get(_HEX)
    pyramid = payload.get(_PYRAMID)
    if hex_sweep is not None:
        checks["v49_hex_sweep_topology_layers_bias_periodic_owner"] = isinstance(hex_sweep, Mapping) and _hex_ok(hex_sweep)
    if pyramid is not None:
        checks["v49_pyramid_apex_orientation_jacobian_boundary_owner"] = isinstance(pyramid, Mapping) and _pyramid_ok(pyramid)
    return _report("cubit_v49_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    journal = payload.get(_JOURNAL)
    netgen = payload.get(_NETGEN)
    if journal is not None:
        checks["v49_journal_checkpoint_allocator_cursor_revision_owner"] = isinstance(journal, Mapping) and _journal_ok(journal)
    if netgen is not None:
        checks["v49_netgen_order_curved_nodes_boundaries_index_digest_owner"] = isinstance(netgen, Mapping) and _netgen_ok(netgen)
    return _report("cubit_v49_source_identity_v1", checks) if checks else {}
