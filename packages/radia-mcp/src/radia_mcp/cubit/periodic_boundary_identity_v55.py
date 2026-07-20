"""Headless periodic, boundary-layer, merge, and Exodus identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


PERIODIC = (
    "periodic_hex_nodeequivalence_transform_highorderjacobian_block_owner_identity"
)
BOUNDARY_LAYER = (
    "boundarylayer_hex_thickness_growth_cornercollapse_quality_owner_identity"
)
MERGE = "merge_tolerance_entityprovenance_idremap_group_owner_identity"
EXODUS = "exodus_block_attribute_truth_table_elementorder_owner_identity"
_ELEMENT_ORDERS = {
    "HEX8",
    "HEX20",
    "HEX27",
    "TET4",
    "TET10",
    "PYRAMID5",
    "PYRAMID13",
}


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _number(value: object, *, positive: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    number = float(value)
    return math.isfinite(number) and (not positive or number > 0.0)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get(
        "accepted_result_sha256"
    ) == row.get("result_sha256")


def _matrix4(value: object) -> list[list[float]] | None:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 4
    ):
        return None
    matrix: list[list[float]] = []
    for row in value:
        if (
            not isinstance(row, Sequence)
            or isinstance(row, (str, bytes))
            or len(row) != 4
            or not all(_number(item) for item in row)
        ):
            return None
        matrix.append([float(item) for item in row])
    return matrix


def _rotation_determinant(matrix: Sequence[Sequence[float]]) -> float:
    a, b, c = matrix[0][:3]
    d, e, f = matrix[1][:3]
    g, h, i = matrix[2][:3]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def _periodic_ok(row: Mapping[str, object]) -> bool:
    pairs = row.get("node_equivalence")
    transform = _matrix4(row.get("periodic_transform_4x4"))
    jacobians = row.get("high_order_jacobian_samples")
    pairs_ok = (
        isinstance(pairs, Sequence)
        and not isinstance(pairs, (str, bytes))
        and bool(pairs)
    )
    pair_rows: list[tuple[int, int]] = []
    if pairs_ok:
        for pair in pairs:
            if (
                not isinstance(pair, Sequence)
                or isinstance(pair, (str, bytes))
                or len(pair) != 2
                or not all(
                    isinstance(node, int) and not isinstance(node, bool) and node > 0
                    for node in pair
                )
            ):
                pairs_ok = False
                break
            pair_rows.append((pair[0], pair[1]))
        pairs_ok = pairs_ok and (
            len({pair[0] for pair in pair_rows}) == len(pair_rows)
            and len({pair[1] for pair in pair_rows}) == len(pair_rows)
            and not ({pair[0] for pair in pair_rows} & {pair[1] for pair in pair_rows})
        )
    transform_ok = (
        transform is not None
        and all(
            math.isclose(value, expected, abs_tol=1.0e-12)
            for value, expected in zip(transform[3], [0.0, 0.0, 0.0, 1.0])
        )
        and math.isclose(abs(_rotation_determinant(transform)), 1.0, rel_tol=1.0e-10)
    )
    return (
        _generation(
            row,
            "equivalence_generation",
            "transform_generation",
            "jacobian_generation",
            "block_generation",
            "owner_generation",
            "result_generation",
        )
        and pairs_ok
        and transform_ok
        and isinstance(jacobians, Sequence)
        and not isinstance(jacobians, (str, bytes))
        and len(jacobians) == len(pair_rows)
        and all(_number(value, positive=True) for value in jacobians)
        and all(
            row.get("result_" + name) == row.get(name)
            for name in (
                "node_equivalence",
                "periodic_transform_4x4",
                "high_order_jacobian_samples",
                "element_block",
                "mesh_owner",
            )
        )
        and str(row.get("element_block") or "").startswith("block:")
        and str(row.get("mesh_owner") or "").startswith("headless:")
        and _result(row)
    )


def _boundary_layer_ok(row: Mapping[str, object]) -> bool:
    topology = row.get("corner_topology")
    layer_count = row.get("layer_count")
    collapsed = row.get("collapsed_layer_count")
    topology_ok = (
        isinstance(topology, Sequence)
        and not isinstance(topology, (str, bytes))
        and bool(topology)
        and len(topology) == len(set(topology))
        and all(isinstance(item, str) and item.startswith("corner:") for item in topology)
    )
    return (
        _generation(
            row,
            "thickness_generation",
            "growth_generation",
            "topology_generation",
            "collapse_generation",
            "quality_generation",
            "owner_generation",
            "result_generation",
        )
        and _number(row.get("first_layer_thickness_m"), positive=True)
        and _number(row.get("growth_ratio"), positive=True)
        and float(row["growth_ratio"]) >= 1.0
        and isinstance(layer_count, int)
        and not isinstance(layer_count, bool)
        and layer_count >= 1
        and topology_ok
        and isinstance(collapsed, int)
        and not isinstance(collapsed, bool)
        and collapsed == 0
        and _number(row.get("minimum_scaled_jacobian"), positive=True)
        and all(
            row.get("result_" + name) == row.get(name)
            for name in (
                "first_layer_thickness_m",
                "growth_ratio",
                "layer_count",
                "corner_topology",
                "collapsed_layer_count",
                "minimum_scaled_jacobian",
                "mesh_owner",
            )
        )
        and str(row.get("mesh_owner") or "").startswith("headless:")
        and _result(row)
    )


def _merge_ok(row: Mapping[str, object]) -> bool:
    provenance = row.get("source_entity_provenance")
    remap = row.get("entity_id_remap")
    groups = row.get("group_membership")
    if not all(isinstance(value, Mapping) and value for value in (provenance, remap, groups)):
        return False
    provenance_ok = all(
        isinstance(entity, str)
        and entity.startswith(("vertex:", "curve:", "surface:", "volume:"))
        and isinstance(owner, str)
        and owner
        for entity, owner in provenance.items()
    )
    remap_ok = all(
        isinstance(source, str)
        and source in provenance
        and isinstance(target, str)
        and target in provenance
        and source != target
        for source, target in remap.items()
    )
    surviving = (set(provenance) - set(remap)) | set(remap.values())
    groups_ok = all(
        isinstance(group, str)
        and group.startswith("group:")
        and isinstance(members, Sequence)
        and not isinstance(members, (str, bytes))
        and bool(members)
        and len(members) == len(set(members))
        and set(members).issubset(surviving)
        for group, members in groups.items()
    )
    return (
        _generation(
            row,
            "tolerance_generation",
            "provenance_generation",
            "remap_generation",
            "group_generation",
            "revision_generation",
            "owner_generation",
            "result_generation",
        )
        and _number(row.get("merge_tolerance_m"), positive=True)
        and provenance_ok
        and remap_ok
        and groups_ok
        and all(
            row.get("replayed_" + name) == row.get(name)
            for name in (
                "merge_tolerance_m",
                "source_entity_provenance",
                "entity_id_remap",
                "group_membership",
                "database_revision",
                "merge_owner",
            )
        )
        and str(row.get("database_revision") or "").startswith("database:")
        and str(row.get("merge_owner") or "").startswith("headless:")
        and _result(row)
    )


def _exodus_ok(row: Mapping[str, object]) -> bool:
    attributes = row.get("block_attributes")
    truth_table = row.get("variable_truth_table")
    orders = row.get("element_order")
    if not all(isinstance(value, Mapping) and value for value in (attributes, truth_table, orders)):
        return False
    blocks = list(attributes)
    attributes_ok = all(
        isinstance(block, str)
        and block.startswith("block:")
        and isinstance(values, Mapping)
        and bool(values)
        and all(isinstance(name, str) and name for name in values)
        for block, values in attributes.items()
    )
    truth_ok = all(
        isinstance(variable, str)
        and variable
        and isinstance(values, Sequence)
        and not isinstance(values, (str, bytes))
        and len(values) == len(blocks)
        and all(value in (0, 1, False, True) for value in values)
        for variable, values in truth_table.items()
    )
    orders_ok = set(orders) == set(blocks) and all(
        isinstance(order, str) and order in _ELEMENT_ORDERS for order in orders.values()
    )
    return (
        _generation(
            row,
            "attribute_generation",
            "truth_table_generation",
            "order_generation",
            "qa_generation",
            "owner_generation",
            "result_generation",
        )
        and attributes_ok
        and truth_ok
        and orders_ok
        and all(
            row.get("replayed_" + name) == row.get(name)
            for name in (
                "block_attributes",
                "variable_truth_table",
                "element_order",
                "qa_revision",
                "file_owner",
            )
        )
        and str(row.get("qa_revision") or "").startswith("qa:")
        and str(row.get("file_owner") or "").startswith("headless:")
        and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {
        "policy": policy,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if payload.get(PERIODIC) is not None:
        checks["v55_periodic_hex_equivalence_transform_jacobian_owner"] = (
            isinstance(payload[PERIODIC], Mapping) and _periodic_ok(payload[PERIODIC])
        )
    if payload.get(BOUNDARY_LAYER) is not None:
        checks["v55_boundary_layer_thickness_growth_collapse_quality_owner"] = (
            isinstance(payload[BOUNDARY_LAYER], Mapping)
            and _boundary_layer_ok(payload[BOUNDARY_LAYER])
        )
    return _report("cubit_v55_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if payload.get(MERGE) is not None:
        checks["v55_merge_tolerance_provenance_remap_group_owner"] = (
            isinstance(payload[MERGE], Mapping) and _merge_ok(payload[MERGE])
        )
    if payload.get(EXODUS) is not None:
        checks["v55_exodus_block_truth_table_order_qa_owner"] = (
            isinstance(payload[EXODUS], Mapping) and _exodus_ok(payload[EXODUS])
        )
    return _report("cubit_v55_source_identity_v1", checks) if checks else {}
