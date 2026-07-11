"""Readable motor-loss handoff gate for lumped and 3D thermal models."""
from __future__ import annotations

import json
from collections import Counter, deque
from typing import Any


HEX_CELL_TYPES = {"hex", "hex8", "hex20", "hex27", "hexahedron"}


def _load_json(value: str, expected: type, name: str) -> Any:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{name} must be valid JSON: {exc.msg}") from exc
    if not isinstance(parsed, expected):
        raise ValueError(f"{name} must decode to {expected.__name__}")
    return parsed


def _close(a: float, b: float, relative_tolerance: float) -> bool:
    scale = max(abs(a), abs(b), 1.0)
    return abs(a - b) <= relative_tolerance * scale


def _reachable_nodes(nodes: set[str], branches: list[dict], ambient: str) -> set[str]:
    graph = {node: set() for node in nodes}
    for branch in branches:
        a = str(branch.get("from") or "")
        b = str(branch.get("to") or "")
        if a in graph and b in graph:
            graph[a].add(b)
            graph[b].add(a)
    seen = {ambient} if ambient in graph else set()
    queue = deque(seen)
    while queue:
        node = queue.popleft()
        for neighbour in graph[node] - seen:
            seen.add(neighbour)
            queue.append(neighbour)
    return seen


def evaluate_motor_thermal_handoff(
    loss_buckets: dict,
    network: dict,
    mesh_regions: list,
    relative_tolerance: float = 1.0e-9,
) -> dict:
    """Check one loss table against an LPTN and an all-hex thermal mesh.

    ``loss_buckets`` maps physical region names to non-negative watts.
    Network nodes map those regions with ``source_regions`` and connect to the
    fixed ``ambient_node`` through positive ``resistance_K_per_W`` branches.
    Mesh rows map the same regions to positive hexahedral cell counts and
    ``loss_W``.  The gate checks topology, one-to-one region ownership, and
    conservation of total heat without pretending to solve either model.
    """
    errors: list[str] = []
    if relative_tolerance <= 0:
        raise ValueError("relative_tolerance must be positive")

    losses: dict[str, float] = {}
    for raw_name, raw_value in loss_buckets.items():
        name = str(raw_name).strip()
        if not name:
            errors.append("loss bucket names must be non-empty")
            continue
        try:
            value = float(raw_value)
        except (TypeError, ValueError):
            errors.append(f"loss bucket {name!r} is not numeric")
            continue
        if value < 0:
            errors.append(f"loss bucket {name!r} must be non-negative")
        losses[name] = value
    if not losses:
        errors.append("at least one loss bucket is required")

    nodes = network.get("nodes")
    branches = network.get("branches")
    ambient = str(network.get("ambient_node") or "").strip()
    if not isinstance(nodes, list) or not nodes:
        errors.append("network.nodes must be a non-empty list")
        nodes = []
    if not isinstance(branches, list) or not branches:
        errors.append("network.branches must be a non-empty list")
        branches = []

    node_ids = [str(row.get("id") or "").strip() for row in nodes if isinstance(row, dict)]
    node_set = {item for item in node_ids if item}
    if len(node_ids) != len(nodes) or len(node_set) != len(node_ids):
        errors.append("network node ids must be present and unique")
    if not ambient or ambient not in node_set:
        errors.append("network.ambient_node must name one network node")

    network_regions: list[str] = []
    for row in nodes:
        if not isinstance(row, dict):
            continue
        node_id = str(row.get("id") or "").strip()
        regions = row.get("source_regions", [])
        if not isinstance(regions, list):
            errors.append(f"node {node_id!r} source_regions must be a list")
            continue
        network_regions.extend(str(region).strip() for region in regions if str(region).strip())
        if node_id != ambient:
            try:
                capacitance = float(row.get("capacitance_J_per_K"))
            except (TypeError, ValueError):
                capacitance = -1.0
            if capacitance <= 0:
                errors.append(f"node {node_id!r} needs positive capacitance_J_per_K")

    for index, branch in enumerate(branches):
        if not isinstance(branch, dict):
            errors.append(f"network branch {index} must be an object")
            continue
        a = str(branch.get("from") or "").strip()
        b = str(branch.get("to") or "").strip()
        if a not in node_set or b not in node_set or a == b:
            errors.append(f"network branch {index} has invalid endpoints")
        try:
            resistance = float(branch.get("resistance_K_per_W"))
        except (TypeError, ValueError):
            resistance = -1.0
        if resistance <= 0:
            errors.append(f"network branch {index} needs positive resistance_K_per_W")

    reachable = _reachable_nodes(node_set, branches, ambient)
    if node_set and reachable != node_set:
        errors.append("every thermal-network node must have a path to ambient")

    network_counts = Counter(network_regions)
    missing_network = sorted(set(losses) - set(network_counts))
    extra_network = sorted(set(network_counts) - set(losses))
    duplicate_network = sorted(name for name, count in network_counts.items() if count != 1)
    if missing_network:
        errors.append(f"network is missing loss regions: {missing_network}")
    if extra_network:
        errors.append(f"network has unknown loss regions: {extra_network}")
    if duplicate_network:
        errors.append(f"network loss regions must be assigned once: {duplicate_network}")

    mesh_counts: Counter[str] = Counter()
    mesh_loss_by_region: dict[str, float] = {}
    for index, row in enumerate(mesh_regions):
        if not isinstance(row, dict):
            errors.append(f"mesh region {index} must be an object")
            continue
        region = str(row.get("region") or "").strip()
        cell_type = str(row.get("cell_type") or "").strip().lower()
        if not region:
            errors.append(f"mesh region {index} needs a region name")
            continue
        mesh_counts[region] += 1
        if cell_type not in HEX_CELL_TYPES:
            errors.append(f"mesh region {region!r} must use a hexahedral cell type")
        try:
            cell_count = int(row.get("cell_count"))
        except (TypeError, ValueError):
            cell_count = 0
        if cell_count <= 0:
            errors.append(f"mesh region {region!r} needs a positive cell_count")
        try:
            mesh_loss = float(row.get("loss_W"))
        except (TypeError, ValueError):
            mesh_loss = -1.0
        if mesh_loss < 0:
            errors.append(f"mesh region {region!r} needs non-negative loss_W")
        mesh_loss_by_region[region] = mesh_loss

    missing_mesh = sorted(set(losses) - set(mesh_counts))
    extra_mesh = sorted(set(mesh_counts) - set(losses))
    duplicate_mesh = sorted(name for name, count in mesh_counts.items() if count != 1)
    if missing_mesh:
        errors.append(f"mesh is missing loss regions: {missing_mesh}")
    if extra_mesh:
        errors.append(f"mesh has unknown loss regions: {extra_mesh}")
    if duplicate_mesh:
        errors.append(f"mesh loss regions must be assigned once: {duplicate_mesh}")

    mismatched_mesh = sorted(
        name for name, value in losses.items()
        if name in mesh_loss_by_region
        and not _close(value, mesh_loss_by_region[name], relative_tolerance)
    )
    if mismatched_mesh:
        errors.append(f"mesh loss_W does not match source buckets: {mismatched_mesh}")

    source_total = sum(losses.values())
    network_total = sum(losses.get(name, 0.0) for name in network_regions)
    mesh_total = sum(value for value in mesh_loss_by_region.values() if value >= 0)
    totals_match = (
        _close(source_total, network_total, relative_tolerance)
        and _close(source_total, mesh_total, relative_tolerance)
    )
    if not totals_match:
        errors.append("source, network, and mesh total heat must agree")

    return {
        "schema": "radia-motor-thermal-handoff/v1",
        "policy": "loss_to_lptn_and_hex_thermal_mesh",
        "status": "ok" if not errors else "needs_attention",
        "source_total_loss_W": source_total,
        "network_total_loss_W": network_total,
        "mesh_total_loss_W": mesh_total,
        "loss_region_count": len(losses),
        "network_node_count": len(node_set),
        "network_branch_count": len(branches),
        "mesh_region_count": len(mesh_counts),
        "checks": {
            "network_regions_one_to_one": not (missing_network or extra_network or duplicate_network),
            "network_connected_to_ambient": bool(node_set) and reachable == node_set,
            "mesh_regions_one_to_one": not (missing_mesh or extra_mesh or duplicate_mesh),
            "mesh_all_hexahedral": all(
                isinstance(row, dict)
                and str(row.get("cell_type") or "").strip().lower() in HEX_CELL_TYPES
                for row in mesh_regions
            ),
            "regional_loss_values_match": not mismatched_mesh,
            "total_heat_conserved": totals_match,
        },
        "errors": errors,
        "next_steps": [
            "Run the LPTN for fast parameter sweeps and controller-oriented studies.",
            "Run the 3D all-hex thermal model for spatial gradients and hotspot resolution.",
            "Compare region-average temperatures and total boundary heat flow before promotion.",
        ],
    }


def motor_thermal_handoff_gate(
    loss_buckets_json: str,
    network_json: str,
    mesh_regions_json: str,
    relative_tolerance: float = 1.0e-9,
) -> str:
    """JSON wrapper used by the MCP tool."""
    losses = _load_json(loss_buckets_json, dict, "loss_buckets_json")
    network = _load_json(network_json, dict, "network_json")
    mesh = _load_json(mesh_regions_json, list, "mesh_regions_json")
    result = evaluate_motor_thermal_handoff(losses, network, mesh, relative_tolerance)
    return json.dumps(result, indent=2, sort_keys=True)
