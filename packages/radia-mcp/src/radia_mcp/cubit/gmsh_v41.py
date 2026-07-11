"""Small, strict inventory helpers for ASCII Gmsh 4.1 mesh exports."""
from __future__ import annotations

from collections import Counter
import shlex
from typing import Mapping


ELEMENT_TYPES = {
    1: ("line", 1, 2),
    2: ("triangle", 1, 3),
    3: ("quad", 1, 4),
    4: ("tet", 1, 4),
    5: ("hex", 1, 8),
    6: ("wedge", 1, 6),
    7: ("pyramid", 1, 5),
    8: ("line", 2, 3),
    9: ("triangle", 2, 6),
    10: ("quad", 2, 9),
    11: ("tet", 2, 10),
    12: ("hex", 2, 27),
    13: ("wedge", 2, 18),
    14: ("pyramid", 2, 14),
    15: ("point", 1, 1),
    16: ("quad", 2, 8),
    17: ("hex", 2, 20),
    18: ("wedge", 2, 15),
    19: ("pyramid", 2, 13),
}


def _section(lines: list[str], name: str) -> list[str]:
    start_marker = f"${name}"
    end_marker = f"$End{name}"
    try:
        start = lines.index(start_marker) + 1
        end = lines.index(end_marker, start)
    except ValueError as exc:
        raise ValueError(f"missing Gmsh section {start_marker}") from exc
    return lines[start:end]


def summarize_gmsh_v41_ascii(text: str, *, source: str = "inline") -> dict:
    """Return topology, order and connectivity inventory from ASCII Gmsh 4.1."""

    lines = [line.strip() for line in str(text).splitlines() if line.strip()]
    mesh_format = _section(lines, "MeshFormat")
    if len(mesh_format) != 1:
        raise ValueError("Gmsh MeshFormat must contain exactly one row")
    format_parts = mesh_format[0].split()
    if len(format_parts) != 3:
        raise ValueError("invalid Gmsh MeshFormat row")
    version, binary_flag, data_size = format_parts
    if version != "4.1":
        raise ValueError(f"expected Gmsh 4.1, got {version}")
    if binary_flag != "0":
        raise ValueError("only ASCII Gmsh 4.1 is supported")

    physical_names = []
    if "$PhysicalNames" in lines:
        section = _section(lines, "PhysicalNames")
        expected_names = int(section[0])
        for row in section[1:]:
            parts = shlex.split(row)
            if len(parts) != 3:
                raise ValueError(f"invalid PhysicalNames row: {row}")
            physical_names.append({"dimension": int(parts[0]), "tag": int(parts[1]), "name": parts[2]})
        if len(physical_names) != expected_names:
            raise ValueError("PhysicalNames count does not match its header")

    nodes = _section(lines, "Nodes")
    node_header = [int(value) for value in nodes[0].split()]
    if len(node_header) != 4:
        raise ValueError("invalid Gmsh Nodes header")
    node_blocks, node_count, min_node_tag, max_node_tag = node_header

    elements = _section(lines, "Elements")
    element_header = [int(value) for value in elements[0].split()]
    if len(element_header) != 4:
        raise ValueError("invalid Gmsh Elements header")
    declared_blocks, declared_elements, min_element_tag, max_element_tag = element_header
    index = 1
    blocks = []
    type_counts: Counter[int] = Counter()
    family_counts: Counter[str] = Counter()
    volume_counts: Counter[str] = Counter()
    surface_counts: Counter[str] = Counter()
    order_counts: Counter[int] = Counter()
    connectivity_mismatches = []
    parsed_elements = 0
    for block_index in range(declared_blocks):
        if index >= len(elements):
            raise ValueError("Elements section ended before all blocks were read")
        header = [int(value) for value in elements[index].split()]
        index += 1
        if len(header) != 4:
            raise ValueError(f"invalid element block header at block {block_index + 1}")
        entity_dim, entity_tag, element_type, block_count = header
        info = ELEMENT_TYPES.get(element_type)
        family, order, expected_nodes = info if info else (f"type_{element_type}", None, None)
        for _ in range(block_count):
            if index >= len(elements):
                raise ValueError("Elements section ended inside a block")
            record = [int(value) for value in elements[index].split()]
            index += 1
            connectivity_count = max(0, len(record) - 1)
            if expected_nodes is not None and connectivity_count != expected_nodes:
                connectivity_mismatches.append({
                    "element_tag": record[0] if record else None,
                    "element_type": element_type,
                    "expected_nodes": expected_nodes,
                    "actual_nodes": connectivity_count,
                })
        type_counts[element_type] += block_count
        family_counts[family] += block_count
        if entity_dim == 3:
            volume_counts[family] += block_count
        if entity_dim == 2:
            surface_counts[family] += block_count
        if order is not None:
            order_counts[order] += block_count
        parsed_elements += block_count
        blocks.append({
            "entity_dimension": entity_dim,
            "entity_tag": entity_tag,
            "element_type": element_type,
            "family": family,
            "order": order,
            "connectivity_nodes": expected_nodes,
            "count": block_count,
        })
    if index != len(elements):
        raise ValueError("unexpected trailing rows in Elements section")

    checks = {
        "mesh_format_is_ascii_v41": version == "4.1" and binary_flag == "0",
        "element_block_count_matches": len(blocks) == declared_blocks,
        "element_count_matches": parsed_elements == declared_elements,
        "known_connectivity_matches": not connectivity_mismatches,
        "node_count_positive": node_count > 0,
    }
    return {
        "source": source,
        "policy": "gmsh_v41_ascii_inventory_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "mesh_format": version,
        "binary": False,
        "data_size": int(data_size),
        "node_block_count": node_blocks,
        "node_count": node_count,
        "min_node_tag": min_node_tag,
        "max_node_tag": max_node_tag,
        "element_block_count": declared_blocks,
        "element_count": declared_elements,
        "min_element_tag": min_element_tag,
        "max_element_tag": max_element_tag,
        "physical_names": physical_names,
        "element_blocks": blocks,
        "element_type_counts": {str(key): value for key, value in sorted(type_counts.items())},
        "element_family_counts": dict(sorted(family_counts.items())),
        "volume_family_counts": dict(sorted(volume_counts.items())),
        "surface_family_counts": dict(sorted(surface_counts.items())),
        "element_order_counts": {str(key): value for key, value in sorted(order_counts.items())},
        "connectivity_mismatches": connectivity_mismatches,
        "checks": checks,
    }


def gmsh_v41_mixed_order_series_gate(
    rows,
    *,
    authoritative_vol_inventory: Mapping | None = None,
    transition_names=("pyram", "pyramid_transition", "transition"),
) -> dict:
    """Check order elevation and metadata fidelity for a mixed Gmsh 4.1 series."""

    series = []
    for raw in rows:
        if not isinstance(raw, Mapping):
            raise ValueError("each order-series row must be a mapping")
        order = int(raw.get("order", 0))
        inventory = raw.get("inventory")
        if not isinstance(inventory, Mapping):
            raise ValueError("each row must contain an inventory mapping")
        series.append((order, inventory))
    series.sort(key=lambda item: item[0])
    orders = [order for order, _ in series]
    if len(set(orders)) != len(orders):
        raise ValueError("order values must be unique")

    signatures = [dict(inv.get("element_family_counts", {})) for _, inv in series]
    physical_signatures = [
        sorted((int(row["dimension"]), int(row["tag"]), str(row["name"])) for row in inv.get("physical_names", []))
        for _, inv in series
    ]
    order_matches = []
    for requested_order, inventory in series:
        blocks = inventory.get("element_blocks", [])
        order_matches.append(bool(blocks) and all(int(block.get("order") or 0) == requested_order for block in blocks))

    first_inventory = series[0][1] if series else {}
    volume = dict(first_inventory.get("volume_family_counts", {}))
    surface = dict(first_inventory.get("surface_family_counts", {}))
    node_counts = [int(inv.get("node_count", 0) or 0) for _, inv in series]
    transition_set = {str(name).lower() for name in transition_names}
    physical_rows = first_inventory.get("physical_names", [])
    transition_physical_rows = [
        row for row in physical_rows if str(row.get("name", "")).lower() in transition_set
    ]
    gmsh_transition_volume_dim = bool(transition_physical_rows) and all(
        int(row.get("dimension", -1)) == 3 for row in transition_physical_rows
    )

    vol = dict(authoritative_vol_inventory or {})
    vol_counts = dict(vol.get("volume_kind_counts", {}))
    vol_surface_counts = dict(vol.get("surface_kind_counts", {}))
    vol_materials = {str(value).lower() for value in dict(vol.get("materials", {})).values()}
    vol_transition_authority = (
        int(vol_counts.get("pyramid", 0) or 0) > 0
        and bool(vol_materials & transition_set)
    )
    checks = {
        "orders_are_one_then_two": orders == [1, 2],
        "all_inventories_are_ascii_v41": bool(series) and all(
            inv.get("mesh_format") == "4.1" and inv.get("binary") is False for _, inv in series
        ),
        "all_connectivity_matches_declared_type": bool(series) and all(
            not inv.get("connectivity_mismatches") for _, inv in series
        ),
        "requested_order_matches_element_types": bool(order_matches) and all(order_matches),
        "topology_counts_invariant": bool(signatures) and all(signature == signatures[0] for signature in signatures),
        "physical_names_invariant": bool(physical_signatures) and all(signature == physical_signatures[0] for signature in physical_signatures),
        "higher_order_adds_nodes": len(node_counts) == 2 and node_counts[1] > node_counts[0],
        "hex_primary_present": int(volume.get("hex", 0) or 0) > 0,
        "pyramid_transition_present": int(volume.get("pyramid", 0) or 0) > 0,
        "tet_compatibility_present": int(volume.get("tet", 0) or 0) > 0,
        "surface_tri_and_quad_present": int(surface.get("triangle", 0) or 0) > 0 and int(surface.get("quad", 0) or 0) > 0,
        "vol_topology_matches_gmsh": not vol or (
            vol_counts == volume and vol_surface_counts == surface
        ),
        "vol_transition_metadata_authoritative": not vol or vol_transition_authority,
        "gmsh_transition_physical_dimension_is_volume": gmsh_transition_volume_dim,
    }
    hard_names = [
        "orders_are_one_then_two",
        "all_inventories_are_ascii_v41",
        "all_connectivity_matches_declared_type",
        "requested_order_matches_element_types",
        "topology_counts_invariant",
        "physical_names_invariant",
        "higher_order_adds_nodes",
        "hex_primary_present",
        "pyramid_transition_present",
        "tet_compatibility_present",
        "surface_tri_and_quad_present",
        "vol_topology_matches_gmsh",
        "vol_transition_metadata_authoritative",
    ]
    hard_ok = all(checks[name] for name in hard_names)
    if hard_ok and checks["gmsh_transition_physical_dimension_is_volume"]:
        status = "ok"
    elif hard_ok and vol_transition_authority:
        status = "ok_with_vol_metadata_authority"
    else:
        status = "needs_attention"
    return {
        "policy": "gmsh_v41_mixed_order_series_gate_v1",
        "status": status,
        "success": hard_ok and (checks["gmsh_transition_physical_dimension_is_volume"] or vol_transition_authority),
        "orders": orders,
        "node_counts": node_counts,
        "element_family_counts": signatures,
        "volume_family_counts": volume,
        "surface_family_counts": surface,
        "transition_physical_names": transition_physical_rows,
        "checks": checks,
        "warnings": [] if checks["gmsh_transition_physical_dimension_is_volume"] else [
            "The transition physical name is not a 3D Gmsh group; keep the labeled .vol inventory authoritative and use .msh for mesh inspection."
        ],
        "lesson": (
            "Parse Gmsh 4.1 by entity blocks, not as v2 element rows. Order elevation must preserve mixed topology and declared connectivity. "
            "When transition-block dimensions are incomplete, retain the labeled .vol as solver metadata authority."
        ),
    }
