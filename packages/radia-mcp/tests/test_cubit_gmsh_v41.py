import json

from radia_mcp.cubit.gmsh_v41 import (
    gmsh_v41_mixed_order_series_gate,
    summarize_gmsh_v41_ascii,
)
from radia_mcp.cubit.server import cubit_gmsh_v41_inventory


def _mesh(element_type, connectivity, node_count, *, transition_dim=0):
    tags = "\n".join(str(value) for value in range(1, node_count + 1))
    coords = "\n".join(f"{value} 0 0" for value in range(node_count))
    element_rows = [
        (3, 1, element_type["hex"], [connectivity["hex"]]),
        (3, 2, element_type["tet"], [connectivity["tet"]]),
        (3, 2, element_type["pyramid"], [connectivity["pyramid"]]),
        (2, 4, element_type["triangle"], [connectivity["triangle"]]),
        (2, 5, element_type["quad"], [connectivity["quad"]]),
    ]
    blocks = []
    element_tag = 1
    for dim, entity, kind, rows in element_rows:
        blocks.append(f"{dim} {entity} {kind} {len(rows)}")
        for count in rows:
            blocks.append(" ".join([str(element_tag)] + [str((i % node_count) + 1) for i in range(count)]))
            element_tag += 1
    return f"""$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
5
3 1 \"map\"
3 2 \"tet\"
{transition_dim} 3 \"pyram\"
2 4 \"boundary_tri\"
2 5 \"boundary_face\"
$EndPhysicalNames
$Nodes
1 {node_count} 1 {node_count}
3 1 0 {node_count}
{tags}
{coords}
$EndNodes
$Elements
5 5 1 5
{chr(10).join(blocks)}
$EndElements
"""


ORDER1_TYPES = {"hex": 5, "tet": 4, "pyramid": 7, "triangle": 2, "quad": 3}
ORDER1_NODES = {"hex": 8, "tet": 4, "pyramid": 5, "triangle": 3, "quad": 4}
ORDER2_TYPES = {"hex": 17, "tet": 11, "pyramid": 19, "triangle": 9, "quad": 16}
ORDER2_NODES = {"hex": 20, "tet": 10, "pyramid": 13, "triangle": 6, "quad": 8}


def _rows(transition_dim=0):
    first = summarize_gmsh_v41_ascii(_mesh(ORDER1_TYPES, ORDER1_NODES, 20, transition_dim=transition_dim))
    second = summarize_gmsh_v41_ascii(_mesh(ORDER2_TYPES, ORDER2_NODES, 60, transition_dim=transition_dim))
    return [{"order": 1, "inventory": first}, {"order": 2, "inventory": second}]


VOL = {
    "volume_kind_counts": {"hex": 1, "pyramid": 1, "tet": 1},
    "surface_kind_counts": {"quad": 1, "triangle": 1},
    "materials": {"1": "map", "2": "tet", "3": "pyram"},
}


def test_gmsh_v41_parser_reads_block_headers_and_high_order_connectivity():
    inventory = _rows()[1]["inventory"]
    assert inventory["status"] == "ok"
    assert inventory["element_type_counts"] == {"9": 1, "11": 1, "16": 1, "17": 1, "19": 1}
    assert inventory["connectivity_mismatches"] == []


def test_mixed_order_gate_uses_vol_when_transition_physical_dimension_is_incomplete():
    result = gmsh_v41_mixed_order_series_gate(_rows(), authoritative_vol_inventory=VOL)
    assert result["status"] == "ok_with_vol_metadata_authority"
    assert result["checks"]["gmsh_transition_physical_dimension_is_volume"] is False
    assert result["checks"]["vol_transition_metadata_authoritative"] is True


def test_mixed_order_gate_rejects_topology_drift():
    rows = _rows(transition_dim=3)
    rows[1]["inventory"]["element_family_counts"]["tet"] = 2
    result = gmsh_v41_mixed_order_series_gate(rows, authoritative_vol_inventory=VOL)
    assert result["status"] == "needs_attention"
    assert result["checks"]["topology_counts_invariant"] is False


def test_gmsh_v41_inventory_mcp_tool_rejects_v2_text():
    result = json.loads(cubit_gmsh_v41_inventory("$MeshFormat\n2.2 0 8\n$EndMeshFormat"))
    assert result["status"] == "invalid_input"
