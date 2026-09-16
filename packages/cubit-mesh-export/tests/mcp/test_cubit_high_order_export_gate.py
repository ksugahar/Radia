import copy
import json

from cubit_mesh_export.mcp.high_order_export_gate import (
    cubit_loft_high_order_vol_series_gate,
)
from cubit_mesh_export.mcp.server import (
    cubit_loft_high_order_vol_series_gate as mcp_series_gate,
)


def _inventory(order):
    return {
        "surface_section": "surfaceelementsuv",
        "surface_elements": 40,
        "surface_kind_counts": {"quad": 40},
        "volume_elements": 24,
        "volume_kind_counts": {"hex": 24},
        "points": 51,
        "materials": {"1": "map"},
        "boundary_names": {str(i): f"Surface_{i}" for i in range(1, 7)},
        "curvedelements_present": order > 1,
        "routing_hint": "cubit_hex_or_mixed_path",
    }


def _rows():
    curved = [0, 118, 1164, 2658, 4952]
    return [
        {
            "order": order,
            "curved_node_count": curved[order - 1],
            "quality_minimum": 0.6565,
            "inventory": _inventory(order),
            "sidecar": {
                "order": order,
                "n_elements": 24,
                "n_points": 51,
                "materials": {"map": 1.980740e-9},
                "boundaries": {f"Surface_{i}": 1.0 for i in range(1, 7)},
            },
        }
        for order in range(1, 6)
    ]


def test_live_shape_high_order_series_passes_and_mcp_dispatches():
    result = cubit_loft_high_order_vol_series_gate(_rows())
    assert result["status"] == "ok"
    assert result["curved_node_counts"] == [0, 118, 1164, 2658, 4952]
    assert json.loads(mcp_series_gate(_rows()))["status"] == "ok"


def test_high_order_series_rejects_topology_and_curving_drift():
    bad = copy.deepcopy(_rows())
    bad[-1]["inventory"]["volume_elements"] = 23
    bad[-1]["curved_node_count"] = 100
    result = cubit_loft_high_order_vol_series_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["topology_and_labels_invariant"] is False
    assert result["checks"]["curved_nodes_start_zero_then_strictly_increase"] is False
