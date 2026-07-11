import copy
import json

from radia_mcp.cubit.high_order_export_gate import (
    cubit_headless_netgen_export_gate,
    cubit_loft_high_order_vol_series_gate,
)
from radia_mcp.cubit.server import (
    cubit_headless_netgen_export_gate as mcp_headless_gate,
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


def _summary():
    return {
        "source_journal": "05_loft.jou",
        "headless": True,
        "persistent_gui_started": False,
        "legacy_command": "radia_export",
        "replay_command": "export netgen",
        "legacy_command_available_headless": False,
        "native_command_available_headless": True,
        "produced_orders": [1, 2, 3, 4, 5],
        "vol_file_count": 5,
        "sidecar_file_count": 5,
        "artifact_set_complete": True,
        "hex_count": 24,
        "quality_minimum": 0.6565,
        "exit_code": 2,
        "known_startup_plugin_path_diagnostics": True,
        "model_or_export_errors": False,
    }


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


def test_headless_native_export_accepts_classified_exit_two_and_dispatches():
    result = cubit_headless_netgen_export_gate(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["exit_code_explained_by_known_startup_diagnostics"] is True
    assert json.loads(mcp_headless_gate(_summary()))["status"] == "ok"


def test_headless_export_rejects_gui_plugin_command_and_missing_order():
    bad = _summary()
    bad["replay_command"] = "radia_export"
    bad["produced_orders"] = [1, 2, 3, 4]
    bad["vol_file_count"] = 4
    result = cubit_headless_netgen_export_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["native_export_netgen_command_used"] is False
    assert result["checks"]["all_expected_orders_and_sidecars_created"] is False
