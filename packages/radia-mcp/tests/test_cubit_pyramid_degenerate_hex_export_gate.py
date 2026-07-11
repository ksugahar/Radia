import copy
import json

from radia_mcp.cubit.server import cubit_pyramid_degenerate_hex_export_gate as mcp_gate
from radia_mcp.cubit.vol_inventory import cubit_pyramid_degenerate_hex_export_gate


def _deck(chexa, cpyram, chexa_nodes, chexa_unique, cpyram_nodes=()):
    return {
        "sha256": (str(chexa) + str(cpyram)) * 32,
        "card_counts": {"CHEXA": chexa, "CPYRAM": cpyram, "CTETRA": 0, "CPENTA": 0},
        "chexa_node_counts": list(chexa_nodes),
        "chexa_unique_node_counts": list(chexa_unique),
        "cpyram_node_counts": list(cpyram_nodes),
        "regular_chexa_count": sum(a == b for a, b in zip(chexa_nodes, chexa_unique)),
        "degenerate_chexa_count": sum(a > b for a, b in zip(chexa_nodes, chexa_unique)),
    }


def _summary():
    decks = [
        _deck(2, 2, [8, 8], [8, 8], [5, 5]),
        _deck(4, 0, [8, 8, 8, 8], [8, 8, 5, 5]),
        _deck(2, 2, [20, 20], [20, 20], [13, 13]),
        _deck(4, 0, [20, 20, 8, 8], [20, 20, 5, 5]),
    ]
    for index, deck in enumerate(decks):
        deck["sha256"] = str(index + 1) * 64
    return {
        "source_journal": "mixed_pyramid_export.jou",
        "source_sha256": "a" * 64,
        "execution_mode": "python_api_headless",
        "headless_flags": ["-nographics", "-batch"],
        "persistent_gui_started": False,
        "batch_wrapper_mode": "single_line_compile_wrapper",
        "direct_multiline_batch_rejected": True,
        "process_exit_code": 2,
        "startup_diagnostics": [
            "Could not open file: <install>/plugins",
            "Could not open file: -commandplugindir",
        ],
        "script_error_lines": [],
        "result_artifact_fresh": True,
        "element_counts": {"hex": 2, "pyramid": 2, "tet": 20, "wedge": 4},
        "block_inventory": {
            "1": {"hex": 2, "pyramid": 0, "tet": 0, "wedge": 0},
            "2": {"hex": 0, "pyramid": 2, "tet": 0, "wedge": 0},
        },
        "export_scope_claim": "registered_blocks_only",
        "quality": {
            "hex": {"count": 2, "minimum": 1.0},
            "pyramid": {
                "api_value_count": 0,
                "geometric_volume_count": 2,
                "geometric_volume_minimum": 0.2357,
            },
        },
        "pyramid_card_deck": decks[0],
        "nopyramid_deck": decks[1],
        "pyramid_card_deck_order2": decks[2],
        "nopyramid_deck_order2": decks[3],
        "pyramid_conversion_order_policy": "linearized_degenerate_chexa8",
        "order2_nopyramid_uniform_order_claimed": False,
        "total_cad_volume": 6.0,
    }


def test_accepts_explicit_order2_pyramid_linearization_and_block_subset():
    result = cubit_pyramid_degenerate_hex_export_gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["order2_nopyramid_chexa_node_counts"] == [20, 20, 8, 8]
    assert json.loads(mcp_gate(_summary()))["status"] == "ok"


def test_rejects_uniform_order_and_full_database_export_claims():
    bad = copy.deepcopy(_summary())
    bad["order2_nopyramid_uniform_order_claimed"] = True
    bad["export_scope_claim"] = "full_database"
    result = cubit_pyramid_degenerate_hex_export_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "export_scope_is_registered_blocks_only",
        "order2_nopyramid_linearization_recorded",
    }


def test_rejects_multiline_batch_and_unrelated_launcher_failure():
    bad = copy.deepcopy(_summary())
    bad["batch_wrapper_mode"] = "direct_multiline_batch"
    bad["direct_multiline_batch_rejected"] = False
    bad["startup_diagnostics"] = ["unknown fatal error"]
    result = cubit_pyramid_degenerate_hex_export_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["launcher_classification"] == "execution_error"
