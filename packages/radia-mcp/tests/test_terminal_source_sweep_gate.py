import json

from radia_mcp.radia_ngsolve.server import cyclic_terminal_source_sweep_gate


def good() -> dict:
    return {
        "formulations": [
            {
                "label": "volume",
                "source_sweep_order": [1, 2, 3, 4, 5],
                "active_terminal_charge_c": [2.0e-12] * 5,
            },
            {
                "label": "boundary",
                "source_sweep_order": [1, 2, 3, 4, 5],
                "active_terminal_charge_c": [1.4e-11] * 5,
            },
        ]
    }


def test_accepts_cyclic_sweeps_without_forcing_equal_absolute_values():
    result = json.loads(cyclic_terminal_source_sweep_gate(json.dumps(good())))
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_to_minimum_formulation_ratio"] == 7


def test_rejects_broken_outer_order_and_cyclic_charge():
    summary = good()
    summary["formulations"][1]["source_sweep_order"] = [1, 2, 3, 4, 4]
    summary["formulations"][1]["active_terminal_charge_c"][2] *= 1.2
    result = json.loads(cyclic_terminal_source_sweep_gate(json.dumps(summary)))
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "formulation_2_order_is_permutation",
        "formulation_2_cyclic_symmetry",
    }


def test_rejects_empty_charge_vectors_instead_of_vacuously_accepting_them():
    summary = good()
    summary["formulations"][0]["active_terminal_charge_c"] = []
    summary["formulations"][0]["source_sweep_order"] = []
    result = json.loads(cyclic_terminal_source_sweep_gate(json.dumps(summary)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["formulation_1_charges_positive_finite"] is False
