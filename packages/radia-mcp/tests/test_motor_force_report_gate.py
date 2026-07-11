import json

from radia_mcp.motor.force_report_gate import (
    evaluate_force_report_method_metadata,
    force_report_method_metadata_gate,
)


def _report():
    return {
        "force_unit": "N",
        "component_frame": "global_cartesian",
        "methods": [
            {
                "family": "maxwell_stress",
                "domain": "closed_air_surface",
                "vector": [100.0, -20.0, 0.0],
            },
            {
                "family": "virtual_work_nodal",
                "domain": "moving_body_nodes",
                "vector": [100.5, -19.8, 0.0],
            },
        ],
        "action_force": [100.0, -20.0, 0.0],
        "reaction_force": [-100.1, 19.9, 0.0],
    }


def test_force_report_gate_accepts_independent_methods_and_balance():
    gate = evaluate_force_report_method_metadata(_report(), 0.02)
    assert gate["status"] == "ok"
    assert all(gate["checks"].values())
    assert json.loads(force_report_method_metadata_gate(json.dumps(_report()), 0.02))["status"] == "ok"


def test_force_report_gate_rejects_same_family_and_wrong_frame():
    report = _report()
    report["methods"][1]["family"] = "maxwell_stress"
    report["component_frame"] = "unspecified"
    gate = evaluate_force_report_method_metadata(report, 0.02)
    assert gate["status"] == "needs_attention"
    assert gate["checks"]["independent_force_families"] is False
    assert gate["checks"]["component_frame_is_explicit"] is False


def test_force_report_gate_rejects_broken_action_reaction():
    report = _report()
    report["reaction_force"] = [-80.0, 10.0, 0.0]
    gate = evaluate_force_report_method_metadata(report, 0.02)
    assert gate["status"] == "needs_attention"
    assert gate["checks"]["action_reaction_closes"] is False
