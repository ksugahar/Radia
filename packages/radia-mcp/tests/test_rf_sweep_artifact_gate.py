import json

from radia_mcp.radia_ngsolve.rf_sweep_artifact_gate import rf_sweep_artifact_summary_gate


def _summary():
    return {"solved": True, "touchstone_fresh": True, "touchstone_suffix": ".s2p",
            "frequency_rows": 1001, "matrix_channels": ["S11", "S12", "S21", "S22"],
            "max_singular_value": 0.9997, "max_reciprocity_abs": 1.2e-4,
            "solver_version": "2026", "run_date_utc": "2026-07-11T00:00:00Z"}


def test_rf_sweep_artifact_summary_accepts_complete_passive_result():
    gate = rf_sweep_artifact_summary_gate(json.dumps(_summary()))
    assert gate["status"] == "ok" and all(gate["checks"].values())


def test_rf_sweep_artifact_summary_rejects_stale_active_result():
    row = _summary(); row["touchstone_fresh"] = False; row["max_singular_value"] = 1.1
    gate = rf_sweep_artifact_summary_gate(json.dumps(row))
    assert gate["status"] == "needs_attention"
    assert not gate["checks"]["touchstone_is_fresh"] and not gate["checks"]["passive"]
