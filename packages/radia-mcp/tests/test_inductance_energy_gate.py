import json

from radia_mcp.radia_ngsolve.inductance_energy_gate import inductance_energy_mutual_gate
from radia_mcp.radia_ngsolve.server import inductance_energy_mutual_gate as mcp_gate


def _args():
    return {
        "self_inductance": 619.447084639947,
        "energy_inductance": 619.4470846399607,
        "mutual_inductance": 1.974711005479,
        "analytic_mutual_inductance": 1.973920879949,
        "inductance_unit": "nH",
    }


def test_accepts_energy_and_analytic_mutual_closure():
    result = inductance_energy_mutual_gate(**_args())
    assert result["status"] == "ok"
    assert json.loads(mcp_gate(**_args()))["status"] == "ok"
    assert result["conventions"]["reciprocity_status"] == "not_tested_by_one_direction_gate"


def test_rejects_energy_or_analytic_drift():
    result = inductance_energy_mutual_gate(
        **{**_args(), "energy_inductance": 500.0, "analytic_mutual_inductance": 1.5}
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["energy_identity_closes"] is False
    assert result["checks"]["mutual_analytic_reference_closes"] is False


def test_rejects_ambiguous_unit():
    result = inductance_energy_mutual_gate(**{**_args(), "inductance_unit": "model unit"})
    assert result["checks"]["inductance_unit_explicit"] is False
