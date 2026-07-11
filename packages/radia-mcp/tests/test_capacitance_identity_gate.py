import json

from radia_mcp.radia_ngsolve.capacitance_identity_gate import (
    two_conductor_capacitance_identity_gate,
)
from radia_mcp.radia_ngsolve.server import (
    two_conductor_capacitance_identity_gate as mcp_gate,
)


VOLTAGES = [0.0, 0.999999999824475]
CHARGES = [-2.266192540666828e-11, 2.266192711409279e-11]
ENERGY = 1.133096376703627e-11


def test_capacitance_identity_accepts_live_two_conductor_result():
    result = two_conductor_capacitance_identity_gate(
        VOLTAGES,
        CHARGES,
        ENERGY,
        planar_depth_m=1.0,
    )
    assert result["status"] == "ok"
    assert result["capacitance_relative_error"] < 2.0e-8
    assert result["charge_balance_relative_error"] < 8.0e-8
    assert result["planar_depth_m"] == 1.0


def test_capacitance_identity_rejects_stale_energy_and_same_sign_charge():
    result = two_conductor_capacitance_identity_gate(
        VOLTAGES,
        [abs(CHARGES[0]), CHARGES[1]],
        ENERGY * 0.5,
        planar_depth_m=1.0,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["charge_energy_capacitance_agree"] is False
    assert result["checks"]["conductor_charge_balance_ok"] is False


def test_capacitance_identity_mcp_dispatches_json():
    result = json.loads(mcp_gate(VOLTAGES, CHARGES, ENERGY, 1, 1.0))
    assert result["status"] == "ok"
    assert result["policy"] == "two_conductor_capacitance_identity_gate_v1"
