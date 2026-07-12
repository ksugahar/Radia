from __future__ import annotations

import copy
import json
import math

from radia_mcp.radia_ngsolve.nonlinear_bh_curve_gate import (
    nonlinear_bh_piecewise_material_gate,
)
from radia_mcp.radia_ngsolve.server import nonlinear_bh_piecewise_material_gate as mcp_gate


def _summary() -> dict:
    mu0 = 4.0e-7 * math.pi
    h = [0.0, 10.0, 100.0, 1000.0]
    interval_mu = [1000.0, 100.0, 1.0]
    b = [0.0]
    for left, right, mu_diff in zip(h[:-1], h[1:], interval_mu, strict=True):
        b.append(b[-1] + mu0 * mu_diff * (right - left))
    return {
        "contract": {
            "interpolation": "piecewise_linear",
            "differential_interval": "left_interval_ending_at_knot",
            "secant_definition": "B/(mu0*H)",
            "differential_definition": "deltaB/(mu0*deltaH)",
            "saturation_tail_expected": True,
        },
        "bh_rows": [
            {"h_a_per_m": h_value, "b_t": b_value}
            for h_value, b_value in zip(h, b, strict=True)
        ],
        "secant_rows": [
            {"h_a_per_m": h_value, "relative_mu": b_value / (mu0 * h_value)}
            for h_value, b_value in zip(h[1:], b[1:], strict=True)
        ],
        "differential_rows": [
            {"h_a_per_m": h_value, "relative_mu_diff": mu_diff}
            for h_value, mu_diff in zip(h[1:], interval_mu, strict=True)
        ],
    }


def test_piecewise_bh_gate_accepts_secant_and_left_interval_tangent() -> None:
    result = nonlinear_bh_piecewise_material_gate(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["secant_permeability_identity"] is True
    assert result["checks"]["differential_permeability_identity"] is True
    assert json.loads(mcp_gate(_summary()))["status"] == "ok"


def test_piecewise_bh_gate_rejects_central_difference_substitution() -> None:
    bad = copy.deepcopy(_summary())
    bad["differential_rows"][1]["relative_mu_diff"] = 0.5 * (1000.0 + 100.0)
    bad["contract"]["differential_interval"] = "central_difference"
    result = nonlinear_bh_piecewise_material_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["piecewise_linear_contract_recorded"] is False
    assert result["checks"]["differential_permeability_identity"] is False


def test_piecewise_bh_gate_rejects_wrong_secant_and_nonmonotone_b() -> None:
    bad = _summary()
    bad["secant_rows"][0]["relative_mu"] *= 1.2
    bad["bh_rows"][2]["b_t"] = -1.0
    result = nonlinear_bh_piecewise_material_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["b_is_monotone_nondecreasing"] is False
    assert result["checks"]["secant_permeability_identity"] is False
