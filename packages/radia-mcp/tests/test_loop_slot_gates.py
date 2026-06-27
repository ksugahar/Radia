import math

import pytest

from radia_mcp.radia_ngsolve.slot_gates import (
    coenergy_torque_periodic_summary,
    parallel_wire_force_per_length,
    two_port_sparameter_health,
)


def test_parallel_wire_force_gate_signed_and_scaled():
    f = parallel_wire_force_per_length(10.0, 20.0, 0.05)
    assert f == pytest.approx(8.0e-4)
    assert parallel_wire_force_per_length(10.0, -20.0, 0.05) == pytest.approx(-f)
    assert parallel_wire_force_per_length(10.0, 20.0, 0.10) == pytest.approx(0.5 * f)
    with pytest.raises(ValueError):
        parallel_wire_force_per_length(1.0, 1.0, 0.0)


def test_coenergy_torque_gate_uses_absolute_tolerance_at_zero_crossings():
    n = 64
    amp = 0.25
    theta = [2.0 * math.pi * i / n for i in range(n)]
    # W' = -A cos(theta); T = dW'/dtheta = A sin(theta).
    coenergy = [-amp * math.cos(t) for t in theta]
    torque = [amp * math.sin(t) for t in theta]

    summary = coenergy_torque_periodic_summary(
        theta,
        coenergy,
        torque,
        rtol=2.0e-3,
        atol=1.0e-12,
    )

    assert summary["status"] == "ok"
    assert summary["max_abs_error"] < 5.0e-4
    zero_rows = [row for row in summary["rows"] if abs(row["reference_torque_nm"]) < 1.0e-12]
    assert zero_rows
    assert max(row["abs_error"] for row in zero_rows) < 1.0e-12


def test_two_port_sparameter_health_checks_passivity_and_reciprocity():
    health = two_port_sparameter_health(0.1, 0.7, s12=0.7, s22=0.1)
    assert health["status"] == "ok"
    assert health["reciprocal"] is True
    assert health["passive"] is True
    assert health["passive_margin"] > 0.0

    nonreciprocal = two_port_sparameter_health(0.1, 0.7, s12=0.6, s22=0.1)
    assert nonreciprocal["reciprocal"] is False
    assert nonreciprocal["status"] == "needs_attention"

    active = two_port_sparameter_health(0.2, 1.1, s12=1.1, s22=0.2)
    assert active["passive"] is False
    assert active["status"] == "needs_attention"
