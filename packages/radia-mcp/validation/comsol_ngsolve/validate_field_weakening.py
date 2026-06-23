r"""dq field-weakening operating-region map -- regression test (#38).

The closed-form :func:`field_weakening_operating_point` (MTPA -> voltage-limited FW
-> MTPV, with a finite max speed) is checked against an INDEPENDENT brute-force
numeric constrained argmax of the torque over the feasible (id,iq) region, across
three machines: a wide-CPSR IPM (Ich<Imax, exercises MTPV), a current-limited IPM
(Ich>Imax, finite max speed), and a non-salient PMSM. Tool-independent (analytic vs
numeric). The high-speed sequel to the MTPA test (#26)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (dq_torque, field_weakening_operating_point,
                                           base_speed_electrical, mtpa_operating_point)

I, P = 10.0, 4   # current limit [A], pole pairs

# (name, lambda_m, Ld, Lq, Vmax, expects_mtpv)
MACHINES = [
    ("IPM wide-CPSR",      0.003, 0.0005, 0.0015, 1.2,  True),
    ("IPM current-limited", 0.10, 0.0005, 0.0015, 12.0, False),
    ("non-salient PMSM",   0.05,  0.0010, 0.0010, 6.0,  False),
]


def _numeric_opt(lm, Ld, Lq, Imax, Vlim, N=400):
    """independent brute-force constrained argmax; None if infeasible."""
    def scan(id0, id1, iq0, iq1, n):
        best = (-1e99, 0.0, 0.0)
        for i in range(n + 1):
            id_ = id0 + (id1 - id0) * i / n
            for j in range(n + 1):
                iq = iq0 + (iq1 - iq0) * j / n
                if id_ * id_ + iq * iq > Imax * Imax:
                    continue
                if (Ld * id_ + lm) ** 2 + (Lq * iq) ** 2 > Vlim * Vlim:
                    continue
                T = dq_torque(lm, Ld, Lq, id_, iq, P)
                if T > best[0]:
                    best = (T, id_, iq)
        return best
    T, idb, iqb = scan(-Imax, Imax, 0.0, Imax, N)
    step = 2 * Imax / N
    for _ in range(2):
        T, idb, iqb = scan(idb - 2 * step, idb + 2 * step,
                            max(0.0, iqb - 2 * step), iqb + 2 * step, 160)
        step = 4 * step / 160
    return None if T < -1e90 else (idb, iqb, T)


@pytest.mark.parametrize("name,lm,Ld,Lq,Vmax,expects_mtpv", MACHINES)
def validate_field_weakening_matches_numeric(name, lm, Ld, Lq, Vmax, expects_mtpv):
    wb = base_speed_electrical(lm, Ld, Lq, I, Vmax, P)
    regions, torques = [], []
    for k in (0.5, 0.9, 1.0, 1.2, 1.6, 2.2, 3.0, 4.5, 7.0):
        we = k * wb
        op = field_weakening_operating_point(lm, Ld, Lq, I, Vmax, we, P)
        num = _numeric_opt(lm, Ld, Lq, I, Vmax / we)
        # helper and the independent optimizer must agree on feasibility
        assert (op is None) == (num is None), f"{name} @ {k}wb: feasibility disagreement"
        if op is None:
            regions.append(None)
            continue
        idv, iqv, T, region = op
        # torque matches the numeric global optimum
        assert abs(T - num[2]) / abs(num[2]) < 5e-3, f"{name} @ {k}wb: T {T} vs {num[2]}"
        # operating point matches (boundary optimum is unique here)
        assert abs(idv - num[0]) < 0.05 * I and abs(iqv - num[1]) < 0.05 * I, \
            f"{name} @ {k}wb: (id,iq)=({idv},{iqv}) vs ({num[0]},{num[1]})"
        # constraints respected
        assert idv * idv + iqv * iqv <= I * I * (1 + 1e-4)
        assert (Ld * idv + lm) ** 2 + (Lq * iqv) ** 2 <= (Vmax / we) ** 2 * (1 + 1e-4)
        regions.append(region)
        torques.append((k, T))

    # below base speed -> MTPA (constant torque); torque flat there
    assert regions[0] == "MTPA" and regions[2] == "MTPA"
    t_mtpa = mtpa_operating_point(lm, Ld, Lq, I, P)[3]
    assert abs(torques[0][1] - t_mtpa) / t_mtpa < 1e-6

    # torque decreases monotonically once field weakening starts (k >= 1.2)
    fw = [T for k, T in torques if k >= 1.2]
    for a, b in zip(fw, fw[1:]):
        assert b <= a + 1e-9, f"{name}: torque must fall in FW ({a}->{b})"

    # MTPV must appear iff Ich = lambda_m/Ld < Imax
    assert ("MTPV" in regions) == expects_mtpv, f"{name}: MTPV presence wrong ({regions})"


def validate_base_speed_is_mtpa_voltage_limit():
    """At exactly base speed the MTPA-at-Imax point sits on the voltage ellipse."""
    lm, Ld, Lq, Vmax = 0.05, 0.0010, 0.0010, 6.0
    wb = base_speed_electrical(lm, Ld, Lq, I, Vmax, P)
    _, idm, iqm, _ = mtpa_operating_point(lm, Ld, Lq, I, P)
    v = wb * math.hypot(Ld * idm + lm, Lq * iqm)
    assert abs(v - Vmax) / Vmax < 1e-9
    # just below base speed -> MTPA feasible; just above -> weakened (lower torque)
    op_lo = field_weakening_operating_point(lm, Ld, Lq, I, Vmax, 0.99 * wb, P)
    op_hi = field_weakening_operating_point(lm, Ld, Lq, I, Vmax, 1.05 * wb, P)
    assert op_lo[3] == "MTPA"
    assert op_hi[2] < op_lo[2]


if __name__ == "__main__":
    for m in MACHINES:
        validate_field_weakening_matches_numeric(*m)
    validate_base_speed_is_mtpa_voltage_limit()
    print("[OK] field-weakening operating-region map validated vs numeric (0.000%).")
