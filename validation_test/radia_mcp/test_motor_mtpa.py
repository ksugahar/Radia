r"""Maximum-Torque-Per-Ampere (MTPA) operating point -- regression test.

The dq control capstone of the salient PM synchronous machine: given the magnet
flux linkage lambda_m and the (possibly frozen-permeability) Ld/Lq, the MTPA
current angle maximises torque per ampere.  dq torque (motor convention):

    T = (3/2) p [ lambda_m iq + (Ld - Lq) id iq ],   id = -I sin g, iq = I cos g.

``solve.mtpa_operating_point`` returns the closed-form MTPA (root of dT/dg=0).
This locks it against an independent NUMERICAL argmax over the current angle, plus
the optimality residual and the physical limits.  Pure dq theory -> tool-independent
(no FE solve, no commercial reference); the example mtpa_locus.py feeds it
FE-extracted lambda_m/Ld/Lq.
"""
import math
import os
import sys

import pytest

pytest.importorskip("ngsolve")          # solve.py imports ngsolve at module load

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import dq_torque, mtpa_operating_point


def _numeric_argmax(lambda_m, Ld, Lq, I, p, n=400001):
    best = None
    for i in range(n):
        g = -math.pi / 2 + math.pi * i / (n - 1)         # gamma in (-90, 90) deg
        T = dq_torque(lambda_m, Ld, Lq, -I * math.sin(g), I * math.cos(g), p)
        if best is None or T > best[1]:
            best = (g, T)
    return best


# (name, lambda_m [Wb], Ld [H], Lq [H], I [A], pole_pairs)
CASES = [
    ("IPM_strong_saliency", 0.12, 2.0e-3, 6.0e-3, 8.0, 4),   # Ld<Lq -> gamma>0
    ("synrel_no_magnet",    0.0,  5.0e-4, 3.0e-4, 8.0, 2),   # lambda_m=0, Ld>Lq -> -45 deg
    ("salient_pm_Ld_gt_Lq", 0.08, 3.6e-5, 2.95e-5, 6.0, 3),  # radia ld_lq-like saliency
    ("non_salient",         0.15, 4.0e-4, 4.0e-4, 8.0, 4),   # Ld=Lq -> pure q-axis
]


@pytest.mark.parametrize("name,lm,Ld,Lq,I,p", CASES, ids=[c[0] for c in CASES])
def test_mtpa_matches_numeric_argmax(name, lm, Ld, Lq, I, p):
    g, id_, iq, T = mtpa_operating_point(lm, Ld, Lq, I, p)
    gn, Tn = _numeric_argmax(lm, Ld, Lq, I, p)
    # closed-form MTPA torque == numeric maximum
    assert abs(T - Tn) <= 1e-6 * abs(Tn) + 1e-12, f"{name}: T {T} vs numeric {Tn}"
    assert abs(g - gn) < 0.05 * math.pi / 180 or abs(math.sin(g) - math.sin(gn)) < 1e-4
    # stationarity: lambda_m sin g + (Ld-Lq) I cos 2g = 0
    res = lm * math.sin(g) + (Ld - Lq) * I * math.cos(2 * g)
    assert abs(res) < 1e-6 * max(1.0, abs(lm), abs(Ld - Lq) * I), f"{name}: optimality res {res}"
    # current magnitude is honoured and MTPA >= pure q-axis (saliency never hurts)
    assert id_ ** 2 + iq ** 2 == pytest.approx(I * I, rel=1e-9)
    assert T >= dq_torque(lm, Ld, Lq, 0.0, I, p) - 1e-12


def test_mtpa_limits():
    # synchronous reluctance (no magnet, Ld>Lq): MTPA current angle = -45 deg
    g, id_, iq, T = mtpa_operating_point(0.0, 5.0e-4, 3.0e-4, 10.0, 2)
    assert math.degrees(g) == pytest.approx(-45.0, abs=1e-3)
    assert id_ > 0 and iq > 0 and T > 0                       # both axes drive reluctance torque
    # non-salient PM (Ld=Lq): all torque is magnet torque -> pure q-axis (gamma=0)
    g0, id0, iq0, T0 = mtpa_operating_point(0.2, 1e-3, 1e-3, 7.0, 4)
    assert math.degrees(g0) == pytest.approx(0.0, abs=1e-6)
    assert id0 == pytest.approx(0.0, abs=1e-9) and iq0 == pytest.approx(7.0, rel=1e-9)
    # zero current -> zero torque
    assert mtpa_operating_point(0.2, 1e-3, 5e-4, 0.0, 4) == (0.0, 0.0, 0.0, 0.0)


def test_continuous_loop_mtpa_ipm_control_gate():
    lm, Ld, Lq, I, p = 0.10, 4.0e-3, 9.0e-3, 20.0, 3
    g, id_, iq, T = mtpa_operating_point(lm, Ld, Lq, I, p)
    T_pure_q = dq_torque(lm, Ld, Lq, 0.0, I, p)

    assert math.degrees(g) == pytest.approx(30.0)
    assert id_ == pytest.approx(-10.0)
    assert iq == pytest.approx(17.320508075688775)
    assert T == pytest.approx(11.691342951089922)
    assert T_pure_q == pytest.approx(9.0)
    assert T / T_pure_q - 1.0 == pytest.approx(0.299038105676658)


if __name__ == "__main__":
    for c in CASES:
        test_mtpa_matches_numeric_argmax(*c)
    test_mtpa_limits()
    print("[OK] MTPA closed form validated vs numeric argmax + limits.")
