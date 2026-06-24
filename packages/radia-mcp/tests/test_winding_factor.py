r"""AC machine winding factor k_w = k_d * k_p (distribution + pitch) -- regression test (#51).

Pure winding geometry -> tool-independent. Checks: k_d == the direct phasor sum of the q slot
EMFs, the concentrated (q=1 -> k_d=1) and full-pitch (beta=1 -> k_p,1=1) limits, the classic
5/6-pitch suppression of the 5th & 7th harmonics, and k_w,1 < 1 for a distributed/short-pitched
winding. Scales the back-EMF / torque constant (#34, F3) and shapes its harmonics."""
import cmath
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (winding_distribution_factor, winding_pitch_factor,
                                           winding_factor, integral_slot_winding_factor)

Q, P2, M = 36, 4, 3
QPPP = Q // (P2 * M)                 # = 3 slots/pole/phase
GAMMA = (P2 // 2) * 360.0 / Q        # = 20 deg elec


def test_kd_matches_phasor_sum():
    g = math.radians(GAMMA)
    for n in (1, 2, 3, 5, 7, 11, 13):
        kd_closed = winding_distribution_factor(n, QPPP, GAMMA)
        kd_phasor = abs(sum(cmath.exp(1j * n * i * g) for i in range(QPPP))) / QPPP
        assert abs(abs(kd_closed) - kd_phasor) < 1e-12, f"n={n}: |kd| {kd_closed} vs phasor {kd_phasor}"


def test_limits():
    # concentrated winding q=1 -> kd = 1 (all harmonics)
    for n in (1, 3, 5):
        assert abs(winding_distribution_factor(n, 1, GAMMA) - 1.0) < 1e-12
    # full pitch beta=1 -> kp1 = 1
    assert abs(winding_pitch_factor(1, 1.0) - 1.0) < 1e-12
    # fundamental winding factor < 1 for a distributed, short-pitched winding
    kw1 = winding_factor(1, QPPP, GAMMA, 5.0 / 6.0)
    assert 0.9 < kw1 < 1.0
    assert math.isclose(kw1, 0.9271, abs_tol=1e-3)


def test_short_pitch_kills_5th_7th():
    # the classic 5/6 short pitch: |kp5| = |kp7| = sin(pi/12) ~ 0.259 (down from 1)
    # (kp5 = sin(5*(5/6)*pi/2) = sin(25 pi/12) = sin(pi/12); kp7 = sin(35 pi/12) = sin(pi/12))
    beta = 5.0 / 6.0
    assert math.isclose(abs(winding_pitch_factor(5, beta)), math.sin(math.pi / 12), rel_tol=1e-9)
    assert math.isclose(abs(winding_pitch_factor(7, beta)), math.sin(math.pi / 12), rel_tol=1e-9)
    assert abs(winding_pitch_factor(5, beta)) < 0.3 and abs(winding_pitch_factor(7, beta)) < 0.3
    # full winding factor: the 5th/7th are strongly suppressed vs the fundamental
    kw1 = abs(winding_factor(1, QPPP, GAMMA, beta))
    assert abs(winding_factor(5, QPPP, GAMMA, beta)) < 0.15 * kw1
    assert abs(winding_factor(7, QPPP, GAMMA, beta)) < 0.15 * kw1


def test_integral_slot_winding_factor_from_slots_and_poles():
    res = integral_slot_winding_factor(Q, P2, harmonic=1, phases=M, pitch_fraction=5.0 / 6.0)

    assert res["pole_pairs"] == P2 // 2
    assert res["slots_per_pole_per_phase"] == QPPP
    assert res["slots_per_phase"] == Q // M
    assert res["slot_angle_elec_deg"] == pytest.approx(GAMMA)
    assert res["distribution_factor"] == pytest.approx(winding_distribution_factor(1, QPPP, GAMMA))
    assert res["pitch_factor"] == pytest.approx(winding_pitch_factor(1, 5.0 / 6.0))
    assert res["winding_factor"] == pytest.approx(winding_factor(1, QPPP, GAMMA, 5.0 / 6.0))
    assert res["winding_factor_abs"] == pytest.approx(abs(res["winding_factor"]))

    with pytest.raises(ValueError, match="integral-slot"):
        integral_slot_winding_factor(27, 8)


if __name__ == "__main__":
    test_kd_matches_phasor_sum()
    test_limits()
    test_short_pitch_kills_5th_7th()
    test_integral_slot_winding_factor_from_slots_and_poles()
    print("[OK] winding factor k_w = k_d k_p validated (phasor sum, limits, harmonic suppression).")
