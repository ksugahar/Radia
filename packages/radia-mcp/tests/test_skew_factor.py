r"""Skew factor (axial skew -> EMF-harmonic & cogging reduction) -- regression test (#58).

Skewing over an electrical angle theta_sk scales the n-th harmonic by k_sk,n =
sin(n theta_sk/2)/(n theta_sk/2). The axial twin of the distribution factor (#51). A harmonic with
n theta_sk = 2 pi is nulled -> one-slot-pitch skew cancels the cogging/slot-passing ripple.
Closed-form helpers (tool-independent), validated against the phasor average integral."""
import cmath
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (skew_factor, slot_pitch_skew_angle,
                                           skewed_winding_factor, winding_factor)


def _phasor_average(n, theta_sk, m=40000):
    return abs(sum(cmath.exp(1j * n * (i + 0.5) * theta_sk / m) for i in range(m)) / m)


def test_skew_factor_closed_form():
    # k_sk == the phasor average (the definition), several (n, theta)
    for (n, th) in [(1, 0.3), (5, 0.5), (7, 0.8), (3, 1.2)]:
        assert math.isclose(abs(skew_factor(n, th)), _phasor_average(n, th), rel_tol=1e-3)
    # theta_sk -> 0: no skew, k_sk -> 1
    assert math.isclose(skew_factor(7, 1e-9), 1.0, rel_tol=1e-9)
    # fundamental cost: 0 < k_sk1 < 1 for a finite skew, and ~1 for a small skew
    assert 0.0 < skew_factor(1, 0.35) < 1.0
    assert skew_factor(1, 0.05) > 0.999
    # NULL when the skew spans one full period of harmonic n (n theta_sk = 2 pi)
    for n in (6, 11, 18):
        assert abs(skew_factor(n, 2 * math.pi / n)) < 1e-9
    # k_sk decreases (in magnitude envelope) as the harmonic rises at fixed skew
    assert abs(skew_factor(1, 0.6)) > abs(skew_factor(3, 0.6)) > abs(skew_factor(5, 0.6))


def test_one_slot_pitch_skew_cancels_cogging():
    Q, p = 36, 2                                   # 36 slots, 4 poles
    theta_slot = slot_pitch_skew_angle(Q, p)
    assert math.isclose(theta_slot, 2 * math.pi * p / Q, rel_tol=1e-12)
    # the slot-passing harmonic (order Q/p) is nulled by one-slot-pitch skew -> cogging cancelled
    assert abs(skew_factor(Q // p, theta_slot)) < 1e-9
    # fundamental barely affected (q=3 integer-slot, 20 deg slot pitch)
    assert skew_factor(1, theta_slot) > 0.99
    # skewed winding factor = winding factor * skew factor
    q = Q // (2 * p * 3)
    for n in (1, 5, 7):
        assert math.isclose(skewed_winding_factor(n, q, 360.0 * p / Q, theta_slot),
                            winding_factor(n, q, 360.0 * p / Q) * skew_factor(n, theta_slot),
                            rel_tol=1e-12)
    # skew shrinks every harmonic below its un-skewed winding factor
    for n in (5, 7, 11, 13):
        assert abs(skewed_winding_factor(n, q, 360.0 * p / Q, theta_slot)) < \
               abs(winding_factor(n, q, 360.0 * p / Q))


if __name__ == "__main__":
    test_skew_factor_closed_form()
    test_one_slot_pitch_skew_cancels_cogging()
    print("[OK] skew factor: sinc == phasor average; one-slot-pitch skew nulls the cogging harmonic.")
