r"""Winding MMF space harmonics & rotating-field synthesis (#72) -- test.

F_phi,n = (4/pi)(k_w,n/n)(N_ph I/2p); the balanced 3-phase combination gives rotating waves of
amplitude (m/2)F_phi,n with the triplen cancelled and n=6k+-1 surviving (forward 6k+1, backward
6k-1). Validated against the direct rotating-field trig sum and the (2/pi) fundamental of #69."""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (winding_factor, single_phase_mmf_harmonic,
                                           rotating_mmf_amplitude, mmf_harmonic_direction)


def _mmf_sum(n, theta, t, m=3):
    return sum(math.cos(t - k * 2 * math.pi / m) * math.cos(n * (theta - k * 2 * math.pi / m))
               for k in range(m))


def test_rotating_field_theorem():
    # the m-phase sum equals (m/2) cos(n theta - dir t) for surviving n, 0 for triplen
    for n in (1, 2, 3, 5, 7, 9, 11, 13):
        d = mmf_harmonic_direction(n, 3)
        for (th, t) in [(0.4, 1.1), (1.7, 0.3), (2.5, 2.0)]:
            pred = rotating_mmf_amplitude(n, 3, 1.0) * (math.cos(n * th - d * t) if d else 0.0)
            assert math.isclose(_mmf_sum(n, th, t), pred, abs_tol=1e-12)


def test_directions_and_cancellation():
    # 3-phase: triplen cancel, 6k+1 forward, 6k-1 backward, even cancel
    assert mmf_harmonic_direction(1, 3) == 1 and mmf_harmonic_direction(7, 3) == 1
    assert mmf_harmonic_direction(13, 3) == 1
    assert mmf_harmonic_direction(5, 3) == -1 and mmf_harmonic_direction(11, 3) == -1
    assert mmf_harmonic_direction(3, 3) == 0 and mmf_harmonic_direction(9, 3) == 0
    assert mmf_harmonic_direction(2, 3) == -1      # theorem: n=2 backward (even absent in real windings)
    # rotating amplitude: triplen -> 0, else (m/2) * single phase
    assert rotating_mmf_amplitude(3, 3, 5.0) == 0.0
    assert math.isclose(rotating_mmf_amplitude(5, 3, 5.0), 7.5, rel_tol=1e-12)
    assert math.isclose(rotating_mmf_amplitude(1, 5, 4.0), 10.0, rel_tol=1e-12)   # 5-phase -> 5/2


def test_amplitude_formula_and_69_link():
    kw1, N_ph, I, p = 0.933, 120, 10.0, 2
    F1 = single_phase_mmf_harmonic(1, kw1, N_ph, I, p)
    # F_phi,1 = (4/pi)(kw1/1)(N I/2p) = (2/pi) kw1 N I / p
    assert math.isclose(F1, (2.0 / math.pi) * kw1 * N_ph * I / p, rel_tol=1e-12)
    # harmonic envelope ~ kw_n/(n): the 5th is well below the fundamental
    F5 = single_phase_mmf_harmonic(5, 0.05, N_ph, I, p)
    assert F5 < 0.1 * F1
    # the rotating fundamental is 3/2 the per-phase (the (m/2) of L_m=(3/2)L_mu, #69)
    assert math.isclose(rotating_mmf_amplitude(1, 3, F1), 1.5 * F1, rel_tol=1e-12)


if __name__ == "__main__":
    test_rotating_field_theorem()
    test_directions_and_cancellation()
    test_amplitude_formula_and_69_link()
    print("[OK] MMF harmonics: rotating-field theorem; triplen cancel / 6k+-1 survive; F1=(2/pi)kw1 N I/p.")
