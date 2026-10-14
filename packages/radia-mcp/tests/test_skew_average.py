r"""Multi-slice SKEW post-processor (ONELAB #5 FE-multislice 2D->3D correction)
cross-validated against the analytic skew_factor.

`skew_average` axially-averages a rotor-angle sweep over the skew; on a pure
harmonic it must reduce the amplitude by EXACTLY skew_factor(n, skew_span) =
sinc(n*skew_span/2) -- the FE <-> analytic cross-check.  It also nulls the
slot-passing / cogging harmonic at one slot-pitch skew, and the N-slice form
converges to the continuous (spectral) one as N grows.  numpy-only, no FE solve.
"""
import math
import os
import sys

import numpy as np

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (skew_average, skew_factor,
                                           slot_pitch_skew_angle)


def _amp(y):
    """Amplitude of a (near-)zero-mean sinusoid from its RMS: A = sqrt(2 <y^2>)."""
    return math.sqrt(2.0 * float(np.mean(np.asarray(y) ** 2)))


def _sweep(n, Ns=720):
    th = np.linspace(0.0, 2.0 * np.pi, Ns, endpoint=False)
    return th, np.sin(n * th)


def test_continuous_skew_matches_skew_factor():
    for n in (1, 2, 6):
        for sk in (0.3, 0.7, 1.1):
            th, v = _sweep(n)
            out = skew_average(th, v, sk, n_slices=0)
            assert math.isclose(_amp(out), abs(skew_factor(n, sk)),
                                rel_tol=1e-7, abs_tol=1e-9), (n, sk)


def test_zero_skew_is_identity():
    th, v = _sweep(3)
    assert np.allclose(skew_average(th, v, 0.0, n_slices=0), v, atol=1e-9)


def test_one_harmonic_period_skew_nulls_it():
    n = 4
    th, v = _sweep(n)
    out = skew_average(th, v, 2.0 * np.pi / n, n_slices=0)   # sinc(pi) = 0
    assert _amp(out) < 1e-9


def test_slot_pitch_skew_kills_cogging_harmonic():
    slots, pole_pairs = 12, 2
    sk = slot_pitch_skew_angle(slots, pole_pairs)            # 2 pi p / Q
    # the slot-passing harmonic (order = slots) sees skew_factor(Q, sk)=sinc(pi p)=0
    assert abs(skew_factor(slots, sk)) < 1e-12
    th, v = _sweep(slots)
    assert _amp(skew_average(th, v, sk, n_slices=0)) < 1e-9


def test_discrete_slices_converge_to_continuous():
    n, sk = 3, 0.8
    th, v = _sweep(n)
    cont = _amp(skew_average(th, v, sk, n_slices=0))
    a2 = _amp(skew_average(th, v, sk, n_slices=2))
    a50 = _amp(skew_average(th, v, sk, n_slices=50))
    assert math.isclose(cont, abs(skew_factor(n, sk)), rel_tol=1e-7)
    assert abs(a50 - cont) < abs(a2 - cont)                 # more slices -> closer
    assert abs(a50 - cont) < 5e-3


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("[OK] skew_average (ONELAB #5 multi-slice) == analytic skew_factor cross-validated")


if __name__ == "__main__":
    main()
