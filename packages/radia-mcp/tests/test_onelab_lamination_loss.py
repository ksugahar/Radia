r"""ONELAB 2D-correction #6 (lamination iron eddy = order-0 homogenisation) PARITY.

ONELAB/GetDP models laminated iron as non-conducting in the field solve and recovers
the classical eddy loss a-posteriori with the homogenisation-order-0 term quoted in
im_3kW_data.geo:

    P = (d^2 / 12) * sigma * <(dB/dt)^2>

This locks that radia-ngsolve's ``coreloss`` reproduces that term EXACTLY (for
sinusoidal flux), i.e. radia already covers ONELAB correction #6 -- and exceeds it,
since the ONELAB bundle has ONLY the classical term (no hysteresis/excess, off by
default) whereas radia carries the full Bertotti separation + per-harmonic path.
Pure arithmetic, no FE solve.  See radia_mcp.motor.onelab_knowledge `twod_corrections`.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.coreloss import (bertotti_loss_density,
                                              classical_eddy_loss_density,
                                              lamination_kc)


def _onelab_classical_eddy(Bpk, f, sigma, d):
    """ONELAB order-0 homogenisation term P = (d^2/12) sigma <(dB/dt)^2>, computed
    independently. For B = Bpk sin(2 pi f t): <(dB/dt)^2> = (1/2)(2 pi f Bpk)^2."""
    mean_sq_dBdt = 0.5 * (2.0 * math.pi * f * Bpk) ** 2     # = 2 pi^2 f^2 Bpk^2
    return (d ** 2 / 12.0) * sigma * mean_sq_dBdt


def test_classical_eddy_equals_onelab_d2_over_12_term():
    sigma, d, Bpk, f = 2.0e6, 0.5e-3, 1.5, 50.0          # VH800-class lamina, 0.5 mm
    onelab = _onelab_classical_eddy(Bpk, f, sigma, d)
    radia = classical_eddy_loss_density(Bpk, f, sigma, d)
    assert math.isclose(radia, onelab, rel_tol=1e-12), (radia, onelab)
    # and via the Bertotti coefficient k_c = pi^2 sigma d^2 / 6
    assert math.isclose(lamination_kc(sigma, d) * (f * Bpk) ** 2, onelab, rel_tol=1e-12)


def test_radia_exceeds_onelab_with_hyst_and_excess():
    """ONELAB bundle = classical only; radia adds hysteresis + excess on top."""
    sigma, d, Bpk, f = 2.0e6, 0.5e-3, 1.0, 50.0
    P_class = classical_eddy_loss_density(Bpk, f, sigma, d)
    k_c = lamination_kc(sigma, d)
    P_full = bertotti_loss_density(Bpk, f, k_h=120.0, k_c=k_c, k_e=0.8)
    assert P_full > P_class                               # hyst + excess strictly add
    # the classical component of the full model is exactly the ONELAB term
    assert math.isclose(bertotti_loss_density(Bpk, f, k_h=0.0, k_c=k_c, k_e=0.0),
                        _onelab_classical_eddy(Bpk, f, sigma, d), rel_tol=1e-12)


def main():
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
    print("[OK] ONELAB lamination-loss correction #6 covered by radia coreloss "
          "(classical-eddy parity exact; radia adds hysteresis + excess)")


if __name__ == "__main__":
    main()
