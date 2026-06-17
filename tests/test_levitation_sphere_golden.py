"""Golden test: isotropic AC levitation force on a conducting sphere.

Locks the verified result of
examples/CLN/scripts/team28_levitation/levitation_sphere_force.py:
  - the analytic sphere response G(x) hits its physical limits
    (DC Re[G]->0, high-freq Re[G]->-1/2),
  - Re[G] < 0 across the band (a lift at every frequency),
  - the CLN/Cauer reduction of the sphere modal Foster ladder converges,
  - the levitation-force coefficient matches the perfect-conductor
    analytic limit (pi a^3 / 2 mu0) |grad B0^2| (sign + normalization).

Pure numpy -- no NGSolve / Cubit, so this always runs in CI.
"""
import math
import os
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
# CLN was absorbed into radia.levitation: examples/CLN/scripts/team28_levitation ->
# examples/levitation/{sphere,ellipsoid}. Add both so sphere + ellipsoid modules resolve.
_LEV = os.path.join(_HERE, "..", "examples", "levitation")
sys.path.insert(0, os.path.join(_LEV, "sphere"))
sys.path.insert(0, os.path.join(_LEV, "ellipsoid"))

import levitation_sphere_force as L  # noqa: E402

TWO_PI = 2 * math.pi


def test_G_limits():
    g_lo = L.G_exact(TWO_PI * 1e-2)
    g_hi = L.G_exact(TWO_PI * 1e9)
    # DC: no eddy response, Re[G] -> 0 from below
    assert abs(g_lo.real) < 1e-3 and g_lo.real <= 0.0
    # high freq: perfect-conductor flux exclusion, Re[G] -> -1/2
    assert abs(g_hi.real + 0.5) < 1e-3


def test_ReG_negative_everywhere_is_lift():
    for k in range(0, 17):           # 1 Hz ... 1e8 Hz
        f = 10 ** (k * 0.5)
        assert L.G_exact(TWO_PI * f).real <= 1e-9


def test_closed_form_matches_foster_modal_sum():
    for f in (50, 500, 5e3, 5e4, 5e5):
        w = TWO_PI * f
        assert abs(L.G_exact(w) - L.G_modal(w)) < 1e-3


def test_cln_cauer_reduction_converges():
    w = TWO_PI * 5e3
    tgt = L.G_modal(w)
    err = [abs(L.G_cln(w, m) - tgt) / abs(tgt) for m in (1, 2, 3, 4, 6, 8)]
    # stage 1 is the eddy-free DC response -> large error
    assert err[0] > 0.1
    # by stage 4 within 0.1%, by stage 6 essentially exact
    assert err[3] < 1e-3, f"stage 4 rel err {err[3]*100:.4f}% (expect <0.1%)"
    assert err[4] < 1e-5, f"stage 6 rel err {err[4]*100:.6f}% (expect <1e-3%)"
    assert err[4] <= err[3]


def test_force_coefficient_pinned_to_perfect_conductor():
    a, mu0 = L.a, L.mu0
    B0, gradB2 = 0.1, -0.2
    F_pc = (math.pi * a**3 / (2 * mu0)) * abs(gradB2)
    # high frequency -> the induced-dipole force must approach the
    # independent perfect-conductor energy-method coefficient
    lift_hi = abs(L.levitation_force(TWO_PI * 1e9, B0, gradB2))
    assert abs(lift_hi / F_pc - 1.0) < 1e-3, (
        f"high-freq lift {lift_hi*1e3:.4f} mN vs perfect-conductor "
        f"{F_pc*1e3:.4f} mN")


def test_force_monotone_rises_with_frequency():
    a, mu0 = L.a, L.mu0
    B0, gradB2 = 0.1, -0.2
    lifts = [abs(L.levitation_force(TWO_PI * f, B0, gradB2))
             for f in (50, 500, 5e3, 5e4, 5e5, 5e6)]
    # the lift rises monotonically from ~0 (DC) to the saturation value
    assert all(lifts[i] < lifts[i + 1] for i in range(len(lifts) - 1))
    F_pc = (math.pi * a**3 / (2 * mu0)) * abs(gradB2)
    assert lifts[0] < 0.01 * F_pc        # 50 Hz: deep in the no-eddy regime
    assert lifts[-1] > 0.98 * F_pc       # 5 MHz: near perfect-conductor
