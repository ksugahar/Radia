r"""Rectangular-cavity quality factor (#68) -- test.

Q_TE101 = (k a d)^3 b eta / (2 pi^2 R_s (2a^3 b + 2 b d^3 + a^3 d + a d^3)), the wall-loss-limited
unloaded Q. Validated against the numerically-integrated TE101 field energy and wall loss, the
Q ~ sqrt(sigma) (1/R_s) scaling, and the surface-resistance / skin-depth relation. Tool-independent."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import (rectangular_cavity_frequency, rectangular_cavity_q,
                                               surface_resistance, skin_depth, MU0, C0)

EPS0 = 1.0 / (MU0 * C0 ** 2)
CAVS = [(0.04, 0.02, 0.03), (0.0229, 0.0102, 0.0229), (0.05, 0.03, 0.05)]


def _q_from_fields(a, b, d, sigma, n=200):
    np = pytest.importorskip("numpy")
    f = rectangular_cavity_frequency(a, b, d, 1, 0, 1)
    w = 2 * math.pi * f
    Rs = surface_resistance(f, sigma)
    Camp = 1.0 / (w * MU0)
    x = np.linspace(0, a, n); z = np.linspace(0, d, n)
    X, Z = np.meshgrid(x, z, indexing="ij")
    U = (EPS0 / 2.0) * np.trapezoid(np.trapezoid((np.sin(math.pi * X / a) * np.sin(math.pi * Z / d)) ** 2,
                                                  z, axis=1), x) * b
    Hx = lambda xx, zz: Camp * (math.pi / d) * np.sin(math.pi * xx / a) * np.cos(math.pi * zz / d)
    Hz = lambda xx, zz: Camp * (math.pi / a) * np.cos(math.pi * xx / a) * np.sin(math.pi * zz / d)
    surf = (np.trapezoid(Hz(0.0, z) ** 2, z) + np.trapezoid(Hz(a, z) ** 2, z)
            + np.trapezoid(Hx(x, 0.0) ** 2, x) + np.trapezoid(Hx(x, d) ** 2, x)) * b \
        + 2 * np.trapezoid(np.trapezoid(Hx(X, Z) ** 2 + Hz(X, Z) ** 2, z, axis=1), x)
    return w * U / ((Rs / 2.0) * surf)


def test_q_closed_form_matches_field_integration():
    pytest.importorskip("numpy")
    sigma = 5.8e7
    for (a, b, d) in CAVS:
        assert math.isclose(rectangular_cavity_q(a, b, d, sigma), _q_from_fields(a, b, d, sigma),
                            rel_tol=1e-3)


def test_q_scales_as_sqrt_sigma():
    a, b, d = 0.04, 0.02, 0.03
    # Q ~ 1/R_s ~ sqrt(sigma): 4x conductivity -> 2x Q
    assert math.isclose(rectangular_cavity_q(a, b, d, 4 * 5.8e7),
                        2.0 * rectangular_cavity_q(a, b, d, 5.8e7), rel_tol=1e-9)
    # realistic copper X-band cavity: Q in the thousands-to-tens-of-thousands
    assert 3e3 < rectangular_cavity_q(a, b, d, 5.8e7) < 5e4


def test_surface_resistance_and_skin_depth():
    f, sigma = 1e10, 5.8e7
    Rs, delta = surface_resistance(f, sigma), skin_depth(f, sigma)
    assert math.isclose(Rs, 1.0 / (sigma * delta), rel_tol=1e-12)       # R_s = 1/(sigma delta)
    assert math.isclose(delta, math.sqrt(2.0 / (2 * math.pi * f * MU0 * sigma)), rel_tol=1e-12)
    assert 0.5e-6 < delta < 0.8e-6                                       # copper @ 10 GHz ~ 0.66 um
    # R_s grows as sqrt(frequency)
    assert math.isclose(surface_resistance(4 * f, sigma), 2.0 * Rs, rel_tol=1e-12)


def test_te101_frequency_anchor():
    # the 40x20x30 cavity is the #59 TE101 mode; the frequency does NOT depend on b (n=0)
    a, b, d = 0.04, 0.02, 0.03
    f = rectangular_cavity_frequency(a, b, d, 1, 0, 1)
    assert math.isclose(f, 0.5 * C0 * math.hypot(1 / a, 1 / d), rel_tol=1e-12)
    assert math.isclose(f, rectangular_cavity_frequency(a, 2 * b, d, 1, 0, 1), rel_tol=1e-12)
    # a taller cavity raises Q (more volume per end-wall loss) but sub-linearly (side walls also grow)
    r = rectangular_cavity_q(a, 2 * b, d, 5.8e7) / rectangular_cavity_q(a, b, d, 5.8e7)
    assert 1.0 < r < 2.0


if __name__ == "__main__":
    test_q_closed_form_matches_field_integration()
    test_q_scales_as_sqrt_sigma()
    test_surface_resistance_and_skin_depth()
    test_te101_frequency_anchor()
    print("[OK] cavity Q: closed form == field-integrated omega U/P_wall; Q ~ sqrt(sigma); R_s=1/(sigma delta).")
