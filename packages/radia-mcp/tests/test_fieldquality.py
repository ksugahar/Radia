"""Unit tests for radia_mcp.radia_ngsolve.fieldquality -- the multipole / field-
quality post-processing. The core extraction is pure numpy/cmath (NO NGSolve),
so these run everywhere and lock the CERN/European convention and the exact
line-current reference. The FEM cross-check (a real quadrupole solve) lives in
validation/force/validate_force_xval.py::validate_multipole_quadrupole_fem.
"""
import math
import cmath
import pytest

from radia_mcp.radia_ngsolve.fieldquality import (
    multipole_coefficients, line_current_multipoles, superpose_multipoles,
    field_errors, MU0)


def _line_field(I, r0, phi0, R):
    """Closed-form planar field of one line current: By + i Bx = (mu0 I)/(2pi)/(z - z0).
    Returns a callable phi -> (Bx, By) on the circle |z| = R."""
    z0 = cmath.rect(r0, phi0)

    def B(phi):
        z = cmath.rect(R, phi)
        f = (MU0 * I) / (2.0 * math.pi) / (z - z0)   # = By + i Bx
        return (f.imag, f.real)                       # (Bx, By)
    return B


def test_extraction_matches_line_current_closed_form():
    """The DFT extractor must reproduce line_current_multipoles to ~machine eps
    when fed that current's exact field (validates convention + DFT bins)."""
    I, r0, phi0, R = 137.0, 0.06, 0.7, 0.025
    num = multipole_coefficients(_line_field(I, r0, phi0, R), R, n_max=8)
    ana = line_current_multipoles(I, r0, phi0, R, n_max=8)
    scale = max(ana["C"][1:])
    for n in range(1, 9):
        assert abs(num["b"][n] - ana["b"][n]) < 1e-9 * scale
        assert abs(num["a"][n] - ana["a"][n]) < 1e-9 * scale


def test_line_current_on_axis_is_purely_normal():
    """A current on the +x axis (phi0=0) gives purely NORMAL multipoles (a_n=0),
    decaying as b_n ~ R^{n-1}/r0^n."""
    I, r0, R = 100.0, 0.05, 0.02
    d = line_current_multipoles(I, r0, 0.0, R, n_max=6)
    for n in range(1, 7):
        assert abs(d["a"][n]) < 1e-15 * d["C"][1]
        expect = MU0 * I / (2 * math.pi) * R ** (n - 1) / r0 ** n
        assert d["b"][n] == pytest.approx(-expect, rel=1e-12)


def test_normal_quadrupole_symmetry():
    """4 wires (+,-,+,-) at 0/90/180/270 deg -> NORMAL quadrupole: main = n2 with
    a2 ~ 0; the only non-vanishing harmonics are the ALLOWED n = 2, 6, 10; the
    12-pole sits at exactly (R_ref/r0)^4 of the main."""
    I, r0, R = 1000.0, 0.05, 0.02
    refs = [line_current_multipoles((+1 if k % 2 == 0 else -1) * I, r0,
                                    math.radians(90 * k), R, 10) for k in range(4)]
    q = superpose_multipoles(*refs)
    fe = field_errors(q)
    assert fe["main_n"] == 2
    assert abs(q["a"][2]) < 1e-9 * q["C"][2]              # normal, not skew
    for n in (1, 3, 4, 5, 7, 8, 9):                        # forbidden
        assert fe["rel_C"][n] < 1e-6
    assert fe["rel_C"][6] == pytest.approx((R / r0) ** 4 * 1e4, rel=1e-9)   # 256 units
    assert fe["rel_C"][10] == pytest.approx((R / r0) ** 8 * 1e4, rel=1e-9)


def test_skew_quadrupole_when_rotated_45deg():
    """The SAME quad rotated 45 deg (wires at 45/135/225/315) is a SKEW quad:
    a2 dominant, b2 ~ 0 -- the convention's normal/skew split tracks orientation."""
    I, r0, R = 1000.0, 0.05, 0.02
    refs = [line_current_multipoles((+1 if k % 2 == 0 else -1) * I, r0,
                                    math.radians(45 + 90 * k), R, 6) for k in range(4)]
    q = superpose_multipoles(*refs)
    assert field_errors(q)["main_n"] == 2
    assert abs(q["b"][2]) < 1e-9 * q["C"][2]              # skew, not normal


def test_field_errors_main_is_unit_scale():
    d = line_current_multipoles(50.0, 0.04, 0.3, 0.015, n_max=6)
    fe = field_errors(d, units=1e4)
    assert fe["rel_C"][fe["main_n"]] == pytest.approx(1e4)


def test_n_samples_guard_rejects_aliasing():
    with pytest.raises(ValueError):
        multipole_coefficients(_line_field(1.0, 0.05, 0.0, 0.02), 0.02,
                               n_max=10, n_samples=16)   # 16 <= 2*10


def test_superpose_rejects_mismatched_reference():
    a = line_current_multipoles(1.0, 0.05, 0.0, 0.02, n_max=6)
    b = line_current_multipoles(1.0, 0.05, 0.0, 0.03, n_max=6)   # different R_ref
    with pytest.raises(ValueError):
        superpose_multipoles(a, b)
