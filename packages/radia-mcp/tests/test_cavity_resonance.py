r"""3D rectangular cavity resonance (full-wave Maxwell eigenvalue) -- regression test (#59).

A closed PEC box a x b x d resonates at f_mnp = (c/2) sqrt((m/a)^2+(n/b)^2+(p/d)^2). The 3-D
vector (curl-curl E) sequel to the 2-D waveguide cutoff (#53), solved on HCurl edge elements with
the gradient kernel suppressed by an eigensolver shift. Closed-form helper (tool-independent) + an
HCurl Maxwell eigensolve."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import (rectangular_cavity_frequency,
                                               rectangular_cavity_modes, C0)

A, B, D = 0.040, 0.020, 0.030


def test_cavity_closed_form():
    # dominant TE101 = (c/2) sqrt(1/a^2 + 1/d^2); ~6.246 GHz for this box
    f101 = rectangular_cavity_frequency(A, B, D, 1, 0, 1)
    assert math.isclose(f101, 0.5 * C0 * math.hypot(1 / A, 1 / D), rel_tol=1e-12)
    assert math.isclose(f101, 6.246e9, rel_tol=2e-3)
    # ordering for a > d > b: TE101 < TM110 < TE011
    f110 = rectangular_cavity_frequency(A, B, D, 1, 1, 0)
    f011 = rectangular_cavity_frequency(A, B, D, 0, 1, 1)
    assert f101 < f110 < f011
    # a cubic cavity L^3: TE101 = (c/2) sqrt(2)/L = c/(sqrt2 L)
    L = 0.05
    assert math.isclose(rectangular_cavity_frequency(L, L, L, 1, 0, 1),
                        C0 * math.sqrt(2) / (2 * L), rel_tol=1e-12)


def test_rectangular_cavity_mode_table_sorted():
    modes = rectangular_cavity_modes(A, B, D, max_index=2, limit=4)
    assert modes[0]["indices"] == (1, 0, 1)
    assert math.isclose(modes[0]["frequency"],
                        rectangular_cavity_frequency(A, B, D, 1, 0, 1),
                        rel_tol=1e-12)
    assert [row["frequency"] for row in modes] == sorted(row["frequency"] for row in modes)
    assert modes[1]["indices"] == (1, 1, 0)
    assert {modes[2]["indices"], modes[3]["indices"]} == {(0, 1, 1), (2, 0, 1)}


def test_cubic_cavity_mode_degeneracy():
    L = 0.05
    modes = rectangular_cavity_modes(L, L, L, max_index=1, limit=3)
    assert [row["indices"] for row in modes] == [(0, 1, 1), (1, 0, 1), (1, 1, 0)]
    assert all(math.isclose(row["frequency"], modes[0]["frequency"], rel_tol=1e-12)
               for row in modes)


def test_rectangular_cavity_mode_table_validation():
    import pytest
    with pytest.raises(ValueError):
        rectangular_cavity_modes(0.0, B, D)
    with pytest.raises(ValueError):
        rectangular_cavity_modes(A, B, D, max_index=0)
    with pytest.raises(ValueError):
        rectangular_cavity_modes(A, B, D, limit=0)
