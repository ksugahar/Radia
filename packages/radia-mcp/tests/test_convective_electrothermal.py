r"""Convective electro-thermal (Joule heat + Newton cooling) -- regression test (#44).

A slab with a uniform Joule source q, convectively cooled (film h) on both faces, has the
centre rise dT = qL^2/(8k) + qL/(2h) = conduction parabola + convective film. radia's Robin
heat solver reproduces it; the closed form is the tool-independent gate. A reference
commercial FE code matches the same closed form to 0.000% live (recorded internally). As
h -> inf the film vanishes and the fixed-T parabola (#19/#41) is recovered."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.multiphysics import (convective_slab_peak_dT,
                                                  fixed_temperature_slab_peak_dT,
                                                  fixed_temperature_slab_temperature_rise)

L, WD = 0.02, 0.10
Q, K = 5.0e5, 2.0


def test_peak_dt_formula():
    # dT = conduction qL^2/8k + film qL/2h
    assert math.isclose(convective_slab_peak_dT(Q, L, K, 500.0),
                        Q*L*L/(8*K) + Q*L/(2*500.0), rel_tol=1e-12)
    # for these values: conduction 12.5 K + film 10 K = 22.5 K
    assert math.isclose(convective_slab_peak_dT(Q, L, K, 500.0), 22.5, rel_tol=1e-12)
    # h -> inf : film vanishes -> the fixed-T parabola qL^2/8k = 12.5 K
    assert math.isclose(convective_slab_peak_dT(Q, L, K, 1e12), Q*L*L/(8*K), rel_tol=1e-6)
    # smaller h (worse cooling) -> larger film -> hotter
    assert convective_slab_peak_dT(Q, L, K, 250.0) > convective_slab_peak_dT(Q, L, K, 500.0)


def test_fixed_temperature_slab_parabola():
    assert math.isclose(fixed_temperature_slab_peak_dT(Q, L, K), Q * L * L / (8 * K),
                        rel_tol=1e-12)
    assert math.isclose(fixed_temperature_slab_peak_dT(Q, L, K), 12.5, rel_tol=1e-12)
    assert math.isclose(fixed_temperature_slab_temperature_rise(Q, L, K, 0.0), 0.0,
                        abs_tol=1e-12)
    assert math.isclose(fixed_temperature_slab_temperature_rise(Q, L, K, L), 0.0,
                        abs_tol=1e-12)
    assert math.isclose(fixed_temperature_slab_temperature_rise(Q, L, K, L / 2),
                        fixed_temperature_slab_peak_dT(Q, L, K), rel_tol=1e-12)
    xs = [0.0, L / 4, L / 2, 3 * L / 4, L]
    profile = fixed_temperature_slab_temperature_rise(Q, L, K, xs)
    assert profile[0] == pytest.approx(0.0)
    assert profile[-1] == pytest.approx(0.0)
    assert profile[1] == pytest.approx(profile[3])
    assert profile[2] == pytest.approx(12.5)
