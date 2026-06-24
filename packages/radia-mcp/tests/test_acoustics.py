"""Closed-form acoustic radiation helpers -- fast regression checks."""

import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.acoustics import (
    helmholtz_green_3d,
    helmholtz_green_low_frequency_series,
    pulsating_sphere_radiation,
)


def test_pulsating_sphere_impedance_and_power_conservation():
    a = 0.1
    c = 343.0
    rho = 1.2041
    ka = 1.25
    f = ka * c / (2.0 * math.pi * a)
    v0 = 0.02
    out = pulsating_sphere_radiation(a, f, v0, rho=rho, c=c, sample_radius=12.0 * a)

    denom = 1.0 + ka * ka
    assert out["ka"] == pytest.approx(ka)
    assert out["specific_resistance"] / (rho * c) == pytest.approx(ka * ka / denom)
    assert out["specific_reactance"] / (rho * c) == pytest.approx(ka / denom)

    expected_power = 0.5 * 4.0 * math.pi * a * a * out["specific_resistance"] * v0 * v0
    assert out["radiated_power"] == pytest.approx(expected_power)
    assert out["sample_power"] == pytest.approx(out["radiated_power"], rel=1e-14)


def test_pulsating_sphere_scaling_and_far_pressure():
    a = 0.05
    c = 343.0
    v0 = 0.01
    low1 = pulsating_sphere_radiation(a, 0.05 * c / (2.0 * math.pi * a), v0)
    low2 = pulsating_sphere_radiation(a, 0.10 * c / (2.0 * math.pi * a), v0)
    assert low2["specific_resistance"] / low1["specific_resistance"] == pytest.approx(4.0, rel=0.01)
    assert low2["specific_reactance"] / low1["specific_reactance"] == pytest.approx(2.0, rel=0.01)

    r10 = pulsating_sphere_radiation(a, 500.0, v0, c=c, sample_radius=10.0 * a)
    r20 = pulsating_sphere_radiation(a, 500.0, v0, c=c, sample_radius=20.0 * a)
    assert abs(r20["sample_pressure"]) / abs(r10["sample_pressure"]) == pytest.approx(0.5)

    high = pulsating_sphere_radiation(a, 10.0 * c / (2.0 * math.pi * a), v0)
    assert high["radiation_efficiency"] == pytest.approx(100.0 / 101.0)


def test_pulsating_sphere_validation():
    with pytest.raises(ValueError):
        pulsating_sphere_radiation(0.0, 100.0, 1.0)
    with pytest.raises(ValueError):
        pulsating_sphere_radiation(0.1, 0.0, 1.0)
    with pytest.raises(ValueError):
        pulsating_sphere_radiation(0.1, 100.0, 1.0, sample_radius=0.05)


def test_helmholtz_green_low_frequency_series_terms():
    r = 2.0
    k = 1.0e-3
    out = helmholtz_green_low_frequency_series(r, k, order=4)

    assert out["laplace_term"].real == pytest.approx(1.0 / (4.0 * math.pi * r))
    assert out["laplace_term"].imag == pytest.approx(0.0)
    assert out["terms"][1].real == pytest.approx(0.0)
    assert out["terms"][1].imag == pytest.approx(-k / (4.0 * math.pi))
    assert out["terms"][2].real == pytest.approx(-(k * k) * r / (8.0 * math.pi))
    assert out["terms"][2].imag == pytest.approx(0.0)
    assert out["abs_error"] < 2.0e-17
    assert out["approx"] == pytest.approx(helmholtz_green_3d(r, k))


def test_helmholtz_green_series_convergence_and_validation():
    r = 0.3
    for kr in (1.0e-4, 1.0e-2, 0.1, 0.5):
        k = kr / r
        err2 = helmholtz_green_low_frequency_series(r, k, order=2)["abs_error"]
        err6 = helmholtz_green_low_frequency_series(r, k, order=6)["abs_error"]
        assert err6 < err2
        assert err6 < 1.0e-6

    with pytest.raises(ValueError):
        helmholtz_green_3d(0.0, 1.0)
    with pytest.raises(ValueError):
        helmholtz_green_low_frequency_series(1.0, 1.0, order=-1)
