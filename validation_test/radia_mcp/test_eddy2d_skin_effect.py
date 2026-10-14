"""Skin-effect AC resistance of an isolated round wire (frequency-domain 2D eddy).

Gates the FEM solve in radia_mcp.radia_ngsolve.eddy2d.ac_resistance_round_wire
against the exact Kelvin-function ratio Rac/Rdc(xi). xi = a*sqrt(omega*mu*sigma)
is the radius in (sqrt2 x) skin depths; the ratio is 1 at DC and ~ xi/(2 sqrt2)+1/4
at high frequency. Both legs are independent (FEM J_z integral vs ber/bei series).
"""
import math
import os
import sys

import pytest

pytest.importorskip("scipy", reason="Kelvin-function reference requires scipy")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.eddy2d import (
    kelvin_ac_resistance_ratio, ac_resistance_round_wire)


@pytest.fixture(scope="module")
def fem_deps():
    pytest.importorskip("numpy", reason="FEM skin-effect path requires numpy")
    pytest.importorskip("ngsolve", reason="FEM skin-effect path requires ngsolve")
    pytest.importorskip("netgen.geom2d", reason="FEM skin-effect path requires netgen")


def test_kelvin_dc_limit():
    # DC: no skin effect, Rac/Rdc -> 1
    assert kelvin_ac_resistance_ratio(0.05) == pytest.approx(1.0, abs=2e-4)


def test_kelvin_high_frequency_asymptote():
    # high xi: Rac/Rdc ~ xi/(2 sqrt2) + 1/4
    xi = 20.0
    approx = xi / (2.0 * math.sqrt(2.0)) + 0.25
    assert kelvin_ac_resistance_ratio(xi) == pytest.approx(approx, rel=0.02)


def test_kelvin_monotonic_increasing():
    vals = [kelvin_ac_resistance_ratio(xi) for xi in (0.5, 1.0, 2.0, 3.0, 4.0, 6.0)]
    assert all(b > a for a, b in zip(vals, vals[1:]))
    assert vals[0] == pytest.approx(1.0, abs=5e-3)   # near-DC


@pytest.mark.parametrize("xi", [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0])
def test_fem_matches_kelvin(fem_deps, xi):
    res = ac_resistance_round_wire(xi=xi, order=4)
    exact = kelvin_ac_resistance_ratio(xi)
    assert res["xi"] == pytest.approx(xi, rel=1e-12)
    assert res["rac_rdc"] == pytest.approx(exact, rel=5e-3)


def test_fem_wire_area(fem_deps):
    res = ac_resistance_round_wire(xi=2.0, a=1.0, order=4)
    assert res["area"] == pytest.approx(math.pi, rel=1e-3)


def test_freq_and_xi_paths_agree(fem_deps):
    # specifying physical freq must reproduce the same xi/ratio as the xi path
    a, sigma, mu = 1.0e-3, 5.8e7, 4e-7 * math.pi
    freq = 2.0e4
    res = ac_resistance_round_wire(a=a, sigma=sigma, freq=freq, order=4)
    xi_expected = a * math.sqrt(2 * math.pi * freq * mu * sigma)
    assert res["xi"] == pytest.approx(xi_expected, rel=1e-9)
    assert res["rac_rdc"] == pytest.approx(
        kelvin_ac_resistance_ratio(xi_expected), rel=1e-2)
