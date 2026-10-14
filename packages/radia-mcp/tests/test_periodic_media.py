"""1-D periodic (Bragg) media band structure -- the photonic-bandgap / distributed-Bragg-reflector
member of the high-frequency family. Gated on the exact transfer-matrix dispersion:
a quarter-wave stack has 1/2 tr(f0) = -1/2(z+1/z) (mid stop-band) and a gap centred on f0;
no index contrast -> no gap. Pure analytic (numpy)."""
import math
import os
import sys

import pytest

np = pytest.importorskip("numpy")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.periodic_media import (
    bragg_half_trace, bragg_bandgaps, quarter_wave_stack_layers, C_LIGHT)

F0 = 1.0e9
N1, N2 = 2.0, 1.0


def test_quarter_wave_midgap_half_trace():
    # at the design f0 both layers are a quarter wave -> 1/2 tr(f0) = -1/2(z+1/z) exactly
    L1, L2 = quarter_wave_stack_layers(F0, N1, N2)
    z = N1 / N2
    ht0 = float(bragg_half_trace(N1, L1, N2, L2, F0))
    assert ht0 == pytest.approx(-0.5 * (z + 1.0 / z), abs=1e-9)
    assert abs(ht0) > 1.0                                  # f0 is inside a stop band (any contrast)


def test_quarter_wave_gap_centered_on_f0():
    L1, L2 = quarter_wave_stack_layers(F0, N1, N2)
    gaps = bragg_bandgaps(N1, L1, N2, L2, 0.05e9, 1.95e9)
    assert gaps, "expected a stop band around f0"
    g = min(gaps, key=lambda ab: abs(0.5 * (ab[0] + ab[1]) - F0))
    center = 0.5 * (g[0] + g[1])
    assert center == pytest.approx(F0, rel=1e-3)           # first gap symmetric about f0
    assert g[0] < F0 < g[1]


def test_no_contrast_no_gap():
    # n1 == n2 is a homogeneous medium -> |1/2 tr| <= 1 everywhere, no stop band
    L1, L2 = quarter_wave_stack_layers(F0, 1.5, 1.5)
    f = np.linspace(0.05e9, 3.0e9, 40000)
    ht = bragg_half_trace(1.5, L1, 1.5, L2, f)
    assert np.max(np.abs(ht)) <= 1.0 + 1e-9
    assert bragg_bandgaps(1.5, L1, 1.5, L2, 0.05e9, 3.0e9) == []


def test_gap_widens_with_contrast():
    # higher index contrast -> wider mid-gap |1/2 tr(f0)| = 1/2(z+1/z) grows with z
    def midgap(n1):
        L1, L2 = quarter_wave_stack_layers(F0, n1, 1.0)
        return abs(float(bragg_half_trace(n1, L1, 1.0, L2, F0)))
    assert midgap(3.0) > midgap(2.0) > midgap(1.2) > 1.0


def test_homogeneous_dispersion_recovers_plane_wave():
    # n1=n2=n: 1/2 tr = cos(omega n (L1+L2)/c) -- the ordinary plane-wave phase
    n = 1.4
    L1, L2 = 0.013, 0.021
    f = 0.7e9
    ht = float(bragg_half_trace(n, L1, n, L2, f))
    phase = 2 * math.pi * f * n * (L1 + L2) / C_LIGHT
    assert ht == pytest.approx(math.cos(phase), abs=1e-9)
