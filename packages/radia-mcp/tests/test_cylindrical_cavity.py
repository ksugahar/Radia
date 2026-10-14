r"""Cylindrical metallic cavity resonances.

The p=0 cylindrical-cavity frequencies reduce to circular-waveguide cutoff frequencies; the
TM010 pillbox mode is independent of cavity length and is the accelerator-cavity anchor.
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import (
    C0,
    circular_waveguide_cutoff,
    cylindrical_cavity_frequency,
)


def test_tm010_pillbox_mode():
    pytest.importorskip("scipy")
    radius = 0.050
    f = cylindrical_cavity_frequency(radius, 0.030, "TM", 0, 1, 0)
    assert math.isclose(f, C0 * 2.4048255577 / (2.0 * math.pi * radius), rel_tol=1e-10)
    # TM010 has no axial variation, so the length does not enter.
    assert math.isclose(f, cylindrical_cavity_frequency(radius, 0.300, "TM", 0, 1, 0), rel_tol=1e-12)
    # Radial scaling: doubling radius halves the frequency.
    assert math.isclose(cylindrical_cavity_frequency(2.0 * radius, 0.030, "TM", 0, 1, 0),
                        0.5 * f, rel_tol=1e-12)


def test_p0_modes_match_circular_waveguide_cutoff():
    pytest.importorskip("scipy")
    radius = 0.0127
    length = 0.040
    for mode, m, n in (("TM", 0, 1), ("TM", 1, 1), ("TE", 1, 1), ("TE", 0, 1)):
        assert math.isclose(cylindrical_cavity_frequency(radius, length, mode, m, n, 0),
                            circular_waveguide_cutoff(radius, mode, m, n), rel_tol=1e-12)


def test_axial_index_adds_standing_wave_term():
    pytest.importorskip("scipy")
    radius = 0.030
    length = 0.080
    f0 = cylindrical_cavity_frequency(radius, length, "TM", 0, 1, 0)
    f1 = cylindrical_cavity_frequency(radius, length, "TM", 0, 1, 1)
    kr = 2.0 * math.pi * f0 / C0
    expected = C0 * math.hypot(kr, math.pi / length) / (2.0 * math.pi)
    assert f1 > f0
    assert math.isclose(f1, expected, rel_tol=1e-12)
    # Shorter cavity raises only the p>0 frequency.
    assert cylindrical_cavity_frequency(radius, 0.5 * length, "TM", 0, 1, 1) > f1
    assert math.isclose(cylindrical_cavity_frequency(radius, 0.5 * length, "TM", 0, 1, 0),
                        f0, rel_tol=1e-12)


def test_validation():
    pytest.importorskip("scipy")
    with pytest.raises(ValueError):
        cylindrical_cavity_frequency(0.0, 0.1, "TM", 0, 1, 0)
    with pytest.raises(ValueError):
        cylindrical_cavity_frequency(0.1, -0.1, "TM", 0, 1, 0)
    with pytest.raises(ValueError):
        cylindrical_cavity_frequency(0.1, 0.1, "XX", 0, 1, 0)
    with pytest.raises(ValueError):
        cylindrical_cavity_frequency(0.1, 0.1, "TM", -1, 1, 0)
    with pytest.raises(ValueError):
        cylindrical_cavity_frequency(0.1, 0.1, "TM", 0, 0, 0)
    with pytest.raises(ValueError):
        cylindrical_cavity_frequency(0.1, 0.1, "TM", 0, 1, -1)
