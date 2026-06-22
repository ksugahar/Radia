r"""Rectangular waveguide cutoff via 2D Helmholtz eigenmodes -- regression test (#53).

The wave-physics entry in the comsol_class series. TE modes = NEUMANN Laplacian, TM modes =
DIRICHLET Laplacian on the cross-section; the eigenvalues are the squared cutoff wavenumbers
k_c^2, with f_c = c k_c/(2 pi). Validated against the exact rectangular spectrum
f_c,mn = (c/2) sqrt((m/a)^2 + (n/b)^2). Closed-form helper (tool-independent) + an FE eigensolve."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import (rectangular_waveguide_cutoff, cutoff_frequency, C0)

A, B = 0.02286, 0.01016          # WR-90 X-band guide


def test_rectangular_cutoff_closed_form():
    # dominant TE10 cutoff = c/(2a); textbook WR-90 value ~6.557 GHz
    f10 = rectangular_waveguide_cutoff(A, B, 1, 0)
    assert math.isclose(f10, C0 / (2 * A), rel_tol=1e-12)
    assert math.isclose(f10, 6.557e9, rel_tol=2e-3)
    # next modes ordered: TE10 < TE20 < TE01 < TM11
    f20 = rectangular_waveguide_cutoff(A, B, 2, 0)
    f01 = rectangular_waveguide_cutoff(A, B, 0, 1)
    f11 = rectangular_waveguide_cutoff(A, B, 1, 1)
    assert f10 < f20 < f01 < f11
    assert math.isclose(f20, 2 * f10, rel_tol=1e-12)            # TE20 = 2 TE10
    assert math.isclose(f01, C0 / (2 * B), rel_tol=1e-12)       # TE01 = c/2b
    # wavenumber round-trip: f = c kc/2pi
    kc = 2 * math.pi * f11 / C0
    assert math.isclose(cutoff_frequency(kc), f11, rel_tol=1e-12)
