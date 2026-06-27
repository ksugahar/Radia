r"""Rectangular waveguide cutoff via 2D Helmholtz eigenmodes -- regression test (#53).

The wave-physics validation entry. TE modes = NEUMANN Laplacian, TM modes =
DIRICHLET Laplacian on the cross-section; the eigenvalues are the squared cutoff wavenumbers
k_c^2, with f_c = c k_c/(2 pi). Validated against the exact rectangular spectrum
f_c,mn = (c/2) sqrt((m/a)^2 + (n/b)^2). Closed-form helper (tool-independent) + an FE eigensolve."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import (rectangular_waveguide_cutoff, cutoff_frequency, C0)

A, B = 0.02286, 0.01016          # WR-90 X-band guide


def validate_rectangular_cutoff_closed_form():
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


def validate_waveguide_eigenmodes_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, TaskManager
    from netgen.occ import OCCGeometry, MoveTo
    from radia_mcp.radia_ngsolve.waveguide import helmholtz_cutoff_wavenumbers_2d

    face = MoveTo(0, 0).Rectangle(A, B).Face()
    face.edges.name = "wall"
    mesh = Mesh(OCCGeometry(face, dim=2).GenerateMesh(maxh=min(A, B) / 16))

    with TaskManager():
        te = helmholtz_cutoff_wavenumbers_2d(mesh, 3, bc="neumann")     # TE10, TE20, TE01
        tm = helmholtz_cutoff_wavenumbers_2d(mesh, 2, bc="dirichlet")   # TM11, TM21

    # TE (Neumann) spectrum matches the closed form; the constant null mode is dropped
    for kc, (m, n) in zip(te, [(1, 0), (2, 0), (0, 1)]):
        assert math.isclose(cutoff_frequency(kc), rectangular_waveguide_cutoff(A, B, m, n),
                            rel_tol=3e-3)
    # TM (Dirichlet) spectrum matches; lowest TM11 lies ABOVE every TE mode here
    for kc, (m, n) in zip(tm, [(1, 1), (2, 1)]):
        assert math.isclose(cutoff_frequency(kc), rectangular_waveguide_cutoff(A, B, m, n),
                            rel_tol=3e-3)
    assert cutoff_frequency(tm[0]) > cutoff_frequency(te[-1])           # TM11 > TE01
    # the BC distinction is the physics: dominant mode is TE10 = c/(2a), no TM below it
    assert math.isclose(cutoff_frequency(te[0]), C0 / (2 * A), rel_tol=3e-3)


if __name__ == "__main__":
    validate_rectangular_cutoff_closed_form()
    validate_waveguide_eigenmodes_fe()
    print("[OK] waveguide cutoff: TE=Neumann / TM=Dirichlet Laplacian eigenvalues == exact spectrum.")
