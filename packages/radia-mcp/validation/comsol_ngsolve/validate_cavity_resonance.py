r"""3D rectangular cavity resonance (full-wave Maxwell eigenvalue) -- regression test (#59).

A closed PEC box a x b x d resonates at f_mnp = (c/2) sqrt((m/a)^2+(n/b)^2+(p/d)^2). The 3-D
vector (curl-curl E) sequel to the 2-D waveguide cutoff (#53), solved on HCurl edge elements with
the gradient kernel suppressed by an eigensolver shift. Closed-form helper (tool-independent) + an
HCurl Maxwell eigensolve."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import rectangular_cavity_frequency, C0

A, B, D = 0.040, 0.020, 0.030


def validate_cavity_closed_form():
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


def validate_cavity_maxwell_fe():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, TaskManager
    from netgen.occ import Box, OCCGeometry, Pnt
    from radia_mcp.radia_ngsolve.waveguide import maxwell_cavity_modes_3d

    box = Box(Pnt(0, 0, 0), Pnt(A, B, D))
    for f in box.faces:
        f.name = "pec"
    mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.006))

    shift = (math.pi / A) ** 2 + (math.pi / D) ** 2
    with TaskManager():
        fe = maxwell_cavity_modes_3d(mesh, 4, shift)

    # the dominant FE mode is TE101, to FE accuracy; gradient kernel suppressed (no near-0 modes)
    assert math.isclose(fe[0], rectangular_cavity_frequency(A, B, D, 1, 0, 1), rel_tol=3e-3)
    assert fe[0] > 1e9                                       # not a spurious zero mode
    # the next physical modes match the box spectrum (TM110, TE011)
    assert math.isclose(fe[1], rectangular_cavity_frequency(A, B, D, 1, 1, 0), rel_tol=3e-3)
    # ascending and all physical
    assert all(b >= a for a, b in zip(fe, fe[1:]))


if __name__ == "__main__":
    validate_cavity_closed_form()
    validate_cavity_maxwell_fe()
    print("[OK] 3D cavity: HCurl Maxwell eigenvalue == exact box spectrum (TE101 dominant).")
