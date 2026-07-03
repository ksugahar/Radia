"""Golden test: HDiv-VIM demag factor on a SMOOTH NON-SPHERICAL body (production #3).

The sphere/cube tests fix the demag factor at the isotropic 1/3.  A prolate spheroid has a non-trivial
ANALYTIC demag factor != 1/3 (e.g. a 2:1 prolate has N_z ~ 0.1736 along the long axis), so it checks
that the production C++ analytic charge-Gram path gets a non-1/3 demag right -- broadening the operator's
validated domain to smooth bodies with anisotropic demag.

The C++ analytic charge Gram lands within ~0.3% of the analytic N_z, confirming the analytic surface
integral is correct for a general (non-equilateral, curved-surface) triangulation.
"""
import os
import sys
from math import log, sqrt

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "vim_legacy"))
from radia.vim import hdiv_demag_solve  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Ellipsoid, Pnt, Vec  # noqa: E402


def _prolate_Nz_analytic(c, a=1.0):
    """Demag factor along the long (z) axis of a prolate spheroid (semi-axes a=b<c)."""
    e = sqrt(1.0 - (a / c) ** 2)
    return (1.0 - e * e) / e ** 2 * ((1.0 / (2.0 * e)) * log((1.0 + e) / (1.0 - e)) - 1.0)


def test_prolate_spheroid_demag_factor():
    c = 2.0                                              # 2:1 prolate, long axis = z
    geo = CSGeometry()
    geo.Add(Ellipsoid(Pnt(0, 0, 0), Vec(1, 0, 0), Vec(0, 1, 0), Vec(0, 0, c)))
    with ng.TaskManager():
        mesh = ng.Mesh(geo.GenerateMesh(maxh=0.6))
        Dw = hdiv_demag_solve(mesh, mu_r=1000.0, H_ext=ng.CoefficientFunction((0, 0, 1.0)))["demag"]
    Na = _prolate_Nz_analytic(c)
    assert 0.15 < Na < 0.20, f"analytic prolate N_z sanity: {Na:.4f}"
    assert abs(Dw - Na) < 1.5e-2 * Na, f"prolate N_z {Dw:.4f} not within 1.5% of analytic {Na:.4f}"
