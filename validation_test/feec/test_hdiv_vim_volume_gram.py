"""Golden test: the analytic VOLUME-charge Gram (the C++ charge Gram, exact near AND far) gives a
PHYSICAL demag factor on sphere AND cube.

For NON-uniform M (nonlinear, div M != 0) the VOLUME-involving Gram blocks (cell-cell, cell-face) matter,
and a crude centroid-monopole leaves a residual.  The production C++ analytic charge Gram (PhiTet /
TriPotential, exact) closes it while keeping the linear demag factor exact (-> 1/3 on sphere AND cube).
The dense Python Gram path was removed, so the demag factor is taken through the production
hdiv_demag_solve (the C++ charge-Gram H-matrix).
"""
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

from radia.vim import Solve  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt  # noqa: E402

_HEXT = ng.CoefficientFunction((0, 0, 1.0))


def _demag(geo_kind):
    g = CSGeometry()
    if geo_kind == "sphere":
        g.Add(Sphere(Pnt(0, 0, 0), 1.0)); h = 0.4
    else:
        g.Add(OrthoBrick(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))); h = 0.3
    with ng.TaskManager():
        mesh = ng.Mesh(g.GenerateMesh(maxh=h))
        return Solve(mesh, mu_r=1000.0, H_ext=_HEXT)["demag"]


def test_volume_gram_demag_sphere():
    D = _demag("sphere")
    assert abs(D - 1.0 / 3) < 5e-3, f"C++ analytic-Gram sphere demag {D:.4f} not within 0.5% of 1/3"


def test_volume_gram_demag_cube():
    D = _demag("cube")
    assert abs(D - 1.0 / 3) < 0.06, f"C++ analytic-Gram cube demag {D:.4f} not ~1/3 (RT0 coarse band)"
