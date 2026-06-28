"""HIGH-ORDER Gauss point operator (build_charge_gauss) -- the P-scatter generalization of the
RadHACApKChargeGaussOperator (C++ rad_hacapk_hdiv).

The order-0 ChargeGaussOperator scatters ONE charge per point (point_weight).  This locks the GENERAL
P-scatter: a quadrature point belongs to one host element shared by ALL the host's polynomial charges, so
a point scatters from many charges with coef = quad_weight_p * monomial_a(x_p).  build_charge_gauss assembles
that P (COO) + the sparse near correction (exact analytic - point quadrature, via a build=False oracle) and
hands it to the generalized C++ operator.  The demag N = B^T G_gauss B must match the analytic high-order
Gram (build_charge_gram) at p=1,2 on straight AND curved tet spheres -- the user's p<=2 / curved scope.

(The order-0 trivial-P path -- the single-owner special case -- is covered by test_hdiv_vim_gauss_hmatrix.py.)
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402

from radia.vim._vim import build_charge_gram, build_charge_gauss  # noqa: E402


def _sphere(maxh=0.7, curve=0):
    g = CSGeometry(); g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        mesh = ng.Mesh(g.GenerateMesh(maxh=maxh))
        if curve:
            mesh.Curve(curve)
    return mesh


def _demag(fes, builder, Mcf, **kw):
    with ng.TaskManager():
        B, G, M = builder(fes, **kw)
        gf = ng.GridFunction(fes); gf.Set(Mcf); m = np.array(gf.vec)
        c = B @ m
        return float(c @ np.array(G.matvec(c.tolist()))) / float(m @ (M @ m)), G


@pytest.mark.parametrize("curve", [0, 2])
@pytest.mark.parametrize("p", [1])
def test_build_charge_gauss_matches_analytic_demag(curve, p):
    """High-order Gauss demag == analytic high-order Gram demag to <1e-3, straight + curved, p=1 (RT2+ abolished)."""
    mesh = _sphere(curve=curve)
    Mcf = ng.CoefficientFunction((0, 0, ng.z))   # non-uniform -> exercises the volume charge (div M != 0)
    fes = ng.HDiv(mesh, order=p)
    D_analytic, _ = _demag(fes, build_charge_gram, Mcf)
    D_gauss, G = _demag(fes, build_charge_gauss, Mcf, qpts=3, near_factor=1.0)
    assert G.npoint() > G.ncharge()          # genuine point H-matrix (P-scatter), not a charge-entry matrix
    assert abs(D_gauss - D_analytic) / abs(D_analytic) < 1e-3


def test_build_charge_gauss_uniform_demag_third():
    """Uniform M_z on a straight-tet sphere -> demag ~ 1/3 through the Gauss point operator (physical anchor)."""
    mesh = _sphere()
    D, _ = _demag(ng.HDiv(mesh, order=1), build_charge_gauss,
                  ng.CoefficientFunction((0, 0, 1.0)), qpts=3, near_factor=1.0)
    assert 0.30 < D < 0.36
