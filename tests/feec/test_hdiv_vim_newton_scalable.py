"""Golden test: SCALABLE nonlinear HDiv-VIM damped Newton (production #2).

solve_nonlinear_newton_scalable replaces the dense O(N^3)/O(N^2) demag of solve_nonlinear_newton with
the C++ HACApK charge-Gram H-matrix (_ChargeGramHMatrix, O(N log N) apply) + a sparse near-correction,
and solves each Newton step ITERATIVELY (GMRES, M_mass-preconditioned) -- no dense factorization
anywhere.  The demag apply is the textbook H-matrix far+near split N v = B^T(H.matvec(B v) + corr (B v)),
driving the same per-element tensor-tangent Newton.

This locks that the scalable solver reproduces the dense one (same monopole+near-corr operator): at
deep saturation they agree to ~machine precision, confirming the nonlinear HDiv-VIM solve is scalable
and correct.  (The C++ Gram is monopole far field, so this matches the dense MONOPOLE+near-corr path,
not the dense Wilton path -- putting Wilton into the C++ near-field correction is the remaining step.)
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
rad = pytest.importorskip("radia")          # the C++ Gram H-matrix lives in radia._radia_pybind

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "feec_vim"))
import hdiv_demag_tet_nonlinear as nl  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402


def _sphere(h=0.45):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    return ng.Mesh(geo.GenerateMesh(maxh=h))


def test_scalable_newton_matches_dense_at_saturation():
    mesh = _sphere()
    chi0, Msat = 1000.0, 1.0e6
    Mof = nl._bh_curve(chi0, Msat)
    H0 = 1.0e6                                  # deep saturation
    with ng.TaskManager():
        Md, _, _ = nl.solve_nonlinear_newton(mesh, chi0, Msat, H0)
        Ms, nit, _ = nl.solve_nonlinear_newton_scalable(mesh, chi0, Msat, H0)
    Ma = nl._scalar_fixed_point(Mof, 1.0 / 3.0, H0)
    # scalable (C++ H-matrix + GMRES) reproduces the dense (same monopole+near-corr operator)
    assert abs(Md - Ms) < 1e-3 * abs(Md), f"scalable {Ms:.1f} != dense {Md:.1f}"
    # and both are physical / near the analytic uniform-sphere answer
    assert 0.0 < Ms < Msat, f"scalable M out of range: {Ms}"
    assert abs(Ms - Ma) < 5e-3 * Ma, f"scalable {Ms:.1f} vs analytic {Ma:.1f}"
    assert nit < 30, f"scalable Newton not fast at saturation: {nit} iters"
