"""Golden: the ngsolve.bem-style HDiv-VIM DemagOperator (radia.vim.DemagOperator).

Mirrors ngsolve.bem: construct from an HDiv FESpace (order from the fes), expose `.mat` -- an H-matrix-
backed NGSolve BaseMatrix N = B^T G B that composes with NGSolve solvers.  order=0 (RT0) and order=p go
through ONE call.  Locks:
  * DemagFactor (the Rayleigh quotient) ~ 1/3 on a cube; order-0 uses the exact analytic Gram (fast),
    order>0 the Sauter-Schwab quadrature, agreeing to within the high-order quad error (the unified API
    gives the same physics at every order; RT0 is just order=0);
  * `.mat` is a real NGSolve BaseMatrix: a matvec works AND it composes into (M + N) which a NGSolve GMRes
    solve runs to convergence -- i.e. the operator plugs into NGSolve's framework like ngsolve.bem's
    SingleLayerPotentialOperator.mat.

NGSolve + Netgen required (importorskip); the C++ charge-Gram H-matrix backs `.mat`.  Meshing is
non-deterministic, so the demag is checked against a band + the order-invariance structure.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt  # noqa: E402
from ngsolve.krylovspace import GMRes  # noqa: E402

from radia.vim import DemagOperator  # noqa: E402  (the ngsolve.bem-style public API)


def _cube(h):
    return ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=h))


def test_demagoperator_factor_rt1():
    """DemagOperator(HDiv(mesh, order=1)).DemagFactor -> ~1/3 on a cube.  HDiv-VIM is RT1-only: RT0 (order=0)
    is retired (per-element inaccurate) and RT2+ is retired, so only RT1 is exercised (and RT0/RT2 must raise)."""
    mesh = _cube(0.8)
    with ng.TaskManager():
        D = DemagOperator(ng.HDiv(mesh, order=1), eps=1e-7).DemagFactor(ng.CF((0, 0, 1)))
        # RT0 + RT2 are retired -> the operator rejects them (No-Fallbacks)
        for bad in (0, 2):
            with pytest.raises(ValueError, match="RT1"):
                DemagOperator(ng.HDiv(mesh, order=bad), eps=1e-7)
    assert 0.31 < D < 0.345, f"RT1 DemagFactor {D:.5f} not ~1/3"


def test_demagoperator_mat_composes_with_ngsolve():
    """`.mat` is an NGSolve BaseMatrix: matvec works AND (M + N) feeds a NGSolve GMRes solve (ngsolve.bem
    idiom -- the operator plugs into NGSolve's solver framework)."""
    mesh = _cube(0.8)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        N = DemagOperator(fes, eps=1e-7)
        assert N.mat.height == fes.ndof and N.mat.width == fes.ndof
        u, v = fes.TnT()
        M = ng.BilinearForm(u * v * ng.dx).Assemble()
        gfu = ng.GridFunction(fes)
        gfu.Set(ng.CF((0, 0, 1)))
        y = N.mat.CreateColVector()
        N.mat.Mult(gfu.vec, y)
        assert np.linalg.norm(y.FV().NumPy()) > 0, "N.mat matvec produced a zero vector"
        # SPD system (M SPD + N PSD) -> GMRes through the composed BaseMatrix
        A = M.mat + N.mat
        rhs = (M.mat * gfu.vec).Evaluate()
        sol = GMRes(A=A, b=rhs, freedofs=fes.FreeDofs(), tol=1e-8, maxsteps=300, printrates=False)
        res = (A * sol - rhs).Evaluate()
        assert np.linalg.norm(res.FV().NumPy()) < 1e-5 * np.linalg.norm(rhs.FV().NumPy()), "GMRes did not converge"


# (test_demagoperator_gauss_backend_matches_analytic removed 2026-06-29: the 'gauss' point-operator backend
#  is retired -- HDiv-VIM (RT1, tet) uses the analytic charge Gram only.  See test_hdiv_vim_rt1_contract.py.)
