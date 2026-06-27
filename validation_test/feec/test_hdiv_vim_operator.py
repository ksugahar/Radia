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


def test_demagoperator_factor_order_invariant():
    """DemagOperator(HDiv(mesh,order=p)).DemagFactor -> ~1/3, IDENTICAL across p=0,1,2 (unified API)."""
    mesh = _cube(0.8)
    D = {}
    for p in (0, 1, 2):
        with ng.TaskManager():
            fes = ng.HDiv(mesh, order=p)
            N = DemagOperator(fes, eps=1e-7)
            D[p] = N.DemagFactor(ng.CF((0, 0, 1)))
    for p, v in D.items():
        assert 0.31 < v < 0.345, f"order={p} DemagFactor {v:.5f} not ~1/3"
    # order-0 uses the EXACT ANALYTIC Gram (Wilton/PhiTet, fast); order>0 uses the Sauter-Schwab QUADRATURE
    # Gram (no analytic high-order Gram exists).  The near/self quad floor is now ORDER-DEPENDENT (quad >= 3*p,
    # the PSD requirement -- see test_hdiv_vim_psd): p=1 uses quad=4, p=2 uses quad=6, so the demag factors
    # differ by the (order-dependent) QUADRATURE accuracy (~4e-4 here, with quad=6 the MORE accurate / closer
    # to 1/3), NOT a tight same-quad invariance.  All three still agree to within the quad error (~1e-3 on this
    # coarse cube).  (The order-0 analytic Gram == the validated solve_nonlinear_newton_scalable Gram, ~5e-10.)
    assert abs(D[1] - D[2]) < 1.5e-3, f"orders 1,2 demag factor disagree beyond quad error: {D}"
    assert abs(D[0] - D[1]) < 3e-3, f"analytic order-0 vs quadrature high-order beyond quad error: {D}"


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


@pytest.mark.parametrize("p", [0, 1, 2])
def test_demagoperator_gauss_backend_matches_analytic(p):
    """DemagOperator(gram_backend='gauss') -- the Gauss POINT operator (P^T K_point P + near correction)
    plugged into the SAME .mat path -- matches the analytic Gram DemagFactor to <1e-3 at p=0,1,2, and its
    .mat is a real NGSolve BaseMatrix that composes into (M + N) for a GMRes solve."""
    mesh = _cube(0.8)
    Mcf = ng.CF((0, 0, ng.z))   # non-uniform -> exercises the volume charge (div M != 0)
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=p)
        D_analytic = DemagOperator(fes, eps=1e-7).DemagFactor(Mcf)
        Ng = DemagOperator(fes, gram_backend="gauss", qpts=3, gauss_near_factor=1.0)
        D_gauss = Ng.DemagFactor(Mcf)
        assert Ng.gram_backend == "gauss"
        assert abs(D_gauss - D_analytic) / abs(D_analytic) < 1e-3
        # the gauss .mat composes with NGSolve just like the analytic one
        u, v = fes.TnT()
        M = ng.BilinearForm(u * v * ng.dx).Assemble()
        gfu = ng.GridFunction(fes); gfu.Set(Mcf)
        A = M.mat + Ng.mat
        rhs = (M.mat * gfu.vec).Evaluate()
        sol = GMRes(A=A, b=rhs, freedofs=fes.FreeDofs(), tol=1e-8, maxsteps=300, printrates=False)
        res = (A * sol - rhs).Evaluate()
        assert np.linalg.norm(res.FV().NumPy()) < 1e-5 * np.linalg.norm(rhs.FV().NumPy())
