"""Golden test: SCALABLE nonlinear HDiv-VIM damped Newton (production #2).

solve_nonlinear_newton_scalable replaces the dense O(N^3)/O(N^2) demag of solve_nonlinear_newton with
the C++ HACApK charge-Gram H-matrix (_ChargeGramHMatrix, O(N log N) apply) and solves each Newton step
ITERATIVELY (GMRES, M_mass-preconditioned) -- no dense factorization anywhere.

The C++ Gram is the ANALYTIC charge Gram (PhiTet/TriPotential, exact near AND far), so the demag apply
is simply N v = B^T(H.matvec(B v)) -- no separate sparse near-correction.  It lands on the analytic
uniform-sphere fixed point AND gives a physical M on a NON-uniform-M body (cube, div M != 0) -- the
C-yoke-class case the old monopole-far + sparse-near split under-resolved.  (The dense Python Gram path
was removed, so the references here are the analytic fixed point + physical-range checks.)
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
rad = pytest.importorskip("radia")          # the C++ Gram H-matrix lives in radia._radia_pybind

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
from radia.vim import _nonlinear as nl  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt  # noqa: E402


def _sphere(h=0.6):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    return ng.Mesh(geo.GenerateMesh(maxh=h))


def _cube(h=0.6):
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5)))
    return ng.Mesh(geo.GenerateMesh(maxh=h))


def test_scalable_newton_matches_analytic_sphere():
    """Uniform-M body: the scalable analytic Newton lands on the analytic uniform-sphere answer at deep
    saturation, fast."""
    mesh = _sphere()
    chi0, Msat, H0 = 1000.0, 1.0e6, 1.0e6        # deep saturation
    Mof = nl._bh_curve(chi0, Msat)
    with ng.TaskManager():
        Ms, nit, _ = nl.solve_nonlinear_newton_scalable(mesh, chi0, Msat, H0)
    Ma = nl._scalar_fixed_point(Mof, 1.0 / 3.0, H0)
    assert 0.0 < Ms < Msat, f"scalable M out of range: {Ms}"
    assert abs(Ms - Ma) < 5e-3 * Ma, f"scalable {Ms:.1f} vs analytic {Ma:.1f}"
    assert nit < 30, f"scalable Newton not fast at saturation: {nit} iters"


def test_scalable_newton_default_gram_eps_is_tight():
    """Fix B (2026-06-09): the DEFAULT gram_eps must stay tight (<=1e-9).  The ACA charge-Gram H-matrix
    asymmetry scales DIRECTLY with eps (1e-6->15%, 1e-8->1.6%, 1e-10->3e-5); the old default 1e-7 gave
    a ~5% asymmetric/INACCURATE operator at scale (n_charge~few-thousand) -> the nonlinear Newton's
    tangent was inconsistent -> slow LINEAR convergence (iters 6->27->37 over a mesh sweep) + a 0.6%
    wrong Mz.  Tightening the default to 1e-10 restores a mesh-INDEPENDENT 6-iter quadratic convergence
    (verified C-yoke h=0.008/0.006/0.005 -> 6/6/6 iters).  Lock the default so it cannot silently
    regress to the too-loose value.  (NOT the loops/-N: B1 and B2 were both red herrings; the root
    cause was the ACA tolerance -- see docs/hdiv_vim/PRODUCTIONIZATION.md.)"""
    import inspect
    eps = inspect.signature(nl.solve_nonlinear_newton_scalable).parameters["gram_eps"].default
    assert eps <= 1e-9, (f"default gram_eps={eps:.0e} too loose -> the nonlinear Newton degrades to "
                         f"non-mesh-independent linear convergence (need <=1e-9, was the buggy 1e-7)")


def test_scalable_newton_fails_loud_not_silent_underconverged():
    """Fix A (2026-06-09): with too few Newton iters the scalable solver MUST raise (fail loud), never
    silently return an under-converged M.  The old code broke on Mavg stagnation (|M_now-M_prev|<1e-8),
    which fires during the slow globalization phase while relF is still O(0.1) -> it returned a WRONG
    (drifted) Mz (the C-yoke "drift to 509k" was this false convergence, NOT an operator error).  The
    sound test is relF = ||F||/||M_mass m|| < newton_tol, with a hard raise if maxit is hit first."""
    mesh = _sphere(0.6)
    with ng.TaskManager():
        with pytest.raises(RuntimeError, match="did NOT converge"):
            nl.solve_nonlinear_newton_scalable(mesh, 1000.0, 1.0e6, 1.0e6, maxit=1, picard_warmstart=1)


def test_scalable_newton_cube_nonuniform_physical():
    """NON-uniform-M body (cube, div M != 0 near corners -> volume charge): the scalable analytic Newton
    gives a PHYSICAL M (in range, partially saturated) and converges.  The analytic C++ Gram (PhiTet
    volume charge) handles the div M != 0 case the old monopole-far + sparse-near split under-resolved."""
    mesh = _cube()
    chi0, Msat, H0 = 1000.0, 1.0e6, 2.0e3        # partial saturation -> genuinely non-uniform M
    with ng.TaskManager():
        Ms, nit, D = nl.solve_nonlinear_newton_scalable(mesh, chi0, Msat, H0)
    # M well inside (0, Msat), clearly responding (not a numerical zero), and consistent with the
    # uniform-body linear-response upper bound M <= chi0/(1+chi0 D) H0 (demag-limited; here ~6000-7000).
    assert 0.0 < Ms < Msat, f"scalable cube M out of range: {Ms}"
    assert Ms > 100.0, f"scalable cube M suspiciously small at H0={H0}: {Ms}"
    assert Ms < 1.5 * chi0 / (1.0 + chi0 * D) * H0, f"scalable cube M too large vs linear bound: {Ms}"
    assert nit < 60, f"scalable cube Newton too slow: {nit} iters"
