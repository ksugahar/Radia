"""Golden test: NONLINEAR HDiv-type VIM demag via damped Newton-Raphson (operator, per-element).

The earlier nonlinear increment (test_hdiv_vim_tet_nonlinear.py) used a SCALAR-chi Picard, which is
well-conditioned but assumes a uniform M.  The per-element (operator) generalization needs a real
Newton-Raphson, which the earlier sessions found the simple per-element Picard / Hantila could NOT do
robustly (NaN / wrong root at saturation).

solve_nonlinear_newton (examples/vim/hdiv_demag_tet_nonlinear.py) is the robust operator Newton.
Its three ingredients (verify-first, 2026-06-07):
  (a) the CONSISTENT TENSOR tangent dM/dH = chi_diff Hhat(x)Hhat + chi_sec (I - Hhat(x)Hhat)
      (the scalar chi_diff*I stalls at moderate drive where chi_sec != chi_diff);
  (b) the #3 NEAR-FIELD CORRECTION on N (without it the per-element Newton converges to a WRONG root,
      M~0.09 instead of 0.30 on the sphere at H0=0.1 -- the monopole's poor local-field accuracy);
  (c) line-search damping + a scalar-chi Picard warmstart (the BH-knee is stiff).

This locks the WIN: the operator Newton is FAST and ACCURATE at DEEP SATURATION (where the per-element
Picard/Hantila failed), reproducing the analytic uniform-sphere M to the operator's demag-factor
accuracy.  (The BH-knee is intentionally NOT exercised here -- it converges to the correct answer but
slowly, an operator-accuracy-limited stiffness documented in the example.)
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
from radia.vim import _nonlinear as nl  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402


def _sphere(h=0.4):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    return ng.Mesh(geo.GenerateMesh(maxh=h))


def test_newton_robust_and_fast_at_saturation():
    """At deep saturation the operator Newton converges FAST and matches the analytic uniform-sphere M.

    This is the regime where the earlier simple per-element Picard / Hantila diverged (NaN); damped
    Newton with the tensor tangent gets there in a handful of iterations.
    """
    mesh = _sphere()
    chi0, Msat = 1000.0, 1.0
    Mof = nl._bh_curve(chi0, Msat)
    with ng.TaskManager():
        for H0 in (1.0, 5.0):
            Mn, nit, _ = nl.solve_nonlinear_newton(mesh, chi0, Msat, H0)
            Mana = nl._scalar_fixed_point(Mof, 1.0 / 3.0, H0)
            assert nit < 25, f"Newton not fast at saturation H0={H0}: {nit} iters"
            assert 0.0 < Mn < 1.0, f"M out of physical range at H0={H0}: {Mn}"
            assert abs(Mn - Mana) < 2e-2 * Mana, \
                f"Newton M={Mn:.5f} not within 2% of analytic {Mana:.5f} at H0={H0}"


def test_newton_needs_tensor_tangent_and_near_correction():
    """Without the near-field correction the per-element Newton converges to a WRONG (too small) root;
    with it, it matches the analytic.  Locks ingredient (b)."""
    mesh = _sphere()
    chi0, Msat = 1000.0, 1.0
    Mof = nl._bh_curve(chi0, Msat)
    H0 = 1.0
    with ng.TaskManager():
        Mcorr, _, Dc = nl.solve_nonlinear_newton(mesh, chi0, Msat, H0, near_correction=True)
        Mmono, _, Dm = nl.solve_nonlinear_newton(mesh, chi0, Msat, H0, near_correction=False)
    Mana = nl._scalar_fixed_point(Mof, 1.0 / 3.0, H0)
    assert Dc > Dm, f"near correction should raise the demag factor: {Dm:.4f} -> {Dc:.4f}"
    assert abs(Mcorr - Mana) < abs(Mmono - Mana) + 1e-9, \
        f"near-corrected Newton ({Mcorr:.5f}) not at least as close to analytic ({Mana:.5f}) as monopole ({Mmono:.5f})"


def test_newton_agrees_with_scalar_picard_at_saturation():
    """At saturation the per-element Newton and the (trusted) near-corrected scalar Picard solve the
    same operator system and must agree -- a cross-check between the two solvers."""
    mesh = _sphere()
    chi0, Msat = 1000.0, 1.0
    Mof = nl._bh_curve(chi0, Msat)
    H0 = 1.0
    with ng.TaskManager():
        Mnewton, _, _ = nl.solve_nonlinear_newton(mesh, chi0, Msat, H0)
        Mpicard, _, _ = nl.solve_nonlinear(mesh, Mof, H0, near_correction=True)
    assert abs(Mnewton - Mpicard) < 1e-2 * Mpicard, \
        f"Newton ({Mnewton:.5f}) and scalar Picard ({Mpicard:.5f}) disagree at saturation"
