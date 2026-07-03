"""Golden test: NONLINEAR HDiv-type VIM demag via damped Newton-Raphson (operator, per-element).

solve_nonlinear_newton routes through the production C++ analytic charge-Gram demag operator (the dense
Python charge-Gram Newton was removed).  This locks the WIN: the operator Newton is FAST and ACCURATE at
DEEP SATURATION (where the per-element Picard/Hantila failed), reproducing the analytic uniform-sphere M
to the operator's demag-factor accuracy.
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "vim_legacy"))
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
            assert nit < 30, f"Newton not fast at saturation H0={H0}: {nit} iters"
            assert 0.0 < Mn < 1.0, f"M out of physical range at H0={H0}: {Mn}"
            assert abs(Mn - Mana) < 2e-2 * Mana, \
                f"Newton M={Mn:.5f} not within 2% of analytic {Mana:.5f} at H0={H0}"
