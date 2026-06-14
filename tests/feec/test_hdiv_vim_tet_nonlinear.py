"""Golden test: NONLINEAR HDiv-type VIM demag (applied-field A+ solve + BH-curve Picard).

The nonlinear foundation (examples/vim/hdiv_demag_tet_nonlinear.py):
  - applied-field system A+ = (1/chi) M_mass + N (PLUS N -- the eigenvalue framing's -N is NOT the
    applied-field solve; A+ reproduces the analytic sphere M/H = chi/(1+chi D));
  - secant-chi Picard with a saturating BH curve converges + saturates (M -> M_sat);
  - the #3 near-field correction raises the demag factor toward the analytic 1/3, bringing the
    nonlinear M much closer to the analytic uniform-sphere answer.

Locks (NGSolve + Netgen required):
  - the Picard converges and the sphere SATURATES (M/Msat rises with H0, bounded < Msat);
  - the near-field-corrected nonlinear M is CLOSER to the analytic (D=1/3) than the monopole M, and
    within ~8% of it at strong drive.
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
    with ng.TaskManager():
        return ng.Mesh(geo.GenerateMesh(maxh=h))


def test_nonlinear_picard_converges_and_saturates():
    mesh = _sphere()
    Mof = nl._bh_curve(1000.0, 1.0)
    Mlow, it_lo, _ = nl.solve_nonlinear(mesh, Mof, 1e-3, near_correction=True)
    Mhigh, it_hi, _ = nl.solve_nonlinear(mesh, Mof, 1e-1, near_correction=True)
    assert it_lo < 250 and it_hi < 250, f"Picard did not converge ({it_lo}, {it_hi} iters)"
    # saturation: M rises with drive but stays bounded below M_sat=1
    assert Mhigh > Mlow, f"M did not increase with drive ({Mlow} -> {Mhigh})"
    assert 0.2 < Mhigh < 1.0, f"high-drive M not in the saturating range: {Mhigh:.4f}"


def test_nonlinear_nearcorrection_closer_to_analytic():
    """At strong drive the near-field-corrected nonlinear M is closer to the analytic (sphere D=1/3)
    uniform-sphere answer than the centroid-monopole M (the monopole under-estimates the demag)."""
    mesh = _sphere()
    Mof = nl._bh_curve(1000.0, 1.0)
    H0 = 1e-1
    Mmono, _, Dm = nl.solve_nonlinear(mesh, Mof, H0, near_correction=False)
    Mcorr, _, Dc = nl.solve_nonlinear(mesh, Mof, H0, near_correction=True)
    Mana = nl._scalar_fixed_point(Mof, 1.0 / 3.0, H0)
    assert Dc > Dm, f"near correction should raise the demag factor: {Dm:.4f} -> {Dc:.4f}"
    assert abs(Mcorr - Mana) < abs(Mmono - Mana), \
        f"near-corrected ({Mcorr:.4f}) not closer to analytic ({Mana:.4f}) than monopole ({Mmono:.4f})"
    assert abs(Mcorr - Mana) < 0.08 * Mana, f"near-corrected {Mcorr:.4f} vs analytic {Mana:.4f} too far"
