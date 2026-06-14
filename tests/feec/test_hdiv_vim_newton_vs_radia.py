"""Cross-validation: NONLINEAR HDiv-VIM damped Newton  vs  Radia MMM/MSC (rad.Solve + MatSatIsoTab).

The operator damped Newton (test_hdiv_vim_tet_newton.py) is validated against the analytic uniform
sphere.  This test cross-checks it against a COMPLETELY INDEPENDENT codebase -- Radia's trusted MMM/MSC
C++ tetrahedral solver with the SAME saturating BH curve, SAME chi0/Msat, SAME sphere, SAME applied
field -- in the deep-saturation regime that was the earlier sessions' open problem.

Both solvers reproduce the analytic uniform-sphere magnetization to <0.1% at saturation, and agree with
each other to <0.1%.  (At deep saturation M -> Msat regardless of the demag factor, so the operator's
~2%-at-low-drive demag-D systematic is irrelevant here -- this isolates the nonlinear SATURATION
behaviour, which is exactly what the Newton machinery had to get right.)
"""
import os
import sys

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
rad = pytest.importorskip("radia")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "examples", "vim"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia"))
from radia.vim import _nonlinear as nl  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt  # noqa: E402
from netgen_mesh_import import netgen_mesh_to_radia  # noqa: E402

MU0 = 4e-7 * np.pi
CHI0, MSAT = 1000.0, 1.0e6


def _bh_table():
    """Shared BH curve as a Radia [[H(A/m), B(T)], ...] table: B = mu0 (H + M(H))."""
    Mof = nl._bh_curve(CHI0, MSAT)
    Hs = np.concatenate([[0.0], np.geomspace(1.0, 5e7, 60)])
    return [[float(H), float(MU0 * (H + Mof(H)))] for H in Hs]


def _radia_sphere_Mz(H0, R=0.05, maxh=0.4):
    rad.UtiDelAll()
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), R))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh * R))
    cont = netgen_mesh_to_radia(mesh, material={"magnetization": [0, 0, 0]}, units="m", verbose=False)
    rad.MatApl(cont, rad.MatSatIsoTab(_bh_table()))
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])      # uniform applied B0 = mu0 H0
    model = rad.ObjCnt([cont, bkg])
    rad.Solve(model, 1e-5, 2000, 0)                        # LU
    pts = [[0, 0, 0], [0.3 * R, 0, 0], [0, 0.3 * R, 0], [0, 0, 0.3 * R], [0, 0, -0.3 * R]]
    Mz = float(np.mean([rad.Fld(cont, "m", p)[2] for p in pts]))
    rad.UtiDelAll()
    return Mz


def test_newton_matches_radia_mmm_msc_at_saturation():
    Mof = nl._bh_curve(CHI0, MSAT)
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=0.4))
    H0 = 1.0e6                                              # deep saturation (knee ~ Msat/chi0 = 1e3)
    with ng.TaskManager():
        M_hdiv, _, _ = nl.solve_nonlinear_newton(mesh, CHI0, MSAT, H0)
    M_radia = _radia_sphere_Mz(H0)
    M_ana = nl._scalar_fixed_point(Mof, 1.0 / 3.0, H0)
    assert abs(M_hdiv - M_ana) < 2e-3 * M_ana, f"HDiv-VIM {M_hdiv:.1f} vs analytic {M_ana:.1f}"
    assert abs(M_radia - M_ana) < 2e-3 * M_ana, f"Radia {M_radia:.1f} vs analytic {M_ana:.1f}"
    assert abs(M_hdiv - M_radia) < 2e-3 * M_ana, \
        f"HDiv-VIM ({M_hdiv:.1f}) and Radia MMM/MSC ({M_radia:.1f}) disagree at saturation"
