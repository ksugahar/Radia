"""Golden: HDiv-VIM CURVED (isoparametric P2) tet NONLINEAR demag solve (curved-nonlinear increment 1).

The curved tet charge Gram (ChargeGram(curve_order=2), CurvedTri/TetPotential -- locked by
test_hdiv_vim_curved_gram.py) is now wired to the symmetric energy-Newton nonlinear solver
(_solve_nonlinear_energy_cpp): Solve(mesh, bh_table=..., curve_order=2) lifts the former
order>0-nonlinear NotImplementedError for the CURVED path.  The FLAT order>0 nonlinear path is now wired too
(verified against the analytic fixed point below).  No new C++ -- the energy-Newton is Gram-AGNOSTIC
(it consumes only the configured C++ Gram apply + mass-Riesz solve, present on every high-order Gram object).

Locks: (1) curve_order=2 + bh_table runs + reports nonlinear=True, solver='energy-newton-cpp', curve_order=2;
(2) a uniform sphere magnetizes to the analytic spheroid fixed point (~1e-2); (3) matches the analytic
nonlinear solve on the same geometry (the demag factor is curving-insensitive); (4) the curved (Duffy) Gram
stays PSD enough for the energy slide-down -- the solve CONVERGES (iters < maxit).
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import Sphere, Pnt, OCCGeometry          # noqa: E402
from radia.vim import Solve                    # noqa: E402
from radia.vim import _nonlinear as nl                    # noqa: E402

# analytic-ish soft iron (chi0=1000, Msat=1e6): synth a [[H,B]] table from the smooth M(H) curve
_CHI0, _MSAT = 1000.0, 1.0e6
_Mof = nl._bh_curve(_CHI0, _MSAT)
_H = np.concatenate([[0.0], np.logspace(-1, 7, 40)])
_B = nl._MU0 * (_H + np.array([_Mof(h) for h in _H]))
_BH = np.column_stack([_H, _B]).tolist()


def _sphere(maxh=0.6):
    return ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=maxh))


def test_curved_tet_nonlinear_runs_and_matches_analytic():
    """curve_order=2 + bh_table at RT1: the energy-Newton runs on the curved Gram and a uniform sphere matches
    the analytic spheroid fixed point AND the flat RT1 nonlinear solve (curving-insensitive demag)."""
    H0 = 5000.0
    Man = nl._scalar_fixed_point(_Mof, 1.0 / 3.0, H0)
    with ng.TaskManager():
        rc = Solve(_sphere(), bh_table=_BH, H_ext=ng.CoefficientFunction((0, 0, H0)),
                              order=1, curve_order=2)
        rf = Solve(_sphere(), bh_table=_BH, H_ext=ng.CoefficientFunction((0, 0, H0)), order=1)
    assert rc["nonlinear"] is True
    assert rc["linear_solver"] == "energy-newton-cpp", rc["linear_solver"]
    assert rc["curve_order"] == 2
    assert rc["iters"] < 300, rc["iters"]                         # converged -> curved Gram is PSD enough
    # a uniform sphere magnetizes to the analytic spheroid fixed point
    assert abs(rc["M_avg"][2] - Man) / abs(Man) < 1e-2, (rc["M_avg"][2], Man)
    # the demag factor is curving-insensitive -> curved RT1 nonlinear ~ flat RT1 nonlinear
    assert abs(rc["M_avg"][2] - rf["M_avg"][2]) / abs(rf["M_avg"][2]) < 1e-2
    # demag factor near the sphere's 1/3 (uniform-M limit)
    assert 0.30 < rc["demag"] < 0.37, rc["demag"]


def test_flat_rt1_nonlinear_matches_analytic():
    """FLAT (non-curved) RT1 nonlinear (energy-Newton on the flat high-order Gram): a uniform sphere magnetizes
    to the analytic spheroid fixed point.  Order 1 is the production nonlinear path and the
    reference is the closed-form fixed point (the former 'flat order>0 blocked' guard was removed; flat RT1
    nonlinear was verified against the analytic fixed point before promotion)."""
    H0 = 5000.0
    Man = nl._scalar_fixed_point(_Mof, 1.0 / 3.0, H0)
    with ng.TaskManager():
        r1 = Solve(_sphere(), bh_table=_BH, H_ext=ng.CoefficientFunction((0, 0, H0)), order=1)
    assert r1["nonlinear"] is True
    assert r1["linear_solver"] == "energy-newton-cpp", r1["linear_solver"]
    assert r1["iters"] < 300, r1["iters"]
    assert abs(r1["M_avg"][2] - Man) / abs(Man) < 1e-2, (r1["M_avg"][2], Man)
    assert 0.30 < r1["demag"] < 0.37, r1["demag"]


def test_curved_tet_nonlinear_jacobi_inner_matches_mass_riesz():
    """The large-run Jacobi inner preconditioner is valid on the curved Gram too."""
    H0 = 5000.0
    mesh = _sphere(maxh=0.8)
    with ng.TaskManager():
        ref = Solve(mesh, bh_table=_BH, H_ext=ng.CoefficientFunction((0, 0, H0)),
                    order=1, curve_order=2)
        tuned = Solve(mesh, bh_table=_BH, H_ext=ng.CoefficientFunction((0, 0, H0)),
                      order=1, curve_order=2, preconditioner="jacobi",
                      newton_continuation=2, newton_reuse_tangent_steps=3,
                      newton_inner_tol="auto")
    rel = np.linalg.norm(tuned["M_avg"] - ref["M_avg"]) / max(np.linalg.norm(ref["M_avg"]), 1.0)
    assert rel < 5e-4, (tuned["M_avg"], ref["M_avg"], rel)
    assert tuned["curve_order"] == 2
    assert tuned["nonlinear_solve_stats"]["nonlinear_inner_preconditioner"] == "jacobi"
