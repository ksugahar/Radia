"""Golden: HDiv-VIM CURVED (isoparametric P2) tet NONLINEAR demag solve (curved-nonlinear increment 1).

The curved tet charge Gram (build_charge_gram(curve_order=2), CurvedTri/TetPotential -- locked by
test_hdiv_vim_curved_gram.py) is now wired to the symmetric energy-Newton nonlinear solver
(_solve_nonlinear_energy_cpp): hdiv_demag_solve(mesh, bh_table=..., curve_order=2) lifts the former
order>0-nonlinear NotImplementedError for the CURVED path (flat order>0 nonlinear stays blocked). No new C++
-- the energy-Newton is Gram-AGNOSTIC (it consumes only H.matvec + H.solve_linear_material_mass_riesz, both
present on the curved m_highorder Gram object).

Locks: (1) curve_order=2 + bh_table runs + reports nonlinear=True, solver='energy-newton-cpp', curve_order=2;
(2) a uniform sphere magnetizes to the analytic spheroid fixed point (~1e-2); (3) matches the RT0 (order=0)
nonlinear solve on the same geometry (the demag factor is curving-insensitive); (4) the curved (Duffy) Gram
stays PSD enough for the energy slide-down -- the solve CONVERGES (iters < maxit).
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import Sphere, Pnt, OCCGeometry          # noqa: E402
from radia.vim import hdiv_demag_solve                    # noqa: E402
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
    """curve_order=2 + bh_table: the energy-Newton runs on the curved Gram and a uniform sphere matches
    the analytic spheroid fixed point AND the RT0 nonlinear solve (curving-insensitive demag)."""
    H0 = 5000.0
    Man = nl._scalar_fixed_point(_Mof, 1.0 / 3.0, H0)
    with ng.TaskManager():
        rc = hdiv_demag_solve(_sphere(), bh_table=_BH, H_ext=ng.CoefficientFunction((0, 0, H0)),
                              order=0, curve_order=2)
        r0 = hdiv_demag_solve(_sphere(), bh_table=_BH, H_ext=ng.CoefficientFunction((0, 0, H0)), order=0)
    assert rc["nonlinear"] is True
    assert rc["linear_solver"] == "energy-newton-cpp", rc["linear_solver"]
    assert rc["curve_order"] == 2
    assert rc["iters"] < 300, rc["iters"]                         # converged -> curved Gram is PSD enough
    # a uniform sphere magnetizes to the analytic spheroid fixed point
    assert abs(rc["M_avg"][2] - Man) / abs(Man) < 1e-2, (rc["M_avg"][2], Man)
    # the demag factor is curving-insensitive -> curved nonlinear ~ RT0 nonlinear
    assert abs(rc["M_avg"][2] - r0["M_avg"][2]) / abs(r0["M_avg"][2]) < 1e-2
    # demag factor near the sphere's 1/3 (uniform-M limit)
    assert 0.30 < rc["demag"] < 0.37, rc["demag"]


def test_flat_highorder_nonlinear_still_blocked():
    """The guard is narrowed, not removed: flat (non-curved) order>0 nonlinear still fails loud."""
    with ng.TaskManager():
        with pytest.raises(NotImplementedError, match="flat order>0"):
            hdiv_demag_solve(_sphere(), bh_table=_BH, H_ext=ng.CoefficientFunction((0, 0, 5000.0)), order=2)
