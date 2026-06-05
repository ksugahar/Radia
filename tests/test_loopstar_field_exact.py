"""Golden test: the loop-star gauge (rad.SetLoopStarGauge) is FIELD-EXACT.

SolveLoopStar solves the MSC system in the tree-cotree (loop/star) split and
KEEPS the loop content via block Gauss-Seidel iterative refinement (star solve
<-> loop solve on the true residual), so the external field matches the direct /
plain-BiCGSTAB solution at every mu_r. This locks two things:

  1. (CI-safe) cube in a uniform background: loop-star field == plain field, and
     the keep-loops block Gauss-Seidel drives the full residual down.  Catches any
     regression that corrupts the loop-star solve.
  2. (LAB-only, skipped if the C-type mesh data is absent) the C-type electromagnet
     -- the one geometry whose topological loops are NOT externally field-silent,
     so REMOVING them (the old behaviour) gave a ~0.5% field error.  This subtest
     fails if the keep-loops recovery is reverted.

Background: examples/c_type_electromagnet/nonlinear/{loopstar_lowmu_sweep,
ctype_loopstar_test}.py and memory/project_loopstar_lowmu_finding_2026_05_31.md.
"""
import os
import sys
import importlib.util

import numpy as np
import pytest

import radia as rad


def _build_cube(nx, mu):
    """nx^3 unit-cell soft-iron cube (MatLin mu) in a uniform 0.5 T background."""
    a = 0.01
    hexes = []
    mat = rad.MatLin(mu)
    for i in range(nx):
        for j in range(nx):
            for k in range(nx):
                x0, y0, z0 = i * a, j * a, k * a
                v = [[x0, y0, z0], [x0 + a, y0, z0], [x0 + a, y0 + a, z0], [x0, y0 + a, z0],
                     [x0, y0, z0 + a], [x0 + a, y0, z0 + a], [x0 + a, y0 + a, z0 + a], [x0, y0 + a, z0 + a]]
                o = rad.ObjHexahedron(v, [0, 0, 0])
                rad.MatApl(o, mat)
                hexes.append(o)
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, 0.5])
    return rad.ObjCnt(hexes + [bkg])


def _solve_field(builder, mu, loopstar, pts):
    rad.SetDeflateNullspace(False, 0.0)
    rad.SetLoopStarGauge(bool(loopstar))
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    model = builder(mu)
    rad.Solve(model, 0.0001, 200, 2)            # method=2 (HACApK)
    B = [np.array(rad.Fld(model, 'b', p)) for p in pts]
    kl = rad.GetKeepLoopStats() if loopstar else None
    rad.UtiDelAll()
    rad.SetLoopStarGauge(False)
    return B, kl


@pytest.mark.parametrize("mu", [2.0, 100.0, 10000.0])
def test_loopstar_cube_field_exact(mu):
    """Cube: loop-star field == plain field; keep-loops GS converges."""
    pts = [[0.05, 0.0, 0.0], [0.0, 0.05, 0.0], [0.0, 0.0, 0.08], [0.06, 0.06, 0.06]]
    builder = lambda m: _build_cube(4, m)
    Bp, _ = _solve_field(builder, mu, False, pts)
    Bl, kl = _solve_field(builder, mu, True, pts)
    dmax = max(np.linalg.norm(bl - bp) / (np.linalg.norm(bp) + 1e-30)
               for bp, bl in zip(Bp, Bl))
    assert kl["n_loop"] > 0, "loop-star built no loop basis (back-sub not exercised)"
    assert kl["res_final_rel"] < 1e-4, f"keep-loops GS did not converge: {kl}"
    assert dmax < 1e-5, f"loop-star field deviates from plain by {dmax:.2e} at mu={mu}"


def _load_ctype():
    """Import the C-type loader from examples; return (load_geometry, build_model) or None."""
    here = os.path.dirname(os.path.abspath(__file__))
    nl = os.path.join(here, "..", "examples", "c_type_electromagnet", "nonlinear")
    db = os.path.join(nl, "deflation_benchmark.py")
    if not os.path.isfile(db):
        return None
    sys.path.insert(0, nl)
    spec = importlib.util.spec_from_file_location("deflation_benchmark", db)
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
        mod.load_geometry("6x6x6")              # probes that the mesh data exists
    except Exception:
        return None
    return mod


@pytest.mark.slow
@pytest.mark.parametrize("mu", [2.0, 100.0])
def test_loopstar_ctype_field_exact(mu):
    """C-type (the discriminating geometry): keep-loops fixes the old ~0.5% error.

    Skipped where the C-type mesh data is unavailable (e.g. CI). On LAB this is the
    regression guard that fails if loop-star reverts to removing the loops.
    """
    mod = _load_ctype()
    if mod is None:
        pytest.skip("C-type mesh data not available")
    nodes, elements = mod.load_geometry("6x6x6")

    def builder(mu_):
        model, _yoke, _ne = mod.build_model(nodes, elements, ("linear", mu_))
        return model

    pts = [[0, 0, 0], [0.0, 0.0, 0.02], [0.03, 0.0, 0.0]]
    Bp, _ = _solve_field(builder, mu, False, pts)
    Bl, kl = _solve_field(builder, mu, True, pts)
    dmax = max(np.linalg.norm(bl - bp) / (np.linalg.norm(bp) + 1e-30)
               for bp, bl in zip(Bp, Bl))
    assert dmax < 1e-3, (
        f"loop-star C-type field deviates from plain by {dmax:.2e} at mu={mu} "
        f"(keep-loops recovery may be broken): {kl}")


# ----------------------------------------------------------------------------
# H-ILU golden: the A_SS H-LU-preconditioned BiCGStab path (SetLoopStarHMatMode 2)
# replaces the dense K-dense LU at large scale (>15^3, where K-dense's ns^2 caps
# at 8 GB).  Two locks:
#   1. (CI-safe) cube: H-ILU reproduces the K-dense loop-star solve -- a broken
#      H-LU factor (e.g. the materialize-free triangular solves htrsm_lln_dense_rhs
#      / htrsm_lUTn_dense_rhs) makes hlu_solve_relerr blow up.
#   2. (LAB-only) C-type: H-ILU is field-exact AND the factor materialize scratch
#      stays LOW.  Pre-fix (2026-06-04) the htrsm mixed fallback densified the
#      internal L/U factor (~2000 Me at 6^3); the materialize-free recursive solves
#      cut it to ~50 Me.  This subtest fails if that de-materialize is reverted.
# Background: memory/project_ass_scalable_solver_2026_06_04.md.
# ----------------------------------------------------------------------------
def _hilu_stats(builder, mu, tol):
    """Run the A_SS H-ILU path (mode 2); return (GetStarHMatStats, HLUMaterializeStats,
    GetKeepLoopStats).  The keep-loops res_final_rel is the FIELD-correctness gate: it is
    the relative residual of the FULL system after the loop recovery, so res_final_rel<1e-4
    means the loop-star+H-ILU field matches the plain solve WITHOUT needing K-dense/ELF
    (the only viable correctness anchor at >15^3).  res>1 means the keep-loops GS diverged
    (as it does on the distorted C-type at the unphysical mu_r=1e5)."""
    rad.SetDeflateNullspace(False, 0.0)
    rad.SetLoopStarGauge(True)
    rad.SetLoopStarHMatMode(2)
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    rad.SetStarHLUTruncTol(tol)
    try:
        rad.HLUSetTruncTol(tol)
    except Exception:
        pass
    model = builder(mu)
    rad.Solve(model, 0.0001, 200, 2)            # method=2 (HACApK); mode 2 -> H-LU precond
    st = rad.GetStarHMatStats()
    try:
        ms = rad.HLUMaterializeStats()
    except Exception:
        ms = {}
    try:
        kl = rad.GetKeepLoopStats()
    except Exception:
        kl = {}
    rad.SetLoopStarHMatMode(0)
    rad.SetLoopStarGauge(False)
    rad.UtiDelAll()
    return st, ms, kl


@pytest.mark.slow
def test_hilu_cube_matches_kdense():
    """H-ILU (mode 2) reproduces the K-dense loop-star solve on a cube (mu=1e5).

    The cube is REGULAR hexes -> its loops ARE ker(N) -> keep-loops converges even at
    mu=1e5 (the divergence is a DISTORTED-hex phenomenon).  This locks the H-LU factor."""
    st, _, _ = _hilu_stats(lambda m: _build_cube(6, m), 1e5, 1e-5)
    assert st["n_star"] > 10, f"A_SS too small to exercise the H-LU tree: {st}"
    assert st["factor_time_s"] >= 0.0, f"H-LU factor failed (factor_time_s<0): {st}"
    assert st["hlu_solve_relerr"] < 1e-5, (
        f"H-ILU solve deviates from K-dense by {st['hlu_solve_relerr']:.2e} "
        f"-- the H-LU factor (materialize-free htrsm) may be broken: {st}")


@pytest.mark.slow
def test_hilu_ctype_dematerialized():
    """C-type (distorted hex) at mu_r=1e4: H-ILU star solve matches K-dense, the
    keep-loops field recovery CONVERGES (field-correct), AND the factor materialize
    stays low.  Three locks in one: (1) relerr<1e-5 = correct H-LU factor; (2)
    keep-loops res_final_rel<1e-4 = field-exact (the >15^3 anchor); (3) materialize
    band <200 Me = the triangular-solve de-materialize (htrsm_lln_dense_rhs /
    htrsm_lUTn_dense_rhs) is not reverted (pre-fix ~2000 Me at 6^3; post-fix ~50 Me).
    Skipped where C-type data is absent.
    """
    mod = _load_ctype()
    if mod is None:
        pytest.skip("C-type mesh data not available")
    nodes, elements = mod.load_geometry("6x6x6")
    builder = lambda m: mod.build_model(nodes, elements, ("linear", m))[0]
    # mu_r=1e4 is the field-correct regime (real soft iron): the keep-loops GS converges,
    # so the loop-star+H-ILU field matches the plain solve.  (At the unphysical mu_r=1e5 the
    # keep-loops DIVERGES on the distorted C-type -- a loop-star limit, not an H-ILU bug;
    # see memory/project_ass_scalable_solver_2026_06_04.md.)
    st, ms, kl = _hilu_stats(builder, 1e4, 1e-4)
    assert st["hlu_solve_relerr"] < 1e-5, (
        f"H-ILU star solve not matching K-dense: relerr={st['hlu_solve_relerr']:.2e}: {st}")
    klr = kl.get("res_final_rel", 1e9)
    assert klr < 1e-4, (
        f"keep-loops did NOT converge (res_final_rel={klr:.2e}) -- the loop-star field is "
        f"wrong; the full loop-star+H-ILU pipeline is only field-exact for mu_r<=1e4: {kl}")
    n_elems = ms.get("n_elems", -1)
    assert 0 <= n_elems < 200e6, (
        f"H-LU materialize scratch {n_elems/1e6:.0f} Me exceeds the de-materialize "
        f"band (pre-fix ~2000 Me; post-fix ~50 Me) -- the htrsm materialize-free "
        f"triangular solve may be reverted: n_elems={n_elems}")
