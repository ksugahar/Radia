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
