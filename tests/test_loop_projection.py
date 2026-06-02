"""Golden test: Helmholtz-Hodge loop removal (rad.SetLoopProjection).

After the (plain) HACApK solve, ProjectOutLoops subtracts the non-physical loop
(circulating surface-charge = ker(N)) component from sigma by a CG projection off
the topological cycle space: c solves (L^T L) c = L^T sigma, sigma -= L c. This
locks three properties:

  1. The loop component is removed: GetLoopProjStats loop_after << loop_before.
  2. The projection is cheap and well-conditioned: cg_iters is small and roughly
     mu_r-independent (L^T L cond ~ O(1)) -> no significant slowdown.
  3. N annihilates loops, so N*sigma (the on-element field) is unchanged -> the
     external field is essentially the same as plain (dB/B small); removing the
     loop changes the magnetization DISTRIBUTION (non-physical circulation), not
     the physical field.

This is the SCALABLE, mu_r-independent way to get a loop-free physical answer
(uses only the O(N) topological cycle basis -- NOT the +3 non-topological
near-null modes, which grow with N and are NOT loops). Background:
memory/project_loopstar_lowmu_finding_2026_05_31.md.
"""
import os
import sys
import importlib.util

import numpy as np
import pytest

import radia as rad


def _build_cube(nx, mu):
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


def _solve(builder, mu, projection, pts):
    rad.SetDeflateNullspace(False, 0.0)
    rad.SetLoopStarGauge(False)
    rad.SetLoopProjection(bool(projection))
    rad.SolverConfig(bicgstab_tol=1e-8, hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)
    model = builder(mu)
    rad.Solve(model, 0.0001, 4000, 2)
    B = [np.array(rad.Fld(model, 'b', p)) for p in pts]
    st = rad.GetLoopProjStats() if projection else None
    rad.UtiDelAll()
    rad.SetLoopProjection(False)
    return B, st


@pytest.mark.parametrize("mu", [2.0, 100.0, 10000.0])
def test_loop_projection_cube(mu):
    """Cube: loop removed, CG cheap, external field essentially unchanged."""
    pts = [[0.05, 0.0, 0.0], [0.0, 0.0, 0.08], [0.06, 0.06, 0.06]]
    builder = lambda m: _build_cube(4, m)
    Bp, _ = _solve(builder, mu, False, pts)
    Bj, st = _solve(builder, mu, True, pts)
    assert st["n_loop"] > 0, "no loop basis built (projection not exercised)"
    # loop component removed by many orders of magnitude
    assert st["loop_after"] < 1e-6 * st["loop_before"], (
        f"loop not removed: before={st['loop_before']:.2e} after={st['loop_after']:.2e}")
    # well-conditioned -> few CG iters (cube cycle basis is orthonormal)
    assert st["cg_iters"] <= 30, f"projection CG too many iters: {st['cg_iters']}"
    # N*sigma unchanged -> external field essentially the same as plain
    dmax = max(np.linalg.norm(bj - bp) / (np.linalg.norm(bp) + 1e-30)
               for bp, bj in zip(Bp, Bj))
    assert dmax < 1e-6, f"loop removal changed the external field by {dmax:.2e} at mu={mu}"


def _load_ctype():
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
        mod.load_geometry("6x6x6")
    except Exception:
        return None
    return mod


@pytest.mark.slow
@pytest.mark.parametrize("mu", [2.0, 100.0])
def test_loop_projection_ctype(mu):
    """C-type: topological cycle basis removes the loop on the 1/4 model too.

    Skipped where the C-type mesh data is unavailable (e.g. CI). cg_iters stays
    mu_r-independent (~54) -> the projection is well-conditioned even on the
    multiply-connected 1/4 geometry.
    """
    mod = _load_ctype()
    if mod is None:
        pytest.skip("C-type mesh data not available")
    nodes, elements = mod.load_geometry("6x6x6")
    SCALE = mod.SCALE

    def builder(mu_):
        mat = rad.MatLin(float(mu_))
        hexes = []
        for nid in elements:
            v = [[nodes[n][0] * SCALE, nodes[n][1] * SCALE, nodes[n][2] * SCALE] for n in nid]
            o = rad.ObjHexahedron(v, [0, 0, 0])
            rad.MatApl(o, mat)
            hexes.append(o)
        coil = rad.ObjRecCur([0, 0, 0], [0.30, 0.30, 0.20], [0, 0, 1.0e6])
        return rad.ObjCnt(hexes + [coil])

    pts = [[0, 0, 0.05], [0.03, 0, 0]]
    Bp, _ = _solve(builder, mu, False, pts)
    Bj, st = _solve(builder, mu, True, pts)
    assert st["n_loop"] > 0
    assert st["loop_after"] < 1e-6 * st["loop_before"], (
        f"loop not removed: before={st['loop_before']:.2e} after={st['loop_after']:.2e}")
    assert st["cg_iters"] <= 120, f"projection CG too many iters: {st['cg_iters']}"
    dmax = max(np.linalg.norm(bj - bp) / (np.linalg.norm(bp) + 1e-30)
               for bp, bj in zip(Bp, Bj))
    assert dmax < 1e-3, f"loop removal changed C-type external field by {dmax:.2e} at mu={mu}"
