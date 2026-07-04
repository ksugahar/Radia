"""Public-API gate: the production entry `radia.vim.Solve` now solves a PURE-HEX order-1 mesh
(LINEAR + NONLINEAR), not only tet.

The hex RT1 charge Gram was wired at `ChargeGram(HDiv(hexmesh, order=1))` (golden
test_hdiv_vim_hex_rt1_wiring), and `_solve_highorder`'s linear (symmetric mass-Riesz CG) and nonlinear
(all-C++ energy-Newton) paths are Gram-AGNOSTIC -- so hex flows through the same production code as tet.
Reviewed + the topology guard relaxed 2026-07-04: pure-tet, pure-hex, and pure-wedge order-1 meshes are
accepted.  The 'auto' dispatch now routes mesh-backed pure HEX/WEDGE soft irons into HDiv-VIM as well;
mixed / pyramid meshes still fail loud or route to the collocation MMMM bridge.

This locks:
  * hex LINEAR: exact cube demag 1/3, transverse ~0, and agreement with collocation MMMM on the SAME mesh;
  * hex NONLINEAR (BH table): converges (bounded iters) and agrees with collocation MMMM;
  * wedge: direct HDiv-VIM solve works; mixed / pyramid remain outside this path.

The known bursty NGSolve GetTrafo first-touch flake in the hex-charge extraction (memory
ngsolve-gettrafo-first-touch-garbage) is retried (the determinism contract fail-LOUD raises; rerun is the
documented remedy) so the gate is CI-robust.
"""
import math

import numpy as np
import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia.vim import Solve, MeshSoftIron  # noqa: E402

pytestmark = pytest.mark.slow

MU0 = 4.0e-7 * math.pi
L = 0.02
H0 = 1000.0
MU_R = 100.0
_mp = lambda x, y, z: (L * (x - 0.5), L * (y - 0.5), L * (z - 0.5))

# saturating BH table (chi0 = 2000) both backends consume
_CHI0, _MSAT = 2000.0, 1.6 / MU0
_Hs = np.concatenate([[0.0], np.logspace(-1, 7, 80)])
_Ms = _CHI0 * _Hs / (1.0 + _CHI0 * _Hs / _MSAT)
BH = [[float(h), float(b)] for h, b in zip(_Hs, MU0 * (_Hs + _Ms))]


def _cube(n):
    return MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n, mapping=_mp)


def _hdiv_hex_retry(mesh, **kw):
    """hdiv_demag_solve with retries on the KNOWN bursty GetTrafo first-touch flake (fail-loud, rerun is the
    documented remedy).  Each attempt re-draws the hex-charge extraction."""
    last = None
    for _ in range(5):
        try:
            with ng.TaskManager():
                return Solve(mesh, **kw)
        except RuntimeError as e:
            if "GetTrafo lattice evaluation unstable" in str(e):
                last = e
                continue
            raise
    raise last


def _collocation_mz(n, nonlinear):
    """Volume-average M_z of the SAME hex cube via the collocation MMMM backend (cross-method reference)."""
    rad.UtiDelAll()
    with ng.TaskManager():
        core = MeshSoftIron(_cube(n), mu_r=MU_R)
    if nonlinear:
        rad.MatApl(core, rad.MatSatIsoTab(BH))
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])
    rad.Solve(rad.ObjCnt([core, bkg]), 1e-6, 3000, 0, demag_backend="collocation_mmmm")
    mz = float(np.mean([m[2] for (_c, m) in rad.ObjM(core)]))
    rad.UtiDelAll()
    return mz


def test_hdiv_demag_solve_hex_linear():
    """hdiv_demag_solve on a pure-hex cube: exact demag 1/3, transverse ~0, agrees with collocation MMMM."""
    res = _hdiv_hex_retry(_cube(6), mu_r=MU_R, H_ext=ng.CoefficientFunction((0, 0, H0)))
    assert res["nonlinear"] is False
    assert abs(res["demag"] - 1.0 / 3.0) < 5e-3, f"hex cube demag {res['demag']} off 1/3"
    mz = res["M_avg"][2]
    assert abs(res["M_avg"][0]) < 1e-2 * mz and abs(res["M_avg"][1]) < 1e-2 * mz, "transverse M not ~0"
    assert res["iters"] < 100, f"hex linear CG not mesh-robust: {res['iters']}"
    mz_c = _collocation_mz(6, False)
    rel = abs(mz - mz_c) / abs(mz_c)
    assert rel < 0.03, f"hex HDiv Mz {mz:.1f} vs collocation MMMM {mz_c:.1f} rel {rel:.2e}"


def test_hdiv_demag_solve_hex_nonlinear():
    """hdiv_demag_solve on a pure-hex cube with a BH table: energy-Newton converges, agrees with MMMM."""
    res = _hdiv_hex_retry(_cube(6), bh_table=BH, H_ext=ng.CoefficientFunction((0, 0, H0)))
    assert res["nonlinear"] is True
    assert res["iters"] < 100, f"hex energy-Newton not bounded: {res['iters']}"
    mz = res["M_avg"][2]
    mz_c = _collocation_mz(6, True)
    rel = abs(mz - mz_c) / abs(mz_c)
    assert rel < 0.03, f"hex HDiv nonlinear Mz {mz:.1f} vs collocation MMMM {mz_c:.1f} rel {rel:.2e}"


def test_hdiv_demag_solve_wedge_linear():
    """A prism/wedge (6-vertex) mesh now SOLVES via the C++ wedge-mode charge Gram (2026-07-04, memory
    hdiv-tet-hex-coupling-pyramid-gated).  Cross-method + cross-mesh check: the prism-cube HDiv-VIM
    volume-average M_z matches the SAME-size HEX-cube collocation MMMM reference (both discretize the same
    cube -> ~1/3 demag -> the same converged M_z, up to discretization/method).  The C++ wedge Gram is
    eig(M_mass^-1 N) in [0,1] (0.992/0.998 @ n=2/3 in the de-risk)."""
    try:
        mesh = MakeStructured3DMesh(prism=True, nx=2, ny=2, nz=2, mapping=_mp)
    except TypeError:
        pytest.skip("this NGSolve MakeStructured3DMesh has no prism= kwarg")
    verts = {len(el.vertices) for el in mesh.Elements(ng.VOL)}
    if verts != {6}:
        pytest.skip(f"prism mesh did not produce pure wedges (got {sorted(verts)})")
    res = _hdiv_hex_retry(mesh, mu_r=MU_R, H_ext=ng.CoefficientFunction((0, 0, H0)))   # H field, A/m
    mz = res["M_avg"][2]
    mz_c = _collocation_mz(2, False)                              # hex-cube collocation MMMM reference
    rel = abs(mz - mz_c) / abs(mz_c)
    assert rel < 0.06, f"wedge HDiv linear Mz {mz:.1f} vs hex collocation MMMM {mz_c:.1f} rel {rel:.2e}"
    assert res["nonlinear"] is False


def test_rad_solve_demag_backend_hdiv_on_hex():
    """The wrapper path: rad.Solve(demag_backend='hdiv') on a HEX soft_iron routes through _solve_via_hdiv
    -> hdiv_demag_solve and SOLVES (no TET-only raise), writing per-element M back to the iron handles."""
    rad.UtiDelAll()
    with ng.TaskManager():
        iron = MeshSoftIron(_cube(4), mu_r=MU_R)
    src = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])
    top = rad.ObjCnt([iron, src])
    last = None
    for _ in range(5):
        try:
            with ng.TaskManager():
                rad.Solve(top, 1e-6, 3000, 0, demag_backend="hdiv")    # explicit hdiv on hex -> now solves
            break
        except RuntimeError as e:
            if "GetTrafo lattice evaluation unstable" in str(e):
                last = e
                continue
            raise
    else:
        raise last
    M = np.array([m for (_c, m) in rad.ObjM(iron)], float)             # hdiv-solved M written back via ObjSetM
    assert M.shape[0] == 64, f"expected 64 hex, got {M.shape[0]}"
    mz = float(M[:, 2].mean())
    assert mz > 1e2, f"hdiv-solved hex Mz not substantial: {mz}"       # mu_r=100, H0=1000 -> Mz ~3400
    assert abs(float(M[:, 0].mean())) < 1e-2 * mz and abs(float(M[:, 1].mean())) < 1e-2 * mz
    rad.UtiDelAll()


def test_rad_solve_auto_on_hex_uses_hdiv():
    """The wrapper path with the production default: rad.Solve(auto) on HEX soft_iron_from_mesh returns the
    HDiv-VIM result dict and writes per-element M back."""
    rad.UtiDelAll()
    rad.set_demag_backend("auto")
    with ng.TaskManager():
        iron = MeshSoftIron(_cube(3), mu_r=MU_R)
    src = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])
    top = rad.ObjCnt([iron, src])
    last = None
    for _ in range(5):
        try:
            with ng.TaskManager():
                res = rad.Solve(top, 1e-6, 3000, 0)
            break
        except RuntimeError as e:
            if "GetTrafo lattice evaluation unstable" in str(e):
                last = e
                continue
            raise
    else:
        raise last

    assert isinstance(res, dict), type(res)
    assert res["n_el"] == 27
    assert res["order"] == 1
    M = np.array([m for (_c, m) in rad.ObjM(iron)], float)
    assert M.shape[0] == 27
    assert float(M[:, 2].mean()) > 1e2
    rad.UtiDelAll()
    rad.set_demag_backend("auto")
