"""Public-API gate: the production entry `radia.vim.Solve` now solves a PURE-HEX order-1 mesh
(LINEAR + NONLINEAR), not only tet.

The hex RT1 charge Gram was wired at `ChargeGram(HDiv(hexmesh, order=1))` (golden
test_hdiv_vim_hex_rt1_wiring), and `_solve_highorder`'s linear (symmetric mass-Riesz CG) and nonlinear
(all-C++ energy-Newton) paths are Gram-AGNOSTIC -- so hex flows through the same production code as tet.
Reviewed + the topology guard relaxed 2026-07-04: pure-tet, pure-hex, and pure-wedge order-1 meshes are
accepted. The 'auto' dispatch now routes mesh-backed pure HEX/WEDGE soft irons into HDiv-VIM as well;
mixed / pyramid meshes fail loud until their HDiv support lands.

This locks:
  * hex LINEAR: exact cube demag 1/3 and transverse ~0;
  * hex NONLINEAR (BH table): converges with bounded iterations;
  * wedge: direct HDiv-VIM solve works; mixed / pyramid remain outside this path.

The affine extraction contains an internal finite/consecutive-sample guard for
NGSolve first-touch state.  This gate deliberately runs once: a failure must
remain visible instead of being masked by a whole-solve retry.
"""
import math

import numpy as np
import pytest

pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402
import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia.vim import Solve, MeshSoftIron  # noqa: E402
from radia.vim._vim import (  # noqa: E402
    _Q2_LATTICE_2D,
    _Q2_LATTICE_3D,
    _hex_q2_lattice_nodes_ngsolve_linear,
    _quad_q2_lattice_nodes_ngsolve_linear,
    _trafo_lattice_nodes,
)

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


def test_ngsolve_vol_linear_hex_lattice_matches_gettrafo():
    """The flat-hex fast path is NGSolve .vol reference order, not Cubit/GMSH vertex order."""
    mesh = _cube(2)
    ir_hex = ng.IntegrationRule(_Q2_LATTICE_3D, [1.0] * 27)
    ir_quad = ng.IntegrationRule(_Q2_LATTICE_2D, [1.0] * 9)
    max_hex = 0.0
    for i in range(mesh.GetNE(ng.VOL)):
        e = ng.ElementId(ng.VOL, i)
        err = np.max(np.abs(
            _hex_q2_lattice_nodes_ngsolve_linear(mesh, e) - _trafo_lattice_nodes(mesh, e, ir_hex)))
        max_hex = max(max_hex, float(err))
    max_quad = 0.0
    for i in range(mesh.GetNE(ng.BND)):
        e = ng.ElementId(ng.BND, i)
        err = np.max(np.abs(
            _quad_q2_lattice_nodes_ngsolve_linear(mesh, e) - _trafo_lattice_nodes(mesh, e, ir_quad)))
        max_quad = max(max_quad, float(err))
    assert max_hex < 1e-14
    assert max_quad < 1e-14


def _hdiv_hex_run(mesh, **kw):
    """Run one production solve; internal geometry stabilization is fail-loud."""
    with ng.TaskManager():
        return Solve(mesh, **kw)


def test_vim_solve_hex_linear():
    """vim.Solve on a pure-hex cube: exact demag 1/3 and transverse ~0."""
    res = _hdiv_hex_run(_cube(6), mu_r=MU_R, H_ext=ng.CoefficientFunction((0, 0, H0)))
    assert res["nonlinear"] is False
    assert abs(res["demag"] - 1.0 / 3.0) < 5e-3, f"hex cube demag {res['demag']} off 1/3"
    mz = res["M_avg"][2]
    assert abs(res["M_avg"][0]) < 1e-2 * mz and abs(res["M_avg"][1]) < 1e-2 * mz, "transverse M not ~0"
    assert res["iters"] < 100, f"hex linear CG not mesh-robust: {res['iters']}"
    assert 2.5e3 < mz < 3.5e3, f"hex HDiv Mz {mz:.1f} outside cube demag band"


def test_vim_solve_hex_nonlinear():
    """vim.Solve on a pure-hex cube with a BH table: energy-Newton converges."""
    res = _hdiv_hex_run(_cube(6), bh_table=BH, H_ext=ng.CoefficientFunction((0, 0, H0)))
    assert res["nonlinear"] is True
    assert res["preconditioner_requested"] == "auto"
    assert res["preconditioner"] == "jacobi"
    assert res["preconditioner_policy"] == "auto:hex-wedge-energy-newton-jacobi"
    assert res["iters"] < 100, f"hex energy-Newton not bounded: {res['iters']}"
    mz = res["M_avg"][2]
    assert 1.0e3 < mz < 1.0e4, f"hex nonlinear HDiv Mz {mz:.1f} outside expected band"


def test_vim_solve_wedge_linear():
    """A prism/wedge (6-vertex) mesh now SOLVES via the C++ wedge-mode charge Gram. The C++ wedge Gram is
    eig(M_mass^-1 N) in [0,1] (0.992/0.998 @ n=2/3 in the de-risk)."""
    try:
        mesh = MakeStructured3DMesh(prism=True, nx=2, ny=2, nz=2, mapping=_mp)
    except TypeError:
        pytest.skip("this NGSolve MakeStructured3DMesh has no prism= kwarg")
    verts = {len(el.vertices) for el in mesh.Elements(ng.VOL)}
    if verts != {6}:
        pytest.skip(f"prism mesh did not produce pure wedges (got {sorted(verts)})")
    res = _hdiv_hex_run(mesh, mu_r=MU_R, H_ext=ng.CoefficientFunction((0, 0, H0)))   # H field, A/m
    mz = res["M_avg"][2]
    assert abs(res["demag"] - 1.0 / 3.0) < 6e-2, f"wedge cube demag {res['demag']} off 1/3"
    assert 2.0e3 < mz < 4.0e3, f"wedge HDiv Mz {mz:.1f} outside cube demag band"
    assert res["nonlinear"] is False


def test_rad_solve_demag_backend_hdiv_on_hex():
    """The wrapper path: rad.Solve(demag_backend='hdiv') on a HEX soft_iron routes through _solve_via_hdiv
    -> vim.Solve and SOLVES (no TET-only raise), writing per-element M back to the iron handles."""
    rad.UtiDelAll()
    with ng.TaskManager():
        iron = MeshSoftIron(_cube(4), mu_r=MU_R)
    src = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])
    top = rad.ObjCnt([iron, src])
    with ng.TaskManager():
        rad.Solve(top, 1e-6, 3000, 0, demag_backend="hdiv")    # explicit hdiv on hex -> now solves
    M = np.array([m for (_c, m) in rad.ObjM(iron)], float)             # hdiv-solved M written back via ObjSetM
    assert M.shape[0] == 64, f"expected 64 hex, got {M.shape[0]}"
    mz = float(M[:, 2].mean())
    assert mz > 1e2, f"hdiv-solved hex Mz not substantial: {mz}"       # mu_r=100, H0=1000 -> Mz ~3400
    assert abs(float(M[:, 0].mean())) < 1e-2 * mz and abs(float(M[:, 1].mean())) < 1e-2 * mz
    rad.UtiDelAll()


def test_rad_solve_auto_on_hex_uses_hdiv():
    """The wrapper path with the production default: rad.Solve(auto) on HEX MeshSoftIron returns the
    HDiv-VIM result dict and writes per-element M back."""
    rad.UtiDelAll()
    rad.set_demag_backend("auto")
    with ng.TaskManager():
        iron = MeshSoftIron(_cube(3), mu_r=MU_R)
    src = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])
    top = rad.ObjCnt([iron, src])
    with ng.TaskManager():
        res = rad.Solve(top, 1e-6, 3000, 0)

    assert isinstance(res, dict), type(res)
    assert res["n_el"] == 27
    assert res["order"] == 1
    M = np.array([m for (_c, m) in rad.ObjM(iron)], float)
    assert M.shape[0] == 27
    assert float(M[:, 2].mean()) > 1e2
    rad.UtiDelAll()
    rad.set_demag_backend("auto")


def test_rad_solve_auto_preserves_mesh_soft_iron_rt2_order():
    """The user-intent bridge must not collapse an RT2 request back to RT1."""
    rad.UtiDelAll()
    rad.set_demag_backend("auto")
    with ng.TaskManager():
        iron = MeshSoftIron(_cube(2), mu_r=MU_R, order=2)
    src = rad.ObjBckg(lambda p: [0.0, 0.0, MU0 * H0])
    top = rad.ObjCnt([iron, src])
    with ng.TaskManager():
        res = rad.Solve(top, 1e-7, 3000, 0)

    assert res["order"] == 2
    assert res["n_el"] == 8
    assert np.isfinite(res["M_avg"]).all()
    assert float(res["M_avg"][2]) > 1e2
    rad.UtiDelAll()
    rad.set_demag_backend("auto")
