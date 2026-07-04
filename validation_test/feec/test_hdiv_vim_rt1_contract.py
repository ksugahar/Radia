"""Golden: the HDiv-VIM TET / RT1-only CONTRACT (Sugahara 2026-06-29).

HDiv-VIM is the RT1 (HDiv order=1) high-order element of the soft-iron demag stack, on a pure-TET or
pure-HEX mesh.  RT0 is retired (per-element INACCURATE -- the demag factor is right ~1/3 but the
per-element M leaks; RT1 is what fixes it); RT2+ is retired (no per-element gain over RT1, slower);
wedge / pyramid / mixed meshes, pm_M mixing, IMA image symmetry, the 'gauss' point Gram and the 'hlu'
system-A solver are all retired from / not in HDiv-VIM -> the collocation MMMM backend.

HEX was UNLOCKED 2026-07-04: the wired hex RT1 charge Gram + the Gram-agnostic energy-Newton already
solve it, so the tet-only guard was relaxed to pure-TET-OR-pure-HEX.  rad.Solve's 'auto' split STILL
routes a HEX mesh-backed iron to collocation MMMM by DEFAULT (KEEP-BOTH; MMMM = coarse tier), but an
explicit demag_backend='hdiv' -- and a direct hdiv_demag_solve(hexmesh, ...) call -- now solve pure hex.
The per-element hex correctness lock lives in test_hdiv_vim_hex_public_solve.py.

This test LOCKS the fail-loud retirement (No-Fallbacks: the raise IS the feature -- it names the backend to
use instead) + the RT1 happy path + the rad.Solve 'auto' split.  It replaces the per-feature golden tests
for the retired paths (gauss / hlu / pm / image / wedge).

NGSolve + Netgen required.  See memory/hdiv_vim_tet_rt1_only.md.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Sphere, OCCGeometry, Pnt  # noqa: E402

import radia as rad  # noqa: E402
from radia.vim import hdiv_demag_solve, soft_iron_from_mesh, DemagOperator, build_charge_gram  # noqa: E402

_HEXT = ng.CoefficientFunction((0, 0, 1e4))


def _sphere(maxh=0.5):
    return ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=maxh))


# ---------------------------------------------------------------- RT1 happy path (the SUPPORTED config)
def test_rt1_is_the_default_order_and_solves():
    """The default order is now 1 (RT1); a linear tet solve gives the sphere demag ~1/3."""
    mesh = _sphere()
    with ng.TaskManager():
        r = hdiv_demag_solve(mesh, mu_r=100.0, H_ext=_HEXT)     # no order= -> default RT1
    assert r.get("order", 1) == 1
    assert 0.31 < r["demag"] < 0.345, r["demag"]


def test_rt1_nonlinear_solves():
    """Flat RT1 nonlinear (the energy-Newton on the high-order Gram) runs at the default order."""
    BH = [[0.0, 0.0], [500.0, 1.3], [5000.0, 1.85], [5e5, 2.4]]
    mesh = _sphere()
    with ng.TaskManager():
        r = hdiv_demag_solve(mesh, bh_table=BH, H_ext=_HEXT)
    assert r["nonlinear"] and r["M_avg"][2] > 1e3


# ---------------------------------------------------------------- retired-from-HDiv-VIM (fail-loud)
def test_rt0_retired():
    """order=0 (RT0) is retired -- per-element inaccurate; the error names the collocation MMMM alternative."""
    mesh = _sphere()
    with pytest.raises(ValueError, match="RT1"):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r=100.0, H_ext=_HEXT, order=0)


def test_rt2_retired():
    """order>1 (RT2+) is retired (no per-element gain over RT1)."""
    mesh = _sphere()
    with pytest.raises(ValueError, match="RT1"):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r=100.0, H_ext=_HEXT, order=2)


def test_pm_mixing_retired():
    """pm_M (permanent-magnet mixing) is retired -- HDiv-VIM is a soft-iron solver."""
    mesh = _sphere()
    with pytest.raises(NotImplementedError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r=100.0, H_ext=_HEXT, pm_M={"default": [0, 0, 1.0]})


def test_image_symmetry_retired():
    """IMA mirror symmetry (image) is retired -- reduced models use collocation MMMM."""
    mesh = _sphere()
    with pytest.raises(NotImplementedError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r=100.0, H_ext=_HEXT, image="+x")


def test_gauss_and_hlu_backends_retired():
    """The 'gauss' point Gram and the 'hlu' system-A solver were RT0 experiments -- both retired."""
    mesh = _sphere()
    with pytest.raises(ValueError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r=100.0, H_ext=_HEXT, gram_backend="gauss")
    with pytest.raises(ValueError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r=100.0, H_ext=_HEXT, linear_solver="hlu")


def test_operator_and_gram_are_rt1_tet_only():
    """The ngsolve.bem-style DemagOperator + build_charge_gram reject order!=1 (RT0/RT2 retired)."""
    mesh = _sphere()
    with ng.TaskManager():
        with pytest.raises(ValueError, match="RT1"):
            DemagOperator(ng.HDiv(mesh, order=0))
        with pytest.raises(ValueError, match="RT1"):
            build_charge_gram(ng.HDiv(mesh, order=0))
        # RT1 is accepted
        D = DemagOperator(ng.HDiv(mesh, order=1)).DemagFactor(ng.CF((0, 0, 1)))
    assert 0.31 < D < 0.345, D


# ---------------------------------------------------------------- non-tet -> collocation MMMM
def test_hex_soft_iron_auto_routes_to_collocation_mmmm():
    """rad.Solve 'auto' routes a HEX mesh-backed soft iron to collocation MMMM by DEFAULT (KEEP-BOTH;
    MMMM = coarse tier) -- it solves, returning the C++ tuple.  (Explicit demag_backend='hdiv' on hex now
    SOLVES via the wired hex RT1 Gram; that + the per-element correctness are locked in
    test_hdiv_vim_hex_public_solve.py, kept out of this contract test to avoid a heavy/flaky hex build.)"""
    from ngsolve.meshes import MakeStructured3DMesh
    rad.UtiDelAll()
    mp = lambda x, y, z: (0.01 * (x - 0.5), 0.01 * (y - 0.5), 0.01 * (z - 0.5))  # noqa: E731
    with ng.TaskManager():
        hexm = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3, mapping=mp)
    iron = soft_iron_from_mesh(hexm, mu_r=1000.0)
    src = rad.ObjBckg(lambda p: [0, 0, 0.1])
    top = rad.ObjCnt([iron, src])
    res = rad.Solve(top, 1e-4, 1000, 0)                          # auto -> hex -> collocation MMMM (C++)
    assert isinstance(res, tuple), type(res)                     # C++ collocation MMMM returns a tuple
    rad.UtiDelAll()


def test_hdiv_demag_solve_accepts_wedge_mesh_directly():
    """hdiv_demag_solve now accepts pure-tet, pure-hex, AND pure-WEDGE/prism (6-vertex) via the C++
    wedge-mode charge Gram (2026-07-04, memory hdiv-tet-hex-coupling-pyramid-gated).  A prism-meshed cube
    solves and returns the ~1/3 cube demag factor (the C++ wedge Gram is eig(M_mass^-1 N) in [0,1]:
    0.992/0.998 @ n=2/3 in the de-risk).  A MIXED / pyramid mesh still routes to collocation MMMM."""
    from ngsolve.meshes import MakeStructured3DMesh
    mp = lambda x, y, z: (0.01 * (x - 0.5), 0.01 * (y - 0.5), 0.01 * (z - 0.5))  # noqa: E731
    try:
        with ng.TaskManager():
            prism = MakeStructured3DMesh(prism=True, nx=2, ny=2, nz=2, mapping=mp)
    except TypeError:
        pytest.skip("this NGSolve MakeStructured3DMesh has no prism= kwarg")
    verts = {len(el.vertices) for el in prism.Elements(ng.VOL)}
    if verts != {6}:
        pytest.skip(f"prism mesh did not produce pure wedges (got {sorted(verts)})")
    with ng.TaskManager():
        out = hdiv_demag_solve(prism, mu_r=100.0, H_ext=_HEXT)
    # cube demag factor ~ 1/3 (geometric, mu-independent); the C++ wedge Gram is PSD + eig<=1
    assert 0.30 <= out["demag"] <= 0.36, out["demag"]
    assert out["n_el"] == 16 and out["nonlinear"] is False
