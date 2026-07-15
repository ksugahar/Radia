"""Golden: the HDiv-VIM RT1/RT2 production contract.

RT1 supports pure-TET, pure-HEX, pure-WEDGE, 2D, IMA, curved geometry, and field evaluation.
RT2 is the production pure-TET higher-order route, including flat/curved IMA and the persistent C++
field evaluator.  Specialized planar/HEX/WEDGE charge kernels remain RT1.  Pyramid / mixed meshes and
The deleted region-dictionary PM API is absent; fixed M uses an independent
MagnetizationSource HDiv space.

HEX/WEDGE were UNLOCKED 2026-07-04: the wired RT1 charge Grams + the Gram-agnostic solve path already
solve them, so the auto guard is now pure-TET / pure-HEX / pure-WEDGE.  The per-element hex correctness
lock lives in test_hdiv_vim_hex_public_solve.py; wedge spectrum/cube locks live in
test_hdiv_vim_wedge_spectrum.py.

This test LOCKS the production fail-loud boundaries (No-Fallbacks: the raise IS
the feature), the RT1/RT2 happy paths, and the rad.Solve 'auto' split.

NGSolve + Netgen required.  See memory/hdiv_rt1_field_production.md.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Box, Sphere, OCCGeometry, Pnt  # noqa: E402

import radia as rad  # noqa: E402
from radia.vim import Solve, MeshSoftIron, DemagOperator, ChargeGram  # noqa: E402

_HEXT = ng.CoefficientFunction((0, 0, 1e4))


def _sphere(maxh=0.5):
    return ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=maxh))


def _half_box(maxh=2.0):
    return ng.Mesh(OCCGeometry(Box(Pnt(0, -1, -1), Pnt(1, 1, 1))).GenerateMesh(maxh=maxh))


# ---------------------------------------------------------------- RT1 happy path (the SUPPORTED config)
def test_rt1_is_the_default_order_and_solves():
    """The default order is now 1 (RT1); a linear tet solve gives the sphere demag ~1/3."""
    mesh = _sphere()
    with ng.TaskManager():
        r = Solve(mesh, mu_r=100.0, H_ext=_HEXT)     # no order= -> default RT1
    assert r.get("order", 1) == 1
    assert 0.31 < r["demag"] < 0.345, r["demag"]


def test_rt1_nonlinear_solves():
    """Flat RT1 nonlinear (the energy-Newton on the high-order Gram) runs at the default order."""
    BH = [[0.0, 0.0], [500.0, 1.3], [5000.0, 1.85], [5e5, 2.4]]
    mesh = _sphere()
    with ng.TaskManager():
        r = Solve(mesh, bh_table=BH, H_ext=_HEXT)
    assert r["nonlinear"] and r["M_avg"][2] > 1e3


def test_public_rt2_pure_tet_solves():
    """RT2 is a supported public material solve on pure tetrahedra."""
    mesh = _sphere()
    with ng.TaskManager():
        result = Solve(mesh, mu_r=100.0, H_ext=_HEXT, order=2)
    assert result["order"] == 2
    assert 0.31 < result["demag"] < 0.345, result["demag"]


def test_deleted_pm_region_api_has_no_compatibility_alias():
    """Fixed M uses MagnetizationSource; deleted pm_M is absent, not a compatibility alias."""
    import inspect

    from radia.vim import MagnetizationSource

    assert "pm_M" not in inspect.signature(Solve).parameters
    assert MagnetizationSource.__name__ == "MagnetizationSource"


def test_curved_image_symmetry_is_wired():
    """Curve(2) and IMA share the configured C++ solve and field evaluator."""
    mesh = _half_box()
    with ng.TaskManager():
        result = Solve(mesh, mu_r=100.0, H_ext=_HEXT, image="+x", curve_order=2)
    assert result["curve_order"] == 2 and result["image"] == "+x"
    assert result["symmetry_constrained_dofs"] > 0
    assert result["_charge_gram"].constraint_count == 0  # removed by NGSolve Compress before assembly
    assert result["field_evaluator_stats"]["image_count"] == 1


def test_operator_and_gram_support_pure_tet_rt1_rt2():
    """The NGSolve-style operator and charge builder expose both production TET orders."""
    mesh = _sphere()
    with ng.TaskManager():
        values = {}
        for order in (1, 2):
            fes = ng.HDiv(mesh, order=order)
            B, G, M = ChargeGram(fes)
            assert B.shape[1] == fes.ndof and M.shape == (fes.ndof, fes.ndof)
            assert G.ndof() == B.shape[0]
            values[order] = DemagOperator(fes).DemagFactor(ng.CF((0, 0, 1)))
    assert all(0.31 < value < 0.345 for value in values.values()), values
    assert abs(values[2] - 1.0 / 3.0) < abs(values[1] - 1.0 / 3.0), values
    assert abs(values[2] - values[1]) < 2e-3, values


def test_rt2_ima_uses_the_high_order_cpp_operator_and_field_evaluator():
    mesh = _half_box()
    with ng.TaskManager():
        result = Solve(mesh, mu_r=100.0, H_ext=_HEXT, order=2, image="+x")
    assert result["order"] == 2 and result["image"] == "+x"
    assert result["_charge_gram"].operator_configured
    assert result["symmetry_constrained_dofs"] > 0
    assert result["field_evaluator_stats"]["image_count"] == 1


def test_rt2_curved_geometry_is_wired():
    mesh = _half_box()
    mesh.Curve(2)
    with ng.TaskManager():
        result = Solve(mesh, mu_r=100.0, H_ext=_HEXT, order=2)
    assert result["order"] == 2 and result["curve_order"] == 2
    assert result["_charge_gram"].operator_configured
    assert result["field_evaluator_stats"]["source_count"] > 0


# ---------------------------------------------------------------- pure HEX/WEDGE -> HDiv-VIM eligibility
def test_hex_soft_iron_auto_is_hdiv_eligible():
    """rad.Solve 'auto' now routes a HEX mesh-backed soft iron to HDiv-VIM.

    Keep this contract cheap: test the read-only dispatch predicate here; the actual solve is locked in
    test_hdiv_vim_hex_public_solve.py.
    """
    from ngsolve.meshes import MakeStructured3DMesh
    from radia.vim import _radsolve

    rad.UtiDelAll()
    mp = lambda x, y, z: (0.01 * (x - 0.5), 0.01 * (y - 0.5), 0.01 * (z - 0.5))  # noqa: E731
    with ng.TaskManager():
        hexm = MakeStructured3DMesh(hexes=True, nx=3, ny=3, nz=3, mapping=mp)
    with pytest.raises(NotImplementedError, match="RT2.*pure-TET"):
        Solve(hexm, mu_r=100.0, H_ext=_HEXT, order=2)
    iron = MeshSoftIron(hexm, mu_r=1000.0)
    top = rad.ObjCnt([iron])
    assert _radsolve.is_hdiv_eligible(top)
    rad.UtiDelAll()


def test_vim_solve_accepts_wedge_mesh_directly():
    """vim.Solve now accepts pure-tet, pure-hex, AND pure-WEDGE/prism (6-vertex) via the C++
    wedge-mode charge Gram (2026-07-04, memory hdiv-tet-hex-coupling-pyramid-gated).  A prism-meshed cube
    solves and returns the ~1/3 cube demag factor (the C++ wedge Gram is eig(M_mass^-1 N) in [0,1]:
    0.992/0.998 @ n=2/3 in the de-risk)."""
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
        out = Solve(prism, mu_r=100.0, H_ext=_HEXT)
    # cube demag factor ~ 1/3 (geometric, mu-independent); the C++ wedge Gram is PSD + eig<=1
    assert 0.30 <= out["demag"] <= 0.36, out["demag"]
    assert out["n_el"] == 16 and out["nonlinear"] is False
