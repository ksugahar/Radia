"""HDiv-VIM production contract: ``rad.Solve`` write-back makes ``rad.Fld`` quantitative.

These tests sit above the pure VIM math tests.  They lock the Radia-facing
contract that matters to applications: solve the mesh-backed soft iron and
evaluate the full RT1 charge field through ``rad.Fld``.  Per-element constant-M
write-back remains metadata/visualization only and is not a field oracle.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

import radia as rad  # noqa: E402
from radia import vim  # noqa: E402

MU0 = 4.0e-7 * math.pi
A = 0.01
MU_R = 1000.0
H0 = 1.0e4
PROBES = np.array([
    [0.0, 0.0, 3.0 * A],
    [1.4 * A, -0.6 * A, 2.5 * A],
    [-1.5 * A, 0.7 * A, -2.4 * A],
], float)


def _hexbox(x0, x1, y0, y1, z0, z1, nx, ny, nz):
    return MakeStructured3DMesh(
        hexes=True, nx=nx, ny=ny, nz=nz,
        mapping=lambda X, Y, Z: (
            x0 + (x1 - x0) * X,
            y0 + (y1 - y0) * Y,
            z0 + (z1 - z0) * Z,
        ),
    )


def _field(obj, pts=PROBES):
    return np.array([rad.Fld(obj, "b", p.tolist()) for p in pts], float)


def _radia_container_from_mesh_M(mesh, M):
    objs = []
    for el, m in zip(mesh.Elements(ng.VOL), M):
        verts = [[float(c) for c in mesh.vertices[v.nr].point] for v in el.vertices]
        mag = [float(c) for c in m]
        if len(verts) == 4:
            objs.append(rad.ObjTetrahedron(verts, mag))
        elif len(verts) == 6:
            objs.append(rad.ObjWedge(verts, mag))
        elif len(verts) == 8:
            objs.append(rad.ObjHexahedron(verts, mag))
        else:
            raise AssertionError(f"unsupported element with {len(verts)} vertices")
    return rad.ObjCnt(objs)


def test_radsolve_hdiv_radfld_uses_full_rt1_solution():
    """rad.Fld on a solved MeshSoftIron is exactly vim.FieldFromSolution, not constant-M write-back."""
    rad.UtiDelAll()
    from radia.vim import _radsolve
    _radsolve.clear_registry()
    with ng.TaskManager():
        mesh = _hexbox(-A, A, -A, A, -A, A, 3, 3, 3)
        iron = vim.MeshSoftIron(mesh, mu_r=MU_R)
        top = rad.ObjCnt([iron, rad.ObjBckg(lambda _p: [0.0, 0.0, MU0 * H0])])
        res = rad.Solve(top, 1e-6, 2000, 0)
        bcollapsed = _field(_radia_container_from_mesh_M(mesh, res["M"]))
        bgot = _field(iron)
        bref = MU0 * vim.FieldFromSolution(res, PROBES)
        m_from_cells = np.mean(np.asarray(res["M"], float), axis=0)
    rel = np.linalg.norm(bgot - bref) / max(np.linalg.norm(bref), 1e-30)
    assert rel < 1e-13, f"rad.Fld(HDiv) != C++ RT1 field (rel {rel:.2e})"
    assert np.all(np.isfinite(bcollapsed))
    assert np.allclose(m_from_cells, res["M_avg"], rtol=1e-13, atol=1e-9)
    assert "field_contract" in res
    rad.UtiDelAll(); _radsolve.clear_registry()


def test_radsolve_hdiv_image_radfld_is_repeatable_rt1_field():
    """An image solve evaluates the reflected RT1 solution and is repeatable to roundoff."""
    rad.UtiDelAll()
    from radia.vim import _radsolve
    _radsolve.clear_registry()
    with ng.TaskManager():
        half_mesh = _hexbox(-A, A, -A, A, 0.0, A, 2, 2, 1)
        half_iron = vim.MeshSoftIron(half_mesh, mu_r=MU_R)
        top = rad.ObjCnt([half_iron, rad.ObjBckg(lambda _p: [0.0, 0.0, MU0 * H0])])
        half_res = rad.Solve(top,
                             1e-6, 2000, 0, image="-z")
        b_half = _field(half_iron)
        b_direct = MU0 * vim.FieldFromSolution(half_res, PROBES)
        again = rad.Solve(top, 1e-6, 2000, 0, image="-z")
        b_again = _field(half_iron)
    assert np.array_equal(b_half, b_direct) or np.allclose(b_half, b_direct, rtol=1e-14, atol=1e-18)
    assert np.array_equal(b_half, b_again) or np.allclose(b_half, b_again, rtol=1e-14, atol=1e-18)
    assert half_res["image"] == "-z"
    assert "M_avg_reduced" in half_res
    assert np.allclose(np.mean(np.asarray(half_res["M"], float), axis=0), half_res["M_avg_reduced"],
                       rtol=1e-13, atol=1e-9)
    rad.UtiDelAll(); _radsolve.clear_registry()


def test_hdiv_image_demag_matches_explicit_full_to_roundoff():
    """The IMA Gram/energy itself equals the explicit full-domain Gram to roundoff on the matching mesh."""
    Hz = ng.CoefficientFunction((0.0, 0.0, H0))
    with ng.TaskManager():
        full = vim.Solve(_hexbox(-A, A, -A, A, -A, A, 1, 1, 2), mu_r=MU_R, H_ext=Hz,
                         gram_eps=1e-12, tol=1e-12)
        half = vim.Solve(_hexbox(-A, A, -A, A, 0.0, A, 1, 1, 1), mu_r=MU_R, H_ext=Hz,
                         image="-z", gram_eps=1e-12, tol=1e-12)
    assert abs(half["demag"] - full["demag"]) < 10.0 * np.finfo(float).eps


def test_hdiv_image_radfld_matches_unconstrained_explicit_full_to_roundoff():
    """Target contract: ``image=`` field and an explicit full solve should agree to ~10 eps.

    The matching full and half meshes use the same affine hex reference rule.  This locks the stronger
    application-level contract rather than accepting a tolerance for a reflection defect in ChargeGram.
    """
    rad.UtiDelAll()
    from radia.vim import _radsolve
    _radsolve.clear_registry()
    with ng.TaskManager():
        full_mesh = _hexbox(-A, A, -A, A, -A, A, 1, 1, 2)
        full_iron = vim.MeshSoftIron(full_mesh, mu_r=MU_R)
        full_top = rad.ObjCnt([full_iron, rad.ObjBckg(lambda _p: [0.0, 0.0, MU0 * H0])])
        rad.Solve(full_top, 1e-6, 2000, 0)
        b_full = _field(full_iron)

        half_mesh = _hexbox(-A, A, -A, A, 0.0, A, 1, 1, 1)
        half_iron = vim.MeshSoftIron(half_mesh, mu_r=MU_R)
        half_top = rad.ObjCnt([half_iron, rad.ObjBckg(lambda _p: [0.0, 0.0, MU0 * H0])])
        rad.Solve(half_top, 1e-6, 2000, 0, image="-z")
        b_half = _field(half_iron)

    rel = np.linalg.norm(b_half - b_full) / max(np.linalg.norm(b_full), 1e-30)
    assert rel < 10.0 * np.finfo(float).eps


def test_hdiv_multicell_hex_image_field_matches_full_to_roundoff():
    """A multi-cell x reflection must not depend on directed FAR quadrature.

    This geometry reaches separated cell/face blocks, unlike the minimal z-split
    contract above.  A one-sided FAR block is matrix-symmetric after upper-triangle
    mirroring but is not reflection-invariant at finite quadrature order.
    """
    H = ng.CoefficientFunction((0.0, 0.0, H0))
    with ng.TaskManager():
        full = vim.Solve(
            _hexbox(-A, A, -A, A, -A, A, 2, 2, 2),
            mu_r=MU_R, H_ext=H, gram_eps=1e-14, leaf=256,
            tol=1e-14, maxit=10000,
        )
        half = vim.Solve(
            _hexbox(0.0, A, -A, A, -A, A, 1, 2, 2),
            mu_r=MU_R, H_ext=H, image="+x", gram_eps=1e-14, leaf=256,
            tol=1e-14, maxit=10000,
        )
    full_field = vim.FieldFromSolution(full, PROBES, algorithm="direct")
    half_field = vim.FieldFromSolution(half, PROBES, algorithm="direct")
    relative = np.linalg.norm(full_field - half_field) / max(np.linalg.norm(full_field), 1e-30)
    assert relative < 10.0 * np.finfo(float).eps
    assert abs(full["demag"] - half["demag"]) < 10.0 * np.finfo(float).eps
    assert half["symmetry_constrained_dofs"] == 16
    assert full["hmat_stats"]["hex_far_one_sided_threshold"] == 0.0
