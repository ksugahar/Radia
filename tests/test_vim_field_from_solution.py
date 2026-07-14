"""Golden tests for the C++ radia.vim.FieldFromSolution RT1 field path.

Locks:
  * exact linear-M coefficient recovery from element-wise moments (BDM1 = full P1);
  * a uniform-M box against Radia's independent analytic cuboid field;
  * the mu_r sphere end-to-end: far field matches the mesh-moment dipole and the
    constant-M-collapse evaluation of the same solution;
  * fail-loud contract (res without gfM / wrong order).
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
rad = pytest.importorskip("radia")

from radia.vim import _field_batch as fb
from radia import vim


def test_linear_coefficient_recovery_exact():
    from netgen.occ import Box, Pnt, OCCGeometry
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Box(Pnt(0.1, -0.2, 0.05), Pnt(0.35, 0.1, 0.3)))
                       .GenerateMesh(maxh=0.12))
        fes = ng.HDiv(mesh, order=1)
        gfM = ng.GridFunction(fes)
        A0 = np.array([3.0, -2.0, 1.5])
        G0 = np.array([[0.5, -1.2, 0.3], [2.0, 0.7, -0.4], [-0.8, 0.9, 1.1]])
        gfM.Set(ng.CoefficientFunction(tuple(
            A0[i] + G0[i, 0] * ng.x + G0[i, 1] * ng.y + G0[i, 2] * ng.z
            for i in range(3))))
        a, G, c, V = fb._linear_M_coefficients(gfM)
    assert np.max(np.abs(a - (A0[None, :] + np.einsum("ij,ej->ei", G0, c)))) < 1e-8
    assert np.max(np.abs(G - G0[None, :, :])) < 1e-7


def test_uniform_box_matches_radia_cpp():
    from netgen.occ import Box, Pnt, OCCGeometry
    rad.UtiDelAll()
    dims = np.array([0.04, 0.05, 0.03])
    center = 0.5 * dims
    M0 = np.array([2.0e5, -1.0e5, 3.0e5])
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(*dims))).GenerateMesh(maxh=0.015))
        gfM = ng.GridFunction(ng.HDiv(mesh, order=1))
        gfM.Set(ng.CoefficientFunction(tuple(M0)))
    box = rad.ObjRecMag(center.tolist(), dims.tolist(), M0.tolist())
    pts = np.array([[0.10, 0.03, 0.02], [-0.05, 0.02, 0.04],
                    [0.02, 0.12, -0.04], [0.08, -0.04, 0.09]])
    H_rad = np.asarray(rad.Fld(box, "h", pts))
    H_rt = vim.FieldFromSolution(
        {"gfM": gfM, "order": 1, "curve_order": None}, pts)
    rel = (np.linalg.norm(H_rt - H_rad, axis=1)
           / np.maximum(np.linalg.norm(H_rad, axis=1), 1e-30))
    mixed = np.array([[0.01, 0.02, 0.01], [0.20, 0.20, 0.20],
                      [0.03, 0.01, 0.02]])
    got_m = fb.magnetization_from_solution(
        {"gfM": gfM, "order": 1, "curve_order": None}, mixed)
    expected_m = np.array([M0, [0.0, 0.0, 0.0], M0])
    rad.UtiDelAll()
    assert rel.max() < 2e-8
    assert np.allclose(got_m, expected_m, rtol=2e-14, atol=2e-10)


def test_sphere_end_to_end_and_fail_loud():
    from netgen.occ import OCCGeometry, Sphere, Pnt
    rng = np.random.default_rng(3)
    rad.UtiDelAll()
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 0.05)).GenerateMesh(maxh=0.02))
    iron = vim.MeshSoftIron(mesh, mu_r=1000.0)
    bkg = rad.ObjBckg(lambda p: [0.0, 0.0, 4e-7 * np.pi * 1000.0])
    top = rad.ObjCnt([iron, bkg])
    res = rad.Solve(top)
    assert "gfM" in res
    assert "_field_evaluator" in res
    evaluator = res["_field_evaluator"]
    assert res["field_evaluator_stats"]["source_kind"] == "analytic-tet"
    with ng.TaskManager():
        V_el = np.asarray(ng.Integrate(ng.CoefficientFunction(1.0), mesh, ng.VOL,
                                       element_wise=True), float)
    m_dip = float((np.asarray(res["M"]) * V_el[:, None]).sum(axis=0)[2])
    u = rng.normal(size=(20, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    far = 0.15 * u
    r = np.linalg.norm(far, axis=1)
    rh = far / r[:, None]
    mz = np.array([0.0, 0.0, m_dip])
    H_an = (3.0 * (rh @ mz)[:, None] * rh - mz[None, :]) / (4.0 * np.pi * r ** 3)[:, None]
    H_lin = vim.FieldFromSolution(res, far)
    assert res["_field_evaluator"] is evaluator
    H_col = np.asarray(rad.Fld(iron, "h", far))
    scale = np.linalg.norm(H_an, axis=1).mean()
    assert np.linalg.norm(H_lin - H_an, axis=1).max() / scale < 2e-2   # facet multipole tail
    assert np.linalg.norm(H_lin - H_col, axis=1).max() / scale < 2e-3  # same solution
    with pytest.raises(ValueError):
        vim.FieldFromSolution({"order": 1}, far)                       # no gfM
    bad = dict(res)
    bad["order"] = 0
    with pytest.raises(NotImplementedError):
        vim.FieldFromSolution(bad, far)
    rad.UtiDelAll()
