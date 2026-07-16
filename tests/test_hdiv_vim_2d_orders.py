"""BDM1/Q2 and BDM2/Q3 contracts for the production planar HDiv-VIM path."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.occ import OCCGeometry, WorkPlane  # noqa: E402
from netgen.meshing import (  # noqa: E402
    Element1D, Element2D, FaceDescriptor, Mesh as NetgenMesh,
    MeshPoint, Pnt,
)

from radia import vim  # noqa: E402


def _ellipse_mesh(curve_order):
    mesh = ng.Mesh(
        OCCGeometry(WorkPlane().Ellipse(0.2, 0.1).Face(), dim=2)
        .GenerateMesh(maxh=0.1))
    mesh.Curve(curve_order)
    return mesh


def _mixed_tri_quad_mesh():
    mesh = NetgenMesh(dim=2)
    mesh.SetMaterial(1, "body")
    mesh.Add(FaceDescriptor(surfnr=1, domin=0, bc=1))
    mesh.SetBCName(0, "outer")
    points = [
        mesh.Add(MeshPoint(Pnt(x, y, 0.0)))
        for x, y in ((0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1))
    ]
    mesh.Add(Element2D(1, [points[i] for i in (0, 1, 4, 3)]))
    mesh.Add(Element2D(1, [points[i] for i in (1, 2, 5)]))
    mesh.Add(Element2D(1, [points[i] for i in (1, 5, 4)]))
    for first, second in ((0, 1), (1, 2), (2, 5), (5, 4), (4, 3), (3, 0)):
        mesh.Add(Element1D([points[first], points[second]], index=1))
    return ng.Mesh(mesh)


def test_planar_rt2_q3_matches_ellipse_demag_and_native_field():
    mesh = _ellipse_mesh(3)
    with ng.TaskManager():
        result = vim.Solve(
            mesh, order=2, mu_r=1000.0,
            H_ext=ng.CoefficientFunction((1.0e5, 0.0)), tol=1e-11)
        field = result["body"].H_at(
            np.array([[0.3, 0.0], [0.0, 0.2], [-0.3, 0.1]]), result["m"])

    assert result["order"] == 2
    assert result["geometry_order"] == 3
    np.testing.assert_allclose(result["demag_factors"], (1.0/3.0, 2.0/3.0), atol=1e-5)
    assert np.isfinite(field).all()


def test_planar_rt1_and_rt2_agree_on_the_same_q2_mesh():
    mesh = _ellipse_mesh(2)
    applied = ng.CoefficientFunction((1.0e5/np.sqrt(2.0), 1.0e5/np.sqrt(2.0)))
    with ng.TaskManager():
        bdm1 = vim.Solve(mesh, order=1, mu_r=1000.0, H_ext=applied, tol=1e-11)
        bdm2 = vim.Solve(mesh, order=2, mu_r=1000.0, H_ext=applied, tol=1e-11)

    np.testing.assert_allclose(bdm1["M_avg"], bdm2["M_avg"], rtol=5e-5)
    assert bdm2["ndof"] > bdm1["ndof"]


def test_planar_bdm1_rejects_q3_geometry():
    with pytest.raises(ValueError, match="does not support geometry order 3 for 2D tri BDM1"):
        vim.PlanarDemagBody(_ellipse_mesh(3), order=1)


@pytest.mark.parametrize("order", [1, 2])
def test_planar_mixed_tri_quad_mesh_runs_the_advertised_order(order):
    mesh = _mixed_tri_quad_mesh()
    with ng.TaskManager():
        result = vim.Solve(
            mesh, order=order, mu_r=10.0,
            H_ext=ng.CoefficientFunction((100.0, 0.0)), tol=1.0e-10)

    assert {len(element.vertices) for element in mesh.Elements(ng.VOL)} == {3, 4}
    assert result["order"] == order
    assert np.isfinite(result["M_avg"]).all()
    assert result["M_avg"][0] > 100.0
