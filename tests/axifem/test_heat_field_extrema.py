"""Focused regressions for physical temperature-range reporting."""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from netgen.meshing import (
    EdgeDescriptor,
    Element1D,
    Element2D,
    FaceDescriptor,
    MeshPoint,
    Pnt,
)
from netgen.meshing import (
    Mesh as NetgenMesh,
)

from radia.panels.calc_heat import _temperature_extrema


@pytest.fixture(scope="module", autouse=True)
def _taskmanager():
    with ng.TaskManager():
        yield


def _single_quad_mesh():
    netgen_mesh = NetgenMesh(dim=2)
    netgen_mesh.SetMaterial(1, "workpiece")
    netgen_mesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    netgen_mesh.SetBCName(0, "outer")
    edge = EdgeDescriptor()
    edge.edgenr = 1
    edge.surfnr = (1, -1)
    edge.domin = 1
    edge.domout = 0
    edge.name = "outer"
    netgen_mesh.Add(edge)
    points = [
        netgen_mesh.Add(MeshPoint(Pnt(x, y, 0.0)))
        for x, y in ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    ]
    netgen_mesh.Add(Element2D(1, points))
    for first, second in ((0, 1), (1, 2), (2, 3), (3, 0)):
        netgen_mesh.Add(Element1D([points[first], points[second]], index=1))
    return ng.Mesh(netgen_mesh)


def test_q2_extrema_include_nonvertex_field_modes():
    mesh = _single_quad_mesh()
    temperature = ng.GridFunction(ng.H1(mesh, order=2))
    temperature.Set(1.0 - (ng.x - 0.37) ** 2 - (ng.y - 0.41) ** 2)

    raw_coefficients = np.asarray(temperature.vec.FV().NumPy())
    minimum, maximum, metadata = _temperature_extrema(temperature, mesh, 2)

    assert maximum > 0.99
    assert maximum > float(np.max(raw_coefficients)) + 0.05
    assert minimum == pytest.approx(1.0 - 0.63**2 - 0.59**2)
    assert metadata["method"] == (
        "vertices-and-volume-boundary-integration-points"
    )
    assert metadata["integration_order"] >= 6
    assert metadata["sample_count"] > len(mesh.vertices)


def test_order1_extrema_use_exact_vertex_range():
    mesh = _single_quad_mesh()
    temperature = ng.GridFunction(ng.H1(mesh, order=1))
    temperature.Set(2.0 + 3.0 * ng.x - ng.y)

    minimum, maximum, metadata = _temperature_extrema(temperature, mesh, 1)

    assert minimum == pytest.approx(1.0)
    assert maximum == pytest.approx(5.0)
    assert metadata == {
        "method": "vertices",
        "integration_order": 0,
        "sample_count": 4,
    }
