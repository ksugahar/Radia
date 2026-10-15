"""Cubit-free load contract for the native exporter's independently authored cube."""
from pathlib import Path

import ngsolve as ng
import pytest

from cubit_mesh_export.check import check_consistency


FIXTURES = Path(__file__).with_name("fixtures")


def test_named_curve_sideset_roundtrips_on_the_declared_cube_edge():
    path = FIXTURES / "bbnd_brick.vol"
    assert path.is_file(), "Required native-export fixture missing; regenerate with bbnd_brick.py"
    mesh = ng.Mesh(str(path))
    names = mesh.GetBBoundaries()
    assert "test_edge" in names
    selected = [element for element in mesh.Elements(ng.BBND)
                if names[element.index] == "test_edge"]
    assert selected
    length = 0.0
    for element in selected:
        a, b = [mesh[vertex].point for vertex in element.vertices]
        for point in (a, b):
            assert point[0] == pytest.approx(-0.5, abs=1e-12)
            assert point[1] == pytest.approx(-0.5, abs=1e-12)
            assert -0.5-1e-12 <= point[2] <= 0.5+1e-12
        length += abs(a[2]-b[2])
    assert length == pytest.approx(1.0, abs=1e-12)
    report = check_consistency(
        path, contract=FIXTURES / "bbnd_brick.contract.json",
        quality=True, require_tetrahedra=True, min_curve_order=2)
    assert report["passed"], report
