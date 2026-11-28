"""Guards: check-vol measures volumes, areas and edge lengths of a curved .vol
with quadrature fine enough for its geometry order, so the CAD gate compares
the mesh, not a quadrature error.

A curved order-5 cylinder's volume read at NGSolve's default integration order
(5) is off by about 1e-6 relative; the reference here is the same integral at
order 30, converged to about 1e-14.

Assumes: NGSolve with netgen.occ is importable.
"""
import math

import pytest

pytest.importorskip("ngsolve")
from netgen.occ import Axes, Cylinder, OCCGeometry, Pnt, X, Z  # noqa: E402
from ngsolve import BBND, BND, CF, Integrate, Mesh, TaskManager  # noqa: E402

from radia.vol_check import check_consistency  # noqa: E402

ORDER = 5
REFERENCE_ORDER = 30


def _curved_cylinder(path):
    body = Cylinder(Axes(Pnt(0, 0, 0), n=Z, h=X), r=0.05, h=0.04)
    body.mat("body")
    body.faces.name = "skin"
    body.edges.name = "rim"
    mesh = Mesh(OCCGeometry(body).GenerateMesh(maxh=0.04))
    mesh.Curve(ORDER)
    mesh.ngmesh.Save(str(path))
    return Mesh(str(path))


def test_curved_measures_match_converged_quadrature(tmp_path):
    vol = tmp_path / "cylinder.vol"
    with TaskManager():
        mesh = _curved_cylinder(vol)
        assert mesh.GetCurveOrder() == ORDER, mesh.GetCurveOrder()
        expected = {
            ("materials", "body", "ng_volume"): Integrate(CF(1), mesh, order=REFERENCE_ORDER),
            ("boundaries", "skin", "ng_area"): Integrate(CF(1), mesh, BND, order=REFERENCE_ORDER),
            ("edges", "rim", "ng_length"): Integrate(CF(1), mesh, BBND, order=REFERENCE_ORDER),
        }
    report = check_consistency(str(vol), quality=False)
    errors = {}
    for (group, name, key), reference in expected.items():
        entry = next(e for e in report[group] if e["name"] == name)
        errors[f"{group}.{name}.{key}"] = (entry[key], reference, entry[key] / reference - 1.0)
    bad = {k: v for k, v in errors.items() if not math.isclose(v[0], v[1], rel_tol=1e-11)}
    assert not bad, f"(actual, expected, relative error): {bad}"
