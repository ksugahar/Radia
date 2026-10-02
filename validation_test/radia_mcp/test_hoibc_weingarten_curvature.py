"""The HOIBC recipe's mean-curvature projection matches analytic geometry.

PEEC HOIBC recipe 6 tells agents to project Trace(Weingarten)/2 and compare
it with 1/R on a sphere and 1/(2R) on a cylinder side, with the sign flipped
when the conductor is meshed as a hole.  These checks run that projection.
"""

import pytest

ngsolve = pytest.importorskip("ngsolve")
occ = pytest.importorskip("netgen.occ")

from radia_mcp.peec.hoibc_knowledge import get_hoibc_documentation

R = 0.5
TOL = 0.02  # relative, curved order-3 mesh, surface maxh = R / 5


def _mean_curvature(shape, label, *, hole=False):
    for face in shape.faces:
        face.maxh = R / 5  # refine only the curved surface under test
    geometry = occ.Box(occ.Pnt(-1, -1, -1.5), occ.Pnt(1, 1, 1.5)) - shape if hole else shape
    mesh = ngsolve.Mesh(occ.OCCGeometry(geometry).GenerateMesh(maxh=0.6))
    mesh.Curve(3)
    region = mesh.Boundaries(label)
    space = ngsolve.SurfaceL2(mesh, order=0, definedon=region)
    gf = ngsolve.GridFunction(space)
    curvature = ngsolve.Trace(ngsolve.specialcf.Weingarten(3)) / 2
    gf.Set(curvature, definedon=region)
    area = ngsolve.Integrate(1.0, mesh, ngsolve.BND, definedon=region)
    return ngsolve.Integrate(gf, mesh, ngsolve.BND, definedon=region) / area


def _sphere():
    sphere = occ.Sphere(occ.Pnt(0, 0, 0), R)
    for face in sphere.faces:
        face.name = "conductor_bnd"
    return sphere


def test_sphere_mean_curvature_is_one_over_r():
    assert _mean_curvature(_sphere(), "conductor_bnd") == pytest.approx(1 / R, rel=TOL)


def test_hole_orientation_reverses_the_sign():
    value = _mean_curvature(_sphere(), "conductor_bnd", hole=True)
    assert value == pytest.approx(-1 / R, rel=TOL)


def test_cylinder_side_mean_curvature_is_one_over_two_r():
    cylinder = occ.Cylinder(occ.Pnt(0, 0, -1), occ.Z, r=R, h=2)
    for face in cylinder.faces:
        face.name = "cap"
    max(cylinder.faces, key=lambda face: face.mass).name = "side"
    assert _mean_curvature(cylinder, "side") == pytest.approx(1 / (2 * R), rel=TOL)


def test_recipe_states_the_checked_values():
    text = get_hoibc_documentation("ngsolve_recipes")
    assert "1/(2R) on a cylinder side" in text
    assert "hole gives -1/R" in text
