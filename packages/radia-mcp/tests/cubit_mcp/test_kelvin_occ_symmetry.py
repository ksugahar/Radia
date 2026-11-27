"""Guards: add_kelvin_occ keeps the positive side of each symmetry plane, so the
Kelvin sphere is a translated copy of the air domain and Netgen meshes the
periodic identification.

Assumes: Netgen OCC is importable; coarse meshes only (a few hundred tets).
A failure here means kelvin_int and kelvin_ext no longer pair, which Netgen
would otherwise mesh silently without any identification.
"""
from __future__ import annotations

import math

import pytest

occ = pytest.importorskip("netgen.occ")
ngsolve = pytest.importorskip("ngsolve")

from radia_mcp.cubit.kelvin import add_kelvin_occ  # noqa: E402

R = 0.2


def _air(lo):
    quad = occ.Box(occ.Pnt(*lo), occ.Pnt(*(c + 2.0 for c in lo)))
    air = occ.Sphere(occ.Pnt(0, 0, 0), R) * quad
    air.name = "air"
    return air


@pytest.mark.parametrize("symmetry, lo", [
    (["z"], (-1.0, -1.0, 0.0)),
    (["x", "y"], (0.0, 0.0, -1.0)),
    (["x"], (0.0, -1.0, -1.0)),
])
def test_kelvin_sphere_is_translated_copy_and_identified(symmetry, lo):
    shape, info = add_kelvin_occ(_air(lo), R, symmetry=symmetry)
    off = info["center"]
    kelvin = [s for s in shape.solids if s.name == "kelvin"]
    air = [s for s in shape.solids if s.name == "air"]
    assert len(kelvin) == 1 and len(air) == 1
    kc, ac = kelvin[0].center, air[0].center
    assert math.dist((kc.x, kc.y, kc.z), (ac.x + off[0], ac.y + off[1], ac.z + off[2])) < 1e-9

    mesh = ngsolve.Mesh(occ.OCCGeometry(shape).GenerateMesh(maxh=0.1))
    pts = [tuple(p.p) for p in mesh.ngmesh.Points()]
    pairs = [(p[0], p[1]) for p in mesh.ngmesh.GetIdentifications()]
    n_int = len({v.nr for el in mesh.Elements(ngsolve.BND)
                 if el.mat == "kelvin_int" for v in el.vertices})
    assert len(pairs) == n_int > 0
    resid = max(min(math.dist(pts[b.nr - 1], tuple(pts[a.nr - 1][k] + off[k] for k in range(3))),
                    math.dist(pts[a.nr - 1], tuple(pts[b.nr - 1][k] + off[k] for k in range(3))))
                for a, b in pairs)
    assert resid < 1e-9


@pytest.mark.parametrize("symmetry, lo", [
    (["x"], (-2.0, -1.0, -1.0)),
    (["x", "y"], (-2.0, -2.0, -1.0)),
])
def test_air_on_negative_side_is_refused(symmetry, lo):
    with pytest.raises(RuntimeError, match="positive side"):
        add_kelvin_occ(_air(lo), R, symmetry=symmetry)


def test_one_eighth_is_refused():
    with pytest.raises(ValueError, match="1/8"):
        add_kelvin_occ(_air((0.0, 0.0, 0.0)), R, symmetry=["x", "y", "z"])
