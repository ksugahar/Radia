# -*- coding: utf-8 -*-
"""Regression test for the 'netgen-scale-after-generate' lint rule.

Rescaling a mesh returned by GenerateMesh moves the nodes but not the OCC
geometry saved with it, so a later Mesh.Curve(k>=2) projects onto the
unscaled CAD. The rule flags that without flagging scaled shapes.
"""
from radia_mcp.radia_ngsolve.rules import ALL_RULES, check_netgen_scale_after_generate


def _hits(lines):
    return check_netgen_scale_after_generate("t.py", lines)


def test_scale_after_generate_flagged():
    f = _hits(["ngm = geo.GenerateMesh(maxh=1.8)", "ngm.Scale(1e-3)", "ngm.Save(out)"])
    assert len(f) == 1
    assert f[0]["rule"] == "netgen-scale-after-generate"
    assert f[0]["severity"] == "HIGH"
    assert f[0]["line"] == 2


def test_occgeometry_chain_flagged():
    f = _hits(["ng = OCCGeometry(shape).GenerateMesh(maxh=h)", "ng.Scale(0.001)"])
    assert len(f) == 1


def test_shape_scale_before_meshing_ok():
    f = _hits(["shape = shape.Scale(Pnt(0, 0, 0), 1e-3)",
               "ngm = OCCGeometry(shape).GenerateMesh(maxh=1.8e-3)"])
    assert f == []


def test_comment_ok():
    f = _hits(["ngm = geo.GenerateMesh(maxh=1)", "# ngm.Scale(1e-3) is wrong"])
    assert f == []


def test_registered():
    assert check_netgen_scale_after_generate in ALL_RULES
