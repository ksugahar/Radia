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


def test_same_line_statements_flagged():
    f = _hits(["ngm = geo.GenerateMesh(maxh=1); ngm.Scale(.001)"])
    assert [x["line"] for x in f] == [1]


def test_spaced_constructor_flagged():
    f = _hits(["ngm = OCCGeometry(Sphere(Pnt(0, 0, 0), 10)).GenerateMesh(maxh=5)",
               "ngm.Scale(.001)"])
    assert len(f) == 1


def test_alias_flagged():
    f = _hits(["ngm = geo.GenerateMesh(maxh=1)", "other = ngm", "other.Scale(.001)"])
    assert [x["line"] for x in f] == [3]


def test_rebinding_clears():
    f = _hits(["ngm = geo.GenerateMesh(maxh=1)", "ngm = shape",
               "ngm.Scale(Pnt(0, 0, 0), .001)"])
    assert f == []


def test_string_ok():
    f = _hits(["ngm = geo.GenerateMesh(maxh=1)", '"""ngm.Scale(.001)"""'])
    assert f == []


def test_function_scopes_are_separate():
    f = _hits(["def a():", "    ngm = geo.GenerateMesh(maxh=1)",
               "def b(ngm):", "    ngm.Scale(.001)"])
    assert f == []


def test_scale_in_rebinding_rhs_flagged():
    f = _hits(["ngm = geo.GenerateMesh()", "ngm = ngm.Scale(.001)"])
    assert [x["line"] for x in f] == [2]


def test_conditional_rebinding_still_flagged():
    f = _hits(["ngm = geo.GenerateMesh()", "if condition:", "    ngm = shape",
               "ngm.Scale(.001)"])
    assert [x["line"] for x in f] == [4]


def test_loop_target_clears():
    f = _hits(["ngm = geo.GenerateMesh()", "for ngm in shapes:",
               "    ngm.Scale(Pnt(0, 0, 0), .001)"])
    assert f == []


def test_loop_carried_generation_flagged():
    f = _hits(["m = None", "for k in range(2):", "    if m is not None:",
               "        m.Scale(.001)", "    m = geo.GenerateMesh()"])
    assert [x["line"] for x in f] == [4]


def test_registered():
    assert check_netgen_scale_after_generate in ALL_RULES
