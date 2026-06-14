# -*- coding: utf-8 -*-
"""Regression test for the 'ngsolve-curve-after-vol-import' lint rule.

A .vol exported with a curve order (e.g. Cubit `export netgen ... order N`) already
stores its high-order node positions. Calling mesh.Curve() after ngsolve.Mesh(".vol")
recomputes curved nodes from CAD geometry the imported mesh does not have, dropping
the baked-in curving -> the surface flattens to facets (~1% floor). The rule flags
that footgun WITHOUT flagging the legitimate OCC/GenerateMesh().Curve() case.
"""
from radia_mcp.radia_ngsolve.rules import (
    check_curve_after_vol_import,
    ALL_RULES,
)


def _hits(lines):
    return check_curve_after_vol_import("t.py", lines)


# ── buggy patterns: must be flagged ───────────────────────────
def test_separate_literal_vol_then_curve():
    f = _hits(['mesh = ng.Mesh("C:/temp/full_hex_o2.vol")', "mesh.Curve(2)"])
    assert len(f) == 1
    assert f[0]["rule"] == "ngsolve-curve-after-vol-import"
    assert f[0]["severity"] == "HIGH"
    assert f[0]["line"] == 2


def test_chained_literal_vol_curve():
    f = _hits(['mesh = ng.Mesh("oct.vol").Curve(3)'])
    assert len(f) == 1
    assert f[0]["line"] == 1


def test_curve_via_path_variable():
    f = _hits(['p = r"C:\\temp\\oct.vol"', "mesh = ngsolve.Mesh(p)", "mesh.Curve(2)"])
    assert len(f) == 1
    assert f[0]["line"] == 3


def test_load_then_curve():
    f = _hits(["ngm = Mesh()", 'ngm.Load("foo.vol")', "ngm.Curve(2)"])
    assert len(f) == 1


# ── legitimate patterns: must NOT be flagged ──────────────────
def test_occ_chained_curve_ok():
    f = _hits([
        "mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0,0,0),1)).GenerateMesh(maxh=0.3)).Curve(3)",
    ])
    assert f == []


def test_occ_separate_curve_ok():
    f = _hits(["mesh = ng.Mesh(geo.GenerateMesh(maxh=0.5))", "mesh.Curve(2)"])
    assert f == []


def test_comment_mentioning_pattern_ok():
    f = _hits(['# do not call mesh.Curve() after Mesh("x.vol")',
               "mesh = ng.Mesh(geo.GenerateMesh(maxh=0.5))"])
    assert f == []


def test_rule_is_registered():
    assert check_curve_after_vol_import in ALL_RULES
