# -*- coding: utf-8 -*-
"""Guards: the 'radia-vol-mesh-into-curved-route' lint rule.

Assumes: a plain Mesh(".vol") of a CAD-free file borrows the process geometry, and
a Cubit export cannot be curved after export; radia.mesh_curve.load_vol is the
supported load (docs/ngsolve_integration/curve_order.md).
"""
from radia_mcp.radia_ngsolve.rules import ALL_RULES, check_vol_mesh_into_curved_route


def _hits(lines):
    return check_vol_mesh_into_curved_route("t.py", lines)


def test_plain_vol_mesh_into_curved_solve_is_flagged():
    found = _hits(['mesh = ng.Mesh("iron.vol")', "res = vim.Solve(mesh, order=2, curve_order=2)"])
    assert [(f["line"], f["rule"], f["severity"]) for f in found] == [
        (2, "radia-vol-mesh-into-curved-route", "HIGH")], found


def test_path_variable_and_str_wrapper_are_followed():
    found = _hits(['p = "meshes/coil.vol.gz"', "m = ngsolve.Mesh(str(p))",
                   "out = solve_fn(m, curve_order=k)"])
    assert [f["line"] for f in found] == [3], found


def test_netgen_load_into_curved_call_is_flagged():
    found = _hits(["ngm = Mesh(dim=3)", 'ngm.Load("yoke.vol")', "run(ngm, curve_order=2)"])
    assert [f["line"] for f in found] == [3], found


def test_load_vol_and_generated_meshes_are_not_flagged():
    assert _hits(["mesh = load_vol('iron.vol')", "vim.Solve(mesh, order=2, curve_order=2)"]) == []
    assert _hits(["mesh = ng.Mesh(geo.GenerateMesh(maxh=0.2))", "vim.Solve(mesh, curve_order=2)"]) == []


def test_straight_requests_and_comments_are_not_flagged():
    assert _hits(['mesh = ng.Mesh("iron.vol")', "vim.Solve(mesh, curve_order=None)",
                  "vim.Solve(mesh, curve_order=1)", "# vim.Solve(mesh, curve_order=2)"]) == []


def test_rebinding_or_detaching_clears_the_mark():
    assert _hits(['mesh = ng.Mesh("a.vol")', "mesh = load_vol('a.vol')",
                  "vim.Solve(mesh, curve_order=2)"]) == []
    assert _hits(['mesh = ng.Mesh("a.vol")', "mesh.ngmesh.SetGeometry(None)",
                  "vim.Solve(mesh, curve_order=2)"]) == []
    assert _hits(["vim.Solve(mesh, curve_order=2)", 'mesh = ng.Mesh("a.vol")']) == []


def test_rule_is_registered():
    assert check_vol_mesh_into_curved_route in ALL_RULES
