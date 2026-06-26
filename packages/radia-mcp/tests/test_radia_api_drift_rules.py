# -*- coding: utf-8 -*-
"""Regression tests for the Radia API-drift / example-bitrot lint rules.

These mirror the bug_patterns catalog entries surfaced by the examples->docs
consolidation (2026-06): magnetization-as-Tesla, the rad.Solve 4-tuple, and the
ObjArcCur / MatSatIsoFrm constructor-arity drifts, plus the cp932 stdout rewrap
that breaks notebook execution.
"""
from radia_mcp.radia_ngsolve.rules import (
    check_solve_result_bad_index,
    check_objarccur_missing_axis,
    check_matsatisofrm_multiple_lists,
    check_magnetization_looks_tesla,
    check_cp932_stdout_rewrap,
    ALL_RULES,
)


def _h(rule, code):
    return rule("t.py", code.split("\n"))


# ── rad.Solve 4-tuple [residual,_,_,iterations]: [4] is out of range ──
def test_solve_index4_flagged():
    f = _h(check_solve_result_bad_index, "res = rad.Solve(c,1e-4,100)\nx = res[4]")
    assert len(f) == 1 and f[0]["rule"] == "solve-result-bad-index"


def test_solve_index3_iteration_read_ok():
    # [3] is the legitimate iteration count -> must NOT be flagged
    assert _h(check_solve_result_bad_index,
              "res = rad.Solve(c,1e-4,100)\nn = int(res[3])") == []


def test_solve_index0_residual_ok():
    assert _h(check_solve_result_bad_index,
              "res = rad.Solve(c,1e-4,100)\nr = res[0]") == []


# ── ObjArcCur needs 8 args incl. man_auto + axis ──
def test_objarccur_6arg_flagged():
    f = _h(check_objarccur_missing_axis, "g = rad.ObjArcCur([0,0,0],[a,b],[c,d],h,n,j)")
    assert len(f) == 1 and f[0]["rule"] == "objarccur-missing-axis"


def test_objarccur_8arg_ok():
    assert _h(check_objarccur_missing_axis,
              "g = rad.ObjArcCur([0,0,0],[a,b],[c,d],h,n,'man','z',j)") == []


def test_objarccur_multiline_open_skipped():
    # multi-line call may carry 'man','z' on a later line -> do not flag the open line
    assert _h(check_objarccur_missing_axis, "coil = rad.ObjArcCur([0, 0, 0],") == []


# ── MatSatIsoFrm takes ONE nested list ──
def test_matsatisofrm_three_positional_flagged():
    f = _h(check_matsatisofrm_multiple_lists,
           "m = rd.MatSatIsoFrm([1596,1.1],[133,0.4],[18,0.5])")
    assert len(f) == 1 and f[0]["rule"] == "matsatisofrm-multiple-lists"


def test_matsatisofrm_single_nested_ok():
    assert _h(check_matsatisofrm_multiple_lists,
              "m = rd.MatSatIsoFrm([[1596,1.1],[133,0.4]])") == []


# ── magnetization that looks like Br (Tesla) instead of A/m ──
def test_magnetization_tesla_flagged():
    f = _h(check_magnetization_looks_tesla,
           "g = rad.ObjCylMag([0,0,0],r,h,16,'z',[0,0,1.0])")
    assert len(f) == 1 and f[0]["rule"] == "magnetization-looks-tesla"


def test_magnetization_apm_ok():
    assert _h(check_magnetization_looks_tesla,
              "g = rad.ObjHexahedron(verts,[0,0,954930])") == []


def test_dims_and_center_not_flagged_as_mag():
    # 3-nonzero dims / all-zero center must not look like a single-axis Br
    assert _h(check_magnetization_looks_tesla,
              "g = rad.ObjRecMag([0,0,0],[0.02,0.02,0.02],[0,0,954930])") == []


def test_non_constructor_vector_ok():
    assert _h(check_magnetization_looks_tesla, "B = rad.Fld(c,'b',[0,0,0.1])") == []


# ── cp932 stdout rewrap (breaks Jupyter) ──
def test_cp932_stdout_rewrap_flagged():
    f = _h(check_cp932_stdout_rewrap,
           "sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer,'strict')")
    assert len(f) == 1 and f[0]["rule"] == "cp932-stdout-rewrap"


# ── all 5 registered ──
def test_rules_registered():
    for r in (check_solve_result_bad_index, check_objarccur_missing_axis,
              check_matsatisofrm_multiple_lists, check_magnetization_looks_tesla,
              check_cp932_stdout_rewrap):
        assert r in ALL_RULES
