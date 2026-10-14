"""Golden: vector-T convex volume inverse design + ACA+TSVD (Stage A / #2).

Locks the verified vector-T claims (sympy 2026-06-11) in a real NGSolve solve:
  * CONVEX inverse: J = curl T linear in T -> ACA+TSVD min-norm fit ~ machine eps;
  * div J = 0 (J = curl T, exactly divergence-free);
  * GAUGE null space: a pure-gradient T = grad(chi) produces ~ZERO field (it is
    in ker(A)), and TSVD truncates the gauge subspace automatically (k_aca <= M
    << ndof) -- no explicit tree-cotree gauge.

See validation_test/feec/vim_legacy/vector_t_inverse.py and streamfunction("clebsch_3d").
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                "../../validation_test/feec/vim_legacy"))


@pytest.fixture(scope="module")
def result():
    import vector_t_inverse as vt
    res, ok = vt.main()
    return res, ok


def test_convex_fit_machine_precision(result):
    """Linear least-norm A T = B_target solves to ~ machine eps."""
    res, _ = result
    assert res["fit_res"] < 1e-8, f"ACA+TSVD fit residual {res['fit_res']:.2e}"


def test_independent_recompute(result):
    res, _ = result
    assert res["check_res"] < 5e-2, f"direct Bz(curl T) residual {res['check_res']:.2e}"


def test_curl_T_divergence_free(result):
    """J = curl T is divergence-free to machine precision."""
    res, _ = result
    assert res["div_rel"] < 1e-8, f"div J / |J| {res['div_rel']:.2e}"


def test_gauge_null_space(result):
    """A pure-gradient T = grad(chi) yields ~zero field (gauge in ker A);
    TSVD truncates the gauge subspace so k_aca <= M << ndof."""
    res, _ = result
    assert res["gauge_field"] < 1e-2, f"gauge field {res['gauge_field']:.2e}"
    assert res["k_aca"] <= res["M"], f"k_aca {res['k_aca']} > M {res['M']}"
    assert res["k_aca"] < res["ndof"], "gauge null space not truncated"


def test_overall_pass(result):
    _, ok = result
    assert ok
