"""Golden: helicity diagnostic (Stage B / B1).

Locks the dimensionless current-helicity diagnostic
    H_rel = |int T . curl T| / (||T|| ||curl T||)  in [0,1]
that decides whether a clean Clebsch / level-set WIRE EXTRACTION is possible:
  * axial T (T.curl T = 0, Clebsch-type)      -> H_rel ~ 0  (extractable),
  * ABC/Beltrami flow (curl T = k T)          -> H_rel ~ 1  (linked, frontier).

See examples/vim/helicity_diagnostic.py and streamfunction("clebsch_3d").
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                "../../examples/vim"))


@pytest.fixture(scope="module")
def result():
    import helicity_diagnostic as hd
    res, ok = hd.main()
    return res, ok


def test_clebsch_type_is_zero_helicity(result):
    """An axial T (T perp J) is Clebsch-representable -> H_rel ~ 0."""
    res, _ = result
    assert res["h_axial"] < 5e-2, f"axial H_rel {res['h_axial']:.2e}"


def test_beltrami_is_maximal_helicity(result):
    """An ABC/Beltrami flow (curl T = k T) is maximally helical -> H_rel ~ 1."""
    res, _ = result
    assert res["h_abc"] > 0.9, f"ABC H_rel {res['h_abc']:.3f}"


def test_diagnostic_separates_regimes(result):
    """The diagnostic must clearly separate extractable vs frontier."""
    res, _ = result
    assert res["h_abc"] - res["h_axial"] > 0.8


def test_overall_pass(result):
    _, ok = result
    assert ok
