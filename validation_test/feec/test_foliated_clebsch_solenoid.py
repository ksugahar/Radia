"""Golden: foliated current-Clebsch volume stream function + ACA+TSVD (Stage A).

Locks that a 3D VOLUME current J = grad(lambda) x grad(mu) with the foliation
mu = r FIXED (so J is LINEAR in lambda) is inverse-designed by the EXISTING
radia.stream_function (ACA+)+TSVD engine to reproduce a uniform-Bz solenoid:
  * ACA+TSVD fit residual small (lambda reproduces target Bz),
  * independent Biot-Savart recompute from the solved lambda matches,
  * the current is divergence-free (weak div, current conservation).

See examples/vim/foliated_clebsch_solenoid.py and streamfunction("clebsch_3d").
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                "../../examples/vim"))


@pytest.fixture(scope="module")
def result():
    import foliated_clebsch_solenoid as fc
    res, ok = fc.main()
    return res, ok


def test_acatsvd_fit_reproduces_target(result):
    res, _ = result
    assert res["fit_res"] < 1e-3, f"ACA+TSVD fit residual {res['fit_res']:.2e}"


def test_independent_biot_savart_recompute(result):
    """Bz recomputed directly from the solved lambda (independent of the
    assembled response matrix) matches the target."""
    res, _ = result
    assert res["check_res"] < 5e-2, f"direct Bz(lambda) residual {res['check_res']:.2e}"


def test_current_is_divergence_free(result):
    """J = grad(lambda) x grad(mu) is div-free (weak divergence ~ 0)."""
    res, _ = result
    assert res["div_rel"] < 1e-3, f"weak div J / |J| {res['div_rel']:.2e}"


def test_overall_pass(result):
    _, ok = result
    assert ok
