"""Golden: Frontier F3 wall -- nonzero-helicity current has NO global Clebsch.

Locks the HONEST WALL (an intentionally partly-negative result): a current with
nonzero helicity cannot be expressed as a Clebsch current J = grad(lambda) x
grad(mu), so it has NO clean level-set / streamline WIRES.  The golden asserts
the WALL behaviour, not a solve:

  * the ARNOLD-BELTRAMI-CHILDRESS target is genuinely helical (H_rel ~ 1),
  * the foliation-INDEPENDENT helicity GAP between the target (~1) and ANY
    foliated-Clebsch fit (~0) is ~ 1 -- a conserved topological invariant
    (Moffatt 1969; Arnold & Khesin 1998) the convex hypothesis class cannot
    carry, REGARDLESS of the foliation chosen,
  * a Clebsch (azimuthal) current line CLOSES (returns to its seed) but the
    helical ABC current line does NOT -- it drifts ergodically (Dombre et al.
    1986, J. Fluid Mech. 167:353).

The only honest outputs for a helical target are the continuous vector-T
distribution (J = curl T, Stage A) or a multi-patch Clebsch atlas with
cohomology cuts -- NOT windable equal-current wires.

See examples/vim/helical_current_no_clebsch.py and streamfunction("clebsch_3d").
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                "../../examples/vim"))


@pytest.fixture(scope="module")
def result():
    import helical_current_no_clebsch as f3
    res, ok = f3.main()
    return res, ok


def test_target_is_genuinely_helical(result):
    """The ABC/Beltrami target current must have H_rel ~ 1 (nonzero helicity)."""
    res, _ = result
    assert res["h_rel"] > 0.9, f"target H_rel {res['h_rel']:.3f} not helical"


def test_clebsch_fit_carries_no_helicity(result):
    """ANY foliated-Clebsch fit is FORCED to zero helicity (T.J = 0 identity)."""
    res, _ = result
    assert res["clebsch_fit_hel"] < 0.1, \
        f"Clebsch fit helicity {res['clebsch_fit_hel']:.3f} should be ~0"


def test_helicity_gap_is_the_wall(result):
    """The foliation-independent helicity GAP (target ~1 - Clebsch ~0) ~ 1.

    This is the fundamental topological wall: a conserved invariant that the
    helicity-zero Clebsch hypothesis class can never reproduce, for ANY mu.
    """
    res, _ = result
    assert res["helicity_gap"] > 0.8, \
        f"helicity gap {res['helicity_gap']:.3f} too small -- wall not shown"


def test_clebsch_streamline_closes(result):
    """A Clebsch (azimuthal) current line CLOSES into a ring (arclen = 2*pi*r)."""
    res, _ = result
    assert res["clebsch_closed"] is True, "Clebsch current line should close"
    assert res["clebsch_close_err"] < 5e-3, \
        f"Clebsch ring close error {res['clebsch_close_err']:.1e} too large"


def test_helical_streamline_does_not_close(result):
    """The helical ABC current line does NOT close -- it wanders ergodically.

    No equal-current loop family can be cut from non-closing current lines.
    """
    res, _ = result
    assert res["abc_closed"] is False, \
        "ABC helical current line must NOT close (frontier wall)"
    assert res["abc_drift_rings"] > 5.0, \
        f"ABC line drift {res['abc_drift_rings']:.1f} ring-lengths -- expected wandering"


def test_overall_pass(result):
    _, ok = result
    assert ok
