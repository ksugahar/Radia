"""Golden: 2D bending-magnet END PACK optimization vs iron saturation.

Locks the two-part finding of docs/clebsch_hodograph/demos/bending_endpack_saturation_opt.py
on a small-gap (24 mm), near-knee end pack:

  (1) the MEDIAN-plane field the beam sees is NATURALLY saturation-invariant -- the
      normalized profile p(s)=B_z/B_body is the same curve linear and saturated (J~1e-3,
      L_eff shift < 0.2 %) even for a hard flat cut (the gap-reluctance-dominated pole
      face stays equipotential); and

  (2) the real saturation problem is the pole-TIP CORNER, which concentrates flux and
      saturates deeply (peak iron |B| past the knee).  The optimization minimizes the
      corner kappa = peak iron |B|/body over the end chamfer (depth, exponent) -> the
      chamfer RELIEVES the corner (kappa ~1.9 -> ~1.0) while the beam field stays put.

ngsolve only (Optuna if present, else the built-in grid + local refine).
"""
import sys
from pathlib import Path

import pytest

EXDIR = Path(__file__).resolve().parents[2] / "docs" / "clebsch_hodograph" / "demos"
sys.path.insert(0, str(EXDIR))


@pytest.mark.slow
def test_bending_endpack_saturation_opt():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import bending_endpack_saturation_opt as be
    r = be.optimize(trials=12, seed=0)
    f, o = r["flat_cut"], r["optimized"]

    # (1) the BEAM-plane field is naturally saturation-invariant -- even the flat cut.
    assert f["J_profile_invariance"] < 5e-3, f            # ~1e-3
    assert abs(f["L_eff_shift_rel"]) < 0.02, f            # < 2% (actually < 0.2%)
    assert o["J_profile_invariance"] < 5e-3, o
    assert r["beam_profile_stays_invariant"] is True, r

    # (2) the pole-tip corner is a real HOT SPOT on the flat cut (kappa well > 1),
    # deeply saturated (peak iron |B| past the knee), and the chamfer RELIEVES it.
    assert 1.4 < f["corner_kappa_sat"] < 2.6, f           # flat corner concentration
    assert f["peak_iron_B_sat"] > be.BK, f               # past the iron knee
    assert o["corner_kappa_sat"] < 1.35, o               # relieved toward 1
    assert r["corner_relieved"] is True, r
    assert r["kappa_reduction_factor"] > 1.3, r          # ~1.9x

    # the optimized chamfer is a real, positive, in-range end profile.
    assert 0.001 < o["depth_m"] < be.G2, o               # positive lift, < half-gap
    assert 0.4 <= o["exponent"] <= 3.0, o
    assert r["n_evals"] >= 12, r
