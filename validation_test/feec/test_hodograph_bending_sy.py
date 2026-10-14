"""Golden: 2D bending-magnet END in the s-y (longitudinal) plane.

Locks docs/clebsch_hodograph/hodograph_bending_sy.{ipynb,py}:

  (1) fringe feasibility law d* = (2/pi) h (same Cauchy-Kovalevskaya bound as the x-y
      cross-section, now bounding the LONGITUDINAL fringe / Enge-edge sharpness);
  (2) box-free manufactured-solution FEM verify of the demanded fringe + pole-face equipotential;
  (3) end-shaping design: a near-square end has a strongly enhanced pole-face field, while the
      FEM-optimized chamfer drives the peak toward B0 (the Rogowski no-enhancement limit).
"""
import math
import sys
from pathlib import Path

import pytest

DOCDIR = Path(__file__).resolve().parents[2] / "docs" / "clebsch_hodograph"
sys.path.insert(0, str(DOCDIR))


def test_feasibility_and_fringe_analytic():
    import numpy as np
    import hodograph_bending_sy as hb
    # d* = (2/pi) h
    assert abs(hb.d_star(1.0) - 2.0 / math.pi) < 1e-12
    tab = {round(r["d"], 3): r["feasible"] for r in
           hb.feasibility_table([1.4, 1.0, 0.8, 0.5, 0.35])}
    assert tab[1.4] and tab[1.0] and tab[0.8]
    assert not tab[0.5] and not tab[0.35]
    # closed form: B_y(s,0) = g(s), B_s(s,0) = 0, phi(s,0) = 0
    ss = np.linspace(-3.0, 3.0, 121)
    By, Bs = hb.field_complex(ss, 0.0 * ss)
    assert np.max(np.abs(By - hb.g(ss))) < 1e-12
    assert np.max(np.abs(Bs)) < 1e-12
    assert np.max(np.abs(hb.phi_np(ss, 0.0 * ss))) < 1e-12
    # pole face deep inside sits at the gap h
    yp = hb.pole_profile(np.array([-3.0]))
    assert abs(yp[0] - hb.H) < 5e-3


@pytest.mark.slow
def test_fem_fringe_verify():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import hodograph_bending_sy as hb
    res = hb.fem_verify_sy(order=4)
    assert res["feasible"] is True, res
    assert res["By_rel_err"] < 1e-3, res       # ~1e-4 (0.01 %)
    assert res["phi_L2_err"] < 1e-5, res       # ~3e-9
    assert res["pole_match_err"] < 5e-4, res   # ~1e-6


@pytest.mark.slow
def test_end_shaping_design_recovers_rogowski():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import hodograph_bending_sy as hb
    # near-square end: strong pole-face enhancement; optimized chamfer: peak -> ~B0
    _, sq = hb.solve_endpeak(hb.ellipse_profile(0.1, 0.1, 60), maxh=0.07)
    _, opt = hb.solve_endpeak(hb.ellipse_profile(*hb.DESIGN_OPT, 60), maxh=0.07)
    assert sq > 1.5, sq                        # ~2.1 (square end -> singular)
    assert opt < 1.12, opt                     # ~1.03 (Rogowski, no enhancement)
    assert opt < sq - 0.3, (opt, sq)
