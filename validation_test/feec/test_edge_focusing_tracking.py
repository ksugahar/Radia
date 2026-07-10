"""Golden: vertical EDGE FOCUSING of a tilted dipole end, measured by particle tracking.

Locks docs/clebsch_hodograph/edge_focusing_tracking.{ipynb,py}: the particle-tracking
Hill integral recovers the hard-edge vertical edge-focusing law |1/f_z| = tan(beta)/rho
on a genuinely Maxwellian (curl-free) tilted fringe -- the CORRECT measurement of edge
focusing (the field-EFB slope cannot; see memory/edge_focusing_efb_slope_negative).

Pure-numpy + deterministic (no ngsolve): CI-friendly.  The claims locked here are:
  (1) the tracked 1/f_z tracks +tan(beta)/rho over a beta sweep (this fringe orientation
      focuses; magnitude tan/rho is the invariant);
  (2) it CONVERGES to the hard-edge law as the fringe width w -> 0 (slope c(w) -> 1);
  (3) the beta=0 baseline is a finite-fringe residual = -0.5 w/rho -> 0 as w -> 0;
  (4) 1/f_z scales as 1/rho (thin lens): 1/f_z * rho collapses onto tan(beta).
"""
import math
import sys
from pathlib import Path

DOCDIR = Path(__file__).resolve().parents[2] / "docs" / "clebsch_hodograph"
sys.path.insert(0, str(DOCDIR))


def test_hard_edge_law_magnitude():
    import edge_focusing_tracking as ef
    # magnitude tan(beta)/rho, thin-lens 1/rho scaling
    assert abs(ef.hard_edge_law(20.0, 1.0) - math.tan(math.radians(20.0))) < 1e-12
    assert abs(ef.hard_edge_law(20.0, 2.0) - 0.5 * ef.hard_edge_law(20.0, 1.0)) < 1e-12


def test_tracker_recovers_edge_focusing():
    import edge_focusing_tracking as ef
    sw = ef.sweep_beta()                         # rho=1, w=0.02
    # sign of this (curl-free entrance-edge) orientation is +, magnitude ~ tan/rho
    for r in sw:
        if r["beta_deg"] >= 5.0:
            assert r["inv_fz"] > 0.0, r
            # finite fringe under-shoots the delta-limit by < 12 %
            assert abs(r["inv_fz"] - r["hard_edge"]) / abs(r["hard_edge"]) < 0.12, r


def test_tracker_matches_enge_fringe_corrected_law():
    import edge_focusing_tracking as ef
    # The FULL classical law tan(beta - psi)/rho, psi = (K1g/rho)(1+sin^2)/cos with
    # K1g = w/2 for the tanh fringe, matches the tracked values to < 1.5 % INCLUDING
    # the finite fringe (measured 0.03-0.71 % at w=0.02; residual ~ the 2nd-order K2 term).
    for w in (0.02, 0.04):
        sw = ef.sweep_beta(w=w)
        for r in sw:
            if abs(r["enge"]) > 1e-9:
                assert abs(r["inv_fz"] - r["enge"]) / abs(r["enge"]) < 0.02, (w, r)
    # beta=0: the tracked baseline IS the Enge correction -K1g/rho^2 exactly
    b0 = ef.edge_focus_integral(ef.edge_field(0.0, w=0.02), 1.0)["inv_fz"]
    assert abs(b0 - ef.scoff_law(0.0, 1.0, 0.01)) < 2e-4, b0


def test_converges_to_hard_edge_as_w_shrinks():
    import edge_focusing_tracking as ef
    wc = ef.w_convergence()                      # w: 0.08 -> 0.005
    slopes = [d["slope_vs_law"] for d in wc]
    assert slopes[0] < slopes[-1]                # monotone approach
    assert slopes[-1] > 0.98                     # ~0.99 at the finest w
    assert slopes[-1] < 1.02


def test_beta0_baseline_is_finite_fringe_residual():
    import edge_focusing_tracking as ef
    # baseline 1/f_z(beta=0) = -0.5 w/rho, an O(w/rho) residual that vanishes as w -> 0
    for w in (0.04, 0.02, 0.01):
        b = ef.edge_focus_integral(ef.edge_field(0.0, w=w), 1.0)["inv_fz"]
        assert abs(b - (-0.5 * w)) < 0.05 * w + 1e-4, (w, b)


def test_rho_collapse():
    import edge_focusing_tracking as ef
    rc = ef.rho_collapse()
    # 1/f_z * rho lands on tan(beta), nearly rho-independent
    for d in rc:
        for tb, y in zip(d["tan_beta"], d["inv_fz_rho"]):
            if tb > 0.05:
                assert abs(y - tb) < 0.05 + 0.03 * tb, (d["rho"], tb, y)


def test_fem_coil_pair_is_clean():
    """Locks the PART B coil construction against two verified pitfalls (2026-07-10):
    CoilBuilder.mirror('xy') emits mirrored straights running the wrong way (spurious
    odd-in-x dBz/dx up to 0.39 T/m), and rad.TrfOrnt-wrapped containers crash
    rad.RadiaField.  The explicit rounded-parallelogram pair must close to machine
    precision and its mid-plane Bz must be x-even (beta=0) / C2-even (beta=20) to
    ~1e-7 of B0."""
    import pytest
    rad = pytest.importorskip("radia")
    import math
    import numpy as np
    import edge_focusing_tracking as ef

    for bdeg in (0.0, 20.0):
        up = ef._fem_coil_path(+1, math.radians(bdeg))
        gap = np.linalg.norm(np.asarray(up._position)
                             - np.asarray(up.segments[0].start_pos))
        assert gap < 1e-12, (bdeg, gap)
        cnt = ef.fem_build_coil(math.radians(bdeg))
        b0 = rad.Fld(cnt, "b", [0.0, 0.0, 0.0])[2]
        assert 0.05 < b0 < 0.2, b0            # both loops contribute (mirror bug gave 0.063->0.103)
        # C2-odd part must vanish for BOTH beta (the loop pair is C2-symmetric);
        # the mirror('xy') bug sat at 7.7e-3 T here -- 1e-6*B0 is 4 orders below it
        for (x, y) in ((0.02, -0.14), (0.02, -0.10), (0.03, 0.0)):
            c2 = 0.5 * (rad.Fld(cnt, "b", [x, y, 0.0])[2]
                        - rad.Fld(cnt, "b", [-x, -y, 0.0])[2])
            assert abs(c2) < 1e-6 * b0, (bdeg, x, y, c2)
        if bdeg == 0.0:
            # x-odd part must vanish at beta=0 (rounded rectangle is x-symmetric)
            for y in (-0.16, -0.12, -0.10):
                odd = 0.5 * (rad.Fld(cnt, "b", [0.02, y, 0.0])[2]
                             - rad.Fld(cnt, "b", [-0.02, y, 0.0])[2])
                assert abs(odd) < 1e-6 * b0, (y, odd)
        rad.UtiDelAll()
