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
