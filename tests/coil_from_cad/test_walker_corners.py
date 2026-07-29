"""Walker sharp-corner turning + bidirectional march locks.

Background (2026-07-29): ``coil_from_step.extract_centerline`` halted
whenever the Frenet tangent could not bend around a sharp spine corner
(the adaptive step-halving follows curvature, not kinks) and only ever
marched FORWARD from the seed.  Consequences on real developer STEPs:
B_rect_sweep traced a single ~97 mm leg of a ~600 mm racetrack and
J_boss_fused ~30% of its loop -- silently, until the 0614ca18d coverage
check turned them into hard errors.  The fix: a minimum-section-area
direction fan turns corners (same principle as the axis-agnostic seed),
and a backward march from the same seed covers the other side of an
open coil whose seed landed mid-wire.

These tests lock both behaviours on synthetic sharp-corner solids:
  1. square frame (4 x 90-deg corners, closed) -> closed walk on the
     true mid-frame centerline;
  2. U-bar (2 x 90-deg corners, open, seeded mid-bottom by the
     axis-agnostic seed) -> the stitched walk reaches BOTH arm ends.
"""

import numpy as np
import pytest

pytest.importorskip("netgen.occ")


def _square_frame():
    """Closed square 'coil': 100 x 100 mm outer, 14 x 8 mm wire.

    Centerline = mid-frame square at +-43 mm, z = 4 mm.  Authored in
    metres.
    """
    from netgen.occ import Box, Pnt
    outer = Box(Pnt(-0.050, -0.050, 0.0), Pnt(0.050, 0.050, 0.008))
    inner = Box(Pnt(-0.036, -0.036, -1.0), Pnt(0.036, 0.036, 1.0))
    return outer - inner


def _u_bar():
    """Open U: bottom leg + two vertical arms, 14 x 8 mm wire."""
    from netgen.occ import Box, Pnt
    left = Box(Pnt(-0.050, -0.050, 0.0), Pnt(-0.036, 0.050, 0.008))
    right = Box(Pnt(0.036, -0.050, 0.0), Pnt(0.050, 0.050, 0.008))
    bottom = Box(Pnt(-0.050, -0.050, 0.0), Pnt(0.050, -0.036, 0.008))
    return left + right + bottom


def test_square_frame_walk_closes_through_corners():
    """4 sharp 90-deg corners must be turned and the loop closed."""
    from radia.coil_from_step import _axis_agnostic_seed, extract_centerline

    solid = _square_frame()
    seed = _axis_agnostic_seed(solid)
    res = extract_centerline(solid, start_hint=seed)

    assert res.closed, (
        f"square-frame walk did not close ({len(res.polyline)} stations)"
    )
    pts = np.asarray(res.polyline)
    assert pts.shape[0] > 20
    assert np.all(np.isfinite(pts))
    # Stations stay on the mid-frame centerline band: |x|,|y| <= 50 mm
    # (inside the solid) and the extreme stations reach every leg.
    assert float(np.max(np.abs(pts[:, :2]))) < 0.050
    for axis in (0, 1):
        assert float(pts[:, axis].max()) > 0.038, f"axis {axis} + leg missed"
        assert float(pts[:, axis].min()) < -0.038, f"axis {axis} - leg missed"
    # Wire mid-height.
    np.testing.assert_allclose(pts[:, 2], 0.004, atol=1.5e-3)
    # Cross-section area continuity: median = 14 x 8 mm = 112 mm^2.
    assert res.areas is not None
    med = float(np.median(res.areas))
    assert 0.7 * 112e-6 < med < 1.5 * 112e-6


def test_u_bar_bidirectional_covers_both_arms():
    """A mid-bottom seed must cover BOTH arms (backward march) and turn
    both corners; the walk stays open (true ends exist)."""
    from radia.coil_from_step import _axis_agnostic_seed, extract_centerline

    solid = _u_bar()
    seed = _axis_agnostic_seed(solid)
    res = extract_centerline(solid, start_hint=seed)

    assert not res.closed
    pts = np.asarray(res.polyline)
    assert pts.shape[0] > 15
    assert np.all(np.isfinite(pts))
    # Both arms reached (pre-fix: forward-only covered ONE arm).
    assert float(pts[:, 0].max()) > 0.038, "right arm missed"
    assert float(pts[:, 0].min()) < -0.038, "left arm missed"
    # Arms walked upward well past the bottom leg.
    assert float(pts[:, 1].max()) > 0.030, "arm tops not reached"
    assert float(pts[:, 1].min()) < -0.038, "bottom leg missed"
    # Chain is a single open polyline whose two ENDS are the two arm
    # tops (stitch orientation correct: ends far apart, both high-y).
    end_a, end_b = pts[0], pts[-1]
    assert float(np.linalg.norm(end_a - end_b)) > 0.05
    assert end_a[1] > 0.030 and end_b[1] > 0.030
