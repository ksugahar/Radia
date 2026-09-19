"""Automatic fin detection on section outlines (no CAD, pure geometry)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from radia.fin_section import (
    analyze_section,
    feature_panel_weights,
    graded_lane_arclengths,
    local_thickness,
    outline_points_at,
    resample_outline,
)


def _beak_outline(n_arc=64):
    """Same outline as validation_test/.../make_beak_fin_step.py, metres."""
    tl = (3.8, -(0.25**2 - 0.05**2) ** 0.5)
    tu = (3.8, +(0.25**2 - 0.05**2) ** 0.5)
    centre, radius = (3.75, 0.0), 0.25
    a0 = math.atan2(tl[1] - centre[1], tl[0] - centre[0])
    a1 = math.atan2(tu[1] - centre[1], tu[0] - centre[0])
    pts = [(-4, -2), (1.6, -2), (1.6, -0.7), tl]
    pts += [(centre[0] + radius * math.cos(t), centre[1] + radius * math.sin(t))
            for t in np.linspace(a0, a1, n_arc)[1:-1]]
    pts += [tu, (1.6, 0.7), (1.6, 2), (-4, 2)]
    return np.array(pts) * 1e-3


def _rect(w, h):
    return np.array([[-w, -h], [w, -h], [w, h], [-w, h]]) * 0.5


def test_local_thickness_of_rectangle():
    p = resample_outline(_rect(6e-3, 2e-3), 400)
    t = local_thickness(p)
    # skip the corner samples whose averaged tangent is diagonal
    interior = (np.abs(np.abs(p[:, 0]) - 3e-3) > 0.05e-3) & \
               (np.abs(np.abs(p[:, 1]) - 1e-3) > 0.05e-3)
    long_faces = np.abs(p[:, 1]) > 0.99e-3
    assert np.allclose(t[long_faces & interior], 2e-3, rtol=1e-6)
    assert np.allclose(t[~long_faces & interior], 6e-3, rtol=1e-6)


@pytest.mark.parametrize("outline", [
    _rect(6e-3, 2e-3), _rect(4e-3, 4e-3),
    3e-3 * np.c_[np.cos(np.linspace(0, 2 * np.pi, 180, endpoint=False)),
                 np.sin(np.linspace(0, 2 * np.pi, 180, endpoint=False))],
])
def test_plain_sections_have_no_fin(outline):
    analysis = analyze_section(resample_outline(outline, 512))
    assert analysis.fins == ()
    assert analysis.primary is None
    # graded lanes degenerate to equal arc length
    s = graded_lane_arclengths(analysis, 64)
    assert np.allclose(np.diff(s), analysis.perimeter / 64)


def test_beak_fixture_geometry_is_recovered():
    analysis = analyze_section(resample_outline(_beak_outline(), 2048))
    assert len(analysis.fins) == 1
    fin = analysis.primary
    assert analysis.body_thickness == pytest.approx(4e-3, rel=1e-3)
    # root chord at the concave corners (1.6, +-0.7) mm
    assert fin.root_mid[0] == pytest.approx(1.6e-3, abs=0.03e-3)
    assert fin.root_thickness == pytest.approx(1.4e-3, abs=0.05e-3)
    assert fin.extremity[0] == pytest.approx(4.0e-3, abs=0.005e-3)
    assert abs(fin.extremity[1]) < 0.01e-3
    assert fin.axis[0] == pytest.approx(1.0, abs=2e-3)
    assert fin.tip_radius == pytest.approx(0.25e-3, rel=0.02)
    # tip arc from (3.8, -0.245) over the nose to (3.8, +0.245) is 0.685 mm;
    # the faces are tangent to the arc, so the fitted-circle test extends a
    # little onto them (quadratic departure) -- bounded, not exact.
    assert 0.65e-3 <= fin.tip_arc_length <= 0.90e-3
    # the previously hard-coded fixture constants fall out of the geometry
    tip_threshold = (fin.root_mid + fin.axis * (fin.length - fin.tip_span))[0]
    assert tip_threshold == pytest.approx(3.5e-3, abs=0.01e-3)
    probe = fin.probe_points(offset=1.0e-3, n=41)
    assert probe[20][0] == pytest.approx(5.0e-3, abs=0.01e-3)
    assert probe[0][1] == pytest.approx(-2.0e-3, abs=0.02e-3)
    assert probe[-1][1] == pytest.approx(2.0e-3, abs=0.02e-3)
    assert fin.probe_points(offset=1.0e-3, n=1).shape == (1, 2)


def test_masks_reproduce_fixture_thresholds_on_lanes():
    analysis = analyze_section(resample_outline(_beak_outline(), 2048))
    fin = analysis.primary
    lanes = resample_outline(analysis.outline, 256)
    assert np.array_equal(fin.fin_mask(lanes), lanes[:, 0] >= 1.6e-3 - 1e-5)
    assert np.array_equal(fin.tip_mask(lanes), lanes[:, 0] >= 3.5e-3 - 1e-5)
    on_arc = fin.on_tip_arc(analysis.arclength_of(lanes), analysis.perimeter)
    assert 7 <= on_arc.sum() <= 11                    # ~0.69-0.9 mm of 22.97 mm


def test_detection_is_rigid_motion_invariant():
    outline = resample_outline(_beak_outline(), 2048)
    base = analyze_section(outline).primary
    theta = 0.65
    rot = np.array([[math.cos(theta), -math.sin(theta)],
                    [math.sin(theta), math.cos(theta)]])
    moved = analyze_section(outline @ rot.T + np.array([0.01, -0.02])).primary
    assert moved.length == pytest.approx(base.length, rel=1e-6)
    assert moved.tip_radius == pytest.approx(base.tip_radius, rel=1e-3)
    assert moved.root_thickness == pytest.approx(base.root_thickness, rel=1e-3)
    angle = math.atan2(moved.axis[1], moved.axis[0])
    assert angle == pytest.approx(theta + math.atan2(base.axis[1], base.axis[0]),
                                  abs=2e-3)


def test_detection_is_stable_under_outline_density():
    tips, roots = [], []
    for n in (256, 512, 1024, 2048):
        fin = analyze_section(resample_outline(_beak_outline(), n)).primary
        tips.append(fin.tip_radius)
        roots.append(fin.root_mid[0])
    assert np.allclose(tips, 0.25e-3, rtol=0.02)
    assert np.all(np.abs(np.asarray(roots) - 1.6e-3) < 0.08e-3)


def test_two_fins_are_separated():
    body = _rect(8e-3, 4e-3)
    # add a second, shorter fin on the left side by mirroring the beak
    beak = _beak_outline()
    left = beak.copy()
    left[:, 0] = -left[:, 0]
    # merge: right half of beak, left half of mirrored beak
    right_part = beak[beak[:, 0] > 1.5e-3]
    left_part = left[left[:, 0] < -1.5e-3]
    outline = np.vstack([
        [(-1.6e-3, -2e-3), (1.6e-3, -2e-3)],
        right_part[np.argsort(np.arctan2(right_part[:, 1], right_part[:, 0] - 3e-3))],
        [(1.6e-3, 2e-3), (-1.6e-3, 2e-3)],
        left_part[np.argsort(-np.arctan2(left_part[:, 1], -(left_part[:, 0] + 3e-3)))],
    ])
    analysis = analyze_section(resample_outline(outline, 2048))
    assert len(analysis.fins) == 2
    axes = sorted(f.axis[0] for f in analysis.fins)
    assert axes[0] == pytest.approx(-1.0, abs=5e-3)
    assert axes[1] == pytest.approx(1.0, abs=5e-3)
    del body


def test_graded_lanes_put_more_cells_on_the_tip():
    analysis = analyze_section(resample_outline(_beak_outline(), 2048))
    fin = analysis.primary
    for n_lanes in (64, 128, 256):
        s = graded_lane_arclengths(analysis, n_lanes, tip_lanes=8)
        assert len(s) == n_lanes
        assert np.all(np.diff(s) > 0)
        lanes = outline_points_at(analysis.outline, s)
        on_arc = fin.on_tip_arc(analysis.arclength_of(lanes), analysis.perimeter)
        assert on_arc.sum() >= 7                      # ~tip_lanes cells on the arc
        uniform = resample_outline(analysis.outline, n_lanes)
        # at 256 lanes the uniform spacing is already finer than the tip target
        assert fin.tip_mask(lanes).sum() >= fin.tip_mask(uniform).sum()
        # spacing shrinks monotonically toward the tip along each face
        spacing = np.diff(np.r_[s, s[0] + analysis.perimeter])
        face = fin.fin_mask(lanes) & ~fin.tip_mask(lanes)
        proj = fin.project(lanes)
        upper = face & (lanes @ fin.normal > 0)
        order = np.argsort(proj[upper])
        assert np.all(np.diff(spacing[upper][order]) <= 1e-9)


def test_graded_lanes_reject_impossible_budget():
    analysis = analyze_section(resample_outline(_beak_outline(), 1024))
    with pytest.raises(ValueError):
        graded_lane_arclengths(analysis, 8, tip_lanes=64)


def test_feature_weights_conserve_arc_and_resist_half_cell_phase_shift():
    analysis = analyze_section(resample_outline(_beak_outline(), 2048))
    fin = analysis.primary
    n_lanes = 128
    ds = analysis.perimeter / n_lanes
    fractions = []
    dense_region_length = (np.count_nonzero(fin.fin_mask(analysis.outline))
                           * analysis.perimeter / len(analysis.outline))
    for phase in (0.0, 0.5):
        s = np.mod((np.arange(n_lanes) + 0.5 + phase) * ds,
                   analysis.perimeter)
        lanes = outline_points_at(analysis.outline, s)
        weight = feature_panel_weights(analysis, fin, lanes)
        assert np.sum(weight) * ds == pytest.approx(
            dense_region_length, rel=0, abs=analysis.perimeter / 2048)
        k = 1.0 + 0.2 * fin.project(lanes) / fin.length
        current = k * ds
        current /= np.sum(current)
        fractions.append(float(np.sum(current * weight)))
    assert abs(fractions[1] / fractions[0] - 1.0) < 0.005


def test_feature_weights_reject_duplicate_lanes():
    analysis = analyze_section(resample_outline(_beak_outline(), 512))
    lanes = resample_outline(analysis.outline, 64)
    lanes[10] = lanes[9]
    with pytest.raises(ValueError, match="ordered and distinct"):
        feature_panel_weights(analysis, analysis.primary, lanes)
