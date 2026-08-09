#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Regression test: CoilBuilder OCC/STEP pose must match the centerline.

Locks the 2026-08-07 fix.  ``CoilSegment.euler_angles`` holds the NEGATED
intrinsic ZXZ angles (a, b, c) of ``orientation``, so the local->world
rotation is ``Rz(ea[2]) Rx(ea[1]) Rz(ea[0])`` -- which is exactly what
``to_radia`` composes (TrfCmbR chains to the right).  The three
``to_occ_shape`` implementations applied the same three angles in the
OPPOSITE order, i.e. ``Rz(ea[0]) Rx(ea[1]) Rz(ea[2])``.  Swapping the two
Z angles is NOT the inverse rotation; it coincides with the correct one
only for special poses (identity, 180-degree), which is why the bug
survived: a planar racetrack looked fine while every TILTED arc came out
mirrored about its own arc centre.

Consequence before the fix, measured on a cylinder-wrapped saddle coil
(R = 40 mm, 2*phi = 70 deg, L = 120 mm, fillet 6 mm): the STEP held 11
solids instead of 8, spanned y = -89.1..64.0 mm instead of -64..64, and
placed end-turn solids at radii 31.5 / 46.4 mm instead of ~40 mm.  The
FIELD was unaffected (``to_radia`` composed correctly), so the defect was
visible only in the exported CAD.

Invariants locked here:
  1. Rz(ea[2]) Rx(ea[1]) Rz(ea[0]) == orientation^T, for random poses;
  2. every segment's OCC solid CONTAINS its own centerline and never
     strays further than the cross-section half-diagonal beyond it, for
     a sweep of arc tilts and for a rotated start orientation;
  3. a saddle coil exports exactly one solid per segment, on-cylinder.
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from radia.coil_builder import CoilBuilder  # noqa: E402

occ = pytest.importorskip("netgen.occ", reason="netgen.occ not installed")

MM = 1e-3
W_CS, H_CS = 10 * MM, 8 * MM
# half-diagonal of the cross-section: the largest distance any point of
# the solid may sit from the centerline
HALF_DIAG = 0.5 * math.hypot(W_CS, H_CS)


def _euler_product(seg):
    from scipy.spatial.transform import Rotation

    ea = seg.euler_angles
    rz = lambda d: Rotation.from_euler("z", d, degrees=True).as_matrix()
    rx = lambda d: Rotation.from_euler("x", d, degrees=True).as_matrix()
    return rz(ea[2]) @ rx(ea[1]) @ rz(ea[0])


def _shape_bbox(shape):
    b = shape.bounding_box
    return np.array([[b[0][i], b[1][i]] for i in range(3)])


def _path_bbox(coil, n_arc=400):
    segs, _ = coil.to_wire_segments(n_arc=n_arc)
    p = np.array([s[0] for s in segs] + [segs[-1][1]])
    return np.array([[p[:, i].min(), p[:, i].max()] for i in range(3)]), p


def _assert_solid_hugs_path(coil, tag):
    """The solid must cover the centerline and not stray beyond the
    cross-section half-diagonal."""
    shape = coil.to_occ()
    ob = _shape_bbox(shape)
    pb, _ = _path_bbox(coil)
    outward = float(np.maximum(pb[:, 0] - ob[:, 0], ob[:, 1] - pb[:, 1]).max())
    missing = float(np.maximum(ob[:, 0] - pb[:, 0], pb[:, 1] - ob[:, 1]).max())
    assert missing <= 1e-9, (
        f"{tag}: solid does not cover the centerline "
        f"(gap {missing * 1e3:.2f} mm); OCC {ob * 1e3}, path {pb * 1e3}")
    assert outward <= HALF_DIAG + 1e-9, (
        f"{tag}: solid strays {outward * 1e3:.2f} mm beyond the centerline, "
        f"more than the {HALF_DIAG * 1e3:.2f} mm cross-section half-diagonal")


def test_euler_product_is_the_inverse_orientation():
    """The three stored angles, composed Z(ea2) X(ea1) Z(ea0), must be the
    local->world rotation.  The pre-fix order Z(ea0) X(ea1) Z(ea2) is a
    different matrix for a general pose."""
    from scipy.spatial.transform import Rotation

    rng = np.random.default_rng(7)
    for k in range(8):
        ori = Rotation.random(random_state=int(rng.integers(1 << 30))
                              ).as_matrix()
        coil = (CoilBuilder(1000.0)
                .set_start([0.0, 0.0, 0.0], ori)
                .set_cross_section(W_CS, H_CS)
                .add_straight(60 * MM))
        seg = coil.segments[-1]
        assert np.allclose(_euler_product(seg), ori.T, atol=1e-9), k


@pytest.mark.parametrize("tilt", [0, 45, 90, 135, 180, -45, -90])
def test_arc_solid_matches_centerline_for_every_tilt(tilt):
    coil = (CoilBuilder(1000.0)
            .set_start([0.0, 0.0, 0.0])
            .set_cross_section(W_CS, H_CS)
            .add_arc(radius=30 * MM, arc_angle=90, tilt=tilt))
    _assert_solid_hugs_path(coil, f"arc tilt={tilt}")


def test_solids_match_centerline_for_a_rotated_start_frame():
    # travel = +z, radial reference = +x
    ori = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]])
    for tilt in (0, 90, -90):
        coil = (CoilBuilder(1000.0)
                .set_start([0.02, -0.01, 0.0], ori)
                .set_cross_section(W_CS, H_CS)
                .add_straight(40 * MM)
                .add_arc(radius=25 * MM, arc_angle=120, tilt=tilt))
        _assert_solid_hugs_path(coil, f"rotated frame, tilt={tilt}")


def _saddle(current, phi_center_deg, *, R=40 * MM, L=120 * MM,
            phi_half=35.0, r_fillet=6 * MM):
    """The tilt-search recipe validated here was promoted to
    radia.coil_builder.saddle_coil (2026-08-09); the pose test now
    exercises the shipped API instead of a private copy."""
    from radia.coil_builder import saddle_coil

    return saddle_coil(current, radius=R, length=L,
                       span_deg=2 * phi_half, bend_radius=r_fillet,
                       width=W_CS, height=H_CS, axis="y",
                       phi_center_deg=phi_center_deg)


def test_saddle_coil_step_is_on_the_cylinder():
    R, L = 40 * MM, 120 * MM
    coil = _saddle(5000.0, 0.0, R=R, L=L)

    # the centerline itself: closes exactly, stays on the cylinder except
    # for the fillets (excursion ~ r_fillet^2 / 2R)
    _, p = _path_bbox(coil, n_arc=64)
    assert np.linalg.norm(p[-1] - p[0]) < 1e-9
    r = np.hypot(p[:, 0], p[:, 2])
    assert r.min() == pytest.approx(R, abs=1e-9)
    assert r.max() - R < 1.0 * MM

    shape = coil.to_occ()
    solids = list(shape.solids)
    assert len(solids) == len(coil.segments) == 8

    half_w = 0.5 * max(W_CS, H_CS)
    for k, s in enumerate(solids):
        b = s.bounding_box
        centre = np.array([(b[0][i] + b[1][i]) / 2 for i in range(3)])
        assert abs(centre[1]) <= L / 2 + half_w, (k, centre)
        r_c = math.hypot(centre[0], centre[2])
        # arc solids report the chord-midpoint radius (R cos(phi_half) .. R)
        assert R * math.cos(math.radians(35.0)) - half_w <= r_c <= R + half_w, (
            k, r_c)

    _assert_solid_hugs_path(coil, "saddle")
