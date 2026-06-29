"""Sprint 1 smoke test for coil_from_step.extract_centerline.

Round-trip: build a known coil with CoilBuilder -> write STEP ->
extract centerline + profiles back -> verify geometry matches.

Test shapes:
  - circular torus: R=30mm, a=3mm circular cross-section
  - racetrack: 2 straights + 2 semicircles, square cross-section

Pass criteria (Sprint 1):
  - centerline arc length matches analytical to 1%
  - mean cross-section area matches analytical to 1%
  - centerline radius (for torus) matches to 1%
  - profile classification picks the right Profile subclass
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_RADIA = os.path.abspath(os.path.join(HERE, '..', '..', '..',
                                         'src', 'radia'))
sys.path.insert(0, SRC_RADIA)

from coil_from_step import extract_centerline, load_step_solid    # noqa: E402
from coil_profile import RectProfile, CircleProfile                # noqa: E402


def _build_torus_step(path, R=0.030, a=0.003):
    from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Axis)
    wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circ = wp.Circle(a).Face()
    torus = circ.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360)
    torus.WriteStep(path)


def _build_square_torus_step(path, R=0.030, w=0.006, h=0.006):
    from netgen.occ import (WorkPlane, Axes, Pnt, Dir, Axis, MoveTo)
    wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    rect = wp.MoveTo(-w / 2, -h / 2).Rectangle(w, h).Face()
    torus = rect.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360)
    torus.WriteStep(path)


def test_circular_torus():
    print("=" * 70)
    print("  Test 1: circular-section torus  R=30mm, a=3mm")
    print("=" * 70)
    R, a = 0.030, 0.003
    step = os.path.join(HERE, "out_smoke_torus_circular.step")
    _build_torus_step(step, R, a)

    # Walk needs a hint (closed loop). Start at (R,0,0) with tangent +y.
    res = extract_centerline(
        step,
        start_hint=(np.array([R, 0, 0]), np.array([0, 1, 0])),
        step_size=0.003,
        verbose=False,
    )

    L = res.arclen[-1] + np.linalg.norm(res.polyline[-1] - res.polyline[0]) \
        if res.closed else res.arclen[-1]
    L_expected = 2 * math.pi * R
    radii = np.linalg.norm(res.polyline[:, :2], axis=1)
    R_mean = float(np.mean(radii))
    z_max = float(np.max(np.abs(res.polyline[:, 2])))
    A_expected = math.pi * a ** 2

    # Profile check: count how many came back as CircleProfile
    n_circle = sum(isinstance(p, CircleProfile) for p in res.profiles)
    n_rect = sum(isinstance(p, RectProfile) for p in res.profiles)

    print(f"  stations  : {len(res.polyline)}  (closed={res.closed})")
    print(f"  arclength : {L:.4f} m  (expected {L_expected:.4f} m, "
          f"err {(L - L_expected) / L_expected * 100:+.2f}%)")
    print(f"  radius    : {R_mean:.4f} m  (expected {R:.4f} m, "
          f"err {(R_mean - R) / R * 100:+.2f}%)")
    print(f"  |z|_max   : {z_max:.2e} m  (expected ~0)")
    print(f"  profile types : Circle={n_circle}, Rect={n_rect}, "
          f"other={len(res.profiles) - n_circle - n_rect}")
    if res.profiles and isinstance(res.profiles[0], CircleProfile):
        r_back = res.profiles[0].r
        print(f"  profile r : {r_back:.4f} m  (expected {a:.4f} m, "
              f"err {(r_back - a) / a * 100:+.2f}%)")


def test_square_torus():
    print()
    print("=" * 70)
    print("  Test 2: square-section torus  R=30mm, 6x6mm")
    print("=" * 70)
    R, w, h = 0.030, 0.006, 0.006
    step = os.path.join(HERE, "out_smoke_torus_square.step")
    _build_square_torus_step(step, R, w, h)

    res = extract_centerline(
        step,
        start_hint=(np.array([R, 0, 0]), np.array([0, 1, 0])),
        step_size=0.003,
        verbose=False,
    )

    L = res.arclen[-1]
    L_expected = 2 * math.pi * R
    radii = np.linalg.norm(res.polyline[:, :2], axis=1)
    R_mean = float(np.mean(radii))
    A_expected = w * h
    n_rect = sum(isinstance(p, RectProfile) for p in res.profiles)
    n_circle = sum(isinstance(p, CircleProfile) for p in res.profiles)

    print(f"  stations  : {len(res.polyline)}  (closed={res.closed})")
    print(f"  arclength : {L:.4f} m  (expected {L_expected:.4f} m, "
          f"err {(L - L_expected) / L_expected * 100:+.2f}%)")
    print(f"  radius    : {R_mean:.4f} m  (expected {R:.4f} m, "
          f"err {(R_mean - R) / R * 100:+.2f}%)")
    print(f"  profile types : Rect={n_rect}, Circle={n_circle}, "
          f"other={len(res.profiles) - n_rect - n_circle}")
    if res.profiles and isinstance(res.profiles[0], RectProfile):
        rect = res.profiles[0]
        print(f"  profile w x h : {rect.w * 1e3:.2f} x {rect.h * 1e3:.2f} mm "
              f"(expected {w * 1e3:.2f} x {h * 1e3:.2f} mm)")


def _build_racetrack_step(path, length=0.040, radius=0.020, w=0.006, h=0.006):
    """Build a planar racetrack (2 straights + 2 semicircles) via CoilBuilder."""
    sys.path.insert(0, SRC_RADIA)
    from coil_builder import CoilBuilder
    cb = (CoilBuilder(current=1.0)
          .set_start([0, -radius, 0],
                     orientation=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
          .set_cross_section(width=w, height=h)
          .add_straight(length)
          .add_arc(radius=radius, arc_angle=180, tilt=0)
          .add_straight(length)
          .add_arc(radius=radius, arc_angle=180, tilt=0))
    cb.write_step(path)


def test_racetrack():
    print()
    print("=" * 70)
    print("  Test 3: racetrack (2 straights + 2 semicircles)")
    print("=" * 70)
    from coil_from_step import (extract_centerline, polyline_to_segments)
    length, radius = 0.040, 0.020
    w, h = 0.006, 0.006
    step = os.path.join(HERE, "out_smoke_racetrack.step")
    try:
        _build_racetrack_step(step, length, radius, w, h)
    except Exception as e:
        print(f"  Build failed: {e}")
        return
    res = extract_centerline(
        step,
        start_hint=(np.array([0, -radius, 0]), np.array([0, 1, 0])),
        step_size=0.003,
    )
    L_expected = 2 * length + 2 * math.pi * radius
    L = res.arclen[-1] + (
        np.linalg.norm(res.polyline[-1] - res.polyline[0])
        if res.closed else 0.0)
    print(f"  stations  : {len(res.polyline)} (closed={res.closed})")
    print(f"  arclength : {L:.4f} m  (expected {L_expected:.4f} m, "
          f"err {(L - L_expected) / L_expected * 100:+.2f}%)")

    segs = polyline_to_segments(res)
    print(f"  Reconstructed {len(segs)} segments (expected 4):")
    for i, s in enumerate(segs):
        if s.kind == 'arc':
            print(f"    [{i}] arc   length={s.length * 1e3:.2f}mm "
                  f"radius={s.radius * 1e3:.2f}mm "
                  f"angle={s.angle_deg:.2f}deg")
        else:
            print(f"    [{i}] straight length={s.length * 1e3:.2f}mm")

    print()
    print("  Building CoilBuilder...")
    from coil_from_step import to_coil_builder
    try:
        builder, _ = to_coil_builder(res, current=1.0)
        print(f"  CoilBuilder OK with {len(builder.segments)} segments")
    except Exception as e:
        print(f"  ERROR: {e}")


def test_sprint2_segments():
    print()
    print("=" * 70)
    print("  Sprint 2: polyline -> CoilBuilder reconstruction")
    print("=" * 70)
    from coil_from_step import (extract_centerline, polyline_to_segments,
                                 to_coil_builder)
    R, a = 0.030, 0.003
    step = os.path.join(HERE, "out_smoke_torus_circular.step")
    res = extract_centerline(
        step,
        start_hint=(np.array([R, 0, 0]), np.array([0, 1, 0])),
        step_size=0.003,
    )
    segs = polyline_to_segments(res)
    print(f"  Reconstructed {len(segs)} segments:")
    for i, s in enumerate(segs):
        if s.kind == 'arc':
            print(f"    [{i}] arc   length={s.length * 1e3:.2f}mm "
                  f"radius={s.radius * 1e3:.2f}mm "
                  f"angle={s.angle_deg:.2f} deg "
                  f"profile={type(s.profile).__name__}")
        else:
            print(f"    [{i}] straight length={s.length * 1e3:.2f}mm "
                  f"profile={type(s.profile).__name__}")
    print(f"  Total length: {sum(s.length for s in segs) * 1e3:.2f}mm "
          f"(expected {2 * math.pi * R * 1e3:.2f}mm)")

    print()
    print("  Building CoilBuilder...")
    try:
        builder, _ = to_coil_builder(res, current=1.0)
        print(f"  CoilBuilder OK with {len(builder.segments)} segments")
    except Exception as e:
        print(f"  ERROR: {e}")


if __name__ == "__main__":
    test_circular_torus()
    test_square_torus()
    test_sprint2_segments()
    test_racetrack()
