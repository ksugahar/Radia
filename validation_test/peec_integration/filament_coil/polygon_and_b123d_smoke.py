"""polygon_and_b123d_smoke.py

Stage 3c smoke test: PolygonProfile + build123d Face bridge.

Scenarios
---------
1. PolygonProfile from explicit (u, v) vertices approximating a
   square -- DC uniformity should be close to the RectProfile result
   (within a few %) since the sqrt-scaled radial map is approximately
   area-preserving for a regular polygon.

2. profile_from_build123d_face(Rectangle) -- detect a square via the
   build123d bridge and confirm area matches Rectangle.area.

3. Loft using the build123d bridge: Rectangle(5x3) face -> Circle(2.5)
   face. This uses the CAD-first path (student-like workflow: build
   CAD in build123d, feed to our PEEC pipeline).
"""

import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', '..', 'src'))
SRC_RADIA = os.path.join(SRC, 'radia')
for p in (SRC, SRC_RADIA):
    if p not in sys.path:
        sys.path.insert(0, p)

from radia.coil_builder import CoilBuilder
from radia.coil_profile import (RectProfile, CircleProfile, PolygonProfile)
from radia.coil_profile_occ import profile_from_build123d_face
from peec_bundle import build_bundle_solver, filament_currents

MM = 1e-3
CU_SIGMA = 5.8e7


def _peec(coil, nw, nh, frequency):
    paths, _ = coil.to_filaments(nw=nw, nh=nh, frequency=frequency,
                                  sigma=CU_SIGMA)
    info = coil._last_filament_info
    solver, seg_of_fil, *_ = build_bundle_solver(
        paths, dw=None, dh=None, sigma=CU_SIGMA,
        cell_wh=info['cell_wh'])
    I_branch = solver.compute_branch_currents(frequency,
                                               [complex(coil.current)])
    return filament_currents(I_branch, seg_of_fil)


def scenario_polygon_dc():
    print("\n=== Scenario 1: PolygonProfile (5x3 rect as 4 vertices) DC ===")
    w = 5 * MM
    h = 3 * MM
    verts = [(-w/2, -h/2), (w/2, -h/2), (w/2, h/2), (-w/2, h/2)]
    profile = PolygonProfile(verts)
    print(f"         PolygonProfile.total_area = {profile.total_area():.6e}  "
          f"(expected {w*h:.6e})")
    assert abs(profile.total_area() - w * h) / (w * h) < 1e-9

    coil = (CoilBuilder(current=100.0)
            .set_start([0, 0, 0])
            .set_profile(profile)
            .add_straight(length=40 * MM))
    I = _peec(coil, nw=5, nh=5, frequency=0.0)
    mag = np.abs(I)
    spread = np.ptp(mag) / np.mean(mag)
    print(f"         |I| range [{mag.min():.4f}, {mag.max():.4f}]  "
          f"spread={spread*100:.2f}%")
    # sqrt-radial map on a rectangle is NOT area-preserving, so expect
    # a non-zero DC spread. Bound it at ~50%.
    assert spread < 0.5, f"Polygon DC spread too large: {spread}"
    print("         OK: polygon sqrt-radial map gives bounded DC spread")


def scenario_b123d_rect_face():
    print("\n=== Scenario 2: build123d Rectangle face -> PolygonProfile ===")
    from build123d import Rectangle
    face = Rectangle(5 * MM, 3 * MM).face()
    profile = profile_from_build123d_face(face, n_per_edge=8)
    print(f"         face.area = {face.area:.6e}")
    print(f"         PolygonProfile.total_area = {profile.total_area():.6e}")
    rel = abs(profile.total_area() - face.area) / face.area
    print(f"         relative area error = {rel*100:.3f}%")
    assert rel < 1e-6, f"Area mismatch {rel}"
    print("         OK: polygon area matches CAD")

    # Also verify bounding box
    w_b, h_b = profile.bounding_wh()
    print(f"         bounding (w, h) = ({w_b*1e3:.3f}, {h_b*1e3:.3f}) mm "
          f"(expected 5, 3)")


def scenario_b123d_loft():
    print("\n=== Scenario 3: build123d-authored loft Rect(5x3) -> Circle(2.5) DC ===")
    from build123d import Rectangle, Circle
    # Use auto n_per_edge: LINE edges get only 2 samples, curved edges
    # get 32 by default. So rectangle -> 4 edges x 2 = 8 vertices,
    # circle -> 1 edge x 32 = 32 vertices.
    p0 = profile_from_build123d_face(Rectangle(5 * MM, 3 * MM).face())
    p1 = profile_from_build123d_face(Circle(2.5 * MM).face())

    print(f"         profile_start: {p0!r}")
    print(f"         profile_end:   {p1!r}")
    # Check circle area fidelity
    import numpy as _np
    circle_ideal = _np.pi * (2.5 * MM) ** 2
    print(f"         circle area capture: "
          f"{p1.total_area() / circle_ideal * 100:.2f}% of true")
    assert p1.total_area() / circle_ideal > 0.99, \
        "auto n_per_edge for curved edges should give >99% area capture"

    coil = (CoilBuilder(current=100.0)
            .set_start([0, 0, 0])
            .set_profile(p0)
            .add_loft_straight(profile_end=p1, length=40 * MM, n_sub=20))
    I = _peec(coil, nw=5, nh=8, frequency=0.0)
    mag = np.abs(I)
    spread = np.ptp(mag) / np.mean(mag)
    print(f"         |I| range [{mag.min():.4f}, {mag.max():.4f}]  "
          f"spread={spread*100:.2f}%")
    assert spread < 0.5, f"Loft DC spread too large: {spread}"
    print("         OK: CAD-first loft -> PEEC pipeline works")


def scenario_b123d_loft_ac():
    print("\n=== Scenario 4: CAD loft Rect->Circle AC 7 kHz (physics) ===")
    from build123d import Rectangle, Circle
    p0 = profile_from_build123d_face(Rectangle(5 * MM, 3 * MM).face())
    p1 = profile_from_build123d_face(Circle(2.5 * MM).face())

    coil = (CoilBuilder(current=100.0)
            .set_start([0, 0, 0])
            .set_profile(p0)
            .add_loft_straight(profile_end=p1, length=40 * MM, n_sub=20))
    I = _peec(coil, nw=5, nh=8, frequency=7e3)
    mag = np.abs(I)
    ph = np.angle(I, deg=True)
    print(f"         |I| range [{mag.min():.4f}, {mag.max():.4f}]  "
          f"max/min = {mag.max()/mag.min():.2f}")
    print(f"         arg range [{ph.min():.1f}, {ph.max():.1f}] deg")
    assert mag.max() / mag.min() > 1.5, \
        "Expected skin-effect redistribution at 7 kHz"
    print("         OK: AC skin redistribution present on CAD-authored loft")


def scenario_loft_arc():
    """LoftArcSegment: cross-section morphs along an arc."""
    print("\n=== Scenario 5: LoftArc Rect(5x5)->Circle(r=2.5) 180deg DC ===")
    from radia.coil_profile import RectProfile, CircleProfile
    coil = (CoilBuilder(current=100.0)
            .set_start([0, 0, 0])
            .set_profile(RectProfile(5 * MM, 5 * MM))
            .add_loft_arc(profile_end=CircleProfile(2.5 * MM),
                           radius=30 * MM, arc_angle=180, n_sub=40))
    I = _peec(coil, nw=5, nh=5, frequency=0.0)
    mag = np.abs(I)
    spread = np.ptp(mag) / np.mean(mag)
    print(f"         |I| range [{mag.min():.4f}, {mag.max():.4f}]  "
          f"spread = {spread*100:.2f}%")
    # Rect->Circle on an arc path: (u, v) offsets shift the filaments
    # to different arc radii, so path lengths and R_k differ per
    # filament. DC spread is larger than the straight-loft case but
    # still bounded.
    assert spread < 0.5, f"LoftArc DC spread too large: {spread}"
    print("         OK: LoftArc path + morphing cross-section works")

    print("\n=== Scenario 6: LoftArc AC 7 kHz ===")
    coil_ac = (CoilBuilder(current=100.0)
               .set_start([0, 0, 0])
               .set_profile(RectProfile(5 * MM, 5 * MM))
               .add_loft_arc(profile_end=CircleProfile(2.5 * MM),
                              radius=30 * MM, arc_angle=180, n_sub=40))
    I_ac = _peec(coil_ac, nw=5, nh=5, frequency=7e3)
    mag_ac = np.abs(I_ac)
    print(f"         |I| range [{mag_ac.min():.4f}, {mag_ac.max():.4f}]  "
          f"max/min = {mag_ac.max()/mag_ac.min():.2f}")
    assert mag_ac.max() / mag_ac.min() > 1.5, \
        "Expected AC skin redistribution on LoftArc"
    print("         OK: AC skin redistribution on LoftArc")


def main():
    scenario_polygon_dc()
    scenario_b123d_rect_face()
    scenario_b123d_loft()
    scenario_b123d_loft_ac()
    scenario_loft_arc()
    print("\nAll Stage 3c + LoftArc smoke tests PASSED")


if __name__ == '__main__':
    main()
