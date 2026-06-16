"""
Racetrack coil for the C-type dipole reproduction.

Matches the C-type electromagnet racetrack-coil geometry (MMB8T model
coil dimensions).

Coil geometry (mm):
    Y-straights centerline x: ±47.5 (radial offset from yoke center)
    Y-straights y range:      [100.0, 162.5]   (y mid = 131.25)
    X-caps y top:             185.0 (between 167.5 and 202.5)
    X-caps y bottom:           77.5 (between 60.0 and 95.0)
    X-cap length:              50  (x in [-25, +25])
    Centerline arc radius:    22.5 (= (5 + 40) / 2)
    Cross-section width  (radial): 35
    Cross-section height (axial):  105
    Current:                  2000 AT
    Plane:                    z = 0  (planar racetrack)

This differs from `em_sample_coil.py` ONLY in Y_CENTER (131.25 here
vs 151.25 there) and the consequent y positions of all segments.
The em_sample coil deliberately sits 20 mm farther from the gap
than this coil -- which is why em_sample's B at origin is much smaller
than the reference -228.1 mT for the same yoke + sym BCs.

Usage:
    Set this file as the "Coil script" on the EM panel for any
    reduction sample that targets the C-type electromagnet reference.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                '..', '..', '..', 'src', 'radia'))

from coil_builder import CoilBuilder


def build_coil():
    """Build the ELF racetrack coil with Y_CENTER = 131.25 mm.

    Returns:
        CoilBuilder instance (full coil, NOT mirrored -- the FEM/MSC
        Biot-Savart integration uses the full source regardless of
        symmetry reduction).
    """
    mm = 1e-3

    NI = 2000                         # Ampere-turns
    COIL_RADIAL = 35 * mm             # cross-section radial width
    COIL_AXIAL = 105 * mm             # cross-section axial height (z-dir)
    STRAIGHT_LEN = 62.5 * mm          # y-straight length (162.5 - 100)
    ARC_R_INNER = 5 * mm              # inner arc radius (ELF .meg)
    ARC_R_OUTER = 40 * mm             # outer arc radius (ELF .meg)
    ARC_R_MEAN = (ARC_R_INNER + ARC_R_OUTER) / 2  # 22.5 mm

    # ELF Y_CENTER = 131.25 mm  (= (100 + 162.5) / 2)
    Y_CENTER = 131.25 * mm
    Y_INNER = Y_CENTER - ARC_R_MEAN   # 108.75 mm

    X_CAP = 50 * mm                   # X-cap length (-25 to +25)
    X_CENTER = 47.5 * mm              # inner Y-straight x position

    coil = (CoilBuilder(current=NI)
        .set_start([X_CENTER, Y_INNER, 0])
        .set_cross_section(width=COIL_RADIAL, height=COIL_AXIAL)
        .add_straight(STRAIGHT_LEN)            # +y up
        .add_arc(radius=ARC_R_MEAN, arc_angle=90)
        .add_straight(X_CAP)                   # x-cap top (-x)
        .add_arc(radius=ARC_R_MEAN, arc_angle=90)
        .add_straight(STRAIGHT_LEN)            # -y down
        .add_arc(radius=ARC_R_MEAN, arc_angle=90)
        .add_straight(X_CAP)                   # x-cap bottom (+x)
        .add_arc(radius=ARC_R_MEAN, arc_angle=90))

    return coil


if __name__ == "__main__":
    coil = build_coil()
    segments, current = coil.to_wire_segments()
    print(f"ELF coil: {len(segments)} wire segments, NI = {current} A")
    print(f"Y_CENTER = 131.25 mm  (cf. em_sample_coil.py: 151.25 mm)")
    print(f"y-extent of straights: [108.75, 153.75] mm  (centerline)")
    print(f"closed: {coil.is_closed},  gap: {coil.gap:.3e} m")
