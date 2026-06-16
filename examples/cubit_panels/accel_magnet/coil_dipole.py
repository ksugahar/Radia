"""
Racetrack coil for C-type dipole magnet.

Coil wraps around the main leg (return path) of the C-type yoke.
Matches the C-type electromagnet racetrack-coil geometry.

Coil geometry (mm):
    Center: y = 151.25 mm from gap face (around main leg)
    Cross-section: 35 mm (radial) x 105 mm (axial)
    Straight sections: 62.5 mm (beam direction)
    Arc radius: 22.5 mm (centerline, inner=5, outer=40)
    Current: 2000 A-turns

Coordinate system (matching yoke.jou):
    y = gap direction (pole tip at y~0, main leg at y~0.15)
    z = beam direction (symmetry at z=0)
    x = depth/thickness (symmetry at x=0)

Usage:
    Set this file as "Coil script" in the Accelerator Magnet panel.
    The panel calls build_coil() to get wire segments for Biot-Savart.
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                '..', '..', '..', 'src', 'radia'))

from coil_builder import CoilBuilder


def build_coil():
    """Build racetrack coil matching the C-type dipole geometry.

    The coil wraps around the main leg of the C-type yoke.
    Two straight sections run along the beam (z) direction,
    connected by semicircular arcs at each end.

    The coil center is at y = 151.25 mm from the gap face,
    which is the center of the main leg in the yoke.

    Returns:
        CoilBuilder instance (upper coil only; mirror for lower)
    """
    mm = 1e-3

    # ELF coil dimensions
    NI = 2000                         # Ampere-turns
    COIL_RADIAL = 35 * mm             # cross-section radial width
    COIL_AXIAL = 105 * mm             # cross-section axial height (x-dir)
    STRAIGHT_LEN = 62.5 * mm          # straight section length (beam, z)
    ARC_R_INNER = 5 * mm              # inner arc radius
    ARC_R_OUTER = 40 * mm             # outer arc radius
    ARC_R_MEAN = (ARC_R_INNER + ARC_R_OUTER) / 2  # 22.5 mm centerline

    # Coil center position (from gap face)
    # ELF: coil center at (0, 131.25, 0) mm in original coords
    # In panel coords: y = 151.25 mm (main leg center)
    Y_CENTER = 151.25 * mm

    # Straight section positions:
    #   Inner straight: y = Y_CENTER - ARC_R_MEAN = 128.75 mm (toward pole)
    #   Outer straight: y = Y_CENTER + ARC_R_MEAN = 173.75 mm (behind leg)
    Y_INNER = Y_CENTER - ARC_R_MEAN

    # Coil start: inner straight, beam start
    # Straight runs along z from -STRAIGHT_LEN/2 to +STRAIGHT_LEN/2
    # For quarter model (z > 0), the full coil still provides correct field
    # Racetrack: 4 straights + 4 quarter-circle arcs (tilt=0 for XY plane)
    X_CAP = 50 * mm              # X-cap length (from ELF: x=-25 to +25)
    X_CENTER = 47.5 * mm         # inner Y-straight x position

    coil = (CoilBuilder(current=NI)
        .set_start([X_CENTER, Y_INNER, 0])
        .set_cross_section(width=COIL_RADIAL, height=COIL_AXIAL)
        .add_straight(STRAIGHT_LEN)            # Y up
        .add_arc(radius=ARC_R_MEAN, arc_angle=90)   # top-right
        .add_straight(X_CAP)                   # X-cap top (-x)
        .add_arc(radius=ARC_R_MEAN, arc_angle=90)   # top-left
        .add_straight(STRAIGHT_LEN)            # Y down
        .add_arc(radius=ARC_R_MEAN, arc_angle=90)   # bottom-left
        .add_straight(X_CAP)                   # X-cap bottom (+x)
        .add_arc(radius=ARC_R_MEAN, arc_angle=90))   # bottom-right (close)

    return coil


if __name__ == "__main__":
    # Quick test: print wire segments and write STEP
    coil = build_coil()
    segments, current = coil.to_wire_segments()
    print(f"Wire segments: {len(segments)}, current: {current} A")
    print(f"Gap: {coil.gap:.6e} m, closed: {coil.is_closed}")

    mm = 1e-3
    print(f"\nCoil center: y = {151.25 * mm * 1000:.2f} mm from gap face")
    print(f"Inner straight: y = {(151.25 - 22.5) * mm * 1000:.2f} mm")
    print(f"Outer straight: y = {(151.25 + 22.5) * mm * 1000:.2f} mm")

    # STEP export for visualization
    try:
        coil.write_step("coil_dipole.step")
        print("Wrote coil_dipole.step")
    except Exception as e:
        print(f"STEP export skipped: {e}")
