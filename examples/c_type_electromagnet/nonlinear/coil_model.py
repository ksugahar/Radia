#!/usr/bin/env python
"""
Common racetrack coil model for electromagnet examples.

Matches ELF geometry from:
S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/mu=1000/ELF_MMB8T_EIEM2_1x1x1

Coil geometry (from MCL8T elements):
- Center: (0, 131.25, 0) mm
- Cross-section: 35 mm (radial) x 105 mm (axial)
- Straight sections: 62.5 mm (y-direction), 50 mm (x-direction)
- Arc radii: 5 mm (inner), 40 mm (outer)
- Height: z = -52.5 to +52.5 mm (105 mm)
- Current: 2000 AT
"""

import sys
import os

# Add radia module path
work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(work_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad
import numpy as np


def create_racetrack_coil(current_at=2000.0, n_seg=12):
    """
    Create racetrack coil matching ELF geometry.

    ELF MCL8T elements have full conductor cross-section (35mm x 105mm).
    Arc centers are at (±25, 100) and (±25, 162.5) mm.

    Parameters:
        current_at: Ampere-turns (default 2000)
        n_seg: Number of segments for arc approximation

    Returns:
        Radia coil object

    Coil geometry (in mm):
        Cross-section: 35 mm (radial) x 105 mm (axial)
        Y-straights: x = 30 to 65, y = 100 to 162.5
        Arc centers: (±25, 100) and (±25, 162.5)
        Inner arc radius: 5 mm, outer arc radius: 40 mm
        X-caps: from x = -25 to +25, y = 167.5/202.5 (top) and 95/60 (bottom)
    """
    # Coil geometry in meters
    width = 0.035   # 35 mm (radial: 65 - 30)
    height = 0.105  # 105 mm (axial: -52.5 to +52.5)

    # Conductor centerline position for Y-straights
    x_center = 0.0475  # 47.5 mm = (30 + 65) / 2

    # Y extent of Y-straights (from ELF nodes)
    y_straight_bottom = 0.100   # 100 mm
    y_straight_top = 0.1625     # 162.5 mm
    y_mid = (y_straight_bottom + y_straight_top) / 2  # 131.25 mm
    straight_y = y_straight_top - y_straight_bottom   # 62.5 mm

    # Arc geometry (ELF specification)
    r_inner = 0.005  # 5 mm inner arc radius
    r_outer = 0.040  # 40 mm outer arc radius
    r_mean = (r_inner + r_outer) / 2  # 22.5 mm (centerline arc radius)

    # Arc centers (from ELF node positions)
    arc_x = 0.025  # 25 mm = 30 - 5 = 65 - 40

    # X-cap positions (from ELF node positions)
    y_cap_top_inner = 0.1675    # 167.5 mm (inner edge of top cap)
    y_cap_top_outer = 0.2025    # 202.5 mm (outer edge of top cap)
    y_cap_bottom_inner = 0.095  # 95 mm (inner edge of bottom cap)
    y_cap_bottom_outer = 0.060  # 60 mm (outer edge of bottom cap)
    y_cap_top = (y_cap_top_inner + y_cap_top_outer) / 2    # 185 mm
    y_cap_bottom = (y_cap_bottom_inner + y_cap_bottom_outer) / 2  # 77.5 mm

    # X-cap length (from ELF: x = -25 to +25)
    straight_x = 0.050  # 50 mm

    # Cross-section area for current density
    cross_section = width * height  # 35 mm x 105 mm = 0.003675 m^2

    # Current density (A/m^2)
    j = current_at / cross_section  # 2000 / 0.003675 = 544217.7 A/m^2

    coil_objects = []

    # For -z field at origin, current flows COUNTER-CLOCKWISE when viewed from +z
    # (ELF MCL8T convention: current flows from nodes 1-4 to nodes 5-8)
    # Path: Y-straight(+x,up) -> arc(TL) -> X-cap(top,-x) -> arc(TR) ->
    #       Y-straight(-x,down) -> arc(BR) -> X-cap(bot,+x) -> arc(BL)

    # =========================================================================
    # Y-direction straights (main coil legs)
    # Full conductor cross-section (35mm x 105mm) at centerline x = +/-47.5
    # =========================================================================

    # Positive x side: current in +y direction (up) - matches ELF MCL8T
    seg_y_pos = rad.ObjRecCur([x_center, y_mid, 0],
                              [width, straight_y, height], [0, j, 0])
    coil_objects.append(seg_y_pos)

    # Negative x side: current in -y direction (down) - matches ELF MCL8T
    seg_y_neg = rad.ObjRecCur([-x_center, y_mid, 0],
                              [width, straight_y, height], [0, -j, 0])
    coil_objects.append(seg_y_neg)

    # =========================================================================
    # Arc sections at corners (approximated by n_seg straight segments)
    # Arc centers are at (±25, 100) and (±25, 162.5)
    # Centerline follows arc of radius 22.5mm
    # =========================================================================

    # Four corners with their arc centers and angle ranges
    # Current flows COUNTER-CLOCKWISE (viewed from +z), so phi INCREASES
    # Path: Y(+x,+y) -> arc(TL) -> X-cap(top,-x) -> arc(TR) -> Y(-x,-y) ->
    #       arc(BR) -> X-cap(bot,+x) -> arc(BL) -> back to Y(+x)
    corner_configs = [
        # (arc_center_x, arc_center_y, start_angle, end_angle)
        # top-right: from Y-straight (47.5,162.5) to X-cap (-25,185) via arc
        (+arc_x, y_straight_top, 0, 90),      # phi: 0 -> 90 (counter-clockwise)
        # top-left: from X-cap (-25,185) to Y-straight (-47.5,162.5)
        (-arc_x, y_straight_top, 90, 180),    # phi: 90 -> 180 (counter-clockwise)
        # bottom-left: from Y-straight (-47.5,100) to X-cap (25,77.5)
        (-arc_x, y_straight_bottom, 180, 270),  # phi: 180 -> 270 (counter-clockwise)
        # bottom-right: from X-cap (25,77.5) to Y-straight (47.5,100)
        (+arc_x, y_straight_bottom, -90, 0),  # phi: -90 -> 0 (counter-clockwise)
    ]

    for arc_cx, arc_cy, phi_start_deg, phi_end_deg in corner_configs:
        phi_start = np.radians(phi_start_deg)
        phi_end = np.radians(phi_end_deg)

        for i in range(n_seg):
            phi1 = phi_start + (phi_end - phi_start) * i / n_seg
            phi2 = phi_start + (phi_end - phi_start) * (i + 1) / n_seg
            phi_mid = (phi1 + phi2) / 2

            # Segment center position (on centerline arc)
            seg_x = arc_cx + r_mean * np.cos(phi_mid)
            seg_y = arc_cy + r_mean * np.sin(phi_mid)

            # Arc length for this segment
            arc_len = r_mean * abs(phi2 - phi1)

            # Current direction: tangent to arc in counter-clockwise direction
            # For counter-clockwise flow (phi increasing), tangent = (-sin(phi), cos(phi))
            # This gives +y at phi=0, -x at phi=90, etc.
            jx = -j * np.sin(phi_mid)
            jy = j * np.cos(phi_mid)

            seg = rad.ObjRecCur([seg_x, seg_y, 0],
                                [width, arc_len, height], [jx, jy, 0])
            coil_objects.append(seg)

    # =========================================================================
    # X-direction caps (at top and bottom of coil)
    # X-cap spans from x = -25 to +25 (50mm), NOT from -47.5 to +47.5
    # =========================================================================

    # Top cap: current in -x direction (left) - matches ELF counter-clockwise
    seg_x_top = rad.ObjRecCur([0, y_cap_top, 0],
                              [straight_x, width, height], [-j, 0, 0])
    coil_objects.append(seg_x_top)

    # Bottom cap: current in +x direction (right) - matches ELF counter-clockwise
    seg_x_bot = rad.ObjRecCur([0, y_cap_bottom, 0],
                              [straight_x, width, height], [j, 0, 0])
    coil_objects.append(seg_x_bot)

    # Arc sections (8 total: 4 inner, 4 outer)
    # Use ObjArcCur for arc current blocks
    # ObjArcCur(center, radii, phi_range, height, n_seg, man_auto, j_phi)

    # Inner arcs (radius 5 mm, from r=25 to r=30)
    # Arc centers are at the corners: (25, 162.5), (-25, 162.5), (-25, 100), (25, 100)

    # Corner positions for inner arcs
    inner_arc_centers = [
        ([0.025, 0.1625, 0], [0, 90]),    # bottom-right inner
        ([0.025, 0.1625, 0], [90, 180]),   # top-right inner
        ([-0.025, 0.1625, 0], [0, 90]),   # top-left inner
        ([-0.025, 0.1625, 0], [270, 360]), # bottom-left inner
    ]

    # For now, skip the arcs and use the simplified model
    # The straight sections capture the main field contribution

    return rad.ObjCnt(coil_objects)


def get_coil_field_at_origin(coil):
    """Get the field from the coil only (no iron) at the origin."""
    B = rad.Fld(coil, 'b', [0, 0, 0])
    return np.array(B)


if __name__ == "__main__":
    # Test the coil model
    rad.UtiDelAll()

    print("=" * 60)
    print("Racetrack Coil Model Test")
    print("=" * 60)

    print("\nCoil geometry:")
    print("  Center: (0, 131.25, 0) mm")
    print("  Cross-section: 35 mm x 105 mm")
    print("  Inner arc radius: 5 mm")
    print("  Outer arc radius: 40 mm")
    print("  X-straights: 50 mm")
    print("  Y-straights: 62.5 mm")

    # Create coil
    coil = create_racetrack_coil(2000.0)

    # Get coil-only field at origin
    B_coil = get_coil_field_at_origin(coil)
    print(f"\nCoil-only field at (0,0,0):")
    print(f"  Bx = {B_coil[0]*1000:.4f} mT")
    print(f"  By = {B_coil[1]*1000:.4f} mT")
    print(f"  Bz = {B_coil[2]*1000:.4f} mT")
    print(f"  |B| = {np.linalg.norm(B_coil)*1000:.4f} mT")

    # Expected: Small field at origin (coil is far from origin)
    # The iron yoke will amplify this significantly
    print("\nNote: The coil produces a small field at origin.")
    print("With mu_r=1000 iron yoke, the field is amplified ~200x.")
    print("Expected total field: ~228 mT (from ELF)")
