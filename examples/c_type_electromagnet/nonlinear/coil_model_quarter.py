#!/usr/bin/env python
"""
Quarter racetrack coil model for C-type electromagnet.

This creates a quarter coil (x >= 0 quadrant only) for use with IMA symmetry.
The coil parameters match the ELF model:
- 20000 AT (Ampere-turns)
- Racetrack shape with straight sections and semicircular ends
"""

import numpy as np
import radia as rad


def create_coil_segment(center, size, j_vec):
    """Create a single rectangular current segment."""
    return rad.ObjRecCur(center, size, j_vec)


def create_racetrack_coil_quarter(current_at=20000.0, n_arc_seg=4):
    """
    Create quarter racetrack coil (x >= 0 only).

    For use with IMA X symmetry - the coil is only in x >= 0 region.
    The IMA will NOT mirror the coil (only magnetic materials are mirrored),
    so this quarter coil produces only 1/4 of the full coil field.

    Parameters
    ----------
    current_at : float
        Total current in Ampere-turns (default: 20000 AT)
    n_arc_seg : int
        Number of segments for arc approximation (default: 4)

    Returns
    -------
    int
        Radia object handle for the coil container
    """
    # Coil dimensions from ELF (in meters)
    x_inner = 0.030   # Inner radius in x
    x_outer = 0.065   # Outer radius in x
    width = x_outer - x_inner
    x_center = (x_inner + x_outer) / 2

    z_bottom = 0.0
    z_top = 0.0525
    height = z_top - z_bottom
    z_center = (z_bottom + z_top) / 2

    y_bottom = 0.100
    y_top = 0.1625
    straight_len = y_top - y_bottom
    y_mid = (y_bottom + y_top) / 2

    r_inner = 0.005
    r_outer = 0.040
    r_mean = (r_inner + r_outer) / 2
    arc_center_x = 0.025

    y_cap_top = 0.185
    y_cap_bottom = 0.0775

    cross_section = width * height
    j = current_at / cross_section  # Current density A/m^2

    coil_objects = []

    # Helper to add segment with z-mirror (coil symmetric in z)
    def add_segment_z(cx, cy, cz, sx, sy, sz, jx, jy, jz):
        # Original (z >= 0)
        coil_objects.append(create_coil_segment([cx, cy, cz], [sx, sy, sz], [jx, jy, jz]))
        # z-mirror (z <= 0): same current direction (coil symmetric)
        coil_objects.append(create_coil_segment([cx, cy, -cz], [sx, sy, sz], [jx, jy, jz]))

    # Y-straight section (only x >= 0)
    add_segment_z(x_center, y_mid, z_center, width, straight_len, height, 0, j, 0)

    # Top arcs (quarter circle, 0 to 90 degrees)
    for i in range(n_arc_seg):
        phi1 = np.radians(0 + 90 * i / n_arc_seg)
        phi2 = np.radians(0 + 90 * (i + 1) / n_arc_seg)
        phi_mid = (phi1 + phi2) / 2

        seg_x = arc_center_x + r_mean * np.cos(phi_mid)
        seg_y = y_top + r_mean * np.sin(phi_mid)
        arc_len = r_mean * abs(phi2 - phi1)

        jx = -j * np.sin(phi_mid)
        jy = j * np.cos(phi_mid)

        add_segment_z(seg_x, seg_y, z_center, width, arc_len, height, jx, jy, 0)

    # Top cap (partial, x >= 0 only)
    cap_x_center = arc_center_x / 2
    cap_x_len = arc_center_x
    add_segment_z(cap_x_center, y_cap_top, z_center, cap_x_len, width, height, -j, 0, 0)

    # Bottom arcs (quarter circle, -90 to 0 degrees)
    for i in range(n_arc_seg):
        phi1 = np.radians(-90 + 90 * i / n_arc_seg)
        phi2 = np.radians(-90 + 90 * (i + 1) / n_arc_seg)
        phi_mid = (phi1 + phi2) / 2

        seg_x = arc_center_x + r_mean * np.cos(phi_mid)
        seg_y = y_bottom + r_mean * np.sin(phi_mid)
        arc_len = r_mean * abs(phi2 - phi1)

        jx = -j * np.sin(phi_mid)
        jy = j * np.cos(phi_mid)

        add_segment_z(seg_x, seg_y, z_center, width, arc_len, height, jx, jy, 0)

    # Bottom cap (partial, x >= 0 only)
    add_segment_z(cap_x_center, y_cap_bottom, z_center, cap_x_len, width, height, j, 0, 0)

    return rad.ObjCnt(coil_objects)


if __name__ == '__main__':
    # Test coil creation
    rad.UtiDelAll()
    rad.FldUnits('m')

    coil = create_racetrack_coil_quarter(20000.0)

    # Check field at origin
    B = rad.Fld(coil, 'b', [0, 0, 0])
    print(f"Quarter coil field at (0,0,0):")
    print(f"  Bx = {B[0]*1000:.4f} mT")
    print(f"  By = {B[1]*1000:.4f} mT")
    print(f"  Bz = {B[2]*1000:.4f} mT")
    print(f"  |B| = {np.linalg.norm(B)*1000:.4f} mT")
