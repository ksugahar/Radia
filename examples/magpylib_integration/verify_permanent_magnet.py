"""
Verify Radia Permanent Magnet vs Magpylib

Compare B field from a permanent magnet at multiple external observation points.
This test verifies that ObjPolyhdr with fixed magnetization produces correct fields.

Requirements:
- pip install magpylib>=5.0
- pip install radia (or use local build)
"""

import sys
import os
import numpy as np

# Add build directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad

try:
    import magpylib as magpy
    print("magpylib version:", magpy.__version__)
except ImportError:
    print("ERROR: magpylib not installed. Run: pip install magpylib")
    sys.exit(1)


def main():
    """Compare Radia permanent magnet with magpylib."""
    print("=" * 70)
    print("Verify Permanent Magnet: Radia vs Magpylib")
    print("=" * 70)

    # Use mm for Radia (default)
    rad.FldUnits('mm')
    rad.UtiDelAll()

    # Magnet dimensions (in mm)
    dx, dy, dz = 40, 40, 60  # 40x40x60 mm

    # Magnetization: Br = 1.2 T -> Mr = Br/mu_0 (A/m)
    MU_0 = 4 * np.pi * 1e-7
    Br = 1.2  # Tesla
    Mr = Br / MU_0  # A/m
    print("\nMagnet parameters:")
    print("  Dimensions: %d x %d x %d mm" % (dx, dy, dz))
    print("  Br = %.2f T" % Br)
    print("  Mr = %.0f A/m" % Mr)

    # =========================================================================
    # Create Radia permanent magnet using ObjPolyhdr
    # =========================================================================
    print("\n1. Creating Radia permanent magnet (ObjPolyhdr)...")

    # Hexahedron face definition (1-indexed vertices)
    HEX_FACES = [
        [1, 4, 3, 2],  # bottom (-z)
        [5, 6, 7, 8],  # top (+z)
        [1, 2, 6, 5],  # front (-y)
        [3, 4, 8, 7],  # back (+y)
        [1, 5, 8, 4],  # left (-x)
        [2, 3, 7, 6],  # right (+x)
    ]

    # Vertices for centered hexahedron
    half_x, half_y, half_z = dx/2, dy/2, dz/2
    vertices = [
        [-half_x, -half_y, -half_z],  # 1
        [+half_x, -half_y, -half_z],  # 2
        [+half_x, +half_y, -half_z],  # 3
        [-half_x, +half_y, -half_z],  # 4
        [-half_x, -half_y, +half_z],  # 5
        [+half_x, -half_y, +half_z],  # 6
        [+half_x, +half_y, +half_z],  # 7
        [-half_x, +half_y, +half_z],  # 8
    ]

    # Create permanent magnet with fixed magnetization [0, 0, Mr]
    radia_magnet = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, Mr])
    print("  Radia magnet created (ObjPolyhdr, ID: %d)" % radia_magnet)

    # =========================================================================
    # Create Magpylib permanent magnet
    # =========================================================================
    print("\n2. Creating magpylib permanent magnet...")

    # Magpylib uses meters and magnetization in A/m
    magpy_magnet = magpy.magnet.Cuboid(
        magnetization=(0, 0, Mr),  # A/m
        dimension=(dx/1000, dy/1000, dz/1000),  # m
        position=(0, 0, 0)  # m
    )
    print("  Magpylib magnet created")

    # =========================================================================
    # Compare B field at multiple external points
    # =========================================================================
    print("\n3. Comparing B field at external observation points...")
    print("   (All points are OUTSIDE the magnet)")
    print("")

    # Observation points (in mm for Radia, convert to m for magpylib)
    # Points are outside the magnet: magnet extends from -20 to +20 in x,y and -30 to +30 in z
    obs_points_mm = [
        # On z-axis (above/below)
        [0, 0, 40],    # 10mm above top face (z=30)
        [0, 0, 50],    # 20mm above top face
        [0, 0, 100],   # 70mm above top face
        [0, 0, -40],   # 10mm below bottom face
        [0, 0, -100],  # 70mm below bottom face
        # On x-axis
        [30, 0, 0],    # 10mm from side face (x=20)
        [50, 0, 0],    # 30mm from side face
        [100, 0, 0],   # 80mm from side face
        # Corner region
        [30, 30, 40],  # diagonal direction
        [50, 50, 50],  # further diagonal
    ]

    print("Point (mm)          |  Radia B (T)                |  Magpylib B (T)             |  Diff (%)")
    print("-" * 105)

    max_diff_percent = 0

    for pt_mm in obs_points_mm:
        # Radia: input in mm
        B_radia = rad.Fld(radia_magnet, 'b', pt_mm)

        # Magpylib: input in m
        pt_m = [p / 1000 for p in pt_mm]
        B_magpy = magpy_magnet.getB(pt_m)

        # Calculate difference
        B_radia_mag = np.sqrt(sum([b**2 for b in B_radia]))
        B_magpy_mag = np.sqrt(sum([b**2 for b in B_magpy]))

        if B_magpy_mag > 1e-10:
            diff_percent = abs(B_radia_mag - B_magpy_mag) / B_magpy_mag * 100
        else:
            diff_percent = 0

        max_diff_percent = max(max_diff_percent, diff_percent)

        print("[%4.0f,%4.0f,%4.0f]  |  [%+.4e, %+.4e, %+.4e]  |  [%+.4e, %+.4e, %+.4e]  |  %.2f%%" % (
            pt_mm[0], pt_mm[1], pt_mm[2],
            B_radia[0], B_radia[1], B_radia[2],
            B_magpy[0], B_magpy[1], B_magpy[2],
            diff_percent
        ))

    print("-" * 105)
    print("")
    print("Maximum difference: %.2f%%" % max_diff_percent)

    # =========================================================================
    # Verification result
    # =========================================================================
    print("")
    print("=" * 70)
    if max_diff_percent < 1.0:
        print("PASS: Radia permanent magnet matches magpylib within 1%%")
    elif max_diff_percent < 5.0:
        print("PASS: Radia permanent magnet matches magpylib within 5%%")
    else:
        print("WARNING: Significant difference between Radia and magpylib (%.1f%%)" % max_diff_percent)
    print("=" * 70)

    return max_diff_percent < 5.0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
