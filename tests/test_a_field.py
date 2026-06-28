#!/usr/bin/env python
"""
Tests for vector potential (A) field computation.

Tests verify that:
1. magnet_box (MMMM rectangular PM) and a hand-built ObjHexahedron give consistent A values off-axis
2. A is computed using face-based integration
3. A field symmetry properties are correct
4. A = 0 on symmetry axes is physically correct for face-based method

Author: Radia Development Team
Date: 2025-12-31
"""
import sys
import os
import pytest
import numpy as np

# Add src/radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'radia'))
import radia as rad

# Constants
MU_0 = 4 * np.pi * 1e-7
Br = 1.2  # Tesla
Mr = Br / MU_0  # A/m


class TestAFieldBasic:
    """Basic tests for A (vector potential) field computation."""

    def setup_method(self):
        """Setup before each test."""
        rad.UtiDelAll()

    def teardown_method(self):
        """Cleanup after each test."""
        rad.UtiDelAll()

    def test_a_off_axis_recmag_vs_hexahedron(self):
        """Test that A matches between magnet_box and a hand-built ObjHexahedron off-axis.

        magnet_box is now the rectangular permanent magnet (MMMM surface-charge ObjHexahedron),
        so this cross-checks the box-helper's vertex winding against a hand-built hexahedron.
        """
        dx, dy, dz = 0.02, 0.02, 0.03

        rec_mag = rad.magnet_box([0, 0, 0], [2*dx, 2*dy, 2*dz], [0, 0, Mr])
        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # Off-axis points (not on x, y, or z axis)
        test_points = [
            [0.05, 0.0, 0.0],   # x-axis
            [0.0, 0.05, 0.0],   # y-axis
            [0.03, 0.03, 0.03], # diagonal
            [0.04, 0.02, 0.05], # general
        ]

        for pt in test_points:
            A_rec = rad.Fld(rec_mag, 'a', pt)
            A_hex = rad.Fld(hex_mag, 'a', pt)

            A_rec_mag = np.sqrt(A_rec[0]**2 + A_rec[1]**2 + A_rec[2]**2)
            A_hex_mag = np.sqrt(A_hex[0]**2 + A_hex[1]**2 + A_hex[2]**2)

            if A_rec_mag > 1e6:  # Only compare where A is significant
                rel_diff = abs(A_rec_mag - A_hex_mag) / A_rec_mag
                assert rel_diff < 0.01, \
                    f"A mismatch at {pt}: rel_diff={rel_diff*100:.2f}%"

    def test_a_nonzero_off_axis(self):
        """Test that A is non-zero at off-axis points."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # A should be non-zero at off-axis points
        test_points = [
            [0.05, 0.0, 0.0],   # x-axis
            [0.0, 0.05, 0.0],   # y-axis
            [0.03, 0.03, 0.03], # diagonal
        ]

        for pt in test_points:
            A = rad.Fld(hex_mag, 'a', pt)
            A_mag = np.sqrt(A[0]**2 + A[1]**2 + A[2]**2)

            assert A_mag > 1e-6, f"|A| should be non-zero at {pt}, got {A_mag}"

    def test_a_decay_with_distance(self):
        """Test that A decays with distance from magnet."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # A should decrease with distance
        distances = [0.05, 0.08, 0.12, 0.16]
        A_values = []
        for d in distances:
            A = rad.Fld(hex_mag, 'a', [d, 0, 0])  # Along x-axis
            A_mag = np.sqrt(A[0]**2 + A[1]**2 + A[2]**2)
            A_values.append(A_mag)

        for i in range(len(A_values) - 1):
            assert A_values[i] > A_values[i+1], \
                f"|A| should decay: {A_values[i]} <= {A_values[i+1]}"


class TestAFieldSymmetryAxis:
    """Tests for A field behavior on symmetry axes."""

    def setup_method(self):
        rad.UtiDelAll()

    def teardown_method(self):
        rad.UtiDelAll()

    def test_a_zero_on_z_axis_face_based(self):
        """Test that A is zero on z-axis for face-based method.

        For a z-magnetized block, the vector potential A on the z-axis
        should be zero due to symmetry. This is because:
        - A = (1/4pi) * M x BufVect
        - On z-axis, BufVect has only z-component by symmetry
        - M = [0, 0, Mz], so M x BufVect = [0, 0, 0]
        """
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # A should be zero (or very small) on z-axis
        test_points = [
            [0.0, 0.0, 0.05],
            [0.0, 0.0, 0.035],
            [0.0, 0.0, 0.10],
        ]

        for pt in test_points:
            A = rad.Fld(hex_mag, 'a', pt)
            A_mag = np.sqrt(A[0]**2 + A[1]**2 + A[2]**2)

            # A should be essentially zero (numerical precision)
            assert A_mag < 1e-6, \
                f"|A| on z-axis should be ~0 for face-based method, got {A_mag} at {pt}"

    def test_a_perpendicular_to_m_direction(self):
        """Test that A is perpendicular to magnetization direction.

        For uniform M, A = (1/4pi) * M x BufVect, so A should be
        perpendicular to M at all points.
        """
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # M is along z, so A should have no z-component (Az = 0)
        test_points = [
            [0.05, 0.0, 0.0],
            [0.0, 0.05, 0.0],
            [0.03, 0.03, 0.03],
        ]

        for pt in test_points:
            A = rad.Fld(hex_mag, 'a', pt)

            # Az should be very small compared to Ax, Ay
            A_xy = np.sqrt(A[0]**2 + A[1]**2)
            if A_xy > 1e6:  # Only check if A is significant
                assert abs(A[2]) < 0.01 * A_xy, \
                    f"Az should be ~0 for z-magnetized block at {pt}: Az={A[2]}, |A_xy|={A_xy}"


class TestAFieldDifferentMagnetizations:
    """Tests for A field with different magnetization directions."""

    def setup_method(self):
        rad.UtiDelAll()

    def teardown_method(self):
        rad.UtiDelAll()

    def test_a_x_magnetized_block(self):
        """Test A field for x-magnetized block."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        # M along x-axis
        hex_mag = rad.ObjHexahedron(vertices, [Mr, 0, 0])

        # A should be perpendicular to M (in y-z plane)
        pt = [0.05, 0.03, 0.03]
        A = rad.Fld(hex_mag, 'a', pt)

        # Ax should be small compared to Ay, Az
        A_yz = np.sqrt(A[1]**2 + A[2]**2)
        if A_yz > 1e6:
            assert abs(A[0]) < 0.01 * A_yz, \
                f"Ax should be ~0 for x-magnetized block: Ax={A[0]}, |A_yz|={A_yz}"

    def test_a_rotational_symmetry(self):
        """Test A field rotational symmetry.

        For a z-magnetized cube (dx=dy), A at (r, 0, z) should have same
        magnitude as A at (0, r, z) due to symmetry.
        """
        d = 0.02  # Same dx, dy for symmetry

        vertices = [
            [-d, -d, -d], [d, -d, -d], [d, d, -d], [-d, d, -d],
            [-d, -d, d], [d, -d, d], [d, d, d], [-d, d, d],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # Points at same distance from z-axis
        r = 0.05
        z = 0.03

        A1 = rad.Fld(hex_mag, 'a', [r, 0, z])
        A2 = rad.Fld(hex_mag, 'a', [0, r, z])

        A1_mag = np.sqrt(A1[0]**2 + A1[1]**2 + A1[2]**2)
        A2_mag = np.sqrt(A2[0]**2 + A2[1]**2 + A2[2]**2)

        rel_diff = abs(A1_mag - A2_mag) / max(A1_mag, 1e-15)
        assert rel_diff < 0.01, \
            f"|A| should be symmetric: |A1|={A1_mag}, |A2|={A2_mag}"


class TestAFieldArbitraryHexahedron:
    """Tests for A field from arbitrary (non-rectangular) hexahedra."""

    def setup_method(self):
        rad.UtiDelAll()

    def teardown_method(self):
        rad.UtiDelAll()

    def test_a_skewed_hexahedron(self):
        """Test A field from a skewed hexahedron."""
        # Skewed hexahedron
        vertices = [
            [-0.02, -0.02, -0.03],
            [0.02, -0.02, -0.03],
            [0.02, 0.02, -0.03],
            [-0.02, 0.02, -0.03],
            [-0.01, -0.02, 0.03],  # Top shifted
            [0.03, -0.02, 0.03],
            [0.03, 0.02, 0.03],
            [-0.01, 0.02, 0.03],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # A should be non-zero at off-axis point
        pt = [0.05, 0.03, 0.03]
        A = rad.Fld(hex_mag, 'a', pt)
        A_mag = np.sqrt(A[0]**2 + A[1]**2 + A[2]**2)

        assert A_mag > 1e-6, f"|A| should be non-zero for skewed hex, got {A_mag}"

    def test_a_truncated_pyramid(self):
        """Test A field from a truncated pyramid."""
        s = 0.5  # Scale factor for top
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz],
            [dx, -dy, -dz],
            [dx, dy, -dz],
            [-dx, dy, -dz],
            [-dx*s, -dy*s, dz],
            [dx*s, -dy*s, dz],
            [dx*s, dy*s, dz],
            [-dx*s, dy*s, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # A should be non-zero
        pt = [0.05, 0.03, 0.03]
        A = rad.Fld(hex_mag, 'a', pt)
        A_mag = np.sqrt(A[0]**2 + A[1]**2 + A[2]**2)

        assert A_mag > 1e-6, f"|A| should be non-zero for truncated pyramid, got {A_mag}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
