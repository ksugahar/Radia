#!/usr/bin/env python
"""
Tests for scalar potential (Phi) field computation.

Tests verify that:
1. magnet_box (rectangular fixed-M PM) and a hand-built ObjHexahedron give consistent Phi values
2. Phi is computed using face-based integration (not dipole approximation)
3. Phi symmetry properties are correct

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


class TestPhiFieldBasic:
    """Basic tests for Phi field computation."""

    def setup_method(self):
        """Setup before each test."""
        rad.UtiDelAll()

    def teardown_method(self):
        """Cleanup after each test."""
        rad.UtiDelAll()

    def test_phi_on_z_axis_recmag_vs_hexahedron(self):
        """Test that Phi matches between magnet_box and a hand-built ObjHexahedron on z-axis."""
        dx, dy, dz = 0.02, 0.02, 0.03

        # Create magnet_box (rectangular fixed-M PM = surface-charge ObjHexahedron)
        rec_mag = rad.magnet_box([0, 0, 0], [2*dx, 2*dy, 2*dz], [0, 0, Mr])

        # Create ObjHexahedron
        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # Test on z-axis (should match exactly)
        test_points_z = [
            [0.0, 0.0, 0.05],
            [0.0, 0.0, 0.035],
            [0.0, 0.0, 0.10],
            [0.0, 0.0, 0.15],
        ]

        for pt in test_points_z:
            Phi_rec = rad.Fld(rec_mag, 'p', pt)
            Phi_hex = rad.Fld(hex_mag, 'p', pt)

            # Should match within numerical precision
            assert abs(Phi_rec - Phi_hex) < 1e-4 * abs(Phi_rec) + 1e-10, \
                f"Phi mismatch at {pt}: rec={Phi_rec}, hex={Phi_hex}"

    def test_phi_diagonal_points(self):
        """Test Phi at diagonal points (not on symmetry axes)."""
        dx, dy, dz = 0.02, 0.02, 0.03

        rec_mag = rad.magnet_box([0, 0, 0], [2*dx, 2*dy, 2*dz], [0, 0, Mr])
        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # Diagonal points - should be close but may have small differences
        test_points = [
            [0.03, 0.03, 0.03],
            [0.04, 0.02, 0.05],
            [0.05, 0.05, 0.05],
        ]

        for pt in test_points:
            Phi_rec = rad.Fld(rec_mag, 'p', pt)
            Phi_hex = rad.Fld(hex_mag, 'p', pt)

            # Allow 2% difference for diagonal points
            rel_diff = abs(Phi_rec - Phi_hex) / max(abs(Phi_rec), 1e-15)
            assert rel_diff < 0.02, \
                f"Phi mismatch at {pt}: rec={Phi_rec}, hex={Phi_hex}, rel_diff={rel_diff*100:.2f}%"

    def test_phi_nonzero_on_z_axis(self):
        """Test that Phi is non-zero on z-axis for z-magnetized block."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # Phi should be non-zero on z-axis
        pt = [0.0, 0.0, 0.05]
        Phi = rad.Fld(hex_mag, 'p', pt)

        assert abs(Phi) > 1.0, f"Phi should be non-zero on z-axis, got {Phi}"


class TestPhiFieldSymmetry:
    """Tests for Phi field symmetry properties."""

    def setup_method(self):
        rad.UtiDelAll()

    def teardown_method(self):
        rad.UtiDelAll()

    def test_phi_symmetry_on_xy_axes(self):
        """Test Phi symmetry on x and y axes for z-magnetized block.

        For a z-magnetized block, Phi on x-axis and y-axis should be equal
        due to symmetry. Both values should be approximately zero because
        the contributions from top and bottom faces cancel.
        """
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # Phi on x-axis and y-axis should both be approximately zero
        # due to symmetry cancellation (top and bottom face contributions cancel)
        distances = [0.05, 0.06, 0.08, 0.10]
        for d in distances:
            Phi_x = rad.Fld(hex_mag, 'p', [d, 0, 0])
            Phi_y = rad.Fld(hex_mag, 'p', [0, d, 0])

            # Both should be very small (essentially zero due to cancellation)
            assert abs(Phi_x) < 1e-6, f"Phi on x-axis should be ~0, got {Phi_x}"
            assert abs(Phi_y) < 1e-6, f"Phi on y-axis should be ~0, got {Phi_y}"

    def test_phi_decay_with_distance(self):
        """Test that Phi decays with distance from magnet."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # Phi should decrease with distance on z-axis
        distances = [0.05, 0.10, 0.15, 0.20]
        Phi_values = []
        for z in distances:
            Phi = rad.Fld(hex_mag, 'p', [0, 0, z])
            Phi_values.append(abs(Phi))

        for i in range(len(Phi_values) - 1):
            assert Phi_values[i] > Phi_values[i+1], \
                f"Phi should decay: |Phi({distances[i]})| = {Phi_values[i]} <= |Phi({distances[i+1]})| = {Phi_values[i+1]}"


# Note: Use ObjTetrahedron and ObjHexahedron for creating elements.
# For mesh import, use netgen_mesh_import.py which creates elements
# using the internal ObjPolyhdr API.


class TestPhiFieldConsistency:
    """Tests for consistency between A, B, and Phi fields."""

    def setup_method(self):
        rad.UtiDelAll()

    def teardown_method(self):
        rad.UtiDelAll()

    def test_a_and_phi_both_nonzero(self):
        """Test that A and Phi are both computed (not zero) at off-axis points."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # At diagonal point, both A and Phi should be non-zero
        pt = [0.03, 0.03, 0.03]

        A = rad.Fld(hex_mag, 'a', pt)
        A_mag = np.sqrt(A[0]**2 + A[1]**2 + A[2]**2)

        Phi = rad.Fld(hex_mag, 'p', pt)

        assert A_mag > 1e-6, f"|A| should be non-zero at diagonal point, got {A_mag}"
        assert abs(Phi) > 1.0, f"|Phi| should be non-zero at diagonal point, got {Phi}"

    def test_b_field_consistency(self):
        """Test that B field is consistent between magnet_box and a hand-built ObjHexahedron."""
        dx, dy, dz = 0.02, 0.02, 0.03

        rec_mag = rad.magnet_box([0, 0, 0], [2*dx, 2*dy, 2*dz], [0, 0, Mr])
        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        test_points = [
            [0.05, 0.0, 0.0],
            [0.0, 0.0, 0.05],
            [0.03, 0.03, 0.03],
        ]

        for pt in test_points:
            B_rec = rad.Fld(rec_mag, 'b', pt)
            B_hex = rad.Fld(hex_mag, 'b', pt)

            B_rec_mag = np.sqrt(B_rec[0]**2 + B_rec[1]**2 + B_rec[2]**2)
            B_hex_mag = np.sqrt(B_hex[0]**2 + B_hex[1]**2 + B_hex[2]**2)

            rel_diff = abs(B_rec_mag - B_hex_mag) / max(B_rec_mag, 1e-15)
            assert rel_diff < 0.01, \
                f"B field mismatch at {pt}: rel_diff={rel_diff*100:.2f}%"


class TestPhiFieldArbitraryHexahedron:
    """Tests for Phi field from arbitrary (non-rectangular) hexahedra."""

    def setup_method(self):
        rad.UtiDelAll()

    def teardown_method(self):
        rad.UtiDelAll()

    def test_phi_skewed_hexahedron(self):
        """Test Phi from a skewed hexahedron."""
        # Skewed hexahedron (top shifted in x)
        vertices = [
            [-0.02, -0.02, -0.03],
            [0.02, -0.02, -0.03],
            [0.02, 0.02, -0.03],
            [-0.02, 0.02, -0.03],
            [-0.01, -0.02, 0.03],  # Shifted
            [0.03, -0.02, 0.03],   # Shifted
            [0.03, 0.02, 0.03],    # Shifted
            [-0.01, 0.02, 0.03],   # Shifted
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # Phi should be non-zero
        pt = [0.05, 0.0, 0.05]
        Phi = rad.Fld(hex_mag, 'p', pt)

        assert abs(Phi) > 1e3, f"Phi should be non-zero for skewed hex, got {Phi}"

    def test_phi_truncated_pyramid(self):
        """Test Phi from a truncated pyramid hexahedron."""
        # Truncated pyramid (top smaller than bottom)
        s = 0.5  # Scale factor for top
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz],
            [dx, -dy, -dz],
            [dx, dy, -dz],
            [-dx, dy, -dz],
            [-dx*s, -dy*s, dz],  # Smaller top
            [dx*s, -dy*s, dz],
            [dx*s, dy*s, dz],
            [-dx*s, dy*s, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # Phi should be non-zero
        pt = [0.05, 0.0, 0.05]
        Phi = rad.Fld(hex_mag, 'p', pt)

        assert abs(Phi) > 1.0, f"Phi should be non-zero for truncated pyramid, got {Phi}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
