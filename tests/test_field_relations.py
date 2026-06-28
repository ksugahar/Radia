#!/usr/bin/env python
"""
Tests for Maxwell field relations: curl(A) = B, H = -grad(Phi)

These tests verify that the A (vector potential) and Phi (scalar potential)
fields are consistent with B and H fields through the Maxwell relations:
  - B = curl(A)  (in appropriate units)
  - H = -grad(Phi) (in appropriate units)

Note on Radia units:
- Radia uses internal units where A = (1/4pi) * M x BufVect
- The curl(A)/B ratio is approximately 1/mu_0 = 7.96e5
- Similarly for H = -grad(Phi)

Author: Radia Development Team
Date: 2025-12-31
"""
import sys
import os
import pytest
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'radia'))
import radia as rad

# Constants
MU_0 = 4 * np.pi * 1e-7
Br = 1.2  # Tesla
Mr = Br / MU_0  # A/m


def numerical_curl(field_func, point, h=1e-6):
    """Compute curl of a vector field numerically using central differences.

    curl(F) = [dFz/dy - dFy/dz, dFx/dz - dFz/dx, dFy/dx - dFx/dy]
    """
    x, y, z = point

    # dF/dx
    Fp_x = field_func([x + h, y, z])
    Fm_x = field_func([x - h, y, z])
    dFx_dx = (Fp_x[0] - Fm_x[0]) / (2 * h)
    dFy_dx = (Fp_x[1] - Fm_x[1]) / (2 * h)
    dFz_dx = (Fp_x[2] - Fm_x[2]) / (2 * h)

    # dF/dy
    Fp_y = field_func([x, y + h, z])
    Fm_y = field_func([x, y - h, z])
    dFx_dy = (Fp_y[0] - Fm_y[0]) / (2 * h)
    dFy_dy = (Fp_y[1] - Fm_y[1]) / (2 * h)
    dFz_dy = (Fp_y[2] - Fm_y[2]) / (2 * h)

    # dF/dz
    Fp_z = field_func([x, y, z + h])
    Fm_z = field_func([x, y, z - h])
    dFx_dz = (Fp_z[0] - Fm_z[0]) / (2 * h)
    dFy_dz = (Fp_z[1] - Fm_z[1]) / (2 * h)
    dFz_dz = (Fp_z[2] - Fm_z[2]) / (2 * h)

    # curl = [dFz/dy - dFy/dz, dFx/dz - dFz/dx, dFy/dx - dFx/dy]
    curl_x = dFz_dy - dFy_dz
    curl_y = dFx_dz - dFz_dx
    curl_z = dFy_dx - dFx_dy

    return np.array([curl_x, curl_y, curl_z])


def numerical_grad(scalar_func, point, h=1e-6):
    """Compute gradient of a scalar field numerically using central differences.

    grad(f) = [df/dx, df/dy, df/dz]
    """
    x, y, z = point

    df_dx = (scalar_func([x + h, y, z]) - scalar_func([x - h, y, z])) / (2 * h)
    df_dy = (scalar_func([x, y + h, z]) - scalar_func([x, y - h, z])) / (2 * h)
    df_dz = (scalar_func([x, y, z + h]) - scalar_func([x, y, z - h])) / (2 * h)

    return np.array([df_dx, df_dy, df_dz])


class TestCurlAEqualsB:
    """Test that curl(A) is proportional to B."""

    def setup_method(self):
        rad.UtiDelAll()

    def teardown_method(self):
        rad.UtiDelAll()

    def test_curl_a_proportional_to_b_hexahedron(self):
        """Test curl(A) / B ratio is consistent for ObjHexahedron."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        def A_func(pt):
            return rad.Fld(hex_mag, 'a', pt)

        # Test at multiple points
        test_points = [
            [0.05, 0.03, 0.02],
            [0.04, 0.04, 0.04],
            [0.06, 0.02, 0.05],
            [0.03, 0.05, 0.03],
        ]

        ratios = []
        for pt in test_points:
            curl_A = numerical_curl(A_func, pt, h=1e-5)
            B = np.array(rad.Fld(hex_mag, 'b', pt))

            curl_A_mag = np.linalg.norm(curl_A)
            B_mag = np.linalg.norm(B)

            if B_mag > 1e-6:  # Only compute ratio where B is significant
                ratio = curl_A_mag / B_mag
                ratios.append(ratio)

        # All ratios should be consistent (same scaling factor)
        if len(ratios) > 1:
            ratio_mean = np.mean(ratios)
            ratio_std = np.std(ratios)
            rel_std = ratio_std / ratio_mean

            assert rel_std < 0.1, \
                f"curl(A)/B ratios should be consistent: mean={ratio_mean:.2e}, std={ratio_std:.2e}, rel_std={rel_std:.2%}"

    def test_curl_a_proportional_to_b_recmag(self):
        """Test curl(A) / B ratio is consistent for a rectangular PM (magnet_box)."""
        dx, dy, dz = 0.02, 0.02, 0.03

        rec_mag = rad.magnet_box([0, 0, 0], [2*dx, 2*dy, 2*dz], [0, 0, Mr])

        def A_func(pt):
            return rad.Fld(rec_mag, 'a', pt)

        test_points = [
            [0.05, 0.03, 0.02],
            [0.04, 0.04, 0.04],
            [0.06, 0.02, 0.05],
        ]

        ratios = []
        for pt in test_points:
            curl_A = numerical_curl(A_func, pt, h=1e-5)
            B = np.array(rad.Fld(rec_mag, 'b', pt))

            curl_A_mag = np.linalg.norm(curl_A)
            B_mag = np.linalg.norm(B)

            if B_mag > 1e-6:
                ratio = curl_A_mag / B_mag
                ratios.append(ratio)

        if len(ratios) > 1:
            ratio_mean = np.mean(ratios)
            ratio_std = np.std(ratios)
            rel_std = ratio_std / ratio_mean

            assert rel_std < 0.1, \
                f"curl(A)/B ratios should be consistent: mean={ratio_mean:.2e}, rel_std={rel_std:.2%}"

    def test_curl_a_direction_matches_b(self):
        """Test that curl(A) direction matches B direction."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        def A_func(pt):
            return rad.Fld(hex_mag, 'a', pt)

        test_points = [
            [0.05, 0.03, 0.02],
            [0.04, 0.04, 0.04],
        ]

        for pt in test_points:
            curl_A = numerical_curl(A_func, pt, h=1e-5)
            B = np.array(rad.Fld(hex_mag, 'b', pt))

            curl_A_mag = np.linalg.norm(curl_A)
            B_mag = np.linalg.norm(B)

            if curl_A_mag > 1e-6 and B_mag > 1e-6:
                # Normalize and compare directions
                curl_A_norm = curl_A / curl_A_mag
                B_norm = B / B_mag

                # Dot product should be close to +1 (same direction)
                dot = np.dot(curl_A_norm, B_norm)
                assert dot > 0.95, \
                    f"curl(A) and B should point in same direction at {pt}: dot={dot}"


class TestGradPhiEqualsMinusH:
    """Test that -grad(Phi) is proportional to H."""

    def setup_method(self):
        rad.UtiDelAll()

    def teardown_method(self):
        rad.UtiDelAll()

    def test_grad_phi_proportional_to_h_hexahedron(self):
        """Test -grad(Phi) / H ratio is consistent for ObjHexahedron."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        def Phi_func(pt):
            return rad.Fld(hex_mag, 'p', pt)

        # Test at multiple points - avoid symmetry axes where Phi ~ 0
        test_points = [
            [0.05, 0.03, 0.05],
            [0.04, 0.04, 0.06],
            [0.06, 0.02, 0.07],
            [0.03, 0.05, 0.05],
        ]

        ratios = []
        for pt in test_points:
            grad_Phi = numerical_grad(Phi_func, pt, h=1e-5)
            H = np.array(rad.Fld(hex_mag, 'h', pt))

            grad_Phi_mag = np.linalg.norm(grad_Phi)
            H_mag = np.linalg.norm(H)

            if H_mag > 1e-6 and grad_Phi_mag > 1e3:  # Only where both are significant
                ratio = grad_Phi_mag / H_mag
                ratios.append(ratio)

        if len(ratios) > 1:
            ratio_mean = np.mean(ratios)
            ratio_std = np.std(ratios)
            rel_std = ratio_std / ratio_mean

            assert rel_std < 0.15, \
                f"-grad(Phi)/H ratios should be consistent: mean={ratio_mean:.2e}, rel_std={rel_std:.2%}"

    def test_grad_phi_proportional_to_h_recmag(self):
        """Test -grad(Phi) / H ratio is consistent for a rectangular PM (magnet_box)."""
        dx, dy, dz = 0.02, 0.02, 0.03

        rec_mag = rad.magnet_box([0, 0, 0], [2*dx, 2*dy, 2*dz], [0, 0, Mr])

        def Phi_func(pt):
            return rad.Fld(rec_mag, 'p', pt)

        test_points = [
            [0.05, 0.03, 0.05],
            [0.04, 0.04, 0.06],
            [0.06, 0.02, 0.07],
        ]

        ratios = []
        for pt in test_points:
            grad_Phi = numerical_grad(Phi_func, pt, h=1e-5)
            H = np.array(rad.Fld(rec_mag, 'h', pt))

            grad_Phi_mag = np.linalg.norm(grad_Phi)
            H_mag = np.linalg.norm(H)

            if H_mag > 1e-6 and grad_Phi_mag > 1e3:
                ratio = grad_Phi_mag / H_mag
                ratios.append(ratio)

        if len(ratios) > 1:
            ratio_mean = np.mean(ratios)
            ratio_std = np.std(ratios)
            rel_std = ratio_std / ratio_mean

            assert rel_std < 0.15, \
                f"-grad(Phi)/H ratios should be consistent: mean={ratio_mean:.2e}, rel_std={rel_std:.2%}"

    def test_grad_phi_direction_opposite_to_h(self):
        """Test that -grad(Phi) direction matches H direction."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        def Phi_func(pt):
            return rad.Fld(hex_mag, 'p', pt)

        # Points off symmetry axes where Phi is significant
        test_points = [
            [0.04, 0.03, 0.06],
            [0.05, 0.04, 0.07],
        ]

        for pt in test_points:
            grad_Phi = numerical_grad(Phi_func, pt, h=1e-5)
            H = np.array(rad.Fld(hex_mag, 'h', pt))

            grad_Phi_mag = np.linalg.norm(grad_Phi)
            H_mag = np.linalg.norm(H)

            if grad_Phi_mag > 1e3 and H_mag > 1e-6:
                # -grad(Phi) should point in same direction as H
                minus_grad_Phi_norm = -grad_Phi / grad_Phi_mag
                H_norm = H / H_mag

                dot = np.dot(minus_grad_Phi_norm, H_norm)
                assert dot > 0.9, \
                    f"-grad(Phi) and H should point in same direction at {pt}: dot={dot}"


class TestBHRelation:
    """Test B = mu_0 * H relation in air region."""

    def setup_method(self):
        rad.UtiDelAll()

    def teardown_method(self):
        rad.UtiDelAll()

    def test_b_equals_mu0_h_in_air(self):
        """Test B = mu_0 * H in air region (outside magnet)."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        # Test at points clearly outside the magnet
        test_points = [
            [0.05, 0.0, 0.0],
            [0.0, 0.05, 0.0],
            [0.0, 0.0, 0.05],
            [0.04, 0.04, 0.05],
        ]

        for pt in test_points:
            B = np.array(rad.Fld(hex_mag, 'b', pt))
            H = np.array(rad.Fld(hex_mag, 'h', pt))

            B_mag = np.linalg.norm(B)
            H_mag = np.linalg.norm(H)

            if H_mag > 1e-6:
                # B = mu_0 * H in air
                expected_B_mag = MU_0 * H_mag
                rel_diff = abs(B_mag - expected_B_mag) / expected_B_mag

                assert rel_diff < 0.01, \
                    f"B != mu_0*H at {pt}: B={B_mag:.6e}, mu_0*H={expected_B_mag:.6e}, diff={rel_diff:.2%}"

    def test_b_h_direction_same_in_air(self):
        """Test B and H point in same direction in air."""
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        test_points = [
            [0.05, 0.03, 0.02],
            [0.04, 0.04, 0.04],
        ]

        for pt in test_points:
            B = np.array(rad.Fld(hex_mag, 'b', pt))
            H = np.array(rad.Fld(hex_mag, 'h', pt))

            B_mag = np.linalg.norm(B)
            H_mag = np.linalg.norm(H)

            if B_mag > 1e-10 and H_mag > 1e-6:
                B_norm = B / B_mag
                H_norm = H / H_mag

                dot = np.dot(B_norm, H_norm)
                assert dot > 0.999, \
                    f"B and H should point in same direction at {pt}: dot={dot}"


class TestFieldRelationsConsistency:
    """Test overall consistency of field relations."""

    def setup_method(self):
        rad.UtiDelAll()

    def teardown_method(self):
        rad.UtiDelAll()

    def test_recmag_hexahedron_ratio_consistency(self):
        """Test that curl(A)/B ratios are same for magnet_box and a hand-built ObjHexahedron."""
        dx, dy, dz = 0.02, 0.02, 0.03

        rec_mag = rad.magnet_box([0, 0, 0], [2*dx, 2*dy, 2*dz], [0, 0, Mr])
        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        def A_rec(pt):
            return rad.Fld(rec_mag, 'a', pt)

        def A_hex(pt):
            return rad.Fld(hex_mag, 'a', pt)

        pt = [0.05, 0.03, 0.04]

        curl_A_rec = numerical_curl(A_rec, pt, h=1e-5)
        curl_A_hex = numerical_curl(A_hex, pt, h=1e-5)

        B_rec = np.array(rad.Fld(rec_mag, 'b', pt))
        B_hex = np.array(rad.Fld(hex_mag, 'b', pt))

        ratio_rec = np.linalg.norm(curl_A_rec) / np.linalg.norm(B_rec)
        ratio_hex = np.linalg.norm(curl_A_hex) / np.linalg.norm(B_hex)

        # Ratios should be similar
        rel_diff = abs(ratio_rec - ratio_hex) / ratio_rec
        assert rel_diff < 0.1, \
            f"curl(A)/B ratios should be similar: rec={ratio_rec:.2e}, hex={ratio_hex:.2e}"

    def test_scaling_factor_consistent(self):
        """Test that curl(A)/B ratio is consistent and positive.

        Note: Radia uses an internal unit system where A = (1/4pi) * M x BufVect.
        The curl(A)/B ratio depends on this internal scaling and is not exactly 1/mu_0.
        This test verifies that the ratio is consistent and the relationship is valid.
        """
        dx, dy, dz = 0.02, 0.02, 0.03

        vertices = [
            [-dx, -dy, -dz], [dx, -dy, -dz], [dx, dy, -dz], [-dx, dy, -dz],
            [-dx, -dy, dz], [dx, -dy, dz], [dx, dy, dz], [-dx, dy, dz],
        ]
        hex_mag = rad.ObjHexahedron(vertices, [0, 0, Mr])

        def A_func(pt):
            return rad.Fld(hex_mag, 'a', pt)

        # Test at multiple points
        test_points = [
            [0.05, 0.03, 0.04],
            [0.04, 0.04, 0.05],
            [0.06, 0.02, 0.03],
        ]

        ratios = []
        for pt in test_points:
            curl_A = numerical_curl(A_func, pt, h=1e-5)
            B = np.array(rad.Fld(hex_mag, 'b', pt))

            curl_A_mag = np.linalg.norm(curl_A)
            B_mag = np.linalg.norm(B)

            if B_mag > 1e-6:
                ratio = curl_A_mag / B_mag
                ratios.append(ratio)

        # Verify:
        # 1. Ratio is positive
        # 2. Ratio is consistent across different points
        assert all(r > 0 for r in ratios), "curl(A)/B ratio should be positive"

        if len(ratios) > 1:
            ratio_mean = np.mean(ratios)
            ratio_std = np.std(ratios)
            rel_std = ratio_std / ratio_mean

            assert rel_std < 0.15, \
                f"curl(A)/B ratio should be consistent: mean={ratio_mean:.2e}, rel_std={rel_std:.2%}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
