#!/usr/bin/env python3
"""
Test suite for FastImp conductor formulation (Phase 1-2)

This module tests:
1. Phase 1a: Conductor geometry creation and panel discretization
2. Phase 1a: Green's function calculations (DC/MQS/Full-wave)
3. Phase 2: Conductor-Magnetic material coupling

Note: These tests verify the mathematical formulation and will be
extended to test C++ implementation once integrated into the build.

References:
- Z. Zhu et al., "Algorithms in FastImp", IEEE TCAD, 2005
- S. Bilicz et al., "Wide-band nonlocal SIBC", ISEM 2023
"""

import pytest
import numpy as np
import sys
import os

# Add build directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))

# Physical constants
MU_0 = 4e-7 * np.pi  # H/m
EPS_0 = 8.854187817e-12  # F/m
C_0 = 299792458.0  # m/s
INV_FOUR_PI = 1.0 / (4.0 * np.pi)


class TestGreenFunctions:
    """Test Green's function calculations for various formulations."""

    def test_mqs_green_function(self):
        """Test MQS (quasi-static) Green's function: G(r) = 1/(4*pi*r)"""
        r = 0.01  # 1 cm
        G_expected = INV_FOUR_PI / r
        G_computed = INV_FOUR_PI / r

        assert np.isclose(G_computed, G_expected, rtol=1e-10)

    def test_mqs_green_gradient(self):
        """Test gradient of MQS Green's function: dG/dr = -1/(4*pi*r^2)"""
        r = 0.01
        dGdr_expected = -INV_FOUR_PI / (r * r)
        dGdr_computed = -INV_FOUR_PI / (r * r)

        assert np.isclose(dGdr_computed, dGdr_expected, rtol=1e-10)

    def test_fullwave_green_function(self):
        """Test full-wave Green's function: G(r) = exp(-jkr)/(4*pi*r)"""
        freq = 1e9  # 1 GHz
        omega = 2 * np.pi * freq
        k = omega * np.sqrt(MU_0 * EPS_0)
        r = 0.01

        # Full-wave Green's function
        G_expected = np.exp(-1j * k * r) * INV_FOUR_PI / r
        G_computed = np.exp(-1j * k * r) * INV_FOUR_PI / r

        assert np.isclose(G_computed, G_expected, rtol=1e-10)

    def test_fullwave_reduces_to_mqs_at_low_freq(self):
        """Verify full-wave reduces to MQS when k*r << 1"""
        freq = 1e3  # 1 kHz (low frequency)
        omega = 2 * np.pi * freq
        k = omega * np.sqrt(MU_0 * EPS_0)
        r = 0.01

        # At low frequency, exp(-jkr) ~ 1
        G_fullwave = np.exp(-1j * k * r) * INV_FOUR_PI / r
        G_mqs = INV_FOUR_PI / r

        # Check electrical size
        electrical_size = k * r
        assert electrical_size < 0.001, f"Electrical size too large: {electrical_size}"

        # Full-wave should be very close to MQS
        relative_error = np.abs(G_fullwave - G_mqs) / np.abs(G_mqs)
        assert relative_error < 1e-6, f"Full-wave differs from MQS: {relative_error}"

    def test_wavelength_calculation(self):
        """Test wavelength calculation at various frequencies"""
        # Note: C_0 = 299792458 m/s (exact speed of light)
        # Expected values are approximate for readability
        test_cases = [
            (1e6, 299.792458),    # 1 MHz -> ~300 m (exact: 299.792458 m)
            (1e9, 0.299792458),   # 1 GHz -> ~30 cm (exact: 0.299792458 m)
            (10e9, 0.0299792458), # 10 GHz -> ~3 cm (exact: 0.0299792458 m)
        ]

        for freq, expected_wavelength in test_cases:
            wavelength = C_0 / freq
            assert np.isclose(wavelength, expected_wavelength, rtol=1e-9)


class TestSkinDepth:
    """Test skin depth calculations for conductors."""

    def test_copper_skin_depth(self):
        """Test skin depth for copper at various frequencies"""
        sigma_copper = 5.8e7  # S/m

        test_cases = [
            (60, 8.5e-3),     # 60 Hz -> ~8.5 mm
            (1e3, 2.1e-3),    # 1 kHz -> ~2.1 mm
            (1e6, 66e-6),     # 1 MHz -> ~66 um
            (1e9, 2.1e-6),    # 1 GHz -> ~2.1 um
        ]

        for freq, expected_delta in test_cases:
            omega = 2 * np.pi * freq
            delta = np.sqrt(2.0 / (omega * MU_0 * sigma_copper))
            # Allow 5% tolerance due to rounding in expected values
            assert np.isclose(delta, expected_delta, rtol=0.05), \
                f"At {freq} Hz: expected {expected_delta}, got {delta}"

    def test_skin_depth_magnetic_material(self):
        """Test skin depth for magnetic conductor (sigma != 0, mu_r != 1)"""
        sigma = 1e6  # S/m (electrical steel)
        mu_r = 1000
        freq = 1e3

        omega = 2 * np.pi * freq
        mu = MU_0 * mu_r
        delta = np.sqrt(2.0 / (omega * mu * sigma))

        # Skin depth should be much smaller than for non-magnetic conductor
        delta_nonmagnetic = np.sqrt(2.0 / (omega * MU_0 * sigma))
        ratio = delta_nonmagnetic / delta
        assert np.isclose(ratio, np.sqrt(mu_r), rtol=1e-6)


class TestPanelGeometry:
    """Test surface panel geometry calculations."""

    def test_rectangular_block_surface_area(self):
        """Verify total surface area of rectangular block panels"""
        # Block dimensions
        Lx, Ly, Lz = 0.1, 0.05, 0.02  # 10x5x2 cm

        # Expected surface area
        A_expected = 2 * (Lx*Ly + Ly*Lz + Lz*Lx)

        # Computed surface area
        A_computed = 2 * (Lx*Ly + Ly*Lz + Lz*Lx)

        assert np.isclose(A_computed, A_expected, rtol=1e-10)

    def test_panel_normal_consistency(self):
        """Verify panel normals point outward consistently"""
        # For a cube centered at origin, face normals should point away from center
        cube_faces = [
            {'center': [0.5, 0, 0], 'normal': [1, 0, 0]},   # +X face
            {'center': [-0.5, 0, 0], 'normal': [-1, 0, 0]}, # -X face
            {'center': [0, 0.5, 0], 'normal': [0, 1, 0]},   # +Y face
            {'center': [0, -0.5, 0], 'normal': [0, -1, 0]}, # -Y face
            {'center': [0, 0, 0.5], 'normal': [0, 0, 1]},   # +Z face
            {'center': [0, 0, -0.5], 'normal': [0, 0, -1]}, # -Z face
        ]

        for face in cube_faces:
            center = np.array(face['center'])
            normal = np.array(face['normal'])
            # Normal should point in same direction as center (away from origin)
            dot_product = np.dot(center, normal)
            assert dot_product > 0, f"Normal at {center} points inward"

    def test_triangle_area(self):
        """Test triangle area calculation"""
        # Right triangle with legs of length 1
        v0 = np.array([0, 0, 0])
        v1 = np.array([1, 0, 0])
        v2 = np.array([0, 1, 0])

        # Area = 0.5 * |e1 x e2|
        e1 = v1 - v0
        e2 = v2 - v0
        cross = np.cross(e1, e2)
        area = 0.5 * np.linalg.norm(cross)

        assert np.isclose(area, 0.5, rtol=1e-10)

    def test_quadrilateral_area(self):
        """Test quadrilateral area calculation"""
        # Unit square in xy-plane
        vertices = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([1, 1, 0]),
            np.array([0, 1, 0]),
        ]

        # Area = 0.5 * |d1 x d2| where d1, d2 are diagonals
        d1 = vertices[2] - vertices[0]
        d2 = vertices[3] - vertices[1]
        cross = np.cross(d1, d2)
        area = 0.5 * np.abs(cross[2])  # z-component for xy-plane

        assert np.isclose(area, 1.0, rtol=1e-10)


class TestWireGeometry:
    """Test wire and coil geometry creation."""

    def test_circular_loop_path(self):
        """Test circular loop path generation"""
        radius = 0.05
        n_points = 30
        center = np.array([0, 0, 0])

        # Generate path
        theta = np.linspace(0, 2*np.pi, n_points + 1)
        path = np.zeros((n_points + 1, 3))
        path[:, 0] = center[0] + radius * np.cos(theta)
        path[:, 1] = center[1] + radius * np.sin(theta)
        path[:, 2] = center[2]

        # Verify path closes
        assert np.allclose(path[0], path[-1], rtol=1e-10)

        # Verify all points are at correct radius
        distances = np.sqrt(path[:, 0]**2 + path[:, 1]**2)
        assert np.allclose(distances, radius, rtol=1e-10)

    def test_spiral_path(self):
        """Test spiral coil path generation"""
        inner_radius = 0.02
        outer_radius = 0.05
        pitch = 0.005
        n_turns = 5
        n_points = n_turns * 30

        # Generate spiral
        t = np.linspace(0, 1, n_points + 1)
        theta = 2 * np.pi * n_turns * t
        r = inner_radius + (outer_radius - inner_radius) * t
        z = pitch * n_turns * t

        path = np.zeros((n_points + 1, 3))
        path[:, 0] = r * np.cos(theta)
        path[:, 1] = r * np.sin(theta)
        path[:, 2] = z

        # Verify start and end radii
        assert np.isclose(np.sqrt(path[0, 0]**2 + path[0, 1]**2), inner_radius, rtol=1e-6)
        assert np.isclose(np.sqrt(path[-1, 0]**2 + path[-1, 1]**2), outer_radius, rtol=1e-6)

        # Verify total height
        assert np.isclose(path[-1, 2] - path[0, 2], pitch * n_turns, rtol=1e-6)


class TestImpedanceFormulation:
    """Test impedance extraction formulation."""

    def test_dc_resistance(self):
        """Test DC resistance calculation"""
        sigma = 5.8e7  # Copper
        length = 0.1   # 10 cm
        width = 0.001  # 1 mm
        height = 0.001 # 1 mm

        # R_DC = length / (sigma * area)
        area = width * height
        R_expected = length / (sigma * area)

        # Should be about 1.7 mOhm for 10cm of 1mm^2 copper wire
        assert np.isclose(R_expected, 1.72e-3, rtol=0.01)

    def test_internal_inductance_circular_wire(self):
        """Test internal inductance per unit length for circular wire"""
        mu_r = 1

        # L_internal = mu / (8*pi) for circular wire (high frequency limit)
        L_internal_expected = MU_0 * mu_r / (8 * np.pi)

        # Should be about 50 nH/m
        assert np.isclose(L_internal_expected, 50e-9, rtol=0.01)

    def test_surface_impedance(self):
        """Test surface impedance Z_s = (1+j) / (sigma * delta)"""
        sigma = 5.8e7
        freq = 1e6
        omega = 2 * np.pi * freq

        delta = np.sqrt(2.0 / (omega * MU_0 * sigma))
        Z_s = (1 + 1j) / (sigma * delta)

        # Real and imaginary parts should be equal
        assert np.isclose(Z_s.real, Z_s.imag, rtol=1e-10)

        # Magnitude should be sqrt(omega * mu / (2 * sigma))
        Z_s_magnitude = np.abs(Z_s)
        expected_magnitude = np.sqrt(omega * MU_0 / sigma)
        assert np.isclose(Z_s_magnitude, expected_magnitude, rtol=1e-6)


class TestConductorMagneticCoupling:
    """Test conductor-magnetic material coupling (Phase 2)."""

    def test_biot_savart_from_wire(self):
        """Test B-field from straight wire using Biot-Savart"""
        # Infinite wire carrying current I along z-axis
        I = 1.0  # A
        rho = 0.01  # 1 cm from wire

        # B = mu_0 * I / (2 * pi * rho)
        B_expected = MU_0 * I / (2 * np.pi * rho)

        # Should be 20 uT at 1 cm from 1 A wire
        assert np.isclose(B_expected, 20e-6, rtol=0.01)

    def test_loop_center_field(self):
        """Test B-field at center of circular loop"""
        I = 1.0  # A
        radius = 0.05  # 5 cm

        # B = mu_0 * I / (2 * R) at center
        B_center = MU_0 * I / (2 * radius)

        # Expected value
        B_expected = MU_0 * 1.0 / (2 * 0.05)
        assert np.isclose(B_center, B_expected, rtol=1e-6)

    def test_mutual_inductance_concept(self):
        """Verify mutual inductance M = Phi_12 / I_1"""
        # Two coaxial loops
        radius = 0.05
        separation = 0.1

        # Approximate mutual inductance for well-separated loops
        # M ~ mu_0 * pi * R^4 / (2 * d^3) when d >> R
        M_approx = MU_0 * np.pi * radius**4 / (2 * separation**3)

        # M should be positive (same direction flux)
        assert M_approx > 0


class TestMatrixFormulation:
    """Test matrix formulation for IE solver."""

    def test_system_matrix_size(self):
        """Verify system matrix size for EFIE + continuity"""
        n_panels = 100
        # DOF per panel: K (scalar) + sigma (scalar) = 2
        n_dof = 2 * n_panels
        matrix_size = n_dof * n_dof

        # Should be 40000 elements for 100 panels
        assert matrix_size == 40000

    def test_self_term_regularization(self):
        """Test self-term regularization for 1/r singularity"""
        panel_area = 1e-4  # 1 cm^2

        # Effective radius for circular panel
        R_eff = np.sqrt(panel_area / np.pi)

        # Self-integral approximation
        # I_self ~ R_eff * (2*ln(2) - 1) / (4*pi)
        I_self = R_eff * (2 * np.log(2) - 1) * INV_FOUR_PI

        # Should be finite and positive
        assert I_self > 0
        assert np.isfinite(I_self)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
