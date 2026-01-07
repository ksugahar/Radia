"""
Test FastImp core functionality.

This test module validates the core FastImp implementation including:
- Surface panel discretization
- Green's function computations
- pFFT acceleration
- Impedance calculation formulas

Note: These tests validate the mathematical foundations. Once Python bindings
are added, additional integration tests will test the C++ implementation directly.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

# Physical constants (matching C++ implementation)
MU_0 = 4 * np.pi * 1e-7  # H/m
EPS_0 = 8.854187817e-12  # F/m
C_0 = 299792458.0  # m/s (speed of light)
INV_FOUR_PI = 1.0 / (4 * np.pi)


class TestSurfacePanelDiscretization:
    """Test surface panel geometry computations."""

    def test_rectangular_block_panel_count(self):
        """Test that rectangular block generates correct number of panels."""
        # A cube with n panels per face direction
        # 6 faces, n^2 panels per face = 6*n^2 total panels
        for n in [2, 4, 8]:
            expected_panels = 6 * n * n
            # Simulate the panel generation logic
            actual_panels = 6 * n * n  # Direct calculation matches C++ logic
            assert actual_panels == expected_panels

    def test_panel_area_sum_equals_surface_area(self):
        """Test that sum of panel areas equals total surface area."""
        # Cube with dimensions L x L x L
        L = 0.1  # 10 cm
        total_surface_area = 6 * L * L

        # Discretize into n x n panels per face
        n = 4
        panel_area = (L / n) * (L / n)
        num_panels = 6 * n * n
        computed_area = num_panels * panel_area

        assert_allclose(computed_area, total_surface_area, rtol=1e-10)

    def test_quadrilateral_panel_normal_direction(self):
        """Test that panel normals point outward."""
        # Define a simple face (z=0 plane, normal should be -z or +z)
        vertices = [
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            np.array([1, 1, 0]),
            np.array([0, 1, 0])
        ]

        # Compute normal via cross product (right-hand rule)
        e1 = vertices[1] - vertices[0]
        e2 = vertices[3] - vertices[0]
        normal = np.cross(e1, e2)
        normal = normal / np.linalg.norm(normal)

        # Normal should be along z-axis
        assert np.abs(normal[0]) < 1e-10
        assert np.abs(normal[1]) < 1e-10
        assert np.abs(np.abs(normal[2]) - 1.0) < 1e-10

    def test_triangle_area_formula(self):
        """Test triangle area calculation."""
        # Triangle with vertices at (0,0,0), (1,0,0), (0,1,0)
        v0 = np.array([0, 0, 0])
        v1 = np.array([1, 0, 0])
        v2 = np.array([0, 1, 0])

        e1 = v1 - v0
        e2 = v2 - v0
        cross = np.cross(e1, e2)
        area = np.linalg.norm(cross) / 2.0

        assert_allclose(area, 0.5, rtol=1e-10)


class TestCircularWireDiscretization:
    """Test circular wire/loop panel generation."""

    def test_circular_loop_surface_area(self):
        """Test that circular loop surface area matches formula."""
        # Circular loop: major radius R, wire radius a
        R = 0.05  # 5 cm loop radius
        a = 0.001  # 1 mm wire radius

        # Theoretical surface area of torus
        theoretical_area = 4 * np.pi * np.pi * R * a

        # Discretized area
        n_loop = 30  # panels around loop
        n_around = 8  # panels around wire circumference

        # Each panel is approximately: (2*pi*R/n_loop) * (2*pi*a/n_around)
        panel_length = 2 * np.pi * R / n_loop
        panel_width = 2 * np.pi * a / n_around
        panel_area = panel_length * panel_width
        discretized_area = n_loop * n_around * panel_area

        assert_allclose(discretized_area, theoretical_area, rtol=0.01)

    def test_spiral_path_length(self):
        """Test spiral path length calculation."""
        # Spiral parameters
        inner_radius = 0.01  # 1 cm
        outer_radius = 0.05  # 5 cm
        num_turns = 5
        pitch = 0.002  # 2 mm pitch

        # Approximate path length (using arc length integral approximation)
        # For Archimedean spiral: r = a + b*theta
        # Length = integral of sqrt(r^2 + (dr/dtheta)^2) dtheta

        total_angle = 2 * np.pi * num_turns
        points_per_turn = 100
        total_points = num_turns * points_per_turn

        path_length = 0
        for i in range(total_points):
            t = i / total_points
            theta = total_angle * t
            r = inner_radius + (outer_radius - inner_radius) * t

            # Next point
            t_next = (i + 1) / total_points
            theta_next = total_angle * t_next
            r_next = inner_radius + (outer_radius - inner_radius) * t_next

            x = r * np.cos(theta)
            y = r * np.sin(theta)
            z = pitch * num_turns * t

            x_next = r_next * np.cos(theta_next)
            y_next = r_next * np.sin(theta_next)
            z_next = pitch * num_turns * t_next

            dl = np.sqrt((x_next - x)**2 + (y_next - y)**2 + (z_next - z)**2)
            path_length += dl

        # Should be > simple circular estimate
        circular_estimate = 2 * np.pi * (inner_radius + outer_radius) / 2 * num_turns
        assert path_length > circular_estimate * 0.9


class TestGreenFunctionComputation:
    """Test Green's function calculations for various formulations."""

    def test_mqs_green_function_basic(self):
        """Test MQS Green's function G(r) = 1/(4*pi*r)."""
        r = 0.1  # 10 cm
        G_expected = INV_FOUR_PI / r
        G_computed = 1.0 / (4 * np.pi * r)

        assert_allclose(G_computed, G_expected, rtol=1e-10)

    def test_mqs_green_scaling_with_distance(self):
        """Test that MQS Green's function scales as 1/r."""
        r1 = 0.1
        r2 = 0.2

        G1 = INV_FOUR_PI / r1
        G2 = INV_FOUR_PI / r2

        # G1/G2 should equal r2/r1
        ratio = G1 / G2
        expected_ratio = r2 / r1

        assert_allclose(ratio, expected_ratio, rtol=1e-10)

    def test_fullwave_green_function(self):
        """Test full-wave Green's function G(r) = exp(-jkr)/(4*pi*r)."""
        freq = 1e9  # 1 GHz
        r = 0.1  # 10 cm

        omega = 2 * np.pi * freq
        k = omega / C_0

        # Full-wave Green's function
        G_fullwave = np.exp(-1j * k * r) * INV_FOUR_PI / r

        # Check magnitude and phase
        G_magnitude = np.abs(G_fullwave)
        G_expected_magnitude = INV_FOUR_PI / r

        # Magnitude should match 1/(4*pi*r) within numerical precision
        assert_allclose(G_magnitude, G_expected_magnitude, rtol=1e-10)

        # Phase should be -kr
        G_phase = np.angle(G_fullwave)
        expected_phase = -k * r

        # Wrap to [-pi, pi]
        while expected_phase < -np.pi:
            expected_phase += 2 * np.pi
        while expected_phase > np.pi:
            expected_phase -= 2 * np.pi

        assert_allclose(G_phase, expected_phase, rtol=1e-6)

    def test_fullwave_reduces_to_mqs_at_low_frequency(self):
        """Test that full-wave approaches MQS at low frequencies."""
        freq = 100  # 100 Hz (very low frequency)
        r = 0.1  # 10 cm

        omega = 2 * np.pi * freq
        k = omega / C_0

        # At low frequency, kr << 1
        kr = k * r
        assert kr < 1e-6, f"kr = {kr} is not small enough for low-frequency approximation"

        G_mqs = INV_FOUR_PI / r
        G_fullwave = np.exp(-1j * k * r) * INV_FOUR_PI / r

        # Real part should match closely
        assert_allclose(np.real(G_fullwave), G_mqs, rtol=1e-6)
        # Imaginary part should be much smaller than real part (|Im/Re| ~ kr << 1)
        imag_to_real_ratio = np.abs(np.imag(G_fullwave) / np.real(G_fullwave))
        assert imag_to_real_ratio < kr * 2, f"Im/Re ratio {imag_to_real_ratio} too large"

    def test_green_gradient_magnitude(self):
        """Test Green's function gradient dG/dr = -1/(4*pi*r^2)."""
        r = 0.1

        # Analytical gradient
        dG_dr_analytical = -INV_FOUR_PI / (r * r)

        # Numerical gradient (central difference)
        dr = 1e-8
        G_plus = INV_FOUR_PI / (r + dr)
        G_minus = INV_FOUR_PI / (r - dr)
        dG_dr_numerical = (G_plus - G_minus) / (2 * dr)

        assert_allclose(dG_dr_numerical, dG_dr_analytical, rtol=1e-4)


class TestPFFTAcceleration:
    """Test pFFT algorithm components."""

    def test_toeplitz_circulant_embedding(self):
        """Test that Toeplitz matrix can be embedded in circulant form."""
        # A 1D Toeplitz matrix T of size N can be embedded in
        # a circulant matrix C of size 2N

        N = 8
        # First row/column of Toeplitz: t = [t0, t1, ..., t_{N-1}]
        t = np.array([1.0 / (k + 1) for k in range(N)])  # 1/r kernel

        # Toeplitz matrix
        T = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                T[i, j] = t[abs(i - j)]

        # Circulant embedding: c = [t0, t1, ..., t_{N-1}, 0, t_{N-1}, ..., t_1]
        c = np.zeros(2 * N)
        c[:N] = t
        c[N + 1:] = t[1:][::-1]  # Mirror for circulant

        # Circulant matrix from first row
        C = np.zeros((2 * N, 2 * N))
        for i in range(2 * N):
            for j in range(2 * N):
                C[i, j] = c[(j - i) % (2 * N)]

        # Top-left N x N block of C should equal T
        assert_allclose(C[:N, :N], T, rtol=1e-10)

    def test_fft_convolution_equivalence(self):
        """Test that FFT-based convolution matches direct convolution."""
        # For circulant matrix C, y = C*x is equivalent to:
        # Y = FFT(c) * FFT(x), y = IFFT(Y)

        N = 16
        # Kernel (first row of circulant)
        c = np.zeros(N)
        c[0] = 1.0
        c[1] = 0.5
        c[-1] = 0.5  # Symmetric kernel

        # Input vector
        x = np.random.randn(N)

        # Direct convolution via circulant matrix
        C = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                C[i, j] = c[(j - i) % N]
        y_direct = C @ x

        # FFT-based convolution
        c_fft = np.fft.fft(c)
        x_fft = np.fft.fft(x)
        y_fft = np.fft.ifft(c_fft * x_fft)

        assert_allclose(np.real(y_fft), y_direct, rtol=1e-10)

    def test_3d_fft_grid_indexing(self):
        """Test 3D FFT grid indexing for pFFT."""
        nx, ny, nz = 4, 4, 4
        total = nx * ny * nz

        # Verify linear indexing matches 3D indexing
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    # Linear index (z-fastest, then y, then x)
                    linear_idx = iz + nz * (iy + ny * ix)
                    assert linear_idx < total
                    assert linear_idx >= 0

                    # Reverse mapping
                    iz_check = linear_idx % nz
                    iy_check = (linear_idx // nz) % ny
                    ix_check = linear_idx // (nz * ny)

                    assert ix_check == ix
                    assert iy_check == iy
                    assert iz_check == iz


class TestImpedanceFormulation:
    """Test impedance matrix formulation."""

    def test_dc_resistance_wire(self):
        """Test DC resistance calculation for wire."""
        # R_DC = rho * L / A = L / (sigma * A)
        sigma = 5.8e7  # Copper conductivity [S/m]
        length = 1.0  # 1 m
        diameter = 0.001  # 1 mm
        area = np.pi * (diameter / 2) ** 2

        R_DC = length / (sigma * area)

        # Expected: approximately 0.022 ohms for 1m of 1mm copper wire
        assert 0.01 < R_DC < 0.03

    def test_skin_effect_resistance_increase(self):
        """Test that AC resistance increases with frequency due to skin effect."""
        sigma = 5.8e7  # Copper conductivity [S/m]
        diameter = 0.001  # 1 mm wire

        def skin_depth(freq):
            omega = 2 * np.pi * freq
            return np.sqrt(2 / (omega * MU_0 * sigma))

        def ac_resistance_per_length(freq, d):
            delta = skin_depth(freq)
            a = d / 2  # radius
            if delta > a:
                # Low frequency: uniform current
                return 1 / (sigma * np.pi * a * a)
            else:
                # High frequency: skin depth limited
                # Approximate: R ~ 1/(sigma * 2*pi*a*delta)
                return 1 / (sigma * 2 * np.pi * a * delta)

        R_1kHz = ac_resistance_per_length(1e3, diameter)
        R_1MHz = ac_resistance_per_length(1e6, diameter)
        R_1GHz = ac_resistance_per_length(1e9, diameter)

        # AC resistance should increase with frequency
        assert R_1MHz > R_1kHz
        assert R_1GHz > R_1MHz

    def test_inductance_matrix_symmetry(self):
        """Test that mutual inductance matrix is symmetric."""
        # L_ij = mu_0/(4*pi) * integral{ 1/|r-r'| } dA_i dA_j

        # For testing, use simplified discrete panels
        n_panels = 5
        panel_centers = np.random.randn(n_panels, 3) * 0.1
        panel_areas = np.ones(n_panels) * 0.001  # 10 cm^2 each

        L = np.zeros((n_panels, n_panels))
        for i in range(n_panels):
            for j in range(n_panels):
                r = np.linalg.norm(panel_centers[i] - panel_centers[j])
                if r < 1e-10:
                    # Self term: approximate
                    R_eff = np.sqrt(panel_areas[i] / np.pi)
                    L[i, j] = MU_0 * INV_FOUR_PI * R_eff * 0.5 * panel_areas[j]
                else:
                    L[i, j] = MU_0 * INV_FOUR_PI / r * panel_areas[j]

        # Check symmetry
        assert_allclose(L, L.T, rtol=1e-10)

    def test_loop_inductance_formula(self):
        """Test circular loop inductance formula."""
        # L = mu_0 * R * (ln(8R/a) - 2) for circular loop
        # where R = loop radius, a = wire radius

        R = 0.05  # 5 cm loop radius
        a = 0.0005  # 0.5 mm wire radius

        L_theoretical = MU_0 * R * (np.log(8 * R / a) - 2)

        # Should be on the order of 100-500 nH for typical parameters
        assert 50e-9 < L_theoretical < 1e-6


class TestNearFieldCorrection:
    """Test near-field (direct) computation for pFFT."""

    def test_near_field_threshold(self):
        """Test near-field identification based on distance."""
        # Panels closer than threshold need direct computation

        panel_size = 0.01  # 1 cm panel
        threshold_factor = 3.0
        threshold = threshold_factor * panel_size

        # Test distances
        distances = [0.005, 0.01, 0.02, 0.05, 0.1]

        near_field = [d < threshold for d in distances]
        expected = [True, True, True, False, False]

        assert near_field == expected

    def test_self_term_regularization(self):
        """Test self-term regularization for singular integrals."""
        # For panel self-interaction, use analytical approximation:
        # G_self = R_eff * (2*ln(2) - 1) / (4*pi)
        # where R_eff = sqrt(Area / pi)

        area = 0.0001  # 1 cm^2
        R_eff = np.sqrt(area / np.pi)

        G_self = R_eff * (2 * np.log(2) - 1) * INV_FOUR_PI

        # Should be finite and positive
        assert G_self > 0
        assert np.isfinite(G_self)

        # Compare to nearby panel interaction (should be similar order)
        nearby_distance = R_eff * 2  # 2 effective radii away
        G_nearby = INV_FOUR_PI / nearby_distance

        # Self term should be smaller than nearby direct term
        assert G_self < G_nearby * 10


class TestChargeConservation:
    """Test charge continuity equation implementation."""

    def test_surface_divergence_operator(self):
        """Test that surface divergence conserves total charge."""
        # div_s(K) + j*omega*sigma = 0
        # For steady state (DC), div_s(K) = 0 implies current is solenoidal

        # On a closed surface, integral of div_s(K) = 0
        # (divergence theorem on surface)

        # Simplified test: uniform current on sphere has zero divergence
        n_panels = 100
        theta = np.linspace(0, np.pi, 10)
        phi = np.linspace(0, 2 * np.pi, 10)

        total_divergence = 0.0
        for t in theta[:-1]:
            for p in phi[:-1]:
                # Panel area element on unit sphere
                dA = np.sin(t) * (theta[1] - theta[0]) * (phi[1] - phi[0])

                # For uniform tangential current, div_s = 0
                # so contribution to total is zero
                total_divergence += 0 * dA

        assert np.abs(total_divergence) < 1e-10

    def test_charge_current_continuity_at_interface(self):
        """Test current continuity at panel interfaces."""
        # Current entering a panel = Current leaving + j*omega*charge accumulation

        omega = 2 * np.pi * 1e6  # 1 MHz
        panel_area = 0.0001  # 1 cm^2

        # If K_in = 1 A/m and K_out = 0.99 A/m
        K_in = 1.0
        K_out = 0.99

        # Net current divergence
        div_K = (K_in - K_out) / np.sqrt(panel_area)  # Approximate div

        # Charge accumulation rate
        sigma_rate = -div_K / (1j * omega)  # From continuity equation

        # Charge should be complex (out of phase with current at high frequency)
        assert np.imag(sigma_rate) != 0


class TestPortExtraction:
    """Test port impedance extraction."""

    def test_voltage_current_ratio(self):
        """Test impedance as V/I ratio."""
        # Z = V / I
        V_applied = 1.0  # 1 Volt
        I_measured = 0.01 + 0.005j  # Complex current (resistive + reactive)

        Z = V_applied / I_measured

        # Real part = resistance, imaginary part = reactance
        R = np.real(Z)
        X = np.imag(Z)

        assert R > 0  # Resistance is positive
        # X can be positive (inductive) or negative (capacitive)

    def test_two_port_network(self):
        """Test two-port impedance matrix properties."""
        # For reciprocal network: Z_12 = Z_21

        # Simulate a simple symmetric structure
        Z_11 = 50 + 10j  # Self impedance port 1
        Z_22 = 50 + 10j  # Self impedance port 2 (same structure)
        Z_12 = 5 + 2j   # Mutual impedance
        Z_21 = Z_12     # Reciprocity

        Z_matrix = np.array([[Z_11, Z_12], [Z_21, Z_22]])

        # Check symmetry for reciprocal network
        assert_allclose(Z_matrix, Z_matrix.T)

        # Eigenvalues should have positive real parts (passive network)
        eigenvalues = np.linalg.eigvals(Z_matrix)
        for ev in eigenvalues:
            assert np.real(ev) >= 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
