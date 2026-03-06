"""
ESIM-based Conductor Model for PEEC.

This module implements the correct approach for using ESIM in PEEC:
1. Compute surface impedance Z_s from ESIM cell problem
2. Combine with bulk (DC) resistance for total impedance

The key formula:
    Z_conductor = R_dc + Z_skin

where:
    R_dc = rho * L / A                    (bulk DC resistance)
    Z_skin = Z_s * L / (effective_depth)  (skin effect contribution)

For planar conductor with half-thickness a:
    R_dc = rho * L / (w * 2a)
    Z_skin contribution depends on geometry

Author: Radia Development Team
Date: 2026-01-11
"""

import numpy as np
from scipy.constants import mu_0
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

from esim_cell_problem import ESIMFiniteSlabSolver, BHCurveInterpolator


class ESIMConductorModel:
    """
    ESIM-based conductor model for accurate skin effect in PEEC.

    This class combines:
    1. DC bulk resistance (always present)
    2. Skin effect contribution from ESIM surface impedance

    The model is valid for:
    - All frequency regimes (DC to high frequency)
    - Linear and nonlinear materials
    - Rectangular cross-section (planar approximation)
    """

    def __init__(self, width, height, length, sigma, mu_r=1.0, bh_curve=None):
        """
        Initialize ESIM conductor model.

        Parameters:
            width: Conductor width [m]
            height: Conductor height [m]
            length: Conductor length [m]
            sigma: Electrical conductivity [S/m]
            mu_r: Relative permeability (for linear material)
            bh_curve: BH curve data (for nonlinear material)
        """
        self.width = width
        self.height = height
        self.length = length
        self.sigma = sigma
        self.rho = 1.0 / sigma
        self.mu_r = mu_r
        self.bh_curve = bh_curve

        # Cross-section properties
        self.area = width * height
        self.perimeter = 2 * (width + height)

        # Use smaller dimension for skin effect calculation
        self.half_thickness = min(width, height) / 2

        # DC resistance
        self.R_dc = self.rho * length / self.area

    def compute_impedance(self, frequency, H0=None, method='hybrid'):
        """
        Compute conductor impedance at given frequency.

        Parameters:
            frequency: Operating frequency [Hz]
            H0: Surface field amplitude [A/m] (for nonlinear iteration)
            method: 'esim', 'dowell', or 'hybrid' (default)

        Returns:
            Z: Complex conductor impedance [Ohm]
            info: dict with additional information
        """
        if frequency < 1e-6:
            return self.R_dc, {'R_ac_Rdc': 1.0, 'xi': 0.0}

        # Compute skin depth and xi
        omega = 2 * np.pi * frequency
        mu = mu_0 * self.mu_r
        delta = np.sqrt(2 * self.rho / (omega * mu))
        xi = self.half_thickness / delta

        if method == 'dowell':
            return self._dowell_impedance(frequency, xi)
        elif method == 'esim':
            return self._esim_impedance(frequency, H0)
        else:  # hybrid
            return self._hybrid_impedance(frequency, xi, H0)

    def _dowell_impedance(self, frequency, xi):
        """Compute impedance using Dowell's formula."""
        R_ratio = self._dowell_ratio(xi)
        R_ac = self.R_dc * R_ratio

        # Inductance contribution (simplified)
        omega = 2 * np.pi * frequency
        delta = self.half_thickness / xi if xi > 0.01 else float('inf')

        # Internal inductance (approximate)
        # For full skin effect: L_int ~ mu * L / (perimeter * delta)
        # For DC: L_int ~ mu * L / (8 * pi) (round wire approx)
        if xi > 3:
            L_int = mu_0 * self.mu_r * self.length / (self.perimeter * delta)
        elif xi < 0.5:
            L_int = mu_0 * self.mu_r * self.length / (4 * np.pi)  # DC limit
        else:
            # Transition: interpolate
            L_high = mu_0 * self.mu_r * self.length / (self.perimeter * delta)
            L_low = mu_0 * self.mu_r * self.length / (4 * np.pi)
            t = (xi - 0.5) / 2.5
            L_int = (1-t) * L_low + t * L_high

        Z = R_ac + 1j * omega * L_int

        return Z, {'R_ac_Rdc': R_ratio, 'xi': xi, 'delta': delta, 'L_int': L_int}

    def _esim_impedance(self, frequency, H0=None):
        """Compute impedance using ESIM solver."""
        if H0 is None:
            H0 = 100.0  # Default surface field

        solver = ESIMFiniteSlabSolver(
            half_thickness=self.half_thickness,
            sigma=self.sigma,
            frequency=frequency,
            mu_r=self.mu_r,
            bh_curve=self.bh_curve,
            n_nodes=100
        )

        result = solver.solve(H0)
        Z_s = result['Z']
        xi = result['xi']

        # Convert surface impedance to conductor impedance
        # Z_conductor = R_dc + (Z_s - R_dc_surface) where R_dc_surface accounts
        # for the DC part already in R_dc

        # For the finite slab model:
        # Z_s gives impedance per unit area of conductor surface
        # The conductor has 2 surfaces (top and bottom) in parallel
        # Effective impedance: Z = Z_s * L / (2 * w) for width w

        # But this is complicated. Let's use the R_ac/R_dc ratio approach:
        # R_ratio = Re(Z_s) * a / rho
        R_ratio_esim = Z_s.real * self.half_thickness / self.rho

        # Ensure R_ac >= R_dc
        R_ratio = max(R_ratio_esim, 1.0)
        R_ac = self.R_dc * R_ratio

        # Reactive part from ESIM
        X_ratio = Z_s.imag * self.half_thickness / self.rho
        omega = 2 * np.pi * frequency
        L_int = X_ratio * self.R_dc / omega if omega > 0 else 0

        Z = R_ac + 1j * omega * L_int

        return Z, {'R_ac_Rdc': R_ratio, 'xi': xi, 'Z_s': Z_s, 'L_int': L_int}

    def _hybrid_impedance(self, frequency, xi, H0=None):
        """
        Hybrid approach: Dowell + ESIM correction for nonlinear.

        For linear materials: Use Dowell (analytical, fast)
        For nonlinear materials: Use ESIM (handles mu(H))
        """
        if self.bh_curve is not None:
            # Nonlinear: use ESIM
            return self._esim_impedance(frequency, H0)
        else:
            # Linear: use Dowell
            return self._dowell_impedance(frequency, xi)

    def _dowell_ratio(self, xi):
        """Dowell's formula for R_ac/R_dc."""
        if xi < 0.01:
            return 1.0 + xi**4 / 45

        try:
            sh2 = np.sinh(2 * xi)
            sn2 = np.sin(2 * xi)
            ch2 = np.cosh(2 * xi)
            cs2 = np.cos(2 * xi)

            denom = ch2 - cs2
            if abs(denom) < 1e-30:
                return 1.0

            return xi * (sh2 + sn2) / denom

        except (OverflowError, FloatingPointError):
            return xi  # High frequency limit


def test_esim_conductor_model():
    """Test the ESIM conductor model."""
    print("=" * 70)
    print("ESIM Conductor Model Test")
    print("=" * 70)

    # Copper conductor parameters
    width = 0.002     # 2 mm
    height = 0.001    # 1 mm
    length = 0.1      # 100 mm
    sigma = 5.8e7     # Copper

    model = ESIMConductorModel(width, height, length, sigma)

    print(f"\nCopper conductor:")
    print(f"  Dimensions: {width*1000:.1f} x {height*1000:.1f} x {length*1000:.1f} mm")
    print(f"  DC resistance: {model.R_dc*1000:.6f} mOhm")
    print()

    frequencies = [0, 100, 1000, 10000, 50000, 100000]

    print(f"{'Freq [Hz]':>12} {'xi':>8} {'R_ac [mOhm]':>14} {'R/Rdc Dow':>12} {'R/Rdc ESIM':>12}")
    print("-" * 64)

    for freq in frequencies:
        if freq == 0:
            print(f"{'DC':>12} {0:>8.3f} {model.R_dc*1000:>14.6f} {1.0:>12.4f} {1.0:>12.4f}")
            continue

        # Dowell method
        Z_dow, info_dow = model.compute_impedance(freq, method='dowell')

        # ESIM method
        Z_esim, info_esim = model.compute_impedance(freq, method='esim')

        print(f"{freq:>12.0f} {info_dow['xi']:>8.3f} {Z_dow.real*1000:>14.6f} "
              f"{info_dow['R_ac_Rdc']:>12.4f} {info_esim['R_ac_Rdc']:>12.4f}")

    print("\n" + "=" * 70)
    print("Nonlinear Material Test")
    print("=" * 70)

    # Steel BH curve
    bh_curve = [
        [0, 0],
        [100, 0.2],
        [500, 0.9],
        [1000, 1.3],
        [5000, 1.7],
        [10000, 1.9],
    ]

    sigma_steel = 2e6
    model_steel = ESIMConductorModel(width, height, length, sigma_steel, bh_curve=bh_curve)

    print(f"\nSteel conductor with BH curve:")
    print(f"  sigma = {sigma_steel/1e6:.1f} MS/m")
    print(f"  DC resistance: {model_steel.R_dc*1000:.6f} mOhm")
    print()

    freq = 50000  # 50 kHz
    H0_values = [10, 100, 1000, 5000]

    print(f"Frequency: {freq/1000:.1f} kHz")
    print(f"\n{'H0 [A/m]':>12} {'R_ac [mOhm]':>14} {'R_ac/R_dc':>12}")
    print("-" * 42)

    for H0 in H0_values:
        Z, info = model_steel.compute_impedance(freq, H0=H0, method='esim')
        print(f"{H0:>12.0f} {Z.real*1000:>14.6f} {info['R_ac_Rdc']:>12.4f}")


def compare_esim_dowell_accuracy():
    """Compare ESIM and Dowell accuracy across frequency range."""
    print("\n" + "=" * 70)
    print("ESIM vs Dowell Accuracy Comparison")
    print("=" * 70)

    sigma = 5.8e7
    half_thickness = 0.001  # 1 mm

    xi_values = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]

    print(f"\n{'xi':>8} {'Dowell':>12} {'ESIM':>12} {'ESIM/Dowell':>12}")
    print("-" * 48)

    for xi in xi_values:
        # Compute frequency for this xi
        omega = 2 * xi**2 / (half_thickness**2 * mu_0 * sigma)
        freq = omega / (2 * np.pi)

        # Dowell
        def dowell_ratio(xi):
            if xi < 0.01:
                return 1.0
            sh2, sn2 = np.sinh(2*xi), np.sin(2*xi)
            ch2, cs2 = np.cosh(2*xi), np.cos(2*xi)
            return xi * (sh2 + sn2) / (ch2 - cs2)

        r_dowell = dowell_ratio(xi)

        # ESIM
        solver = ESIMFiniteSlabSolver(
            half_thickness=half_thickness,
            sigma=sigma,
            frequency=max(freq, 1),
            mu_r=1.0,
            n_nodes=200
        )
        result = solver.solve(1.0)
        Z_s = result['Z']
        r_esim = max(Z_s.real * half_thickness * sigma, 1.0)

        ratio = r_esim / r_dowell if r_dowell > 0 else 0

        print(f"{xi:>8.2f} {r_dowell:>12.4f} {r_esim:>12.4f} {ratio:>12.4f}")

    print("""
    Note: ESIM gives slightly different values than Dowell because:
    1. Dowell is analytical (exact for planar, linear)
    2. ESIM is numerical (finite difference discretization)
    3. The formulas are mathematically equivalent for linear materials

    For nonlinear materials, ESIM is the only option.
    """)


if __name__ == '__main__':
    test_esim_conductor_model()
    compare_esim_dowell_accuracy()

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print("""
    ESIM Conductor Model Usage:

    1. For LINEAR materials:
       - Use Dowell's formula (fast, analytical)
       - ESIM gives same results (slower, numerical)

    2. For NONLINEAR materials (mu(H) saturation):
       - Use ESIM (handles H-dependent permeability)
       - Iterate with surface field H0 until converged

    3. Key formula:
       R_ac = R_dc * R_ratio
       where R_ratio comes from Dowell (linear) or ESIM (nonlinear)

    4. The hybrid approach:
       - Check if material is nonlinear
       - Use appropriate method automatically
    """)
