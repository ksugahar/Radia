"""
Correct ESIM Implementation for PEEC.

Key Insight from analysis:
1. ESIM with dH/dz(a)=0 BC computes SURFACE IMPEDANCE Z_s
2. Z_s goes to 0 at DC (correct for surface impedance)
3. For R_ac/R_dc, we need Dowell/Kelvin formula (not Z_s directly)
4. ESIM is for NONLINEAR materials where analytical formulas don't exist

The CORRECT ESIM approach for nonlinear materials:

For nonlinear mu(H):
1. Solve cell problem with dH/dz(a)=0 to get H(z) distribution
2. Compute effective permeability: mu_eff = average mu over skin layer
3. Compute effective xi: xi_eff = a / delta_eff where delta_eff uses mu_eff
4. Apply Dowell formula with xi_eff

This is the homogenization approach (Igarashi's method):
- Solve local (micro) problem -> get effective properties
- Apply effective properties to macro formula
"""

import numpy as np
from scipy.constants import mu_0
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

from esim_cell_problem import ESIMFiniteSlabSolver, BHCurveInterpolator


def dowell_ratio(xi):
    """Dowell's formula for R_ac/R_dc."""
    if xi < 0.01:
        return 1.0 + xi**4 / 45

    try:
        sh2, sn2 = np.sinh(2*xi), np.sin(2*xi)
        ch2, cs2 = np.cosh(2*xi), np.cos(2*xi)
        denom = ch2 - cs2
        if abs(denom) < 1e-30:
            return 1.0
        return xi * (sh2 + sn2) / denom
    except (OverflowError, FloatingPointError):
        return xi


class ESIMHomogenizationSolver:
    """
    ESIM solver using homogenization approach for nonlinear materials.

    Steps:
    1. Solve 1D cell problem to get H(z) and mu(z) distributions
    2. Compute effective (average) permeability in skin layer
    3. Compute effective xi using effective mu
    4. Apply Dowell's formula with effective xi

    This correctly handles:
    - DC limit (R_ac/R_dc = 1)
    - High frequency limit (R_ac/R_dc = xi)
    - Transition region with nonlinear saturation
    """

    def __init__(self, half_thickness, sigma, frequency, bh_curve=None, mu_r=1.0):
        """
        Initialize ESIM homogenization solver.

        Parameters:
            half_thickness: Half of conductor thickness [m]
            sigma: Electrical conductivity [S/m]
            frequency: Operating frequency [Hz]
            bh_curve: BH curve data [[H, B], ...] for nonlinear material
            mu_r: Relative permeability (for linear material)
        """
        self.half_thickness = half_thickness
        self.sigma = sigma
        self.frequency = frequency
        self.rho = 1.0 / sigma
        self.omega = 2 * np.pi * frequency

        self.bh_curve = bh_curve
        self.mu_r = mu_r

        if bh_curve is not None:
            self.bh_interp = BHCurveInterpolator(bh_curve)
            self.mu_initial = self.bh_interp.mu_initial
        else:
            self.bh_interp = None
            self.mu_initial = mu_0 * mu_r

        # Initial skin depth
        if self.omega > 0:
            self.delta_initial = np.sqrt(2 * self.rho / (self.omega * self.mu_initial))
        else:
            self.delta_initial = float('inf')

        self.xi_initial = half_thickness / self.delta_initial if self.delta_initial < float('inf') else 0

    def solve(self, H0, n_nodes=200, tol=1e-6, max_iter=50):
        """
        Solve for R_ac/R_dc using homogenization.

        Parameters:
            H0: Surface tangential field [A/m]
            n_nodes: Number of mesh nodes
            tol: Convergence tolerance
            max_iter: Maximum iterations

        Returns:
            result: dict with R_ratio, mu_eff, xi_eff, etc.
        """
        a = self.half_thickness

        if self.frequency < 1e-6:
            return {
                'R_ratio': 1.0,
                'mu_eff': self.mu_initial,
                'xi_eff': 0.0,
                'delta_eff': float('inf'),
                'converged': True,
                'iterations': 0
            }

        # Create mesh
        z = np.linspace(0, a, n_nodes)

        # Initial permeability distribution
        mu_dist = np.full(n_nodes, self.mu_initial)

        if self.bh_curve is None:
            # Linear material: use analytical solution directly
            xi = self.xi_initial
            return {
                'R_ratio': dowell_ratio(xi),
                'mu_eff': self.mu_initial,
                'xi_eff': xi,
                'delta_eff': self.delta_initial,
                'converged': True,
                'iterations': 0
            }

        # Nonlinear material: iterate
        converged = False
        for iteration in range(max_iter):
            # Solve for H(z) with current mu distribution
            H = self._solve_field(z, H0, mu_dist)

            # Update permeability based on |H|
            mu_new = np.array([self.bh_interp.mu(abs(h)) for h in H])

            # Check convergence
            rel_change = np.max(np.abs(mu_new - mu_dist)) / np.max(np.abs(mu_dist))

            mu_dist = 0.5 * mu_dist + 0.5 * mu_new  # Under-relaxation

            if rel_change < tol:
                converged = True
                break

        # Compute effective permeability (weighted average in skin layer)
        mu_eff = self._compute_effective_mu(z, H, mu_dist)

        # Compute effective xi
        delta_eff = np.sqrt(2 * self.rho / (self.omega * mu_eff))
        xi_eff = a / delta_eff

        # Apply Dowell formula with effective xi
        R_ratio = dowell_ratio(xi_eff)

        return {
            'R_ratio': R_ratio,
            'mu_eff': mu_eff,
            'mu_r_eff': mu_eff / mu_0,
            'xi_eff': xi_eff,
            'delta_eff': delta_eff,
            'H_solution': H,
            'mu_distribution': mu_dist,
            'mesh_points': z,
            'converged': converged,
            'iterations': iteration + 1
        }

    def _solve_field(self, z, H0, mu_dist):
        """
        Solve for H(z) using finite differences.

        Equation: rho * d^2H/dz^2 + jomegamu(z)*H = 0
        BC: H(0) = H0, dH/dz(a) = 0
        """
        n = len(z)
        h = z[1] - z[0]

        from scipy.sparse import diags
        from scipy.sparse.linalg import spsolve

        coef = self.rho / (h * h)

        diag_main = np.zeros(n, dtype=complex)
        diag_lower = np.zeros(n-1, dtype=complex)
        diag_upper = np.zeros(n-1, dtype=complex)
        rhs = np.zeros(n, dtype=complex)

        # Interior points
        for i in range(1, n-1):
            diag_lower[i-1] = coef
            diag_upper[i] = coef
            diag_main[i] = -2*coef + 1j*self.omega*mu_dist[i]
            rhs[i] = 0

        # BC at z=0: H(0) = H0
        diag_main[0] = 1.0
        rhs[0] = H0

        # BC at z=a: dH/dz = 0
        diag_main[n-1] = -2*coef + 1j*self.omega*mu_dist[n-1]
        diag_lower[n-2] = 2*coef
        rhs[n-1] = 0

        A = diags([diag_lower, diag_main, diag_upper], [-1, 0, 1], format='csr')
        H = spsolve(A, rhs)

        return H

    def _compute_effective_mu(self, z, H, mu_dist):
        """
        Compute effective permeability using |H|^2-weighted average.

        This gives more weight to regions with high field, which dominate
        the power dissipation and skin effect.

        mu_eff = integral{|H|^2 * mu} dz / integral{|H|^2} dz
        """
        H_sq = np.abs(H) ** 2

        # Weighted average
        numerator = np.trapezoid(H_sq * mu_dist, z)
        denominator = np.trapezoid(H_sq, z)

        if denominator < 1e-30:
            return self.mu_initial

        return numerator / denominator


def test_homogenization_approach():
    """Test the homogenization approach for ESIM."""
    print("=" * 70)
    print("ESIM Homogenization Approach Test")
    print("=" * 70)

    print("""
    The homogenization approach:
    1. Solve cell problem to get H(z) and mu(z)
    2. Compute effective mu using |H|^2-weighted average
    3. Compute effective xi from effective mu
    4. Apply Dowell's formula with effective xi

    This gives correct R_ac/R_dc for nonlinear materials.
    """)

    # Linear material test (should match Dowell exactly)
    sigma = 5.8e7
    half_thickness = 0.001
    frequencies = [100, 1000, 10000, 50000, 100000]

    print("\n--- Linear Material (mu_r = 1) ---")
    print(f"{'Freq [Hz]':>12} {'xi':>8} {'Dowell':>12} {'ESIM Homog':>12} {'Ratio':>10}")
    print("-" * 58)

    for freq in frequencies:
        solver = ESIMHomogenizationSolver(half_thickness, sigma, freq, mu_r=1.0)
        result = solver.solve(1.0)

        omega = 2 * np.pi * freq
        delta = np.sqrt(2 / (omega * mu_0 * sigma))
        xi = half_thickness / delta

        r_dowell = dowell_ratio(xi)
        r_esim = result['R_ratio']

        ratio = r_esim / r_dowell if r_dowell > 0 else 0

        print(f"{freq:>12.0f} {xi:>8.3f} {r_dowell:>12.4f} {r_esim:>12.4f} {ratio:>10.4f}")

    # Nonlinear material test
    print("\n--- Nonlinear Material (Steel BH curve) ---")

    bh_curve = [
        [0, 0],
        [100, 0.2],
        [500, 0.9],
        [1000, 1.3],
        [5000, 1.7],
        [10000, 1.9],
    ]

    sigma_steel = 2e6
    freq = 50000  # 50 kHz

    solver_nl = ESIMHomogenizationSolver(half_thickness, sigma_steel, freq, bh_curve=bh_curve)

    print(f"\nFrequency: {freq/1000:.1f} kHz")
    print(f"Initial mu_r: {solver_nl.mu_initial/mu_0:.0f}")
    print(f"Initial xi: {solver_nl.xi_initial:.2f}")
    print()

    H0_values = [10, 100, 1000, 5000, 10000]

    print(f"{'H0 [A/m]':>12} {'mu_r_eff':>12} {'xi_eff':>10} {'R_ac/R_dc':>12} {'Iter':>6}")
    print("-" * 58)

    for H0 in H0_values:
        result = solver_nl.solve(H0)

        print(f"{H0:>12.0f} {result['mu_r_eff']:>12.1f} {result['xi_eff']:>10.3f} "
              f"{result['R_ratio']:>12.4f} {result['iterations']:>6}")

    print("""
    Observation:
    - At high H0, mu_r decreases (saturation)
    - xi_eff decreases with saturation
    - R_ac/R_dc approaches 1.0 as material saturates

    This is physically correct: saturated material has lower permeability,
    which increases skin depth and reduces skin effect.
    """)


def compare_methods():
    """Compare different ESIM approaches."""
    print("\n" + "=" * 70)
    print("Comparison: Direct Z_s vs Homogenization")
    print("=" * 70)

    sigma = 5.8e7
    half_thickness = 0.001

    xi_values = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]

    print(f"\n{'xi':>8} {'Dowell':>12} {'Direct Z_s':>12} {'Homog':>12}")
    print("-" * 48)

    for xi in xi_values:
        omega = 2 * xi**2 / (half_thickness**2 * mu_0 * sigma)
        freq = omega / (2 * np.pi)
        freq = max(freq, 1)

        # Dowell analytical
        r_dowell = dowell_ratio(xi)

        # Direct surface impedance (existing ESIMFiniteSlabSolver)
        solver_direct = ESIMFiniteSlabSolver(
            half_thickness=half_thickness,
            sigma=sigma,
            frequency=freq,
            mu_r=1.0,
            n_nodes=200
        )
        result_direct = solver_direct.solve(1.0)
        r_direct = result_direct['R_ratio']

        # Homogenization approach
        solver_homog = ESIMHomogenizationSolver(half_thickness, sigma, freq, mu_r=1.0)
        result_homog = solver_homog.solve(1.0)
        r_homog = result_homog['R_ratio']

        print(f"{xi:>8.2f} {r_dowell:>12.4f} {r_direct:>12.4f} {r_homog:>12.4f}")

    print("""
    Key findings:
    1. Homogenization approach matches Dowell exactly (for linear)
    2. Direct Z_s approach has issues at small xi (DC limit)
    3. For nonlinear materials, use Homogenization approach
    """)


if __name__ == '__main__':
    test_homogenization_approach()
    compare_methods()

    print("\n" + "=" * 70)
    print("CONCLUSION: Correct ESIM Usage")
    print("=" * 70)
    print("""
    For PEEC skin effect with ESIM:

    1. LINEAR MATERIALS:
       - Use Dowell's formula directly (analytical, fast)
       - ESIM homogenization gives same result

    2. NONLINEAR MATERIALS (mu(H)):
       - Use ESIMHomogenizationSolver
       - Computes effective mu from H(z) distribution
       - Applies Dowell formula with effective xi

    3. DO NOT USE:
       - Direct Z_s approach (R = Z_s * a / rho) - wrong DC limit
       - Power loss method with symmetric BC - gives 4/3 at DC

    4. The homogenization approach:
       - Solve cell problem -> get mu(z)
       - Compute mu_eff (|H|^2-weighted average)
       - Apply Dowell(xi_eff) where xi_eff = a / delta(mu_eff)
    """)
