"""
ESIM Cell Problem Solver using NGSolve

This module implements the 1D Cell Problem solver for the Effective Surface
Impedance Method (ESIM) based on Karl Hollaus's paper.

Reference:
    K. Hollaus, M. Kaltenbacher, J. Schoberl, "A Nonlinear Effective Surface
    Impedance in a Magnetic Scalar Potential Formulation," IEEE Trans. Magnetics,
    2025, DOI: 10.1109/TMAG.2025.3613932

The Cell Problem:
    rho * d^2H/ds^2 + j*omega*mu(|H|)*H = 0    for s in [0, infinity)

    Boundary conditions:
        H(0) = H0         (surface tangential field)
        H(infinity) = 0   (field vanishes at infinity)

From the solution, we compute:
    P'(H0) = (1/2) * integral{ E * J* } ds    (active power loss per unit area)
    Q'(H0) = (omega/2) * integral{ H * B* } ds (reactive power per unit area)
    Z(H0) = 2*(P' + j*Q') / |H0|^2            (effective surface impedance)

Author: Radia Development Team
Date: 2026-01-08
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.constants import mu_0

try:
    from ngsolve import *
    import netgen.meshing as ngm
    NGSOLVE_AVAILABLE = True
except ImportError:
    NGSOLVE_AVAILABLE = False

# Alternative: Pure numpy/scipy implementation for environments without NGSolve
SCIPY_AVAILABLE = True
try:
    from scipy.linalg import solve_banded
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
except ImportError:
    SCIPY_AVAILABLE = False


class BHCurveInterpolator:
    """
    Interpolator for BH curve data.

    Provides mu(|H|) = B(|H|) / |H| with proper handling of H -> 0 limit.
    """

    def __init__(self, bh_curve):
        """
        Initialize BH curve interpolator.

        Parameters:
            bh_curve: list of [H, B] pairs, where H is in A/m and B is in Tesla.
                      Must be sorted by increasing H.
        """
        self.bh_data = np.array(bh_curve)
        self.H_data = self.bh_data[:, 0]
        self.B_data = self.bh_data[:, 1]

        # Ensure data starts from H=0
        if self.H_data[0] != 0:
            self.H_data = np.insert(self.H_data, 0, 0)
            self.B_data = np.insert(self.B_data, 0, 0)

        # Create B(H) interpolator with cubic spline
        self._B_interp = interp1d(
            self.H_data, self.B_data,
            kind='cubic',
            fill_value='extrapolate',
            bounds_error=False
        )

        # Compute initial permeability (slope at H=0)
        if len(self.H_data) >= 2 and self.H_data[1] > 0:
            self.mu_initial = self.B_data[1] / self.H_data[1]
        else:
            self.mu_initial = mu_0 * 1000  # Default: mu_r = 1000

    def B(self, H_abs):
        """
        Get B for given |H|.

        Parameters:
            H_abs: Absolute value of magnetic field [A/m]

        Returns:
            B: Magnetic flux density [T]
        """
        return float(self._B_interp(H_abs))

    def mu(self, H_abs):
        """
        Get permeability mu = B/H for given |H|.

        Parameters:
            H_abs: Absolute value of magnetic field [A/m]

        Returns:
            mu: Permeability [H/m]
        """
        if H_abs < 1e-10:
            return self.mu_initial

        B = self.B(H_abs)
        return B / H_abs

    def mu_r(self, H_abs):
        """
        Get relative permeability mu_r = mu/mu_0.

        Parameters:
            H_abs: Absolute value of magnetic field [A/m]

        Returns:
            mu_r: Relative permeability (dimensionless)
        """
        return self.mu(H_abs) / mu_0

    def differential_mu(self, H_abs):
        """
        Get differential permeability dB/dH.

        Parameters:
            H_abs: Absolute value of magnetic field [A/m]

        Returns:
            mu_diff: Differential permeability [H/m]
        """
        # Use finite difference
        dH = max(H_abs * 0.001, 1.0)
        B1 = self.B(H_abs - dH/2)
        B2 = self.B(H_abs + dH/2)
        return (B2 - B1) / dH


class ESIMCellProblemSolver:
    """
    Solves the 1D Cell Problem for ESIM using finite differences.

    The Cell Problem is a 1D boundary value problem on a semi-infinite domain
    that computes the effective surface impedance Z(H0) for a given surface
    tangential field amplitude H0.

    This implementation uses a finite difference method with geometric mesh
    grading for accuracy near the surface, where field gradients are largest.
    """

    def __init__(self, bh_curve, sigma, frequency,
                 num_skin_depths=10, n_nodes=100):
        """
        Initialize the Cell Problem solver.

        Parameters:
            bh_curve: BH curve data as [[H1, B1], [H2, B2], ...]
                      H in A/m, B in Tesla
            sigma: Electrical conductivity [S/m]
            frequency: Operating frequency [Hz]
            num_skin_depths: Domain depth in skin depths (default: 10)
            n_nodes: Number of mesh nodes (default: 100)
        """
        if not SCIPY_AVAILABLE:
            raise ImportError("Scipy is required for ESIM Cell Problem solver")

        self.bh_interp = BHCurveInterpolator(bh_curve)
        self.sigma = sigma
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency
        self.rho = 1.0 / sigma  # Resistivity [Ohm*m]

        self.num_skin_depths = num_skin_depths
        self.n_nodes = n_nodes

        # Estimate initial skin depth using initial permeability
        self.delta_initial = np.sqrt(2 * self.rho / (self.omega * self.bh_interp.mu_initial))

        # Domain length
        self.L = self.num_skin_depths * self.delta_initial

        # Create mesh with geometric grading
        self._create_mesh()

    def _create_mesh(self):
        """
        Create 1D mesh with geometric grading near the surface (s=0).

        Uses a geometric progression to cluster nodes near s=0 where
        field gradients are largest (within skin depth).
        """
        # Geometric grading: finer near s=0
        # s_i = L * (r^i - 1) / (r^N - 1), where r is the grading ratio

        grading = 1.1  # Grading ratio (>1 means finer near s=0)
        N = self.n_nodes - 1

        if abs(grading - 1.0) < 1e-10:
            # Uniform mesh
            self.mesh_points = np.linspace(0, self.L, self.n_nodes)
        else:
            # Geometric progression
            i = np.arange(self.n_nodes)
            r = grading
            self.mesh_points = self.L * (r ** i - 1) / (r ** N - 1)

        self.n_elements = self.n_nodes - 1

    def _linear_skin_depth(self, mu):
        """Calculate skin depth for given permeability."""
        return np.sqrt(2 * self.rho / (self.omega * mu))

    def solve(self, H0, tol=1e-6, max_iter=50, relaxation=0.5):
        """
        Solve the Cell Problem for given surface field H0.

        Uses Picard iteration for the nonlinear problem with finite difference
        discretization.

        The governing equation:
            rho * d^2H/ds^2 + j*omega*mu(|H|)*H = 0

        Discretized using central differences:
            rho * (H_{i+1} - 2*H_i + H_{i-1}) / h^2 + j*omega*mu_i*H_i = 0

        Parameters:
            H0: Surface tangential field amplitude [A/m] (real, positive)
            tol: Convergence tolerance
            max_iter: Maximum Picard iterations
            relaxation: Under-relaxation parameter (0 < alpha <= 1)

        Returns:
            result: dict with keys:
                'Z': Complex effective surface impedance [Ohm]
                'P_prime': Active power loss per unit area [W/m^2]
                'Q_prime': Reactive power per unit area [var/m^2]
                'H_solution': Array with H(s) distribution
                'converged': True if converged
                'iterations': Number of iterations
        """
        # Initial guess: exponential decay with initial skin depth
        delta = self.delta_initial
        H = H0 * np.exp(-self.mesh_points / delta) * np.exp(-1j * self.mesh_points / delta)

        # Initial permeability distribution (constant)
        mu_dist = np.full(self.n_nodes, self.bh_interp.mu_initial)

        converged = False

        for iteration in range(max_iter):
            # Build and solve the linear system with current mu distribution
            H_new = self._solve_linear_system(H0, mu_dist)

            # Update permeability based on |H| distribution
            mu_new = np.array([self.bh_interp.mu(abs(h)) for h in H_new])

            # Check convergence (relative change in mu)
            rel_change = np.max(np.abs(mu_new - mu_dist)) / np.max(np.abs(mu_dist))

            # Under-relaxation
            mu_dist = (1 - relaxation) * mu_dist + relaxation * mu_new
            H = H_new

            if rel_change < tol:
                converged = True
                break

        # Compute power losses and impedance
        P_prime, Q_prime = self._compute_power_losses(H, mu_dist)

        # Effective surface impedance: Z = 2*(P' + j*Q') / |H0|^2
        Z = 2 * (P_prime + 1j * Q_prime) / (H0 ** 2)

        # Average final permeability (for reporting)
        mu_final = np.mean(mu_dist[:10])  # Average near surface

        return {
            'Z': Z,
            'P_prime': P_prime,
            'Q_prime': Q_prime,
            'H_solution': H,
            'mesh_points': self.mesh_points,
            'converged': converged,
            'iterations': iteration + 1,
            'mu_final': mu_final
        }

    def _solve_linear_system(self, H0, mu_dist):
        """
        Solve the linearized Cell Problem using finite differences.

        The equation:
            rho * d^2H/ds^2 + j*omega*mu*H = 0

        With boundary conditions:
            H(0) = H0
            H(L) = 0 (approximation for H(infinity) = 0)

        Parameters:
            H0: Surface field value
            mu_dist: Permeability distribution at mesh points

        Returns:
            H: Solution array
        """
        n = self.n_nodes
        s = self.mesh_points

        # Build tridiagonal system: A * H = b
        # For non-uniform mesh, use general finite difference formula

        # Coefficient arrays for tridiagonal matrix
        diag_main = np.zeros(n, dtype=complex)
        diag_lower = np.zeros(n - 1, dtype=complex)
        diag_upper = np.zeros(n - 1, dtype=complex)
        rhs = np.zeros(n, dtype=complex)

        # Interior points (i = 1, ..., n-2)
        for i in range(1, n - 1):
            h_minus = s[i] - s[i - 1]  # h_{i-1/2}
            h_plus = s[i + 1] - s[i]   # h_{i+1/2}
            h_avg = 0.5 * (h_minus + h_plus)

            # Second derivative: d^2H/ds^2 ≈ (H_{i+1} - H_i)/h_plus - (H_i - H_{i-1})/h_minus) / h_avg
            coef_lower = self.rho / (h_minus * h_avg)
            coef_upper = self.rho / (h_plus * h_avg)
            coef_main = -coef_lower - coef_upper + 1j * self.omega * mu_dist[i]

            diag_lower[i - 1] = coef_lower
            diag_main[i] = coef_main
            diag_upper[i] = coef_upper
            rhs[i] = 0.0

        # Boundary condition at s=0: H(0) = H0
        diag_main[0] = 1.0
        rhs[0] = H0

        # Boundary condition at s=L: H(L) = 0
        diag_main[n - 1] = 1.0
        rhs[n - 1] = 0.0

        # Solve tridiagonal system using scipy
        # Construct sparse matrix
        from scipy.sparse import diags as sp_diags
        from scipy.sparse.linalg import spsolve as sp_solve

        A = sp_diags(
            [diag_lower, diag_main, diag_upper],
            offsets=[-1, 0, 1],
            format='csr'
        )

        H = sp_solve(A, rhs)

        return H

    def _compute_power_losses(self, H, mu_dist):
        """
        Compute specific power losses from the Cell Problem solution.

        P' = (1/2) * integral{ rho * |dH/ds|^2 } ds
             (from E = -rho * dH/ds, J = -dH/ds, P = (1/2)*Re(E*J*))

        Q' = (omega/2) * integral{ mu * |H|^2 } ds
             (from reactive power in magnetic field)

        Parameters:
            H: Solution array H(s)
            mu_dist: Permeability distribution

        Returns:
            P_prime: Active power loss per unit area [W/m^2]
            Q_prime: Reactive power per unit area [var/m^2]
        """
        s = self.mesh_points
        n = self.n_nodes

        P_prime = 0.0
        Q_prime = 0.0

        # Trapezoidal integration
        for i in range(n - 1):
            ds = s[i + 1] - s[i]

            # dH/ds at midpoint (central difference)
            dHds = (H[i + 1] - H[i]) / ds

            # H and mu at midpoint (average)
            H_mid = 0.5 * (H[i] + H[i + 1])
            mu_mid = 0.5 * (mu_dist[i] + mu_dist[i + 1])

            # Active power: P' = (1/2) * rho * |dH/ds|^2
            P_prime += 0.5 * self.rho * np.abs(dHds) ** 2 * ds

            # Reactive power: Q' = (omega/2) * mu * |H|^2
            Q_prime += 0.5 * self.omega * mu_mid * np.abs(H_mid) ** 2 * ds

        return float(P_prime), float(Q_prime)

    def generate_esi_table(self, H0_values, tol=1e-6, max_iter=50):
        """
        Generate ESI table for a range of H0 values.

        Parameters:
            H0_values: List or array of surface field amplitudes [A/m]
            tol: Convergence tolerance for each solve
            max_iter: Maximum iterations for each solve

        Returns:
            table: numpy array with columns:
                   [H0, Re(Z), Im(Z), P', Q', mu_final, converged]
        """
        table = []

        for H0 in H0_values:
            result = self.solve(H0, tol=tol, max_iter=max_iter)

            row = [
                H0,
                result['Z'].real,
                result['Z'].imag,
                result['P_prime'],
                result['Q_prime'],
                result['mu_final'],
                1.0 if result['converged'] else 0.0
            ]
            table.append(row)

        return np.array(table)

    def get_linear_sibc(self, mu_r=None):
        """
        Get classical local SIBC impedance for comparison.

        For linear materials: Z = (1+j) / (sigma * delta)
        where delta = sqrt(2*rho / (omega*mu))

        Parameters:
            mu_r: Relative permeability (optional, uses initial if not given)

        Returns:
            Z_local: Complex local surface impedance [Ohm]
        """
        if mu_r is None:
            mu = self.bh_interp.mu_initial
        else:
            mu = mu_0 * mu_r

        delta = np.sqrt(2 * self.rho / (self.omega * mu))
        Z_local = (1 + 1j) / (self.sigma * delta)

        return Z_local


class ESITable:
    """
    ESI (Effective Surface Impedance) lookup table with interpolation.

    Stores pre-computed Z(H0) values and provides fast interpolation
    for the 3D solver's fixed-point iteration.
    """

    def __init__(self, table_data=None, bh_curve=None, sigma=None, frequency=None):
        """
        Initialize ESI table.

        Either provide pre-computed table_data, or (bh_curve, sigma, frequency)
        to generate the table.

        Parameters:
            table_data: Pre-computed table as numpy array
            bh_curve: BH curve data (for generation)
            sigma: Conductivity [S/m] (for generation)
            frequency: Frequency [Hz] (for generation)
        """
        if table_data is not None:
            self._load_table(table_data)
        elif bh_curve is not None and sigma is not None and frequency is not None:
            self._generate_table(bh_curve, sigma, frequency)
        else:
            raise ValueError("Provide either table_data or (bh_curve, sigma, frequency)")

    def _load_table(self, table_data):
        """Load pre-computed table."""
        self.table = np.array(table_data)
        self.H0_values = self.table[:, 0]
        self.Z_real = self.table[:, 1]
        self.Z_imag = self.table[:, 2]
        self.P_prime = self.table[:, 3]
        self.Q_prime = self.table[:, 4]

        # Ensure P_prime and Q_prime are positive for log interpolation
        # Use absolute value and store sign for Q_prime (can be negative in some cases)
        P_prime_safe = np.maximum(self.P_prime, 1e-20)
        Q_prime_safe = np.maximum(np.abs(self.Q_prime), 1e-20)

        # Create interpolators (log-H for better accuracy over wide range)
        log_H0 = np.log10(self.H0_values + 1e-20)

        self._Z_real_interp = interp1d(
            log_H0, self.Z_real,
            kind='cubic', fill_value='extrapolate'
        )
        self._Z_imag_interp = interp1d(
            log_H0, self.Z_imag,
            kind='cubic', fill_value='extrapolate'
        )
        # Use log-log interpolation for power (scales as H^2)
        self._P_interp = interp1d(
            log_H0, np.log10(P_prime_safe),
            kind='cubic', fill_value='extrapolate'
        )
        # Use linear interpolation for Q_prime (can have sign changes)
        self._Q_interp = interp1d(
            log_H0, self.Q_prime,
            kind='cubic', fill_value='extrapolate'
        )

    def _generate_table(self, bh_curve, sigma, frequency, n_points=50):
        """Generate ESI table from BH curve."""
        solver = ESIMCellProblemSolver(bh_curve, sigma, frequency)

        # H0 range: from 1 A/m to 100,000 A/m (log-spaced)
        H0_values = np.logspace(0, 5, n_points)

        table = solver.generate_esi_table(H0_values)
        self._load_table(table)

        # Store parameters
        self.bh_curve = bh_curve
        self.sigma = sigma
        self.frequency = frequency

    def get_impedance(self, H0):
        """
        Get effective surface impedance for given H0.

        Parameters:
            H0: Surface tangential field amplitude [A/m]

        Returns:
            Z: Complex surface impedance [Ohm]
        """
        log_H0 = np.log10(max(H0, 1e-20))
        Z_real = float(self._Z_real_interp(log_H0))
        Z_imag = float(self._Z_imag_interp(log_H0))
        return Z_real + 1j * Z_imag

    def get_power_loss(self, H0):
        """
        Get specific power losses for given H0.

        Parameters:
            H0: Surface tangential field amplitude [A/m]

        Returns:
            P_prime: Active power loss per unit area [W/m^2]
            Q_prime: Reactive power per unit area [var/m^2]
        """
        log_H0 = np.log10(max(H0, 1e-20))
        P_prime = 10 ** float(self._P_interp(log_H0))
        Q_prime = float(self._Q_interp(log_H0))  # Linear interpolation, not log
        return P_prime, Q_prime

    def save(self, filename):
        """Save ESI table to file."""
        header = "# ESIM ESI Table\n"
        header += f"# Sigma: {getattr(self, 'sigma', 'N/A')} S/m\n"
        header += f"# Frequency: {getattr(self, 'frequency', 'N/A')} Hz\n"
        header += "# H0 [A/m]  Re(Z) [Ohm]  Im(Z) [Ohm]  P' [W/m^2]  Q' [var/m^2]\n"

        np.savetxt(filename, self.table[:, :5], header=header,
                   fmt=['%.6e', '%.6e', '%.6e', '%.6e', '%.6e'])

    @classmethod
    def load(cls, filename):
        """Load ESI table from file."""
        data = np.loadtxt(filename)
        return cls(table_data=data)


# Convenience function for Python API
def generate_esi_table_from_bh_curve(bh_curve, sigma, frequency, n_points=50):
    """
    Generate ESI table from BH curve data.

    This is the main entry point for ESIM table generation.

    Parameters:
        bh_curve: BH curve as [[H1, B1], [H2, B2], ...]
                  H in A/m, B in Tesla
        sigma: Conductivity [S/m]
        frequency: Operating frequency [Hz]
        n_points: Number of H0 sampling points

    Returns:
        esi_table: ESITable object
    """
    # Create solver
    solver = ESIMCellProblemSolver(bh_curve, sigma, frequency)

    # Generate H0 values (log-spaced from 1 to 100,000 A/m)
    H0_values = np.logspace(0, 5, n_points)

    # Generate table
    table = solver.generate_esi_table(H0_values)

    return ESITable(table_data=table)


# Example usage and test
if __name__ == "__main__":
    # Test with typical steel BH curve
    bh_curve_steel = [
        [0, 0],
        [100, 0.2],
        [250, 0.5],
        [500, 0.9],
        [1000, 1.3],
        [2500, 1.6],
        [5000, 1.8],
        [10000, 1.95],
        [50000, 2.1],
    ]

    sigma_steel = 2e6  # S/m (hot steel)
    freq = 50000  # 50 kHz

    print("ESIM Cell Problem Solver Test")
    print("=" * 50)
    print(f"Material: Steel (sigma = {sigma_steel/1e6:.1f} MS/m)")
    print(f"Frequency: {freq/1000:.1f} kHz")
    print()

    if NGSOLVE_AVAILABLE:
        # Create solver
        solver = ESIMCellProblemSolver(bh_curve_steel, sigma_steel, freq)

        print(f"Initial skin depth: {solver.delta_initial*1e3:.3f} mm")
        print(f"Domain length: {solver.L*1e3:.1f} mm")
        print(f"Number of elements: {solver.n_elements}")
        print()

        # Solve for different H0 values
        H0_test = [10, 100, 1000, 10000]

        print("Cell Problem Solutions:")
        print("-" * 70)
        print(f"{'H0 [A/m]':>12} {'Re(Z) [mOhm]':>14} {'Im(Z) [mOhm]':>14} {'P [kW/m^2]':>14}")
        print("-" * 70)

        for H0 in H0_test:
            result = solver.solve(H0)
            Z = result['Z']
            P = result['P_prime']

            print(f"{H0:>12.0f} {Z.real*1e3:>14.4f} {Z.imag*1e3:>14.4f} {P/1e3:>14.2f}")

        print()

        # Compare with linear SIBC
        Z_linear = solver.get_linear_sibc()
        print(f"Linear SIBC (mu_r_initial): Z = {Z_linear.real*1e3:.4f} + j{Z_linear.imag*1e3:.4f} mOhm")
    else:
        print("NGSolve not available. Skipping test.")
