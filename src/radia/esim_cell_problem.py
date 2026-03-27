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

        # Create B(H) interpolator
        # Use cubic if enough points, otherwise linear
        interp_kind = 'cubic' if len(self.H_data) >= 4 else 'linear'
        self._B_interp = interp1d(
            self.H_data, self.B_data,
            kind=interp_kind,
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


class ComplexPermeabilityInterpolator:
    """
    Interpolator for complex permeability mu = mu' - j*mu".

    This class handles materials with magnetic losses where the permeability
    is complex: mu = mu' - j*mu"
    - mu' (real part): Energy storage (reactive)
    - mu" (imaginary part): Energy loss (hysteresis, eddy current in grains, etc.)

    The loss tangent is defined as: tan(delta_m) = mu" / mu'

    Complex permeability affects:
    1. Skin depth: delta = sqrt(2*rho / (omega * |mu|))
    2. Surface impedance: Z_s = sqrt(j*omega*mu / sigma)
    3. Power loss: Both ohmic (J^2/sigma) and magnetic (omega*mu"*H^2)
    """

    def __init__(self, mu_data, frequency=None):
        """
        Initialize complex permeability interpolator.

        Parameters:
            mu_data: Permeability data in one of these formats:
                Format 1: [[H, mu'_r, mu"_r], ...] - H-dependent complex permeability
                Format 2: [[f, mu'_r, mu"_r], ...] - Frequency-dependent (for linear materials)
                Format 3: (mu'_r, mu"_r) - Constant complex permeability
                Format 4: dict with 'mu_prime' and 'mu_double_prime' arrays

            frequency: Operating frequency [Hz] (required for Format 2)

        Units:
            - H in A/m
            - f in Hz
            - mu'_r, mu"_r are relative (dimensionless)
        """
        self.frequency = frequency

        if isinstance(mu_data, tuple) and len(mu_data) == 2:
            # Format 3: Constant permeability (mu'_r, mu"_r)
            self.mode = 'constant'
            self.mu_prime_r = float(mu_data[0])
            self.mu_double_prime_r = float(mu_data[1])
            self.mu_initial = mu_0 * complex(self.mu_prime_r, -self.mu_double_prime_r)

        elif isinstance(mu_data, dict):
            # Format 4: Dictionary with arrays
            self.mode = 'H_dependent'
            self._setup_H_dependent(
                mu_data.get('H', np.array([0, 1e6])),
                mu_data.get('mu_prime', np.array([1000, 100])),
                mu_data.get('mu_double_prime', np.array([100, 10]))
            )

        elif isinstance(mu_data, (list, np.ndarray)):
            mu_array = np.array(mu_data)

            if mu_array.shape[1] == 3:
                # Format 1 or 2: [x, mu'_r, mu"_r]
                # Assume H-dependent if first column spans typical H range
                x_data = mu_array[:, 0]
                if x_data[-1] > 1e4:  # Likely H values (A/m)
                    self.mode = 'H_dependent'
                    self._setup_H_dependent(x_data, mu_array[:, 1], mu_array[:, 2])
                else:  # Likely frequency values (Hz)
                    self.mode = 'f_dependent'
                    self._setup_f_dependent(x_data, mu_array[:, 1], mu_array[:, 2])
            else:
                raise ValueError("mu_data must have 3 columns: [x, mu'_r, mu\"_r]")
        else:
            raise ValueError("Invalid mu_data format")

    def _setup_H_dependent(self, H_data, mu_prime_r, mu_double_prime_r):
        """Setup H-dependent permeability interpolation."""
        self.H_data = np.array(H_data)
        self.mu_prime_r_data = np.array(mu_prime_r)
        self.mu_double_prime_r_data = np.array(mu_double_prime_r)

        # Ensure data starts from H=0
        if self.H_data[0] != 0:
            self.H_data = np.insert(self.H_data, 0, 0)
            self.mu_prime_r_data = np.insert(self.mu_prime_r_data, 0, self.mu_prime_r_data[0])
            self.mu_double_prime_r_data = np.insert(self.mu_double_prime_r_data, 0, self.mu_double_prime_r_data[0])

        # Create interpolators
        self._mu_prime_interp = interp1d(
            self.H_data, self.mu_prime_r_data,
            kind='linear', fill_value='extrapolate', bounds_error=False
        )
        self._mu_double_prime_interp = interp1d(
            self.H_data, self.mu_double_prime_r_data,
            kind='linear', fill_value='extrapolate', bounds_error=False
        )

        # Initial permeability at H=0
        self.mu_initial = mu_0 * complex(self.mu_prime_r_data[0], -self.mu_double_prime_r_data[0])

    def _setup_f_dependent(self, f_data, mu_prime_r, mu_double_prime_r):
        """Setup frequency-dependent permeability interpolation."""
        self.f_data = np.array(f_data)
        self.mu_prime_r_f_data = np.array(mu_prime_r)
        self.mu_double_prime_r_f_data = np.array(mu_double_prime_r)

        # Use log-frequency interpolation
        self._mu_prime_f_interp = interp1d(
            np.log10(self.f_data + 1), self.mu_prime_r_f_data,
            kind='linear', fill_value='extrapolate', bounds_error=False
        )
        self._mu_double_prime_f_interp = interp1d(
            np.log10(self.f_data + 1), self.mu_double_prime_r_f_data,
            kind='linear', fill_value='extrapolate', bounds_error=False
        )

        # Initial permeability at given frequency
        if self.frequency is not None:
            mu_p = float(self._mu_prime_f_interp(np.log10(self.frequency + 1)))
            mu_pp = float(self._mu_double_prime_f_interp(np.log10(self.frequency + 1)))
            self.mu_initial = mu_0 * complex(mu_p, -mu_pp)
        else:
            self.mu_initial = mu_0 * complex(self.mu_prime_r_f_data[0], -self.mu_double_prime_r_f_data[0])

    def mu(self, H_abs):
        """
        Get complex permeability mu = mu' - j*mu" for given |H|.

        Parameters:
            H_abs: Absolute value of magnetic field [A/m]

        Returns:
            mu: Complex permeability [H/m]
        """
        if self.mode == 'constant':
            return mu_0 * complex(self.mu_prime_r, -self.mu_double_prime_r)

        elif self.mode == 'H_dependent':
            mu_p = float(self._mu_prime_interp(H_abs))
            mu_pp = float(self._mu_double_prime_interp(H_abs))
            return mu_0 * complex(mu_p, -mu_pp)

        elif self.mode == 'f_dependent':
            if self.frequency is None:
                raise ValueError("Frequency must be set for f-dependent permeability")
            mu_p = float(self._mu_prime_f_interp(np.log10(self.frequency + 1)))
            mu_pp = float(self._mu_double_prime_f_interp(np.log10(self.frequency + 1)))
            return mu_0 * complex(mu_p, -mu_pp)

    def mu_r(self, H_abs):
        """
        Get complex relative permeability mu_r = mu/mu_0.

        Parameters:
            H_abs: Absolute value of magnetic field [A/m]

        Returns:
            mu_r: Complex relative permeability (dimensionless)
        """
        return self.mu(H_abs) / mu_0

    def mu_prime(self, H_abs):
        """Get real part of permeability mu' [H/m]."""
        return self.mu(H_abs).real

    def mu_double_prime(self, H_abs):
        """Get imaginary part magnitude |mu"| [H/m] (returned as positive value)."""
        return abs(self.mu(H_abs).imag)

    def loss_tangent(self, H_abs):
        """
        Get magnetic loss tangent tan(delta_m) = mu" / mu'.

        Parameters:
            H_abs: Absolute value of magnetic field [A/m]

        Returns:
            tan_delta: Magnetic loss tangent (dimensionless)
        """
        mu_complex = self.mu(H_abs)
        mu_p = mu_complex.real
        mu_pp = abs(mu_complex.imag)

        if abs(mu_p) < 1e-20:
            return 0.0

        return mu_pp / mu_p

    def set_frequency(self, frequency):
        """Set operating frequency [Hz]."""
        self.frequency = frequency
        if self.mode == 'f_dependent':
            mu_p = float(self._mu_prime_f_interp(np.log10(frequency + 1)))
            mu_pp = float(self._mu_double_prime_f_interp(np.log10(frequency + 1)))
            self.mu_initial = mu_0 * complex(mu_p, -mu_pp)


class ESIMFiniteSlabSolver:
    """
    Solves the 1D Cell Problem for a FINITE THICKNESS slab conductor.

    This solver handles the skin effect in conductors with finite thickness,
    which is important for accurate R_ac/R_dc calculation in the transition
    region (0.5 < xi < 3) where xi = half_thickness / skin_depth.

    The geometry is a slab of thickness 2*a with:
    - Surface at z = 0 with tangential H = H0
    - Center at z = a with dH/dz = 0 (symmetry)

    Governing equation:
        rho * d^2H/dz^2 + j*omega*mu(|H|)*H = 0    for z in [0, a]

    Boundary conditions:
        H(0) = H0           (surface tangential field)
        dH/dz(a) = 0        (symmetry at center)

    For linear materials, the analytical solution is:
        H(z) = H0 * cosh(gamma*(a-z)) / cosh(gamma*a)
        where gamma = (1+j)/delta

    This gives Dowell's formula for R_ac/R_dc.
    """

    def __init__(self, half_thickness, bh_curve=None, sigma=None, frequency=None,
                 n_nodes=100, complex_mu=None, mu_r=None, geometry='slab'):
        """
        Initialize the Finite Slab/Cylinder Cell Problem solver.

        Parameters:
            half_thickness: Half of conductor thickness [m] (slab) or radius [m] (cylinder)
            bh_curve: BH curve data as [[H1, B1], [H2, B2], ...]
            sigma: Electrical conductivity [S/m]
            frequency: Operating frequency [Hz]
            n_nodes: Number of mesh nodes (default: 100)
            complex_mu: Complex permeability data (alternative to bh_curve)
            mu_r: Relative permeability for linear material (alternative to bh_curve)
            geometry: 'slab' (default) or 'cylinder' (r-weighted for cylindrical workpiece)
        """
        if not SCIPY_AVAILABLE:
            raise ImportError("Scipy is required for ESIM Cell Problem solver")

        if sigma is None or frequency is None:
            raise ValueError("sigma and frequency are required parameters")

        self.half_thickness = half_thickness  # a (slab) or R (cylinder)
        self.sigma = sigma
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency
        self.rho = 1.0 / sigma
        self.geometry = geometry  # 'slab' or 'cylinder'

        self.n_nodes = n_nodes

        # Setup permeability model
        self.use_complex_mu = complex_mu is not None
        self.linear_mu_r = mu_r

        if self.use_complex_mu:
            self.mu_interp = ComplexPermeabilityInterpolator(complex_mu, frequency)
            self.bh_interp = None
            mu_initial = self.mu_interp.mu_initial
        elif mu_r is not None:
            # Linear material with constant mu_r
            self.mu_interp = None
            self.bh_interp = None
            mu_initial = mu_0 * mu_r
        elif bh_curve is not None:
            self.bh_interp = BHCurveInterpolator(bh_curve)
            self.mu_interp = None
            mu_initial = self.bh_interp.mu_initial
        else:
            raise ValueError("Provide bh_curve, complex_mu, or mu_r")

        self.mu_initial = mu_initial
        mu_abs = abs(mu_initial)

        # Skin depth
        if self.omega > 0:
            self.delta = np.sqrt(2 * self.rho / (self.omega * mu_abs))
        else:
            self.delta = float('inf')

        # xi parameter
        self.xi = half_thickness / self.delta if self.delta < float('inf') else 0

        # Create uniform mesh from surface (z=0) to center (z=a)
        self.mesh_points = np.linspace(0, half_thickness, n_nodes)
        self.n_elements = n_nodes - 1

    def _get_mu(self, H_abs):
        """Get permeability for given |H|."""
        if self.use_complex_mu:
            return self.mu_interp.mu(H_abs)
        elif self.linear_mu_r is not None:
            return mu_0 * self.linear_mu_r
        else:
            return self.bh_interp.mu(H_abs)

    def solve(self, H0, tol=1e-6, max_iter=50, relaxation=0.5):
        """
        Solve the finite slab Cell Problem for given surface field H0.

        Parameters:
            H0: Surface tangential field amplitude [A/m]
            tol: Convergence tolerance
            max_iter: Maximum Picard iterations
            relaxation: Under-relaxation parameter

        Returns:
            result: dict with solution data
        """
        n = self.n_nodes
        z = self.mesh_points
        a = self.half_thickness

        # Initial guess using analytical solution for linear material
        if self.delta < float('inf') and self.delta > 0:
            gamma = complex(1, 1) / self.delta
            if self.geometry == 'cylinder':
                # Cylinder: H(r) = H0 * I0(gamma*r) / I0(gamma*R)
                from scipy.special import iv as bessel_iv
                try:
                    H = H0 * bessel_iv(0, gamma * z) / bessel_iv(0, gamma * a)
                except (OverflowError, FloatingPointError):
                    H = H0 * np.exp(-(a - z) / self.delta)
            else:
                # Slab: H(z) = H0 * cosh(gamma*(a-z)) / cosh(gamma*a)
                try:
                    H = H0 * np.cosh(gamma * (a - z)) / np.cosh(gamma * a)
                except (OverflowError, FloatingPointError):
                    H = H0 * np.exp(-z / self.delta) * np.exp(-1j * z / self.delta)
        else:
            H = np.full(n, H0, dtype=complex)  # DC: uniform current

        # Initial permeability distribution
        mu_dist = np.full(n, self.mu_initial, dtype=complex)

        converged = False

        for iteration in range(max_iter):
            # Build and solve linear system
            H_new = self._solve_linear_system(H0, mu_dist)

            # Update permeability
            mu_new = np.array([self._get_mu(abs(h)) for h in H_new], dtype=complex)

            # Check convergence
            rel_change = np.max(np.abs(mu_new - mu_dist)) / np.max(np.abs(mu_dist))

            # Under-relaxation
            mu_dist = (1 - relaxation) * mu_dist + relaxation * mu_new
            H = H_new

            if rel_change < tol:
                converged = True
                break

        # Compute resistance ratio R_ac/R_dc
        R_ratio = self._compute_resistance_ratio(H, mu_dist, H0)

        # Compute power losses
        P_prime, Q_prime, P_magnetic = self._compute_power_losses(H, mu_dist)

        # Effective surface impedance
        Z = 2 * (P_prime + 1j * Q_prime) / (H0 ** 2)

        return {
            'R_ratio': R_ratio,  # R_ac/R_dc
            'Z': Z,
            'P_prime': P_prime,
            'Q_prime': Q_prime,
            'P_magnetic': P_magnetic,
            'H_solution': H,
            'mesh_points': self.mesh_points,
            'converged': converged,
            'iterations': iteration + 1,
            'xi': self.xi,
            'delta': self.delta
        }

    def _solve_linear_system(self, H0, mu_dist):
        """
        Solve the linearized Cell Problem.

        Slab (geometry='slab'):
            Domain: z from 0 (surface) to a (center)
            PDE: rho * d^2H/dz^2 + jw*mu*H = 0
            BC: H(0) = H0, dH/dz(a) = 0
            Solution: H(z) = H0 * cosh(gamma*(a-z)) / cosh(gamma*a)

        Cylinder (geometry='cylinder'):
            Domain: r from 0 (center) to R (surface)
            PDE: rho/r * d/dr(r * dH/dr) + jw*mu*H = 0
            BC: dH/dr(0) = 0 (symmetry), H(R) = H0
            Solution: H(r) = H0 * I0(gamma*r) / I0(gamma*R)
        """
        n = self.n_nodes
        z = self.mesh_points  # z for slab, r for cylinder
        h = z[1] - z[0]
        a = self.half_thickness  # half-thickness (slab) or radius (cylinder)

        diag_main = np.zeros(n, dtype=complex)
        diag_lower = np.zeros(n - 1, dtype=complex)
        diag_upper = np.zeros(n - 1, dtype=complex)
        rhs = np.zeros(n, dtype=complex)

        if self.geometry == 'cylinder':
            # Cylindrical: (1/r)*d/dr[r*dH/dr] + jw*mu*sigma*H = 0
            # FD on uniform mesh r_i = i*h:
            #   [r_{i+1/2}*(H_{i+1}-H_i) - r_{i-1/2}*(H_i-H_{i-1})] / (r_i*h^2)
            #   + jw*mu*sigma*H_i = 0
            # Multiply by rho: rho/r_i * [...] + jw*mu*H_i = 0

            for i in range(1, n - 1):
                r_i = z[i]
                r_ip = z[i] + h / 2  # r_{i+1/2}
                r_im = z[i] - h / 2  # r_{i-1/2}
                coef_p = self.rho * r_ip / (r_i * h * h)
                coef_m = self.rho * r_im / (r_i * h * h)
                diag_lower[i - 1] = coef_m
                diag_upper[i] = coef_p
                diag_main[i] = -(coef_p + coef_m) + 1j * self.omega * mu_dist[i]
                rhs[i] = 0.0

            # BC at r=0: dH/dr = 0 (symmetry)
            # L'Hopital: (1/r)*d/dr(r*dH/dr)|_{r=0} = 2*d^2H/dr^2
            # With symmetry H[-1]=H[1]: d^2H/dr^2 = 2*(H[1]-H[0])/h^2
            # So equation: rho*4*(H[1]-H[0])/h^2 + jw*mu*H[0] = 0
            coef_0 = 4 * self.rho / (h * h)
            diag_main[0] = -coef_0 + 1j * self.omega * mu_dist[0]
            diag_upper[0] = coef_0
            rhs[0] = 0.0

            # BC at r=R: H(R) = H0 (Dirichlet at surface)
            diag_main[n - 1] = 1.0
            if n >= 2:
                diag_lower[n - 2] = 0.0
            rhs[n - 1] = H0

        else:
            # Slab: rho * d^2H/dz^2 + jw*mu*H = 0
            coef = self.rho / (h * h)

            for i in range(1, n - 1):
                diag_lower[i - 1] = coef
                diag_upper[i] = coef
                diag_main[i] = -2 * coef + 1j * self.omega * mu_dist[i]
                rhs[i] = 0.0

            # BC at z=0: H(0) = H0 (Dirichlet at surface)
            diag_main[0] = 1.0
            rhs[0] = H0

            # BC at z=a: dH/dz = 0 (Neumann - symmetry)
            diag_main[n - 1] = -2 * coef + 1j * self.omega * mu_dist[n - 1]
            diag_lower[n - 2] = 2 * coef
            rhs[n - 1] = 0.0

        # Solve tridiagonal system
        from scipy.sparse import diags as sp_diags
        from scipy.sparse.linalg import spsolve as sp_solve

        A = sp_diags(
            [diag_lower, diag_main, diag_upper],
            offsets=[-1, 0, 1],
            format='csr'
        )

        H = sp_solve(A, rhs)
        return H

    def _compute_resistance_ratio(self, H, mu_dist, H0):
        """
        Compute R_ac/R_dc using power loss method for better numerical stability.

        Method: Compare AC power loss to DC power loss for the same total current.

        Current density: J = -dH/dz (Ampere's law)
        Total current (half slab): I = H(0) - H(a)

        AC power (per unit length/width, half slab):
            P_ac = integral_0^a { |J|^2 * rho } dz

        DC power for same current:
            J_dc = I / a  (uniform)
            P_dc = |J_dc|^2 * rho * a = |I|^2 * rho / a

        R_ac/R_dc = P_ac / P_dc = a * integral{|J|^2} dz / |I|^2

        This is the Cauchy-Schwarz inequality form:
            For uniform J: integral{|J|^2} dz * a = |integral{J} dz|^2 = |I|^2
            So R_ac/R_dc = 1 for DC (uniform J)

        For AC with skin effect: J is concentrated near surface, so R_ac > R_dc.
        """
        z = self.mesh_points
        n = self.n_nodes
        a = self.half_thickness

        if self.geometry == 'cylinder':
            # Cylinder: I = integral_0^R J*2pi*r dr, J = -dH/dr
            # Using Ampere: I_total = 2*pi*R*H(R) (for H_z with symmetry)
            # R_ac/R_dc = (R/2) * integral_0^R |J|^2*r dr / (|I/(pi*R^2)|^2 * pi*R^2)
            # Simpler: use P_ac/P_dc method
            I_total = H[n-1]  # H(R) = H0, total current per unit length = H0
            I_abs_sq = np.abs(I_total) ** 2
            if I_abs_sq < 1e-30:
                return self._compute_resistance_ratio_surface_impedance(H)
            J_sq_r_integral = 0.0
            for i in range(n - 1):
                dr = z[i + 1] - z[i]
                J_mid = -(H[i + 1] - H[i]) / dr
                r_mid = 0.5 * (z[i] + z[i + 1])
                J_sq_r_integral += np.abs(J_mid) ** 2 * r_mid * dr
            # P_dc for uniform J in cylinder: J_dc = I/(pi*R^2), P_dc = rho*|J_dc|^2*pi*R^2
            # P_ac = rho * 2*pi * integral |J|^2 * r dr
            # R_ac/R_dc = P_ac/P_dc = (R^2/2) * integral|J|^2*r dr / |I|^2
            # But I = H0 (Ampere, per unit length), I_total = 2*pi*R*H0
            # R_dc = rho * L / A = rho / (pi*R^2) per unit length
            # P_dc = |I_total|^2 * R_dc = 4*pi^2*R^2*|H0|^2 * rho/(pi*R^2) = 4*pi*rho*|H0|^2
            # P_ac = rho * 2*pi * integral|J|^2*r dr
            # R_ratio = 2*pi*rho*integral|J|^2*r dr / (4*pi*rho*|H0|^2)
            #         = integral|J|^2*r dr / (2*|H0|^2)
            R_ratio = J_sq_r_integral / (2 * I_abs_sq) * a
        else:
            # Slab: I = H(0) - H(a)
            I_total = H[0] - H[n-1]
            I_abs_sq = np.abs(I_total) ** 2
            if I_abs_sq < 1e-30:
                return self._compute_resistance_ratio_surface_impedance(H)
            J_sq_integral = 0.0
            for i in range(n - 1):
                dz = z[i + 1] - z[i]
                J_mid = -(H[i + 1] - H[i]) / dz
                J_sq_integral += np.abs(J_mid) ** 2 * dz
            R_ratio = a * J_sq_integral / I_abs_sq

        return max(R_ratio, 1.0)

    def _compute_resistance_ratio_surface_impedance(self, H):
        """
        Fallback method using surface impedance when current is very small.

        Z_s = -rho * dH/dz|_0 / H(0)
        R_ac/R_dc = Re(Z_s) * a / rho
        """
        z = self.mesh_points
        n = self.n_nodes
        a = self.half_thickness
        h = z[1] - z[0]

        # 3-point forward formula for dH/dz at z=0
        if n >= 3:
            dHdz_surface = (-3*H[0] + 4*H[1] - H[2]) / (2*h)
        else:
            dHdz_surface = (H[1] - H[0]) / h

        if abs(H[0]) < 1e-30:
            return 1.0

        Z_s = -self.rho * dHdz_surface / H[0]
        R_ratio = Z_s.real * a / self.rho

        return max(R_ratio, 1.0)

    def _compute_power_losses(self, H, mu_dist):
        """Compute power losses per unit surface area.

        Slab: P' = integral_0^a (1/2)*rho*|dH/dz|^2 dz  [W/m^2]
        Cylinder: P' = (1/R)*integral_0^R (1/2)*rho*|dH/dr|^2 * r dr  [W/m^2]
        """
        z = self.mesh_points
        n = self.n_nodes
        R = self.half_thickness

        P_ohmic = 0.0
        P_magnetic = 0.0
        Q_prime = 0.0

        for i in range(n - 1):
            dz = z[i + 1] - z[i]
            dHdz = (H[i + 1] - H[i]) / dz
            H_mid = 0.5 * (H[i] + H[i + 1])
            mu_mid = 0.5 * (mu_dist[i] + mu_dist[i + 1])
            H_sq = np.abs(H_mid) ** 2

            if self.geometry == 'cylinder':
                # r-weighted, normalized by R (per unit surface area)
                r_mid = 0.5 * (z[i] + z[i + 1])
                w = r_mid / R  # cylindrical weight / surface normalization
            else:
                w = 1.0  # slab: no weight

            P_ohmic += 0.5 * self.rho * np.abs(dHdz) ** 2 * w * dz

            if self.use_complex_mu:
                mu_prime = mu_mid.real
                mu_double_prime = -mu_mid.imag
                P_magnetic += 0.5 * self.omega * mu_double_prime * H_sq * w * dz
                Q_prime += 0.5 * self.omega * mu_prime * H_sq * w * dz
            else:
                Q_prime += 0.5 * self.omega * mu_mid.real * H_sq * w * dz

        P_prime = P_ohmic + P_magnetic
        return float(P_prime), float(Q_prime), float(P_magnetic)

    def get_dowell_ratio(self):
        """
        Get analytical R_ac/R_dc using Dowell's formula (linear material).

        R_ac/R_dc = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))

        where xi = a / delta
        """
        xi = self.xi

        if xi < 0.01:
            return 1.0

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


class ESIMCellProblemSolver:
    """
    Solves the 1D Cell Problem for ESIM using finite differences.

    The Cell Problem is a 1D boundary value problem on a semi-infinite domain
    that computes the effective surface impedance Z(H0) for a given surface
    tangential field amplitude H0.

    This implementation uses a finite difference method with geometric mesh
    grading for accuracy near the surface, where field gradients are largest.

    Supports both real and complex permeability:
    - Real mu: From BH curve (nonlinear saturation)
    - Complex mu = mu' - j*mu": Includes magnetic losses (hysteresis, etc.)

    NOTE: For finite-thickness conductors (transition region 0.5 < xi < 3),
    use ESIMFiniteSlabSolver instead for accurate R_ac/R_dc calculation.
    """

    def __init__(self, bh_curve=None, sigma=None, frequency=None,
                 num_skin_depths=10, n_nodes=100, complex_mu=None):
        """
        Initialize the Cell Problem solver.

        Parameters:
            bh_curve: BH curve data as [[H1, B1], [H2, B2], ...]
                      H in A/m, B in Tesla (for real permeability)
            sigma: Electrical conductivity [S/m]
            frequency: Operating frequency [Hz]
            num_skin_depths: Domain depth in skin depths (default: 10)
            n_nodes: Number of mesh nodes (default: 100)
            complex_mu: Complex permeability data (alternative to bh_curve)
                       Can be:
                       - (mu'_r, mu"_r) tuple for constant complex permeability
                       - [[H, mu'_r, mu"_r], ...] for H-dependent complex permeability
                       - dict with 'H', 'mu_prime', 'mu_double_prime' arrays
        """
        if not SCIPY_AVAILABLE:
            raise ImportError("Scipy is required for ESIM Cell Problem solver")

        if sigma is None or frequency is None:
            raise ValueError("sigma and frequency are required parameters")

        self.sigma = sigma
        self.frequency = frequency
        self.omega = 2 * np.pi * frequency
        self.rho = 1.0 / sigma  # Resistivity [Ohm*m]

        self.num_skin_depths = num_skin_depths
        self.n_nodes = n_nodes

        # Setup permeability model
        self.use_complex_mu = complex_mu is not None

        if self.use_complex_mu:
            # Use complex permeability model
            self.mu_interp = ComplexPermeabilityInterpolator(complex_mu, frequency)
            self.bh_interp = None  # Not used
            mu_initial = self.mu_interp.mu_initial
        else:
            # Use real BH curve model
            if bh_curve is None:
                raise ValueError("Either bh_curve or complex_mu must be provided")
            self.bh_interp = BHCurveInterpolator(bh_curve)
            self.mu_interp = None  # Not used
            mu_initial = self.bh_interp.mu_initial

        # Estimate initial skin depth using initial permeability
        # For complex mu, use |mu| for skin depth estimation
        mu_abs = abs(mu_initial)
        self.delta_initial = np.sqrt(2 * self.rho / (self.omega * mu_abs))

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

    def _get_mu(self, H_abs):
        """
        Get permeability for given |H|, handling both real and complex mu.

        Parameters:
            H_abs: Absolute value of magnetic field [A/m]

        Returns:
            mu: Permeability [H/m] (real or complex)
        """
        if self.use_complex_mu:
            return self.mu_interp.mu(H_abs)
        else:
            return self.bh_interp.mu(H_abs)

    def _get_mu_initial(self):
        """Get initial permeability (at H=0)."""
        if self.use_complex_mu:
            return self.mu_interp.mu_initial
        else:
            return self.bh_interp.mu_initial

    def solve(self, H0, tol=1e-6, max_iter=50, relaxation=0.5):
        """
        Solve the Cell Problem for given surface field H0.

        Uses Picard iteration for the nonlinear problem with finite difference
        discretization.

        The governing equation:
            rho * d^2H/ds^2 + j*omega*mu(|H|)*H = 0

        For complex mu = mu' - j*mu":
            rho * d^2H/ds^2 + j*omega*(mu' - j*mu")*H = 0
            rho * d^2H/ds^2 + (j*omega*mu' + omega*mu")*H = 0

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
                'P_magnetic': Magnetic power loss per unit area [W/m^2] (for complex mu)
                'H_solution': Array with H(s) distribution
                'converged': True if converged
                'iterations': Number of iterations
        """
        # Initial guess: exponential decay with initial skin depth
        delta = self.delta_initial
        H = H0 * np.exp(-self.mesh_points / delta) * np.exp(-1j * self.mesh_points / delta)

        # Initial permeability distribution (constant, may be complex)
        mu_initial = self._get_mu_initial()
        mu_dist = np.full(self.n_nodes, mu_initial, dtype=complex)

        converged = False

        for iteration in range(max_iter):
            # Build and solve the linear system with current mu distribution
            H_new = self._solve_linear_system(H0, mu_dist)

            # Update permeability based on |H| distribution
            mu_new = np.array([self._get_mu(abs(h)) for h in H_new], dtype=complex)

            # Check convergence (relative change in |mu|)
            rel_change = np.max(np.abs(mu_new - mu_dist)) / np.max(np.abs(mu_dist))

            # Under-relaxation (for complex mu, apply to both real and imag parts)
            mu_dist = (1 - relaxation) * mu_dist + relaxation * mu_new
            H = H_new

            if rel_change < tol:
                converged = True
                break

        # Compute power losses and impedance
        P_prime, Q_prime, P_magnetic = self._compute_power_losses(H, mu_dist)

        # Effective surface impedance: Z = 2*(P' + j*Q') / |H0|^2
        # For complex mu, P' includes both ohmic and magnetic losses
        Z = 2 * (P_prime + 1j * Q_prime) / (H0 ** 2)

        # Average final permeability (for reporting)
        mu_final = np.mean(mu_dist[:10])  # Average near surface (complex)

        return {
            'Z': Z,
            'P_prime': P_prime,
            'Q_prime': Q_prime,
            'P_magnetic': P_magnetic,
            'H_solution': H,
            'mesh_points': self.mesh_points,
            'converged': converged,
            'iterations': iteration + 1,
            'mu_final': mu_final,
            'use_complex_mu': self.use_complex_mu
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

        For real permeability:
            P' = (1/2) * integral{ rho * |dH/ds|^2 } ds  (ohmic loss)
            Q' = (omega/2) * integral{ mu * |H|^2 } ds   (reactive power)

        For complex permeability mu = mu' - j*mu":
            P_ohmic = (1/2) * integral{ rho * |dH/ds|^2 } ds
            P_magnetic = (omega/2) * integral{ mu" * |H|^2 } ds  (magnetic loss from mu")
            P' = P_ohmic + P_magnetic  (total active power)
            Q' = (omega/2) * integral{ mu' * |H|^2 } ds  (reactive power from mu')

        Parameters:
            H: Solution array H(s)
            mu_dist: Permeability distribution (may be complex)

        Returns:
            P_prime: Total active power loss per unit area [W/m^2]
            Q_prime: Reactive power per unit area [var/m^2]
            P_magnetic: Magnetic power loss per unit area [W/m^2] (from mu")
        """
        s = self.mesh_points
        n = self.n_nodes

        P_ohmic = 0.0
        P_magnetic = 0.0
        Q_prime = 0.0

        # Trapezoidal integration
        for i in range(n - 1):
            ds = s[i + 1] - s[i]

            # dH/ds at midpoint (central difference)
            dHds = (H[i + 1] - H[i]) / ds

            # H and mu at midpoint (average)
            H_mid = 0.5 * (H[i] + H[i + 1])
            mu_mid = 0.5 * (mu_dist[i] + mu_dist[i + 1])

            # |H|^2 at midpoint
            H_sq = np.abs(H_mid) ** 2

            # Ohmic power loss: P_ohmic = (1/2) * rho * |dH/ds|^2
            P_ohmic += 0.5 * self.rho * np.abs(dHds) ** 2 * ds

            if self.use_complex_mu:
                # For complex mu = mu' - j*mu":
                # mu' = Re(mu), mu" = -Im(mu) (note: mu = mu' - j*mu")
                mu_prime = mu_mid.real
                mu_double_prime = -mu_mid.imag  # mu" is positive

                # Magnetic power loss: P_mag = (omega/2) * mu" * |H|^2
                P_magnetic += 0.5 * self.omega * mu_double_prime * H_sq * ds

                # Reactive power: Q' = (omega/2) * mu' * |H|^2
                Q_prime += 0.5 * self.omega * mu_prime * H_sq * ds
            else:
                # For real mu: no magnetic loss, full mu contributes to Q'
                # Reactive power: Q' = (omega/2) * mu * |H|^2
                Q_prime += 0.5 * self.omega * mu_mid.real * H_sq * ds

        # Total active power = ohmic + magnetic
        P_prime = P_ohmic + P_magnetic

        return float(P_prime), float(Q_prime), float(P_magnetic)

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

    def get_linear_sibc(self, mu_r=None, complex_mu_r=None):
        """
        Get classical local SIBC impedance for comparison.

        For real linear materials:
            Z = (1+j) / (sigma * delta)
            where delta = sqrt(2*rho / (omega*mu))

        For complex permeability mu = mu' - j*mu":
            Z = sqrt(j*omega*mu / sigma)
            This gives both resistive and reactive components that
            depend on the complex mu.

        Parameters:
            mu_r: Real relative permeability (optional, uses initial if not given)
            complex_mu_r: Complex relative permeability (mu'_r, mu"_r) tuple

        Returns:
            Z_local: Complex local surface impedance [Ohm]
        """
        if complex_mu_r is not None:
            # Complex permeability case
            mu_p_r, mu_pp_r = complex_mu_r
            mu = mu_0 * complex(mu_p_r, -mu_pp_r)  # mu = mu' - j*mu"
            # Z = sqrt(j*omega*mu / sigma)
            Z_local = np.sqrt(1j * self.omega * mu / self.sigma)
            return Z_local

        elif self.use_complex_mu:
            # Use the complex mu from initialization
            mu = self._get_mu_initial()
            Z_local = np.sqrt(1j * self.omega * mu / self.sigma)
            return Z_local

        else:
            # Real permeability case
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
        # Store parameters (may be None if using pre-computed table)
        self.bh_curve = bh_curve
        self.sigma = sigma
        self.frequency = frequency

        if table_data is not None:
            self._load_table(table_data)
        elif bh_curve is not None and sigma is not None and frequency is not None:
            self._generate_table(bh_curve, sigma, frequency)
        else:
            raise ValueError("Provide either table_data or (bh_curve, sigma, frequency)")

    def _load_table(self, table_data):
        """Load pre-computed table."""
        # Convert to real values (discard any imaginary parts from numerical noise)
        table_array = np.array(table_data)
        if np.iscomplexobj(table_array):
            table_array = table_array.real
        self.table = np.ascontiguousarray(table_array, dtype=np.float64)
        self.H0_values = np.ascontiguousarray(self.table[:, 0].copy())
        self.Z_real = np.ascontiguousarray(self.table[:, 1].copy())
        self.Z_imag = np.ascontiguousarray(self.table[:, 2].copy())
        self.P_prime = np.ascontiguousarray(self.table[:, 3].copy())
        self.Q_prime = np.ascontiguousarray(self.table[:, 4].copy())

        # Ensure P_prime and Q_prime are positive for log interpolation
        # Use absolute value and store sign for Q_prime (can be negative in some cases)
        P_prime_safe = np.maximum(self.P_prime, 1e-20)
        Q_prime_safe = np.maximum(np.abs(self.Q_prime), 1e-20)

        # Create interpolators (log-H for better accuracy over wide range)
        log_H0 = np.ascontiguousarray(np.log10(self.H0_values + 1e-20))

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

    # Create ESITable with parameters stored
    esi_table = ESITable(table_data=table)
    esi_table.bh_curve = bh_curve
    esi_table.sigma = sigma
    esi_table.frequency = frequency

    return esi_table


# Example usage and test
if __name__ == "__main__":
    # Test with typical steel BH curve (real permeability)
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
    print("=" * 60)

    # ==========================================================
    # Test 1: Real permeability (BH curve)
    # ==========================================================
    print("\n" + "=" * 60)
    print("Test 1: Real Permeability (BH Curve)")
    print("=" * 60)
    print(f"Material: Steel (sigma = {sigma_steel/1e6:.1f} MS/m)")
    print(f"Frequency: {freq/1000:.1f} kHz")
    print()

    # Create solver with BH curve
    solver = ESIMCellProblemSolver(bh_curve=bh_curve_steel, sigma=sigma_steel, frequency=freq)

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

    # ==========================================================
    # Test 2: Complex permeability (constant)
    # ==========================================================
    print("\n" + "=" * 60)
    print("Test 2: Complex Permeability (Constant)")
    print("=" * 60)

    # Typical ferrite at 50 kHz: mu'_r = 2000, mu"_r = 200 (loss tangent = 0.1)
    mu_prime_r = 2000
    mu_double_prime_r = 200
    sigma_ferrite = 1e-2  # S/m (ferrite is nearly insulating)

    print(f"Material: Ferrite (mu'_r = {mu_prime_r}, mu\"_r = {mu_double_prime_r})")
    print(f"Loss tangent: tan(delta_m) = {mu_double_prime_r/mu_prime_r:.3f}")
    print(f"Conductivity: sigma = {sigma_ferrite} S/m")
    print(f"Frequency: {freq/1000:.1f} kHz")
    print()

    # Create solver with complex permeability
    solver_complex = ESIMCellProblemSolver(
        sigma=sigma_ferrite,
        frequency=freq,
        complex_mu=(mu_prime_r, mu_double_prime_r)
    )

    print(f"Initial skin depth: {solver_complex.delta_initial*1e3:.3f} mm")
    print(f"Domain length: {solver_complex.L*1e3:.1f} mm")
    print()

    # Solve for different H0 values
    print("Cell Problem Solutions with Complex mu:")
    print("-" * 85)
    print(f"{'H0 [A/m]':>12} {'Re(Z) [mOhm]':>14} {'Im(Z) [mOhm]':>14} {'P_total [W/m^2]':>16} {'P_mag [W/m^2]':>14}")
    print("-" * 85)

    for H0 in H0_test:
        result = solver_complex.solve(H0)
        Z = result['Z']
        P = result['P_prime']
        P_mag = result['P_magnetic']

        print(f"{H0:>12.0f} {Z.real*1e3:>14.4f} {Z.imag*1e3:>14.4f} {P:>16.2f} {P_mag:>14.2f}")

    print()

    # Compare with analytical SIBC for complex mu
    Z_linear_complex = solver_complex.get_linear_sibc()
    print(f"Linear SIBC (complex mu): Z = {Z_linear_complex.real*1e3:.4f} + j{Z_linear_complex.imag*1e3:.4f} mOhm")

    # ==========================================================
    # Test 3: H-dependent complex permeability
    # ==========================================================
    print("\n" + "=" * 60)
    print("Test 3: H-Dependent Complex Permeability")
    print("=" * 60)

    # H-dependent complex permeability (saturation and loss variation)
    # [H, mu'_r, mu"_r]
    complex_mu_data = [
        [0, 2000, 200],        # Low H: high permeability
        [100, 1800, 180],
        [500, 1500, 150],
        [1000, 1000, 100],
        [5000, 500, 50],
        [10000, 200, 20],
        [50000, 50, 5],        # High H: saturated
    ]

    sigma_test = 1e6  # S/m

    print(f"Material: H-dependent ferromagnetic")
    print(f"Conductivity: sigma = {sigma_test/1e6:.1f} MS/m")
    print(f"Frequency: {freq/1000:.1f} kHz")
    print()
    print("H-dependent mu' and mu\":")
    print(f"  H [A/m]    mu'_r     mu\"_r")
    for row in complex_mu_data:
        print(f"  {row[0]:>8.0f}  {row[1]:>6.0f}  {row[2]:>6.0f}")
    print()

    # Create solver with H-dependent complex permeability
    solver_Hdep = ESIMCellProblemSolver(
        sigma=sigma_test,
        frequency=freq,
        complex_mu=complex_mu_data
    )

    print(f"Initial skin depth: {solver_Hdep.delta_initial*1e3:.3f} mm")
    print()

    # Solve for different H0 values
    print("Cell Problem Solutions (H-dependent complex mu):")
    print("-" * 100)
    print(f"{'H0 [A/m]':>12} {'Re(Z) [mOhm]':>14} {'Im(Z) [mOhm]':>14} {'P_total [kW/m^2]':>18} {'P_mag [kW/m^2]':>16} {'Iter':>6}")
    print("-" * 100)

    for H0 in H0_test:
        result = solver_Hdep.solve(H0)
        Z = result['Z']
        P = result['P_prime']
        P_mag = result['P_magnetic']
        iters = result['iterations']

        print(f"{H0:>12.0f} {Z.real*1e3:>14.4f} {Z.imag*1e3:>14.4f} {P/1e3:>18.2f} {P_mag/1e3:>16.2f} {iters:>6}")

    print()
    print("Test completed!")
