"""
Unified Analysis Framework for Radia

This module provides a unified interface for electromagnetic analysis:
- Static analysis (DC)
- Frequency response analysis (AC sweep)
- Transient analysis (CLN time-domain)

The framework supports:
- PEEC conductor models with skin effect
- MMM magnetic material models
- Coupled PEEC-MMM (CplMag) systems
- CLN (Cauer Ladder Network) model order reduction

References:
    [1] A.E. Ruehli, "Equivalent Circuit Models for Three-Dimensional
        Multiconductor Systems", IEEE Trans. MTT, Vol. 22, No. 3, 1974.
    [2] P. Feldmann, R.W. Freund, "Efficient Linear Circuit Analysis by
        Pade Approximation via the Lanczos Process", IEEE TCAD, Vol. 14, 1995.
    [3] K. Hollaus, "Transient Analysis of Electromagnetic Fields",
        TU Wien, 2003.

Author: Radia Development Team
License: LGPL-2.1
"""

import numpy as np
from enum import Enum, auto
from dataclasses import dataclass, field
from typing import Optional, List, Tuple, Union, Callable
from abc import ABC, abstractmethod

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7  # H/m
EPS_0 = 8.854187817e-12     # F/m
INV_FOUR_PI = 1.0 / (4.0 * np.pi)


class AnalysisType(Enum):
    """Enumeration of analysis types."""
    STATIC = auto()      # DC analysis
    FREQUENCY = auto()   # AC frequency sweep
    TRANSIENT = auto()   # Time-domain transient


class SolverType(Enum):
    """Enumeration of solver types."""
    PEEC = auto()        # Pure PEEC conductor
    MMM = auto()         # Pure MMM magnetic
    CPLMAG = auto()      # Coupled PEEC-MMM


@dataclass
class AnalysisResult:
    """Base class for analysis results."""
    analysis_type: AnalysisType
    solver_type: SolverType
    success: bool
    message: str = ""
    computation_time: float = 0.0


@dataclass
class StaticResult(AnalysisResult):
    """Result of static (DC) analysis."""
    resistance: float = 0.0           # DC resistance [Ohm]
    inductance: float = 0.0           # DC inductance [H]
    current: Optional[np.ndarray] = None  # Current distribution

    def __post_init__(self):
        self.analysis_type = AnalysisType.STATIC


@dataclass
class FrequencyResult(AnalysisResult):
    """Result of frequency response analysis."""
    frequencies: np.ndarray = field(default_factory=lambda: np.array([]))
    impedance: np.ndarray = field(default_factory=lambda: np.array([]))  # Complex Z(f)
    resistance: np.ndarray = field(default_factory=lambda: np.array([]))  # Re(Z)
    reactance: np.ndarray = field(default_factory=lambda: np.array([]))   # Im(Z)
    inductance: np.ndarray = field(default_factory=lambda: np.array([]))  # L = Im(Z)/(2*pi*f)

    def __post_init__(self):
        self.analysis_type = AnalysisType.FREQUENCY

    @property
    def quality_factor(self) -> np.ndarray:
        """Q factor = X_L / R = omega*L / R."""
        with np.errstate(divide='ignore', invalid='ignore'):
            Q = np.abs(self.reactance) / np.abs(self.resistance)
            Q[~np.isfinite(Q)] = 0.0
        return Q


@dataclass
class TransientResult(AnalysisResult):
    """Result of transient analysis."""
    time: np.ndarray = field(default_factory=lambda: np.array([]))
    current: np.ndarray = field(default_factory=lambda: np.array([]))     # i(t)
    voltage: np.ndarray = field(default_factory=lambda: np.array([]))     # v(t)
    flux_linkage: np.ndarray = field(default_factory=lambda: np.array([]))  # psi(t)
    power: np.ndarray = field(default_factory=lambda: np.array([]))       # p(t) = v*i

    def __post_init__(self):
        self.analysis_type = AnalysisType.TRANSIENT

    @property
    def energy(self) -> np.ndarray:
        """Cumulative energy E(t) = integral(p*dt)."""
        if len(self.time) < 2:
            return np.array([0.0])
        dt = np.diff(self.time)
        power_avg = 0.5 * (self.power[:-1] + self.power[1:])
        return np.concatenate([[0.0], np.cumsum(power_avg * dt)])


class AnalysisSolver(ABC):
    """Abstract base class for analysis solvers."""

    def __init__(self):
        self._is_built = False

    @abstractmethod
    def build(self) -> bool:
        """Build the system matrices."""
        pass

    @abstractmethod
    def solve_static(self) -> StaticResult:
        """Perform static (DC) analysis."""
        pass

    @abstractmethod
    def solve_frequency(self, frequencies: np.ndarray) -> FrequencyResult:
        """Perform frequency response analysis."""
        pass

    @abstractmethod
    def solve_transient(self, time: np.ndarray,
                        excitation: Callable[[float], float]) -> TransientResult:
        """Perform transient analysis."""
        pass


class PEECAnalysisSolver(AnalysisSolver):
    """
    PEEC-based analysis solver with CLN model order reduction.

    Supports:
    - Static DC analysis
    - Frequency response with skin effect (Dowell or SIBC)
    - Transient analysis using CLN (Cauer Ladder Network)

    The CLN method represents the frequency-dependent impedance as a
    continued fraction expansion, enabling efficient time-domain simulation.
    """

    def __init__(self):
        super().__init__()

        # PEEC matrices
        self._L = None  # Inductance matrix
        self._R = None  # Resistance matrix
        self._P = None  # Potential coefficient matrix
        self._M_LS = None  # Loop-Star coupling matrix

        # CLN parameters
        self._cln_order = 5  # Lanczos iterations
        self._cln_result = None

        # Skin effect parameters
        self._sigma = 5.8e7  # Conductivity [S/m] (copper default)
        self._use_dowell = True
        self._conductor_width = None
        self._conductor_height = None

    def set_peec_matrices(self, L: np.ndarray, R: np.ndarray,
                          P: Optional[np.ndarray] = None,
                          M_LS: Optional[np.ndarray] = None):
        """
        Set the PEEC system matrices directly.

        Parameters:
            L: Inductance matrix [H]
            R: Resistance matrix [Ohm]
            P: Potential coefficient matrix [1/F] (optional)
            M_LS: Loop-Star coupling matrix (optional)
        """
        self._L = np.asarray(L)
        self._R = np.asarray(R)
        self._P = np.asarray(P) if P is not None else None
        self._M_LS = np.asarray(M_LS) if M_LS is not None else None
        self._is_built = False

    def set_cln_order(self, order: int):
        """Set the CLN model order (number of Lanczos iterations)."""
        if order < 1:
            raise ValueError("CLN order must be >= 1")
        self._cln_order = order
        self._is_built = False

    def set_conductivity(self, sigma: float):
        """Set conductor conductivity [S/m]."""
        self._sigma = sigma

    def set_skin_effect_model(self, use_dowell: bool = True,
                               width: Optional[float] = None,
                               height: Optional[float] = None):
        """
        Configure skin effect model.

        Parameters:
            use_dowell: Use Dowell formula (True) or SIBC (False)
            width: Conductor width for Dowell model [m]
            height: Conductor height for Dowell model [m]
        """
        self._use_dowell = use_dowell
        self._conductor_width = width
        self._conductor_height = height

    def build(self) -> bool:
        """Build CLN model from PEEC matrices."""
        if self._L is None or self._R is None:
            raise ValueError("PEEC matrices not set. Call set_peec_matrices() first.")

        try:
            # Import CLN module
            try:
                from . import cln_core
            except ImportError:
                import cln_core

            # Perform Lanczos reduction
            # Ensure matrices are contiguous float64
            L = np.ascontiguousarray(self._L, dtype=np.float64)
            R = np.ascontiguousarray(self._R, dtype=np.float64)

            # Handle diagonal R matrix
            if R.ndim == 1:
                R_diag = R
            else:
                R_diag = np.diag(R)

            self._cln_result = cln_core.lanczos(R_diag, L, self._cln_order)
            self._is_built = True
            return True

        except Exception as e:
            # Fallback: use Python implementation
            self._cln_result = self._lanczos_python(self._R, self._L, self._cln_order)
            self._is_built = True
            return True

    def _lanczos_python(self, R: np.ndarray, L: np.ndarray,
                        n_iter: int) -> dict:
        """
        Pure Python Lanczos algorithm for CLN reduction.

        The Lanczos process transforms (R, L) into tridiagonal form:
            R -> R_diag (diagonal)
            L -> L_tridiag (symmetric tridiagonal)

        This enables efficient frequency sweep via continued fractions.
        """
        n = L.shape[0]

        # Ensure R is diagonal
        if R.ndim == 2:
            R_diag = np.diag(R)
        else:
            R_diag = R.copy()

        # Initialize Lanczos vectors
        Q = np.zeros((n, n_iter + 1))
        alpha = np.zeros(n_iter)
        beta = np.zeros(n_iter + 1)

        # Starting vector (normalized)
        q = np.ones(n) / np.sqrt(n)
        Q[:, 0] = q

        # Lanczos iteration with K-orthogonalization
        for j in range(n_iter):
            # w = L * q_j
            w = L @ Q[:, j]

            # alpha_j = q_j^T * K * w where K = R (mass matrix)
            alpha[j] = Q[:, j] @ (R_diag * w)

            # Orthogonalize: w = w - alpha_j * q_j - beta_j * q_{j-1}
            w = w - alpha[j] * Q[:, j]
            if j > 0:
                w = w - beta[j] * Q[:, j-1]

            # beta_{j+1} = ||w||_K
            beta[j+1] = np.sqrt(w @ (R_diag * w))

            if beta[j+1] < 1e-14:
                # Early termination (invariant subspace found)
                n_iter = j + 1
                break

            Q[:, j+1] = w / beta[j+1]

        # Build tridiagonal L matrix
        L_tridiag = np.diag(alpha[:n_iter])
        for j in range(n_iter - 1):
            L_tridiag[j, j+1] = beta[j+1]
            L_tridiag[j+1, j] = beta[j+1]

        # R remains diagonal in reduced space
        R_reduced = np.ones(n_iter)  # Identity in K-orthogonal basis

        return {
            'R_diag': R_reduced,
            'L_tridiag': L_tridiag,
            'Q': Q[:, :n_iter],
            'alpha': alpha[:n_iter],
            'beta': beta[1:n_iter+1],
            'n_iter': n_iter
        }

    def solve_static(self) -> StaticResult:
        """
        Perform static (DC) analysis.

        At DC:
            Z_DC = R_DC (resistance only)
            L_DC = L (full inductance)
        """
        import time
        t_start = time.time()

        if self._L is None or self._R is None:
            return StaticResult(
                analysis_type=AnalysisType.STATIC,
                solver_type=SolverType.PEEC,
                success=False,
                message="PEEC matrices not set"
            )

        # DC resistance: sum of diagonal elements (series connection)
        if self._R.ndim == 2:
            R_dc = np.sum(np.diag(self._R))
        else:
            R_dc = np.sum(self._R)

        # DC inductance: sum of all elements (total flux linkage)
        L_dc = np.sum(self._L)

        t_elapsed = time.time() - t_start

        return StaticResult(
            analysis_type=AnalysisType.STATIC,
            solver_type=SolverType.PEEC,
            success=True,
            message="DC analysis completed",
            computation_time=t_elapsed,
            resistance=R_dc,
            inductance=L_dc
        )

    def solve_frequency(self, frequencies: np.ndarray) -> FrequencyResult:
        """
        Perform frequency response analysis using CLN.

        The CLN impedance is computed as a continued fraction:
            Z(s) = R_0 + s*L_0 / (1 + s*L_1/R_1 / (1 + ...))

        With Dowell correction for skin effect:
            Z(s) = R_dc * F_R(delta) + j*omega*L_dc * F_L(delta)

        where delta = sqrt(omega*mu*sigma) * thickness / sqrt(2)
        """
        import time
        t_start = time.time()

        if not self._is_built:
            self.build()

        frequencies = np.asarray(frequencies)
        n_freq = len(frequencies)

        # Get CLN parameters
        if isinstance(self._cln_result, dict):
            R_diag = self._cln_result['R_diag']
            L_tridiag = self._cln_result['L_tridiag']
        else:
            R_diag = self._cln_result.R_diag
            L_tridiag = self._cln_result.L_tridiag

        # Compute impedance at each frequency
        Z = np.zeros(n_freq, dtype=complex)

        for i, f in enumerate(frequencies):
            omega = 2.0 * np.pi * f
            s = 1j * omega

            # CLN continued fraction evaluation
            Z[i] = self._evaluate_cln_impedance(s, R_diag, L_tridiag)

            # Apply Dowell correction if enabled
            if self._use_dowell and f > 0 and self._conductor_height is not None:
                Z[i] = self._apply_dowell_correction(Z[i], f)

        R_array = np.real(Z)
        X_array = np.imag(Z)

        # Compute inductance L = X / (2*pi*f)
        L_array = np.zeros(n_freq)
        for i, f in enumerate(frequencies):
            if f > 0:
                L_array[i] = X_array[i] / (2.0 * np.pi * f)

        t_elapsed = time.time() - t_start

        return FrequencyResult(
            analysis_type=AnalysisType.FREQUENCY,
            solver_type=SolverType.PEEC,
            success=True,
            message=f"Frequency sweep completed ({n_freq} points)",
            computation_time=t_elapsed,
            frequencies=frequencies,
            impedance=Z,
            resistance=R_array,
            reactance=X_array,
            inductance=L_array
        )

    def _evaluate_cln_impedance(self, s: complex,
                                 R_diag: np.ndarray,
                                 L_tridiag: np.ndarray) -> complex:
        """
        Evaluate CLN impedance at complex frequency s.

        Uses backward recursion for continued fraction:
            Z_n = R_n + s*L_nn
            Z_k = R_k + s*L_kk + (s*L_k,k+1)^2 / Z_{k+1}
        """
        n = len(R_diag)

        if n == 1:
            return R_diag[0] + s * L_tridiag[0, 0]

        # Backward recursion
        Z = R_diag[n-1] + s * L_tridiag[n-1, n-1]

        for k in range(n-2, -1, -1):
            L_kk = L_tridiag[k, k]
            L_k_kp1 = L_tridiag[k, k+1] if k < n-1 else 0.0

            Z_self = R_diag[k] + s * L_kk
            coupling = (s * L_k_kp1) ** 2 / Z if abs(Z) > 1e-30 else 0.0
            Z = Z_self + coupling

        return Z

    def _apply_dowell_correction(self, Z_dc: complex, freq: float) -> complex:
        """
        Apply Dowell skin effect correction.

        F_R(xi) = xi * (sinh(2*xi) + sin(2*xi)) / (cosh(2*xi) - cos(2*xi))
        F_L(xi) = (3/(2*xi)) * (sinh(2*xi) - sin(2*xi)) / (cosh(2*xi) - cos(2*xi))

        where xi = h / delta, delta = sqrt(2 / (omega*mu*sigma))
        """
        omega = 2.0 * np.pi * freq
        delta = np.sqrt(2.0 / (omega * MU_0 * self._sigma))

        h = self._conductor_height
        xi = h / delta

        if xi < 0.01:
            # Small argument: F_R -> 1, F_L -> 1
            return Z_dc
        elif xi > 50:
            # Large argument: F_R -> xi, F_L -> 3/(2*xi)
            F_R = xi
            F_L = 1.5 / xi
        else:
            # General case
            sinh_2xi = np.sinh(2.0 * xi)
            sin_2xi = np.sin(2.0 * xi)
            cosh_2xi = np.cosh(2.0 * xi)
            cos_2xi = np.cos(2.0 * xi)

            denom = cosh_2xi - cos_2xi
            if abs(denom) < 1e-30:
                denom = 1e-30

            F_R = xi * (sinh_2xi + sin_2xi) / denom
            F_L = (1.5 / xi) * (sinh_2xi - sin_2xi) / denom

        R_dc = np.real(Z_dc)
        L_dc = np.imag(Z_dc) / omega if omega > 0 else 0.0

        R_ac = R_dc * F_R
        L_ac = L_dc * F_L

        return R_ac + 1j * omega * L_ac

    def solve_transient(self, time: np.ndarray,
                        excitation: Callable[[float], float],
                        v_or_i: str = 'v') -> TransientResult:
        """
        Perform transient analysis using CLN state-space model.

        The CLN model transforms the frequency-domain impedance into
        a state-space representation suitable for time-domain simulation:

            dx/dt = A*x + B*u
            y = C*x + D*u

        where x is the internal state vector, u is the excitation,
        and y is the response (current or voltage).

        Parameters:
            time: Time points [s]
            excitation: Function u(t) returning excitation value
            v_or_i: 'v' for voltage excitation, 'i' for current excitation
        """
        import time as time_module
        t_start = time_module.time()

        if not self._is_built:
            self.build()

        time_array = np.asarray(time)
        n_time = len(time_array)

        # Get CLN parameters
        if isinstance(self._cln_result, dict):
            R_diag = self._cln_result['R_diag']
            L_tridiag = self._cln_result['L_tridiag']
            n_states = self._cln_result['n_iter']
        else:
            R_diag = np.asarray(self._cln_result.R_diag)
            L_tridiag = np.asarray(self._cln_result.L_tridiag)
            n_states = len(R_diag)

        # Build state-space matrices for CLN
        # State: x = [i_1, i_2, ..., i_n] (branch currents)
        # For voltage excitation: v = Z(s) * i -> L*di/dt + R*i = v

        # A = -L^{-1} * R
        # B = L^{-1}
        # C = [1, 0, ..., 0]^T (output is total current)
        # D = 0

        try:
            L_inv = np.linalg.inv(L_tridiag)
        except np.linalg.LinAlgError:
            L_inv = np.linalg.pinv(L_tridiag)

        R_mat = np.diag(R_diag)
        A = -L_inv @ R_mat
        B = L_inv[:, 0]  # First column (excitation enters at first node)
        C = np.zeros(n_states)
        C[0] = 1.0  # Output is first state

        # Initialize state and output arrays
        x = np.zeros(n_states)
        current = np.zeros(n_time)
        voltage = np.zeros(n_time)
        flux = np.zeros(n_time)

        # Time integration using Backward Euler (implicit, stable)
        for i, t in enumerate(time_array):
            u = excitation(t)

            if v_or_i == 'v':
                voltage[i] = u
            else:
                current[i] = u

            if i > 0:
                dt = time_array[i] - time_array[i-1]

                # Backward Euler: (I - dt*A)*x_new = x_old + dt*B*u
                I_minus_dtA = np.eye(n_states) - dt * A
                rhs = x + dt * B * u

                try:
                    x = np.linalg.solve(I_minus_dtA, rhs)
                except np.linalg.LinAlgError:
                    x = np.linalg.lstsq(I_minus_dtA, rhs, rcond=None)[0]

            if v_or_i == 'v':
                current[i] = C @ x
            else:
                # For current excitation, compute voltage
                voltage[i] = R_diag[0] * current[i] + L_tridiag[0, 0] * (
                    (current[i] - current[i-1]) / (time_array[i] - time_array[i-1])
                    if i > 0 else 0.0
                )

            # Flux linkage: psi = L * i (total flux)
            flux[i] = L_tridiag[0, 0] * current[i]

        # Power
        power = voltage * current

        t_elapsed = time_module.time() - t_start

        return TransientResult(
            analysis_type=AnalysisType.TRANSIENT,
            solver_type=SolverType.PEEC,
            success=True,
            message=f"Transient analysis completed ({n_time} time points)",
            computation_time=t_elapsed,
            time=time_array,
            current=current,
            voltage=voltage,
            flux_linkage=flux,
            power=power
        )


class UnifiedAnalysis:
    """
    High-level unified analysis interface.

    Provides a simple API for common analysis tasks:

    Example:
        >>> from radia.analysis import UnifiedAnalysis
        >>>
        >>> # Create analysis object
        >>> analysis = UnifiedAnalysis()
        >>>
        >>> # Set PEEC model
        >>> analysis.set_peec_model(L, R)
        >>>
        >>> # Static analysis
        >>> result = analysis.static()
        >>> print(f"DC Resistance: {result.resistance:.3e} Ohm")
        >>>
        >>> # Frequency sweep
        >>> freqs = np.logspace(0, 6, 100)  # 1 Hz to 1 MHz
        >>> result = analysis.frequency_sweep(freqs)
        >>>
        >>> # Transient analysis
        >>> time = np.linspace(0, 1e-3, 1000)  # 1 ms
        >>> voltage = lambda t: 1.0 if t > 0 else 0.0  # Step function
        >>> result = analysis.transient(time, voltage)
    """

    def __init__(self):
        self._solver = PEECAnalysisSolver()
        self._is_configured = False

    def set_peec_model(self, L: np.ndarray, R: np.ndarray,
                       P: Optional[np.ndarray] = None,
                       M_LS: Optional[np.ndarray] = None,
                       cln_order: int = 5):
        """
        Configure PEEC model for analysis.

        Parameters:
            L: Inductance matrix [H]
            R: Resistance matrix or diagonal [Ohm]
            P: Potential coefficient matrix [1/F] (optional)
            M_LS: Loop-Star coupling matrix (optional)
            cln_order: CLN model order (default: 5)
        """
        self._solver.set_peec_matrices(L, R, P, M_LS)
        self._solver.set_cln_order(cln_order)
        self._is_configured = True

    def set_skin_effect(self, sigma: float = 5.8e7,
                        use_dowell: bool = True,
                        conductor_height: Optional[float] = None):
        """
        Configure skin effect model.

        Parameters:
            sigma: Conductivity [S/m] (default: copper 5.8e7)
            use_dowell: Use Dowell formula (default: True)
            conductor_height: Conductor height for Dowell [m]
        """
        self._solver.set_conductivity(sigma)
        self._solver.set_skin_effect_model(use_dowell, height=conductor_height)

    def static(self) -> StaticResult:
        """Perform static (DC) analysis."""
        if not self._is_configured:
            raise RuntimeError("Model not configured. Call set_peec_model() first.")
        return self._solver.solve_static()

    def frequency_sweep(self, frequencies: np.ndarray) -> FrequencyResult:
        """
        Perform frequency response analysis.

        Parameters:
            frequencies: Array of frequencies [Hz]

        Returns:
            FrequencyResult with impedance, resistance, reactance, inductance
        """
        if not self._is_configured:
            raise RuntimeError("Model not configured. Call set_peec_model() first.")
        return self._solver.solve_frequency(frequencies)

    def transient(self, time: np.ndarray,
                  excitation: Callable[[float], float],
                  excitation_type: str = 'voltage') -> TransientResult:
        """
        Perform transient analysis using CLN method.

        Parameters:
            time: Time points [s]
            excitation: Function u(t) returning excitation value
            excitation_type: 'voltage' or 'current'

        Returns:
            TransientResult with time, current, voltage, flux, power
        """
        if not self._is_configured:
            raise RuntimeError("Model not configured. Call set_peec_model() first.")

        v_or_i = 'v' if excitation_type.lower().startswith('v') else 'i'
        return self._solver.solve_transient(time, excitation, v_or_i)


# Convenience functions for common waveforms
def step_voltage(amplitude: float = 1.0, t_start: float = 0.0) -> Callable[[float], float]:
    """Create step voltage excitation."""
    return lambda t: amplitude if t >= t_start else 0.0


def pulse_voltage(amplitude: float = 1.0, t_start: float = 0.0,
                  duration: float = 1e-3) -> Callable[[float], float]:
    """Create pulse voltage excitation."""
    return lambda t: amplitude if t_start <= t < t_start + duration else 0.0


def sinusoidal_voltage(amplitude: float = 1.0, frequency: float = 1000.0,
                       phase: float = 0.0) -> Callable[[float], float]:
    """Create sinusoidal voltage excitation."""
    omega = 2.0 * np.pi * frequency
    return lambda t: amplitude * np.sin(omega * t + phase)


def ramp_voltage(slope: float = 1.0, t_start: float = 0.0,
                 max_value: Optional[float] = None) -> Callable[[float], float]:
    """Create ramp voltage excitation."""
    def ramp(t):
        if t < t_start:
            return 0.0
        v = slope * (t - t_start)
        if max_value is not None:
            v = min(v, max_value)
        return v
    return ramp
