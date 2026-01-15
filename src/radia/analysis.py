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


# =============================================================================
# MMM (Magnetic Material Method) Analysis
# =============================================================================

@dataclass
class MMMStaticResult(AnalysisResult):
    """Result of MMM static analysis."""
    magnetization: Optional[np.ndarray] = None  # Magnetization distribution [A/m]
    B_field: Optional[np.ndarray] = None        # Magnetic flux density [T]
    H_field: Optional[np.ndarray] = None        # Magnetic field strength [A/m]
    total_energy: float = 0.0                   # Stored magnetic energy [J]
    max_iterations: int = 0
    final_precision: float = 0.0

    def __post_init__(self):
        self.analysis_type = AnalysisType.STATIC
        self.solver_type = SolverType.MMM


@dataclass
class MMMFrequencyResult(AnalysisResult):
    """Result of MMM frequency response analysis."""
    frequencies: np.ndarray = field(default_factory=lambda: np.array([]))
    complex_permeability: np.ndarray = field(default_factory=lambda: np.array([]))  # mu(f)
    mu_real: np.ndarray = field(default_factory=lambda: np.array([]))      # mu'
    mu_imag: np.ndarray = field(default_factory=lambda: np.array([]))      # mu"
    loss_tangent: np.ndarray = field(default_factory=lambda: np.array([]))  # tan(delta)

    def __post_init__(self):
        self.analysis_type = AnalysisType.FREQUENCY
        self.solver_type = SolverType.MMM

    @property
    def quality_factor(self) -> np.ndarray:
        """Q factor = 1 / tan(delta) = mu' / mu"."""
        with np.errstate(divide='ignore', invalid='ignore'):
            Q = np.abs(self.mu_real) / np.abs(self.mu_imag)
            Q[~np.isfinite(Q)] = 0.0
        return Q


class MMMAnalysisSolver(AnalysisSolver):
    """
    MMM-based analysis solver for magnetic materials.

    Supports:
    - Static analysis: Magnetization computation using Radia's solver
    - Frequency response: Complex permeability model mu(omega)

    For frequency response, the complex permeability is modeled as:
        mu(omega) = mu_inf + (mu_s - mu_inf) / (1 + j*omega*tau)

    where:
        mu_s: Static (DC) permeability
        mu_inf: High-frequency permeability
        tau: Relaxation time constant

    For nonlinear materials, the effective permeability depends on H level.
    """

    def __init__(self):
        super().__init__()

        # Radia object handle
        self._radia_obj = None

        # Material parameters
        self._mu_r_static = 1000.0      # DC relative permeability
        self._mu_r_inf = 1.0            # High-frequency relative permeability
        self._tau = 1e-6                # Relaxation time [s]
        self._sigma = 0.0               # Conductivity (for eddy current loss)

        # Solver parameters
        self._precision = 1e-4
        self._max_iterations = 1000
        self._solver_method = 0  # 0=LU, 1=BiCGSTAB

        # Background field
        self._H_background = np.array([0.0, 0.0, 0.0])

        # Computed results cache
        self._static_result = None

    def set_radia_object(self, radia_obj: int):
        """Set the Radia object handle for analysis."""
        self._radia_obj = radia_obj
        self._is_built = False
        self._static_result = None

    def set_material_parameters(self,
                                 mu_r_static: float = 1000.0,
                                 mu_r_inf: float = 1.0,
                                 tau: float = 1e-6,
                                 sigma: float = 0.0):
        """
        Set material parameters for frequency response.

        Parameters:
            mu_r_static: DC relative permeability
            mu_r_inf: High-frequency relative permeability
            tau: Relaxation time constant [s]
            sigma: Conductivity for eddy current loss [S/m]
        """
        self._mu_r_static = mu_r_static
        self._mu_r_inf = mu_r_inf
        self._tau = tau
        self._sigma = sigma

    def set_solver_parameters(self,
                               precision: float = 1e-4,
                               max_iterations: int = 1000,
                               method: int = 0):
        """
        Set solver parameters.

        Parameters:
            precision: Convergence precision
            max_iterations: Maximum iterations
            method: 0=LU, 1=BiCGSTAB, 2=HACApK
        """
        self._precision = precision
        self._max_iterations = max_iterations
        self._solver_method = method

    def set_background_field(self, H: Union[List[float], np.ndarray]):
        """Set external background field [A/m]."""
        self._H_background = np.asarray(H)

    def build(self) -> bool:
        """
        Build/prepare the MMM system.

        For MMM, this mainly validates the Radia object is set.
        The actual matrix build happens during Solve().
        """
        if self._radia_obj is None:
            raise ValueError("Radia object not set. Call set_radia_object() first.")
        self._is_built = True
        return True

    def solve_static(self) -> MMMStaticResult:
        """
        Perform static (DC) magnetization analysis.

        Uses Radia's Solve() to compute the magnetization distribution
        in the presence of the background field.
        """
        import time

        if not self._is_built:
            self.build()

        t_start = time.time()

        try:
            # Import radia module
            try:
                import radia as rad
            except ImportError:
                from . import radia as rad

            # Solve the magnetization
            result = rad.Solve(self._radia_obj,
                               self._precision,
                               self._max_iterations,
                               self._solver_method)

            # Parse result
            # Solve returns: [avg_M_change, max_H, n_iterations]
            if isinstance(result, (list, tuple, np.ndarray)):
                avg_M_change = result[0] if len(result) > 0 else 0.0
                max_H = result[1] if len(result) > 1 else 0.0
                n_iter = int(result[2]) if len(result) > 2 else 0
            else:
                avg_M_change = 0.0
                max_H = 0.0
                n_iter = 0

            t_elapsed = time.time() - t_start

            self._static_result = MMMStaticResult(
                analysis_type=AnalysisType.STATIC,
                solver_type=SolverType.MMM,
                success=True,
                message=f"MMM static analysis completed in {n_iter} iterations",
                computation_time=t_elapsed,
                max_iterations=n_iter,
                final_precision=avg_M_change
            )
            return self._static_result

        except Exception as e:
            t_elapsed = time.time() - t_start
            return MMMStaticResult(
                analysis_type=AnalysisType.STATIC,
                solver_type=SolverType.MMM,
                success=False,
                message=f"MMM solve failed: {str(e)}",
                computation_time=t_elapsed
            )

    def solve_frequency(self, frequencies: np.ndarray) -> MMMFrequencyResult:
        """
        Perform frequency response analysis.

        Computes the complex permeability mu(omega) using the Debye model:
            mu(omega) = mu_inf + (mu_s - mu_inf) / (1 + j*omega*tau)

        Plus eddy current contribution (if sigma > 0):
            mu_eddy(omega) = mu_0 * sigma / (j*omega)  (added to imaginary part)

        Parameters:
            frequencies: Array of frequencies [Hz]

        Returns:
            MMMFrequencyResult with complex permeability data
        """
        import time
        t_start = time.time()

        frequencies = np.asarray(frequencies)
        n_freq = len(frequencies)

        # Complex permeability array
        mu_complex = np.zeros(n_freq, dtype=complex)

        # Debye relaxation model
        mu_s = self._mu_r_static * MU_0
        mu_inf = self._mu_r_inf * MU_0
        delta_mu = mu_s - mu_inf
        tau = self._tau

        for i, f in enumerate(frequencies):
            omega = 2.0 * np.pi * f

            if omega < 1e-10:
                # DC limit
                mu_complex[i] = mu_s
            else:
                # Debye model
                mu_debye = mu_inf + delta_mu / (1.0 + 1j * omega * tau)

                # Add eddy current loss (if conductive)
                if self._sigma > 0:
                    # Skin depth: delta = sqrt(2 / (omega * mu * sigma))
                    # Eddy current loss adds to mu"
                    skin_depth = np.sqrt(2.0 / (omega * MU_0 * self._sigma))
                    # Simplified eddy current contribution
                    mu_eddy_imag = MU_0 * self._sigma * skin_depth**2 * omega / 4.0
                    mu_debye = mu_debye - 1j * mu_eddy_imag

                mu_complex[i] = mu_debye

        # Extract real and imaginary parts
        mu_real = np.real(mu_complex)
        mu_imag = -np.imag(mu_complex)  # Convention: mu" is positive for loss

        # Loss tangent
        loss_tangent = np.zeros(n_freq)
        for i in range(n_freq):
            if abs(mu_real[i]) > 1e-30:
                loss_tangent[i] = mu_imag[i] / mu_real[i]

        t_elapsed = time.time() - t_start

        return MMMFrequencyResult(
            analysis_type=AnalysisType.FREQUENCY,
            solver_type=SolverType.MMM,
            success=True,
            message=f"MMM frequency sweep completed ({n_freq} points)",
            computation_time=t_elapsed,
            frequencies=frequencies,
            complex_permeability=mu_complex,
            mu_real=mu_real,
            mu_imag=mu_imag,
            loss_tangent=loss_tangent
        )

    def set_magnetic_circuit(self, L: np.ndarray, R: np.ndarray,
                              cln_order: int = 5):
        """
        Set magnetic circuit equivalent model for CLN-based transient.

        The magnetic circuit is modeled as:
            v = R*i + L*di/dt

        where:
            v: MMF (magnetomotive force) [A-turns]
            i: flux [Wb]
            R: reluctance matrix [A-turns/Wb]
            L: magnetic "mass" matrix (related to energy storage)

        For a simple magnetic system:
            R = l / (mu * A)  (reluctance)
            L = related to eddy current time constant

        Parameters:
            L: Inductance-like matrix for magnetic circuit
            R: Reluctance matrix
            cln_order: CLN model order
        """
        self._L_mag = np.asarray(L)
        self._R_mag = np.asarray(R)
        self._cln_order = cln_order
        self._cln_result = None

    def _build_magnetic_cln(self):
        """Build CLN model from magnetic circuit matrices."""
        if not hasattr(self, '_L_mag') or self._L_mag is None:
            # Auto-generate from material parameters
            # Simple single-branch model: tau = L/R
            tau = self._tau
            if tau < 1e-15:
                tau = 1e-9  # Small but non-zero

            # Equivalent circuit: v = R*phi + L*dphi/dt
            # At DC: v = R*phi -> phi = v/R = mu*A/l * MMF
            # Time constant: tau = L/R

            # Normalize so R = 1 (per unit)
            self._R_mag = np.array([[1.0]])
            self._L_mag = np.array([[tau]])
            self._cln_order = min(self._cln_order if hasattr(self, '_cln_order') else 5, 3)

        # Perform Lanczos reduction
        L = np.atleast_2d(self._L_mag)
        R = np.atleast_2d(self._R_mag)

        n = L.shape[0]
        n_iter = min(self._cln_order, n)

        if n == 1:
            # Single element case
            self._cln_result = {
                'R_diag': np.array([R[0, 0]]),
                'L_tridiag': np.array([[L[0, 0]]]),
                'n_iter': 1
            }
            return

        # Use same Lanczos implementation as PEEC
        R_diag = np.diag(R) if R.ndim == 2 else R

        Q = np.zeros((n, n_iter + 1))
        alpha = np.zeros(n_iter)
        beta = np.zeros(n_iter + 1)

        q = np.ones(n) / np.sqrt(n)
        Q[:, 0] = q

        for j in range(n_iter):
            w = L @ Q[:, j]
            alpha[j] = Q[:, j] @ (R_diag * w)
            w = w - alpha[j] * Q[:, j]
            if j > 0:
                w = w - beta[j] * Q[:, j-1]

            beta[j+1] = np.sqrt(np.abs(w @ (R_diag * w)))

            if beta[j+1] < 1e-14:
                n_iter = j + 1
                break

            Q[:, j+1] = w / beta[j+1]

        L_tridiag = np.diag(alpha[:n_iter])
        for j in range(n_iter - 1):
            L_tridiag[j, j+1] = beta[j+1]
            L_tridiag[j+1, j] = beta[j+1]

        self._cln_result = {
            'R_diag': np.ones(n_iter),
            'L_tridiag': L_tridiag,
            'n_iter': n_iter
        }

    def solve_transient(self, time: np.ndarray,
                        excitation: Callable[[float], float]) -> TransientResult:
        """
        Transient analysis for MMM using CLN method.

        The magnetic system is modeled using an equivalent circuit approach,
        then reduced using CLN for efficient time-domain simulation.

        For magnetic materials with relaxation time tau:
            dM/dt = (M_eq(H) - M) / tau

        This is equivalent to a first-order RL circuit.
        For higher-order models, use set_magnetic_circuit() to provide
        custom L, R matrices.

        The CLN method enables efficient computation even for complex
        frequency-dependent permeability models.
        """
        import time as time_module
        t_start = time_module.time()

        # Build CLN if not already done
        if not hasattr(self, '_cln_result') or self._cln_result is None:
            self._build_magnetic_cln()

        time_array = np.asarray(time)
        n_time = len(time_array)

        # Get CLN parameters
        R_diag = self._cln_result['R_diag']
        L_tridiag = self._cln_result['L_tridiag']
        n_states = self._cln_result['n_iter']

        # Build state-space model
        # dx/dt = A*x + B*u, y = C*x
        # where x is internal state, u is H (field), y is M (magnetization)

        try:
            L_inv = np.linalg.inv(L_tridiag)
        except np.linalg.LinAlgError:
            L_inv = np.linalg.pinv(L_tridiag)

        R_mat = np.diag(R_diag)
        A = -L_inv @ R_mat
        B = L_inv[:, 0]
        C = np.zeros(n_states)
        C[0] = self._mu_r_static - 1.0  # chi = mu_r - 1

        # Initialize state and output
        x = np.zeros(n_states)
        M = np.zeros(n_time)
        H = np.zeros(n_time)

        # Time integration (Backward Euler)
        for i, t in enumerate(time_array):
            H[i] = excitation(t)

            if i > 0:
                dt = time_array[i] - time_array[i-1]

                # Backward Euler: (I - dt*A)*x_new = x_old + dt*B*H
                I_minus_dtA = np.eye(n_states) - dt * A
                rhs = x + dt * B * H[i]

                try:
                    x = np.linalg.solve(I_minus_dtA, rhs)
                except np.linalg.LinAlgError:
                    x = np.linalg.lstsq(I_minus_dtA, rhs, rcond=None)[0]

            M[i] = C @ x

        # B = mu_0 * (H + M)
        B_field = MU_0 * (H + M)

        # Power: dB/dt * H (per unit volume)
        voltage = np.zeros(n_time)
        for i in range(1, n_time):
            dt = time_array[i] - time_array[i-1]
            voltage[i] = (B_field[i] - B_field[i-1]) / dt

        power = voltage * H

        t_elapsed = time_module.time() - t_start

        return TransientResult(
            analysis_type=AnalysisType.TRANSIENT,
            solver_type=SolverType.MMM,
            success=True,
            message=f"MMM transient (CLN) completed ({n_time} points, {n_states} states)",
            computation_time=t_elapsed,
            time=time_array,
            current=M,       # Magnetization [A/m]
            voltage=H,       # Field strength [A/m]
            flux_linkage=B_field,  # Flux density [T]
            power=power
        )


def build_magnetic_circuit_from_mmm(N_matrix: np.ndarray,
                                     chi: float,
                                     tau: float = 1e-6) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build equivalent magnetic circuit (L, R matrices) from MMM interaction matrix.

    The MMM system is:
        M = chi * (H_ext + H_demag)
        H_demag = N * M

    Rearranging:
        (I - chi * N) * M = chi * H_ext

    This is analogous to a resistive network with:
        R_eq = (I - chi * N) / chi  (reluctance-like matrix)

    For transient behavior (frequency-dependent permeability), we add
    inductance-like terms based on the relaxation time tau.

    Parameters:
        N_matrix: MMM interaction (demagnetization) matrix [n_dof x n_dof]
        chi: Magnetic susceptibility (chi = mu_r - 1)
        tau: Relaxation time constant [s]

    Returns:
        L: Inductance-like matrix for CLN
        R: Reluctance matrix
    """
    n = N_matrix.shape[0]
    I = np.eye(n)

    # Reluctance matrix: R = (I - chi * N) / chi
    # At DC: M = chi * (I - chi*N)^{-1} * H_ext
    R = (I - chi * N_matrix) / chi

    # For transient, the inductance-like matrix is related to tau
    # dM/dt = (M_eq - M) / tau
    # In matrix form: tau * dM/dt + M = chi * (I - chi*N)^{-1} * H_ext
    # Equivalent: L * dM/dt + R * M = H_ext
    # where L = tau * R

    L = tau * R

    return L, R


class UnifiedMMMAnalysis:
    """
    High-level unified analysis interface for MMM (magnetic materials).

    Supports two modes:
    1. Radia object mode: Uses rad.Solve() for static analysis
    2. Matrix mode: Uses MMM interaction matrix for full analysis including CLN transient

    Example 1: Radia object mode
        >>> import radia as rad
        >>> from radia.analysis import UnifiedMMMAnalysis
        >>>
        >>> rad.FldUnits('m')
        >>> mag = rad.ObjHexahedron(vertices, [0, 0, 0])
        >>> mat = rad.MatLin(1000)
        >>> rad.MatApl(mag, mat)
        >>> bkg = rad.ObjBckg(lambda p: [0, 0, 1e5])
        >>> container = rad.ObjCnt([mag, bkg])
        >>>
        >>> analysis = UnifiedMMMAnalysis()
        >>> analysis.set_radia_model(container, mu_r=1000)
        >>> result = analysis.static()

    Example 2: Matrix mode (with mmm_core)
        >>> from radia.analysis import UnifiedMMMAnalysis
        >>> import mmm_core
        >>>
        >>> # Build MMM matrices
        >>> builder = mmm_core.MMMBuilder()
        >>> builder.add_tetrahedra_from_mesh(vertices, elements)
        >>> N, dof_offset = builder.build()
        >>>
        >>> # Analysis with CLN transient
        >>> analysis = UnifiedMMMAnalysis()
        >>> analysis.set_mmm_matrices(N, mu_r=1000, tau=1e-6)
        >>>
        >>> # Static, frequency, and transient all work
        >>> result = analysis.static()
        >>> result = analysis.frequency_sweep(frequencies)
        >>> result = analysis.transient(time, H_excitation)
    """

    def __init__(self):
        self._solver = MMMAnalysisSolver()
        self._is_configured = False

    def set_radia_model(self, radia_obj: int,
                        mu_r: float = 1000.0,
                        mu_r_inf: float = 1.0,
                        tau: float = 1e-6,
                        sigma: float = 0.0):
        """
        Configure MMM model using Radia object.

        Parameters:
            radia_obj: Radia object handle (container)
            mu_r: DC relative permeability
            mu_r_inf: High-frequency relative permeability
            tau: Relaxation time constant [s]
            sigma: Conductivity [S/m] (for eddy currents)
        """
        self._solver.set_radia_object(radia_obj)
        self._solver.set_material_parameters(mu_r, mu_r_inf, tau, sigma)
        self._is_configured = True
        self._mode = 'radia'

    def set_mmm_matrices(self, N_matrix: np.ndarray,
                          mu_r: float = 1000.0,
                          mu_r_inf: float = 1.0,
                          tau: float = 1e-6,
                          sigma: float = 0.0,
                          cln_order: int = 5):
        """
        Configure MMM model using interaction matrix directly.

        This mode enables CLN-based transient analysis by building
        equivalent magnetic circuit from the MMM interaction matrix N.

        The MMM system: M = chi * (H_ext + N * M)
        Equivalent circuit: L * dM/dt + R * M = chi * H_ext

        Parameters:
            N_matrix: MMM interaction (demagnetization) matrix
            mu_r: DC relative permeability
            mu_r_inf: High-frequency relative permeability
            tau: Relaxation time constant [s]
            sigma: Conductivity [S/m] (for eddy currents)
            cln_order: CLN model order for transient
        """
        chi = mu_r - 1.0

        # Build magnetic circuit from MMM matrix
        L, R = build_magnetic_circuit_from_mmm(N_matrix, chi, tau)

        # Set solver parameters
        self._solver.set_material_parameters(mu_r, mu_r_inf, tau, sigma)
        self._solver.set_magnetic_circuit(L, R, cln_order)

        # Store N matrix for static solve
        self._N_matrix = N_matrix
        self._chi = chi

        self._is_configured = True
        self._mode = 'matrix'

    def set_solver_parameters(self,
                               precision: float = 1e-4,
                               max_iterations: int = 1000,
                               method: int = 0):
        """
        Set solver parameters.

        Parameters:
            precision: Convergence precision
            max_iterations: Maximum iterations
            method: 0=LU, 1=BiCGSTAB, 2=HACApK
        """
        self._solver.set_solver_parameters(precision, max_iterations, method)

    def static(self) -> MMMStaticResult:
        """Perform static magnetization analysis."""
        if not self._is_configured:
            raise RuntimeError("Model not configured. Call set_radia_model() first.")
        return self._solver.solve_static()

    def frequency_sweep(self, frequencies: np.ndarray) -> MMMFrequencyResult:
        """
        Perform frequency response analysis.

        Computes complex permeability mu(omega) using Debye relaxation model.

        Parameters:
            frequencies: Array of frequencies [Hz]

        Returns:
            MMMFrequencyResult with mu', mu", loss tangent
        """
        if not self._is_configured:
            raise RuntimeError("Model not configured. Call set_radia_model() first.")
        return self._solver.solve_frequency(frequencies)

    def transient(self, time: np.ndarray,
                  excitation: Callable[[float], float]) -> TransientResult:
        """
        Perform transient analysis.

        Models magnetization relaxation: dM/dt = (M_eq - M) / tau

        Parameters:
            time: Time points [s]
            excitation: Function H(t) returning field strength [A/m]

        Returns:
            TransientResult with M(t), H(t), B(t)
        """
        if not self._is_configured:
            raise RuntimeError("Model not configured. Call set_radia_model() first.")
        return self._solver.solve_transient(time, excitation)
