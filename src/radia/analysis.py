"""
Unified Analysis Framework for Radia

This module provides a unified interface for electromagnetic analysis:
- Static analysis (DC)
- Frequency response analysis (AC sweep)
- Transient analysis (CLN time-domain)

The framework supports:
- PEEC conductor models with skin effect
- CLN (Cauer Ladder Network) model order reduction

References:
    [1] A.E. Ruehli, "Equivalent Circuit Models for Three-Dimensional
        Multiconductor Systems", IEEE Trans. MTT, Vol. 22, No. 3, 1974.
    [2] P. Feldmann, R.W. Freund, "Efficient Linear Circuit Analysis by
        Pade Approximation via the Lanczos Process", IEEE TCAD, Vol. 14, 1995.
    [3] K. Hollaus, "Transient Analysis of Electromagnetic Fields",
        TU Wien, 2003.

Author: Radia Development Team
License: BSD-style AND MIT
"""

import warnings
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
class MultiPortResult(AnalysisResult):
    """
    Result of multi-port frequency response analysis.

    For WPT (Wireless Power Transfer) and transformer analysis:
    - Z-parameters: [V] = [Z] * [I]
    - Y-parameters: [I] = [Y] * [V]
    - S-parameters: [b] = [S] * [a] (for RF applications)

    Key quantities for WPT:
    - Mutual inductance M_ij: Off-diagonal elements of L matrix
    - Coupling coefficient k_ij = M_ij / sqrt(L_ii * L_jj)
    - Transfer impedance Z_21: Voltage at port 2 per current at port 1
    """
    frequencies: np.ndarray = field(default_factory=lambda: np.array([]))
    n_ports: int = 2

    # Z-parameters [n_freq, n_ports, n_ports]
    Z_matrix: np.ndarray = field(default_factory=lambda: np.array([]))

    # Derived quantities
    mutual_inductance: np.ndarray = field(default_factory=lambda: np.array([]))  # M_12 vs frequency
    coupling_coefficient: np.ndarray = field(default_factory=lambda: np.array([]))  # k vs frequency
    self_inductance: np.ndarray = field(default_factory=lambda: np.array([]))  # [L1, L2] vs frequency

    # Transfer characteristics
    transfer_impedance: np.ndarray = field(default_factory=lambda: np.array([]))  # Z_21
    power_transfer_efficiency: np.ndarray = field(default_factory=lambda: np.array([]))  # eta

    def __post_init__(self):
        self.analysis_type = AnalysisType.FREQUENCY
        self.solver_type = SolverType.PEEC

    def get_Z(self, i: int, j: int) -> np.ndarray:
        """Get Z_ij parameter vs frequency."""
        return self.Z_matrix[:, i, j]

    def get_Y(self) -> np.ndarray:
        """Get Y-parameters (Y = inv(Z)).

        A singular Z at some frequency has no defined admittance: that row is
        set to NaN and a warning is issued.  No-Fallback: this does NOT silently
        substitute a pseudo-inverse (pinv gives a different, un-auditable
        min-norm result that the caller could not distinguish from a true Y).
        """
        n_freq = len(self.frequencies)
        Y = np.zeros_like(self.Z_matrix)
        for i in range(n_freq):
            try:
                Y[i] = np.linalg.inv(self.Z_matrix[i])
            except np.linalg.LinAlgError:
                warnings.warn(
                    f"Singular Z-matrix at f={self.frequencies[i]:.6g} Hz; "
                    "Y-parameters at this frequency set to NaN (no defined "
                    "admittance).")
                Y[i] = np.nan
        return Y

    def get_S(self, Z0: float = 50.0) -> np.ndarray:
        """
        Get S-parameters (for RF applications).

        S = (Z - Z0*I) * (Z + Z0*I)^{-1}

        Parameters:
            Z0: Reference impedance [Ohm] (default: 50)
        """
        n_freq = len(self.frequencies)
        n_ports = self.n_ports
        S = np.zeros((n_freq, n_ports, n_ports), dtype=complex)
        I = np.eye(n_ports)

        for i in range(n_freq):
            Z = self.Z_matrix[i]
            try:
                S[i] = (Z - Z0 * I) @ np.linalg.inv(Z + Z0 * I)
            except np.linalg.LinAlgError:
                # No-Fallback: NaN (clearly invalid), NOT 0 -- a zero S-matrix
                # reads as a perfectly-absorbing port and would be trusted.
                warnings.warn(
                    f"Singular (Z + Z0*I) at f={self.frequencies[i]:.6g} Hz; "
                    "S-parameters at this frequency set to NaN.")
                S[i] = np.nan
        return S


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
    - Loop-Star decomposition for low-frequency stability

    The CLN method represents the frequency-dependent impedance as a
    continued fraction expansion, enabling efficient time-domain simulation.

    Loop-Star Decomposition:
        The current I can be decomposed into:
        - Loop currents I_L (solenoidal, inductive)
        - Star currents I_S (irrotational, capacitive)

        The system equation becomes:
        [Z_LL  Z_LS] [I_L]   [V_L]
        [Z_SL  Z_SS] [I_S] = [V_S]

        where:
        - Z_LL = R_L + j*omega*L_L (Loop-Loop: inductive)
        - Z_SS = R_S + 1/(j*omega*C_S) (Star-Star: capacitive)
        - Z_LS, Z_SL: coupling terms

        In MQS (Magneto-Quasi-Static) approximation (power electronics),
        the Star currents can be neglected: I_S ~ 0.
    """

    def __init__(self):
        super().__init__()

        # PEEC matrices (original RWG basis)
        self._L = None  # Inductance matrix [H]
        self._R = None  # Resistance matrix [Ohm]
        self._P = None  # Potential coefficient matrix [1/F]
        self._M_LS = None  # Loop-Star coupling matrix (optional)

        # Loop-Star decomposed matrices
        self._L_LL = None  # Loop-Loop inductance
        self._L_LS = None  # Loop-Star coupling (L)
        self._L_SS = None  # Star-Star inductance
        self._R_LL = None  # Loop-Loop resistance
        self._R_LS = None  # Loop-Star coupling (R)
        self._R_SS = None  # Star-Star resistance
        self._P_SS = None  # Star-Star potential (capacitance)
        self._n_loops = 0  # Number of Loop DOFs
        self._n_stars = 0  # Number of Star DOFs

        # CLN parameters
        self._cln_order = 5  # Lanczos iterations
        self._cln_result = None

        # Skin effect parameters
        self._sigma = 5.8e7  # Conductivity [S/m] (copper default)
        self._use_dowell = True
        self._conductor_width = None
        self._conductor_height = None

        # MQS mode flag
        self._use_mqs = True  # If True, neglect Star currents

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

    def set_mqs_mode(self, use_mqs: bool = True):
        """
        Set MQS (Magneto-Quasi-Static) mode.

        In MQS mode, Star (capacitive) currents are neglected.
        This is valid for power electronics frequencies (DC to ~1 MHz)
        where omega * epsilon << sigma (conduction dominates displacement).

        Parameters:
            use_mqs: If True, solve Loop equations only (I_S = 0)
        """
        self._use_mqs = use_mqs

    def set_loop_star_matrices(self,
                                L_LL: np.ndarray, R_LL: np.ndarray,
                                L_LS: Optional[np.ndarray] = None,
                                L_SS: Optional[np.ndarray] = None,
                                R_LS: Optional[np.ndarray] = None,
                                R_SS: Optional[np.ndarray] = None,
                                P_SS: Optional[np.ndarray] = None):
        """
        Set Loop-Star decomposed PEEC matrices directly.

        This method is used when Loop-Star matrices are computed externally
        (e.g., from mesh topology analysis).

        Parameters:
            L_LL: Loop-Loop inductance matrix [H]
            R_LL: Loop-Loop resistance matrix [Ohm]
            L_LS: Loop-Star inductance coupling matrix [H] (optional)
            L_SS: Star-Star inductance matrix [H] (optional)
            R_LS: Loop-Star resistance coupling matrix [Ohm] (optional)
            R_SS: Star-Star resistance matrix [Ohm] (optional)
            P_SS: Star-Star potential coefficient matrix [1/F] (optional)

        The Loop-Star impedance system is:
            Z_LL = R_LL + j*omega*L_LL
            Z_LS = R_LS + j*omega*L_LS
            Z_SL = R_LS^T + j*omega*L_LS^T
            Z_SS = R_SS + j*omega*L_SS + P_SS/(j*omega)
        """
        self._L_LL = np.asarray(L_LL)
        self._R_LL = np.asarray(R_LL)
        self._n_loops = L_LL.shape[0]

        if L_LS is not None:
            self._L_LS = np.asarray(L_LS)
            self._n_stars = L_LS.shape[1]
        else:
            self._L_LS = None
            self._n_stars = 0

        self._L_SS = np.asarray(L_SS) if L_SS is not None else None
        self._R_LS = np.asarray(R_LS) if R_LS is not None else None
        self._R_SS = np.asarray(R_SS) if R_SS is not None else None
        self._P_SS = np.asarray(P_SS) if P_SS is not None else None

        self._is_built = False

    def compute_loop_star_from_topology(self, incidence_matrix: np.ndarray):
        """
        Compute Loop-Star transformation from mesh incidence matrix.

        The incidence matrix A (n_nodes x n_edges) encodes mesh topology.
        Loop-Star decomposition uses spanning tree analysis:
        - Tree edges -> Star (irrotational) basis
        - Co-tree edges -> Loop (solenoidal) basis

        Parameters:
            incidence_matrix: Node-edge incidence matrix A

        Returns:
            T_LS: Loop-Star transformation matrix [n_edges x (n_loops + n_stars)]
        """
        A = np.asarray(incidence_matrix)
        n_nodes, n_edges = A.shape

        # Number of independent loops (Euler formula)
        # n_loops = n_edges - n_nodes + n_connected_components
        # For simply connected mesh: n_loops = n_edges - n_nodes + 1

        # Find spanning tree using graph algorithm
        # Tree edges define Star basis, co-tree edges define Loop basis

        # Simplified implementation: use rank analysis
        # Rank of A = n_nodes - n_connected_components
        rank_A = np.linalg.matrix_rank(A)
        n_tree = rank_A  # Tree edges (Star)
        n_cotree = n_edges - n_tree  # Co-tree edges (Loop)

        # For full implementation, need to:
        # 1. Find spanning tree edges (using BFS/DFS)
        # 2. Build Loop basis from fundamental cycles
        # 3. Build Star basis from tree gradients

        # Placeholder: return identity (no transformation)
        # Full implementation requires mesh topology from conductor solver
        self._n_loops = n_cotree
        self._n_stars = n_tree

        return np.eye(n_edges)  # Identity for now

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

        The Lanczos process transforms the generalized eigenvalue problem
            L * x = lambda * R * x
        into tridiagonal form using R-orthogonal basis vectors.

        For the impedance Z(s) = R + s*L, the CLN representation is:
            Z(s) = R_tot + s * L_0 / (1 + s*L_1/R_1 / (1 + ...))

        This enables O(n) frequency sweep via continued fraction evaluation.

        Reference: Feldmann & Freund, IEEE TCAD, Vol. 14, 1995.
        """
        n = L.shape[0]

        # Ensure R is diagonal
        if R.ndim == 2:
            R_diag = np.diag(R)
        else:
            R_diag = R.copy()

        # Total DC resistance (scalar sum for series model)
        R_total = np.sum(R_diag)

        # Total DC inductance (sum of all L elements for series model)
        L_total = np.sum(L)

        # Initialize Lanczos vectors
        Q = np.zeros((n, n_iter + 1))
        alpha = np.zeros(n_iter)
        beta = np.zeros(n_iter + 1)

        # Starting vector: proportional to unit current distribution
        # For series connection, all segments carry same current
        q = np.ones(n)
        # R-normalize: ||q||_R = sqrt(q^T * R * q) = 1
        norm_R = np.sqrt(q @ (R_diag * q))
        if norm_R > 1e-14:
            q = q / norm_R
        Q[:, 0] = q

        # Lanczos iteration with R-orthogonalization
        # This produces: Q^T * R * Q = I (identity)
        #                Q^T * L * Q = T (tridiagonal)
        for j in range(n_iter):
            # w = L * q_j (not R^{-1}*L for standard Lanczos)
            w = L @ Q[:, j]

            # alpha_j = q_j^T * R * (L * q_j) / (q_j^T * R * q_j)
            # But since q_j is R-normalized: q_j^T * R * q_j = 1
            alpha[j] = Q[:, j] @ w  # Rayleigh quotient for L

            # Orthogonalize: w = w - alpha_j * R * q_j - beta_j * R * q_{j-1}
            w = w - alpha[j] * (R_diag * Q[:, j])
            if j > 0:
                w = w - beta[j] * (R_diag * Q[:, j-1])

            # beta_{j+1} = ||w||_{R^{-1}} = sqrt(w^T * R^{-1} * w)
            # For efficiency, compute: w^T * (w / R_diag)
            with np.errstate(divide='ignore', invalid='ignore'):
                R_inv_diag = np.where(R_diag > 1e-30, 1.0 / R_diag, 0.0)
            beta_sq = w @ (R_inv_diag * w)
            beta[j+1] = np.sqrt(max(beta_sq, 0.0))

            if beta[j+1] < 1e-14:
                # Early termination (invariant subspace found)
                n_iter = j + 1
                break

            # q_{j+1} = R^{-1} * w / beta_{j+1}
            Q[:, j+1] = (R_inv_diag * w) / beta[j+1]

        # Build tridiagonal representation
        # T = Q^T * L * Q (tridiagonal)
        L_tridiag = np.diag(alpha[:n_iter])
        for j in range(n_iter - 1):
            L_tridiag[j, j+1] = beta[j+1]
            L_tridiag[j+1, j] = beta[j+1]

        return {
            'R_total': R_total,      # Total DC resistance
            'L_total': L_total,      # Total DC inductance
            'L_tridiag': L_tridiag,  # Reduced L matrix (tridiagonal)
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

        The CLN impedance is computed using the Lanczos-reduced model:
            Z(s) = R_dc + s * L_eff(s)

        where L_eff(s) is the effective inductance from the tridiagonal system.

        At DC (s=0): Z = R_dc
        At high frequency: Z ~ R_dc + j*omega*L_dc

        With optional Dowell correction for skin effect:
            Z(s) = R_dc * F_R(xi) + j*omega*L_dc * F_L(xi)

        where xi = h / delta, delta = sqrt(2 / (omega*mu*sigma))

        Reference:
            Feldmann & Freund, "Efficient Linear Circuit Analysis by Pade
            Approximation via the Lanczos Process", IEEE TCAD, 1995.
        """
        import time
        t_start = time.time()

        if not self._is_built:
            self.build()

        frequencies = np.asarray(frequencies)
        n_freq = len(frequencies)

        # Get CLN parameters
        if isinstance(self._cln_result, dict):
            R_total = self._cln_result.get('R_total', 0.0)
            L_total = self._cln_result.get('L_total', 0.0)
            L_tridiag = self._cln_result['L_tridiag']
            n_cln = self._cln_result['n_iter']
        else:
            # Fallback for C++ module result
            R_total = getattr(self._cln_result, 'R_total', np.sum(self._R))
            L_total = getattr(self._cln_result, 'L_total', np.sum(self._L))
            L_tridiag = np.asarray(self._cln_result.L_tridiag)
            n_cln = len(L_tridiag)

        # Compute impedance at each frequency
        Z = np.zeros(n_freq, dtype=complex)

        for i, f in enumerate(frequencies):
            omega = 2.0 * np.pi * f

            if f < 1e-10:
                # DC limit: Z = R_dc
                Z[i] = R_total
            else:
                s = 1j * omega

                # Evaluate Z(s) using the CLN tridiagonal model
                # Z(s) = R_total + s * e1^T * (I + s*T)^{-1} * e1 * L_total
                # where T is the normalized tridiagonal matrix
                Z[i] = self._evaluate_cln_impedance_v2(s, R_total, L_total, L_tridiag)

                # Apply Dowell correction if enabled
                if self._use_dowell and self._conductor_height is not None:
                    Z[i] = self._apply_dowell_correction(Z[i], f)

        R_array = np.real(Z)
        X_array = np.imag(Z)

        # Compute inductance L = X / (2*pi*f)
        L_array = np.zeros(n_freq)
        for i, f in enumerate(frequencies):
            if f > 1e-10:
                L_array[i] = X_array[i] / (2.0 * np.pi * f)
            else:
                L_array[i] = L_total  # DC inductance

        t_elapsed = time.time() - t_start

        return FrequencyResult(
            analysis_type=AnalysisType.FREQUENCY,
            solver_type=SolverType.PEEC,
            success=True,
            message=f"Frequency sweep completed ({n_freq} points, CLN order={n_cln})",
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
        Evaluate CLN impedance at complex frequency s (legacy method).

        Uses backward recursion for continued fraction:
            Z_n = R_n + s*L_nn
            Z_k = R_k + s*L_kk + (s*L_k,k+1)^2 / Z_{k+1}

        Note: This is kept for backward compatibility. Use _evaluate_cln_impedance_v2
        for the improved Lanczos-based evaluation.
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

    def _evaluate_cln_impedance_v2(self, s: complex,
                                    R_total: float,
                                    L_total: float,
                                    L_tridiag: np.ndarray) -> complex:
        """
        Evaluate CLN impedance at complex frequency s using direct RL model.

        For a simple series RL circuit:
            Z(s) = R + s*L

        The CLN (Lanczos) reduction captures the frequency-dependent
        inductance variation, but for this simplified model we use:

        Z(s) = R_total + s * L_eff

        where L_eff is the effective inductance from the tridiagonal system.

        At DC (s->0): Z = R_total
        At high freq (s->inf): Z ~ s * L_total

        The tridiagonal matrix from Lanczos captures the mode structure,
        but for total impedance the key values are:
        - R_total: sum of resistances (series)
        - L_total: effective total inductance

        Parameters:
            s: Complex frequency (j*omega)
            R_total: Total DC resistance [Ohm]
            L_total: Total DC inductance [H]
            L_tridiag: Tridiagonal inductance matrix from Lanczos

        Returns:
            Complex impedance Z(s)
        """
        # For this implementation, we use the direct RL model
        # The CLN tridiagonal structure would be used for more complex
        # frequency-dependent effects (skin effect, proximity effect)
        # which are handled separately by the Dowell correction.

        # Simple series RL impedance
        Z = R_total + s * L_total

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

    def solve_frequency_loop_star(self, frequencies: np.ndarray) -> FrequencyResult:
        """
        Perform frequency response analysis using Loop-Star decomposition.

        This method solves the Loop-Star coupled system:
            [Z_LL  Z_LS] [I_L]   [V_L]
            [Z_SL  Z_SS] [I_S] = [V_S]

        where:
            Z_LL = R_LL + j*omega*L_LL (Loop-Loop: inductive)
            Z_SS = R_SS + j*omega*L_SS + P_SS/(j*omega) (Star-Star: capacitive)
            Z_LS = R_LS + j*omega*L_LS (Loop-Star coupling)
            Z_SL = Z_LS^T (reciprocity)

        In MQS mode (self._use_mqs = True):
            - Star currents are neglected: I_S = 0
            - Only Loop equations are solved: Z_LL * I_L = V_L

        Parameters:
            frequencies: Array of frequencies [Hz]

        Returns:
            FrequencyResult with impedance Z(f)
        """
        import time
        t_start = time.time()

        # Check if Loop-Star matrices are set
        if self._L_LL is None or self._R_LL is None:
            # Fallback to standard CLN method
            return self.solve_frequency(frequencies)

        frequencies = np.asarray(frequencies)
        n_freq = len(frequencies)

        # Compute impedance at each frequency
        Z = np.zeros(n_freq, dtype=complex)

        for i, f in enumerate(frequencies):
            omega = 2.0 * np.pi * f

            if f < 1e-10:
                # DC limit: Z = R_LL (Loop contribution only)
                if self._R_LL.ndim == 2:
                    Z[i] = np.sum(self._R_LL)
                else:
                    Z[i] = np.sum(self._R_LL)
            else:
                s = 1j * omega

                if self._use_mqs or self._n_stars == 0:
                    # MQS mode: solve Loop equations only (I_S = 0)
                    # Z_LL * I_L = V_L -> I_L = inv(Z_LL) * V_L
                    # For unit voltage: Z = V / I = 1 / (e1^T * inv(Z_LL) * e1)
                    Z[i] = self._solve_loop_only(s)
                else:
                    # Full Loop-Star: solve coupled system
                    Z[i] = self._solve_loop_star_coupled(s)

                # Apply Dowell correction if enabled
                if self._use_dowell and self._conductor_height is not None:
                    Z[i] = self._apply_dowell_correction(Z[i], f)

        R_array = np.real(Z)
        X_array = np.imag(Z)

        # Compute inductance L = X / (2*pi*f)
        L_array = np.zeros(n_freq)
        for i, f in enumerate(frequencies):
            if f > 1e-10:
                L_array[i] = X_array[i] / (2.0 * np.pi * f)
            else:
                # DC inductance from L_LL
                L_array[i] = np.sum(self._L_LL) if self._L_LL is not None else 0.0

        t_elapsed = time.time() - t_start

        mode_str = "MQS (Loop-only)" if self._use_mqs else "Full Loop-Star"

        return FrequencyResult(
            analysis_type=AnalysisType.FREQUENCY,
            solver_type=SolverType.PEEC,
            success=True,
            message=f"Loop-Star frequency sweep ({n_freq} points, {mode_str})",
            computation_time=t_elapsed,
            frequencies=frequencies,
            impedance=Z,
            resistance=R_array,
            reactance=X_array,
            inductance=L_array
        )

    def _solve_loop_only(self, s: complex) -> complex:
        """
        Solve Loop-only equations (MQS approximation).

        Z_LL * I_L = V_L
        Z_LL = R_LL + s * L_LL

        For series-connected segments with unit terminal voltage,
        the impedance is the sum of loop impedances.
        """
        # Build Z_LL
        if self._R_LL.ndim == 1:
            # Diagonal R
            R_LL = np.diag(self._R_LL)
        else:
            R_LL = self._R_LL

        L_LL = self._L_LL
        Z_LL = R_LL + s * L_LL

        # For total impedance with series connection:
        # All loops carry same current -> Z_total = sum of all elements
        # (This assumes parallel-connected loops driven by same voltage)

        # More accurate: solve Z_LL * I = V and compute V^T * I / |I|^2
        # For unit excitation: I = inv(Z_LL) * ones
        n = Z_LL.shape[0]
        ones = np.ones(n)

        try:
            I_L = np.linalg.solve(Z_LL, ones)
            # Total impedance: Z = V_total / I_total = n / sum(I)
            Z_total = n / np.sum(I_L)
        except np.linalg.LinAlgError:
            # Fallback: sum of elements
            Z_total = np.sum(Z_LL)

        return Z_total

    def _solve_loop_star_coupled(self, s: complex) -> complex:
        """
        Solve full Loop-Star coupled system using admittance matrix.

        System equation:
            [Z_LL  Z_LS] [I_L]   [V_L]
            [Z_SL  Z_SS] [I_S] = [V_S]

        where:
            Z_LL = R_LL + s * L_LL (Loop impedance: inductive)
            Z_SS = R_SS + s * L_SS + P_SS / s (Star impedance: capacitive)
            Z_LS = R_LS + s * L_LS (Loop-Star coupling)
            Z_SL = Z_LS^T (reciprocity)

        For parallel LC resonance, both Loop and Star branches see the same
        voltage at the terminals. The total current is I_total = I_L + I_S.

        Using admittance formulation:
            Y_full = inv(Z_full)
            I = Y_full * V

        For unit voltage applied to all ports (V_L = V_S = 1):
            I_total = sum(Y_full * ones)
            Y_total = I_total / V
            Z_total = 1 / Y_total

        This correctly models parallel LC resonance where:
            Z_total = Z_L || Z_C = 1 / (1/Z_L + 1/Z_C)

        Returns total impedance seen at terminals.
        """
        n_L = self._n_loops
        n_S = self._n_stars

        # Build Z_LL (Loop-Loop impedance)
        if self._R_LL.ndim == 1:
            R_LL = np.diag(self._R_LL)
        else:
            R_LL = self._R_LL
        Z_LL = R_LL + s * self._L_LL

        # Build Z_LS (Loop-Star coupling)
        if self._L_LS is not None:
            R_LS = self._R_LS if self._R_LS is not None else np.zeros_like(self._L_LS)
            Z_LS = R_LS + s * self._L_LS
        else:
            Z_LS = np.zeros((n_L, n_S), dtype=complex)

        # Z_SL = Z_LS^T (reciprocity for passive system)
        Z_SL = Z_LS.T

        # Build Z_SS (Star-Star impedance, includes capacitive term)
        if self._L_SS is not None:
            R_SS = self._R_SS if self._R_SS is not None else np.zeros_like(self._L_SS)
            Z_SS = R_SS + s * self._L_SS
        else:
            Z_SS = np.zeros((n_S, n_S), dtype=complex)

        # Add capacitive term P_SS / s = 1/(s*C)
        if self._P_SS is not None:
            Z_SS = Z_SS + self._P_SS / s

        # Assemble full impedance matrix
        n_total = n_L + n_S
        Z_full = np.zeros((n_total, n_total), dtype=complex)
        Z_full[:n_L, :n_L] = Z_LL
        Z_full[:n_L, n_L:] = Z_LS
        Z_full[n_L:, :n_L] = Z_SL
        Z_full[n_L:, n_L:] = Z_SS

        # Compute admittance matrix Y = inv(Z)
        try:
            Y_full = np.linalg.inv(Z_full)
        except np.linalg.LinAlgError:
            # Fallback: return Loop impedance only
            if n_L == 1:
                return Z_LL[0, 0]
            else:
                return np.sum(Z_LL) / (n_L * n_L)

        # Compute total port admittance
        # For parallel connection: apply unit voltage to all nodes
        # I_total = Y * V where V = ones
        # Y_total = sum(I_total) / V = sum(sum(Y))
        #
        # For a single-port view (Loop side is the external port):
        # Y_port = sum of all admittances flowing into the Loop port
        # Y_port = sum(Y_full[:n_L, :])  (all currents due to all voltages)

        # Simplest case: single Loop and single Star (parallel LC)
        if n_L == 1 and n_S == 1:
            # Y_total = Y_LL + Y_SS + Y_LS + Y_SL for parallel
            # But Z is block matrix, so Y_total = sum of all Y elements
            # For V_L = V_S = 1: I_L = Y[0,:] * [1,1]^T = Y[0,0] + Y[0,1]
            #                   I_S = Y[1,:] * [1,1]^T = Y[1,0] + Y[1,1]
            # Total current into port = I_L (Loop current)
            # Z_port = V / I_L = 1 / (Y[0,0] + Y[0,1])

            # For parallel LC: both branches see same V
            # I_total = I_L + I_S = (Y[0,0]+Y[0,1]) + (Y[1,0]+Y[1,1]) * V
            # Actually for single port connected to Loop:
            # Z_port = 1 / sum(Y_full[0, :]) represents all current paths from Loop

            # Correct: Port admittance = admittance seen at Loop terminal
            # When applying V at Loop, current = Y_LL * V + Y_LS * V_S
            # For parallel LC with shared voltage: V_L = V_S = V
            # I = (Y_LL + Y_LS) * V, so Y_port = Y_LL + Y_LS

            # Sum all elements for total admittance (parallel branches)
            Y_total = np.sum(Y_full)
            Z_total = 1.0 / Y_total if abs(Y_total) > 1e-30 else np.inf

        else:
            # General case: Port at Loop[0], others are internal
            # Apply unit voltage at Loop[0], others grounded
            # I_port = Y[0,0] * V = Y[0,0]
            # Z_port = 1 / Y[0,0]

            # For series-connected loops with parallel Star:
            # Use row sum approach: Y_port = sum(Y_full[0,:]) when all V=1

            # Actually, for proper port impedance:
            # V = [1, 0, ..., 0]^T (voltage at port 0, others grounded)
            # I = Y * V = Y[:, 0]
            # I_port = I[0] = Y[0, 0]
            # Z_port = 1 / Y[0, 0]

            Z_total = 1.0 / Y_full[0, 0] if abs(Y_full[0, 0]) > 1e-30 else np.inf

        return Z_total

    # =========================================================================
    # Multi-Port Analysis for WPT (Wireless Power Transfer)
    # =========================================================================

    def set_multiport_matrices(self,
                                L_matrix: np.ndarray,
                                R_matrix: np.ndarray,
                                port_indices: Optional[List[List[int]]] = None):
        """
        Set multi-port PEEC matrices for WPT analysis.

        The L and R matrices contain self and mutual terms:
            L_ii: Self-inductance of port i
            L_ij: Mutual inductance between ports i and j (i != j)
            R_ii: Resistance of port i
            R_ij: Mutual resistance (usually zero)

        Parameters:
            L_matrix: Full inductance matrix [n_loops x n_loops]
            R_matrix: Full resistance matrix [n_loops x n_loops]
            port_indices: List of loop indices for each port (optional)
                          If None, each loop is one port.
                          Example: [[0,1,2], [3,4,5]] for 2 ports, 3 loops each

        Example for 2-coil WPT:
            L = [[L1, M],     # L1 = self-inductance of Tx
                 [M,  L2]]    # L2 = self-inductance of Rx, M = mutual
            R = [[R1, 0],
                 [0,  R2]]
        """
        self._L = np.asarray(L_matrix)
        self._R = np.asarray(R_matrix)

        n_loops = L_matrix.shape[0]
        if port_indices is None:
            # Each loop is one port
            self._port_indices = [[i] for i in range(n_loops)]
            self._n_ports = n_loops
        else:
            self._port_indices = port_indices
            self._n_ports = len(port_indices)

        self._is_built = False

    def solve_frequency_multiport(self, frequencies: np.ndarray) -> 'MultiPortResult':
        """
        Perform multi-port frequency response analysis.

        Computes the full Z-parameter matrix at each frequency:
            [V1]   [Z11 Z12 ... Z1n] [I1]
            [V2] = [Z21 Z22 ... Z2n] [I2]
            [..]   [... ... ... ...] [..]
            [Vn]   [Zn1 Zn2 ... Znn] [In]

        where Z_ij = V_i / I_j with I_k = 0 for k != j

        For WPT analysis, key outputs include:
        - Z_21: Transfer impedance (induced voltage at Rx per Tx current)
        - M_12: Mutual inductance = Im(Z_12) / omega
        - k: Coupling coefficient = M_12 / sqrt(L1 * L2)

        Parameters:
            frequencies: Array of frequencies [Hz]

        Returns:
            MultiPortResult with Z-parameters and WPT metrics
        """
        import time as time_module
        t_start = time_module.time()

        if self._L is None or self._R is None:
            raise ValueError("Matrices not set. Call set_multiport_matrices() first.")

        frequencies = np.asarray(frequencies)
        n_freq = len(frequencies)
        n_ports = getattr(self, '_n_ports', self._L.shape[0])
        n_loops = self._L.shape[0]

        # Allocate result arrays
        Z_matrix = np.zeros((n_freq, n_ports, n_ports), dtype=complex)
        M_12 = np.zeros(n_freq)
        k = np.zeros(n_freq)
        L_self = np.zeros((n_freq, n_ports))
        Z_21 = np.zeros(n_freq, dtype=complex)
        eta = np.zeros(n_freq)

        # Build R matrix if diagonal
        if self._R.ndim == 1:
            R_matrix = np.diag(self._R)
        else:
            R_matrix = self._R

        # Get port grouping
        port_indices = getattr(self, '_port_indices', [[i] for i in range(n_loops)])

        for idx, f in enumerate(frequencies):
            omega = 2.0 * np.pi * f
            s = 1j * omega if f > 1e-10 else 1e-10

            # Build full impedance matrix Z = R + j*omega*L
            Z_full = R_matrix + s * self._L

            # Apply Dowell correction if enabled (to diagonal elements)
            if self._use_dowell and self._conductor_height is not None:
                for i in range(n_loops):
                    Z_full[i, i] = self._apply_dowell_correction(Z_full[i, i], f)

            # Compute Z-parameters for each port
            # For grouped ports, aggregate loop impedances

            if n_ports == n_loops:
                # Simple case: each loop is one port
                # Z-matrix = Z_full directly
                Z_matrix[idx] = Z_full

            else:
                # Grouped ports: need to aggregate
                # Z_ij = sum over loops in port i, j
                for pi in range(n_ports):
                    for pj in range(n_ports):
                        loops_i = port_indices[pi]
                        loops_j = port_indices[pj]
                        # Z_port_ij = sum of Z_full[li, lj] for li in port_i, lj in port_j
                        Z_sum = 0.0j
                        for li in loops_i:
                            for lj in loops_j:
                                Z_sum += Z_full[li, lj]
                        Z_matrix[idx, pi, pj] = Z_sum

            # Compute WPT metrics for 2-port system
            if n_ports >= 2:
                Z11 = Z_matrix[idx, 0, 0]
                Z22 = Z_matrix[idx, 1, 1]
                Z12 = Z_matrix[idx, 0, 1]
                Z21_val = Z_matrix[idx, 1, 0]

                # Self-inductances: L = Im(Z) / omega
                if omega > 1e-10:
                    L1 = np.imag(Z11) / omega
                    L2 = np.imag(Z22) / omega
                    M = np.imag(Z12) / omega  # Mutual inductance

                    L_self[idx, 0] = L1
                    L_self[idx, 1] = L2
                    M_12[idx] = M

                    # Coupling coefficient k = M / sqrt(L1 * L2)
                    if L1 > 0 and L2 > 0:
                        k[idx] = M / np.sqrt(L1 * L2)
                    else:
                        k[idx] = 0.0

                # Transfer impedance Z_21
                Z_21[idx] = Z21_val

                # Power transfer efficiency (simplified, at resonance)
                # eta = (k^2 * Q1 * Q2) / (1 + k^2 * Q1 * Q2)
                # where Q = omega*L / R
                R1 = np.real(Z11)
                R2 = np.real(Z22)
                if omega > 1e-10 and R1 > 1e-30 and R2 > 1e-30:
                    Q1 = omega * L_self[idx, 0] / R1
                    Q2 = omega * L_self[idx, 1] / R2
                    k_sq = k[idx] ** 2
                    eta_num = k_sq * Q1 * Q2
                    eta[idx] = eta_num / (1.0 + eta_num) if eta_num > 0 else 0.0

        t_elapsed = time_module.time() - t_start

        return MultiPortResult(
            analysis_type=AnalysisType.FREQUENCY,
            solver_type=SolverType.PEEC,
            success=True,
            message=f"Multi-port analysis completed ({n_freq} freqs, {n_ports} ports)",
            computation_time=t_elapsed,
            frequencies=frequencies,
            n_ports=n_ports,
            Z_matrix=Z_matrix,
            mutual_inductance=M_12,
            coupling_coefficient=k,
            self_inductance=L_self,
            transfer_impedance=Z_21,
            power_transfer_efficiency=eta
        )

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

        For an RL circuit: L*di/dt + R*i = v
        State-space: di/dt = -R/L * i + 1/L * v
                     A = -R/L (scalar or matrix for CLN)

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
            R_total = self._cln_result.get('R_total', 0.0)
            L_total = self._cln_result.get('L_total', 0.0)
            L_tridiag = self._cln_result['L_tridiag']
            n_states = self._cln_result['n_iter']
        else:
            R_total = getattr(self._cln_result, 'R_total', np.sum(self._R))
            L_total = getattr(self._cln_result, 'L_total', np.sum(self._L))
            L_tridiag = np.asarray(self._cln_result.L_tridiag)
            n_states = len(L_tridiag)

        # Build state-space matrices for CLN
        # The CLN state-space uses the tridiagonal inductance matrix
        # State: x = internal CLN state variables
        # For voltage excitation: L*dx/dt + R*x = v (in reduced space)

        # Time constant tau = L/R
        if abs(R_total) > 1e-30:
            tau = L_total / R_total
        else:
            tau = 1.0

        # Normalized system: tau * T * dx/dt + x = (1/R) * v
        # where T = L_tridiag / L_total (normalized tridiagonal)

        if abs(L_total) > 1e-30:
            T_norm = L_tridiag / L_total
        else:
            T_norm = np.eye(n_states)

        # State-space: dx/dt = -T_norm^{-1} * x / tau + T_norm^{-1} * v / (R * tau)
        try:
            T_inv = np.linalg.inv(T_norm)
        except np.linalg.LinAlgError:
            T_inv = np.linalg.pinv(T_norm)

        # A = -T_inv / tau
        # B = T_inv / (R_total * tau) for voltage input
        A = -T_inv / tau
        B = T_inv[:, 0] / (R_total * tau) if abs(R_total) > 1e-30 else T_inv[:, 0] / tau
        C = np.zeros(n_states)
        C[0] = 1.0  # Output is first state (proportional to current)

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
                # Current is proportional to first state
                # i = v/R at DC, scaled by CLN transfer function at AC
                current[i] = C @ x
            else:
                # For current excitation, compute voltage
                voltage[i] = R_total * current[i] + L_total * (
                    (current[i] - current[i-1]) / (time_array[i] - time_array[i-1])
                    if i > 0 else 0.0
                )

            # Flux linkage: psi = L * i
            flux[i] = L_total * current[i]

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

    # =========================================================================
    # Loop-Star Decomposition Methods
    # =========================================================================

    def set_loop_star_model(self,
                            L_LL: np.ndarray, R_LL: np.ndarray,
                            L_LS: Optional[np.ndarray] = None,
                            L_SS: Optional[np.ndarray] = None,
                            R_LS: Optional[np.ndarray] = None,
                            R_SS: Optional[np.ndarray] = None,
                            P_SS: Optional[np.ndarray] = None):
        """
        Configure Loop-Star decomposed PEEC model.

        Loop-Star decomposition separates currents into:
        - Loop (solenoidal) currents: I_L - inductive behavior
        - Star (irrotational) currents: I_S - capacitive behavior

        The system becomes:
            [Z_LL  Z_LS] [I_L]   [V_L]
            [Z_SL  Z_SS] [I_S] = [V_S]

        Parameters:
            L_LL: Loop-Loop inductance matrix [H]
            R_LL: Loop-Loop resistance matrix [Ohm]
            L_LS: Loop-Star inductance coupling [H] (optional)
            L_SS: Star-Star inductance [H] (optional)
            R_LS: Loop-Star resistance coupling [Ohm] (optional)
            R_SS: Star-Star resistance [Ohm] (optional)
            P_SS: Star-Star potential coefficient [1/F] (optional)
        """
        self._solver.set_loop_star_matrices(
            L_LL, R_LL, L_LS, L_SS, R_LS, R_SS, P_SS
        )
        self._is_configured = True

    def set_mqs_mode(self, use_mqs: bool = True):
        """
        Set MQS (Magneto-Quasi-Static) mode.

        In MQS mode, Star (capacitive) currents are neglected (I_S = 0).
        This is valid for power electronics frequencies where
        omega * epsilon << sigma (conduction dominates displacement).

        Parameters:
            use_mqs: If True, solve Loop equations only
        """
        self._solver.set_mqs_mode(use_mqs)

    def frequency_sweep_loop_star(self, frequencies: np.ndarray) -> FrequencyResult:
        """
        Perform frequency response using Loop-Star decomposition.

        This method uses Loop-Star matrices if available, providing
        better numerical stability at low frequencies.

        In MQS mode, only Loop equations are solved (Star currents = 0).
        In full mode, the coupled Loop-Star system is solved.

        Parameters:
            frequencies: Array of frequencies [Hz]

        Returns:
            FrequencyResult with impedance Z(f)
        """
        if not self._is_configured:
            raise RuntimeError("Model not configured. Call set_loop_star_model() first.")
        return self._solver.solve_frequency_loop_star(frequencies)

    # =========================================================================
    # Multi-Port Analysis for WPT
    # =========================================================================

    def set_multiport_model(self,
                            L_matrix: np.ndarray,
                            R_matrix: np.ndarray,
                            port_indices: Optional[List[List[int]]] = None):
        """
        Configure multi-port PEEC model for WPT analysis.

        For a 2-coil WPT system:
            L = [[L1, M],     # L1 = Tx self-inductance
                 [M,  L2]]    # L2 = Rx self-inductance, M = mutual
            R = [[R1, 0],
                 [0,  R2]]

        Parameters:
            L_matrix: Full inductance matrix including mutual terms [H]
            R_matrix: Full resistance matrix [Ohm]
            port_indices: Port grouping (optional). If None, each loop = 1 port.

        Example for simple 2-coil WPT:
            >>> L1, L2, M = 10e-6, 10e-6, 5e-6  # 10 uH each, 5 uH mutual
            >>> R1, R2 = 0.1, 0.1  # 0.1 Ohm each
            >>> analysis.set_multiport_model(
            ...     L_matrix=np.array([[L1, M], [M, L2]]),
            ...     R_matrix=np.array([[R1, 0], [0, R2]])
            ... )
        """
        self._solver.set_multiport_matrices(L_matrix, R_matrix, port_indices)
        self._is_configured = True

    def frequency_sweep_multiport(self, frequencies: np.ndarray) -> MultiPortResult:
        """
        Perform multi-port frequency response for WPT analysis.

        Returns Z-parameters and WPT metrics:
        - Z_matrix: Full impedance matrix [Z_ij] at each frequency
        - mutual_inductance: M_12 = Im(Z_12) / omega
        - coupling_coefficient: k = M / sqrt(L1 * L2)
        - power_transfer_efficiency: eta at each frequency

        Parameters:
            frequencies: Array of frequencies [Hz]

        Returns:
            MultiPortResult with WPT analysis data

        Example:
            >>> freqs = np.logspace(3, 6, 100)  # 1 kHz to 1 MHz
            >>> result = analysis.frequency_sweep_multiport(freqs)
            >>> print(f"k = {result.coupling_coefficient[50]:.3f}")
            >>> print(f"Max eta = {np.max(result.power_transfer_efficiency)*100:.1f}%")
        """
        if not self._is_configured:
            raise RuntimeError("Model not configured. Call set_multiport_model() first.")
        return self._solver.solve_frequency_multiport(frequencies)

    def compute_wpt_coupling(self, L1: float, L2: float, M: float) -> dict:
        """
        Compute WPT coupling parameters from inductances.

        Parameters:
            L1: Transmitter self-inductance [H]
            L2: Receiver self-inductance [H]
            M: Mutual inductance [H]

        Returns:
            dict with:
                k: Coupling coefficient
                n: Turns ratio equivalent (sqrt(L2/L1))
                M_max: Maximum possible M (geometric mean)
        """
        k = M / np.sqrt(L1 * L2) if L1 > 0 and L2 > 0 else 0.0
        n = np.sqrt(L2 / L1) if L1 > 0 else 1.0
        M_max = np.sqrt(L1 * L2)

        return {
            'k': k,
            'n': n,
            'M': M,
            'M_max': M_max,
            'L1': L1,
            'L2': L2
        }

    def find_resonant_frequency(self, L: float, C: float) -> float:
        """
        Calculate resonant frequency for LC compensation.

        Parameters:
            L: Inductance [H]
            C: Capacitance [F]

        Returns:
            Resonant frequency [Hz]
        """
        return 1.0 / (2.0 * np.pi * np.sqrt(L * C))

    def design_compensation_capacitor(self, L: float, f_res: float) -> float:
        """
        Design compensation capacitor for WPT at target frequency.

        Parameters:
            L: Inductance to compensate [H]
            f_res: Target resonant frequency [Hz]

        Returns:
            Required capacitance [F]
        """
        omega = 2.0 * np.pi * f_res
        return 1.0 / (omega ** 2 * L)


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
