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
        """Get Y-parameters (Y = inv(Z))."""
        n_freq = len(self.frequencies)
        Y = np.zeros_like(self.Z_matrix)
        for i in range(n_freq):
            try:
                Y[i] = np.linalg.inv(self.Z_matrix[i])
            except np.linalg.LinAlgError:
                Y[i] = np.linalg.pinv(self.Z_matrix[i])
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
                warnings.warn(
                    f"Singular Z-matrix at f={self.frequencies[i]:.6g} Hz, "
                    "S-parameters set to zero"
                )
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

        # Material parameters (constant model)
        self._mu_r_static = 1000.0      # DC relative permeability
        self._mu_r_inf = 1.0            # High-frequency relative permeability
        self._tau = 1e-6                # Relaxation time [s]
        self._sigma = 0.0               # Conductivity (for eddy current loss)

        # H-dependent permeability table (optional)
        # Format: [[H1, mu_r1], [H2, mu_r2], ...] where H is in A/m
        self._mu_H_table = None
        self._H_operating = 0.0         # Operating H level [A/m]

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
        Set material parameters for frequency response (constant model).

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
        self._mu_H_table = None  # Clear H-dependent table

    def set_h_dependent_permeability(self,
                                      mu_H_table: np.ndarray,
                                      tau: float = 1e-6,
                                      sigma: float = 0.0):
        """
        Set H-dependent permeability for nonlinear frequency response.

        For magnetic materials, the permeability decreases with increasing
        H field due to saturation. This method allows specifying a
        mu_r(H) curve that will be used for frequency response analysis.

        The frequency response uses the Debye model with H-dependent mu:
            mu(omega, H) = mu_inf + (mu(H) - mu_inf) / (1 + j*omega*tau)

        where mu(H) is interpolated from the table.

        Parameters:
            mu_H_table: Array of [H, mu_r] pairs, e.g.
                        [[0, 5000], [100, 4500], [1000, 2000], [10000, 100]]
                        H in A/m, mu_r is relative permeability
            tau: Relaxation time constant [s]
            sigma: Conductivity for eddy current loss [S/m]

        Example:
            # Typical soft iron saturation curve
            mu_H_data = [
                [0, 5000],      # Initial permeability at low H
                [100, 4500],
                [500, 3000],
                [1000, 2000],
                [5000, 500],
                [10000, 100],   # Deep saturation
            ]
            solver.set_h_dependent_permeability(mu_H_data, tau=1e-6)
        """
        self._mu_H_table = np.asarray(mu_H_table)
        self._tau = tau
        self._sigma = sigma

        # Set static mu_r from first point (H=0 or minimum H)
        if len(self._mu_H_table) > 0:
            # Sort by H and get mu at minimum H
            sorted_table = self._mu_H_table[np.argsort(self._mu_H_table[:, 0])]
            self._mu_r_static = sorted_table[0, 1]
            # High-frequency permeability (assumed 1.0 for ferromagnetic)
            self._mu_r_inf = 1.0

    def set_operating_field(self, H: float):
        """
        Set the operating H field level for H-dependent analysis.

        For frequency response with H-dependent permeability, this sets
        the bias point around which small-signal behavior is computed.

        Parameters:
            H: Operating H field magnitude [A/m]
        """
        self._H_operating = abs(H)

    def _interpolate_mu_r(self, H: float) -> float:
        """
        Interpolate relative permeability from H-dependent table.

        Parameters:
            H: Magnetic field strength [A/m]

        Returns:
            Relative permeability mu_r at given H
        """
        if self._mu_H_table is None or len(self._mu_H_table) == 0:
            return self._mu_r_static

        H_vals = self._mu_H_table[:, 0]
        mu_vals = self._mu_H_table[:, 1]

        # Sort by H
        sort_idx = np.argsort(H_vals)
        H_sorted = H_vals[sort_idx]
        mu_sorted = mu_vals[sort_idx]

        # Interpolate (linear)
        if H <= H_sorted[0]:
            return mu_sorted[0]
        elif H >= H_sorted[-1]:
            return mu_sorted[-1]
        else:
            return np.interp(H, H_sorted, mu_sorted)

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

    def solve_frequency(self, frequencies: np.ndarray,
                        H_levels: Optional[np.ndarray] = None) -> MMMFrequencyResult:
        """
        Perform frequency response analysis.

        Computes the complex permeability mu(omega) using the Debye model:
            mu(omega) = mu_inf + (mu_s - mu_inf) / (1 + j*omega*tau)

        For H-dependent permeability (if set_h_dependent_permeability was called):
            mu(omega, H) = mu_inf + (mu(H) - mu_inf) / (1 + j*omega*tau)

        Plus eddy current contribution (if sigma > 0):
            mu_eddy(omega) = mu_0 * sigma / (j*omega)  (added to imaginary part)

        Parameters:
            frequencies: Array of frequencies [Hz]
            H_levels: Optional array of H field levels [A/m] for H-dependent analysis.
                      If None, uses self._H_operating (set by set_operating_field()).
                      If a single value, uses that for all frequencies.
                      If same length as frequencies, uses corresponding H for each f.

        Returns:
            MMMFrequencyResult with complex permeability data
        """
        import time
        t_start = time.time()

        frequencies = np.asarray(frequencies)
        n_freq = len(frequencies)

        # Complex permeability array
        mu_complex = np.zeros(n_freq, dtype=complex)

        # Determine H level(s) for analysis
        if H_levels is not None:
            H_levels = np.asarray(H_levels)
            if H_levels.size == 1:
                H_array = np.full(n_freq, float(H_levels))
            else:
                H_array = H_levels
        else:
            H_array = np.full(n_freq, self._H_operating)

        # High-frequency permeability
        mu_inf = self._mu_r_inf * MU_0
        tau = self._tau

        for i, f in enumerate(frequencies):
            omega = 2.0 * np.pi * f

            # Get mu_r at operating H level
            if self._mu_H_table is not None:
                mu_r_H = self._interpolate_mu_r(H_array[i])
            else:
                mu_r_H = self._mu_r_static

            mu_s = mu_r_H * MU_0
            delta_mu = mu_s - mu_inf

            if omega < 1e-10:
                # DC limit
                mu_complex[i] = mu_s
            else:
                # Debye model with H-dependent static permeability
                mu_debye = mu_inf + delta_mu / (1.0 + 1j * omega * tau)

                # Add eddy current loss (if conductive)
                if self._sigma > 0:
                    # Skin depth: delta = sqrt(2 / (omega * mu * sigma))
                    # Use local mu for skin depth calculation
                    mu_eff = np.real(mu_debye)
                    if mu_eff < MU_0:
                        mu_eff = MU_0
                    skin_depth = np.sqrt(2.0 / (omega * mu_eff * self._sigma))
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

        # Build message
        if self._mu_H_table is not None:
            msg = f"MMM frequency sweep completed ({n_freq} points, H-dependent mu)"
        else:
            msg = f"MMM frequency sweep completed ({n_freq} points)"

        return MMMFrequencyResult(
            analysis_type=AnalysisType.FREQUENCY,
            solver_type=SolverType.MMM,
            success=True,
            message=msg,
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

    def set_h_dependent_permeability(self,
                                      mu_H_table: np.ndarray,
                                      tau: float = 1e-6,
                                      sigma: float = 0.0):
        """
        Set H-dependent permeability for nonlinear frequency response.

        For magnetic materials, the permeability decreases with increasing
        H field due to saturation.

        Parameters:
            mu_H_table: Array of [H, mu_r] pairs, e.g.
                        [[0, 5000], [100, 4500], [1000, 2000], [10000, 100]]
                        H in A/m, mu_r is relative permeability
            tau: Relaxation time constant [s]
            sigma: Conductivity for eddy current loss [S/m]

        Example:
            >>> analysis = UnifiedMMMAnalysis()
            >>> analysis.set_radia_model(container, mu_r=5000)
            >>>
            >>> # Set saturation curve
            >>> mu_H_data = [
            ...     [0, 5000],
            ...     [1000, 2000],
            ...     [10000, 100],
            ... ]
            >>> analysis.set_h_dependent_permeability(mu_H_data)
            >>>
            >>> # Frequency sweep at different operating points
            >>> result_low_H = analysis.frequency_sweep(freqs, H_operating=100)
            >>> result_high_H = analysis.frequency_sweep(freqs, H_operating=5000)
        """
        self._solver.set_h_dependent_permeability(mu_H_table, tau, sigma)

    def set_operating_field(self, H: float):
        """
        Set the operating H field level for H-dependent analysis.

        Parameters:
            H: Operating H field magnitude [A/m]
        """
        self._solver.set_operating_field(H)

    def static(self) -> MMMStaticResult:
        """Perform static magnetization analysis."""
        if not self._is_configured:
            raise RuntimeError("Model not configured. Call set_radia_model() first.")
        return self._solver.solve_static()

    def frequency_sweep(self, frequencies: np.ndarray,
                        H_operating: Optional[float] = None) -> MMMFrequencyResult:
        """
        Perform frequency response analysis.

        Computes complex permeability mu(omega) using Debye relaxation model.
        If H-dependent permeability is set, uses the operating H level.

        Parameters:
            frequencies: Array of frequencies [Hz]
            H_operating: Operating H field [A/m] for H-dependent analysis.
                         If None, uses the value set by set_operating_field().

        Returns:
            MMMFrequencyResult with mu', mu", loss tangent
        """
        if not self._is_configured:
            raise RuntimeError("Model not configured. Call set_radia_model() first.")

        if H_operating is not None:
            self._solver.set_operating_field(H_operating)

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
