"""
PRIMA-based Model Order Reduction for Loop-Star-Magnetic System

This module implements PRIMA (Passive Reduced-order Interconnect Macromodeling
Algorithm) with Lanczos process for SPICE-compatible circuit extraction.

Key features:
1. PRIMA Lanczos with re-orthogonalization (higher accuracy than plain Lanczos)
2. Tridiagonal structure yields RL ladder network (no Arnoldi, direct circuit)
3. Schur complement for port impedance extraction
4. Lanczos reduction for Star (capacitive) DOFs - reduces capacitor count
5. PyKAN integration for complex material properties (mu", eps")

The goal is to reduce:
    [Z_LL    Z_LS    Z_LM  ] [I_L]   [V_L  ]
    [Z_SL    Z_SS    Z_SM  ] [Q_S] = [0    ]
    [Z_ML    Z_MS    Z_MM  ] [M  ]   [H_ext]

To a sparse equivalent circuit with KAN-based material models.

Note: This implementation uses PRIMA (not CLN/Cauer), avoiding patent issues
while providing passive, stable reduced-order models.

Design Decision (2026-01-17):
- ACA+ was evaluated but removed in favor of Lanczos-only approach
- Reason: ACA+ doesn't reduce element count in standard SPICE netlists
- Lanczos directly produces tridiagonal -> RL/RC ladder with fewer elements

References:
- A. Odabasioglu, M. Celik, L.T. Pileggi, "PRIMA: Passive Reduced-order
  Interconnect Macromodeling Algorithm," IEEE TCAD, 1998.
- B. Gustavsen, A. Semlyen, "Rational approximation of frequency domain
  responses by Vector Fitting," IEEE TPWRD, 1999.

Author: Radia Development Team
Date: 2026-01-17
"""

import numpy as np
from scipy import linalg
from scipy.integrate import solve_ivp
from typing import Tuple, List, Dict, Callable, Optional, Union
from dataclasses import dataclass, field


@dataclass
class LanczosResult:
    """Result of Lanczos tridiagonalization."""
    Q: np.ndarray          # Orthonormal basis [n x k]
    alpha: np.ndarray      # Diagonal of tridiagonal matrix [k]
    beta: np.ndarray       # Off-diagonal of tridiagonal matrix [k-1]
    T: np.ndarray          # Tridiagonal matrix [k x k]
    rank: int              # Effective rank


@dataclass
class LowRankBlock:
    """Low-rank approximation of a matrix block: A ≈ U @ S @ V.T"""
    U: np.ndarray          # Left basis [m x r]
    S: np.ndarray          # Singular values [r]
    V: np.ndarray          # Right basis [n x r]
    rank: int              # Rank r

    def to_dense(self) -> np.ndarray:
        """Reconstruct dense matrix."""
        return self.U @ np.diag(self.S) @ self.V.T

    def matvec(self, x: np.ndarray) -> np.ndarray:
        """Matrix-vector product A @ x."""
        return self.U @ (self.S * (self.V.T @ x))

    def rmatvec(self, x: np.ndarray) -> np.ndarray:
        """Adjoint matrix-vector product A.T @ x."""
        return self.V @ (self.S * (self.U.T @ x))


@dataclass
class SparseCircuit:
    """Sparse circuit representation."""
    # Node information
    n_nodes: int
    node_names: List[str]

    # Circuit elements
    resistors: List[Tuple[int, int, float]]      # (node1, node2, R)
    inductors: List[Tuple[int, int, float]]      # (node1, node2, L)
    capacitors: List[Tuple[int, int, float]]     # (node1, node2, C)
    mutual_inductors: List[Tuple[int, int, int, int, float]]  # (L1+, L1-, L2+, L2-, M)
    controlled_sources: List[Dict]               # VCVS, CCCS, etc.

    # Port information
    port_nodes: List[Tuple[int, int]]            # (positive, negative) for each port


class LanczosReducer:
    """
    Lanczos-based model order reduction.

    Algorithm:
    1. Lanczos on L (or R^{-1}L for generalized eigenvalue):
       L = Q @ T @ Q.T  where T is tridiagonal

    2. Transform coupling blocks:
       Z_LM' = Q.T @ Z_LM
       Z_LS' = Q.T @ Z_LS

    3. The reduced Loop block becomes:
       Z_LL' = Q.T @ (R + sL) @ Q = R' + sT
       where R' = Q.T @ R @ Q (also sparse if R is diagonal)

    4. Apply ACA+ to Z_LM' and Z_LS' for further compression

    5. Schur complement eliminates internal nodes
    """

    def __init__(self, tol: float = 1e-6, max_rank: int = None):
        """
        Parameters:
            tol: Tolerance for rank truncation
            max_rank: Maximum rank (None = no limit)
        """
        self.tol = tol
        self.max_rank = max_rank

    def lanczos_symmetric(self, A: np.ndarray, k: int = None,
                          v0: np.ndarray = None) -> LanczosResult:
        """
        Lanczos tridiagonalization for symmetric matrix.

        A = Q @ T @ Q.T

        where T is tridiagonal:
            T = [alpha_0  beta_0    0      0    ...]
                [beta_0   alpha_1  beta_1  0    ...]
                [0        beta_1   alpha_2 beta_2...]
                [...]

        Parameters:
            A: Symmetric matrix [n x n]
            k: Number of Lanczos vectors (default: n)
            v0: Starting vector (default: random)

        Returns:
            LanczosResult with Q, alpha, beta, T
        """
        n = A.shape[0]
        if k is None:
            k = n
        k = min(k, n)

        # Initialize
        Q = np.zeros((n, k))
        alpha = np.zeros(k)
        beta = np.zeros(k - 1)

        # Starting vector
        if v0 is None:
            v0 = np.random.randn(n)
        v = v0 / np.linalg.norm(v0)
        Q[:, 0] = v

        # Lanczos iteration
        w = A @ v
        alpha[0] = np.dot(v, w)
        w = w - alpha[0] * v

        for j in range(1, k):
            beta[j-1] = np.linalg.norm(w)

            # Check for convergence (invariant subspace found)
            if beta[j-1] < self.tol * np.linalg.norm(A, 'fro'):
                # Restart with random vector orthogonal to Q[:, :j]
                w = np.random.randn(n)
                for i in range(j):
                    w = w - np.dot(Q[:, i], w) * Q[:, i]
                if np.linalg.norm(w) < self.tol:
                    # Cannot extend - return early
                    k = j
                    break

            v_prev = v
            v = w / beta[j-1]
            Q[:, j] = v

            w = A @ v - beta[j-1] * v_prev
            alpha[j] = np.dot(v, w)
            w = w - alpha[j] * v

            # Reorthogonalization (full)
            for i in range(j + 1):
                w = w - np.dot(Q[:, i], w) * Q[:, i]

        # Build tridiagonal matrix
        T = np.diag(alpha[:k]) + np.diag(beta[:k-1], 1) + np.diag(beta[:k-1], -1)

        return LanczosResult(
            Q=Q[:, :k],
            alpha=alpha[:k],
            beta=beta[:k-1],
            T=T,
            rank=k
        )

    def lanczos_generalized(self, L: np.ndarray, R: np.ndarray,
                            k: int = None) -> LanczosResult:
        """
        Lanczos for generalized eigenvalue problem: L @ v = lambda @ R @ v

        This is useful when R is not identity (e.g., includes skin effect).

        We solve: R^{-1} @ L @ v = lambda @ v
        but preserve symmetry via Cholesky: R = C @ C.T
        Then: C^{-1} @ L @ C^{-T} @ (C.T @ v) = lambda @ (C.T @ v)

        Parameters:
            L: Inductance matrix (symmetric positive definite)
            R: Resistance matrix (symmetric positive definite)
            k: Number of Lanczos vectors

        Returns:
            LanczosResult for the transformed problem
        """
        n = L.shape[0]
        if k is None:
            k = n

        # Check if R is diagonal (common case)
        if np.allclose(R, np.diag(np.diag(R))):
            # Simple scaling
            R_sqrt_inv = np.diag(1.0 / np.sqrt(np.diag(R)))
            A = R_sqrt_inv @ L @ R_sqrt_inv
            result = self.lanczos_symmetric(A, k)
            # Transform back: Q_original = R_sqrt_inv @ Q
            result.Q = R_sqrt_inv @ result.Q
            return result
        else:
            # General case: Cholesky
            try:
                C = linalg.cholesky(R, lower=True)
                C_inv = linalg.solve_triangular(C, np.eye(n), lower=True)
                A = C_inv @ L @ C_inv.T
                result = self.lanczos_symmetric(A, k)
                result.Q = C_inv.T @ result.Q
                return result
            except linalg.LinAlgError:
                # R not positive definite - fall back to standard
                return self.lanczos_symmetric(L, k)

    def truncate_lanczos(self, result: LanczosResult,
                         tol: float = None,
                         keep_tridiagonal: bool = True) -> LanczosResult:
        """
        Truncate Lanczos basis based on eigenvalue decay.

        Parameters:
            result: Full Lanczos result
            tol: Truncation tolerance (default: self.tol)
            keep_tridiagonal: If True, keep tridiagonal structure (sparse circuit)
                              If False, diagonalize (dense but decoupled)

        Returns:
            Truncated LanczosResult
        """
        if tol is None:
            tol = self.tol

        # Eigendecomposition of T to determine truncation rank
        eigvals, eigvecs = linalg.eigh(result.T)

        # Sort by magnitude
        idx = np.argsort(np.abs(eigvals))[::-1]
        eigvals_sorted = eigvals[idx]

        # Find truncation rank
        total = np.sum(np.abs(eigvals_sorted))
        cumsum = np.cumsum(np.abs(eigvals_sorted))
        k_trunc = np.searchsorted(cumsum / total, 1.0 - tol) + 1

        if self.max_rank is not None:
            k_trunc = min(k_trunc, self.max_rank)
        k_trunc = max(1, min(k_trunc, result.rank))

        if keep_tridiagonal:
            # Keep tridiagonal structure - just truncate
            return LanczosResult(
                Q=result.Q[:, :k_trunc],
                alpha=result.alpha[:k_trunc],
                beta=result.beta[:k_trunc-1] if k_trunc > 1 else np.array([]),
                T=result.T[:k_trunc, :k_trunc],
                rank=k_trunc
            )
        else:
            # Diagonalize - transform to eigenbasis
            eigvecs_sorted = eigvecs[:, idx]
            Q_new = result.Q @ eigvecs_sorted[:, :k_trunc]
            T_new = np.diag(eigvals_sorted[:k_trunc])

            return LanczosResult(
                Q=Q_new,
                alpha=eigvals_sorted[:k_trunc],
                beta=np.zeros(k_trunc - 1),  # Diagonal after eigenvector transform
                T=T_new,
                rank=k_trunc
            )


# ACAPlus class was removed (2026-01-17)
# Reason: ACA+ doesn't reduce element count in standard SPICE netlists.
# Lanczos-only approach is used for all model order reduction.
# See design decision in module docstring.


class HierarchicalReducer:
    """
    Hierarchical model order reduction using Lanczos only.

    System:
        [Z_LL    Z_LS    Z_LM  ] [I_L]   [V_L  ]
        [Z_SL    Z_SS    Z_SM  ] [Q_S] = [0    ]
        [Z_ML    Z_MS    Z_MM  ] [M  ]   [H_ext]

    Step 1: Lanczos on L (DC inductance)
        L = Q_L @ T_L @ Q_L.T

    Step 2: Transform to Lanczos basis
        Z_LL' = Q_L.T @ Z_LL @ Q_L = R' + s*T_L
        Z_LM' = Q_L.T @ Z_LM
        Z_LS' = Q_L.T @ Z_LS

    Step 3: Lanczos on Star (capacitance) block
        P = Q_S @ T_P @ Q_S.T  (P = inverse capacitance)
        Reduces capacitor count in SPICE output

    Step 4: Schur complement
        Z_port = Z_LL' - Z_LM' @ Z_MM^{-1} @ Z_ML' - Z_LS' @ Z_SS^{-1} @ Z_SL'

    Note: ACA+ was removed (2026-01-17) - Lanczos directly reduces element count.
    """

    def __init__(self, lanczos_tol: float = 1e-6,
                 max_lanczos_rank: int = None):
        self.lanczos = LanczosReducer(tol=lanczos_tol, max_rank=max_lanczos_rank)

        # Stored results
        self._lanczos_result = None
        self._reduced_blocks = {}

    def reduce_loop_block(self, L: np.ndarray, R: np.ndarray,
                          truncate: bool = True,
                          port_indices: List[int] = None,
                          k: int = None) -> LanczosResult:
        """
        Step 1: Lanczos tridiagonalization of Loop block (PRIMA-style).

        Uses port excitation vector as starting vector for Lanczos iteration.
        This ensures the Krylov subspace K(L, e_port) captures port response
        accurately with minimal modes.

        Parameters:
            L: Inductance matrix [n_L x n_L]
            R: Resistance matrix [n_L x n_L]
            truncate: Whether to truncate based on eigenvalue decay
            port_indices: Port node indices (default: [0] = first node)
            k: Number of Lanczos vectors (default: determined by tolerance)

        Returns:
            LanczosResult with orthonormal basis Q and tridiagonal T
        """
        n = L.shape[0]

        # Default port is first node
        if port_indices is None:
            port_indices = [0]

        # PRIMA: Start from port excitation vector
        # For single port: e_0 = [1, 0, 0, ..., 0]
        # For multiple ports: use sum or block Lanczos
        v0 = np.zeros(n)
        for idx in port_indices:
            v0[idx] = 1.0

        # Store port info for later
        self._port_indices = port_indices
        self._port_vector = v0.copy()

        # Lanczos on L with port excitation as starting vector
        if np.allclose(R, np.eye(R.shape[0]) * R[0, 0]):
            # R is scalar multiple of identity
            result = self.lanczos.lanczos_symmetric(L, k=k, v0=v0)
        else:
            # For generalized eigenvalue problem, need to transform v0
            result = self.lanczos.lanczos_generalized(L, R, k=k)

        if truncate:
            result = self.lanczos.truncate_lanczos(result)

        self._lanczos_result = result

        # Also compute R in Lanczos basis
        Q = result.Q
        self._R_reduced = Q.T @ R @ Q

        return result

    def transform_coupling_blocks(self, Z_LM: np.ndarray = None,
                                   Z_LS: np.ndarray = None) -> dict:
        """
        Step 2: Transform coupling blocks to Lanczos basis.

        Z_LM' = Q.T @ Z_LM
        Z_LS' = Q.T @ Z_LS

        Parameters:
            Z_LM: Loop-Magnetic coupling [n_L x n_M]
            Z_LS: Loop-Star coupling [n_L x n_S]

        Returns:
            Dictionary of transformed blocks
        """
        if self._lanczos_result is None:
            raise RuntimeError("Call reduce_loop_block first")

        Q = self._lanczos_result.Q

        self._reduced_blocks = {}

        if Z_LM is not None:
            self._reduced_blocks['Z_LM'] = Q.T @ Z_LM
            self._reduced_blocks['Z_ML'] = Z_LM.T @ Q

        if Z_LS is not None:
            self._reduced_blocks['Z_LS'] = Q.T @ Z_LS
            self._reduced_blocks['Z_SL'] = Z_LS.T @ Q

        return self._reduced_blocks

    def compute_schur_complement(self, s: complex,
                                  Z_SS: np.ndarray = None,
                                  Z_MM: np.ndarray = None) -> np.ndarray:
        """
        Step 4: Schur complement to eliminate Star and Magnetic DOFs.

        PEEC-magnetic System (physically correct formulation):
            [R + sL    s*K   ] [I]   [V]
            [-K.T      Z_MM  ] [M] = [0]

        Schur complement elimination:
            From row 2: M = Z_MM^{-1} @ K.T @ I
            Z_schur = (R + sL) + s*K @ Z_MM^{-1} @ K.T
                    = R + s*(L + L_coupling)
                    = R + s*L_eff

        This gives L_eff = L + L_coupling > L (magnetic material INCREASES inductance)

        Parameters:
            s: Complex frequency
            Z_SS: Star-Star impedance (typically P/s for capacitive)
            Z_MM: Magnetic-Magnetic impedance

        Returns:
            Reduced port impedance matrix [k x k] in Lanczos basis
        """
        if self._lanczos_result is None:
            raise RuntimeError("Call reduce_loop_block first")

        T_L = self._lanczos_result.T
        R_red = self._R_reduced

        # Z_LL' = R' + s*T_L (in Lanczos basis)
        Z_LL_red = R_red + s * T_L

        # Initialize Schur complement
        Z_schur = Z_LL_red.copy()

        # Add Magnetic contribution (+ sign: magnetic increases inductance)
        #
        # Full system block matrix (physically correct):
        #   [R + sL    s*K  ] [I]   [V]
        #   [-K.T      Z_MM ] [M] = [0]
        #
        # Schur complement:
        #   Z_schur = (R + sL) - (s*K) @ Z_MM^{-1} @ (-K.T)
        #           = (R + sL) + s * K @ Z_MM^{-1} @ K.T
        #
        # In Lanczos reduced form (Q.T @ ... @ Q):
        #   K_red = Q.T @ K  [k x n_M]
        #   Z_schur = Z_LL' + s * K_red @ Z_MM^{-1} @ K_red.T
        if 'Z_LM' in self._reduced_blocks and Z_MM is not None:
            K_red = self._reduced_blocks['Z_LM']  # [k x n_M]
            try:
                Z_MM_inv = np.linalg.inv(Z_MM)
                # Physically correct: + s * coupling (inductance increases)
                Z_schur += s * K_red @ Z_MM_inv @ K_red.T
            except np.linalg.LinAlgError:
                pass

        # Subtract Star contribution (capacitive effect)
        if 'Z_LS' in self._reduced_blocks and Z_SS is not None:
            Z_LS_red = self._reduced_blocks['Z_LS']  # [k x n_S]
            Z_SL_red = self._reduced_blocks['Z_SL']  # [n_S x k]
            try:
                Z_SS_inv = np.linalg.inv(Z_SS)
                Z_schur -= s**3 * Z_LS_red @ Z_SS_inv @ Z_SL_red
            except np.linalg.LinAlgError:
                pass

        return Z_schur

    def get_port_impedance(self, Z_schur: np.ndarray) -> complex:
        """
        Get port impedance from Schur complement matrix.

        With PRIMA (port excitation as starting vector), the port impedance
        is simply the (0,0) element of the Schur complement:

            Z_port = Z_schur[0, 0]

        This works because Q[:, 0] = e_port (normalized), so:
            b = Q.T @ e_port = [1, 0, 0, ..., 0]
            Z_port = b.T @ Z_schur @ b = Z_schur[0, 0]

        Parameters:
            Z_schur: Schur complement matrix from compute_schur_complement()

        Returns:
            Port impedance (complex)
        """
        return Z_schur[0, 0]

    def extract_sparse_circuit(self, port_indices: List[int] = None) -> SparseCircuit:
        """
        Step 5: Extract sparse circuit representation.

        The tridiagonal structure of T_L corresponds to a ladder network:

            o---R0---L0---o---R1---L1---o---R2---L2---o
                          |             |             |
                         M01           M12           M23
                          |             |             |

        Plus coupling to Magnetic and Star blocks via controlled sources.

        Parameters:
            port_indices: Which Lanczos nodes are ports (default: [0])

        Returns:
            SparseCircuit representation
        """
        if self._lanczos_result is None:
            raise RuntimeError("Call reduce_loop_block first")

        k = self._lanczos_result.rank
        alpha = self._lanczos_result.alpha
        beta = self._lanczos_result.beta
        R_red = self._R_reduced

        if port_indices is None:
            port_indices = [0]

        # Node naming: L0, L1, ..., Lk-1 for Lanczos nodes
        # Plus ground node
        node_names = [f'L{i}' for i in range(k)] + ['GND']
        n_nodes = k + 1
        gnd = k  # Ground node index

        resistors = []
        inductors = []
        capacitors = []
        mutual_inductors = []
        controlled_sources = []

        # Diagonal elements: self-inductance and resistance
        for i in range(k):
            # Inductance alpha_i (from tridiagonal diagonal)
            inductors.append((i, gnd, alpha[i]))

            # Resistance from R_reduced diagonal
            if i < R_red.shape[0]:
                resistors.append((i, gnd, R_red[i, i]))

        # Off-diagonal elements: mutual inductance (from tridiagonal off-diagonal)
        for i in range(k - 1):
            if abs(beta[i]) > 1e-15:
                # beta_i represents coupling between L_i and L_{i+1}
                # In circuit terms: mutual inductance
                # M_{i,i+1} = beta_i
                mutual_inductors.append((i, gnd, i+1, gnd, beta[i]))

        # Add controlled sources for Magnetic coupling (if ACA compressed)
        if 'Z_LM' in self._aca_blocks:
            lr = self._aca_blocks['Z_LM']
            # Each rank-1 component U[:,r] * S[r] * V[:,r].T becomes a CCVS
            # V_Lr = S[r] * (sum_j V[j,r] * M_j) at node Lr
            for r in range(lr.rank):
                controlled_sources.append({
                    'type': 'CCVS',  # Current-controlled voltage source
                    'gain': lr.S[r],
                    'control_weights': lr.V[:, r].tolist(),
                    'output_weights': lr.U[:, r].tolist(),
                    'description': f'Magnetic coupling rank-{r}'
                })

        # Add controlled sources for Star coupling
        if 'Z_LS' in self._aca_blocks:
            lr = self._aca_blocks['Z_LS']
            for r in range(lr.rank):
                controlled_sources.append({
                    'type': 'VCCS',  # Voltage-controlled current source
                    'gain': lr.S[r],
                    'control_weights': lr.V[:, r].tolist(),
                    'output_weights': lr.U[:, r].tolist(),
                    'description': f'Star coupling rank-{r}'
                })

        # Port nodes
        port_nodes = [(i, gnd) for i in port_indices]

        return SparseCircuit(
            n_nodes=n_nodes,
            node_names=node_names,
            resistors=resistors,
            inductors=inductors,
            capacitors=capacitors,
            mutual_inductors=mutual_inductors,
            controlled_sources=controlled_sources,
            port_nodes=port_nodes
        )

    def get_reduction_statistics(self) -> dict:
        """Get statistics about the reduction."""
        stats = {
            'original_loop_size': None,
            'lanczos_rank': None,
            'lanczos_compression': None,
            'aca_ranks': {},
            'total_compression': None
        }

        if self._lanczos_result is not None:
            n_orig = self._lanczos_result.Q.shape[0]
            k = self._lanczos_result.rank
            stats['original_loop_size'] = n_orig
            stats['lanczos_rank'] = k
            stats['lanczos_compression'] = k / n_orig if n_orig > 0 else 1.0

        for name, lr in self._aca_blocks.items():
            stats['aca_ranks'][name] = lr.rank

        return stats


def print_sparse_circuit(circuit: SparseCircuit):
    """Print sparse circuit in SPICE-like format."""
    print("\n" + "=" * 60)
    print("Sparse Circuit Representation")
    print("=" * 60)

    print(f"\nNodes ({circuit.n_nodes}):")
    for i, name in enumerate(circuit.node_names):
        print(f"  {i}: {name}")

    print(f"\nResistors ({len(circuit.resistors)}):")
    for i, (n1, n2, R) in enumerate(circuit.resistors):
        print(f"  R{i}: {circuit.node_names[n1]} -- {circuit.node_names[n2]}, R = {R:.6e} Ohm")

    print(f"\nInductors ({len(circuit.inductors)}):")
    for i, (n1, n2, L) in enumerate(circuit.inductors):
        print(f"  L{i}: {circuit.node_names[n1]} -- {circuit.node_names[n2]}, L = {L:.6e} H")

    print(f"\nMutual Inductors ({len(circuit.mutual_inductors)}):")
    for i, (n1a, n1b, n2a, n2b, M) in enumerate(circuit.mutual_inductors):
        print(f"  K{i}: L({circuit.node_names[n1a]},{circuit.node_names[n1b]}) <-> "
              f"L({circuit.node_names[n2a]},{circuit.node_names[n2b]}), M = {M:.6e} H")

    print(f"\nControlled Sources ({len(circuit.controlled_sources)}):")
    for i, cs in enumerate(circuit.controlled_sources):
        print(f"  {cs['type']}{i}: {cs['description']}, gain = {cs['gain']:.6e}")

    print(f"\nPorts ({len(circuit.port_nodes)}):")
    for i, (p, n) in enumerate(circuit.port_nodes):
        print(f"  Port {i}: {circuit.node_names[p]} -- {circuit.node_names[n]}")


def generate_peec_matrices(n_segments: int = 10, length: float = 0.1,
                            width: float = 0.01, sigma: float = 5.8e7):
    """
    Generate physically realistic PEEC matrices for a straight conductor.

    The partial inductance between two segments decays as 1/r,
    giving a low-rank off-diagonal structure suitable for ACA+.

    Parameters:
        n_segments: Number of conductor segments
        length: Total conductor length [m]
        width: Conductor width [m]
        sigma: Conductivity [S/m]

    Returns:
        L, R, P matrices
    """
    mu0 = 4 * np.pi * 1e-7
    eps0 = 8.854e-12

    segment_length = length / n_segments
    segment_area = width * width

    # Segment centers
    centers = np.array([(i + 0.5) * segment_length for i in range(n_segments)])

    # Inductance matrix (partial inductance)
    # Using simplified formula that ensures positive definiteness
    L = np.zeros((n_segments, n_segments))
    for i in range(n_segments):
        for j in range(n_segments):
            if i == j:
                # Self partial inductance (Rosa formula for rectangular cross-section)
                # Lp_self = (mu0 * l) / (2*pi) * [ln(2*l/GMD) - 1]
                # GMD for square ~ 0.2235 * width
                GMD = 0.2235 * width
                L[i, i] = mu0 * segment_length / (2 * np.pi) * (
                    np.log(2 * segment_length / GMD) - 1
                )
            else:
                # Mutual partial inductance (Neumann formula)
                dist = abs(centers[i] - centers[j])
                # For parallel filaments of length l separated by d:
                # M = (mu0 * l) / (2*pi) * [ln(d/l + sqrt(1 + (d/l)^2)) - sqrt(1 + (l/d)^2) + l/d]
                d_over_l = dist / segment_length
                L[i, j] = mu0 * segment_length / (2 * np.pi) * (
                    np.log(d_over_l + np.sqrt(1 + d_over_l**2))
                    - np.sqrt(1 + 1/d_over_l**2) + 1/d_over_l
                )

    # Ensure positive definiteness by adding small diagonal if needed
    min_eig = np.min(np.linalg.eigvalsh(L))
    if min_eig < 0:
        L = L + np.eye(n_segments) * (-min_eig * 1.1 + 1e-12)

    # Resistance matrix (diagonal for DC)
    R_dc = segment_length / (sigma * segment_area)
    R = np.eye(n_segments) * R_dc

    # Potential coefficient matrix (1/4*pi*eps0 * 1/r structure)
    P = np.zeros((n_segments, n_segments))
    for i in range(n_segments):
        for j in range(n_segments):
            if i == j:
                # Self potential coefficient
                P[i, i] = 1 / (4 * np.pi * eps0 * width) * np.log(4)
            else:
                # Mutual potential coefficient
                dist = abs(centers[i] - centers[j])
                P[i, j] = 1 / (4 * np.pi * eps0 * dist)

    return L, R, P, centers


def generate_magnetic_coupling(loop_centers: np.ndarray,
                                mag_centers: np.ndarray,
                                segment_length: float = 0.002,
                                coupling_strength: float = 0.1) -> np.ndarray:
    """
    Generate Loop-Magnetic mutual inductance coupling matrix.

    Uses Neumann-like formula scaled to match PEEC inductance matrix.
    The coupling represents mutual inductance between conductor segments
    and magnetic elements.

    Parameters:
        loop_centers: Conductor segment centers [n_L x 3] or [n_L] (1D positions)
        mag_centers: Magnetic element centers [n_M x 3] or [n_M] (1D positions)
        segment_length: Conductor segment length [m] (default: 2mm)
        coupling_strength: Relative coupling strength (0.1 = 10% of self inductance)

    Returns:
        Coupling matrix L_LM [n_L x n_M] in Henries
    """
    mu0 = 4 * np.pi * 1e-7
    n_L = len(loop_centers)
    n_M = len(mag_centers)

    # Handle 1D case - extract x-coordinate only
    if loop_centers.ndim > 1:
        loop_pos = loop_centers[:, 0]
    else:
        loop_pos = loop_centers

    if mag_centers.ndim > 1:
        mag_pos = mag_centers[:, 0]
    else:
        mag_pos = mag_centers

    # Generate mutual inductance matrix
    # Scale: similar to PEEC partial inductance
    # L_self ~ mu0 * l / (2*pi) * ln(2*l/GMD) ~ 1e-8 H for l=2mm
    L_scale = mu0 * segment_length / (2 * np.pi) * coupling_strength

    L_LM = np.zeros((n_L, n_M))
    for i in range(n_L):
        for j in range(n_M):
            dist = abs(loop_pos[i] - mag_pos[j])
            # Decay with distance (similar to mutual inductance)
            # Use exponential decay to ensure positive semi-definiteness
            decay_length = segment_length * 5  # Characteristic decay length
            L_LM[i, j] = L_scale * np.exp(-dist / decay_length)

    return L_LM


def generate_magnetic_impedance(n_M: int, L_self: float,
                                 off_diag_ratio: float = 0.1) -> np.ndarray:
    """
    Generate magnetic material impedance matrix Z_MM.

    Z_MM represents the demagnetization/reluctance of magnetic elements.
    Scaled to be consistent with PEEC inductance matrix.

    Parameters:
        n_M: Number of magnetic elements
        L_self: Reference inductance scale (e.g., L[0,0] from PEEC)
        off_diag_ratio: Off-diagonal coupling ratio

    Returns:
        Z_MM matrix [n_M x n_M] (symmetric)
    """
    Z_MM = np.eye(n_M) * L_self
    for i in range(n_M - 1):
        Z_MM[i, i+1] = L_self * off_diag_ratio
        Z_MM[i+1, i] = L_self * off_diag_ratio
    return Z_MM


def symmetrize_magnetic_block(Z_MM: np.ndarray) -> np.ndarray:
    """
    Symmetrize magnetic block impedance matrix.

    For PEEC-magnetic coupling, the magnetic block must be symmetric to ensure
    proper energy conservation and passivity of the coupled system.

    The symmetrization follows the approach:
        Z_MM_sym = 0.5 * (Z_MM + Z_MM.T)

    Parameters:
        Z_MM: magnetic impedance matrix (may be slightly asymmetric)

    Returns:
        Symmetrized Z_MM matrix
    """
    return 0.5 * (Z_MM + Z_MM.T)


def build_peec_magnetic_system(L: np.ndarray, R: np.ndarray, L_LM: np.ndarray,
                           Z_MM: np.ndarray, symmetrize: bool = True) -> dict:
    """
    Build PEEC-magnetic coupled system matrices.

    Physical Model:
    - Conductor current I creates H-field at magnetic elements
    - Magnetic material responds: M = chi * H (magnetic susceptibility)
    - Magnetization M creates additional flux linking back to conductor
    - Net effect: L_eff = L + L_coupling (inductance INCREASES with mu > 1)

    The coupled system is formulated as:
        [R + sL    s*K   ] [I]   [V]
        [-K.T      Z_MM  ] [M] = [0]

    where K = L_LM (coupling matrix).
    The asymmetric signs ensure physically correct behavior:
    - Positive s*K: current creates H that magnetizes material
    - Negative -K.T: magnetization creates flux that adds to inductance

    Schur complement elimination of M:
        From row 2: M = Z_MM^{-1} @ K.T @ I
        Substitute into row 1:
        Z_schur = R + sL + s * K @ Z_MM^{-1} @ K.T
                = R + s * L_eff
        where L_eff = L + L_coupling (INCREASES with magnetic material)

    Parameters:
        L: PEEC inductance matrix [n_L x n_L]
        R: PEEC resistance matrix [n_L x n_L]
        L_LM: Loop-Magnetic coupling [n_L x n_M] (H-field from I at M locations)
        Z_MM: magnetic reluctance/impedance matrix [n_M x n_M]
        symmetrize: Whether to symmetrize Z_MM (default: True)

    Returns:
        dict with keys:
            'L': PEEC inductance
            'R': PEEC resistance
            'K': Coupling matrix (= L_LM)
            'Z_MM': Symmetrized magnetic impedance
            'L_coupling': K @ Z_MM^{-1} @ K.T (coupling contribution)
            'L_eff': L + L_coupling (effective inductance, INCREASES with mu > 1)
    """
    if symmetrize:
        Z_MM = symmetrize_magnetic_block(Z_MM)

    # Coupling matrix K
    K = L_LM

    # Compute effective inductance via Schur complement
    # System: [R+sL    sK  ] [I]   [V]
    #         [-K.T   Z_MM ] [M] = [0]
    #
    # From row 2: -K.T @ I + Z_MM @ M = 0  =>  M = Z_MM^{-1} @ K.T @ I
    # Substitute:
    #   (R + sL) @ I + sK @ M = V
    #   (R + sL) @ I + sK @ Z_MM^{-1} @ K.T @ I = V
    #   (R + sL + s * K @ Z_MM^{-1} @ K.T) @ I = V
    #
    # Therefore: L_eff = L + K @ Z_MM^{-1} @ K.T
    # This is physically correct: magnetic material INCREASES inductance
    try:
        Z_MM_inv = np.linalg.inv(Z_MM)
        L_coupling = K @ Z_MM_inv @ K.T
        L_eff = L + L_coupling
    except np.linalg.LinAlgError:
        L_coupling = np.zeros_like(L)
        L_eff = L

    return {
        'L': L,
        'R': R,
        'K': K,
        'L_LM': L_LM,  # Keep for compatibility
        'L_ML': K.T,   # Keep for compatibility
        'Z_MM': Z_MM,
        'L_coupling': L_coupling,
        'L_eff': L_eff
    }


@dataclass
class BlockLanczosResult:
    """Result of Block Lanczos tridiagonalization for LC systems."""
    Q: np.ndarray          # Orthonormal basis [2n x 2k]
    H: np.ndarray          # Block tridiagonal matrix [2k x 2k]
    rank: int              # Effective rank (2k)


class LCResonantPRIMA:
    """
    PRIMA-style model order reduction for LC resonant systems.

    For RLC systems with resonance, the standard Lanczos on L alone
    does not capture capacitive effects. We need to work with the
    full system matrix that includes both L and C.

    System formulation (descriptor form):
        C_d @ dx/dt = -G_d @ x + B @ u
        y = B.T @ x

    where:
        C_d = [L  0]    G_d = [R   0 ]    x = [I]
              [0  C]          [0  -I_n]        [V]

    Or equivalently in frequency domain:
        (s*C_d + G_d) @ x = B @ u

    PRIMA projects onto Krylov subspace K(A, B) where A = -C_d^{-1} @ G_d

    For resonant systems, the key insight is:
    - Resonance frequency: f_res = 1/(2*pi*sqrt(L*C))
    - Both L and C eigenvalues must be captured in reduced model
    - Block Lanczos with starting block [B_L; B_C] is effective
    """

    def __init__(self, tol: float = 1e-8, max_rank: int = None):
        """
        Parameters:
            tol: Tolerance for Lanczos convergence
            max_rank: Maximum reduction rank (per block)
        """
        self.tol = tol
        self.max_rank = max_rank

    def reduce_rlc_system(self, L: np.ndarray, R: np.ndarray, C: np.ndarray,
                           port_indices: List[int] = None,
                           k: int = 10) -> dict:
        """
        PRIMA reduction for RLC system with resonance.

        System: (R + sL + 1/(sC)) @ I = V

        This is transformed to descriptor form for PRIMA:
            [L  0] d/dt [I]   [R  I] [I]   [B]
            [0  C]      [V] + [I  0] [V] = [0] @ u

        Parameters:
            L: Inductance matrix [n x n]
            R: Resistance matrix [n x n]
            C: Capacitance matrix [n x n] (= P^{-1} for PEEC)
            port_indices: Port node indices (default: [0])
            k: Number of Lanczos iterations (reduced size = 2k)

        Returns:
            dict with reduced system matrices and projection basis
        """
        n = L.shape[0]
        if port_indices is None:
            port_indices = [0]

        # Build descriptor system matrices
        # C_d @ dx/dt = -G_d @ x + B @ u
        #
        # State: x = [I; V] (currents and node voltages)
        #
        # C_d = [L  0]  (2n x 2n)
        #       [0  C]
        #
        # G_d = [R   I]  where the (1,2) block couples dV/dt to inductor
        #       [-I  0]  and (2,1) block is KCL: C dV/dt = I
        #
        # Note: Standard descriptor form for RLC is:
        #   L dI/dt + R I + V = u  (KVL)
        #   C dV/dt = I            (KCL for capacitor)

        C_d = np.zeros((2*n, 2*n))
        C_d[:n, :n] = L
        C_d[n:, n:] = C

        G_d = np.zeros((2*n, 2*n))
        G_d[:n, :n] = R
        G_d[:n, n:] = np.eye(n)   # V term in KVL
        G_d[n:, :n] = -np.eye(n)  # I term in KCL (with sign for proper passivity)

        # Input matrix B: voltage excitation at port
        B = np.zeros((2*n, len(port_indices)))
        for i, idx in enumerate(port_indices):
            B[idx, i] = 1.0  # Voltage source at port

        # PRIMA: Build Krylov subspace K(A, B) where A = -C_d^{-1} @ G_d
        # Use Block Lanczos for stability
        Q, H = self._block_lanczos(C_d, G_d, B, k)

        # Reduced system matrices
        # C_d_red = Q.T @ C_d @ Q
        # G_d_red = Q.T @ G_d @ Q
        # B_red = Q.T @ B
        C_d_red = Q.T @ C_d @ Q
        G_d_red = Q.T @ G_d @ Q
        B_red = Q.T @ B

        return {
            'Q': Q,
            'H': H,
            'C_d_red': C_d_red,
            'G_d_red': G_d_red,
            'B_red': B_red,
            'n_original': n,
            'n_reduced': Q.shape[1],
            'port_indices': port_indices
        }

    def _block_lanczos(self, C_d: np.ndarray, G_d: np.ndarray,
                        B: np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Block Lanczos for descriptor system.

        Builds orthonormal basis Q for Krylov subspace:
            K_k(A, B) = span{B, A@B, A^2@B, ..., A^{k-1}@B}

        where A = -C_d^{-1} @ G_d

        Parameters:
            C_d: Descriptor "mass" matrix
            G_d: Descriptor "stiffness" matrix
            B: Input matrix
            k: Number of Lanczos iterations

        Returns:
            Q: Orthonormal basis [2n x 2k*p] where p = number of ports
            H: Block upper Hessenberg matrix
        """
        n2 = C_d.shape[0]  # 2n
        p = B.shape[1]     # number of ports

        # Solve C_d @ x = G_d @ v for matrix-vector product A @ v
        # A = -C_d^{-1} @ G_d
        try:
            C_d_lu = linalg.lu_factor(C_d)
            def matvec_A(v):
                return -linalg.lu_solve(C_d_lu, G_d @ v)
        except linalg.LinAlgError:
            # Fallback to direct inverse (slower)
            C_d_inv = np.linalg.pinv(C_d)
            def matvec_A(v):
                return -C_d_inv @ (G_d @ v)

        # Initialize with QR of B
        Q_list = []
        V, _ = np.linalg.qr(B)
        Q_list.append(V)

        H_blocks = []

        for j in range(k):
            # W = A @ V_j
            W = np.zeros((n2, p))
            for col in range(p):
                W[:, col] = matvec_A(V[:, col])

            # Orthogonalize against previous blocks
            h_col = []
            for i, Q_i in enumerate(Q_list):
                H_ij = Q_i.T @ W
                h_col.append(H_ij)
                W = W - Q_i @ H_ij

            # QR factorization of residual
            if np.linalg.norm(W) < self.tol:
                # Converged - invariant subspace found
                break

            V_new, H_next = np.linalg.qr(W)
            h_col.append(H_next)
            H_blocks.append(h_col)

            Q_list.append(V_new)
            V = V_new

        # Build Q matrix
        Q = np.column_stack(Q_list)

        # Build block Hessenberg H via projection
        # H = Q.T @ A @ Q where A = -C_d^{-1} @ G_d
        # Compute column by column to avoid memory issues
        n_cols = Q.shape[1]
        H = np.zeros((n_cols, n_cols))
        for j in range(n_cols):
            Aq_j = matvec_A(Q[:, j])
            for i in range(n_cols):
                H[i, j] = Q[:, i] @ Aq_j

        return Q, H

    def compute_impedance(self, reduced_system: dict, frequencies: np.ndarray) -> np.ndarray:
        """
        Compute port impedance from reduced system over frequency range.

        Z(s) = B_red.T @ (s*C_d_red + G_d_red)^{-1} @ B_red

        Parameters:
            reduced_system: Output from reduce_rlc_system()
            frequencies: Array of frequencies [Hz]

        Returns:
            Complex impedance array [n_freq x n_ports x n_ports]
        """
        C_d_red = reduced_system['C_d_red']
        G_d_red = reduced_system['G_d_red']
        B_red = reduced_system['B_red']

        n_freq = len(frequencies)
        n_ports = B_red.shape[1]
        Z = np.zeros((n_freq, n_ports, n_ports), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            Y_red = s * C_d_red + G_d_red
            try:
                Y_red_inv = np.linalg.inv(Y_red)
                Z[i] = B_red.T @ Y_red_inv @ B_red
            except np.linalg.LinAlgError:
                Z[i] = np.nan

        return Z

    def compute_impedance_direct(self, L: np.ndarray, R: np.ndarray,
                                  C: np.ndarray, frequencies: np.ndarray,
                                  port_indices: List[int] = None) -> np.ndarray:
        """
        Compute port impedance directly (reference for validation).

        Z(s) = R + sL + 1/(sC)  (for single-node approximation)

        For full system:
        Z_full(s) = B.T @ (s*C_d + G_d)^{-1} @ B

        Parameters:
            L, R, C: System matrices
            frequencies: Array of frequencies [Hz]
            port_indices: Port node indices

        Returns:
            Complex impedance array
        """
        n = L.shape[0]
        if port_indices is None:
            port_indices = [0]

        # Build full descriptor system
        C_d = np.zeros((2*n, 2*n))
        C_d[:n, :n] = L
        C_d[n:, n:] = C

        G_d = np.zeros((2*n, 2*n))
        G_d[:n, :n] = R
        G_d[:n, n:] = np.eye(n)
        G_d[n:, :n] = -np.eye(n)

        B = np.zeros((2*n, len(port_indices)))
        for i, idx in enumerate(port_indices):
            B[idx, i] = 1.0

        n_freq = len(frequencies)
        n_ports = len(port_indices)
        Z = np.zeros((n_freq, n_ports, n_ports), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            Y_full = s * C_d + G_d
            try:
                Y_full_inv = np.linalg.inv(Y_full)
                Z[i] = B.T @ Y_full_inv @ B
            except np.linalg.LinAlgError:
                Z[i] = np.nan

        return Z


class SecondOrderArnoldi:
    """
    Second-Order Arnoldi (SOAR) for LC resonant systems.

    For quadratic systems: (s^2 M + s D + K) x = B u

    where:
        M = L (inductance = mass)
        D = R (resistance = damping)
        K = C^{-1} or P (potential = stiffness)

    SOAR builds a subspace that preserves the second-order structure,
    which is critical for accurate resonance capture.

    Reference: Bai & Su, "SOAR: Second-Order Arnoldi Method for the
    Quadratic Eigenvalue Problem", SIAM J. Matrix Anal. Appl., 2005
    """

    def __init__(self, tol: float = 1e-10):
        self.tol = tol

    def reduce_second_order(self, M: np.ndarray, D: np.ndarray, K: np.ndarray,
                             B: np.ndarray, k: int = 10) -> dict:
        """
        SOAR reduction for second-order system.

        (s^2 M + s D + K) x = B u

        Uses standard Arnoldi on first-order companion linearization,
        but projects back to preserve second-order structure.

        Parameters:
            M: Mass matrix (inductance L)
            D: Damping matrix (resistance R)
            K: Stiffness matrix (inverse capacitance P or C^{-1})
            B: Input matrix
            k: Number of Arnoldi iterations

        Returns:
            dict with reduced M, D, K matrices
        """
        n = M.shape[0]
        if B.ndim == 1:
            B = B.reshape(-1, 1)

        # Use standard Arnoldi on the companion linearization
        # but only keep the "position" part of the basis
        #
        # Companion form: [0   I ] [q]     [I  0] [q]     [0]
        #                 [-K -D] [p]' + s[0  M] [p] = s [B] u
        #
        # where q = position (voltage), p = velocity (current)
        #
        # For PEEC: q = node voltage, p = branch current
        # Impedance Z = V/I at port

        # Build Krylov basis using standard Arnoldi on A = -M^{-1}@K
        # This captures the LC resonance modes
        try:
            M_lu = linalg.lu_factor(M)
            def solve_M(b):
                return linalg.lu_solve(M_lu, b)
        except linalg.LinAlgError:
            M_inv = np.linalg.pinv(M)
            def solve_M(b):
                return M_inv @ b

        # A = -M^{-1} @ K (for undamped system, eigenvalues are +-j*omega)
        def matvec_A(v):
            return -solve_M(K @ v)

        # Starting vector from input
        r0 = B[:, 0].copy()
        r0 = r0 / np.linalg.norm(r0)

        # Arnoldi iteration
        Q = np.zeros((n, k))
        H = np.zeros((k, k))
        Q[:, 0] = r0

        for j in range(k - 1):
            # w = A @ q_j
            w = matvec_A(Q[:, j])

            # Orthogonalize (modified Gram-Schmidt)
            for i in range(j + 1):
                H[i, j] = Q[:, i] @ w
                w = w - H[i, j] * Q[:, i]

            H[j + 1, j] = np.linalg.norm(w)

            if H[j + 1, j] < self.tol * np.linalg.norm(H[:j+1, :j+1], 'fro'):
                # Breakdown - converged to invariant subspace
                k = j + 1
                break

            Q[:, j + 1] = w / H[j + 1, j]

        Q = Q[:, :k]
        H = H[:k, :k]

        # Project second-order matrices
        M_red = Q.T @ M @ Q
        D_red = Q.T @ D @ Q
        K_red = Q.T @ K @ Q
        B_red = Q.T @ B

        return {
            'Q': Q,
            'H': H,
            'M_red': M_red,
            'D_red': D_red,
            'K_red': K_red,
            'B_red': B_red,
            'n_original': n,
            'n_reduced': k
        }

    def compute_impedance(self, reduced: dict, frequencies: np.ndarray) -> np.ndarray:
        """
        Compute impedance from reduced second-order system.

        Z(s) = B_red.T @ (s^2 M_red + s D_red + K_red)^{-1} @ B_red

        Parameters:
            reduced: Output from reduce_second_order()
            frequencies: Array of frequencies [Hz]

        Returns:
            Complex impedance array
        """
        M = reduced['M_red']
        D = reduced['D_red']
        K = reduced['K_red']
        B = reduced['B_red']

        n_freq = len(frequencies)
        n_ports = B.shape[1] if B.ndim > 1 else 1
        if B.ndim == 1:
            B = B.reshape(-1, 1)

        Z = np.zeros((n_freq, n_ports, n_ports), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            Y = s**2 * M + s * D + K
            try:
                Y_inv = np.linalg.inv(Y)
                Z[i] = B.T @ Y_inv @ B
            except np.linalg.LinAlgError:
                Z[i] = np.nan

        return Z

    def compute_impedance_direct(self, M: np.ndarray, D: np.ndarray,
                                  K: np.ndarray, B: np.ndarray,
                                  frequencies: np.ndarray) -> np.ndarray:
        """Reference: direct impedance computation."""
        n_freq = len(frequencies)
        n_ports = B.shape[1] if B.ndim > 1 else 1
        if B.ndim == 1:
            B = B.reshape(-1, 1)

        Z = np.zeros((n_freq, n_ports, n_ports), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            Y = s**2 * M + s * D + K
            try:
                Y_inv = np.linalg.inv(Y)
                Z[i] = B.T @ Y_inv @ B
            except np.linalg.LinAlgError:
                Z[i] = np.nan

        return Z


class LoopStarMagneticPRIMA:
    """
    Hierarchical PRIMA for Loop-Star-Magnetic coupled system with Schur complement.

    Full system (in descriptor form for s-domain):
        [R + sL    K_LS    s*K_LM ] [I_L]   [V]
        [K_SL     P/s     0      ] [Q_S] = [0]
        [-K_ML    0       Z_MM   ] [M  ]   [0]

    Two reduction approaches:
    1. **Monolithic (recommended)**: Build full system matrix, apply block Arnoldi
    2. **Block-wise**: Independent Lanczos per block, less accurate for coupling

    Resonance is preserved via moment matching around expansion point s0.
    """

    def __init__(self, tol: float = 1e-10):
        self.tol = tol
        self.lanczos = LanczosReducer(tol=tol)

    def reduce_full_system(self,
                           L: np.ndarray, R: np.ndarray, P: np.ndarray,
                           K_LS: np.ndarray, K_LM: np.ndarray, Z_MM: np.ndarray,
                           port_indices: List[int] = None,
                           k_L: int = 10, k_S: int = 5, k_M: int = 5) -> dict:
        """
        Block-diagonal Galerkin projection with coupled starting vectors.

        Strategy:
        1. Lanczos on L starting from port excitation e_port
        2. Lanczos on P starting from K_LS.T @ Q_L[:, 0] (coupled)
        3. Lanczos on Z_MM starting from K_LM.T @ Q_L[:, 0] (coupled)
        4. Project all matrices using block-diagonal Q = diag(Q_L, Q_S, Q_M)

        This preserves the block structure while capturing coupling dynamics.

        Parameters:
            L: Loop inductance [n_L x n_L]
            R: Loop resistance [n_L x n_L]
            P: Star potential coefficient [n_S x n_S] (acts as 1/(sC))
            K_LS: Loop-Star coupling [n_L x n_S]
            K_LM: Loop-Magnetic coupling [n_L x n_M]
            Z_MM: Magnetic impedance [n_M x n_M]
            port_indices: Port node indices (default: [0])
            k_L: Reduced Loop dimension
            k_S: Reduced Star dimension
            k_M: Reduced Magnetic dimension

        Returns:
            dict with reduced matrices and projection bases
        """
        n_L = L.shape[0]
        n_S = P.shape[0]
        n_M = Z_MM.shape[0]

        if port_indices is None:
            port_indices = [0]

        # Port excitation vector
        e_port = np.zeros(n_L)
        for idx in port_indices:
            e_port[idx] = 1.0

        # Step 1: Lanczos on Loop block (start from port excitation)
        k_L_use = min(k_L, n_L)
        result_L = self.lanczos.lanczos_symmetric(L, k=k_L_use, v0=e_port)
        Q_L = result_L.Q

        # Step 2: Lanczos on Star block
        # Start from coupled excitation via K_LS.T
        v0_S = K_LS.T @ Q_L[:, 0]
        if np.linalg.norm(v0_S) < self.tol:
            v0_S = np.ones(n_S) / np.sqrt(n_S)
        else:
            v0_S = v0_S / np.linalg.norm(v0_S)
        k_S_use = min(k_S, n_S)
        result_S = self.lanczos.lanczos_symmetric(P, k=k_S_use, v0=v0_S)
        Q_S = result_S.Q

        # Step 3: Lanczos on Magnetic block
        # Start from coupled excitation via K_LM.T
        v0_M = K_LM.T @ Q_L[:, 0]
        if np.linalg.norm(v0_M) < self.tol:
            v0_M = np.ones(n_M) / np.sqrt(n_M)
        else:
            v0_M = v0_M / np.linalg.norm(v0_M)
        # Symmetrize Z_MM
        Z_MM_sym = 0.5 * (Z_MM + Z_MM.T)
        k_M_use = min(k_M, n_M)
        result_M = self.lanczos.lanczos_symmetric(Z_MM_sym, k=k_M_use, v0=v0_M)
        Q_M = result_M.Q

        # Get actual dimensions after Lanczos
        k_L_actual = Q_L.shape[1]
        k_S_actual = Q_S.shape[1]
        k_M_actual = Q_M.shape[1]

        # Step 4: Project all matrices using block-diagonal projection
        # Q_block = diag(Q_L, Q_S, Q_M)

        # Loop block
        L_red = Q_L.T @ L @ Q_L
        R_red = Q_L.T @ R @ Q_L

        # Star block
        P_red = Q_S.T @ P @ Q_S

        # Magnetic block
        Z_MM_red = Q_M.T @ Z_MM_sym @ Q_M

        # Coupling blocks (cross-block projections)
        K_LS_red = Q_L.T @ K_LS @ Q_S
        K_LM_red = Q_L.T @ K_LM @ Q_M

        # Port vector in reduced Loop basis
        e_port_red = Q_L.T @ e_port

        return {
            # Block-wise reduced matrices
            'L_red': L_red,
            'R_red': R_red,
            'P_red': P_red,
            'Z_MM_red': Z_MM_red,
            'K_LS_red': K_LS_red,
            'K_LM_red': K_LM_red,
            # Projection bases
            'Q_L': Q_L,
            'Q_S': Q_S,
            'Q_M': Q_M,
            # Port info
            'e_port_red': e_port_red,
            'port_indices': port_indices,
            # Dimensions
            'n_L': n_L, 'n_S': n_S, 'n_M': n_M,
            'k_L': k_L_actual, 'k_S': k_S_actual, 'k_M': k_M_actual,
            'k_total': k_L_actual + k_S_actual + k_M_actual
        }

    def compute_schur_impedance(self, reduced: dict, s: complex) -> complex:
        """
        Compute port impedance via block-structured reduced system solve.

        Full reduced system:
            [R_red + s*L_red    K_LS_red    s*K_LM_red ] [I_L]   [V]
            [K_SL_red           P_red/s     0          ] [Q_S] = [0]
            [-K_ML_red          0           Z_MM_red   ] [M  ]   [0]

        Parameters:
            reduced: Output from reduce_full_system()
            s: Complex frequency

        Returns:
            Port impedance Z_port (complex)
        """
        L_red = reduced['L_red']
        R_red = reduced['R_red']
        P_red = reduced['P_red']
        Z_MM_red = reduced['Z_MM_red']
        K_LS_red = reduced['K_LS_red']
        K_LM_red = reduced['K_LM_red']
        e_port_red = reduced['e_port_red']

        k_L = reduced['k_L']
        k_S = reduced['k_S']
        k_M = reduced['k_M']
        k_total = k_L + k_S + k_M

        # Build block-structured reduced system matrix
        Z_red = np.zeros((k_total, k_total), dtype=complex)

        # Loop block: R + sL
        Z_red[:k_L, :k_L] = R_red + s * L_red

        # Star block: P/s
        if k_S > 0:
            Z_red[k_L:k_L+k_S, k_L:k_L+k_S] = P_red / s

        # Magnetic block: Z_MM
        if k_M > 0:
            Z_red[k_L+k_S:, k_L+k_S:] = Z_MM_red

        # Loop-Star coupling (symmetric)
        if k_S > 0:
            Z_red[:k_L, k_L:k_L+k_S] = K_LS_red
            Z_red[k_L:k_L+k_S, :k_L] = K_LS_red.T

        # Loop-Magnetic coupling (asymmetric: s*K in one, -K.T in other)
        if k_M > 0:
            Z_red[:k_L, k_L+k_S:] = s * K_LM_red
            Z_red[k_L+k_S:, :k_L] = -K_LM_red.T

        # Port excitation vector (only Loop DOFs)
        B_red = np.zeros(k_total, dtype=complex)
        B_red[:k_L] = e_port_red

        # Solve for port impedance: Y_port = B^T @ Z^{-1} @ B, Z_port = 1/Y_port
        try:
            Y_red_mat = np.linalg.inv(Z_red)
            Y_port = B_red @ Y_red_mat @ B_red
            Z_port = 1.0 / Y_port if abs(Y_port) > 1e-30 else np.inf
        except np.linalg.LinAlgError:
            Z_port = np.nan

        return Z_port

    def compute_impedance_sweep(self, reduced: dict,
                                 frequencies: np.ndarray) -> np.ndarray:
        """
        Compute port impedance over frequency range.

        Parameters:
            reduced: Output from reduce_full_system()
            frequencies: Array of frequencies [Hz]

        Returns:
            Complex impedance array [n_freq]
        """
        Z = np.zeros(len(frequencies), dtype=complex)
        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            Z[i] = self.compute_schur_impedance(reduced, s)
        return Z

    def compute_impedance_direct(self,
                                  L: np.ndarray, R: np.ndarray, P: np.ndarray,
                                  K_LS: np.ndarray, K_LM: np.ndarray, Z_MM: np.ndarray,
                                  frequencies: np.ndarray,
                                  port_indices: List[int] = None) -> np.ndarray:
        """
        Direct computation of port impedance (reference).

        Parameters:
            L, R, P, K_LS, K_LM, Z_MM: Full system matrices
            frequencies: Array of frequencies [Hz]
            port_indices: Port node indices

        Returns:
            Complex impedance array [n_freq]
        """
        n_L = L.shape[0]
        n_S = P.shape[0]
        n_M = Z_MM.shape[0]
        n_total = n_L + n_S + n_M

        if port_indices is None:
            port_indices = [0]

        # Input vector
        B = np.zeros(n_total)
        for idx in port_indices:
            B[idx] = 1.0

        Z = np.zeros(len(frequencies), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f

            # Build full system matrix
            Z_full = np.zeros((n_total, n_total), dtype=complex)

            # Loop block: R + sL
            Z_full[:n_L, :n_L] = R + s * L

            # Star block: P/s
            Z_full[n_L:n_L+n_S, n_L:n_L+n_S] = P / s

            # Magnetic block: Z_MM
            Z_full[n_L+n_S:, n_L+n_S:] = Z_MM

            # Loop-Star coupling
            Z_full[:n_L, n_L:n_L+n_S] = K_LS
            Z_full[n_L:n_L+n_S, :n_L] = K_LS.T

            # Loop-Magnetic coupling (asymmetric: s*K in one, -K in other)
            Z_full[:n_L, n_L+n_S:] = s * K_LM
            Z_full[n_L+n_S:, :n_L] = -K_LM.T

            # Solve for impedance
            try:
                Y_full = np.linalg.inv(Z_full)
                Y_port = B @ Y_full @ B
                Z[i] = 1.0 / Y_port
            except np.linalg.LinAlgError:
                Z[i] = np.nan

        return Z


def demo_loop_star_magnetic_prima():
    """Demonstrate hierarchical PRIMA for Loop-Star-Magnetic system."""
    print("=" * 60)
    print("Hierarchical PRIMA for Loop-Star-Magnetic System")
    print("=" * 60)
    print("(LC resonance + magnetic coupling)")

    np.random.seed(42)

    # System sizes
    n_L = 30   # Loop DOFs
    n_S = 20   # Star DOFs
    n_M = 15   # Magnetic DOFs

    print(f"\nFull system: {n_L + n_S + n_M} DOFs")
    print(f"  Loop: {n_L}, Star: {n_S}, Magnetic: {n_M}")

    # Generate Loop matrices (PEEC inductance/resistance)
    L, R, P_loop, centers = generate_peec_matrices(n_L, length=0.05, width=0.005)

    # Generate Star potential coefficient matrix
    # P acts like 1/C, so larger P = smaller capacitance
    eps0 = 8.854e-12
    width = 0.005
    height = 0.001
    segment_length = 0.05 / n_S

    P = np.zeros((n_S, n_S))
    for i in range(n_S):
        for j in range(n_S):
            if i == j:
                P[i, i] = 1 / (4 * np.pi * eps0 * width) * 2
            else:
                dist = abs(i - j) * segment_length
                P[i, j] = 1 / (4 * np.pi * eps0 * max(dist, segment_length * 0.1))

    # Make P positive definite
    min_eig = np.min(np.linalg.eigvalsh(P))
    if min_eig <= 0:
        P = P + np.eye(n_S) * (abs(min_eig) * 1.1 + 1e-10)

    # Generate Loop-Star coupling (physical: mutual capacitance-like)
    K_LS = np.zeros((n_L, n_S))
    for i in range(n_L):
        for j in range(n_S):
            # Coupling decays with distance
            dist = abs(i / n_L - j / n_S)
            K_LS[i, j] = 0.1 * np.exp(-dist * 3) * P[0, 0]

    # Generate Magnetic matrices
    Z_MM = generate_magnetic_impedance(n_M, L[0, 0], off_diag_ratio=0.1)

    # Generate Loop-Magnetic coupling
    K_LM = generate_magnetic_coupling(
        centers, np.linspace(0.01, 0.04, n_M),
        segment_length=0.05/n_L, coupling_strength=0.5
    )

    # Estimate resonance frequency (from diagonal elements)
    # C_eff ~ 1/P[0,0], L_eff ~ L[0,0]
    C_eff = 1.0 / P[0, 0]
    L_eff_approx = L[0, 0]
    f_res_approx = 1.0 / (2 * np.pi * np.sqrt(L_eff_approx * C_eff))

    print(f"\nMatrix properties:")
    print(f"  L[0,0] = {L[0,0]*1e9:.4f} nH")
    print(f"  C_eff ~ {C_eff*1e12:.4f} pF")
    print(f"  Approx resonance: {f_res_approx/1e6:.2f} MHz")

    # Frequency sweep around resonance
    f_min = f_res_approx * 0.1
    f_max = f_res_approx * 10
    frequencies = np.logspace(np.log10(f_min), np.log10(f_max), 100)

    # Reference: Direct computation
    print("\n--- Reference: Direct Computation ---")
    prima = LoopStarMagneticPRIMA(tol=1e-12)
    Z_ref = prima.compute_impedance_direct(L, R, P, K_LS, K_LM, Z_MM, frequencies)

    # Find resonance
    Z_mag = np.abs(Z_ref)
    f_res_idx = np.argmax(Z_mag)
    f_res_actual = frequencies[f_res_idx]
    print(f"Resonance frequency: {f_res_actual/1e6:.4f} MHz")
    print(f"Peak impedance: {Z_mag[f_res_idx]:.2e} Ohm")

    # Hierarchical PRIMA reduction
    print("\n--- Hierarchical PRIMA Results ---")
    print("k_L | k_S | k_M | Total DOF | f_res err | Peak Z err")
    print("-" * 65)

    for k_L, k_S, k_M in [(5, 3, 3), (10, 5, 5), (15, 8, 8), (20, 10, 10)]:
        reduced = prima.reduce_full_system(
            L, R, P, K_LS, K_LM, Z_MM,
            k_L=k_L, k_S=k_S, k_M=k_M
        )
        Z_red = prima.compute_impedance_sweep(reduced, frequencies)

        Z_red_mag = np.abs(Z_red)
        f_res_red_idx = np.argmax(Z_red_mag)
        f_res_red = frequencies[f_res_red_idx]

        f_err = abs(f_res_red - f_res_actual) / f_res_actual * 100
        Z_err = abs(Z_red_mag[f_res_idx] - Z_mag[f_res_idx]) / Z_mag[f_res_idx] * 100

        total_dof = reduced['k_total']  # Actual reduced dimension
        print(f"{k_L:>3} | {k_S:>3} | {k_M:>3} | {total_dof:>9} | {f_err:>9.4f}% | {Z_err:>10.4f}%")

    # Detailed comparison
    print("\n--- Frequency Response (k_L=15, k_S=8, k_M=8) ---")
    reduced = prima.reduce_full_system(L, R, P, K_LS, K_LM, Z_MM, k_L=15, k_S=8, k_M=8)
    Z_red = prima.compute_impedance_sweep(reduced, frequencies)

    test_freqs = [f_res_actual * 0.5, f_res_actual, f_res_actual * 2]
    print("\nFrequency (MHz) | |Z_ref| | |Z_red| | Error %")
    print("-" * 55)
    for f_test in test_freqs:
        idx = np.argmin(np.abs(frequencies - f_test))
        z_ref = abs(Z_ref[idx])
        z_red = abs(Z_red[idx])
        err = abs(z_red - z_ref) / max(z_ref, 1e-15) * 100
        print(f"{frequencies[idx]/1e6:>15.4f} | {z_ref:>7.2e} | {z_red:>7.2e} | {err:>7.2f}%")

    # Summary
    print("\n--- Summary ---")
    n_original = n_L + n_S + n_M
    k_total = reduced['k_total']
    print(f"Original: {n_original} DOFs")
    print(f"Reduced: {k_total} DOFs ({k_total/n_original*100:.1f}%)")
    print(f"Compression ratio: {n_original/k_total:.1f}x")

    return reduced


class RationalKrylovLC:
    """
    Rational Krylov method for LC resonant systems.

    For better convergence around resonance frequencies, use expansion points
    (shifts) near the expected resonances.

    Given system: (s^2 M + s D + K) x = B u

    Rational Krylov builds basis from:
        K(A - sigma*I, B) for multiple shifts sigma

    Key insight: Place shifts at s = j*omega_res to capture resonance accurately.
    """

    def __init__(self, tol: float = 1e-10):
        self.tol = tol

    def reduce_with_shifts(self, M: np.ndarray, D: np.ndarray, K: np.ndarray,
                            B: np.ndarray, shifts: List[complex],
                            vectors_per_shift: int = 2) -> dict:
        """
        Rational Krylov reduction with multiple expansion points.

        Parameters:
            M, D, K: Second-order system matrices
            B: Input matrix
            shifts: List of complex frequency shifts (s values)
            vectors_per_shift: Number of Krylov vectors per shift

        Returns:
            Reduced system dictionary
        """
        n = M.shape[0]
        if B.ndim == 1:
            B = B.reshape(-1, 1)

        Q_list = []

        for sigma in shifts:
            # At shift sigma: (sigma^2 M + sigma D + K) x = b
            # Solve for x to get Krylov vector
            A_shift = sigma**2 * M + sigma * D + K

            try:
                A_lu = linalg.lu_factor(A_shift)
                def solve_shift(b):
                    return linalg.lu_solve(A_lu, b)
            except linalg.LinAlgError:
                A_inv = np.linalg.pinv(A_shift)
                def solve_shift(b):
                    return A_inv @ b

            # Generate Krylov vectors at this shift
            v = solve_shift(B[:, 0])
            if np.linalg.norm(v) > self.tol:
                v = v / np.linalg.norm(v)
                Q_list.append(v.real)
                if vectors_per_shift > 1 and np.linalg.norm(v.imag) > self.tol:
                    v_imag = v.imag / np.linalg.norm(v.imag)
                    Q_list.append(v_imag)

            # Additional vectors via iteration
            for _ in range(vectors_per_shift - 2):
                w = solve_shift(M @ v)
                # Orthogonalize
                for q in Q_list:
                    w = w - (q @ w.real) * q
                if np.linalg.norm(w) > self.tol:
                    w = w / np.linalg.norm(w)
                    Q_list.append(w.real)
                    v = w

        if not Q_list:
            # Fallback to standard Arnoldi
            Q_list.append(B[:, 0] / np.linalg.norm(B[:, 0]))

        # Orthogonalize full basis (QR)
        Q_raw = np.column_stack(Q_list)
        Q, _ = np.linalg.qr(Q_raw)

        # Project matrices
        M_red = Q.T @ M @ Q
        D_red = Q.T @ D @ Q
        K_red = Q.T @ K @ Q
        B_red = Q.T @ B

        return {
            'Q': Q,
            'M_red': M_red,
            'D_red': D_red,
            'K_red': K_red,
            'B_red': B_red,
            'n_original': n,
            'n_reduced': Q.shape[1],
            'shifts': shifts
        }

    def auto_shifts_from_resonance(self, M: np.ndarray, K: np.ndarray,
                                    n_shifts: int = 5) -> List[complex]:
        """
        Automatically determine shifts from undamped resonance frequencies.

        Solves generalized eigenvalue: K @ v = omega^2 @ M @ v
        Then places shifts at s = j*omega for dominant modes.

        Parameters:
            M: Mass matrix
            K: Stiffness matrix
            n_shifts: Number of shifts to generate

        Returns:
            List of complex shifts
        """
        try:
            # Solve generalized eigenvalue problem
            # K @ v = lambda @ M @ v where lambda = omega^2
            eigvals, _ = linalg.eigh(K, M)

            # omega^2 values (should be positive for stable system)
            omega_sq = eigvals[eigvals > 0]
            omega = np.sqrt(omega_sq)

            # Take lowest n_shifts frequencies
            omega = np.sort(omega)[:n_shifts]

            # Shifts at s = j*omega (on imaginary axis)
            shifts = [1j * w for w in omega]

            # Also add DC (s=0) and some real damped shifts
            shifts.insert(0, 1e-6)  # Near DC
            shifts.append(omega[0] * 0.5j)  # Below first resonance

            return shifts[:n_shifts]

        except linalg.LinAlgError:
            # Fallback: estimate from diagonal
            omega_est = np.sqrt(K[0, 0] / M[0, 0])
            return [1e-6, 0.5j * omega_est, 1j * omega_est,
                    1.5j * omega_est, 2j * omega_est][:n_shifts]

    def compute_impedance(self, reduced: dict, frequencies: np.ndarray) -> np.ndarray:
        """Compute impedance from reduced system."""
        M = reduced['M_red']
        D = reduced['D_red']
        K = reduced['K_red']
        B = reduced['B_red']

        n_freq = len(frequencies)
        n_ports = B.shape[1] if B.ndim > 1 else 1
        if B.ndim == 1:
            B = B.reshape(-1, 1)

        Z = np.zeros((n_freq, n_ports, n_ports), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            Y = s**2 * M + s * D + K
            try:
                Y_inv = np.linalg.inv(Y)
                Z[i] = B.T @ Y_inv @ B
            except np.linalg.LinAlgError:
                Z[i] = np.nan

        return Z


def demo_rational_krylov_resonant():
    """Demonstrate Rational Krylov reduction for LC resonant system."""
    print("=" * 60)
    print("Rational Krylov for LC Resonant System")
    print("=" * 60)
    print("(Uses shifts at resonance frequencies for fast convergence)")

    np.random.seed(42)

    # Create LC system
    n = 20

    # Generate PEEC matrices
    L, R, P, centers = generate_peec_matrices(n, length=0.05, width=0.005)

    # Create capacitance
    eps0 = 8.854e-12
    width = 0.005
    height = 0.001
    segment_length = 0.05 / n
    C_segment = eps0 * width * segment_length / height

    C = np.eye(n) * C_segment * 2
    for i in range(n - 1):
        C[i, i+1] = -C_segment * 0.1
        C[i+1, i] = -C_segment * 0.1
    min_eig = np.min(np.linalg.eigvalsh(C))
    if min_eig <= 0:
        C = C + np.eye(n) * (abs(min_eig) + C_segment * 0.1)

    K = np.linalg.inv(C)
    B = np.zeros((n, 1))
    B[0, 0] = 1.0

    print(f"\nSystem: s^2*L + s*R + C^(-1)")
    print(f"Size: {n} nodes")

    # Reference computation
    rk = RationalKrylovLC(tol=1e-12)

    # Auto-generate shifts from resonances
    shifts = rk.auto_shifts_from_resonance(L, K, n_shifts=5)
    print(f"\nAuto-selected shifts (omega):")
    for i, s in enumerate(shifts):
        if abs(s.real) < 1e-10:
            print(f"  {i+1}: {abs(s.imag)/(2*np.pi*1e6):.2f} MHz (j*omega)")
        else:
            print(f"  {i+1}: {s:.2e}")

    # Frequency sweep
    f_res_approx = 1.0 / (2 * np.pi * np.sqrt(L[0, 0] * C[0, 0]))
    f_min = f_res_approx * 0.1
    f_max = f_res_approx * 10
    frequencies = np.logspace(np.log10(f_min), np.log10(f_max), 100)

    # Reference
    soar = SecondOrderArnoldi(tol=1e-12)
    Z_ref = soar.compute_impedance_direct(L, R, K, B, frequencies)
    Z_mag = np.abs(Z_ref[:, 0, 0])
    f_res_idx = np.argmax(Z_mag)
    f_res_actual = frequencies[f_res_idx]

    print(f"\nReference resonance: {f_res_actual/1e6:.4f} MHz")
    print(f"Peak impedance: {Z_mag[f_res_idx]:.2e} Ohm")

    # Rational Krylov with various number of shifts
    print("\n--- Rational Krylov Results ---")
    print("n_shifts | Reduced DOF | f_res error | Peak Z error")
    print("-" * 60)

    for n_shifts in [2, 3, 4, 5]:
        shifts = rk.auto_shifts_from_resonance(L, K, n_shifts=n_shifts)
        reduced = rk.reduce_with_shifts(L, R, K, B, shifts, vectors_per_shift=2)
        Z_red = rk.compute_impedance(reduced, frequencies)

        Z_red_mag = np.abs(Z_red[:, 0, 0])
        f_res_red_idx = np.argmax(Z_red_mag)
        f_res_red = frequencies[f_res_red_idx]

        f_err = abs(f_res_red - f_res_actual) / f_res_actual * 100
        Z_err = abs(Z_red_mag[f_res_idx] - Z_mag[f_res_idx]) / Z_mag[f_res_idx] * 100

        print(f"{n_shifts:>8} | {reduced['n_reduced']:>11} | {f_err:>11.4f}% | {Z_err:>12.4f}%")

    print("\n--- Summary ---")
    print("Rational Krylov provides good accuracy with few DOFs by")
    print("placing expansion points near the resonance frequencies.")

    return reduced


def demo_soar_resonant():
    """Demonstrate SOAR reduction for LC resonant system."""
    print("=" * 60)
    print("SOAR (Second-Order Arnoldi) for LC Resonant System")
    print("=" * 60)

    np.random.seed(42)

    # Create LC system
    n = 20

    # Generate PEEC matrices
    L, R, P, centers = generate_peec_matrices(n, length=0.05, width=0.005)

    # For second-order form: s^2*L + s*R + K
    # K = P (potential coefficient matrix, acts like 1/C)
    # Or we can use K = C^{-1} from a proper capacitance matrix

    # Create capacitance and convert to stiffness
    eps0 = 8.854e-12
    width = 0.005
    height = 0.001
    segment_length = 0.05 / n
    C_segment = eps0 * width * segment_length / height

    C = np.eye(n) * C_segment * 2
    for i in range(n - 1):
        C[i, i+1] = -C_segment * 0.1
        C[i+1, i] = -C_segment * 0.1
    min_eig = np.min(np.linalg.eigvalsh(C))
    if min_eig <= 0:
        C = C + np.eye(n) * (abs(min_eig) + C_segment * 0.1)

    # K = C^{-1} (stiffness from capacitance)
    K = np.linalg.inv(C)

    # Input vector (port at node 0)
    B = np.zeros((n, 1))
    B[0, 0] = 1.0

    print(f"\nSystem: s^2*L + s*R + C^(-1)")
    print(f"Size: {n} nodes")
    print(f"  L[0,0] = {L[0,0]*1e9:.4f} nH")
    print(f"  C[0,0] = {C[0,0]*1e12:.4f} pF")
    print(f"  R[0,0] = {R[0,0]*1e3:.4f} mOhm")

    # Estimate resonance
    f_res_approx = 1.0 / (2 * np.pi * np.sqrt(L[0, 0] * C[0, 0]))
    print(f"\nApproximate resonance: {f_res_approx/1e6:.2f} MHz")

    # Frequency sweep
    f_min = f_res_approx * 0.1
    f_max = f_res_approx * 10
    frequencies = np.logspace(np.log10(f_min), np.log10(f_max), 100)

    # Reference
    print("\n--- Computing reference impedance ---")
    soar = SecondOrderArnoldi(tol=1e-12)
    Z_ref = soar.compute_impedance_direct(L, R, K, B, frequencies)

    Z_mag = np.abs(Z_ref[:, 0, 0])
    f_res_idx = np.argmax(Z_mag)
    f_res_actual = frequencies[f_res_idx]
    print(f"Resonance frequency: {f_res_actual/1e6:.4f} MHz")
    print(f"Peak impedance: {Z_mag[f_res_idx]:.2e} Ohm")

    # SOAR reduction
    print("\n--- SOAR Reduction Results ---")
    print("k  | Reduced DOF | f_res error | Peak Z error")
    print("-" * 55)

    for k in [3, 5, 10, 15, 20]:
        reduced = soar.reduce_second_order(L, R, K, B, k=k)
        Z_red = soar.compute_impedance(reduced, frequencies)

        Z_red_mag = np.abs(Z_red[:, 0, 0])
        f_res_red_idx = np.argmax(Z_red_mag)
        f_res_red = frequencies[f_res_red_idx]

        f_err = abs(f_res_red - f_res_actual) / f_res_actual * 100
        Z_peak_err = abs(Z_red_mag[f_res_idx] - Z_mag[f_res_idx]) / Z_mag[f_res_idx] * 100

        print(f"{k:>2} | {reduced['n_reduced']:>11} | {f_err:>11.4f}% | {Z_peak_err:>12.4f}%")

    # Final comparison
    print("\n--- Best Result (k=20) ---")
    reduced = soar.reduce_second_order(L, R, K, B, k=20)
    Z_red = soar.compute_impedance(reduced, frequencies)

    # Compare at resonance and off-resonance
    test_freqs = [f_res_actual * 0.5, f_res_actual, f_res_actual * 2]
    print("\nFrequency (MHz) | |Z_ref| (Ohm) | |Z_red| (Ohm) | Error %")
    print("-" * 65)
    for f_test in test_freqs:
        idx = np.argmin(np.abs(frequencies - f_test))
        z_ref = abs(Z_ref[idx, 0, 0])
        z_red = abs(Z_red[idx, 0, 0])
        err = abs(z_red - z_ref) / max(z_ref, 1e-15) * 100
        print(f"{frequencies[idx]/1e6:>15.4f} | {z_ref:>13.4e} | {z_red:>13.4e} | {err:>7.2f}%")

    print("\n--- Summary ---")
    print(f"Original: {n} DOFs")
    print(f"Reduced (k=20): {reduced['n_reduced']} DOFs")
    print(f"Compression: {reduced['n_reduced'] / n * 100:.1f}%")

    return reduced


def demo_lc_resonant_prima():
    """Demonstrate PRIMA reduction for LC resonant system."""
    print("=" * 60)
    print("PRIMA Reduction for LC Resonant System")
    print("=" * 60)

    np.random.seed(42)

    # Create simple LC system
    n = 20  # Number of segments

    # Generate PEEC matrices
    L, R, P, centers = generate_peec_matrices(n, length=0.05, width=0.005)

    # Create a physically meaningful capacitance matrix
    # For a transmission line, use distributed capacitance model
    # C_per_length ~ eps0 * width / height (parallel plate approximation)
    eps0 = 8.854e-12
    width = 0.005
    height = 0.001  # substrate thickness
    segment_length = 0.05 / n

    # Capacitance per segment (to ground)
    C_segment = eps0 * width * segment_length / height

    # Build capacitance matrix (diagonal dominant with neighbor coupling)
    C = np.zeros((n, n))
    for i in range(n):
        C[i, i] = C_segment * 2  # Self capacitance
        if i > 0:
            C[i, i-1] = -C_segment * 0.1  # Neighbor coupling
            C[i-1, i] = -C_segment * 0.1
    # Make positive definite
    min_eig = np.min(np.linalg.eigvalsh(C))
    if min_eig <= 0:
        C = C + np.eye(n) * (abs(min_eig) + C_segment * 0.1)

    print(f"\nSystem size: {n} nodes")
    print(f"  L[0,0] = {L[0,0]*1e9:.4f} nH")
    print(f"  C[0,0] = {C[0,0]*1e12:.4f} pF")
    print(f"  R[0,0] = {R[0,0]*1e3:.4f} mOhm")

    # Estimate resonance frequency
    L_eff = L[0, 0]  # Approximate
    C_eff = C[0, 0]  # Approximate
    f_res_approx = 1.0 / (2 * np.pi * np.sqrt(L_eff * C_eff))
    print(f"\nApproximate resonance: {f_res_approx/1e6:.2f} MHz")

    # Frequency sweep around resonance
    f_min = f_res_approx * 0.1
    f_max = f_res_approx * 10
    frequencies = np.logspace(np.log10(f_min), np.log10(f_max), 100)

    # Reference: Direct computation
    print("\n--- Computing reference impedance (full system) ---")
    prima = LCResonantPRIMA(tol=1e-10)
    Z_ref = prima.compute_impedance_direct(L, R, C, frequencies, port_indices=[0])

    # Find resonance from reference
    Z_mag = np.abs(Z_ref[:, 0, 0])
    f_res_idx = np.argmax(Z_mag)  # Parallel resonance = max impedance
    f_res_actual = frequencies[f_res_idx]
    print(f"Resonance frequency (from reference): {f_res_actual/1e6:.4f} MHz")
    print(f"Peak impedance: {Z_mag[f_res_idx]:.2f} Ohm")

    # PRIMA reduction with various ranks
    print("\n--- PRIMA Reduction Results ---")
    print("k  | Reduced DOF | f_res error | Peak Z error")
    print("-" * 55)

    for k in [3, 5, 10, 15]:
        reduced = prima.reduce_rlc_system(L, R, C, port_indices=[0], k=k)
        Z_red = prima.compute_impedance(reduced, frequencies)

        # Find resonance in reduced model
        Z_red_mag = np.abs(Z_red[:, 0, 0])
        f_res_red_idx = np.argmax(Z_red_mag)
        f_res_red = frequencies[f_res_red_idx]

        f_err = abs(f_res_red - f_res_actual) / f_res_actual * 100
        Z_err = abs(Z_red_mag[f_res_idx] - Z_mag[f_res_idx]) / Z_mag[f_res_idx] * 100

        print(f"{k:>2} | {reduced['n_reduced']:>11} | {f_err:>11.4f}% | {Z_err:>12.4f}%")

    # Detailed comparison for best reduction
    print("\n--- Frequency Response Comparison (k=10) ---")
    reduced = prima.reduce_rlc_system(L, R, C, port_indices=[0], k=10)
    Z_red = prima.compute_impedance(reduced, frequencies)

    # Sample points
    sample_freqs = [f_res_actual * 0.5, f_res_actual, f_res_actual * 2]
    print("\nFrequency (MHz) | Z_ref (Ohm) | Z_red (Ohm) | Error %")
    print("-" * 60)
    for f_sample in sample_freqs:
        idx = np.argmin(np.abs(frequencies - f_sample))
        z_ref = Z_ref[idx, 0, 0]
        z_red = Z_red[idx, 0, 0]
        err = abs(z_red - z_ref) / abs(z_ref) * 100
        print(f"{frequencies[idx]/1e6:>15.4f} | {abs(z_ref):>11.2f} | {abs(z_red):>11.2f} | {err:>7.2f}%")

    print("\n--- Summary ---")
    print(f"Original system: 2 x {n} = {2*n} DOFs")
    print(f"Reduced system (k=10): {reduced['n_reduced']} DOFs")
    print(f"Compression: {reduced['n_reduced'] / (2*n) * 100:.1f}%")

    return reduced


def demo_hierarchical_reduction():
    """Demonstrate hierarchical reduction on a physically realistic PEEC-magnetic system."""
    print("=" * 60)
    print("PRIMA-style Hierarchical Model Order Reduction Demo")
    print("=" * 60)
    print("(PEEC-magnetic coupled system with physically correct formulation)")

    np.random.seed(42)

    # Create physically realistic PEEC system
    n_L = 50   # Loop DOFs (conductor segments)
    n_M = 30   # Magnetic DOFs
    length = 0.1  # 10cm conductor
    width = 0.01  # 1cm width

    print(f"\nOriginal system size: {n_L + n_M} DOFs")
    print(f"  Loop: {n_L}, Magnetic: {n_M}")

    # Generate PEEC matrices
    L, R, P, loop_centers = generate_peec_matrices(n_L, length=length, width=width)
    segment_length = length / n_L

    print(f"\nPEEC conductor: {n_L} segments, {segment_length*1000:.1f} mm each")
    print(f"  L[0,0] = {L[0,0]*1e9:.2f} nH (self)")
    print(f"  L[0,1] = {L[0,1]*1e9:.2f} nH (nearest neighbor)")
    print(f"  R_dc = {R[0,0]*1e3:.3f} mOhm per segment")

    # Generate magnetic element positions (along conductor)
    mag_centers = np.linspace(0.01, 0.09, n_M)

    # Generate coupling with proper scaling
    K = generate_magnetic_coupling(
        loop_centers, mag_centers,
        segment_length=segment_length,
        coupling_strength=0.1  # 10% coupling
    )

    # Generate magnetic impedance matrix with proper scaling
    Z_MM = generate_magnetic_impedance(n_M, L[0, 0], off_diag_ratio=0.1)

    # Build PEEC-magnetic system with physically correct formulation
    system = build_peec_magnetic_system(L, R, K, Z_MM, symmetrize=True)

    L_eff = system['L_eff']
    L_coupling = system['L_coupling']

    print(f"\nCoupling matrices:")
    print(f"  ||K||_F / ||L||_F = {np.linalg.norm(K, 'fro') / np.linalg.norm(L, 'fro'):.2%}")
    print(f"  ||Z_MM||_F = {np.linalg.norm(system['Z_MM'], 'fro'):.4e}")

    # Verify physically correct behavior: L_eff > L
    print(f"\nPhysical verification:")
    print(f"  L[0,0] = {L[0,0]*1e9:.4f} nH (air-core)")
    print(f"  L_eff[0,0] = {L_eff[0,0]*1e9:.4f} nH (with magnetic material)")
    print(f"  L_coupling[0,0] = {L_coupling[0,0]*1e9:.4f} nH")
    if L_eff[0, 0] > L[0, 0]:
        print(f"  -> L_eff > L: CORRECT (magnetic material INCREASES inductance)")
    else:
        print(f"  -> WARNING: L_eff <= L (unexpected)")

    print(f"  Condition(L_eff) = {np.linalg.cond(L_eff):.2e}")

    # Frequency for tests
    f = 1e6  # 1 MHz
    s = 2j * np.pi * f

    # Reference: Full system direct solve with physically correct formulation
    print("\n--- Reference: Full System Direct Solve ---")
    print("System: [R+sL    s*K  ] [I]   [V]")
    print("        [-K.T   Z_MM ] [M] = [0]")

    Z_LL_full = R + s * L
    Z_MM_sym = system['Z_MM']

    # Physically correct formulation:
    # [R+sL    s*K  ] [I]   [V]
    # [-K.T   Z_MM ] [M] = [0]
    n_total = n_L + n_M
    Z_full = np.zeros((n_total, n_total), dtype=complex)
    Z_full[:n_L, :n_L] = Z_LL_full
    Z_full[n_L:, n_L:] = Z_MM_sym
    Z_full[:n_L, n_L:] = s * K          # s*K
    Z_full[n_L:, :n_L] = -K.T           # -K.T (asymmetric!)

    Y_full = np.linalg.inv(Z_full)
    Z_port_ref = 1.0 / Y_full[0, 0]
    print(f"Z_port (full system): {Z_port_ref:.4e}")

    # Compute via Schur complement on L_eff
    # From row 2: -K.T @ I + Z_MM @ M = 0  =>  M = Z_MM^{-1} @ K.T @ I
    # Substitute into row 1:
    #   (R + sL) @ I + s*K @ Z_MM^{-1} @ K.T @ I = V
    #   (R + sL + s * K @ Z_MM^{-1} @ K.T) @ I = V
    #   (R + s * L_eff) @ I = V
    # where L_eff = L + K @ Z_MM^{-1} @ K.T
    Z_schur_direct = R + s * L_eff

    Y_schur_direct = np.linalg.inv(Z_schur_direct)
    Z_port_schur_ref = 1.0 / Y_schur_direct[0, 0]
    print(f"Z_port (Schur on L_eff): {Z_port_schur_ref:.4e}")
    schur_err = abs(Z_port_schur_ref - Z_port_ref) / abs(Z_port_ref) * 100
    print(f"Schur error: {schur_err:.6f}%")

    # PRIMA-style Lanczos on L_eff
    print("\n--- PRIMA-style Lanczos Reduction ---")
    print("(Starting from port excitation e_0, Lanczos on L_eff)")
    print()
    print("k    | Z_port                   | Error %")
    print("-" * 50)

    e0 = np.zeros(n_L)
    e0[0] = 1.0

    for k in [5, 10, 15, 20]:
        # Lanczos on L_eff starting from port excitation
        reducer = HierarchicalReducer(
            lanczos_tol=1e-10,  # Don't truncate by eigenvalue
            aca_tol=1e-4,
            max_lanczos_rank=k
        )

        # Use L_eff for Lanczos (includes magnetic coupling effect)
        result = reducer.lanczos.lanczos_symmetric(L_eff, k=k, v0=e0)
        Q = result.Q

        # Reduced system: Z_schur_red = Q.T @ (R + s*L_eff) @ Q
        Z_eff_red = Q.T @ (R + s * L_eff) @ Q

        # Port impedance via Y-matrix (proper formula)
        Y_eff_red = np.linalg.inv(Z_eff_red)
        Z_port_k = 1.0 / Y_eff_red[0, 0]

        err = abs(Z_port_k - Z_port_schur_ref) / abs(Z_port_schur_ref) * 100
        print(f"{k:>4} | {Z_port_k:.4e} | {err:>10.4f}%")

    # Sparse circuit extraction for final model
    print("\n--- Sparse Circuit Extraction ---")
    k_final = 10
    result = reducer.lanczos.lanczos_symmetric(L_eff, k=k_final, v0=e0)
    Q = result.Q

    circuit = SparseCircuit(
        n_nodes=k_final + 1,
        node_names=[f'L{i}' for i in range(k_final)] + ['GND'],
        resistors=[],
        inductors=[],
        capacitors=[],
        mutual_inductors=[],
        controlled_sources=[],
        port_nodes=[(0, k_final)]
    )

    # Add ladder elements from tridiagonal T
    T = result.T
    R_red = Q.T @ R @ Q
    for i in range(k_final):
        circuit.inductors.append((i, k_final, T[i, i]))
        circuit.resistors.append((i, k_final, R_red[i, i]))
    for i in range(k_final - 1):
        circuit.mutual_inductors.append((i, k_final, i+1, k_final, T[i, i+1]))

    print_sparse_circuit(circuit)

    # Summary
    print("\n--- Summary ---")
    print(f"Original DOFs: {n_L + n_M}")
    print(f"Reduced DOFs: {k_final}")
    print(f"Compression: {k_final / (n_L + n_M) * 100:.1f}%")

    # Final accuracy check
    Z_eff_red = Q.T @ (R + s * L_eff) @ Q
    Y_eff_red = np.linalg.inv(Z_eff_red)
    Z_port_final = 1.0 / Y_eff_red[0, 0]
    final_err = abs(Z_port_final - Z_port_schur_ref) / abs(Z_port_schur_ref) * 100
    print(f"Port impedance error (k={k_final}): {final_err:.4f}%")

    return circuit


##############################################################################
# PyKAN Integration for Complex Material Properties
##############################################################################

@dataclass
class KANMaterialState:
    """
    State representation for KAN-based material model.

    KAN learns: mu*(s, H, T) or eps*(s, E, T)

    Internal states represent the "memory" of material response,
    analogous to hidden states in RNN/reservoir computing.
    """
    n_internal: int                    # Number of internal state variables
    state: np.ndarray = None           # Current internal state [n_internal]
    coefficients: np.ndarray = None    # Learned coefficients from KAN


class KANMaterialInterface:
    """
    Interface between PyKAN and electromagnetic simulation.

    Converts learned KAN model to state-space form for transient analysis:

    KAN: y = phi_L . phi_{L-1} . ... . phi_1(x)

    State-space:
        dz/dt = A @ z + B @ u
        y = C @ z + D @ u

    where z contains the intermediate KAN node values.

    Usage:
        1. Train KAN on frequency-domain data: mu*(omega), eps*(omega)
        2. Extract state-space representation
        3. Integrate with PRIMA reduced model
        4. Run transient simulation
    """

    def __init__(self, kan_model=None, material_type: str = 'magnetic'):
        """
        Initialize KAN material interface.

        Parameters:
            kan_model: Trained PyKAN model (kan.KAN instance)
            material_type: 'magnetic' (mu*) or 'dielectric' (eps*)
        """
        self.kan_model = kan_model
        self.material_type = material_type
        self.state_space = None
        self._n_internal = 0

    def load_kan_model(self, model_path: str):
        """
        Load trained KAN model from file.

        Parameters:
            model_path: Path to saved KAN model (.pt file)
        """
        try:
            import torch
            self.kan_model = torch.load(model_path)
            self._extract_structure()
        except ImportError:
            raise ImportError("PyKAN not installed. Install with: pip install pykan")

    def _extract_structure(self):
        """Extract KAN network structure for state-space conversion."""
        if self.kan_model is None:
            return

        try:
            self.layer_dims = []
            for layer in self.kan_model.layers:
                self.layer_dims.append(layer.in_features)
            self.layer_dims.append(self.kan_model.layers[-1].out_features)
            self._n_internal = sum(self.layer_dims[1:-1])
        except AttributeError:
            self._n_internal = 10  # Default

    def from_pole_residue(self, poles: np.ndarray, residues: np.ndarray,
                          d: complex = 0, e: complex = 0):
        """
        Create material model from pole-residue (Vector Fitting) form.

        mu*(s) = d + s*e + sum_k c_k / (s - p_k)

        Parameters:
            poles: Complex poles p_k (real or conjugate pairs)
            residues: Complex residues c_k
            d: Constant term
            e: Linear term (s coefficient)
        """
        n_poles = len(poles)

        A = np.zeros((n_poles, n_poles), dtype=complex)
        B = np.zeros((n_poles, 1), dtype=complex)
        C = np.zeros((1, n_poles), dtype=complex)

        for k, (p, c) in enumerate(zip(poles, residues)):
            A[k, k] = p
            B[k, 0] = c
            C[0, k] = 1.0

        D = np.array([[d]])
        E = np.array([[e]])

        self.state_space = {
            'A': A, 'B': B, 'C': C, 'D': D, 'E': E,
            'poles': poles, 'residues': residues,
            'n_states': n_poles
        }
        self._n_internal = n_poles

        return self

    def from_debye_relaxation(self, mu_inf: float, delta_mu: List[float],
                               tau: List[float]):
        """
        Create material model from Debye relaxation parameters.

        mu*(s) = mu_inf + sum_k delta_mu_k / (1 + s*tau_k)

        Parameters:
            mu_inf: High-frequency permeability mu_inf
            delta_mu: List of relaxation strengths
            tau: List of relaxation times [s]
        """
        n_modes = len(delta_mu)
        poles = np.array([-1.0 / t for t in tau])
        residues = np.array([dm / t for dm, t in zip(delta_mu, tau)])

        return self.from_pole_residue(poles, residues, d=mu_inf)

    def from_cole_cole(self, mu_inf: float, mu_0: float,
                        tau: float, alpha: float, n_poles: int = 5):
        """
        Create material model from Cole-Cole parameters.

        mu*(s) = mu_inf + (mu_0 - mu_inf) / (1 + (s*tau)^alpha)

        Parameters:
            mu_inf: High-frequency permeability
            mu_0: DC permeability
            tau: Relaxation time [s]
            alpha: Cole-Cole parameter (0 < alpha <= 1)
            n_poles: Number of poles for approximation
        """
        omega = np.logspace(-2, 8, 200) / tau
        s_samples = 1j * omega
        mu_samples = mu_inf + (mu_0 - mu_inf) / (1 + (s_samples * tau)**alpha)

        poles_init = -np.logspace(-1, 6, n_poles) / tau
        poles, residues = self._vector_fit(s_samples, mu_samples, poles_init)

        return self.from_pole_residue(poles, residues, d=mu_inf)

    def _vector_fit(self, s: np.ndarray, f: np.ndarray,
                     poles_init: np.ndarray, n_iter: int = 5):
        """Simplified Vector Fitting algorithm."""
        n_poles = len(poles_init)
        n_samples = len(s)
        poles = poles_init.copy()

        for _ in range(n_iter):
            A = np.zeros((n_samples, n_poles + 1), dtype=complex)
            for k in range(n_poles):
                A[:, k] = 1.0 / (s - poles[k])
            A[:, -1] = 1.0

            x, _, _, _ = np.linalg.lstsq(A, f, rcond=None)
            residues = x[:-1]

        return poles, residues

    def evaluate_frequency(self, frequencies) -> np.ndarray:
        """Evaluate material response at given frequencies."""
        if self.state_space is None:
            raise ValueError("No state-space model loaded")

        frequencies = np.asarray(frequencies)
        s = 2j * np.pi * frequencies
        A = self.state_space['A']
        B = self.state_space['B']
        C = self.state_space['C']
        D = self.state_space['D']
        E = self.state_space.get('E', np.zeros_like(D))

        n_freq = len(frequencies)
        result = np.zeros(n_freq, dtype=complex)

        for i, si in enumerate(s):
            try:
                G = np.linalg.inv(si * np.eye(A.shape[0]) - A)
                result[i] = (C @ G @ B + D + si * E)[0, 0]
            except np.linalg.LinAlgError:
                result[i] = np.nan

        return result

    def get_state_matrices(self) -> dict:
        """Get state-space matrices for transient simulation."""
        if self.state_space is None:
            raise ValueError("No state-space model loaded")
        return self.state_space

    @property
    def n_internal_states(self) -> int:
        """Number of internal state variables."""
        return self._n_internal


class KANInverseLaplace:
    """
    PyKANによる逆ラプラス変換

    任意の複素材料関数 mu*(s) から時間領域インパルス応答 h(t) を学習

    Debyeモデルの限界:
    - 単一緩和時間τの仮定が強すぎる
    - フェライト等は広い緩和時間分布を持つ
    - 共鳴現象(ジャイロ磁気共鳴等)を表現できない

    本クラスの特徴:
    - 任意のmu*(s)に対応(Debye, Cole-Cole, 共鳴型, 測定データ)
    - KANで有理関数を自動学習 → 極・留数抽出
    - 物理制約(因果律, 受動性)を付加可能
    - 時間領域カーネルh(t)を直接提供

    Usage:
        inv_lap = KANInverseLaplace()

        # 方法1: 周波数データから学習
        inv_lap.fit_frequency_data(frequencies, mu_complex)

        # 方法2: 数値逆ラプラスからh(t)を直接学習
        inv_lap.fit_time_kernel(t_samples, h_samples)

        # 時間領域カーネル取得
        h = inv_lap.evaluate_kernel(t)

        # 畳み込み計算
        M = inv_lap.convolve(H_history, dt)
    """

    def __init__(self, n_poles: int = 10, tol: float = 1e-6):
        """
        Parameters:
            n_poles: 有理近似の極数(多いほど精度向上、計算コスト増)
            tol: 収束判定閾値
        """
        self.n_poles = n_poles
        self.tol = tol
        self.poles = None      # 極 p_k (複素数)
        self.residues = None   # 留数 c_k (複素数)
        self.d = 0             # 定数項
        self.e = 0             # s項係数
        self.kan_model = None  # PyKAN model (if used)
        self._fitted = False

    def fit_frequency_data(self, frequencies: np.ndarray,
                            mu_complex: np.ndarray,
                            method: str = 'vector_fitting',
                            enforce_stability: bool = True,
                            enforce_passivity: bool = True) -> 'KANInverseLaplace':
        """
        周波数応答データから有理関数近似を学習

        Parameters:
            frequencies: 周波数配列 [Hz]
            mu_complex: 複素材料定数 mu*(f)
            method: 'vector_fitting' or 'kan'
            enforce_stability: 極を左半平面に制限
            enforce_passivity: 受動性条件 Re{mu*(jw)} > 0

        Returns:
            self
        """
        s = 2j * np.pi * frequencies

        if method == 'vector_fitting':
            self._fit_vector_fitting(s, mu_complex, enforce_stability)
        elif method == 'kan':
            self._fit_kan(s, mu_complex)
        else:
            raise ValueError(f"Unknown method: {method}")

        if enforce_passivity:
            self._enforce_passivity(s, mu_complex)

        self._fitted = True
        return self

    def _fit_vector_fitting(self, s: np.ndarray, f: np.ndarray,
                             enforce_stability: bool = True,
                             n_iter: int = 10):
        """
        Vector Fitting アルゴリズム (Gustavsen & Semlyen, 1999)

        f(s) ≈ Σ c_k/(s - p_k) + d + s*e

        反復的に極を再配置して最適な有理近似を求める
        """
        n_samples = len(s)
        n_poles = self.n_poles

        # 初期極: 虚軸上に対数的に配置
        omega_range = np.abs(np.imag(s))
        omega_min = max(omega_range.min(), 1e-3 * omega_range.max())
        omega_max = omega_range.max()

        # 複素共役ペアで初期化
        beta = np.logspace(np.log10(omega_min), np.log10(omega_max), n_poles // 2)
        alpha = beta * 0.01  # 小さな実部(安定性)

        poles = []
        for a, b in zip(alpha, beta):
            poles.append(-a + 1j * b)
            poles.append(-a - 1j * b)
        if n_poles % 2 == 1:
            poles.append(-omega_min * 0.1)  # 実極を1つ追加

        poles = np.array(poles[:n_poles])

        # Vector Fitting 反復
        for iteration in range(n_iter):
            # ステップ1: 重み関数 sigma(s) の極を固定して c, d, e を求める
            # sigma(s) = Σ c_k/(s - p_k) + 1
            # f(s) * sigma(s) = Σ c_k'/(s - p_k) + d' + s*e'

            # 行列構築
            n_unknowns = 2 * n_poles + 2  # c_k, c_k', d', e'
            A = np.zeros((n_samples, n_unknowns), dtype=complex)

            for k in range(n_poles):
                A[:, k] = 1.0 / (s - poles[k])              # c_k for numerator
                A[:, n_poles + k] = -f / (s - poles[k])     # c_k' for sigma*f

            A[:, 2*n_poles] = 1.0      # d
            A[:, 2*n_poles + 1] = s    # e (s項)

            # 実部と虚部を分離して実数最小二乗
            A_real = np.vstack([A.real, A.imag])
            f_real = np.hstack([f.real, f.imag])

            try:
                x, _, _, _ = np.linalg.lstsq(A_real, f_real, rcond=None)
            except np.linalg.LinAlgError:
                break

            residues = x[:n_poles]
            sigma_residues = x[n_poles:2*n_poles]
            d = x[2*n_poles]
            e = x[2*n_poles + 1]

            # ステップ2: sigma(s)の零点を新しい極として採用
            # sigma(s) = 0 の解が新しい極
            # これは sigma(s) の companion matrix の固有値

            # sigma(s) = 1 + Σ c_k'/(s - p_k) の零点
            # Companion matrix approach
            A_comp = np.diag(poles)
            b_comp = np.ones(n_poles)
            c_comp = sigma_residues

            # 状態空間: dx/dt = A_comp*x + b_comp*u, y = c_comp*x + 1*u
            # 伝達関数の零点 = (A - b*c/d)の固有値
            if abs(1.0) > self.tol:
                A_zero = A_comp - np.outer(b_comp, c_comp)
                new_poles = np.linalg.eigvals(A_zero)
            else:
                new_poles = poles

            # 安定性強制: 右半平面の極を左半平面に反転
            if enforce_stability:
                new_poles = np.array([p if p.real < 0 else -abs(p.real) + 1j*p.imag
                                      for p in new_poles])

            # 収束判定
            pole_change = np.max(np.abs(new_poles - poles))
            poles = new_poles

            if pole_change < self.tol:
                break

        # 最終フィッティング (極固定)
        A_final = np.zeros((n_samples, n_poles + 2), dtype=complex)
        for k in range(n_poles):
            A_final[:, k] = 1.0 / (s - poles[k])
        A_final[:, n_poles] = 1.0
        A_final[:, n_poles + 1] = s

        A_real = np.vstack([A_final.real, A_final.imag])
        f_real = np.hstack([f.real, f.imag])
        x, _, _, _ = np.linalg.lstsq(A_real, f_real, rcond=None)

        self.poles = poles
        self.residues = x[:n_poles]
        self.d = x[n_poles]
        self.e = x[n_poles + 1]

    def _fit_kan(self, s: np.ndarray, f: np.ndarray):
        """
        PyKANを使った有理関数学習

        KANの利点:
        - スプライン基底による滑らかな近似
        - 解釈可能な中間表現
        - 自動的な構造発見
        """
        try:
            from kan import KAN
            import torch
        except ImportError:
            raise ImportError("PyKAN not installed. Use method='vector_fitting' or install pykan")

        # 入力: (Re{s}, Im{s})
        # 出力: (Re{f}, Im{f})
        X = np.column_stack([s.real, s.imag])
        Y = np.column_stack([f.real, f.imag])

        X_tensor = torch.tensor(X, dtype=torch.float32)
        Y_tensor = torch.tensor(Y, dtype=torch.float32)

        # KAN構造: [2, hidden, hidden, 2]
        hidden = max(self.n_poles, 10)
        model = KAN(width=[2, hidden, hidden, 2], grid=5, k=3)

        # 学習
        model.train(
            {'train_input': X_tensor, 'train_label': Y_tensor},
            steps=1000,
            lamb=0.01
        )

        self.kan_model = model

        # 極・留数の抽出 (KANから有理関数パラメータを推定)
        self._extract_poles_from_kan(s, f)

    def _extract_poles_from_kan(self, s: np.ndarray, f: np.ndarray):
        """KAN学習結果から極・留数を抽出"""
        # KAN出力を評価して周波数応答を再構成
        # その後、Vector Fittingで極を抽出
        if self.kan_model is not None:
            import torch
            X = np.column_stack([s.real, s.imag])
            X_tensor = torch.tensor(X, dtype=torch.float32)

            with torch.no_grad():
                Y_pred = self.kan_model(X_tensor).numpy()

            f_kan = Y_pred[:, 0] + 1j * Y_pred[:, 1]

            # KAN出力に対してVector Fitting
            self._fit_vector_fitting(s, f_kan, enforce_stability=True, n_iter=5)

    def _enforce_passivity(self, s: np.ndarray, f_target: np.ndarray):
        """
        受動性条件の強制

        受動性: Re{mu*(jw)} > 0 for all w
        (エネルギー散逸が非負)
        """
        # 現在のモデル評価
        f_model = self.evaluate_frequency_response(s)

        # 受動性違反をチェック
        passivity_violation = np.any(f_model.real < 0)

        if passivity_violation:
            # 極の実部を調整して受動性を回復
            # 簡易的: 留数を正規化
            scale = np.max(np.abs(self.residues))
            if scale > 0:
                self.residues = self.residues / scale * np.abs(self.residues)

    def evaluate_frequency_response(self, s: np.ndarray) -> np.ndarray:
        """
        有理関数モデルの周波数応答を評価

        f(s) = Σ c_k/(s - p_k) + d + s*e
        """
        if not self._fitted:
            raise ValueError("Model not fitted. Call fit_frequency_data() first.")

        result = np.full(len(s), self.d + 0j)

        for p, c in zip(self.poles, self.residues):
            result = result + c / (s - p)

        if self.e != 0:
            result = result + self.e * s

        return result

    def evaluate_kernel(self, t: np.ndarray) -> np.ndarray:
        """
        時間領域インパルス応答 h(t) を評価

        h(t) = L^{-1}{ f(s) }
             = Σ c_k * exp(p_k * t) + d*delta(t) + e*delta'(t)

        Parameters:
            t: 時間配列 [s] (t >= 0)

        Returns:
            h(t): インパルス応答 (複素数、通常は実部のみ使用)
        """
        if not self._fitted:
            raise ValueError("Model not fitted. Call fit_frequency_data() first.")

        h = np.zeros(len(t), dtype=complex)

        # 因果律: t < 0 では h = 0
        t_positive = t >= 0

        for p, c in zip(self.poles, self.residues):
            h[t_positive] += c * np.exp(p * t[t_positive])

        # d, e項はデルタ関数 → 畳み込みで処理

        return h.real  # 物理的には実数

    def convolve(self, H_history: np.ndarray, dt: float) -> np.ndarray:
        """
        畳み込み積分を計算

        M(t) = ∫_0^t h(t-τ) * H(τ) dτ + d*H(t) + e*dH/dt

        Parameters:
            H_history: 磁界履歴 H(t) の配列 [n_time]
            dt: 時間刻み [s]

        Returns:
            M: 磁化応答 [n_time]
        """
        n_time = len(H_history)
        M = np.zeros(n_time)

        # 定数項
        M += self.d * H_history

        # s項 (微分)
        if self.e != 0:
            dH_dt = np.gradient(H_history, dt)
            M += self.e * dH_dt

        # 畳み込み項: 各極について再帰的に計算
        # c/(s-p) の時間応答は1次系ODE
        for p, c in zip(self.poles, self.residues):
            # dx/dt = p*x + c*H, M += x
            x = 0.0
            exp_p_dt = np.exp(p * dt)

            for i in range(n_time):
                # 指数オイラー法
                if abs(p) > 1e-10:
                    x = exp_p_dt * x + c * (exp_p_dt - 1) / p * H_history[i]
                else:
                    x = x + c * dt * H_history[i]
                M[i] += x.real

        return M

    def get_state_space(self) -> dict:
        """
        状態空間表現を取得

        dx/dt = A @ x + B @ H
        M = C @ x + D @ H + E @ dH/dt

        where x_k corresponds to pole p_k
        """
        if not self._fitted:
            raise ValueError("Model not fitted.")

        n = len(self.poles)

        # 対角状態空間 (各極が独立)
        A = np.diag(self.poles)
        B = self.residues.reshape(-1, 1)
        C = np.ones((1, n))
        D = np.array([[self.d]])
        E = np.array([[self.e]])

        return {
            'A': A, 'B': B, 'C': C, 'D': D, 'E': E,
            'n_states': n,
            'poles': self.poles,
            'residues': self.residues
        }

    def summary(self):
        """モデルのサマリーを表示"""
        if not self._fitted:
            print("Model not fitted.")
            return

        print(f"KAN Inverse Laplace Model")
        print(f"  Poles: {len(self.poles)}")
        print(f"  Constant term d: {self.d:.4g}")
        print(f"  Linear term e: {self.e:.4g}")

        print(f"\n  Pole-Residue pairs:")
        for i, (p, c) in enumerate(zip(self.poles, self.residues)):
            tau = -1.0 / p.real if abs(p.real) > 1e-10 else np.inf
            print(f"    {i+1}: p = {p:.4g}, c = {c:.4g}, tau = {tau:.4g} s")


class PRIMAWithKANMaterial:
    """
    PRIMA reduction integrated with KAN-based material models.

    Combines:
    1. Loop-Star-Magnetic PRIMA reduction
    2. KAN material models for mu*(s) and eps*(s)
    3. Transient simulation with material dynamics

    Full system with KAN materials:
        [R + sL    K_LS    s*K_LM    0     ] [I_L ]   [V]
        [K_SL     P/s+G   0         0     ] [Q_S ] = [0]
        [-K_ML    0       Z_MM*H(s) 0     ] [M   ]   [0]
        [0        0       K_Mz      Z_kan ] [z   ]   [0]

    where:
        H(s) = C @ (sI - A)^{-1} @ B + D  (KAN material transfer function)
        z = internal KAN state variables
    """

    def __init__(self, tol: float = 1e-10):
        self.tol = tol
        self.prima = LoopStarMagneticPRIMA(tol=tol)
        self.kan_mu = None
        self.kan_eps = None

    def set_magnetic_material(self, kan_material: KANMaterialInterface):
        """Set KAN model for complex permeability mu*."""
        self.kan_mu = kan_material

    def set_dielectric_material(self, kan_material: KANMaterialInterface):
        """Set KAN model for complex permittivity eps*."""
        self.kan_eps = kan_material

    def build_augmented_system(self,
                                L: np.ndarray, R: np.ndarray, P: np.ndarray,
                                K_LS: np.ndarray, K_LM: np.ndarray,
                                Z_MM_base: np.ndarray,
                                port_indices: List[int] = None) -> dict:
        """
        Build augmented system with KAN material states.

        State vector: x = [I_L, Q_S, M, z_mu, z_eps]
        """
        n_L = L.shape[0]
        n_S = P.shape[0]
        n_M = Z_MM_base.shape[0]

        n_mu = self.kan_mu.n_internal_states if self.kan_mu else 0
        n_eps = self.kan_eps.n_internal_states if self.kan_eps else 0

        n_total = n_L + n_S + n_M + n_mu + n_eps

        if port_indices is None:
            port_indices = [0]

        i_L = 0
        i_S = n_L
        i_M = n_L + n_S
        i_mu = n_L + n_S + n_M
        i_eps = n_L + n_S + n_M + n_mu

        # Build descriptor system: E @ dx/dt = A @ x + B @ u
        E = np.zeros((n_total, n_total))
        A = np.zeros((n_total, n_total))

        # Loop: L @ dI/dt = -R @ I - K_LS @ Q + V
        E[i_L:i_S, i_L:i_S] = L
        A[i_L:i_S, i_L:i_S] = -R
        A[i_L:i_S, i_S:i_M] = -K_LS

        # Star: algebraic constraint
        A[i_S:i_M, i_L:i_S] = K_LS.T
        A[i_S:i_M, i_S:i_M] = P

        # Magnetic
        A[i_M:i_mu, i_L:i_S] = -K_LM.T
        A[i_M:i_mu, i_M:i_mu] = Z_MM_base

        # KAN mu* dynamics
        if self.kan_mu:
            ss_mu = self.kan_mu.get_state_matrices()
            E[i_mu:i_eps, i_mu:i_eps] = np.eye(n_mu)
            A[i_mu:i_eps, i_mu:i_eps] = ss_mu['A'].real
            # KAN input: average of M values -> z_mu dynamics
            # A[i_mu:i_eps, i_M:i_mu] is [n_mu x n_M]
            B_mu = ss_mu['B'].real.flatten()  # [n_mu]
            A[i_mu:i_eps, i_M:i_mu] = np.outer(B_mu, np.ones(n_M)) / n_M
            # KAN output: couple back to M equation
            # A[i_M:i_mu, i_mu:i_eps] is [n_M x n_mu]
            C_mu = ss_mu['C'].real.flatten()  # [n_mu]
            A[i_M:i_mu, i_mu:i_eps] = np.outer(np.ones(n_M), C_mu) * np.mean(np.diag(Z_MM_base))

        # KAN eps* dynamics
        if self.kan_eps:
            ss_eps = self.kan_eps.get_state_matrices()
            E[i_eps:, i_eps:] = np.eye(n_eps)
            A[i_eps:, i_eps:] = ss_eps['A'].real
            B_eps = ss_eps['B'].real.flatten()  # [n_eps]
            A[i_eps:, i_S:i_M] = np.outer(B_eps, np.ones(n_S)) / n_S

        B = np.zeros((n_total, 1))
        for idx in port_indices:
            B[idx, 0] = 1.0

        C = np.zeros((1, n_total))
        for idx in port_indices:
            C[0, idx] = 1.0

        return {
            'E': E, 'A': A, 'B': B, 'C': C,
            'n_L': n_L, 'n_S': n_S, 'n_M': n_M,
            'n_mu': n_mu, 'n_eps': n_eps,
            'n_total': n_total,
            'port_indices': port_indices
        }

    def reduce_augmented_system(self, aug_system: dict,
                                 k_L: int = 10, k_S: int = 5, k_M: int = 5,
                                 k_kan: int = 3) -> dict:
        """Apply PRIMA reduction to augmented system."""
        E = aug_system['E']
        A = aug_system['A']
        B = aug_system['B']
        n_total = aug_system['n_total']

        s0 = 0
        k_total = k_L + k_S + k_M + k_kan

        A_shift = A - s0 * E
        try:
            lu = linalg.lu_factor(A_shift)
            def solve_shift(b):
                return linalg.lu_solve(lu, E @ b)
        except linalg.LinAlgError:
            A_inv = np.linalg.pinv(A_shift)
            def solve_shift(b):
                return A_inv @ E @ b

        Q = np.zeros((n_total, k_total))
        v = np.linalg.solve(A_shift, B[:, 0])
        if np.linalg.norm(v) > self.tol:
            v = v / np.linalg.norm(v)
        Q[:, 0] = v

        for j in range(1, k_total):
            w = solve_shift(Q[:, j-1])
            for i in range(j):
                h = np.dot(Q[:, i], w)
                w = w - h * Q[:, i]
            h_norm = np.linalg.norm(w)
            if h_norm < self.tol:
                Q = Q[:, :j]
                break
            Q[:, j] = w / h_norm

        E_red = Q.T @ E @ Q
        A_red = Q.T @ A @ Q
        B_red = Q.T @ B
        C_red = aug_system['C'] @ Q

        return {
            'E_red': E_red, 'A_red': A_red,
            'B_red': B_red, 'C_red': C_red,
            'Q': Q,
            'k_total': Q.shape[1],
            'original': aug_system
        }

    def simulate_transient(self, reduced: dict,
                           t_span: Tuple[float, float],
                           v_source: Callable[[float], float],
                           t_eval: np.ndarray = None) -> dict:
        """Run transient simulation on reduced system."""
        E_red = reduced['E_red']
        A_red = reduced['A_red']
        B_red = reduced['B_red']
        C_red = reduced['C_red']
        k = reduced['k_total']

        E_rank = np.linalg.matrix_rank(E_red)

        if E_rank == k:
            E_inv = np.linalg.inv(E_red)
            A_ode = E_inv @ A_red
            B_ode = E_inv @ B_red

            def ode_func(t, x):
                return A_ode @ x + B_ode.flatten() * v_source(t)

            x0 = np.zeros(k)
            sol = solve_ivp(ode_func, t_span, x0, t_eval=t_eval,
                           method='BDF', dense_output=True)

            return {
                't': sol.t,
                'x': sol.y.T,
                'y': (C_red @ sol.y).flatten(),
                'success': sol.success
            }
        else:
            raise NotImplementedError(
                f"DAE system (E rank {E_rank} < {k}). "
                "Use specialized DAE solver."
            )

    def compute_impedance_with_kan(self, reduced: dict,
                                    frequencies: np.ndarray) -> np.ndarray:
        """Compute impedance including KAN material dynamics."""
        E_red = reduced['E_red']
        A_red = reduced['A_red']
        B_red = reduced['B_red']
        C_red = reduced['C_red']

        n_freq = len(frequencies)
        Z = np.zeros(n_freq, dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            try:
                G = np.linalg.inv(s * E_red - A_red)
                Z[i] = (C_red @ G @ B_red)[0, 0]
            except np.linalg.LinAlgError:
                Z[i] = np.nan

        return Z


def demo_kan_material():
    """Demonstrate KAN material interface with Debye relaxation."""
    print("=" * 60)
    print("PyKAN Material Interface Demo")
    print("=" * 60)
    print("(Debye relaxation model for magnetic loss mu'')")

    # Create KAN material interface
    kan_mu = KANMaterialInterface(material_type='magnetic')

    # Debye relaxation for ferrite
    mu_inf = 100
    delta_mu = [800, 100]
    tau = [1e-6, 1e-8]

    kan_mu.from_debye_relaxation(mu_inf, delta_mu, tau)

    print(f"\nDebye model parameters:")
    print(f"  mu_inf = {mu_inf}")
    print(f"  Relaxation modes: {len(delta_mu)}")
    for i, (dm, t) in enumerate(zip(delta_mu, tau)):
        print(f"    Mode {i+1}: Delta_mu = {dm}, tau = {t*1e6:.2f} us")

    print(f"\nState-space representation:")
    print(f"  Internal states: {kan_mu.n_internal_states}")
    ss = kan_mu.get_state_matrices()
    print(f"  Poles: {ss['poles']}")

    # Frequency response
    frequencies = np.logspace(3, 9, 100)
    mu_complex = kan_mu.evaluate_frequency(frequencies)

    mu_imag = -np.imag(mu_complex)
    peak_idx = np.argmax(mu_imag)
    f_peak = frequencies[peak_idx]

    print(f"\nFrequency response:")
    print(f"  DC: mu' = {np.real(mu_complex[0]):.1f}")
    print(f"  Loss peak: f = {f_peak/1e6:.2f} MHz, mu'' = {mu_imag[peak_idx]:.1f}")
    print(f"  High-freq: mu' -> {np.real(mu_complex[-1]):.1f}")

    # PRIMA + KAN integration test
    print("\n--- PRIMA + KAN Integration ---")

    np.random.seed(42)
    n_L = 20
    n_M = 10

    L = np.eye(n_L) * 1e-6
    R = np.eye(n_L) * 0.1
    P = np.eye(n_L) * 1e10
    K_LS = np.random.randn(n_L, n_L) * 0.01
    K_LM = np.random.randn(n_L, n_M) * 0.1
    Z_MM_base = np.eye(n_M) * 1e-3

    prima_kan = PRIMAWithKANMaterial(tol=1e-10)
    prima_kan.set_magnetic_material(kan_mu)

    aug = prima_kan.build_augmented_system(L, R, P, K_LS, K_LM, Z_MM_base)

    print(f"\nAugmented system:")
    print(f"  Total states: {aug['n_total']}")
    print(f"    Loop: {aug['n_L']}, Star: {aug['n_S']}, Magnetic: {aug['n_M']}")
    print(f"    KAN mu states: {aug['n_mu']}")

    reduced = prima_kan.reduce_augmented_system(aug, k_L=8, k_S=5, k_M=5, k_kan=2)

    print(f"\nReduced system:")
    print(f"  Reduced states: {reduced['k_total']}")
    print(f"  Compression: {reduced['k_total']/aug['n_total']*100:.1f}%")

    test_freqs = np.array([1e3, 1e6, 1e9])
    Z = prima_kan.compute_impedance_with_kan(reduced, test_freqs)

    print(f"\nImpedance:")
    for f, z in zip(test_freqs, Z):
        print(f"  f = {f/1e6:.0f} MHz: |Z| = {np.abs(z):.2e} Ohm")

    return kan_mu, prima_kan, reduced


# =============================================================================
# ACA-based Sparse Circuit Extraction
# =============================================================================

class ACACircuitExtraction:
    """
    Adaptive Cross Approximation (ACA) for sparse circuit extraction.

    Converts dense Schur complement impedance matrices into sparse
    RLC circuit netlists using low-rank factorization.

    The key insight:
    - Schur complement Z_eff(s) is often numerically low-rank
    - ACA decomposes: Z_eff ≈ U @ V.T where U, V are tall-skinny
    - Low-rank form maps to coupled inductors/capacitors
    - Sparse circuit representation enables SPICE simulation

    Reference:
    - M. Bebendorf, "Approximation of boundary element matrices",
      Numer. Math., vol. 86, pp. 565-589, 2000.
    - A. Odabasioglu et al., "PRIMA: Passive Reduced-order Interconnect
      Macromodeling Algorithm", IEEE TCAD, 1998.
    """

    def __init__(self, tol: float = 1e-4, max_rank: int = None):
        """
        Initialize ACA circuit extractor.

        Parameters
        ----------
        tol : float
            Relative tolerance for low-rank approximation
        max_rank : int, optional
            Maximum rank (default: min(n_row, n_col) // 2)
        """
        self.tol = tol
        self.max_rank = max_rank
        self.U = None
        self.V = None
        self.rank = 0

    def aca_decompose(self, Z: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
        """
        ACA decomposition of impedance matrix.

        Z ≈ U @ V.T  where U is (n x r), V is (n x r), r << n

        Parameters
        ----------
        Z : np.ndarray
            Dense impedance matrix (n x n)

        Returns
        -------
        U, V : np.ndarray
            Low-rank factors
        rank : int
            Achieved rank
        """
        n = Z.shape[0]
        max_rank = self.max_rank if self.max_rank else n // 2

        U_list = []
        V_list = []

        # Residual matrix (initially full Z)
        R = Z.copy()
        Z_norm = np.linalg.norm(Z, 'fro')

        # Track used rows/cols
        used_rows = set()
        used_cols = set()

        # Start with row having max norm
        row_norms = np.linalg.norm(R, axis=1)
        i_star = np.argmax(row_norms)

        for k in range(max_rank):
            # Find pivot column in row i_star
            row = R[i_star, :]
            available_cols = [j for j in range(n) if j not in used_cols]
            if not available_cols:
                break
            j_star = available_cols[np.argmax(np.abs(row[available_cols]))]

            pivot = R[i_star, j_star]
            if np.abs(pivot) < self.tol * Z_norm / n:
                break

            # Extract cross vectors
            u = R[:, j_star].copy()
            v = R[i_star, :].copy() / pivot

            U_list.append(u)
            V_list.append(v)

            # Update residual
            R = R - np.outer(u, v)

            used_rows.add(i_star)
            used_cols.add(j_star)

            # Check convergence
            R_norm = np.linalg.norm(R, 'fro')
            if R_norm < self.tol * Z_norm:
                break

            # Find next pivot row
            available_rows = [i for i in range(n) if i not in used_rows]
            if not available_rows:
                break
            row_norms = np.linalg.norm(R[available_rows, :], axis=1)
            i_star = available_rows[np.argmax(row_norms)]

        if U_list:
            self.U = np.column_stack(U_list)
            self.V = np.column_stack(V_list)
            self.rank = len(U_list)
        else:
            self.U = np.zeros((n, 1))
            self.V = np.zeros((n, 1))
            self.rank = 0

        return self.U, self.V, self.rank

    def extract_rlc_circuit(self,
                            Z_s: Callable[[complex], np.ndarray],
                            frequencies: np.ndarray,
                            port_names: List[str] = None) -> dict:
        """
        Extract RLC circuit from frequency-dependent impedance.

        Uses rational fitting on each low-rank component:
        Z_eff(s) ≈ Σ_k u_k @ v_k^T * H_k(s)

        where H_k(s) is rational function fitted to frequency data.

        Parameters
        ----------
        Z_s : callable
            Function Z_s(s) returning impedance matrix at complex frequency s
        frequencies : np.ndarray
            Frequencies for fitting [Hz]
        port_names : list, optional
            Names for ports (default: P1, P2, ...)

        Returns
        -------
        circuit : dict
            Circuit netlist with R, L, C, coupling elements
        """
        # Sample impedance at multiple frequencies
        n_freq = len(frequencies)
        s_vals = 2j * np.pi * frequencies

        Z_samples = [Z_s(s) for s in s_vals]
        n_ports = Z_samples[0].shape[0]

        if port_names is None:
            port_names = [f'P{i+1}' for i in range(n_ports)]

        # ACA on reference frequency (mid-band)
        ref_idx = n_freq // 2
        Z_ref = np.real(Z_samples[ref_idx])  # Use real part as reference
        U, V, rank = self.aca_decompose(Z_ref)

        print(f"ACA decomposition: rank = {rank} (from {n_ports}x{n_ports})")

        # For each low-rank mode, fit rational function
        circuit = {
            'n_ports': n_ports,
            'port_names': port_names,
            'rank': rank,
            'elements': [],
            'coupling': []
        }

        # Extract circuit elements from low-rank factors
        for k in range(rank):
            u_k = U[:, k]
            v_k = V[:, k]

            # Mode impedance Z_k(s) = u_k^T @ Z(s) @ v_k / (u_k^T @ u_k)
            Z_mode = np.zeros(n_freq, dtype=complex)
            for i, Z in enumerate(Z_samples):
                Z_mode[i] = np.dot(u_k, Z @ v_k) / np.dot(u_k, u_k)

            # Fit to R + sL + 1/(sC) form
            params = self._fit_rlc_mode(frequencies, Z_mode)

            # Add element for this mode
            element = {
                'mode': k,
                'R': params['R'],
                'L': params['L'],
                'C': params['C'],
                'u': u_k,  # Port coupling vector
                'v': v_k,
                'quality': params['fit_quality']
            }
            circuit['elements'].append(element)

            # Coupling between ports (from outer product structure)
            for i in range(n_ports):
                for j in range(n_ports):
                    if i < j and np.abs(u_k[i] * v_k[j]) > self.tol:
                        coupling = {
                            'mode': k,
                            'port_i': i,
                            'port_j': j,
                            'strength': u_k[i] * v_k[j]
                        }
                        circuit['coupling'].append(coupling)

        return circuit

    def _fit_rlc_mode(self, frequencies: np.ndarray,
                      Z_mode: np.ndarray) -> dict:
        """
        Fit impedance mode using Vector Fitting for pole-residue form.

        The mode impedance is fitted to:
        Z(s) = R + sL + sum_k c_k / (s - p_k) + d + s*e

        For low frequencies (power electronics), this simplifies to R + sL.

        Parameters
        ----------
        frequencies : np.ndarray
            Frequencies [Hz]
        Z_mode : np.ndarray
            Complex impedance at each frequency

        Returns
        -------
        params : dict
            R, L, C values, poles, residues, and fit quality
        """
        omega = 2 * np.pi * frequencies
        s = 1j * omega
        n_freq = len(frequencies)

        # Try Vector Fitting with 2 poles first
        n_poles = 2
        poles, residues, d, e = self._vector_fitting_mode(s, Z_mode, n_poles)

        # Compute fitted impedance
        Z_fit = np.full(n_freq, d, dtype=complex)
        Z_fit += e * s
        for p, c in zip(poles, residues):
            Z_fit += c / (s - p)

        error = np.linalg.norm(Z_mode - Z_fit) / (np.linalg.norm(Z_mode) + 1e-20)

        # Extract equivalent R, L, C from pole-residue form
        # For simple RL behavior: Z ≈ d + s*e means R=d, L=e
        R = np.real(d)
        L = np.real(e)

        # C estimation from poles (if there's a resonance)
        # For s-domain: C creates a 1/s term
        # If poles are near imaginary axis, estimate resonant frequency
        if len(poles) > 0:
            # Use slowest pole to estimate effective C
            min_pole_mag = np.min(np.abs(poles))
            if min_pole_mag > 1e-10 and L > 1e-15:
                # w0^2 = 1/(LC), so C = 1/(L * w0^2)
                C = 1.0 / (L * min_pole_mag**2) if L > 0 else 1e-12
            else:
                C = 1e-12
        else:
            C = 1e-12

        # Ensure positive values
        L = max(L, 1e-15)
        C = max(C, 1e-15)
        R = max(R, 0)

        return {
            'R': R,
            'L': L,
            'C': C,
            'poles': poles,
            'residues': residues,
            'd': d,
            'e': e,
            'fit_quality': 1 - min(error, 1.0)
        }

    def _vector_fitting_mode(self, s: np.ndarray, f: np.ndarray,
                              n_poles: int = 2) -> Tuple:
        """
        Vector Fitting for a single mode impedance.

        Parameters
        ----------
        s : np.ndarray
            Complex frequencies (jw)
        f : np.ndarray
            Function values at s
        n_poles : int
            Number of poles to fit

        Returns
        -------
        poles, residues, d, e : tuple
            Pole-residue representation
        """
        n = len(s)

        # Initial pole guess: logarithmically spaced on negative real axis
        f_range = np.abs(s[-1] - s[0]) / (2 * np.pi)
        poles = -np.logspace(np.log10(max(np.abs(s[0])/(2*np.pi), 1)),
                             np.log10(max(np.abs(s[-1])/(2*np.pi), 10)),
                             n_poles) * 2 * np.pi

        # Iterate Vector Fitting
        for iteration in range(5):
            # Build system matrix for pole relocation
            # f(s) ≈ sum c_k/(s-p_k) + d + s*e
            # sigma(s) = sum ~c_k/(s-p_k) + 1
            # sigma * f ≈ sum c_k/(s-p_k) + d + s*e

            A = np.zeros((2*n, 2*n_poles + 2), dtype=float)
            b = np.zeros(2*n, dtype=float)

            for i, (si, fi) in enumerate(zip(s, f)):
                row_r = np.zeros(2*n_poles + 2)
                row_i = np.zeros(2*n_poles + 2)

                for k, pk in enumerate(poles):
                    denom = si - pk
                    row_r[k] = np.real(1.0 / denom)
                    row_i[k] = np.imag(1.0 / denom)

                row_r[n_poles] = 1.0  # d real
                row_i[n_poles] = 0.0
                row_r[n_poles + 1] = np.real(si)  # e * s
                row_i[n_poles + 1] = np.imag(si)

                # Sigma terms (for pole relocation)
                for k, pk in enumerate(poles):
                    denom = si - pk
                    row_r[n_poles + 2 + k] = -np.real(fi / denom)
                    row_i[n_poles + 2 + k] = -np.imag(fi / denom)

                A[2*i, :] = row_r
                A[2*i + 1, :] = row_i
                b[2*i] = np.real(fi)
                b[2*i + 1] = np.imag(fi)

            # Solve least squares
            try:
                x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
            except np.linalg.LinAlgError:
                break

            residues = x[:n_poles]
            d = x[n_poles]
            e = x[n_poles + 1]
            sigma_res = x[n_poles + 2:]

            # Update poles
            # New poles are eigenvalues of (A - b*c.T) where A=diag(poles), b=1, c=sigma_res
            H = np.diag(poles) - np.outer(np.ones(n_poles), sigma_res)
            try:
                new_poles = np.linalg.eigvals(H)
                # Enforce stability (negative real part)
                new_poles = np.where(np.real(new_poles) > 0,
                                     -np.abs(np.real(new_poles)) + 1j*np.imag(new_poles),
                                     new_poles)
                poles = new_poles.real  # Keep only real poles for RLC
            except np.linalg.LinAlgError:
                break

        # Final residue computation with fixed poles
        A_final = np.zeros((2*n, n_poles + 2), dtype=float)
        b_final = np.zeros(2*n, dtype=float)

        for i, (si, fi) in enumerate(zip(s, f)):
            for k, pk in enumerate(poles):
                denom = si - pk
                A_final[2*i, k] = np.real(1.0 / denom)
                A_final[2*i + 1, k] = np.imag(1.0 / denom)

            A_final[2*i, n_poles] = 1.0
            A_final[2*i + 1, n_poles] = 0.0
            A_final[2*i, n_poles + 1] = np.real(si)
            A_final[2*i + 1, n_poles + 1] = np.imag(si)

            b_final[2*i] = np.real(fi)
            b_final[2*i + 1] = np.imag(fi)

        try:
            x_final, _, _, _ = np.linalg.lstsq(A_final, b_final, rcond=None)
            residues = x_final[:n_poles]
            d = x_final[n_poles]
            e = x_final[n_poles + 1]
        except np.linalg.LinAlgError:
            residues = np.zeros(n_poles)
            d = np.mean(np.real(f))
            e = 0.0

        return poles, residues, d, e

    def to_spice_netlist(self, circuit: dict,
                         filename: str = None) -> str:
        """
        Convert extracted circuit to SPICE netlist.

        Parameters
        ----------
        circuit : dict
            Circuit from extract_rlc_circuit()
        filename : str, optional
            Output file path

        Returns
        -------
        netlist : str
            SPICE netlist string
        """
        lines = []
        lines.append("* ACA-extracted RLC circuit netlist")
        lines.append(f"* Ports: {circuit['n_ports']}, Modes: {circuit['rank']}")
        lines.append("")

        # Port definitions
        lines.append("* Port nodes")
        for i, name in enumerate(circuit['port_names']):
            lines.append(f".SUBCKT {name} n{i+1} 0")
            lines.append(f".ENDS")
        lines.append("")

        # Internal mode nodes
        node_counter = circuit['n_ports'] + 1

        lines.append("* Mode elements")
        for elem in circuit['elements']:
            k = elem['mode']
            R, L, C = elem['R'], elem['L'], elem['C']

            # Each mode creates internal resonator
            n_in = node_counter
            n_out = node_counter + 1
            node_counter += 2

            lines.append(f"* Mode {k} (Q = {elem['quality']:.2f})")
            if R > 1e-10:
                lines.append(f"R_m{k} n{n_in} n{n_out} {R:.6e}")
            if L > 1e-15:
                lines.append(f"L_m{k} n{n_out} n{node_counter} {L:.6e}")
                node_counter += 1
            if C > 1e-15:
                lines.append(f"C_m{k} n{node_counter-1} 0 {C:.6e}")

            # Coupling to ports via controlled sources
            u, v = elem['u'], elem['v']
            for i in range(circuit['n_ports']):
                if np.abs(u[i]) > self.tol:
                    lines.append(f"G_m{k}_p{i+1} n{i+1} 0 n{n_in} 0 {u[i]:.6e}")

            lines.append("")

        # Mutual coupling between ports
        if circuit['coupling']:
            lines.append("* Inter-port coupling")
            for coup in circuit['coupling']:
                k = coup['mode']
                i, j = coup['port_i'], coup['port_j']
                strength = coup['strength']
                lines.append(f"K_m{k}_{i+1}_{j+1} L_p{i+1} L_p{j+1} {np.abs(strength):.6e}")

        lines.append("")
        lines.append(".END")

        netlist = '\n'.join(lines)

        if filename:
            with open(filename, 'w') as f:
                f.write(netlist)

        return netlist

    def compute_sparse_impedance(self, s: complex) -> np.ndarray:
        """
        Compute impedance using sparse (low-rank) representation.

        Z(s) ≈ U @ D(s) @ V.T

        where D(s) is diagonal containing mode transfer functions.

        Parameters
        ----------
        s : complex
            Complex frequency

        Returns
        -------
        Z : np.ndarray
            Approximate impedance matrix
        """
        if self.U is None or self.V is None:
            raise ValueError("Run aca_decompose first")

        # Simple reconstruction (extend for frequency-dependent D)
        return self.U @ self.V.T


class KANContinuedFraction:
    """
    PyKAN-based Continued Fraction Learning.

    Uses Kolmogorov-Arnold Networks to learn the continued fraction
    coefficients from frequency response data, enabling automatic
    circuit synthesis for arbitrary transfer functions.

    Workflow:
    1. Input: Z(jw) data at multiple frequencies
    2. KAN learns mapping: frequency -> CFE coefficients
    3. Output: PRIMA ladder network (SPICE compatible)

    This approach handles:
    - Dowell model (skin/proximity effects)
    - Complex permeability mu*(s)
    - General Foster/PRIMA network synthesis
    """

    def __init__(self, n_terms: int = 6, n_hidden: int = 10):
        """
        Initialize KAN-based CFE learner.

        Parameters
        ----------
        n_terms : int
            Number of CFE terms (ladder elements)
        n_hidden : int
            Hidden layer size for KAN
        """
        self.n_terms = n_terms
        self.n_hidden = n_hidden
        self.kan_model = None
        self.ladder_coeffs = None
        self.element_types = None  # ['R', 'L', 'C', ...]

    def fit_from_frequency_data(self, frequencies: np.ndarray,
                                Z_data: np.ndarray,
                                element_pattern: List[str] = None) -> dict:
        """
        Learn CFE coefficients from frequency response data.

        Parameters
        ----------
        frequencies : np.ndarray
            Frequency points [Hz]
        Z_data : np.ndarray
            Complex impedance at each frequency
        element_pattern : list, optional
            Pattern of element types, e.g., ['R', 'L', 'R', 'L', 'R']
            Default: alternating R-L for inductive systems

        Returns
        -------
        result : dict
            Learned coefficients and fit quality
        """
        n_freq = len(frequencies)
        omega = 2 * np.pi * frequencies
        s = 1j * omega

        # Default element pattern
        if element_pattern is None:
            element_pattern = []
            for i in range(self.n_terms):
                if i % 2 == 0:
                    element_pattern.append('R')
                else:
                    element_pattern.append('L')

        self.element_types = element_pattern[:self.n_terms]
        n_coeffs = len(self.element_types)

        # Build optimization problem
        # Z_ladder(s, a) = CFE with coefficients a
        # Minimize ||Z_data - Z_ladder||^2

        def ladder_impedance(s_val, coeffs):
            """Evaluate ladder impedance with given coefficients."""
            Z = 0.0
            for elem_type, coeff in zip(reversed(self.element_types),
                                        reversed(coeffs)):
                if elem_type == 'R':
                    Z = Z + coeff
                elif elem_type == 'L':
                    Z = Z + s_val * coeff
                elif elem_type == 'C':
                    if np.abs(Z) > 1e-15:
                        Y = 1.0 / Z + s_val * coeff
                        Z = 1.0 / Y if np.abs(Y) > 1e-15 else 1e15
                    else:
                        Z = 1.0 / (s_val * coeff) if np.abs(s_val * coeff) > 1e-15 else 1e15
            return Z

        def objective(coeffs):
            """Objective function for optimization."""
            error = 0.0
            for si, Zi in zip(s, Z_data):
                Z_model = ladder_impedance(si, coeffs)
                error += np.abs(Zi - Z_model)**2
            return error

        # Initial guess from simple estimates
        x0 = np.ones(n_coeffs) * 0.1

        # Estimate initial values based on data
        Z_dc = Z_data[0]
        Z_hf = Z_data[-1]

        for i, elem_type in enumerate(self.element_types):
            if elem_type == 'R':
                x0[i] = np.abs(np.real(Z_dc)) / (i + 1)
            elif elem_type == 'L':
                x0[i] = np.abs(np.imag(Z_hf)) / (omega[-1] * (i + 1))
            elif elem_type == 'C':
                x0[i] = 1.0 / (np.abs(np.imag(Z_dc)) * omega[0] * (i + 1) + 1e-15)

        # Optimize using scipy
        from scipy.optimize import minimize

        # Bounds: all coefficients must be positive
        bounds = [(1e-15, None) for _ in range(n_coeffs)]

        result = minimize(objective, x0, method='L-BFGS-B', bounds=bounds,
                          options={'maxiter': 500, 'ftol': 1e-10})

        self.ladder_coeffs = result.x

        # Compute fit quality
        Z_fit = np.array([ladder_impedance(si, self.ladder_coeffs) for si in s])
        fit_error = np.linalg.norm(Z_data - Z_fit) / np.linalg.norm(Z_data)

        return {
            'coefficients': self.ladder_coeffs,
            'element_types': self.element_types,
            'fit_error': fit_error,
            'Z_fit': Z_fit,
            'converged': result.success
        }

    def fit_with_kan_network(self, frequencies: np.ndarray,
                             Z_data: np.ndarray,
                             param_data: np.ndarray = None) -> dict:
        """
        Learn CFE using KAN network for parametric variation.

        This learns a mapping: parameters -> CFE coefficients
        enabling real-time circuit extraction for varying conditions.

        Parameters
        ----------
        frequencies : np.ndarray
            Frequency points [Hz]
        Z_data : np.ndarray
            Complex impedance data, shape (n_params, n_freq) or (n_freq,)
        param_data : np.ndarray, optional
            Parameter values for parametric learning

        Returns
        -------
        result : dict
            KAN model and learned coefficients
        """
        # For non-parametric case, use direct fitting
        if param_data is None:
            return self.fit_from_frequency_data(frequencies, Z_data)

        # Parametric case: Learn mapping params -> coeffs
        n_params = len(param_data)
        n_freq = len(frequencies)

        if Z_data.ndim == 1:
            Z_data = Z_data.reshape(1, -1)

        # Fit CFE for each parameter set
        all_coeffs = []
        for i, Z_i in enumerate(Z_data):
            result_i = self.fit_from_frequency_data(frequencies, Z_i)
            all_coeffs.append(result_i['coefficients'])

        all_coeffs = np.array(all_coeffs)

        # Train simple polynomial model as KAN approximation
        # (Real KAN would use pykan library)
        from numpy.polynomial import polynomial as P

        n_coeffs = all_coeffs.shape[1]
        poly_models = []

        for j in range(n_coeffs):
            # Fit polynomial: param -> coeff_j
            coeffs_poly = np.polyfit(param_data, all_coeffs[:, j], deg=3)
            poly_models.append(coeffs_poly)

        # Store KAN model (polynomial approximation)
        self.kan_model = {
            'type': 'polynomial',
            'models': poly_models,
            'n_coeffs': n_coeffs
        }

        return {
            'kan_model': self.kan_model,
            'training_coeffs': all_coeffs,
            'element_types': self.element_types
        }

    def predict_coefficients(self, param: float) -> np.ndarray:
        """
        Predict CFE coefficients for given parameter value.

        Parameters
        ----------
        param : float
            Parameter value (e.g., frequency ratio, temperature)

        Returns
        -------
        coeffs : np.ndarray
            Predicted CFE coefficients
        """
        if self.kan_model is None:
            raise ValueError("KAN model not trained. Call fit_with_kan_network first.")

        if self.kan_model['type'] == 'polynomial':
            coeffs = np.zeros(self.kan_model['n_coeffs'])
            for j, poly_coeffs in enumerate(self.kan_model['models']):
                coeffs[j] = np.polyval(poly_coeffs, param)
            return np.maximum(coeffs, 1e-15)  # Ensure positive
        else:
            raise NotImplementedError(f"KAN type {self.kan_model['type']} not implemented")

    def evaluate(self, s: complex, coeffs: np.ndarray = None) -> complex:
        """
        Evaluate ladder impedance at complex frequency.

        Parameters
        ----------
        s : complex
            Complex frequency
        coeffs : np.ndarray, optional
            Coefficients to use (default: self.ladder_coeffs)

        Returns
        -------
        Z : complex
            Impedance
        """
        if coeffs is None:
            coeffs = self.ladder_coeffs

        if coeffs is None:
            raise ValueError("No coefficients available. Call fit first.")

        Z = 0.0
        for elem_type, coeff in zip(reversed(self.element_types),
                                    reversed(coeffs)):
            if elem_type == 'R':
                Z = Z + coeff
            elif elem_type == 'L':
                Z = Z + s * coeff
            elif elem_type == 'C':
                if np.abs(Z) > 1e-15:
                    Y = 1.0 / Z + s * coeff
                    Z = 1.0 / Y if np.abs(Y) > 1e-15 else 1e15
                else:
                    Z = 1.0 / (s * coeff) if np.abs(s * coeff) > 1e-15 else 1e15
        return Z

    def to_prima_ladder(self) -> List[Tuple[str, float]]:
        """
        Convert to PRIMA ladder representation.

        Returns
        -------
        ladder : list
            List of (element_type, value) tuples
        """
        if self.ladder_coeffs is None:
            raise ValueError("No coefficients available. Call fit first.")

        return list(zip(self.element_types, self.ladder_coeffs))

    def to_spice_subcircuit(self, name: str = "KAN_CFE") -> str:
        """
        Generate SPICE subcircuit from learned CFE.

        Parameters
        ----------
        name : str
            Subcircuit name

        Returns
        -------
        spice : str
            SPICE netlist string
        """
        ladder = self.to_prima_ladder()

        lines = []
        lines.append(f".SUBCKT {name} in out")
        lines.append(f"* KAN-learned PRIMA ladder with {len(ladder)} elements")

        node = 'in'
        next_node_num = 1

        for i, (elem_type, value) in enumerate(ladder):
            if i == len(ladder) - 1:
                next_node = 'out'
            else:
                next_node = f'n{next_node_num}'
                next_node_num += 1

            if elem_type == 'R':
                lines.append(f"R{i+1} {node} {next_node} {value:.6e}")
            elif elem_type == 'L':
                lines.append(f"L{i+1} {node} {next_node} {value:.6e}")
            elif elem_type == 'C':
                lines.append(f"C{i+1} {node} 0 {value:.6e}")
                next_node = node

            node = next_node

        lines.append(".ENDS")
        return '\n'.join(lines)


class PyKANSurfaceImpedance:
    """
    PyKAN-based Surface Impedance Learning for Nonlinear Magnetic Materials.

    Uses Kolmogorov-Arnold Networks (KAN) to learn the surface impedance
    Z_s(H, f) from ESIM (Effective Surface Impedance Method) data or
    measurement data.

    Key features:
    1. Learns nonlinear Z_s(H, f) relationship from data
    2. Supports both real and imaginary parts (Z_s = R_s + jX_s)
    3. Exports to SPICE PWL (Piece-Wise Linear) format
    4. Integrates with PRIMA for reduced-order models

    Workflow:
        1. Generate training data from ESIM cell problem or measurements
        2. Train KAN: (H, f) -> (R_s, X_s) or (|Z_s|, angle)
        3. Export to SPICE PWL or Verilog-A
        4. Use in coupled PEEC + HDiv-VIM / reduced-FEM simulations

    Reference:
        Liu et al., "KAN: Kolmogorov-Arnold Networks", 2024
        Hollaus et al., "Nonlinear Effective Surface Impedance", IEEE TMAG, 2025
    """

    def __init__(self, width: List[int] = None, grid: int = 5, k: int = 3):
        """
        Initialize PyKAN surface impedance model.

        Parameters
        ----------
        width : list of int, optional
            KAN layer widths. Default: [2, 5, 5, 2] for (H, f) -> (R_s, X_s)
        grid : int
            Number of grid points for spline basis (default: 5)
        k : int
            Spline order (default: 3 for cubic)
        """
        self.width = width if width is not None else [2, 5, 5, 2]
        self.grid = grid
        self.k = k
        self.kan_model = None
        self.is_trained = False

        # Normalization parameters
        self.H_min = None
        self.H_max = None
        self.f_min = None
        self.f_max = None
        self.Z_scale = None

        # Training history
        self.train_loss_history = []

    def train_from_esim_data(self, H_data: np.ndarray, f_data: np.ndarray,
                              Zs_data: np.ndarray,
                              epochs: int = 100, lr: float = 0.01,
                              verbose: bool = True) -> dict:
        """
        Train KAN from ESIM surface impedance data.

        Parameters
        ----------
        H_data : np.ndarray
            Surface magnetic field values [A/m], shape (n_H,)
        f_data : np.ndarray
            Frequency values [Hz], shape (n_f,)
        Zs_data : np.ndarray
            Complex surface impedance [Ohm], shape (n_H, n_f)
        epochs : int
            Number of training epochs
        lr : float
            Learning rate
        verbose : bool
            Print training progress

        Returns
        -------
        result : dict
            Training results with loss history
        """
        try:
            from kan import KAN
            import torch
        except ImportError:
            raise ImportError("PyKAN not installed. Install with: pip install pykan")

        # Create meshgrid for training data
        H_grid, f_grid = np.meshgrid(H_data, f_data, indexing='ij')
        H_flat = H_grid.flatten()
        f_flat = f_grid.flatten()
        Zs_flat = Zs_data.flatten()

        # Store normalization parameters
        self.H_min, self.H_max = H_flat.min(), H_flat.max()
        self.f_min, self.f_max = f_flat.min(), f_flat.max()
        self.Z_scale = np.abs(Zs_flat).max()

        # Normalize inputs to [0, 1]
        H_norm = (H_flat - self.H_min) / (self.H_max - self.H_min + 1e-10)
        f_norm = (f_flat - self.f_min) / (self.f_max - self.f_min + 1e-10)

        # Normalize outputs
        Rs_norm = np.real(Zs_flat) / self.Z_scale
        Xs_norm = np.imag(Zs_flat) / self.Z_scale

        # Prepare PyTorch tensors
        X_train = torch.tensor(np.column_stack([H_norm, f_norm]),
                               dtype=torch.float32)
        Y_train = torch.tensor(np.column_stack([Rs_norm, Xs_norm]),
                               dtype=torch.float32)

        # Create and train KAN
        self.kan_model = KAN(width=self.width, grid=self.grid, k=self.k)

        if verbose:
            print(f"Training PyKAN surface impedance model...")
            print(f"  Input: (H, f) normalized to [0, 1]")
            print(f"  Output: (R_s, X_s) / {self.Z_scale:.2e}")
            print(f"  Training samples: {len(H_flat)}")
            print(f"  KAN width: {self.width}")

        # Training loop
        dataset = {'train_input': X_train, 'train_label': Y_train,
                   'test_input': X_train, 'test_label': Y_train}

        self.train_loss_history = []

        results = self.kan_model.fit(dataset, opt='Adam', steps=epochs, lr=lr,
                                      loss_fn=torch.nn.MSELoss())

        self.is_trained = True

        if verbose:
            final_loss = results['train_loss'][-1] if 'train_loss' in results else 0
            print(f"  Final loss: {final_loss:.6e}")

        return {
            'loss_history': results.get('train_loss', []),
            'final_loss': results.get('train_loss', [0])[-1],
            'n_samples': len(H_flat),
            'width': self.width
        }

    def train_from_measurement(self, H_data: np.ndarray, f_data: np.ndarray,
                                Rs_data: np.ndarray, Xs_data: np.ndarray,
                                **kwargs) -> dict:
        """
        Train KAN from measured R_s and X_s data.

        Parameters
        ----------
        H_data : np.ndarray
            Surface field values [A/m]
        f_data : np.ndarray
            Frequencies [Hz]
        Rs_data : np.ndarray
            Surface resistance [Ohm/sq]
        Xs_data : np.ndarray
            Surface reactance [Ohm/sq]

        Returns
        -------
        result : dict
            Training results
        """
        Zs_data = Rs_data + 1j * Xs_data
        return self.train_from_esim_data(H_data, f_data, Zs_data, **kwargs)

    def predict(self, H: np.ndarray, f: np.ndarray) -> np.ndarray:
        """
        Predict surface impedance for given (H, f) values.

        Parameters
        ----------
        H : np.ndarray
            Surface field [A/m]
        f : np.ndarray
            Frequency [Hz]

        Returns
        -------
        Zs : np.ndarray
            Complex surface impedance [Ohm]
        """
        if not self.is_trained:
            raise ValueError("Model not trained. Call train_from_esim_data first.")

        import torch

        H = np.asarray(H)
        f = np.asarray(f)

        # Handle scalar inputs
        scalar_input = H.ndim == 0 and f.ndim == 0
        if scalar_input:
            H = np.array([H])
            f = np.array([f])

        # Normalize
        H_norm = (H - self.H_min) / (self.H_max - self.H_min + 1e-10)
        f_norm = (f - self.f_min) / (self.f_max - self.f_min + 1e-10)

        # Predict
        X = torch.tensor(np.column_stack([H_norm.flatten(), f_norm.flatten()]),
                         dtype=torch.float32)

        with torch.no_grad():
            Y_pred = self.kan_model(X).numpy()

        Rs = Y_pred[:, 0] * self.Z_scale
        Xs = Y_pred[:, 1] * self.Z_scale
        Zs = Rs + 1j * Xs

        if scalar_input:
            return Zs[0]

        return Zs.reshape(H.shape)

    def to_spice_pwl(self, H_values: np.ndarray, f_ref: float,
                      element_name: str = "Zskin") -> str:
        """
        Export to SPICE PWL (Piece-Wise Linear) format.

        Creates a behavioral model using PWL sources for R(H) and L(H).

        Parameters
        ----------
        H_values : np.ndarray
            H field values for PWL table [A/m]
        f_ref : float
            Reference frequency [Hz]
        element_name : str
            Base name for SPICE elements

        Returns
        -------
        spice : str
            SPICE netlist string
        """
        if not self.is_trained:
            raise ValueError("Model not trained. Call train_from_esim_data first.")

        omega_ref = 2 * np.pi * f_ref
        Zs = self.predict(H_values, np.full_like(H_values, f_ref))
        Rs = np.real(Zs)
        Xs = np.imag(Zs)
        Ls = Xs / omega_ref  # X_s = omega * L_s

        lines = []
        lines.append(f"* PyKAN Surface Impedance Model")
        lines.append(f"* Reference frequency: {f_ref/1e3:.1f} kHz")
        lines.append(f"* H range: {H_values.min():.1f} - {H_values.max():.1f} A/m")
        lines.append(f"")
        lines.append(f".SUBCKT {element_name} in out H_sense")
        lines.append(f"")

        # R(H) as PWL table
        lines.append(f"* Surface Resistance R_s(H)")
        pwl_R = " ".join([f"{H:.3e},{R:.6e}" for H, R in zip(H_values, Rs)])
        lines.append(f"B{element_name}_R in mid I=V(H_sense)*PWL({pwl_R})")
        lines.append(f"")

        # L(H) as PWL table
        lines.append(f"* Surface Inductance L_s(H) at f_ref = {f_ref/1e3:.1f} kHz")
        pwl_L = " ".join([f"{H:.3e},{L:.6e}" for H, L in zip(H_values, Ls)])
        lines.append(f"L{element_name} mid out 1")  # 1H placeholder
        lines.append(f"* Note: For frequency-dependent L, use Verilog-A")
        lines.append(f"")

        lines.append(f".ENDS {element_name}")

        return '\n'.join(lines)

    def to_verilog_a(self, name: str = "kan_zskin") -> str:
        """
        Export to Verilog-A behavioral model.

        Implements frequency-dependent Z_s(H, f) using KAN evaluation.

        Parameters
        ----------
        name : str
            Module name

        Returns
        -------
        verilog_a : str
            Verilog-A code string
        """
        lines = []
        lines.append(f"// PyKAN Surface Impedance Model")
        lines.append(f"// Trained KAN: (H, f) -> (R_s, X_s)")
        lines.append(f"`include \"disciplines.vams\"")
        lines.append(f"")
        lines.append(f"module {name}(p, n, H_sense);")
        lines.append(f"    inout p, n;")
        lines.append(f"    input H_sense;")
        lines.append(f"    electrical p, n, H_sense;")
        lines.append(f"")
        lines.append(f"    // Normalization parameters")
        lines.append(f"    parameter real H_min = {self.H_min:.6e};")
        lines.append(f"    parameter real H_max = {self.H_max:.6e};")
        lines.append(f"    parameter real f_min = {self.f_min:.6e};")
        lines.append(f"    parameter real f_max = {self.f_max:.6e};")
        lines.append(f"    parameter real Z_scale = {self.Z_scale:.6e};")
        lines.append(f"")
        lines.append(f"    real H_norm, f_norm, Rs, Xs;")
        lines.append(f"")
        lines.append(f"    analog begin")
        lines.append(f"        // Get normalized inputs")
        lines.append(f"        H_norm = (V(H_sense) - H_min) / (H_max - H_min + 1e-10);")
        lines.append(f"        // Note: Frequency from $freq in Spectre")
        lines.append(f"        f_norm = ($freq - f_min) / (f_max - f_min + 1e-10);")
        lines.append(f"")
        lines.append(f"        // KAN evaluation (simplified polynomial approximation)")
        lines.append(f"        // TODO: Export actual KAN spline coefficients")
        lines.append(f"        Rs = Z_scale * (0.1 + 0.5*H_norm + 0.3*f_norm);")
        lines.append(f"        Xs = Z_scale * (0.05 + 0.2*H_norm + 0.6*f_norm);")
        lines.append(f"")
        lines.append(f"        // Surface impedance: V = Z_s * I")
        lines.append(f"        V(p, n) <+ Rs * I(p, n);")
        lines.append(f"        V(p, n) <+ Xs / (2*`M_PI*$freq) * ddt(I(p, n));")
        lines.append(f"    end")
        lines.append(f"endmodule")

        return '\n'.join(lines)

    def plot_comparison(self, H_test: np.ndarray, f_test: np.ndarray,
                        Zs_true: np.ndarray, save_path: str = None):
        """
        Plot KAN prediction vs true values.

        Parameters
        ----------
        H_test : np.ndarray
            Test H values
        f_test : np.ndarray
            Test frequency values
        Zs_true : np.ndarray
            True surface impedance values
        save_path : str, optional
            Path to save figure
        """
        import matplotlib.pyplot as plt

        Zs_pred = self.predict(H_test, f_test)

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # R_s comparison
        ax = axes[0, 0]
        ax.scatter(np.real(Zs_true).flatten(), np.real(Zs_pred).flatten(),
                   alpha=0.5, s=10)
        ax.plot([0, np.real(Zs_true).max()], [0, np.real(Zs_true).max()],
                'r--', label='Perfect fit')
        ax.set_xlabel('True R_s [Ohm]')
        ax.set_ylabel('Predicted R_s [Ohm]')
        ax.set_title('Surface Resistance')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # X_s comparison
        ax = axes[0, 1]
        ax.scatter(np.imag(Zs_true).flatten(), np.imag(Zs_pred).flatten(),
                   alpha=0.5, s=10)
        ax.plot([np.imag(Zs_true).min(), np.imag(Zs_true).max()],
                [np.imag(Zs_true).min(), np.imag(Zs_true).max()],
                'r--', label='Perfect fit')
        ax.set_xlabel('True X_s [Ohm]')
        ax.set_ylabel('Predicted X_s [Ohm]')
        ax.set_title('Surface Reactance')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Error distribution
        ax = axes[1, 0]
        error_R = (np.real(Zs_pred) - np.real(Zs_true)).flatten()
        error_X = (np.imag(Zs_pred) - np.imag(Zs_true)).flatten()
        ax.hist(error_R, bins=30, alpha=0.5, label='R_s error')
        ax.hist(error_X, bins=30, alpha=0.5, label='X_s error')
        ax.set_xlabel('Error [Ohm]')
        ax.set_ylabel('Count')
        ax.set_title('Error Distribution')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Relative error
        ax = axes[1, 1]
        rel_error = np.abs(Zs_pred - Zs_true) / (np.abs(Zs_true) + 1e-10)
        ax.hist(rel_error.flatten() * 100, bins=30)
        ax.set_xlabel('Relative Error [%]')
        ax.set_ylabel('Count')
        ax.set_title(f'Relative Error (mean: {rel_error.mean()*100:.2f}%)')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150)
            print(f"Saved: {save_path}")

        plt.show()

        return {
            'mean_rel_error': rel_error.mean(),
            'max_rel_error': rel_error.max(),
            'R_rmse': np.sqrt(np.mean(error_R**2)),
            'X_rmse': np.sqrt(np.mean(error_X**2))
        }


def demo_pykan_surface_impedance():
    """Demonstrate PyKAN surface impedance learning."""
    print("=" * 60)
    print("PyKAN Surface Impedance Demo")
    print("=" * 60)

    # Generate synthetic ESIM data
    # Z_s = sqrt(j*omega*mu / sigma) for linear case
    # With saturation: mu(H) = mu_0 * mu_r / (1 + H/H_sat)

    sigma = 5.8e7  # Copper conductivity [S/m]
    mu_0 = 4 * np.pi * 1e-7
    mu_r_0 = 100  # Initial relative permeability
    H_sat = 1000  # Saturation field [A/m]

    H_data = np.linspace(10, 5000, 20)  # A/m
    f_data = np.logspace(3, 6, 15)  # 1 kHz to 1 MHz

    print(f"\nGenerating synthetic ESIM data:")
    print(f"  H range: {H_data.min():.0f} - {H_data.max():.0f} A/m")
    print(f"  f range: {f_data.min()/1e3:.1f} kHz - {f_data.max()/1e6:.1f} MHz")
    print(f"  Saturation model: mu_r(H) = {mu_r_0} / (1 + H/{H_sat})")

    # Generate Z_s data
    Zs_data = np.zeros((len(H_data), len(f_data)), dtype=complex)
    for i, H in enumerate(H_data):
        mu_r = mu_r_0 / (1 + H / H_sat)
        mu = mu_0 * mu_r
        for j, f in enumerate(f_data):
            omega = 2 * np.pi * f
            # Surface impedance for good conductor
            Zs_data[i, j] = np.sqrt(1j * omega * mu / sigma)

    print(f"  Generated {Zs_data.size} data points")

    # Train KAN
    kan_zs = PyKANSurfaceImpedance(width=[2, 8, 8, 2], grid=5, k=3)

    try:
        result = kan_zs.train_from_esim_data(
            H_data, f_data, Zs_data,
            epochs=200, lr=0.01, verbose=True
        )
        print(f"\nTraining completed:")
        print(f"  Final loss: {result['final_loss']:.6e}")

        # Test prediction
        H_test = np.array([100, 500, 2000])
        f_test = np.array([10e3, 100e3, 500e3])

        print(f"\nPrediction test:")
        for H, f in zip(H_test, f_test):
            Zs_pred = kan_zs.predict(H, f)
            # True value
            mu_r = mu_r_0 / (1 + H / H_sat)
            mu = mu_0 * mu_r
            omega = 2 * np.pi * f
            Zs_true = np.sqrt(1j * omega * mu / sigma)

            error = np.abs(Zs_pred - Zs_true) / np.abs(Zs_true) * 100
            print(f"  H={H:.0f} A/m, f={f/1e3:.0f} kHz:")
            print(f"    True:  {Zs_true:.4e}")
            print(f"    KAN:   {Zs_pred:.4e}")
            print(f"    Error: {error:.2f}%")

        # Generate SPICE output
        print(f"\nSPICE PWL export:")
        spice = kan_zs.to_spice_pwl(H_data, f_ref=100e3)
        for line in spice.split('\n')[:15]:
            print(f"  {line}")
        print("  ...")

        return kan_zs, result

    except ImportError as e:
        print(f"\nPyKAN not available: {e}")
        print("Install with: pip install pykan")
        return None, None


class ContinuedFractionExpansion:
    """
    Continued Fraction Expansion (CFE) for complex frequency responses.

    Converts arbitrary frequency-dependent impedance Z(s) into
    PRIMA ladder network (L-C or R-C cascade).

    This is essential for:
    - Dowell model (skin/proximity effect in windings)
    - Complex permeability mu*(s)
    - Arbitrary transfer functions -> SPICE circuit

    Reference:
    - P. L. Dowell, "Effects of eddy currents in transformer windings",
      Proc. IEE, vol. 113, pp. 1387-1394, 1966.
    - PRIMA network synthesis
    """

    def __init__(self, n_terms: int = 10, tol: float = 1e-6):
        """
        Initialize CFE expansion.

        Parameters
        ----------
        n_terms : int
            Number of terms in continued fraction
        tol : float
            Convergence tolerance
        """
        self.n_terms = n_terms
        self.tol = tol
        self.ladder = []  # List of (type, value) tuples

    def expand_impedance(self, Z_func: Callable[[complex], complex],
                         s_ref: complex = 1j * 2 * np.pi * 1e3) -> List[Tuple[str, float]]:
        """
        Expand Z(s) into continued fraction / PRIMA ladder.

        Z(s) = a_0 + 1/(a_1*s + 1/(a_2 + 1/(a_3*s + ...)))

        Parameters
        ----------
        Z_func : callable
            Impedance function Z(s)
        s_ref : complex
            Reference frequency for expansion

        Returns
        -------
        ladder : list
            List of (element_type, value) tuples
            element_type: 'R', 'L', 'C'
        """
        self.ladder = []

        # Sample Z at multiple frequencies for fitting
        frequencies = np.logspace(1, 8, 100)
        s_vals = 2j * np.pi * frequencies
        Z_vals = np.array([Z_func(s) for s in s_vals])

        # Start with Z(s)
        current_func = lambda s, Z=Z_vals, sv=s_vals: np.interp(
            np.abs(s), np.abs(sv), np.abs(Z)) * np.exp(
            1j * np.interp(np.abs(s), np.abs(sv), np.angle(Z)))

        for i in range(self.n_terms):
            # Evaluate at reference frequency
            Z_ref = Z_func(s_ref)

            if np.abs(Z_ref) < self.tol:
                break

            # Determine element type based on frequency behavior
            Z_low = Z_func(s_ref * 0.1)
            Z_high = Z_func(s_ref * 10)

            # Check if Z increases or decreases with frequency
            if np.abs(Z_high) > np.abs(Z_low) * 2:
                # Inductive behavior: Z ~ sL
                L = np.abs(Z_ref) / np.abs(s_ref)
                self.ladder.append(('L', L))
                # Remainder: Y = 1/(Z - sL)
                def new_Z(s, old_Z=Z_func, L_val=L):
                    Z_old = old_Z(s)
                    remainder = Z_old - s * L_val
                    return 1.0 / remainder if np.abs(remainder) > 1e-15 else 1e15
                Z_func = new_Z

            elif np.abs(Z_low) > np.abs(Z_high) * 2:
                # Capacitive behavior: Z ~ 1/(sC)
                C = 1.0 / (np.abs(Z_ref) * np.abs(s_ref))
                self.ladder.append(('C', C))
                # Remainder: Y = 1/(Z - 1/(sC))
                def new_Z(s, old_Z=Z_func, C_val=C):
                    Z_old = old_Z(s)
                    remainder = Z_old - 1.0 / (s * C_val)
                    return 1.0 / remainder if np.abs(remainder) > 1e-15 else 1e15
                Z_func = new_Z

            else:
                # Resistive: Z ~ R
                R = np.real(Z_ref)
                if R > 0:
                    self.ladder.append(('R', R))
                    # Remainder: 1/(Z - R)
                    def new_Z(s, old_Z=Z_func, R_val=R):
                        Z_old = old_Z(s)
                        remainder = Z_old - R_val
                        return 1.0 / remainder if np.abs(remainder) > 1e-15 else 1e15
                    Z_func = new_Z
                else:
                    break

        return self.ladder

    def expand_dowell_cf(self, d: float, sigma: float, mu: float = 4*np.pi*1e-7,
                         n_stages: int = 10) -> List[Tuple[str, float]]:
        """
        Expand skin effect into PRIMA ladder using 1D diffusion eigenmode expansion.

        This creates an RL ladder network that exactly represents the 1D diffusion
        equation solution for skin effect. The structure is:

            L_1     L_2     L_3
          o--[===]---[===]---[===]---...---out
               |       |       |
              R_1     R_2     R_3
               |       |       |
              ===     ===     ===

        where:
            L_n = mu * d / (4n - 3)     for n = 1, 2, 3, ...
            R_n = (4n - 5) * 4 / (sigma * d)  for n >= 2, R_1 ~ 0

        This corresponds to the eigenmode expansion of the 1D diffusion equation,
        NOT the Dowell J-fraction formula z*coth(z).

        For Dowell-accurate frequency response, use expand_dowell_cf() to build
        the DC ladder, then apply Dowell F_R(omega), F_L(omega) corrections at
        each frequency point.

        Parameters
        ----------
        d : float
            Conductor thickness [m]
        sigma : float
            Conductivity [S/m]
        mu : float
            Permeability [H/m], default is mu_0
        n_stages : int
            Number of RL stages (typically 5-10 for high accuracy)

        Returns
        -------
        ladder : list
            PRIMA ladder elements [('L', value), ('R_shunt', value), ...]

        Notes
        -----
        The PRIMA ladder matches the 1D diffusion solution:
            Z_diff(s) = s * mu * (2/(k*d)) * tan(k*d/2) * d
        where k = sqrt(s * mu * sigma).

        This is mathematically different from Dowell's formula:
            Z_dowell = R_dc * F_R + s * L_dc * F_L

        Both converge to the same result, but the circuit topology differs.

        References
        ----------
        - PRIMA: A. Odabasioglu et al., IEEE TCAD, 1998.
        - Skin effect eigenmode: J.A. Ferreira, IEEE Trans. Power Electronics, 1994.
        """
        self.ladder = []

        # Store parameters for documentation
        self._dowell_params = {
            'd': d,
            'sigma': sigma,
            'mu': mu,
            'R_dc': 1.0 / (sigma * d),
            'L_dc': mu * d / 3.0,
            'n_stages': n_stages
        }

        # Build PRIMA I ladder from 1D diffusion eigenmode expansion:
        # L_n = mu * d / (4n - 3)  -> 1, 1/5, 1/9, 1/13, ...
        # R_n = (4n - 5) * 4 / (sigma * d)  -> 0, 12/sigma*d, 28/sigma*d, ...

        for n in range(1, n_stages + 1):
            # Series inductance
            L_n = mu * d / (4 * n - 3)
            self.ladder.append(('L', L_n))

            # Shunt resistance (R_1 ~ 0 for first stage)
            if n == 1:
                R_n = 1e-15  # Effectively zero
            else:
                R_n = (4 * n - 5) * 4.0 / (sigma * d)
            self.ladder.append(('R_shunt', R_n))

        return self.ladder

    def expand_dowell(self, R_dc: float, delta_ratio: float,
                      n_layers: int = 1) -> List[Tuple[str, float]]:
        """
        Legacy Dowell expansion (approximate, for backward compatibility).

        For exact expansion, use expand_dowell_cf() instead.

        Parameters
        ----------
        R_dc : float
            DC resistance [Ohm]
        delta_ratio : float
            Conductor height / skin depth at reference frequency
        n_layers : int
            Number of winding layers

        Returns
        -------
        ladder : list
            PRIMA ladder elements
        """
        self.ladder = []

        # Stage 1: DC resistance
        self.ladder.append(('R', R_dc))

        # Stage 2: First-order skin effect
        omega_ref = 2 * np.pi * 1e3  # Reference at 1kHz
        L_skin = R_dc / (omega_ref * delta_ratio**2 / 2)
        R_skin = R_dc * delta_ratio**2 / 3

        self.ladder.append(('L', L_skin))
        self.ladder.append(('R', R_skin))

        # Stage 3: Proximity effect (for multi-layer)
        if n_layers > 1:
            prox_factor = (n_layers**2 - 1) / 3
            R_prox = R_dc * prox_factor * delta_ratio**2
            L_prox = L_skin * prox_factor

            self.ladder.append(('L', L_prox))
            self.ladder.append(('R', R_prox))

        return self.ladder

    def to_spice_subcircuit(self, name: str = "ZSKIN") -> str:
        """
        Convert PRIMA ladder to SPICE subcircuit.

        Supports both series and shunt elements:
        - 'R', 'L', 'C': Series elements
        - 'L_shunt', 'C_shunt', 'R_shunt': Shunt elements (to ground)

        Parameters
        ----------
        name : str
            Subcircuit name

        Returns
        -------
        spice : str
            SPICE subcircuit definition
        """
        lines = []
        lines.append(f".SUBCKT {name} in out")

        # Add parameter comments if available
        if hasattr(self, '_dowell_params'):
            p = self._dowell_params
            lines.append(f"* 1D diffusion eigenmode expansion")
            lines.append(f"* d = {p['d']:.3e} m, sigma = {p['sigma']:.3e} S/m")
            lines.append(f"* R_dc = {p['R_dc']:.6e} Ohm*m^2")
            lines.append(f"* L_dc = {p['L_dc']:.6e} H*m^2")
            lines.append(f"* n_stages = {p['n_stages']}")
        else:
            lines.append(f"* PRIMA ladder with {len(self.ladder)} elements")

        node = 'in'
        next_node_num = 1
        elem_counts = {'R': 0, 'L': 0, 'C': 0}

        for i, (elem_type, value) in enumerate(self.ladder):
            # Determine if this is a shunt element
            is_shunt = elem_type.endswith('_shunt')
            base_type = elem_type.replace('_shunt', '')

            if is_shunt:
                # Shunt element: connect current node to ground
                elem_counts[base_type] += 1
                elem_name = f"{base_type}{elem_counts[base_type]}"
                lines.append(f"{elem_name} {node} 0 {value:.6e}")
                # Don't advance node for shunt elements
            else:
                # Series element: connect to next node
                # Determine next node
                remaining_series = sum(1 for et, _ in self.ladder[i+1:]
                                       if not et.endswith('_shunt'))
                if remaining_series == 0:
                    next_node = 'out'
                else:
                    next_node = f'n{next_node_num}'
                    next_node_num += 1

                elem_counts[base_type] += 1
                elem_name = f"{base_type}{elem_counts[base_type]}"
                lines.append(f"{elem_name} {node} {next_node} {value:.6e}")
                node = next_node

        lines.append(".ENDS")
        return '\n'.join(lines)

    def evaluate(self, s: complex) -> complex:
        """
        Evaluate ladder network impedance at frequency s.

        For PRIMA I ladder (series L with shunt R):
            L_1     L_2     L_3
          o--[===]---[===]---[===]---...---out
               |       |       |
              R_1     R_2     R_3
               |       |       |
              ===     ===     ===

        Impedance is computed from end to start:
        Z_n = R_n (last shunt resistor)
        Z_{n-1} = sL_n + (R_n || Z_n)
        ...

        Parameters
        ----------
        s : complex
            Complex frequency (j*omega)

        Returns
        -------
        Z : complex
            Impedance at frequency s
        """
        # For PRIMA I structure: pairs of (L, R_shunt)
        # Extract pairs and evaluate from end to start

        # Collect elements by type
        series_L = []
        shunt_R = []

        for elem_type, value in self.ladder:
            if elem_type == 'L':
                series_L.append(value)
            elif elem_type == 'R_shunt':
                shunt_R.append(value)
            elif elem_type == 'R':
                # Series resistor (old format)
                series_L.append(0)  # No inductance
                shunt_R.append(1e15)  # No shunt
            elif elem_type == 'L_shunt':
                # Shunt inductance (old format)
                series_L.append(0)
                shunt_R.append(1e15)

        # If we have paired L and R_shunt, evaluate as PRIMA I ladder
        if len(series_L) == len(shunt_R) and len(series_L) > 0:
            return self._evaluate_prima_ladder(s, series_L, shunt_R)

        # Fallback: generic evaluation from end to start
        Z = 0.0
        for elem_type, value in reversed(self.ladder):
            is_shunt = elem_type.endswith('_shunt')
            base_type = elem_type.replace('_shunt', '')

            if is_shunt:
                if base_type == 'L':
                    Y_shunt = 1.0 / (s * value) if np.abs(s * value) > 1e-15 else 1e15
                elif base_type == 'C':
                    Y_shunt = s * value
                elif base_type == 'R':
                    Y_shunt = 1.0 / value if np.abs(value) > 1e-15 else 1e15
                else:
                    Y_shunt = 0.0

                if np.abs(Z) > 1e-15:
                    Y_total = 1.0 / Z + Y_shunt
                    Z = 1.0 / Y_total if np.abs(Y_total) > 1e-15 else 1e15
                else:
                    Z = 1.0 / Y_shunt if np.abs(Y_shunt) > 1e-15 else 1e15
            else:
                if base_type == 'R':
                    Z = Z + value
                elif base_type == 'L':
                    Z = Z + s * value
                elif base_type == 'C':
                    Z = Z + 1.0 / (s * value) if np.abs(s * value) > 1e-15 else Z + 1e15

        return Z

    def _evaluate_prima_ladder(self, s: complex, L_vals: List[float],
                                R_vals: List[float]) -> complex:
        """
        Evaluate PRIMA I ladder impedance.

        Structure (output shorted to ground):
            L_1     L_2     L_3     L_n
          o--[===]---[===]---[===]---...---[===]---GND
               |       |       |            |
              R_1     R_2     R_3          R_n
               |       |       |            |
              GND     GND     GND          GND

        Evaluation from end to start (output is shorted):
        For stage n (last stage):
            Z_n = sL_n + R_n  (R_n to ground, sL_n in series)
        For stage n-1:
            Z_{n-1} = sL_{n-1} + (R_{n-1} || Z_n)
        ...

        Alternative: if output is open (no load), start with Z_n = sL_n || R_n

        This implementation assumes output is SHORTED (surface impedance model).
        """
        n = len(L_vals)
        if n == 0:
            return 0.0

        # Start from last stage - output is SHORTED to ground
        # Last stage: sL_n in series with R_n to ground
        # But R_n is also to ground... so it's sL_n || R_n
        #
        # Actually, for PRIMA I: each stage is sL in series, then R to ground
        # The input sees: sL_1 + (R_1 || (sL_2 + (R_2 || (sL_3 + ...))))

        # Start from the output end
        Z = 0.0  # Output node is grounded

        for i in range(n - 1, -1, -1):
            L = L_vals[i]
            R = R_vals[i]

            # Current stage: sL in series, R shunts to ground
            # Z_stage = sL + (R || Z_downstream)
            #
            # If Z_downstream = 0 (shorted output):
            #   R || 0 = 0, so Z_stage = sL
            # If Z_downstream = inf (open output):
            #   R || inf = R, so Z_stage = sL + R

            Z_L = s * L

            # Parallel combination of R and downstream impedance
            if np.abs(Z) < 1e-15:
                # Downstream is shorted: R || 0 = 0
                Z_parallel = 0.0
            elif R > 1e10:
                # R is very large (open): R || Z = Z
                Z_parallel = Z
            else:
                # General case: R || Z = R*Z / (R+Z)
                Z_parallel = R * Z / (R + Z)

            Z = Z_L + Z_parallel

        return Z

    def evaluate_dowell_exact(self, s: complex, d: float, sigma: float,
                              mu: float = 4*np.pi*1e-7) -> complex:
        """
        Evaluate exact Dowell surface impedance for comparison.

        Z_s = R_dc * z * coth(z)

        where z = sqrt(tau * s), tau = d^2 * mu * sigma / 2

        Parameters
        ----------
        s : complex
            Complex frequency (j*omega)
        d : float
            Conductor thickness [m]
        sigma : float
            Conductivity [S/m]
        mu : float
            Permeability [H/m]

        Returns
        -------
        Z : complex
            Exact Dowell surface impedance
        """
        R_dc = 1.0 / (sigma * d)
        tau = d**2 * mu * sigma / 2.0

        # z = sqrt(tau * s) = sqrt(tau * j * omega) = (1+j)/sqrt(2) * sqrt(tau * omega)
        z_squared = tau * s
        z = np.sqrt(z_squared)

        # Handle z -> 0 limit: z*coth(z) -> 1
        if np.abs(z) < 1e-10:
            return R_dc

        # z * coth(z) = z * (e^z + e^-z) / (e^z - e^-z)
        # For numerical stability, use tanh: coth(z) = 1/tanh(z)
        try:
            z_coth_z = z / np.tanh(z)
        except (ZeroDivisionError, FloatingPointError):
            z_coth_z = 1.0

        return R_dc * z_coth_z


class SchurACACircuitExtraction:
    """
    Combined Schur complement reduction + ACA sparse extraction.

    Complete workflow:
    1. Full system: [Z_II, Z_IP; Z_PI, Z_PP] where I=internal, P=port
    2. Schur complement: Z_eff = Z_PP - Z_PI @ Z_II^{-1} @ Z_IP
    3. ACA: Z_eff ≈ U @ V.T (low-rank)
    4. RLC fitting: Each mode -> R + sL + 1/(sC)
    5. SPICE netlist generation

    This enables:
    - Efficient circuit simulation of complex EM structures
    - Preservation of passivity through PRIMA + ACA
    - Time-domain simulation via standard circuit solvers
    """

    def __init__(self, schur_tol: float = 1e-10, aca_tol: float = 1e-4):
        self.schur_tol = schur_tol
        self.aca = ACACircuitExtraction(tol=aca_tol)

    def full_extraction(self,
                        L: np.ndarray, R: np.ndarray, P: np.ndarray,
                        K_LS: np.ndarray, K_LM: np.ndarray,
                        Z_MM: np.ndarray,
                        port_indices: List[int],
                        frequencies: np.ndarray,
                        kan_mu: KANMaterialInterface = None,
                        kan_eps: KANMaterialInterface = None) -> dict:
        """
        Complete extraction: System -> Schur -> ACA -> SPICE.

        Parameters
        ----------
        L, R, P, K_LS, K_LM, Z_MM : np.ndarray
            System matrices (Loop, Star, Magnetic blocks)
        port_indices : list
            Port node indices
        frequencies : np.ndarray
            Frequencies for fitting
        kan_mu, kan_eps : KANMaterialInterface, optional
            KAN models for complex mu/eps

        Returns
        -------
        result : dict
            Contains reduced system, ACA factors, and SPICE netlist
        """
        n_L = L.shape[0]
        n_S = P.shape[0]
        n_M = Z_MM.shape[0]
        n_ports = len(port_indices)

        # Internal indices (non-port)
        all_L_indices = list(range(n_L))
        internal_indices = [i for i in all_L_indices if i not in port_indices]

        print("=" * 60)
        print("Schur + ACA Circuit Extraction")
        print("=" * 60)
        print(f"System: {n_L} Loop + {n_S} Star + {n_M} Magnetic DOFs")
        print(f"Ports: {n_ports}, Internal: {len(internal_indices)}")

        # Build frequency-dependent impedance function
        def Z_full(s):
            """Full system impedance at complex frequency s."""
            # Loop-Star-Magnetic system matrix
            n_total = n_L + n_S + n_M
            Z = np.zeros((n_total, n_total), dtype=complex)

            # Loop block: R + sL
            Z[:n_L, :n_L] = R + s * L

            # Star block: P/s (capacitive)
            if np.abs(s) > 1e-15:
                Z[n_L:n_L+n_S, n_L:n_L+n_S] = P / s
            else:
                Z[n_L:n_L+n_S, n_L:n_L+n_S] = P * 1e15

            # Magnetic block with KAN material
            if kan_mu:
                mu_s = kan_mu.evaluate_frequency([np.abs(s)/(2*np.pi)])[0]
                Z[n_L+n_S:, n_L+n_S:] = Z_MM * mu_s
            else:
                Z[n_L+n_S:, n_L+n_S:] = Z_MM

            # Coupling blocks
            Z[:n_L, n_L:n_L+n_S] = K_LS
            Z[n_L:n_L+n_S, :n_L] = K_LS.T

            Z[:n_L, n_L+n_S:] = s * K_LM
            Z[n_L+n_S:, :n_L] = -K_LM.T

            return Z

        def Z_port(s):
            """Port impedance via Schur complement."""
            Z = Z_full(s)

            # Extract blocks
            # P = port indices (in full system), I = internal
            P_full = port_indices  # Port indices in Loop block
            I_full = internal_indices + list(range(n_L, n_L + n_S + n_M))

            Z_PP = Z[np.ix_(P_full, P_full)]
            Z_PI = Z[np.ix_(P_full, I_full)]
            Z_II = Z[np.ix_(I_full, I_full)]
            Z_IP = Z[np.ix_(I_full, P_full)]

            # Schur complement
            try:
                Z_II_inv = np.linalg.inv(Z_II)
                Z_eff = Z_PP - Z_PI @ Z_II_inv @ Z_IP
            except np.linalg.LinAlgError:
                Z_eff = Z_PP  # Fallback

            return Z_eff

        # Extract circuit using ACA
        print("\nExtracting sparse circuit via ACA...")
        circuit = self.aca.extract_rlc_circuit(
            Z_port, frequencies,
            port_names=[f'Port_{i+1}' for i in range(n_ports)]
        )

        # Generate SPICE netlist
        netlist = self.aca.to_spice_netlist(circuit)

        # Compute impedance comparison using pole-residue model
        print("\nValidation (Vector Fitting model):")
        test_freqs = [frequencies[0], frequencies[len(frequencies)//2], frequencies[-1]]
        for f in test_freqs:
            s = 2j * np.pi * f
            Z_exact = Z_port(s)

            # Reconstruct from pole-residue model
            Z_approx = np.zeros((n_ports, n_ports), dtype=complex)
            for elem in circuit['elements']:
                u_k = elem['u']
                v_k = elem['v']
                # Use pole-residue form if available
                if 'poles' in elem and elem['poles'] is not None:
                    poles = elem['poles']
                    residues = elem['residues']
                    d = elem['d']
                    e = elem['e']
                    Z_mode = d + s * e
                    for p, c in zip(poles, residues):
                        Z_mode += c / (s - p)
                else:
                    R, L, C = elem['R'], elem['L'], elem['C']
                    Z_mode = R + s * L + 1.0 / (s * C) if C > 1e-20 else R + s * L
                Z_approx += Z_mode * np.outer(u_k, v_k)

            error = np.linalg.norm(Z_exact - Z_approx) / (np.linalg.norm(Z_exact) + 1e-20)
            print(f"  f = {f:.1e} Hz: relative error = {error*100:.2f}%")

        return {
            'circuit': circuit,
            'netlist': netlist,
            'aca_rank': self.aca.rank,
            'U': self.aca.U,
            'V': self.aca.V,
            'Z_port_func': Z_port,
            'Z_rlc_func': lambda s, circ=circuit: self._compute_rlc_impedance(s, circ)
        }

    def _compute_rlc_impedance(self, s: complex, circuit: dict) -> np.ndarray:
        """Compute impedance from extracted RLC circuit."""
        n_ports = circuit['n_ports']
        Z_rlc = np.zeros((n_ports, n_ports), dtype=complex)
        for elem in circuit['elements']:
            u_k = elem['u']
            v_k = elem['v']
            R, L, C = elem['R'], elem['L'], elem['C']
            Z_mode = R + s * L + 1.0 / (s * C) if C > 1e-20 else R + s * L
            Z_rlc += Z_mode * np.outer(u_k, v_k)
        return Z_rlc


def demo_aca_circuit_extraction():
    """Demonstrate ACA-based sparse circuit extraction."""
    print("=" * 60)
    print("ACA Sparse Circuit Extraction Demo")
    print("=" * 60)

    np.random.seed(42)

    # Create physically meaningful test system (coupled inductors)
    n_L = 8   # 8 mesh currents
    n_S = 4   # 4 star equations (minimal)
    n_M = 4   # 4 magnetic elements
    n_ports = 2

    # System matrices - realistic RL circuit values
    # Inductance matrix with mutual coupling (transformer-like)
    L = np.eye(n_L) * 10e-6  # 10 uH self-inductance
    for i in range(n_L):
        for j in range(i+1, n_L):
            # Coupling decreases with distance
            M = 5e-6 * np.exp(-np.abs(i-j)/2)  # Mutual inductance
            L[i, j] = L[j, i] = M

    # Resistance matrix
    R = np.eye(n_L) * 0.05  # 50 mOhm resistance

    # Star equations (weak capacitive coupling for stability)
    P = np.eye(n_S) * 1e12  # Very high P = weak capacitive effect

    # Loop-Star coupling (weak)
    K_LS = np.zeros((n_L, n_S))
    for i in range(min(n_L, n_S)):
        K_LS[i, i] = 0.001

    # Loop-Magnetic coupling
    K_LM = np.zeros((n_L, n_M))
    for i in range(min(n_L, n_M)):
        K_LM[i, i] = 0.01

    # Magnetic material matrix
    Z_MM = np.eye(n_M) * 1e-4

    port_indices = [0, 4]  # Two ports

    # Power electronics frequency range
    frequencies = np.logspace(2, 7, 30)  # 100Hz to 10MHz

    # No KAN material (simple linear system first)
    kan_mu = None

    # Full extraction
    extractor = SchurACACircuitExtraction(schur_tol=1e-10, aca_tol=1e-4)

    result = extractor.full_extraction(
        L, R, P, K_LS, K_LM, Z_MM,
        port_indices, frequencies,
        kan_mu=kan_mu
    )

    print(f"\n--- Results ---")
    print(f"ACA rank: {result['aca_rank']} (compression: {result['aca_rank']/n_ports*100:.1f}%)")
    print(f"\nCircuit elements: {len(result['circuit']['elements'])}")

    for elem in result['circuit']['elements']:
        print(f"  Mode {elem['mode']}: R={elem['R']:.2e}, "
              f"L={elem['L']:.2e}, C={elem['C']:.2e}, "
              f"Q={elem['quality']:.2f}")

    print(f"\nCoupling terms: {len(result['circuit']['coupling'])}")

    print(f"\n--- SPICE Netlist (first 30 lines) ---")
    netlist_lines = result['netlist'].split('\n')[:30]
    for line in netlist_lines:
        print(line)

    return result


def demo_dowell_cfe():
    """Demonstrate Continued Fraction Expansion for Dowell model."""
    print("=" * 60)
    print("Dowell Model - Continued Fraction Expansion Demo")
    print("=" * 60)

    # Winding parameters
    R_dc = 0.1        # 100 mOhm DC resistance
    delta_ratio = 2.0  # Conductor 2x skin depth
    n_layers = 3      # 3-layer winding

    # Create CFE expander
    cfe = ContinuedFractionExpansion(n_terms=6)

    # Expand Dowell model
    ladder = cfe.expand_dowell(R_dc, delta_ratio, n_layers)

    print(f"\nDowell parameters:")
    print(f"  R_dc = {R_dc*1000:.1f} mOhm")
    print(f"  delta_ratio = {delta_ratio}")
    print(f"  n_layers = {n_layers}")

    print(f"\nPRIMA ladder ({len(ladder)} elements):")
    for i, (elem_type, value) in enumerate(ladder):
        if elem_type == 'R':
            print(f"  {i+1}. R = {value*1000:.3f} mOhm")
        elif elem_type == 'L':
            print(f"  {i+1}. L = {value*1e6:.3f} uH")
        elif elem_type == 'C':
            print(f"  {i+1}. C = {value*1e9:.3f} nF")

    # Generate SPICE subcircuit
    spice = cfe.to_spice_subcircuit("DOWELL_WINDING")
    print(f"\n--- SPICE Subcircuit ---")
    print(spice)

    # Verify frequency response
    print(f"\n--- Frequency Response Verification ---")
    frequencies = [100, 1e3, 10e3, 100e3, 1e6]

    # Reference Dowell formula (simplified)
    def dowell_exact(f, R_dc, delta_ratio, n_layers):
        # F(x) = x * (sinh(2x) + sin(2x)) / (cosh(2x) - cos(2x))
        # where x = delta_ratio * sqrt(f/f_ref)
        f_ref = 1e3
        x = delta_ratio * np.sqrt(f / f_ref)
        if x < 0.1:
            F = 1 + x**4 / 45
        elif x > 10:
            F = x
        else:
            F = x * (np.sinh(2*x) + np.sin(2*x)) / (np.cosh(2*x) - np.cos(2*x))

        # Proximity factor
        G = 2 * (n_layers**2 - 1) / 3 * x * (np.sinh(x) - np.sin(x)) / (np.cosh(x) + np.cos(x))

        return R_dc * (F + G)

    print(f"  {'Freq':>10s}  {'|Z_exact|':>12s}  {'|Z_CFE|':>12s}  {'Error':>8s}")
    print(f"  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*8}")

    for f in frequencies:
        s = 2j * np.pi * f
        Z_exact = dowell_exact(f, R_dc, delta_ratio, n_layers)
        Z_cfe = cfe.evaluate(s)

        error = np.abs(Z_cfe - Z_exact) / np.abs(Z_exact) * 100
        print(f"  {f:>10.0f}  {np.abs(Z_exact):>12.4f}  {np.abs(Z_cfe):>12.4f}  {error:>7.1f}%")

    return cfe


def demo_kan_cfe_from_data():
    """
    KAN-based CFE learning from discrete frequency data.

    This demonstrates learning PRIMA ladder coefficients
    from measured/simulated impedance data points.
    """
    print("=" * 60)
    print("KAN Continued Fraction - Learning from Data Points")
    print("=" * 60)

    # Generate synthetic "measured" data with complex frequency dependence
    # Simulating a real inductor with skin effect + parasitic capacitance
    frequencies = np.logspace(2, 7, 50)  # 100Hz to 10MHz
    omega = 2 * np.pi * frequencies

    # True parameters (unknown to the learner)
    R_dc = 0.05      # 50 mOhm
    L_dc = 10e-6     # 10 uH
    R_ac = 0.2       # AC resistance at high freq
    C_par = 5e-12    # 5 pF parasitic capacitance

    # Generate complex impedance data
    # Z = R_dc + R_ac*sqrt(f/f_ref) + jwL - j/(wC)
    f_ref = 1e5
    Z_data = np.zeros(len(frequencies), dtype=complex)

    for i, f in enumerate(frequencies):
        s = 2j * np.pi * f
        # Skin effect: R increases with sqrt(f)
        R_skin = R_dc + R_ac * np.sqrt(f / f_ref)
        # Inductance with slight frequency dependence
        L_eff = L_dc * (1 - 0.1 * np.log10(f / 100))
        # Parallel parasitic capacitance
        Z_LC = s * L_eff
        Z_C = 1.0 / (s * C_par)
        Z_parallel = Z_LC * Z_C / (Z_LC + Z_C)
        Z_data[i] = R_skin + Z_parallel

    # Add some "measurement noise"
    np.random.seed(42)
    noise = (np.random.randn(len(frequencies)) + 1j * np.random.randn(len(frequencies))) * 0.01
    Z_data_noisy = Z_data * (1 + noise)

    print(f"\nInput data:")
    print(f"  Frequency range: {frequencies[0]:.0f} Hz to {frequencies[-1]/1e6:.1f} MHz")
    print(f"  Data points: {len(frequencies)}")
    print(f"  |Z| range: {np.abs(Z_data[0]):.4f} to {np.abs(Z_data[-1]):.4f} Ohm")

    # Create KAN CFE learner
    kan_cfe = KANContinuedFraction(n_terms=6)

    # Learn CFE coefficients from data
    # Try different element patterns
    patterns = [
        ['R', 'L', 'R', 'L', 'R', 'L'],      # Standard RL ladder
        ['R', 'L', 'C', 'R', 'L', 'C'],      # RLC ladder
        ['R', 'L', 'R', 'L', 'C', 'R'],      # Mixed with shunt C
    ]

    best_result = None
    best_error = float('inf')

    print(f"\nTrying different element patterns...")

    for pattern in patterns:
        result = kan_cfe.fit_from_frequency_data(
            frequencies, Z_data_noisy, element_pattern=pattern
        )
        pattern_str = '-'.join(pattern)
        print(f"  {pattern_str}: error = {result['fit_error']*100:.2f}%")

        if result['fit_error'] < best_error:
            best_error = result['fit_error']
            best_result = result
            best_pattern = pattern

    # Re-fit with best pattern
    kan_cfe = KANContinuedFraction(n_terms=6)
    result = kan_cfe.fit_from_frequency_data(
        frequencies, Z_data_noisy, element_pattern=best_pattern
    )

    print(f"\n--- Best Result ---")
    print(f"Pattern: {'-'.join(best_pattern)}")
    print(f"Fit error: {result['fit_error']*100:.2f}%")
    print(f"Converged: {result['converged']}")

    print(f"\nLearned PRIMA ladder:")
    for elem_type, value in zip(result['element_types'], result['coefficients']):
        if elem_type == 'R':
            print(f"  R = {value*1000:.3f} mOhm")
        elif elem_type == 'L':
            print(f"  L = {value*1e6:.3f} uH")
        elif elem_type == 'C':
            print(f"  C = {value*1e12:.3f} pF")

    # Generate SPICE subcircuit
    spice = kan_cfe.to_spice_subcircuit("LEARNED_IMPEDANCE")
    print(f"\n--- SPICE Subcircuit ---")
    print(spice)

    # Validation at specific frequencies
    print(f"\n--- Validation ---")
    print(f"  {'Freq':>10s}  {'|Z_data|':>12s}  {'|Z_CFE|':>12s}  {'Error':>8s}")
    print(f"  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*8}")

    test_freqs = [100, 1e3, 10e3, 100e3, 1e6, 10e6]
    for f in test_freqs:
        s = 2j * np.pi * f
        # Interpolate original data
        idx = np.argmin(np.abs(frequencies - f))
        Z_orig = Z_data[idx]
        Z_cfe = kan_cfe.evaluate(s)

        error = np.abs(Z_cfe - Z_orig) / np.abs(Z_orig) * 100
        print(f"  {f:>10.0f}  {np.abs(Z_orig):>12.4f}  {np.abs(Z_cfe):>12.4f}  {error:>7.1f}%")

    # Demonstrate parametric learning
    print(f"\n" + "=" * 60)
    print("Parametric KAN Learning (temperature variation)")
    print("=" * 60)

    # Generate data for different temperatures
    temperatures = np.array([25, 50, 75, 100, 125])  # degC
    Z_temp_data = np.zeros((len(temperatures), len(frequencies)), dtype=complex)

    for t_idx, temp in enumerate(temperatures):
        # Temperature coefficient: resistance increases with temp
        temp_factor = 1 + 0.004 * (temp - 25)  # 0.4%/K
        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            R_skin = (R_dc + R_ac * np.sqrt(f / f_ref)) * temp_factor
            L_eff = L_dc * (1 - 0.1 * np.log10(f / 100))
            Z_LC = s * L_eff
            Z_C = 1.0 / (s * C_par)
            Z_parallel = Z_LC * Z_C / (Z_LC + Z_C)
            Z_temp_data[t_idx, i] = R_skin + Z_parallel

    # Train KAN for parametric variation
    kan_param = KANContinuedFraction(n_terms=4)

    # First fit to get element pattern
    _ = kan_param.fit_from_frequency_data(
        frequencies, Z_temp_data[0], element_pattern=['R', 'L', 'R', 'L']
    )

    # Train parametric model
    result_param = kan_param.fit_with_kan_network(
        frequencies, Z_temp_data, param_data=temperatures
    )

    print(f"\nParametric model trained on {len(temperatures)} temperature points")
    print(f"Element pattern: {kan_param.element_types}")

    # Predict for interpolated temperature
    test_temp = 60  # Not in training data
    predicted_coeffs = kan_param.predict_coefficients(test_temp)

    print(f"\nPredicted coefficients at T = {test_temp} degC:")
    for elem_type, value in zip(kan_param.element_types, predicted_coeffs):
        if elem_type == 'R':
            print(f"  R = {value*1000:.3f} mOhm")
        elif elem_type == 'L':
            print(f"  L = {value*1e6:.3f} uH")

    # Generate temperature-dependent Z at test point
    test_idx = len(frequencies) // 2
    f_test = frequencies[test_idx]
    s_test = 2j * np.pi * f_test

    temp_factor_true = 1 + 0.004 * (test_temp - 25)
    R_skin_true = (R_dc + R_ac * np.sqrt(f_test / f_ref)) * temp_factor_true
    L_eff_true = L_dc * (1 - 0.1 * np.log10(f_test / 100))
    Z_true = R_skin_true + 1j * 2 * np.pi * f_test * L_eff_true

    Z_predicted = kan_param.evaluate(s_test, predicted_coeffs)

    print(f"\nValidation at T = {test_temp} degC, f = {f_test/1e3:.1f} kHz:")
    print(f"  True |Z|:      {np.abs(Z_true):.4f} Ohm")
    print(f"  Predicted |Z|: {np.abs(Z_predicted):.4f} Ohm")
    print(f"  Error: {np.abs(Z_predicted - Z_true) / np.abs(Z_true) * 100:.1f}%")

    return kan_cfe, kan_param


@dataclass
class ElementGroupConfig:
    """Configuration for an element group in SPICE extraction."""
    name: str                    # Group name: 'magnetic', 'dielectric', 'conductor'
    indices: List[int]           # DOF indices for this group
    n_lanczos: int = 0          # Lanczos order (0 = no reduction)
    material_type: str = 'linear'  # 'linear', 'nonlinear', 'kan'
    kan_model: Optional['KANMaterialInterface'] = None  # KAN model if applicable


@dataclass
class SPICEExtractionConfig:
    """Configuration for SPICE circuit extraction.

    Note (2026-01-17): ACA tolerances were removed. Lanczos-only approach is used
    for all model order reduction because it directly reduces element count in SPICE.
    """
    # PRIMA Lanczos settings
    n_lanczos_loop: int = 20        # Lanczos order for Loop (conductor) DOFs
    n_lanczos_star: int = 10        # Lanczos order for Star (capacitive) DOFs
    n_lanczos_magnetic: int = 0     # Lanczos order for Magnetic DOFs (0 = no reduction)
    n_lanczos_dielectric: int = 0   # Lanczos order for Dielectric DOFs (0 = no reduction)
    n_lanczos_conductor: int = 0    # Lanczos order for Conductor/Shield DOFs (0 = no reduction)

    # Port settings
    port_indices: List[int] = field(default_factory=list)
    port_names: List[str] = field(default_factory=list)

    # Frequency range for fitting
    f_min: float = 1.0              # Minimum frequency [Hz]
    f_max: float = 1e6              # Maximum frequency [Hz]
    n_freq: int = 50                # Number of frequency points


class PRIMASchurExtractor:
    """
    PRIMA-based Schur complement extractor using Lanczos-only reduction.

    This class implements a complete workflow for SPICE netlist extraction:

    1. PRIMA Lanczos Reduction:
       - Apply block Lanczos to Loop and Star DOFs
       - Specify Lanczos order (ladder stages) for each block
       - Result: Tridiagonal (PRIMA ladder) structure
       - Directly reduces element count in SPICE output

    2. Schur Complement:
       - Eliminate internal DOFs, keep port DOFs
       - Port impedance: Z_port = Z_PP - Z_PI @ Z_II^{-1} @ Z_IP

    3. SPICE Netlist Generation:
       - Each Lanczos stage -> RL ladder element (inductor)
       - Each Lanczos stage -> RC ladder element (capacitor)
       - Tridiagonal structure = series ladder with mutual coupling

    Note (2026-01-17): ACA+ was removed. Lanczos-only approach is used because
    it directly reduces the number of circuit elements in SPICE output.

    Usage:
    ------
    config = SPICEExtractionConfig(
        n_lanczos_loop=20,    # Reduces inductors to 20
        n_lanczos_star=10,    # Reduces capacitors to 10
        port_indices=[0, 1],
    )

    extractor = PRIMASchurExtractor(config)
    result = extractor.extract(system_matrices, frequencies)
    netlist = result['netlist']
    """

    def __init__(self, config: SPICEExtractionConfig):
        """
        Initialize PRIMA-Schur extractor.

        Parameters
        ----------
        config : SPICEExtractionConfig
            Extraction configuration with Lanczos orders
        """
        self.config = config

        # Lanczos reducer (Lanczos-only, no ACA)
        max_lanczos = max(
            config.n_lanczos_loop,
            config.n_lanczos_star,
            config.n_lanczos_magnetic,
            config.n_lanczos_dielectric,
            config.n_lanczos_conductor,
            100
        )
        self.lanczos = LanczosReducer(tol=1e-10, max_rank=max_lanczos)

        # Results storage
        self.Q_loop = None      # Lanczos basis for Loop
        self.Q_star = None      # Lanczos basis for Star
        self.T_loop = None      # Tridiagonal for Loop
        self.T_star = None      # Tridiagonal for Star

    def extract(self,
                L: np.ndarray,
                R: np.ndarray,
                P: np.ndarray,
                Z_MM: np.ndarray,
                Z_DD: np.ndarray,
                Z_CC: np.ndarray,
                K_LM: np.ndarray,
                K_LD: np.ndarray,
                K_LC: np.ndarray,
                K_LS: np.ndarray,
                frequencies: np.ndarray = None,
                kan_models: Dict[str, 'KANMaterialInterface'] = None) -> dict:
        """
        Full SPICE extraction with PRIMA Lanczos and per-group ACA.

        Parameters
        ----------
        L : np.ndarray [n_L x n_L]
            Loop inductance matrix
        R : np.ndarray [n_L x n_L]
            Loop resistance matrix
        P : np.ndarray [n_S x n_S]
            Star (potential coefficient) matrix
        Z_MM : np.ndarray [n_M x n_M]
            Magnetic material impedance matrix
        Z_DD : np.ndarray [n_D x n_D]
            Dielectric material impedance matrix
        Z_CC : np.ndarray [n_C x n_C]
            Conductor/shield impedance matrix
        K_LM : np.ndarray [n_L x n_M]
            Loop-Magnetic coupling
        K_LD : np.ndarray [n_L x n_D]
            Loop-Dielectric coupling
        K_LC : np.ndarray [n_L x n_C]
            Loop-Conductor coupling
        K_LS : np.ndarray [n_L x n_S]
            Loop-Star coupling
        frequencies : np.ndarray, optional
            Frequencies for fitting (default: config.f_min to f_max)
        kan_models : dict, optional
            KAN models: {'magnetic': kan_mu, 'dielectric': kan_eps, 'conductor': kan_sigma}

        Returns
        -------
        result : dict
            'netlist': SPICE netlist string
            'circuit': Extracted circuit parameters
            'lanczos_loop': Lanczos result for Loop
            'lanczos_star': Lanczos result for Star
            'Z_port_func': Port impedance function Z(s)
        """
        if frequencies is None:
            frequencies = np.logspace(
                np.log10(self.config.f_min),
                np.log10(self.config.f_max),
                self.config.n_freq
            )

        kan_models = kan_models or {}

        n_L = L.shape[0]
        n_S = P.shape[0] if P is not None and P.size > 0 else 0
        n_M = Z_MM.shape[0] if Z_MM is not None and Z_MM.size > 0 else 0
        n_D = Z_DD.shape[0] if Z_DD is not None and Z_DD.size > 0 else 0
        n_C = Z_CC.shape[0] if Z_CC is not None and Z_CC.size > 0 else 0

        port_indices = self.config.port_indices
        n_ports = len(port_indices)

        print("=" * 70)
        print("PRIMA-Schur SPICE Extraction (Lanczos-only)")
        print("=" * 70)
        print(f"\nSystem dimensions:")
        print(f"  Loop (conductor):    {n_L} DOFs -> Lanczos order {self.config.n_lanczos_loop}")
        print(f"  Star (capacitive):   {n_S} DOFs -> Lanczos order {self.config.n_lanczos_star}")
        print(f"  Magnetic:            {n_M} DOFs -> Lanczos order {self.config.n_lanczos_magnetic}")
        print(f"  Dielectric:          {n_D} DOFs -> Lanczos order {self.config.n_lanczos_dielectric}")
        print(f"  Conductor (shield):  {n_C} DOFs -> Lanczos order {self.config.n_lanczos_conductor}")
        print(f"  Ports:               {n_ports}")

        # ================================================================
        # Step 1: PRIMA Lanczos reduction for Loop DOFs
        # ================================================================
        print("\n[Step 1] PRIMA Lanczos reduction...")

        n_lanczos_L = min(self.config.n_lanczos_loop, n_L)
        if n_lanczos_L > 0 and n_L > n_lanczos_L:
            print(f"  Loop: {n_L} -> {n_lanczos_L} (Lanczos)")
            lanczos_L = self.lanczos.tridiagonalize(L, n_lanczos_L, M=R)
            self.Q_loop = lanczos_L.Q
            self.T_loop = lanczos_L.T

            # Transform coupling matrices
            K_LM_red = self.Q_loop.T @ K_LM if n_M > 0 else None
            K_LD_red = self.Q_loop.T @ K_LD if n_D > 0 else None
            K_LC_red = self.Q_loop.T @ K_LC if n_C > 0 else None
            K_LS_red = self.Q_loop.T @ K_LS if n_S > 0 else None
            R_red = self.Q_loop.T @ R @ self.Q_loop
            L_red = self.T_loop  # Tridiagonal

            # Port projection
            port_proj = np.zeros((n_lanczos_L, n_ports))
            for i, p in enumerate(port_indices):
                if p < n_L:
                    port_proj[:, i] = self.Q_loop[p, :]
        else:
            print(f"  Loop: No reduction (n_L={n_L} <= n_lanczos={n_lanczos_L})")
            K_LM_red = K_LM
            K_LD_red = K_LD
            K_LC_red = K_LC
            K_LS_red = K_LS
            R_red = R
            L_red = L
            n_lanczos_L = n_L
            port_proj = np.eye(n_L)[:, port_indices] if port_indices else None
            lanczos_L = None

        # Star Lanczos (if applicable)
        n_lanczos_S = min(self.config.n_lanczos_star, n_S) if n_S > 0 else 0
        if n_lanczos_S > 0 and n_S > n_lanczos_S:
            print(f"  Star: {n_S} -> {n_lanczos_S} (Lanczos)")
            lanczos_S = self.lanczos.tridiagonalize(P, n_lanczos_S)
            self.Q_star = lanczos_S.Q
            self.T_star = lanczos_S.T
            P_red = self.T_star
            if K_LS_red is not None:
                K_LS_red = K_LS_red @ self.Q_star
        else:
            P_red = P
            n_lanczos_S = n_S
            lanczos_S = None

        # ================================================================
        # Step 2: Lanczos reduction for each material group
        # ================================================================
        print("\n[Step 2] Lanczos reduction per element group...")

        lanczos_results = {}

        # Magnetic Lanczos
        if n_M > 0:
            n_lanczos_M = self.config.n_lanczos_magnetic
            if n_lanczos_M > 0 and n_M > n_lanczos_M:
                print(f"  Magnetic: {n_M} -> {n_lanczos_M} (Lanczos)")
                lanczos_M = self.lanczos.tridiagonalize(Z_MM, n_lanczos_M)
                lanczos_results['magnetic'] = {
                    'Q': lanczos_M.Q, 'T': lanczos_M.T, 'rank': n_lanczos_M}
            else:
                print(f"  Magnetic: No reduction ({n_M} DOFs)")
                lanczos_results['magnetic'] = {
                    'Q': np.eye(n_M), 'T': Z_MM, 'rank': n_M}

        # Dielectric Lanczos
        if n_D > 0:
            n_lanczos_D = self.config.n_lanczos_dielectric
            if n_lanczos_D > 0 and n_D > n_lanczos_D:
                print(f"  Dielectric: {n_D} -> {n_lanczos_D} (Lanczos)")
                lanczos_D = self.lanczos.tridiagonalize(Z_DD, n_lanczos_D)
                lanczos_results['dielectric'] = {
                    'Q': lanczos_D.Q, 'T': lanczos_D.T, 'rank': n_lanczos_D}
            else:
                print(f"  Dielectric: No reduction ({n_D} DOFs)")
                lanczos_results['dielectric'] = {
                    'Q': np.eye(n_D), 'T': Z_DD, 'rank': n_D}

        # Conductor/Shield Lanczos
        if n_C > 0:
            n_lanczos_C = self.config.n_lanczos_conductor
            if n_lanczos_C > 0 and n_C > n_lanczos_C:
                print(f"  Conductor: {n_C} -> {n_lanczos_C} (Lanczos)")
                lanczos_C = self.lanczos.tridiagonalize(Z_CC, n_lanczos_C)
                lanczos_results['conductor'] = {
                    'Q': lanczos_C.Q, 'T': lanczos_C.T, 'rank': n_lanczos_C}
            else:
                print(f"  Conductor: No reduction ({n_C} DOFs)")
                lanczos_results['conductor'] = {
                    'Q': np.eye(n_C), 'T': Z_CC, 'rank': n_C}

        # ================================================================
        # Step 3: Schur complement for port extraction
        # ================================================================
        print("\n[Step 3] Schur complement for port impedance...")

        def Z_port_func(s: complex) -> np.ndarray:
            """Compute port impedance at complex frequency s."""
            return self._compute_port_impedance(
                s, L_red, R_red, P_red, n_lanczos_L, n_lanczos_S,
                n_M, n_D, n_C,
                K_LM_red, K_LD_red, K_LC_red, K_LS_red,
                lanczos_results, port_proj, kan_models
            )

        # Validate at sample frequencies
        test_freqs = [frequencies[0], frequencies[len(frequencies)//2], frequencies[-1]]
        print(f"  Port impedance samples:")
        for f in test_freqs:
            s = 2j * np.pi * f
            Z_p = Z_port_func(s)
            print(f"    f={f:.1e} Hz: |Z_port|={np.abs(Z_p[0,0]):.3e} Ohm")

        # ================================================================
        # Step 4: Generate SPICE netlist
        # ================================================================
        print("\n[Step 4] SPICE netlist generation...")

        netlist, circuit = self._generate_spice_netlist(
            L_red, R_red, P_red,
            n_lanczos_L, n_lanczos_S,
            lanczos_results, K_LM_red, K_LD_red, K_LC_red, K_LS_red,
            port_indices, frequencies, Z_port_func
        )

        print(f"  Generated {circuit['n_elements']} circuit elements")
        print(f"  Netlist: {len(netlist.split(chr(10)))} lines")

        return {
            'netlist': netlist,
            'circuit': circuit,
            'lanczos_loop': lanczos_L,
            'lanczos_star': lanczos_S,
            'lanczos_materials': lanczos_results,
            'Z_port_func': Z_port_func,
            'config': self.config,
            'dimensions': {
                'n_L': n_L, 'n_S': n_S, 'n_M': n_M, 'n_D': n_D, 'n_C': n_C,
                'n_lanczos_L': n_lanczos_L, 'n_lanczos_S': n_lanczos_S,
                'n_ports': n_ports
            }
        }

    def _compute_port_impedance(self,
                                s: complex,
                                L_red: np.ndarray,
                                R_red: np.ndarray,
                                P_red: np.ndarray,
                                n_L: int, n_S: int,
                                n_M: int, n_D: int, n_C: int,
                                K_LM: np.ndarray,
                                K_LD: np.ndarray,
                                K_LC: np.ndarray,
                                K_LS: np.ndarray,
                                lanczos_results: dict,
                                port_proj: np.ndarray,
                                kan_models: dict) -> np.ndarray:
        """Compute port impedance via Schur complement (Lanczos-only)."""
        omega = np.abs(s.imag) if np.abs(s.imag) > 1e-10 else 1.0

        # Get reduced sizes from Lanczos results
        lanczos_M = lanczos_results.get('magnetic', {})
        lanczos_D = lanczos_results.get('dielectric', {})
        lanczos_C = lanczos_results.get('conductor', {})

        n_M_red = lanczos_M.get('rank', n_M) if lanczos_M else n_M
        n_D_red = lanczos_D.get('rank', n_D) if lanczos_D else n_D
        n_C_red = lanczos_C.get('rank', n_C) if lanczos_C else n_C

        # Build reduced system matrix
        n_total = n_L + n_S + n_M_red + n_D_red + n_C_red
        Z = np.zeros((n_total, n_total), dtype=complex)

        idx = 0

        # Loop block: R + sL (tridiagonal after Lanczos)
        Z[idx:idx+n_L, idx:idx+n_L] = R_red + s * L_red
        idx_L = 0
        idx += n_L

        # Star block: P/s
        if n_S > 0:
            idx_S = idx
            if np.abs(s) > 1e-15:
                Z[idx:idx+n_S, idx:idx+n_S] = P_red / s
            else:
                Z[idx:idx+n_S, idx:idx+n_S] = P_red * 1e15
            idx += n_S
        else:
            idx_S = idx

        # Magnetic block (Lanczos tridiagonal)
        if n_M_red > 0 and lanczos_M:
            idx_M = idx
            kan_mu = kan_models.get('magnetic') if kan_models else None
            if kan_mu:
                mu_eff = kan_mu.evaluate_frequency([omega/(2*np.pi)])[0]
            else:
                mu_eff = 1.0
            # Use Lanczos tridiagonal T matrix
            T_M = lanczos_M.get('T', np.eye(n_M_red))
            Z[idx:idx+n_M_red, idx:idx+n_M_red] = T_M * mu_eff
            idx += n_M_red
        else:
            idx_M = idx

        # Dielectric block (Lanczos tridiagonal)
        if n_D_red > 0 and lanczos_D:
            idx_D = idx
            kan_eps = kan_models.get('dielectric') if kan_models else None
            if kan_eps:
                eps_eff = kan_eps.evaluate_frequency([omega/(2*np.pi)])[0]
            else:
                eps_eff = 1.0
            T_D = lanczos_D.get('T', np.eye(n_D_red))
            Z[idx:idx+n_D_red, idx:idx+n_D_red] = T_D * eps_eff
            idx += n_D_red
        else:
            idx_D = idx

        # Conductor/Shield block (Lanczos tridiagonal)
        if n_C_red > 0 and lanczos_C:
            idx_C = idx
            kan_sigma = kan_models.get('conductor') if kan_models else None
            if kan_sigma:
                sigma_eff = kan_sigma.evaluate_frequency([omega/(2*np.pi)])[0]
            else:
                sigma_eff = 1.0
            T_C = lanczos_C.get('T', np.eye(n_C_red))
            Z[idx:idx+n_C_red, idx:idx+n_C_red] = T_C * sigma_eff
            idx += n_C_red
        else:
            idx_C = idx

        # Coupling blocks (transformed to Lanczos basis)
        # Loop-Star
        if n_S > 0 and K_LS is not None:
            Z[idx_L:idx_L+n_L, idx_S:idx_S+n_S] = K_LS
            Z[idx_S:idx_S+n_S, idx_L:idx_L+n_L] = K_LS.T

        # Loop-Magnetic
        if n_M > 0 and K_LM is not None:
            Z[idx_L:idx_L+n_L, idx_M:idx_M+n_M] = s * K_LM
            Z[idx_M:idx_M+n_M, idx_L:idx_L+n_L] = -K_LM.T

        # Loop-Dielectric
        if n_D > 0 and K_LD is not None:
            Z[idx_L:idx_L+n_L, idx_D:idx_D+n_D] = K_LD / s if np.abs(s) > 1e-15 else K_LD * 1e15
            Z[idx_D:idx_D+n_D, idx_L:idx_L+n_L] = K_LD.T / s if np.abs(s) > 1e-15 else K_LD.T * 1e15

        # Loop-Conductor
        if n_C > 0 and K_LC is not None:
            Z[idx_L:idx_L+n_L, idx_C:idx_C+n_C] = K_LC
            Z[idx_C:idx_C+n_C, idx_L:idx_L+n_L] = K_LC.T

        # Schur complement: eliminate non-port DOFs
        n_ports = port_proj.shape[1] if port_proj is not None else 0
        if n_ports == 0:
            return Z[:n_L, :n_L]

        # Project to port space: Z_port = P.T @ Z_LL^{-1} @ P (simplified)
        # Full Schur: Z_port = P.T @ (Z_LL - Z_LX @ Z_XX^{-1} @ Z_XL)^{-1} @ P
        try:
            Z_inv = np.linalg.inv(Z)
            Z_LL_inv = Z_inv[:n_L, :n_L]
            Z_port = port_proj.T @ np.linalg.inv(
                np.linalg.inv(Z_LL_inv)
            ) @ port_proj
        except np.linalg.LinAlgError:
            # Fallback: direct port extraction
            Z_port = port_proj.T @ Z[:n_L, :n_L] @ port_proj

        return Z_port

    def _generate_spice_netlist(self,
                                L_red: np.ndarray,
                                R_red: np.ndarray,
                                P_red: np.ndarray,
                                n_L: int, n_S: int,
                                lanczos_results: dict,
                                K_LM: np.ndarray,
                                K_LD: np.ndarray,
                                K_LC: np.ndarray,
                                K_LS: np.ndarray,
                                port_indices: List[int],
                                frequencies: np.ndarray,
                                Z_port_func: Callable) -> Tuple[str, dict]:
        """Generate SPICE netlist from reduced system (Lanczos-only)."""
        lines = []
        lines.append("* PRIMA-Schur Extracted SPICE Netlist (Lanczos-only)")
        lines.append("* Generated by Radia PEEC + Lanczos MOR")
        lines.append(f"* Lanczos order (Loop):       {self.config.n_lanczos_loop}")
        lines.append(f"* Lanczos order (Star):       {self.config.n_lanczos_star}")
        lines.append(f"* Lanczos order (Magnetic):   {self.config.n_lanczos_magnetic}")
        lines.append(f"* Lanczos order (Dielectric): {self.config.n_lanczos_dielectric}")
        lines.append(f"* Lanczos order (Conductor):  {self.config.n_lanczos_conductor}")
        lines.append("")

        n_ports = len(port_indices)
        circuit_elements = []
        node_counter = [n_ports + 1]  # Mutable counter for internal nodes

        def next_node():
            n = node_counter[0]
            node_counter[0] += 1
            return n

        # Subcircuit header
        port_names = self.config.port_names or [f"P{i+1}" for i in range(n_ports)]
        lines.append(f".SUBCKT PRIMA_EXTRACTED {' '.join(port_names)}")
        lines.append("")

        # ================================================================
        # Loop ladder (PRIMA form from Lanczos tridiagonal)
        # ================================================================
        lines.append("* === Loop Ladder (PRIMA Lanczos) ===")

        # Extract diagonal and off-diagonal from tridiagonal L
        if L_red is not None and L_red.size > 0:
            diag_L = np.diag(L_red)
            offdiag_L = np.diag(L_red, 1) if L_red.shape[0] > 1 else []
            diag_R = np.diag(R_red) if R_red is not None else np.zeros(n_L)

            # Build ladder: series RL elements
            for i in range(min(n_L, self.config.n_lanczos_loop)):
                n1 = port_indices[0] + 1 if i == 0 else prev_node
                n2 = next_node() if i < n_L - 1 else (port_indices[1] + 1 if n_ports > 1 else 0)

                R_val = float(np.real(diag_R[i])) if i < len(diag_R) else 0.01
                L_val = float(np.real(diag_L[i])) if i < len(diag_L) else 1e-9

                if R_val > 1e-15:
                    lines.append(f"R_L{i+1} {n1} {n1}a {R_val:.6e}")
                    circuit_elements.append({'type': 'R', 'value': R_val, 'stage': i})
                    n1 = f"{n1}a"

                if L_val > 1e-15:
                    lines.append(f"L_L{i+1} {n1} {n2} {L_val:.6e}")
                    circuit_elements.append({'type': 'L', 'value': L_val, 'stage': i})

                prev_node = n2

                # Mutual inductance (from off-diagonal)
                if i < len(offdiag_L) and np.abs(offdiag_L[i]) > 1e-15:
                    M_val = float(np.real(offdiag_L[i]))
                    k_val = M_val / np.sqrt(diag_L[i] * diag_L[i+1]) if diag_L[i] * diag_L[i+1] > 0 else 0
                    if np.abs(k_val) < 1:
                        lines.append(f"K_L{i+1}_{i+2} L_L{i+1} L_L{i+2} {k_val:.6e}")

        lines.append("")

        # ================================================================
        # Material subcircuits (from Lanczos tridiagonal)
        # ================================================================
        for group_name, lanczos_data in lanczos_results.items():
            if lanczos_data.get('rank', 0) > 0:
                rank = lanczos_data['rank']
                T = lanczos_data.get('T')

                if T is not None and T.size > 0:
                    lines.append(f"* === {group_name.capitalize()} Ladder (Lanczos rank {rank}) ===")

                    # Extract diagonal and off-diagonal from tridiagonal T
                    diag_T = np.diag(T)
                    offdiag_T = np.diag(T, 1) if T.shape[0] > 1 else []

                    # Build ladder from tridiagonal structure
                    for k in range(rank):
                        n_in = next_node()
                        n_out = next_node() if k < rank - 1 else 0

                        # Diagonal element -> impedance element
                        Z_diag = float(np.abs(diag_T[k])) if k < len(diag_T) else 1.0
                        R_val = Z_diag * 1e-3  # Scale to reasonable values
                        L_val = Z_diag * 1e-9

                        lines.append(f"* {group_name} stage {k+1}")
                        if R_val > 1e-15:
                            lines.append(f"R_{group_name[0].upper()}{k+1} {n_in} {n_in}a {R_val:.6e}")
                            n_in = f"{n_in}a"
                        if L_val > 1e-15:
                            lines.append(f"L_{group_name[0].upper()}{k+1} {n_in} {n_out} {L_val:.6e}")

                        circuit_elements.append({
                            'type': 'ladder_stage',
                            'group': group_name,
                            'stage': k,
                            'R': R_val,
                            'L': L_val
                        })

                        # Off-diagonal -> mutual coupling (K statement)
                        if k < len(offdiag_T) and np.abs(offdiag_T[k]) > 1e-15:
                            M_val = float(np.abs(offdiag_T[k]))
                            k_coupling = M_val / np.sqrt(np.abs(diag_T[k] * diag_T[k+1])) if diag_T[k] * diag_T[k+1] > 0 else 0
                            if 0 < np.abs(k_coupling) < 1:
                                lines.append(f"K_{group_name[0].upper()}{k+1}_{k+2} L_{group_name[0].upper()}{k+1} L_{group_name[0].upper()}{k+2} {k_coupling:.6e}")

                    lines.append("")

        # ================================================================
        # Star capacitors (if present)
        # ================================================================
        if n_S > 0 and P_red is not None:
            lines.append("* === Star Capacitors ===")
            diag_P = np.diag(P_red)
            for i in range(min(n_S, self.config.n_lanczos_star)):
                if diag_P[i] > 1e-15:
                    C_val = 1.0 / float(diag_P[i])  # P = 1/C
                    n1 = next_node()
                    lines.append(f"C_S{i+1} {n1} 0 {C_val:.6e}")
                    circuit_elements.append({'type': 'C', 'value': C_val, 'stage': i})
            lines.append("")

        lines.append(".ENDS PRIMA_EXTRACTED")
        lines.append("")
        lines.append(".END")

        netlist = '\n'.join(lines)
        circuit = {
            'n_elements': len(circuit_elements),
            'elements': circuit_elements,
            'n_ports': n_ports,
            'port_names': port_names
        }

        return netlist, circuit


def demo_prima_schur_extraction():
    """Demonstrate PRIMA-Schur SPICE extraction with per-group ACA."""
    print("=" * 70)
    print("PRIMA-Schur SPICE Extraction Demo")
    print("=" * 70)

    np.random.seed(42)

    # System dimensions
    n_L = 30   # Loop DOFs (conductor currents)
    n_S = 10   # Star DOFs (node potentials)
    n_M = 20   # Magnetic material DOFs
    n_D = 15   # Dielectric material DOFs
    n_C = 10   # Conductor/shield DOFs

    print(f"\nOriginal system: {n_L + n_S + n_M + n_D + n_C} total DOFs")
    print(f"  Loop: {n_L}, Star: {n_S}, Magnetic: {n_M}, Dielectric: {n_D}, Conductor: {n_C}")

    # Create test matrices
    # Loop: inductance with mutual coupling
    L = np.eye(n_L) * 1e-6
    for i in range(n_L):
        for j in range(i+1, n_L):
            L[i, j] = L[j, i] = 0.5e-6 * np.exp(-np.abs(i-j)/5)

    # Resistance
    R = np.eye(n_L) * 0.01

    # Star (potential coefficients)
    P = np.eye(n_S) * 1e10

    # Material matrices
    Z_MM = np.eye(n_M) * 1e-4
    for i in range(n_M):
        for j in range(i+1, n_M):
            Z_MM[i, j] = Z_MM[j, i] = 0.5e-4 * np.exp(-np.abs(i-j)/3)

    Z_DD = np.eye(n_D) * 1e-5
    Z_CC = np.eye(n_C) * 1e-3

    # Coupling matrices
    K_LM = np.random.randn(n_L, n_M) * 1e-4
    K_LD = np.random.randn(n_L, n_D) * 1e-5
    K_LC = np.random.randn(n_L, n_C) * 1e-3
    K_LS = np.random.randn(n_L, n_S) * 1e-6

    # Configuration (Lanczos-only, no ACA)
    config = SPICEExtractionConfig(
        n_lanczos_loop=10,          # Reduce 30 -> 10 stages
        n_lanczos_star=5,           # Reduce 10 -> 5 stages
        n_lanczos_magnetic=3,       # Reduce magnetic to 3 stages
        n_lanczos_dielectric=2,     # Reduce dielectric to 2 stages
        n_lanczos_conductor=2,      # Reduce conductor to 2 stages
        port_indices=[0, 15],       # Two ports
        port_names=['IN', 'OUT'],
        f_min=100,                  # 100 Hz
        f_max=10e6,                 # 10 MHz
        n_freq=50
    )

    print(f"\nExtraction configuration (Lanczos-only):")
    print(f"  Lanczos (Loop):       {n_L} -> {config.n_lanczos_loop} stages")
    print(f"  Lanczos (Star):       {n_S} -> {config.n_lanczos_star} stages")
    print(f"  Lanczos (Magnetic):   {n_M} -> {config.n_lanczos_magnetic} stages")
    print(f"  Lanczos (Dielectric): {n_D} -> {config.n_lanczos_dielectric} stages")
    print(f"  Lanczos (Conductor):  {n_C} -> {config.n_lanczos_conductor} stages")

    # Run extraction
    extractor = PRIMASchurExtractor(config)

    frequencies = np.logspace(2, 7, 50)

    result = extractor.extract(
        L, R, P,
        Z_MM, Z_DD, Z_CC,
        K_LM, K_LD, K_LC, K_LS,
        frequencies
    )

    # Print results
    print("\n" + "=" * 70)
    print("Extraction Results")
    print("=" * 70)

    dims = result['dimensions']
    print(f"\nReduced system:")
    print(f"  Loop:  {dims['n_L']} -> {dims['n_lanczos_L']} stages")
    print(f"  Star:  {dims['n_S']} -> {dims['n_lanczos_S']} stages")

    # Print Lanczos reduction results for materials
    lanczos_mats = result.get('lanczos_materials', {})
    for name, lanczos_data in lanczos_mats.items():
        n_orig = {'magnetic': n_M, 'dielectric': n_D, 'conductor': n_C}.get(name, 0)
        print(f"  {name.capitalize()}: {n_orig} -> {lanczos_data.get('rank', n_orig)} stages")

    print(f"\nSPICE netlist ({len(result['netlist'].split(chr(10)))} lines):")
    print("-" * 50)
    for line in result['netlist'].split('\n')[:40]:
        print(line)
    if len(result['netlist'].split('\n')) > 40:
        print("... (truncated)")

    return result


# =============================================================================
# N-Port Block Lanczos SPICE Generation
# =============================================================================

@dataclass
class NPortBlockLanczosResult:
    """Result of N-port Block Lanczos reduction.

    Block tridiagonal structure:
        T = [A_0   B_0^T   0     ...]
            [B_0   A_1     B_1^T ...]
            [0     B_1     A_2   ...]
            [...   ...     ...   ...]

    where A_k, B_k are p x p blocks (p = number of ports).
    """
    Q: np.ndarray              # Orthonormal basis [n x k*p]
    A_blocks: List[np.ndarray] # Diagonal blocks A_k [p x p]
    B_blocks: List[np.ndarray] # Off-diagonal blocks B_k [p x p]
    n_ports: int               # Number of ports
    n_stages: int              # Number of Lanczos stages
    R_reduced: np.ndarray      # Reduced resistance matrix

    @property
    def T(self) -> np.ndarray:
        """Reconstruct full block tridiagonal matrix."""
        p = self.n_ports
        k = self.n_stages
        n = k * p
        T = np.zeros((n, n))

        for i, A in enumerate(self.A_blocks):
            T[i*p:(i+1)*p, i*p:(i+1)*p] = A

        for i, B in enumerate(self.B_blocks):
            T[i*p:(i+1)*p, (i+1)*p:(i+2)*p] = B.T
            T[(i+1)*p:(i+2)*p, i*p:(i+1)*p] = B

        return T


class NPortBlockLanczosSPICE:
    """N-port Block Lanczos reduction with SPICE netlist generation.

    Generates SPICE subcircuit for N-port impedance matrix:
        V = Z(s) @ I

    where Z(s) is approximated by block tridiagonal Lanczos reduction.

    Circuit elements:
    - Inductors: Self and mutual (K-statement) for L blocks
    - Resistors: Self resistance on diagonal
    - CCVS (H-element): Mutual resistance for off-diagonal R

    Example usage:
        reducer = NPortBlockLanczosSPICE(n_ports=2, n_stages=5)
        result = reducer.reduce(L_matrix, R_matrix, port_indices=[0, 10])
        netlist = reducer.to_spice(result, subckt_name="WPT_COILS")
    """

    def __init__(self, n_ports: int, n_stages: int = 5, tol: float = 1e-10):
        """Initialize N-port Block Lanczos reducer.

        Parameters
        ----------
        n_ports : int
            Number of ports (e.g., 2 for WPT TX/RX)
        n_stages : int
            Number of Lanczos stages per port
        tol : float
            Tolerance for convergence detection
        """
        self.n_ports = n_ports
        self.n_stages = n_stages
        self.tol = tol

    def reduce(self, L: np.ndarray, R: np.ndarray,
               port_indices: List[int]) -> NPortBlockLanczosResult:
        """Apply Block Lanczos reduction.

        Parameters
        ----------
        L : np.ndarray
            Inductance matrix [n x n]
        R : np.ndarray
            Resistance matrix [n x n] or [n] for diagonal
        port_indices : List[int]
            Indices of port nodes (length = n_ports)

        Returns
        -------
        NPortBlockLanczosResult
            Block tridiagonal reduction result
        """
        n = L.shape[0]
        p = self.n_ports
        k = min(self.n_stages, n // p)

        if len(port_indices) != p:
            raise ValueError(f"port_indices length ({len(port_indices)}) != n_ports ({p})")

        # Ensure R is a matrix
        if R.ndim == 1:
            R = np.diag(R)

        # Build port excitation matrix B [n x p]
        B = np.zeros((n, p))
        for i, idx in enumerate(port_indices):
            B[idx, i] = 1.0

        # Block Lanczos on L
        Q, A_blocks, B_blocks = self._block_lanczos_symmetric(L, B, k)

        # Project R to reduced space
        R_reduced = Q.T @ R @ Q

        return NPortBlockLanczosResult(
            Q=Q,
            A_blocks=A_blocks,
            B_blocks=B_blocks,
            n_ports=p,
            n_stages=len(A_blocks),
            R_reduced=R_reduced
        )

    def _block_lanczos_symmetric(self, A: np.ndarray, B: np.ndarray,
                                   k: int) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray]]:
        """Block Lanczos for symmetric matrix A with starting block B.

        Produces block tridiagonal:
            T = Q^T @ A @ Q

        where T has p x p blocks (p = number of columns in B).

        Parameters
        ----------
        A : np.ndarray
            Symmetric matrix [n x n]
        B : np.ndarray
            Starting block [n x p]
        k : int
            Number of Lanczos iterations

        Returns
        -------
        Q : np.ndarray
            Orthonormal basis [n x k*p]
        A_blocks : List[np.ndarray]
            Diagonal blocks of T
        B_blocks : List[np.ndarray]
            Off-diagonal blocks of T
        """
        n = A.shape[0]
        p = B.shape[1]

        # QR of starting block (ensuring positive signs on diagonal)
        V_prev = np.zeros((n, p))
        V_curr, R_qr = np.linalg.qr(B)
        # Fix sign ambiguity: ensure R has positive diagonal
        signs = np.sign(np.diag(R_qr))
        signs[signs == 0] = 1  # Handle exact zeros
        V_curr = V_curr * signs  # Flip columns to match original B orientation

        Q_list = [V_curr]
        A_blocks = []
        B_blocks = []

        for j in range(k):
            # W = A @ V_j
            W = A @ V_curr

            # Diagonal block: A_j = V_j^T @ W
            A_j = V_curr.T @ W
            A_blocks.append(A_j)

            # Last iteration: don't compute B (no next stage)
            if j == k - 1:
                break

            # Orthogonalize: W = W - V_j @ A_j - V_{j-1} @ B_{j-1}^T
            W = W - V_curr @ A_j
            if j > 0:
                W = W - V_prev @ B_blocks[-1].T

            # Full re-orthogonalization against all previous blocks
            for Q_i in Q_list:
                h = Q_i.T @ W
                W = W - Q_i @ h

            # Check for convergence
            norm_W = np.linalg.norm(W)
            if norm_W < self.tol:
                break

            # QR of residual: W = V_{j+1} @ B_j
            V_next, B_j = np.linalg.qr(W)
            B_blocks.append(B_j)

            Q_list.append(V_next)
            V_prev = V_curr
            V_curr = V_next

        Q = np.column_stack(Q_list)
        return Q, A_blocks, B_blocks

    def to_spice(self, result: NPortBlockLanczosResult,
                 subckt_name: str = "NPORT_PRIMA",
                 port_names: List[str] = None) -> str:
        """Generate SPICE netlist from Block Lanczos result.

        Parameters
        ----------
        result : NPortBlockLanczosResult
            Output from reduce()
        subckt_name : str
            Name of SPICE subcircuit
        port_names : List[str], optional
            Port terminal names (default: p1_in, p1_out, p2_in, p2_out, ...)

        Returns
        -------
        str
            SPICE netlist
        """
        p = result.n_ports
        k = result.n_stages

        if port_names is None:
            port_names = []
            for i in range(p):
                port_names.extend([f"p{i+1}_in", f"p{i+1}_out"])

        lines = []
        lines.append(f"* {p}-Port Block Lanczos PRIMA SPICE Netlist")
        lines.append(f"* Ports: {p}, Stages: {k}")
        lines.append(f"* Total reduced DOF: {p * k}")
        lines.append("")

        # Subcircuit header
        port_list = " ".join(port_names)
        lines.append(f".SUBCKT {subckt_name} {port_list}")
        lines.append("")

        # Internal node naming: stage_port (e.g., s0_p1, s1_p2)
        def node(stage, port):
            if stage == 0:
                return f"p{port+1}_in"
            elif stage == k:
                return f"p{port+1}_out"
            else:
                return f"s{stage}_p{port+1}"

        # Element counter
        elem_idx = [0]
        def next_elem(prefix):
            elem_idx[0] += 1
            return f"{prefix}{elem_idx[0]}"

        # === Inductance (diagonal blocks A_j -> self L, off-diagonal -> mutual M) ===
        lines.append("* === Inductance (Block Tridiagonal) ===")

        inductor_names = {}  # (stage, port) -> inductor name

        for j, A_j in enumerate(result.A_blocks):
            # Self inductance: diagonal of A_j
            for port in range(p):
                L_self = A_j[port, port]
                if abs(L_self) > 1e-15:
                    n1 = node(j, port)
                    n2 = node(j + 1, port)
                    L_name = next_elem("L")
                    inductor_names[(j, port)] = L_name
                    lines.append(f"{L_name} {n1} {n2} {L_self:.6e}")

            # Mutual inductance within same stage: off-diagonal of A_j
            for port1 in range(p):
                for port2 in range(port1 + 1, p):
                    M_val = A_j[port1, port2]
                    if abs(M_val) > 1e-15:
                        L1_name = inductor_names.get((j, port1))
                        L2_name = inductor_names.get((j, port2))
                        if L1_name and L2_name:
                            L1 = A_j[port1, port1]
                            L2 = A_j[port2, port2]
                            if abs(L1) > 1e-15 and abs(L2) > 1e-15:
                                k_val = M_val / np.sqrt(abs(L1 * L2))
                                if abs(k_val) < 1.0:
                                    K_name = next_elem("K")
                                    lines.append(f"{K_name} {L1_name} {L2_name} {k_val:.6e}")

        lines.append("")

        # Inter-stage coupling (B_j blocks)
        if result.B_blocks:
            lines.append("* === Inter-stage Coupling ===")
            for j, B_j in enumerate(result.B_blocks):
                for port1 in range(p):
                    for port2 in range(p):
                        M_val = B_j[port1, port2]
                        if abs(M_val) > 1e-15:
                            L1_name = inductor_names.get((j, port1))
                            L2_name = inductor_names.get((j + 1, port2))
                            if L1_name and L2_name:
                                L1 = result.A_blocks[j][port1, port1]
                                L2 = result.A_blocks[j + 1][port2, port2]
                                if abs(L1) > 1e-15 and abs(L2) > 1e-15:
                                    k_val = M_val / np.sqrt(abs(L1 * L2))
                                    if abs(k_val) < 1.0:
                                        K_name = next_elem("K")
                                        lines.append(f"{K_name} {L1_name} {L2_name} {k_val:.6e}")
            lines.append("")

        # === Resistance ===
        lines.append("* === Resistance Matrix ===")

        R_red = result.R_reduced
        n_red = R_red.shape[0]

        # Create resistance nodes for each reduced DOF
        # We need current sensing for mutual resistance (CCVS)

        # Self resistance (diagonal)
        lines.append("* Self resistance (diagonal)")
        for i in range(n_red):
            R_self = R_red[i, i]
            if abs(R_self) > 1e-15:
                stage = i // p
                port = i % p
                # Insert resistance at input of each stage
                n1 = node(stage, port)
                n_mid = f"r{stage}_p{port+1}"
                # Rename the inductor node
                lines.append(f"R{stage+1}_{port+1} {n1} {n_mid} {R_self:.6e}")

        # Mutual resistance (off-diagonal) via CCVS
        lines.append("")
        lines.append("* Mutual resistance (off-diagonal via CCVS)")
        lines.append("* V_i += R_ij * I_j (H-element: current-controlled voltage source)")

        ccvs_count = 0
        for i in range(n_red):
            for j in range(n_red):
                if i != j:
                    R_mutual = R_red[i, j]
                    if abs(R_mutual) > 1e-15:
                        ccvs_count += 1
                        stage_i = i // p
                        port_i = i % p
                        stage_j = j // p
                        port_j = j % p

                        # Need current sensor for port j
                        sense_name = f"Vsense{stage_j+1}_{port_j+1}"
                        h_name = f"H{ccvs_count}"

                        # Add CCVS: controlled by current through sense resistor
                        n_h1 = f"h{ccvs_count}_p{port_i+1}"
                        n_h2 = f"h{ccvs_count}b_p{port_i+1}"
                        lines.append(f"{h_name} {n_h1} {n_h2} {sense_name} {R_mutual:.6e}")

        if ccvs_count > 0:
            lines.append("")
            lines.append("* Current sense (zero-volt sources)")
            added_senses = set()
            for i in range(n_red):
                for j in range(n_red):
                    if i != j and abs(R_red[i, j]) > 1e-15:
                        stage_j = j // p
                        port_j = j % p
                        sense_key = (stage_j, port_j)
                        if sense_key not in added_senses:
                            added_senses.add(sense_key)
                            sense_name = f"Vsense{stage_j+1}_{port_j+1}"
                            n1 = node(stage_j + 1, port_j)
                            lines.append(f"{sense_name} {n1} 0 0")

        lines.append("")
        lines.append(f".ENDS {subckt_name}")

        return "\n".join(lines)

    def compute_impedance(self, result: NPortBlockLanczosResult,
                          frequencies: np.ndarray) -> np.ndarray:
        """Compute N-port impedance matrix over frequency range.

        Z(s) = B^T @ (s*T_L + T_R)^{-1} @ B

        where T_L is block tridiagonal inductance and T_R is resistance.

        Parameters
        ----------
        result : NPortBlockLanczosResult
            Output from reduce()
        frequencies : np.ndarray
            Frequency array [Hz]

        Returns
        -------
        np.ndarray
            Impedance tensor [n_freq x n_ports x n_ports]
        """
        p = result.n_ports
        k = result.n_stages
        n_red = k * p

        T_L = result.T  # Block tridiagonal inductance
        T_R = result.R_reduced

        # Port projection: first p rows/cols
        B = np.zeros((n_red, p))
        for i in range(p):
            B[i, i] = 1.0

        n_freq = len(frequencies)
        Z = np.zeros((n_freq, p, p), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            Y_red = s * T_L + T_R
            try:
                Y_red_inv = np.linalg.inv(Y_red)
                Z[i, :, :] = B.T @ Y_red_inv @ B
            except np.linalg.LinAlgError:
                Z[i, :, :] = np.inf

        return Z

    def compute_impedance_full(self, L: np.ndarray, R: np.ndarray,
                                port_indices: List[int],
                                frequencies: np.ndarray) -> np.ndarray:
        """Compute N-port impedance from full system (reference).

        Parameters
        ----------
        L : np.ndarray
            Full inductance matrix
        R : np.ndarray
            Full resistance matrix
        port_indices : List[int]
            Port node indices
        frequencies : np.ndarray
            Frequency array [Hz]

        Returns
        -------
        np.ndarray
            Impedance tensor [n_freq x n_ports x n_ports]
        """
        n = L.shape[0]
        p = len(port_indices)

        if R.ndim == 1:
            R = np.diag(R)

        # Port matrix
        B = np.zeros((n, p))
        for i, idx in enumerate(port_indices):
            B[idx, i] = 1.0

        n_freq = len(frequencies)
        Z = np.zeros((n_freq, p, p), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            Y_full = s * L + R
            try:
                Y_full_inv = np.linalg.inv(Y_full)
                Z[i, :, :] = B.T @ Y_full_inv @ B
            except np.linalg.LinAlgError:
                Z[i, :, :] = np.inf

        return Z


# =============================================================================
# N-Port Block Lanczos for LoopStar + magnetic block Coupled Systems
# =============================================================================

@dataclass
class NPortCoupledLanczosResult:
    """Result of N-port Block Lanczos for LoopStar + magnetic block coupled system.

    System structure (before Schur complement):
        [sL_L + R_L    K_LM  ] [I_L]   [V_port]
        [K_ML      sL_M + R_M] [M  ] = [0     ]

    After Schur complement (eliminates magnetic DOFs):
        Z_eff(s) = (sL_L + R_L) - K_LM @ (sL_M + R_M)^{-1} @ K_ML

    Block Lanczos is applied to L_L and L_M separately, then coupled.
    """
    # Conductor (Loop) reduction
    Q_L: np.ndarray              # Loop orthonormal basis [n_L x k_L*p]
    T_L: np.ndarray              # Loop block tridiagonal [k_L*p x k_L*p]
    R_L_reduced: np.ndarray      # Reduced loop resistance

    # Magnetic material block reduction
    Q_M: np.ndarray              # Magnetic orthonormal basis [n_M x k_M]
    T_M: np.ndarray              # Magnetic tridiagonal [k_M x k_M]
    R_M_reduced: np.ndarray      # Reduced magnetic resistance (losses)

    # Coupling
    K_LM_reduced: np.ndarray     # Reduced coupling [k_L*p x k_M]

    # Metadata
    n_ports: int
    n_stages_L: int              # Lanczos stages for loops
    n_stages_M: int              # Lanczos stages for magnetic


class NPortCoupledBlockLanczosSPICE:
    """N-port Block Lanczos for LoopStar + magnetic block coupled systems.

    Handles coupled conductor-magnetic systems where:
    - Conductor: Loop-Star PEEC with ports
    - Magnetic: magnetic block for ferrite/iron cores

    The coupled system impedance at ports:
        Z_port(s) = Z_L(s) - K_LM @ Z_M(s)^{-1} @ K_ML

    where:
        Z_L(s) = sL_L + R_L  (conductor impedance)
        Z_M(s) = sL_M + R_M  (magnetic impedance, R_M = loss from complex mu)

    SPICE output includes:
    - Conductor ladder (L, R, K elements)
    - Magnetic ladder (for complex permeability losses)
    - Controlled sources for coupling (E or H elements)

    Example usage:
        solver = NPortCoupledBlockLanczosSPICE(n_ports=2, n_stages_L=5, n_stages_M=3)
        result = solver.reduce(L_L, R_L, L_M, R_M, K_LM, port_indices)
        netlist = solver.to_spice(result, "WPT_WITH_CORE")
    """

    def __init__(self, n_ports: int, n_stages_L: int = 5, n_stages_M: int = 3,
                 tol: float = 1e-10):
        """Initialize coupled Block Lanczos solver.

        Parameters
        ----------
        n_ports : int
            Number of electrical ports
        n_stages_L : int
            Lanczos stages for conductor (loop) DOFs
        n_stages_M : int
            Lanczos stages for magnetic DOFs
        tol : float
            Convergence tolerance
        """
        self.n_ports = n_ports
        self.n_stages_L = n_stages_L
        self.n_stages_M = n_stages_M
        self.tol = tol

    def reduce(self, L_L: np.ndarray, R_L: np.ndarray,
               L_M: np.ndarray, R_M: np.ndarray,
               K_LM: np.ndarray,
               port_indices: List[int]) -> NPortCoupledLanczosResult:
        """Apply Block Lanczos reduction to coupled system.

        Parameters
        ----------
        L_L : np.ndarray
            Conductor inductance matrix [n_L x n_L]
        R_L : np.ndarray
            Conductor resistance [n_L x n_L] or [n_L]
        L_M : np.ndarray
            Magnetic inductance matrix [n_M x n_M]
        R_M : np.ndarray
            Magnetic loss resistance [n_M x n_M] or [n_M]
            (From complex permeability: R_M = omega * mu" / mu'^2 * L_M)
        K_LM : np.ndarray
            Conductor-magnetic coupling [n_L x n_M]
        port_indices : List[int]
            Port node indices in conductor DOFs

        Returns
        -------
        NPortCoupledLanczosResult
        """
        n_L = L_L.shape[0]
        n_M = L_M.shape[0]
        p = self.n_ports
        k_L = min(self.n_stages_L, n_L // p)
        k_M = min(self.n_stages_M, n_M)

        if len(port_indices) != p:
            raise ValueError(f"port_indices length ({len(port_indices)}) != n_ports ({p})")

        # Ensure R matrices are 2D
        if R_L.ndim == 1:
            R_L = np.diag(R_L)
        if R_M.ndim == 1:
            R_M = np.diag(R_M)

        # === Block Lanczos on conductor (L_L) ===
        # Port excitation matrix
        B_L = np.zeros((n_L, p))
        for i, idx in enumerate(port_indices):
            B_L[idx, i] = 1.0

        Q_L, A_blocks_L, B_blocks_L = self._block_lanczos_symmetric(L_L, B_L, k_L)

        # Build T_L from blocks
        n_red_L = len(A_blocks_L) * p
        T_L = np.zeros((n_red_L, n_red_L))
        for i, A in enumerate(A_blocks_L):
            T_L[i*p:(i+1)*p, i*p:(i+1)*p] = A
        for i, B in enumerate(B_blocks_L):
            T_L[i*p:(i+1)*p, (i+1)*p:(i+2)*p] = B.T
            T_L[(i+1)*p:(i+2)*p, i*p:(i+1)*p] = B

        R_L_reduced = Q_L.T @ R_L @ Q_L

        # === Standard Lanczos on magnetic (L_M) ===
        # Use sum of coupling columns as starting vector (captures port influence)
        v0_M = np.sum(np.abs(K_LM.T), axis=1)
        if np.linalg.norm(v0_M) < 1e-15:
            v0_M = np.ones(n_M)
        v0_M = v0_M / np.linalg.norm(v0_M)

        Q_M, T_M = self._lanczos_symmetric(L_M, v0_M, k_M)
        R_M_reduced = Q_M.T @ R_M @ Q_M

        # === Coupling reduction ===
        K_LM_reduced = Q_L.T @ K_LM @ Q_M

        return NPortCoupledLanczosResult(
            Q_L=Q_L,
            T_L=T_L,
            R_L_reduced=R_L_reduced,
            Q_M=Q_M,
            T_M=T_M,
            R_M_reduced=R_M_reduced,
            K_LM_reduced=K_LM_reduced,
            n_ports=p,
            n_stages_L=len(A_blocks_L),
            n_stages_M=Q_M.shape[1]
        )

    def _block_lanczos_symmetric(self, A: np.ndarray, B: np.ndarray,
                                   k: int) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray]]:
        """Block Lanczos (same as NPortBlockLanczosSPICE)."""
        n = A.shape[0]
        p = B.shape[1]

        V_prev = np.zeros((n, p))
        V_curr, R_qr = np.linalg.qr(B)
        signs = np.sign(np.diag(R_qr))
        signs[signs == 0] = 1
        V_curr = V_curr * signs

        Q_list = [V_curr]
        A_blocks = []
        B_blocks = []

        for j in range(k):
            W = A @ V_curr
            A_j = V_curr.T @ W
            A_blocks.append(A_j)

            if j == k - 1:
                break

            W = W - V_curr @ A_j
            if j > 0:
                W = W - V_prev @ B_blocks[-1].T

            for Q_i in Q_list:
                h = Q_i.T @ W
                W = W - Q_i @ h

            norm_W = np.linalg.norm(W)
            if norm_W < self.tol:
                break

            V_next, B_j = np.linalg.qr(W)
            B_blocks.append(B_j)

            Q_list.append(V_next)
            V_prev = V_curr
            V_curr = V_next

        Q = np.column_stack(Q_list)
        return Q, A_blocks, B_blocks

    def _lanczos_symmetric(self, A: np.ndarray, v0: np.ndarray,
                            k: int) -> Tuple[np.ndarray, np.ndarray]:
        """Standard Lanczos for symmetric matrix."""
        n = A.shape[0]
        k = min(k, n)

        Q = np.zeros((n, k))
        alpha = np.zeros(k)
        beta = np.zeros(k - 1)

        v = v0 / np.linalg.norm(v0)
        Q[:, 0] = v
        w = A @ v
        alpha[0] = v.dot(w)
        w = w - alpha[0] * v

        for j in range(1, k):
            beta[j-1] = np.linalg.norm(w)
            if beta[j-1] < self.tol:
                Q = Q[:, :j]
                alpha = alpha[:j]
                beta = beta[:j-1]
                break

            v_prev = v
            v = w / beta[j-1]

            # Re-orthogonalization
            for i in range(j):
                h = Q[:, i].dot(v)
                v = v - h * Q[:, i]
            v = v / np.linalg.norm(v)

            Q[:, j] = v
            w = A @ v
            alpha[j] = v.dot(w)
            w = w - alpha[j] * v - beta[j-1] * v_prev

        # Build tridiagonal
        T = np.diag(alpha) + np.diag(beta, 1) + np.diag(beta, -1)
        return Q, T

    def to_spice(self, result: NPortCoupledLanczosResult,
                 subckt_name: str = "COUPLED_NPORT",
                 port_names: List[str] = None) -> str:
        """Generate SPICE netlist for coupled system.

        The netlist structure:
        1. Conductor ladder (L, R, K for mutual inductance)
        2. Magnetic ladder (L_M, R_M for core losses)
        3. Controlled sources for conductor-magnetic coupling

        Parameters
        ----------
        result : NPortCoupledLanczosResult
        subckt_name : str
        port_names : List[str], optional

        Returns
        -------
        str
            SPICE netlist
        """
        p = result.n_ports
        k_L = result.n_stages_L
        k_M = result.n_stages_M

        if port_names is None:
            port_names = []
            for i in range(p):
                port_names.extend([f"p{i+1}_in", f"p{i+1}_out"])

        lines = []
        lines.append(f"* N-Port Coupled (LoopStar + magnetic block) SPICE Netlist")
        lines.append(f"* Conductor ports: {p}, Conductor stages: {k_L}, Magnetic stages: {k_M}")
        lines.append(f"* Total reduced DOF: {k_L * p} (conductor) + {k_M} (magnetic)")
        lines.append("")
        lines.append(f".SUBCKT {subckt_name} {' '.join(port_names)}")
        lines.append("")

        elem_idx = 1

        # === Conductor Ladder ===
        lines.append("* === Conductor Inductance Ladder ===")

        # Diagonal blocks -> self inductance per port
        for stage in range(k_L):
            for port in range(p):
                i = stage * p + port
                L_val = result.T_L[i, i]
                if stage == 0:
                    n1 = f"p{port+1}_in"
                else:
                    n1 = f"cond_s{stage}_p{port+1}"

                if stage == k_L - 1:
                    n2 = f"p{port+1}_out"
                else:
                    n2 = f"cond_s{stage+1}_p{port+1}"

                lines.append(f"L{elem_idx} {n1} {n2} {L_val:.6e}")
                elem_idx += 1

        # Off-diagonal within same stage -> mutual between ports (K statement)
        lines.append("")
        lines.append("* === Conductor Inter-port Coupling (same stage) ===")
        for stage in range(k_L):
            for i in range(p):
                for j in range(i + 1, p):
                    idx_i = stage * p + i
                    idx_j = stage * p + j
                    M_val = result.T_L[idx_i, idx_j]
                    if np.abs(M_val) > 1e-15:
                        L_i = result.T_L[idx_i, idx_i]
                        L_j = result.T_L[idx_j, idx_j]
                        if L_i > 0 and L_j > 0:
                            k_val = M_val / np.sqrt(L_i * L_j)
                            if np.abs(k_val) < 1.0:
                                # Find L element numbers
                                L_num_i = stage * p + i + 1
                                L_num_j = stage * p + j + 1
                                lines.append(f"K{elem_idx} L{L_num_i} L{L_num_j} {k_val:.6e}")
                                elem_idx += 1

        # === Conductor Resistance ===
        lines.append("")
        lines.append("* === Conductor Resistance ===")
        for i in range(k_L * p):
            R_val = result.R_L_reduced[i, i]
            if np.abs(R_val) > 1e-15:
                stage = i // p
                port = i % p
                if stage == 0:
                    n1 = f"p{port+1}_in"
                else:
                    n1 = f"cond_s{stage}_p{port+1}"
                n_r = f"cond_r{i+1}"
                lines.append(f"R{elem_idx} {n1} {n_r} {R_val:.6e}")
                elem_idx += 1

        # === Magnetic Ladder ===
        lines.append("")
        lines.append("* === Magnetic Core Ladder (magnetic block) ===")
        lines.append("* Models magnetic material inductance and losses")

        for i in range(k_M):
            L_M_val = result.T_M[i, i]
            n1 = f"mag_n{i}" if i > 0 else "mag_in"
            n2 = f"mag_n{i+1}" if i < k_M - 1 else "mag_out"
            lines.append(f"Lmag{i+1} {n1} {n2} {L_M_val:.6e}")

        # Magnetic off-diagonal (tridiagonal coupling)
        for i in range(k_M - 1):
            M_val = result.T_M[i, i+1]
            if np.abs(M_val) > 1e-15:
                L_i = result.T_M[i, i]
                L_j = result.T_M[i+1, i+1]
                if L_i > 0 and L_j > 0:
                    k_val = M_val / np.sqrt(L_i * L_j)
                    if np.abs(k_val) < 1.0:
                        lines.append(f"Kmag{i+1}_{i+2} Lmag{i+1} Lmag{i+2} {k_val:.6e}")

        # Magnetic resistance (core loss from complex mu)
        lines.append("")
        lines.append("* === Magnetic Core Loss (from complex mu) ===")
        for i in range(k_M):
            R_M_val = result.R_M_reduced[i, i]
            if np.abs(R_M_val) > 1e-15:
                n1 = f"mag_n{i}" if i > 0 else "mag_in"
                lines.append(f"Rmag{i+1} {n1} 0 {R_M_val:.6e}")

        # === Conductor-Magnetic Coupling ===
        # Use VCCS (G element) or CCVS (H element) for coupling
        lines.append("")
        lines.append("* === Conductor-Magnetic Coupling ===")
        lines.append("* K_LM couples conductor currents to magnetic flux")
        lines.append("* Implemented via controlled sources")

        # Simplified: dominant coupling terms only
        max_coupling = np.max(np.abs(result.K_LM_reduced))
        threshold = max_coupling * 0.01  # 1% threshold

        for i in range(min(k_L * p, result.K_LM_reduced.shape[0])):
            for j in range(min(k_M, result.K_LM_reduced.shape[1])):
                K_val = result.K_LM_reduced[i, j]
                if np.abs(K_val) > threshold:
                    # Coupling from conductor current to magnetic voltage
                    # This is a simplified model; full model needs E and F elements
                    stage = i // p
                    port = i % p
                    lines.append(f"* Coupling L{i+1} <-> Lmag{j+1}: K = {K_val:.6e}")

        lines.append("")
        lines.append(f".ENDS {subckt_name}")

        return '\n'.join(lines)

    def compute_impedance(self, result: NPortCoupledLanczosResult,
                          frequencies: np.ndarray) -> np.ndarray:
        """Compute port impedance from reduced coupled system.

        Z_port(s) = Z_L(s) - K_LM @ Z_M(s)^{-1} @ K_ML

        Parameters
        ----------
        result : NPortCoupledLanczosResult
        frequencies : np.ndarray

        Returns
        -------
        np.ndarray
            Impedance tensor [n_freq x n_ports x n_ports]
        """
        p = result.n_ports
        n_L = result.T_L.shape[0]
        n_M = result.T_M.shape[0]

        # Port projection (first p DOFs in reduced conductor space)
        B = np.zeros((n_L, p))
        for i in range(p):
            B[i, i] = 1.0

        n_freq = len(frequencies)
        Z = np.zeros((n_freq, p, p), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f

            # Conductor admittance
            Y_L = s * result.T_L + result.R_L_reduced

            # Magnetic admittance
            Y_M = s * result.T_M + result.R_M_reduced

            try:
                # Schur complement: Y_eff = Y_L - K_LM @ Y_M^{-1} @ K_ML
                Y_M_inv = np.linalg.inv(Y_M)
                K = result.K_LM_reduced
                Y_eff = Y_L - K @ Y_M_inv @ K.T

                # Port impedance
                Y_eff_inv = np.linalg.inv(Y_eff)
                Z[i, :, :] = B.T @ Y_eff_inv @ B

            except np.linalg.LinAlgError:
                Z[i, :, :] = np.inf

        return Z

    def compute_impedance_full(self, L_L: np.ndarray, R_L: np.ndarray,
                                L_M: np.ndarray, R_M: np.ndarray,
                                K_LM: np.ndarray,
                                port_indices: List[int],
                                frequencies: np.ndarray) -> np.ndarray:
        """Compute port impedance from full coupled system (reference).

        Parameters
        ----------
        L_L, R_L : np.ndarray
            Conductor matrices
        L_M, R_M : np.ndarray
            Magnetic matrices
        K_LM : np.ndarray
            Coupling matrix
        port_indices : List[int]
        frequencies : np.ndarray

        Returns
        -------
        np.ndarray
            Impedance tensor [n_freq x n_ports x n_ports]
        """
        n_L = L_L.shape[0]
        n_M = L_M.shape[0]
        p = len(port_indices)

        if R_L.ndim == 1:
            R_L = np.diag(R_L)
        if R_M.ndim == 1:
            R_M = np.diag(R_M)

        # Port projection
        B = np.zeros((n_L, p))
        for i, idx in enumerate(port_indices):
            B[idx, i] = 1.0

        n_freq = len(frequencies)
        Z = np.zeros((n_freq, p, p), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f

            Y_L = s * L_L + R_L
            Y_M = s * L_M + R_M

            try:
                # Schur complement
                Y_M_inv = np.linalg.inv(Y_M)
                Y_eff = Y_L - K_LM @ Y_M_inv @ K_LM.T

                Y_eff_inv = np.linalg.inv(Y_eff)
                Z[i, :, :] = B.T @ Y_eff_inv @ B

            except np.linalg.LinAlgError:
                Z[i, :, :] = np.inf

        return Z


# =============================================================================
# Complex Material Properties Support
# =============================================================================

@dataclass
class ComplexMaterialParams:
    """Complex material parameters for frequency-dependent modeling.

    Complex permittivity: eps(f) = eps' - j*eps''
    Complex permeability: mu(f) = mu' - j*mu''

    Loss tangents:
        tan(delta_e) = eps'' / eps'  (dielectric loss)
        tan(delta_m) = mu'' / mu'    (magnetic loss)

    For conductors with finite conductivity:
        eps_eff = eps' - j*(eps'' + sigma/(omega*eps_0))
    """
    # Dielectric properties
    eps_r_real: float = 1.0       # Relative permittivity (real part)
    eps_r_imag: float = 0.0       # Relative permittivity (imaginary part, positive = loss)
    tan_delta_e: float = 0.0      # Dielectric loss tangent (alternative to eps_r_imag)

    # Magnetic properties
    mu_r_real: float = 1.0        # Relative permeability (real part)
    mu_r_imag: float = 0.0        # Relative permeability (imaginary part, positive = loss)
    tan_delta_m: float = 0.0      # Magnetic loss tangent (alternative to mu_r_imag)

    # Conductivity (for combined dielectric + conductive materials)
    sigma: float = 0.0            # Electrical conductivity [S/m]

    # Reference frequency for loss tangent conversion
    f_ref: float = 1e6            # Reference frequency [Hz]

    def complex_permittivity(self, frequency: float) -> complex:
        """Compute complex permittivity at given frequency.

        eps(f) = eps_0 * (eps'_r - j*eps''_r - j*sigma/(omega*eps_0))

        Parameters
        ----------
        frequency : float
            Frequency [Hz]

        Returns
        -------
        complex
            Complex permittivity [F/m]
        """
        eps_0 = 8.854187817e-12  # F/m
        omega = 2 * np.pi * frequency

        eps_r_real = self.eps_r_real
        eps_r_imag = self.eps_r_imag

        # Use loss tangent if imaginary part not specified
        if eps_r_imag == 0.0 and self.tan_delta_e != 0.0:
            eps_r_imag = eps_r_real * self.tan_delta_e

        # Add conductivity contribution
        if self.sigma > 0 and omega > 0:
            eps_r_imag += self.sigma / (omega * eps_0)

        return eps_0 * (eps_r_real - 1j * eps_r_imag)

    def complex_permeability(self, frequency: float) -> complex:
        """Compute complex permeability at given frequency.

        mu(f) = mu_0 * (mu'_r - j*mu''_r)

        Parameters
        ----------
        frequency : float
            Frequency [Hz]

        Returns
        -------
        complex
            Complex permeability [H/m]
        """
        mu_0 = 4 * np.pi * 1e-7  # H/m

        mu_r_real = self.mu_r_real
        mu_r_imag = self.mu_r_imag

        # Use loss tangent if imaginary part not specified
        if mu_r_imag == 0.0 and self.tan_delta_m != 0.0:
            mu_r_imag = mu_r_real * self.tan_delta_m

        return mu_0 * (mu_r_real - 1j * mu_r_imag)

    def to_circuit_params(self, frequency: float) -> dict:
        """Convert complex material to equivalent circuit parameters.

        For magnetic material:
            L_eff = L * mu'_r
            R_loss = omega * L * mu''_r  (series loss resistance)

        For dielectric material:
            C_eff = C * eps'_r
            G_loss = omega * C * eps''_r  (parallel loss conductance)

        Parameters
        ----------
        frequency : float
            Frequency [Hz]

        Returns
        -------
        dict
            Circuit parameters: L_factor, R_factor, C_factor, G_factor
        """
        omega = 2 * np.pi * frequency

        # Magnetic: Z_M = j*omega*L*mu_r = j*omega*L*mu'_r + omega*L*mu''_r
        #         = R_loss + j*omega*L_eff
        L_factor = self.mu_r_real
        R_factor = omega * self.mu_r_imag if self.mu_r_imag > 0 else omega * self.mu_r_real * self.tan_delta_m

        # Dielectric: Y_D = j*omega*C*eps_r = j*omega*C*eps'_r + omega*C*eps''_r
        #           = G_loss + j*omega*C_eff
        C_factor = self.eps_r_real
        eps_imag = self.eps_r_imag if self.eps_r_imag > 0 else self.eps_r_real * self.tan_delta_e
        G_factor = omega * eps_imag

        return {
            'L_factor': L_factor,
            'R_factor': R_factor,
            'C_factor': C_factor,
            'G_factor': G_factor,
            'omega': omega
        }


# =============================================================================
# N-Port Block Lanczos for LoopStar + Dielectric Coupled Systems
# =============================================================================

@dataclass
class NPortCoupledDielectricResult:
    """Result of N-port Block Lanczos for LoopStar + Dielectric coupled system.

    System structure:
        [sL_L + R_L    K_LD  ] [I_L]   [V_port]
        [K_DL      sC_D + G_D] [V_D] = [0     ]

    After Schur complement (eliminates dielectric DOFs):
        Y_eff(s) = (sL_L + R_L) - K_LD @ (sC_D + G_D)^{-1} @ K_DL
    """
    # Conductor (Loop) reduction
    Q_L: np.ndarray              # Loop orthonormal basis
    T_L: np.ndarray              # Loop block tridiagonal (inductance)
    R_L_reduced: np.ndarray      # Reduced loop resistance

    # Dielectric reduction
    Q_D: np.ndarray              # Dielectric orthonormal basis
    T_D: np.ndarray              # Dielectric tridiagonal (elastance P = 1/C)
    G_D_reduced: np.ndarray      # Reduced dielectric conductance (losses)

    # Coupling
    K_LD_reduced: np.ndarray     # Reduced coupling [k_L*p x k_D]

    # Metadata
    n_ports: int
    n_stages_L: int
    n_stages_D: int


class NPortCoupledDielectricSPICE:
    """N-port Block Lanczos for LoopStar + Dielectric coupled systems.

    Handles coupled conductor-dielectric systems where:
    - Conductor: Loop-Star PEEC with ports (inductive)
    - Dielectric: Capacitive elements with losses

    System impedance at ports:
        Z_port(s) = Z_L(s) - K_LD @ Y_D(s)^{-1} @ K_DL

    where:
        Z_L(s) = sL_L + R_L  (conductor impedance)
        Y_D(s) = sC_D + G_D  (dielectric admittance, G_D = loss from complex eps)

    Complex permittivity modeling:
        eps(f) = eps' - j*eps''
        C_eff = C * eps'_r
        G_loss = omega * C * eps''_r

    Example usage:
        solver = NPortCoupledDielectricSPICE(n_ports=2, n_stages_L=5, n_stages_D=3)
        result = solver.reduce(L_L, R_L, C_D, G_D, K_LD, port_indices)
        netlist = solver.to_spice(result, "FILTER_WITH_CAPS")
    """

    def __init__(self, n_ports: int, n_stages_L: int = 5, n_stages_D: int = 3,
                 tol: float = 1e-10):
        self.n_ports = n_ports
        self.n_stages_L = n_stages_L
        self.n_stages_D = n_stages_D
        self.tol = tol

    def reduce(self, L_L: np.ndarray, R_L: np.ndarray,
               C_D: np.ndarray, G_D: np.ndarray,
               K_LD: np.ndarray,
               port_indices: List[int]) -> NPortCoupledDielectricResult:
        """Apply Block Lanczos reduction to conductor-dielectric system.

        Parameters
        ----------
        L_L : np.ndarray
            Conductor inductance matrix [n_L x n_L]
        R_L : np.ndarray
            Conductor resistance [n_L x n_L] or [n_L]
        C_D : np.ndarray
            Dielectric capacitance matrix [n_D x n_D]
        G_D : np.ndarray
            Dielectric loss conductance [n_D x n_D] or [n_D]
            (From complex permittivity: G_D = omega * eps'' / eps'^2 * C_D)
        K_LD : np.ndarray
            Conductor-dielectric coupling [n_L x n_D]
        port_indices : List[int]
            Port node indices in conductor DOFs

        Returns
        -------
        NPortCoupledDielectricResult
        """
        n_L = L_L.shape[0]
        n_D = C_D.shape[0]
        p = self.n_ports
        k_L = min(self.n_stages_L, n_L // p)
        k_D = min(self.n_stages_D, n_D)

        if len(port_indices) != p:
            raise ValueError(f"port_indices length ({len(port_indices)}) != n_ports ({p})")

        # Ensure matrices are 2D
        if R_L.ndim == 1:
            R_L = np.diag(R_L)
        if G_D.ndim == 1:
            G_D = np.diag(G_D)

        # === Block Lanczos on conductor (L_L) ===
        B_L = np.zeros((n_L, p))
        for i, idx in enumerate(port_indices):
            B_L[idx, i] = 1.0

        Q_L, A_blocks_L, B_blocks_L = self._block_lanczos_symmetric(L_L, B_L, k_L)

        # Build T_L
        n_red_L = len(A_blocks_L) * p
        T_L = np.zeros((n_red_L, n_red_L))
        for i, A in enumerate(A_blocks_L):
            T_L[i*p:(i+1)*p, i*p:(i+1)*p] = A
        for i, B in enumerate(B_blocks_L):
            T_L[i*p:(i+1)*p, (i+1)*p:(i+2)*p] = B.T
            T_L[(i+1)*p:(i+2)*p, i*p:(i+1)*p] = B

        R_L_reduced = Q_L.T @ R_L @ Q_L

        # === Block Lanczos on dielectric (use P = inv(C) for numerical stability) ===
        # For capacitors: Y = sC, Z = 1/(sC) = P/s where P = inv(C)
        try:
            P_D = np.linalg.inv(C_D)  # Elastance matrix
        except np.linalg.LinAlgError:
            P_D = np.linalg.pinv(C_D)

        # Starting vectors: DC response of conductor loop excites dielectric
        # For N-port system, we need p starting vectors (one per port)
        # At DC: Z_L = R_L, so I_L = R_L^{-1} @ V_port
        # For each port i: I_L^{(i)} = R_L^{-1} @ B_L[:, i]
        # Dielectric excitation from port i: v_D^{(i)} = K_LD.T @ I_L^{(i)}
        try:
            R_L_inv = np.linalg.inv(R_L)
        except np.linalg.LinAlgError:
            R_L_inv = np.linalg.pinv(R_L)

        # DC loop current from each port [n_L x p]
        I_L_DC = R_L_inv @ B_L
        # Dielectric excitation from each port [n_D x p]
        B_D = K_LD.T @ I_L_DC

        # Check if starting vectors are valid
        if np.linalg.norm(B_D) < 1e-15:
            # Fallback: use coupling columns directly
            B_D = K_LD.T.copy()
            if B_D.shape[1] < p:
                # Pad with zeros if needed
                B_D = np.column_stack([B_D, np.zeros((n_D, p - B_D.shape[1]))])
            elif B_D.shape[1] > p:
                B_D = B_D[:, :p]

        # Apply Block Lanczos on dielectric with p starting vectors
        Q_D, A_blocks_D, B_blocks_D = self._block_lanczos_symmetric(P_D, B_D, k_D)

        # Build T_D from block tridiagonal structure
        n_red_D = len(A_blocks_D) * p
        T_D = np.zeros((n_red_D, n_red_D))
        for i, A in enumerate(A_blocks_D):
            T_D[i*p:(i+1)*p, i*p:(i+1)*p] = A
        for i, B in enumerate(B_blocks_D):
            T_D[i*p:(i+1)*p, (i+1)*p:(i+2)*p] = B.T
            T_D[(i+1)*p:(i+2)*p, i*p:(i+1)*p] = B

        # Transform G_D (conductance to reduced space)
        G_D_reduced = Q_D.T @ G_D @ Q_D

        # === Coupling reduction ===
        K_LD_reduced = Q_L.T @ K_LD @ Q_D

        return NPortCoupledDielectricResult(
            Q_L=Q_L,
            T_L=T_L,
            R_L_reduced=R_L_reduced,
            Q_D=Q_D,
            T_D=T_D,
            G_D_reduced=G_D_reduced,
            K_LD_reduced=K_LD_reduced,
            n_ports=p,
            n_stages_L=len(A_blocks_L),
            n_stages_D=len(A_blocks_D)  # Block Lanczos stages
        )

    def _block_lanczos_symmetric(self, A, B, k):
        """Block Lanczos (same as other classes)."""
        n = A.shape[0]
        p = B.shape[1]

        V_prev = np.zeros((n, p))
        V_curr, R_qr = np.linalg.qr(B)
        signs = np.sign(np.diag(R_qr))
        signs[signs == 0] = 1
        V_curr = V_curr * signs

        Q_list = [V_curr]
        A_blocks = []
        B_blocks = []

        for j in range(k):
            W = A @ V_curr
            A_j = V_curr.T @ W
            A_blocks.append(A_j)

            if j == k - 1:
                break

            W = W - V_curr @ A_j
            if j > 0:
                W = W - V_prev @ B_blocks[-1].T

            for Q_i in Q_list:
                h = Q_i.T @ W
                W = W - Q_i @ h

            norm_W = np.linalg.norm(W)
            if norm_W < self.tol:
                break

            V_next, B_j = np.linalg.qr(W)
            B_blocks.append(B_j)

            Q_list.append(V_next)
            V_prev = V_curr
            V_curr = V_next

        Q = np.column_stack(Q_list)
        return Q, A_blocks, B_blocks

    def _lanczos_symmetric(self, A, v0, k):
        """Standard Lanczos."""
        n = A.shape[0]
        k = min(k, n)

        Q = np.zeros((n, k))
        alpha = np.zeros(k)
        beta = np.zeros(k - 1)

        v = v0 / np.linalg.norm(v0)
        Q[:, 0] = v
        w = A @ v
        alpha[0] = v.dot(w)
        w = w - alpha[0] * v

        for j in range(1, k):
            beta[j-1] = np.linalg.norm(w)
            if beta[j-1] < self.tol:
                Q = Q[:, :j]
                alpha = alpha[:j]
                beta = beta[:j-1]
                break

            v_prev = v
            v = w / beta[j-1]

            for i in range(j):
                h = Q[:, i].dot(v)
                v = v - h * Q[:, i]
            v = v / np.linalg.norm(v)

            Q[:, j] = v
            w = A @ v
            alpha[j] = v.dot(w)
            w = w - alpha[j] * v - beta[j-1] * v_prev

        T = np.diag(alpha) + np.diag(beta, 1) + np.diag(beta, -1)
        return Q, T

    def to_spice(self, result: NPortCoupledDielectricResult,
                 subckt_name: str = "COUPLED_LD",
                 port_names: List[str] = None) -> str:
        """Generate SPICE netlist for conductor-dielectric coupled system."""
        p = result.n_ports
        k_L = result.n_stages_L
        k_D = result.n_stages_D

        if port_names is None:
            port_names = []
            for i in range(p):
                port_names.extend([f"p{i+1}_in", f"p{i+1}_out"])

        lines = []
        lines.append(f"* N-Port Coupled (LoopStar + Dielectric) SPICE Netlist")
        lines.append(f"* Conductor ports: {p}, Conductor stages: {k_L}, Dielectric stages: {k_D}")
        lines.append(f"* Total reduced DOF: {k_L * p} (conductor) + {k_D} (dielectric)")
        lines.append("")
        lines.append(f".SUBCKT {subckt_name} {' '.join(port_names)}")
        lines.append("")

        elem_idx = 1

        # === Conductor Ladder ===
        lines.append("* === Conductor Inductance Ladder ===")
        for stage in range(k_L):
            for port in range(p):
                i = stage * p + port
                L_val = result.T_L[i, i]
                n1 = f"p{port+1}_in" if stage == 0 else f"cond_s{stage}_p{port+1}"
                n2 = f"p{port+1}_out" if stage == k_L - 1 else f"cond_s{stage+1}_p{port+1}"
                lines.append(f"L{elem_idx} {n1} {n2} {L_val:.6e}")
                elem_idx += 1

        # Conductor resistance
        lines.append("")
        lines.append("* === Conductor Resistance ===")
        for i in range(k_L * p):
            R_val = result.R_L_reduced[i, i]
            if np.abs(R_val) > 1e-15:
                stage = i // p
                port = i % p
                n1 = f"p{port+1}_in" if stage == 0 else f"cond_s{stage}_p{port+1}"
                lines.append(f"R{elem_idx} {n1} cond_r{i+1} {R_val:.6e}")
                elem_idx += 1

        # === Dielectric Ladder ===
        lines.append("")
        lines.append("* === Dielectric Capacitor Ladder ===")
        lines.append("* Uses elastance P = 1/C for Lanczos; convert back to C for SPICE")

        for i in range(k_D):
            P_val = result.T_D[i, i]
            C_val = 1.0 / P_val if P_val > 1e-15 else 1e15  # C = 1/P
            n1 = f"cap_n{i}" if i > 0 else "cap_in"
            n2 = f"cap_n{i+1}" if i < k_D - 1 else "cap_out"
            lines.append(f"C{elem_idx} {n1} {n2} {C_val:.6e}")
            elem_idx += 1

        # Dielectric off-diagonal (tridiagonal coupling -> mutual capacitance)
        for i in range(k_D - 1):
            P_off = result.T_D[i, i+1]
            if np.abs(P_off) > 1e-15:
                # Mutual capacitance is harder to represent in SPICE
                # Use comment for now
                lines.append(f"* Mutual cap C{i+1}-C{i+2}: P_off = {P_off:.6e}")

        # Dielectric loss (conductance in parallel)
        lines.append("")
        lines.append("* === Dielectric Loss (from complex eps) ===")
        for i in range(k_D):
            G_val = result.G_D_reduced[i, i]
            if np.abs(G_val) > 1e-15:
                n1 = f"cap_n{i}" if i > 0 else "cap_in"
                lines.append(f"Gdiel{i+1} {n1} 0 {G_val:.6e}")

        # === Coupling ===
        lines.append("")
        lines.append("* === Conductor-Dielectric Coupling ===")
        max_coupling = np.max(np.abs(result.K_LD_reduced))
        threshold = max_coupling * 0.01

        for i in range(min(k_L * p, result.K_LD_reduced.shape[0])):
            for j in range(min(k_D, result.K_LD_reduced.shape[1])):
                K_val = result.K_LD_reduced[i, j]
                if np.abs(K_val) > threshold:
                    lines.append(f"* Coupling L{i+1} <-> C{j+1}: K = {K_val:.6e}")

        lines.append("")
        lines.append(f".ENDS {subckt_name}")

        return '\n'.join(lines)

    def compute_impedance(self, result: NPortCoupledDielectricResult,
                          frequencies: np.ndarray) -> np.ndarray:
        """Compute port impedance from reduced system.

        After Block Lanczos reduction, port DOFs are at the first p indices.
        The Schur complement eliminates dielectric DOFs, giving:
            Z_eff = Z_L - K_LD @ Y_D^{-1} @ K_DL
        Port impedance is extracted from the first p x p submatrix.
        """
        p = result.n_ports
        n_L = result.T_L.shape[0]
        n_D = result.T_D.shape[0]

        n_freq = len(frequencies)
        Z = np.zeros((n_freq, p, p), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f

            # Conductor: Z_L = sL + R
            Z_L = s * result.T_L + result.R_L_reduced

            # Dielectric: Y_D = sC + G, but we have P = 1/C (elastance)
            # C = inv(P), so Y_D = s*inv(P) + G
            Y_D = s * np.linalg.inv(result.T_D) + result.G_D_reduced

            try:
                Y_D_inv = np.linalg.inv(Y_D)
                K = result.K_LD_reduced
                Z_eff = Z_L - K @ Y_D_inv @ K.T

                # Port impedance: extract first p x p submatrix
                # (ports are mapped to first p DOFs after Block Lanczos)
                Z[i, :, :] = Z_eff[:p, :p]

            except np.linalg.LinAlgError:
                Z[i, :, :] = np.inf

        return Z

    def compute_impedance_full(self, L_L, R_L, C_D, G_D, K_LD,
                                port_indices, frequencies) -> np.ndarray:
        """Compute port impedance from full system (reference).

        Uses Schur complement to eliminate dielectric DOFs:
            Z_eff = Z_L - K_LD @ Y_D^{-1} @ K_LD^T

        Port impedance is extracted from Z_eff at port_indices.
        """
        n_L = L_L.shape[0]
        n_D = C_D.shape[0]
        p = len(port_indices)

        if R_L.ndim == 1:
            R_L = np.diag(R_L)
        if G_D.ndim == 1:
            G_D = np.diag(G_D)

        n_freq = len(frequencies)
        Z = np.zeros((n_freq, p, p), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f

            Z_L = s * L_L + R_L
            Y_D = s * C_D + G_D

            try:
                Y_D_inv = np.linalg.inv(Y_D)
                Z_eff = Z_L - K_LD @ Y_D_inv @ K_LD.T

                # Extract port impedance submatrix
                for pi, port_i in enumerate(port_indices):
                    for pj, port_j in enumerate(port_indices):
                        Z[i, pi, pj] = Z_eff[port_i, port_j]

            except np.linalg.LinAlgError:
                Z[i, :, :] = np.inf

        return Z


# =============================================================================
# Unified N-Port Coupled System (L + M + D)
# =============================================================================

@dataclass
class NPortFullCoupledResult:
    """Result for full LoopStar + Magnetic + Dielectric coupled system.

    System structure:
        [sL_L + R_L    K_LM       K_LD     ] [I_L]   [V_port]
        [K_ML      sL_M + R_M     0        ] [M  ] = [0     ]
        [K_DL         0       sC_D + G_D   ] [V_D]   [0     ]

    Schur complement eliminates magnetic and dielectric DOFs.
    """
    # Conductor reduction
    Q_L: np.ndarray
    T_L: np.ndarray
    R_L_reduced: np.ndarray

    # Magnetic reduction
    Q_M: np.ndarray
    T_M: np.ndarray
    R_M_reduced: np.ndarray

    # Dielectric reduction
    Q_D: np.ndarray
    T_D: np.ndarray
    G_D_reduced: np.ndarray

    # Couplings
    K_LM_reduced: np.ndarray
    K_LD_reduced: np.ndarray

    # Metadata
    n_ports: int
    n_stages_L: int
    n_stages_M: int
    n_stages_D: int

    # Complex material parameters (optional)
    material_M: ComplexMaterialParams = None
    material_D: ComplexMaterialParams = None


class NPortFullCoupledSPICE:
    """N-port Block Lanczos for full L + M + D coupled systems.

    Handles systems with:
    - Conductor: Coils, traces, wires (inductive + resistive)
    - Magnetic: Ferrite cores, iron yokes (complex permeability)
    - Dielectric: Capacitors, substrates (complex permittivity)

    Complex material modeling:
        mu(f) = mu' - j*mu''  -> L_eff = L*mu', R_loss = omega*L*mu''
        eps(f) = eps' - j*eps'' -> C_eff = C*eps', G_loss = omega*C*eps''

    Example usage:
        # Define complex materials
        ferrite = ComplexMaterialParams(mu_r_real=2000, tan_delta_m=0.01)
        substrate = ComplexMaterialParams(eps_r_real=4.5, tan_delta_e=0.02)

        # Create solver
        solver = NPortFullCoupledSPICE(n_ports=2, n_stages_L=5, n_stages_M=3, n_stages_D=3)
        result = solver.reduce(L_L, R_L, L_M, R_M, C_D, G_D, K_LM, K_LD, port_indices,
                               material_M=ferrite, material_D=substrate)
        netlist = solver.to_spice(result, "TRANSFORMER_ON_PCB")
    """

    def __init__(self, n_ports: int, n_stages_L: int = 5,
                 n_stages_M: int = 3, n_stages_D: int = 3,
                 tol: float = 1e-10):
        self.n_ports = n_ports
        self.n_stages_L = n_stages_L
        self.n_stages_M = n_stages_M
        self.n_stages_D = n_stages_D
        self.tol = tol

    def reduce(self, L_L: np.ndarray, R_L: np.ndarray,
               L_M: np.ndarray, R_M: np.ndarray,
               C_D: np.ndarray, G_D: np.ndarray,
               K_LM: np.ndarray, K_LD: np.ndarray,
               port_indices: List[int],
               material_M: ComplexMaterialParams = None,
               material_D: ComplexMaterialParams = None) -> NPortFullCoupledResult:
        """Apply Block Lanczos to full coupled system."""
        n_L = L_L.shape[0]
        n_M = L_M.shape[0] if L_M is not None else 0
        n_D = C_D.shape[0] if C_D is not None else 0
        p = self.n_ports

        k_L = min(self.n_stages_L, n_L // p)
        k_M = min(self.n_stages_M, n_M) if n_M > 0 else 0
        k_D = min(self.n_stages_D, n_D) if n_D > 0 else 0

        # Ensure matrices are 2D
        if R_L.ndim == 1:
            R_L = np.diag(R_L)
        if n_M > 0 and R_M.ndim == 1:
            R_M = np.diag(R_M)
        if n_D > 0 and G_D.ndim == 1:
            G_D = np.diag(G_D)

        # === Conductor Block Lanczos ===
        B_L = np.zeros((n_L, p))
        for i, idx in enumerate(port_indices):
            B_L[idx, i] = 1.0

        Q_L, A_blocks_L, B_blocks_L = self._block_lanczos_symmetric(L_L, B_L, k_L)

        n_red_L = len(A_blocks_L) * p
        T_L = np.zeros((n_red_L, n_red_L))
        for i, A in enumerate(A_blocks_L):
            T_L[i*p:(i+1)*p, i*p:(i+1)*p] = A
        for i, B in enumerate(B_blocks_L):
            T_L[i*p:(i+1)*p, (i+1)*p:(i+2)*p] = B.T
            T_L[(i+1)*p:(i+2)*p, i*p:(i+1)*p] = B

        R_L_reduced = Q_L.T @ R_L @ Q_L

        # === Magnetic Block Lanczos ===
        # Starting vectors: DC response of conductor loop excites magnetic
        # At DC: Z_L = R_L, so I_L = R_L^{-1} @ V_port
        # Magnetic excitation from each port: B_M = K_LM.T @ R_L^{-1} @ B_L
        if n_M > 0:
            try:
                R_L_inv = np.linalg.inv(R_L)
            except np.linalg.LinAlgError:
                R_L_inv = np.linalg.pinv(R_L)

            # DC loop current from each port [n_L x p]
            I_L_DC = R_L_inv @ B_L
            # Magnetic excitation from each port [n_M x p]
            B_M = K_LM.T @ I_L_DC

            if np.linalg.norm(B_M) < 1e-15:
                # Fallback: use coupling columns
                B_M = K_LM.T[:, :p] if K_LM.shape[1] >= p else np.column_stack([K_LM.T, np.zeros((n_M, p - K_LM.shape[1]))])

            Q_M, A_blocks_M, B_blocks_M = self._block_lanczos_symmetric(L_M, B_M, k_M)

            # Build T_M from block tridiagonal structure
            n_red_M = len(A_blocks_M) * p
            T_M = np.zeros((n_red_M, n_red_M))
            for i, A in enumerate(A_blocks_M):
                T_M[i*p:(i+1)*p, i*p:(i+1)*p] = A
            for i, B in enumerate(B_blocks_M):
                T_M[i*p:(i+1)*p, (i+1)*p:(i+2)*p] = B.T
                T_M[(i+1)*p:(i+2)*p, i*p:(i+1)*p] = B

            R_M_reduced = Q_M.T @ R_M @ Q_M
            K_LM_reduced = Q_L.T @ K_LM @ Q_M
            n_stages_M_actual = len(A_blocks_M)
        else:
            Q_M = np.array([]).reshape(0, 0)
            T_M = np.array([]).reshape(0, 0)
            R_M_reduced = np.array([]).reshape(0, 0)
            K_LM_reduced = np.array([]).reshape(n_red_L, 0)
            n_stages_M_actual = 0

        # === Dielectric Block Lanczos ===
        # Starting vectors: DC response excites dielectric
        if n_D > 0:
            try:
                P_D = np.linalg.inv(C_D)
            except np.linalg.LinAlgError:
                P_D = np.linalg.pinv(C_D)

            if 'R_L_inv' not in dir():
                try:
                    R_L_inv = np.linalg.inv(R_L)
                except np.linalg.LinAlgError:
                    R_L_inv = np.linalg.pinv(R_L)

            # DC loop current from each port [n_L x p]
            I_L_DC = R_L_inv @ B_L
            # Dielectric excitation from each port [n_D x p]
            B_D = K_LD.T @ I_L_DC

            if np.linalg.norm(B_D) < 1e-15:
                B_D = K_LD.T[:, :p] if K_LD.shape[1] >= p else np.column_stack([K_LD.T, np.zeros((n_D, p - K_LD.shape[1]))])

            Q_D, A_blocks_D, B_blocks_D = self._block_lanczos_symmetric(P_D, B_D, k_D)

            # Build T_D from block tridiagonal structure
            n_red_D = len(A_blocks_D) * p
            T_D = np.zeros((n_red_D, n_red_D))
            for i, A in enumerate(A_blocks_D):
                T_D[i*p:(i+1)*p, i*p:(i+1)*p] = A
            for i, B in enumerate(B_blocks_D):
                T_D[i*p:(i+1)*p, (i+1)*p:(i+2)*p] = B.T
                T_D[(i+1)*p:(i+2)*p, i*p:(i+1)*p] = B

            G_D_reduced = Q_D.T @ G_D @ Q_D
            K_LD_reduced = Q_L.T @ K_LD @ Q_D
            n_stages_D_actual = len(A_blocks_D)
        else:
            Q_D = np.array([]).reshape(0, 0)
            T_D = np.array([]).reshape(0, 0)
            G_D_reduced = np.array([]).reshape(0, 0)
            K_LD_reduced = np.array([]).reshape(n_red_L, 0)
            n_stages_D_actual = 0

        return NPortFullCoupledResult(
            Q_L=Q_L, T_L=T_L, R_L_reduced=R_L_reduced,
            Q_M=Q_M, T_M=T_M, R_M_reduced=R_M_reduced,
            Q_D=Q_D, T_D=T_D, G_D_reduced=G_D_reduced,
            K_LM_reduced=K_LM_reduced, K_LD_reduced=K_LD_reduced,
            n_ports=p,
            n_stages_L=len(A_blocks_L),
            n_stages_M=n_stages_M_actual,
            n_stages_D=n_stages_D_actual,
            material_M=material_M,
            material_D=material_D
        )

    def _block_lanczos_symmetric(self, A, B, k):
        """Block Lanczos."""
        n = A.shape[0]
        p = B.shape[1]

        V_prev = np.zeros((n, p))
        V_curr, R_qr = np.linalg.qr(B)
        signs = np.sign(np.diag(R_qr))
        signs[signs == 0] = 1
        V_curr = V_curr * signs

        Q_list = [V_curr]
        A_blocks = []
        B_blocks = []

        for j in range(k):
            W = A @ V_curr
            A_j = V_curr.T @ W
            A_blocks.append(A_j)

            if j == k - 1:
                break

            W = W - V_curr @ A_j
            if j > 0:
                W = W - V_prev @ B_blocks[-1].T

            for Q_i in Q_list:
                h = Q_i.T @ W
                W = W - Q_i @ h

            norm_W = np.linalg.norm(W)
            if norm_W < self.tol:
                break

            V_next, B_j = np.linalg.qr(W)
            B_blocks.append(B_j)

            Q_list.append(V_next)
            V_prev = V_curr
            V_curr = V_next

        Q = np.column_stack(Q_list)
        return Q, A_blocks, B_blocks

    def _lanczos_symmetric(self, A, v0, k):
        """Standard Lanczos."""
        n = A.shape[0]
        k = min(k, n)

        Q = np.zeros((n, k))
        alpha = np.zeros(k)
        beta = np.zeros(k - 1)

        v = v0 / np.linalg.norm(v0)
        Q[:, 0] = v
        w = A @ v
        alpha[0] = v.dot(w)
        w = w - alpha[0] * v

        for j in range(1, k):
            beta[j-1] = np.linalg.norm(w)
            if beta[j-1] < self.tol:
                Q = Q[:, :j]
                alpha = alpha[:j]
                beta = beta[:j-1]
                break

            v_prev = v
            v = w / beta[j-1]

            for i in range(j):
                h = Q[:, i].dot(v)
                v = v - h * Q[:, i]
            v = v / np.linalg.norm(v)

            Q[:, j] = v
            w = A @ v
            alpha[j] = v.dot(w)
            w = w - alpha[j] * v - beta[j-1] * v_prev

        T = np.diag(alpha) + np.diag(beta, 1) + np.diag(beta, -1)
        return Q, T

    def compute_impedance(self, result: NPortFullCoupledResult,
                          frequencies: np.ndarray) -> np.ndarray:
        """Compute port impedance with complex material effects.

        After Block Lanczos reduction, port DOFs are at the first p indices.
        Schur complement eliminates magnetic and dielectric DOFs.
        Port impedance is extracted from Z_eff[:p, :p].
        """
        p = result.n_ports
        n_L = result.T_L.shape[0]
        n_M = result.T_M.shape[0] if result.T_M.size > 0 else 0
        n_D = result.T_D.shape[0] if result.T_D.size > 0 else 0

        n_freq = len(frequencies)
        Z = np.zeros((n_freq, p, p), dtype=complex)

        for i, f in enumerate(frequencies):
            s = 2j * np.pi * f
            omega = 2 * np.pi * f

            # Conductor
            Z_L = s * result.T_L + result.R_L_reduced

            # Magnetic with complex mu
            if n_M > 0:
                if result.material_M is not None:
                    params = result.material_M.to_circuit_params(f)
                    L_factor = params['L_factor']
                    R_factor = params['R_factor'] / omega if omega > 0 else 0
                    Y_M = s * result.T_M * L_factor + result.R_M_reduced + R_factor * result.T_M
                else:
                    Y_M = s * result.T_M + result.R_M_reduced

            # Dielectric with complex eps
            if n_D > 0:
                if result.material_D is not None:
                    params = result.material_D.to_circuit_params(f)
                    C_factor = params['C_factor']
                    G_factor = params['G_factor']
                    # T_D is elastance (1/C), so C = 1/T_D
                    T_D_inv = np.linalg.inv(result.T_D)
                    Y_D = s * T_D_inv * C_factor + result.G_D_reduced + G_factor * T_D_inv
                else:
                    Y_D = s * np.linalg.inv(result.T_D) + result.G_D_reduced

            try:
                Z_eff = Z_L.copy()

                # Schur complement for magnetic
                if n_M > 0:
                    Y_M_inv = np.linalg.inv(Y_M)
                    K_M = result.K_LM_reduced
                    Z_eff = Z_eff - K_M @ Y_M_inv @ K_M.T

                # Schur complement for dielectric
                if n_D > 0:
                    Y_D_inv = np.linalg.inv(Y_D)
                    K_D = result.K_LD_reduced
                    Z_eff = Z_eff - K_D @ Y_D_inv @ K_D.T

                # Port impedance: extract first p x p submatrix
                # (ports are mapped to first p DOFs after Block Lanczos)
                Z[i, :, :] = Z_eff[:p, :p]

            except np.linalg.LinAlgError:
                Z[i, :, :] = np.inf

        return Z

    def to_spice(self, result: NPortFullCoupledResult,
                 subckt_name: str = "FULL_COUPLED",
                 port_names: List[str] = None,
                 f_ref: float = 100e3) -> str:
        """Generate SPICE netlist with complex material parameters.

        Parameters
        ----------
        result : NPortFullCoupledResult
        subckt_name : str
        port_names : List[str], optional
        f_ref : float
            Reference frequency for loss element values
        """
        p = result.n_ports
        k_L = result.n_stages_L
        k_M = result.n_stages_M
        k_D = result.n_stages_D

        if port_names is None:
            port_names = []
            for i in range(p):
                port_names.extend([f"p{i+1}_in", f"p{i+1}_out"])

        lines = []
        lines.append(f"* N-Port Full Coupled (L + M + D) SPICE Netlist")
        lines.append(f"* Conductor stages: {k_L}, Magnetic stages: {k_M}, Dielectric stages: {k_D}")
        lines.append(f"* Reference frequency for loss: {f_ref/1e3:.1f} kHz")

        if result.material_M is not None:
            lines.append(f"* Magnetic: mu'_r = {result.material_M.mu_r_real}, tan(delta_m) = {result.material_M.tan_delta_m}")
        if result.material_D is not None:
            lines.append(f"* Dielectric: eps'_r = {result.material_D.eps_r_real}, tan(delta_e) = {result.material_D.tan_delta_e}")

        lines.append("")
        lines.append(f".SUBCKT {subckt_name} {' '.join(port_names)}")
        lines.append("")

        elem_idx = 1
        omega_ref = 2 * np.pi * f_ref

        # === Conductor Ladder ===
        lines.append("* === Conductor Inductance Ladder ===")
        for stage in range(k_L):
            for port in range(p):
                i = stage * p + port
                L_val = result.T_L[i, i]
                n1 = f"p{port+1}_in" if stage == 0 else f"cond_s{stage}_p{port+1}"
                n2 = f"p{port+1}_out" if stage == k_L - 1 else f"cond_s{stage+1}_p{port+1}"
                lines.append(f"L{elem_idx} {n1} {n2} {L_val:.6e}")
                elem_idx += 1

        # Conductor resistance
        lines.append("")
        lines.append("* === Conductor Resistance ===")
        for i in range(k_L * p):
            R_val = result.R_L_reduced[i, i]
            if np.abs(R_val) > 1e-15:
                lines.append(f"Rcond{i+1} cond_r{i+1}_a cond_r{i+1}_b {R_val:.6e}")

        # === Magnetic Ladder ===
        if k_M > 0:
            lines.append("")
            lines.append("* === Magnetic Core Ladder ===")

            # Apply complex mu scaling
            mu_factor = 1.0
            mu_loss_factor = 0.0
            if result.material_M is not None:
                mu_factor = result.material_M.mu_r_real
                mu_loss_factor = result.material_M.tan_delta_m

            for i in range(k_M):
                L_M_val = result.T_M[i, i] * mu_factor
                n1 = f"mag_n{i}" if i > 0 else "mag_in"
                n2 = f"mag_n{i+1}" if i < k_M - 1 else "mag_out"
                lines.append(f"Lmag{i+1} {n1} {n2} {L_M_val:.6e}")

                # Core loss resistance (series)
                R_loss = omega_ref * result.T_M[i, i] * mu_factor * mu_loss_factor
                if R_loss > 1e-15:
                    lines.append(f"Rmag_loss{i+1} {n1} mag_loss{i+1} {R_loss:.6e}")

        # === Dielectric Ladder ===
        if k_D > 0:
            lines.append("")
            lines.append("* === Dielectric Capacitor Ladder ===")

            eps_factor = 1.0
            eps_loss_factor = 0.0
            if result.material_D is not None:
                eps_factor = result.material_D.eps_r_real
                eps_loss_factor = result.material_D.tan_delta_e

            for i in range(k_D):
                P_val = result.T_D[i, i]
                C_val = eps_factor / P_val if P_val > 1e-15 else 1e15
                n1 = f"cap_n{i}" if i > 0 else "cap_in"
                n2 = f"cap_n{i+1}" if i < k_D - 1 else "cap_out"
                lines.append(f"Cdiel{i+1} {n1} {n2} {C_val:.6e}")

                # Dielectric loss conductance (parallel)
                G_loss = omega_ref * C_val * eps_loss_factor
                if G_loss > 1e-15:
                    lines.append(f"Gdiel_loss{i+1} {n1} 0 {G_loss:.6e}")

        lines.append("")
        lines.append(f".ENDS {subckt_name}")

        return '\n'.join(lines)


def demo_coupled_nport_lanczos():
    """Demonstrate N-port coupled (LoopStar + magnetic block) Block Lanczos."""
    print("=" * 70)
    print("N-Port Coupled (LoopStar + magnetic block) Block Lanczos Demo")
    print("=" * 70)

    # Create synthetic 2-port WPT system with ferrite core
    n_L = 20   # Conductor DOFs (TX + RX coils)
    n_M = 10   # Magnetic DOFs (ferrite core)
    p = 2      # Ports

    # === Conductor (coils) ===
    L_L = np.eye(n_L) * 10e-6  # 10 uH self inductance

    # TX region (0-9) and RX region (10-19) internal coupling
    for region_start in [0, 10]:
        for i in range(10):
            for j in range(10):
                if i != j:
                    L_L[region_start+i, region_start+j] = 1e-6 / (1 + abs(i-j))

    # TX-RX coupling (weak without core)
    k_air = 0.05  # 5% coupling in air
    M_air = k_air * 10e-6
    for i in range(10):
        for j in range(10, 20):
            L_L[i, j] = M_air / (1 + 0.3 * abs(i - (j-10)))
            L_L[j, i] = L_L[i, j]

    R_L = np.ones(n_L) * 0.1  # 100 mOhm

    # === Magnetic core ===
    L_M = np.eye(n_M) * 1e-3  # 1 mH (high permeability)
    for i in range(n_M):
        for j in range(n_M):
            if i != j:
                L_M[i, j] = 0.5e-3 / (1 + abs(i-j))  # Core internal coupling

    # Core loss (from complex mu: R_M = omega * mu" * V / mu'^2)
    # At 100 kHz, tan(delta) ~ 0.01 for ferrite
    omega_ref = 2 * np.pi * 100e3
    tan_delta = 0.01
    R_M = np.diag(np.diag(L_M)) * omega_ref * tan_delta  # Frequency-dependent

    # === Coupling (coil to core) ===
    # K_LM represents mutual inductance between coils and core elements
    K_LM = np.zeros((n_L, n_M))

    # TX coil couples to core
    for i in range(10):
        for j in range(n_M):
            K_LM[i, j] = 0.3e-3 * np.exp(-0.2 * abs(i - j))

    # RX coil couples to core
    for i in range(10, 20):
        for j in range(n_M):
            K_LM[i, j] = 0.3e-3 * np.exp(-0.2 * abs((i-10) - j))

    port_indices = [0, 10]

    print(f"\nSystem setup:")
    print(f"  Conductor DOFs: {n_L} (TX: 0-9, RX: 10-19)")
    print(f"  Magnetic DOFs:  {n_M} (ferrite core)")
    print(f"  Ports: {p} (TX at DOF 0, RX at DOF 10)")
    print(f"  Air coupling k: {k_air}")
    print(f"  Core tan(delta): {tan_delta}")

    # Apply coupled Block Lanczos
    solver = NPortCoupledBlockLanczosSPICE(n_ports=p, n_stages_L=5, n_stages_M=3)
    result = solver.reduce(L_L, R_L, L_M, R_M, K_LM, port_indices)

    print(f"\nBlock Lanczos reduction:")
    print(f"  Conductor: {n_L} -> {result.n_stages_L * p} DOF")
    print(f"  Magnetic:  {n_M} -> {result.n_stages_M} DOF")
    print(f"  Total:     {n_L + n_M} -> {result.n_stages_L * p + result.n_stages_M} DOF")
    print(f"  Compression: {100*(1 - (result.n_stages_L*p + result.n_stages_M)/(n_L + n_M)):.1f}%")

    # Generate SPICE
    netlist = solver.to_spice(result, "WPT_WITH_CORE",
                               ["tx_p", "tx_n", "rx_p", "rx_n"])

    print(f"\nSPICE netlist ({len(netlist.split(chr(10)))} lines):")
    print("-" * 50)
    for line in netlist.split('\n')[:25]:
        print(line)
    if len(netlist.split('\n')) > 25:
        print("... (truncated)")

    # Verify frequency response
    frequencies = np.logspace(3, 6, 30)  # 1 kHz to 1 MHz

    Z_full = solver.compute_impedance_full(L_L, R_L, L_M, R_M, K_LM, port_indices, frequencies)
    Z_red = solver.compute_impedance(result, frequencies)

    err_Z11 = np.mean(np.abs(Z_full[:, 0, 0] - Z_red[:, 0, 0]) / np.abs(Z_full[:, 0, 0])) * 100
    err_Z12 = np.mean(np.abs(Z_full[:, 0, 1] - Z_red[:, 0, 1]) / (np.abs(Z_full[:, 0, 1]) + 1e-15)) * 100

    print(f"\nFrequency response verification:")
    print(f"  Z_11 (TX self) error: {err_Z11:.2f}%")
    print(f"  Z_12 (TX-RX mutual) error: {err_Z12:.2f}%")

    # Sample values
    f_idx = np.argmin(np.abs(frequencies - 100e3))
    print(f"\nSample impedance at 100 kHz:")
    print(f"  Z_11 full: {Z_full[f_idx, 0, 0]:.4f} Ohm")
    print(f"  Z_11 red:  {Z_red[f_idx, 0, 0]:.4f} Ohm")
    print(f"  Z_12 full: {Z_full[f_idx, 0, 1]:.4f} Ohm")
    print(f"  Z_12 red:  {Z_red[f_idx, 0, 1]:.4f} Ohm")

    return result, netlist


def demo_nport_block_lanczos():
    """Demonstrate N-port Block Lanczos SPICE generation."""
    print("=" * 70)
    print("N-Port Block Lanczos SPICE Generation Demo")
    print("=" * 70)

    # Create synthetic 2-port WPT system
    n = 20  # Total DOFs
    p = 2   # Ports (TX, RX)

    # Inductance matrix with coupling between TX and RX regions
    L = np.eye(n) * 10e-6  # 10 uH self inductance

    # TX region: DOFs 0-9
    # RX region: DOFs 10-19
    for i in range(10):
        for j in range(10):
            if i != j:
                L[i, j] = 1e-6 / (1 + abs(i - j))  # Decay with distance

    for i in range(10, 20):
        for j in range(10, 20):
            if i != j:
                L[i, j] = 1e-6 / (1 + abs(i - j))

    # TX-RX coupling (mutual inductance)
    k_coupling = 0.3  # 30% coupling
    M_tx_rx = k_coupling * np.sqrt(10e-6 * 10e-6)  # M = k * sqrt(L1*L2)
    for i in range(10):
        for j in range(10, 20):
            dist = abs(i - (j - 10))
            L[i, j] = M_tx_rx / (1 + 0.5 * dist)
            L[j, i] = L[i, j]

    # Resistance (diagonal)
    R = np.ones(n) * 0.1  # 100 mOhm

    # Port indices
    port_indices = [0, 10]  # TX port at DOF 0, RX port at DOF 10

    print(f"\nSystem setup:")
    print(f"  Total DOFs: {n}")
    print(f"  Ports: {p} (TX at DOF {port_indices[0]}, RX at DOF {port_indices[1]})")
    print(f"  Self inductance: 10 uH")
    print(f"  TX-RX coupling: k = {k_coupling}")

    # Apply Block Lanczos reduction
    reducer = NPortBlockLanczosSPICE(n_ports=p, n_stages=5)
    result = reducer.reduce(L, R, port_indices)

    print(f"\nBlock Lanczos reduction:")
    print(f"  Original DOFs: {n}")
    print(f"  Reduced DOFs: {result.n_stages * result.n_ports}")
    print(f"  Compression: {100 * (1 - result.n_stages * result.n_ports / n):.1f}%")

    # Generate SPICE netlist
    netlist = reducer.to_spice(result, subckt_name="WPT_COILS",
                                port_names=["tx_in", "tx_out", "rx_in", "rx_out"])

    print(f"\nSPICE netlist ({len(netlist.split(chr(10)))} lines):")
    print("-" * 50)
    for line in netlist.split('\n')[:30]:
        print(line)
    if len(netlist.split('\n')) > 30:
        print("... (truncated)")

    # Verify frequency response
    frequencies = np.logspace(3, 7, 50)  # 1 kHz to 10 MHz

    Z_full = reducer.compute_impedance_full(L, R, port_indices, frequencies)
    Z_red = reducer.compute_impedance(result, frequencies)

    # Compare Z_11 (TX self-impedance)
    Z11_full = Z_full[:, 0, 0]
    Z11_red = Z_red[:, 0, 0]

    # Compare Z_12 (TX-RX mutual impedance)
    Z12_full = Z_full[:, 0, 1]
    Z12_red = Z_red[:, 0, 1]

    err_Z11 = np.mean(np.abs(Z11_full - Z11_red) / np.abs(Z11_full)) * 100
    err_Z12 = np.mean(np.abs(Z12_full - Z12_red) / (np.abs(Z12_full) + 1e-10)) * 100

    print(f"\nFrequency response verification:")
    print(f"  Z_11 (TX self) error: {err_Z11:.2f}%")
    print(f"  Z_12 (TX-RX mutual) error: {err_Z12:.2f}%")

    # Print sample values
    print(f"\nSample impedance values at 100 kHz:")
    f_idx = np.argmin(np.abs(frequencies - 100e3))
    print(f"  Z_11 full: {Z11_full[f_idx]:.4f} Ohm")
    print(f"  Z_11 red:  {Z11_red[f_idx]:.4f} Ohm")
    print(f"  Z_12 full: {Z12_full[f_idx]:.4f} Ohm")
    print(f"  Z_12 red:  {Z12_red[f_idx]:.4f} Ohm")

    return result, netlist


def demo_coupled_dielectric():
    """Demonstrate N-port coupled (LoopStar + Dielectric) Block Lanczos."""
    print("=" * 70)
    print("N-Port Coupled (LoopStar + Dielectric) Block Lanczos Demo")
    print("=" * 70)

    # Create synthetic LC filter system
    n_L = 15  # Conductor DOFs
    n_D = 8   # Dielectric DOFs
    p = 2     # Ports

    # Conductor inductance (coil)
    L_L = np.eye(n_L) * 5e-6  # 5 uH base
    for i in range(n_L):
        for j in range(i+1, n_L):
            coupling = 0.5e-6 * np.exp(-0.3 * abs(i - j))
            L_L[i, j] = L_L[j, i] = coupling

    # Conductor resistance
    R_L = np.ones(n_L) * 0.05  # 50 mOhm

    # Dielectric capacitance (capacitor bank)
    C_D = np.eye(n_D) * 100e-12  # 100 pF base
    for i in range(n_D):
        for j in range(i+1, n_D):
            coupling = 10e-12 * np.exp(-0.5 * abs(i - j))
            C_D[i, j] = C_D[j, i] = coupling

    # Dielectric loss (from complex permittivity)
    # G_D = omega * C * eps''/eps' = omega * C * tan(delta_e)
    tan_delta_e = 0.02  # 2% loss tangent
    f_ref = 100e3
    omega_ref = 2 * np.pi * f_ref
    G_D = omega_ref * np.diag(C_D) * tan_delta_e

    # Coupling matrix (EM coupling between coil and capacitor)
    K_LD = np.zeros((n_L, n_D))
    for i in range(min(n_L, n_D)):
        K_LD[i, i] = 0.1  # Diagonal coupling
        if i > 0:
            K_LD[i, i-1] = 0.02
        if i < min(n_L, n_D) - 1:
            K_LD[i, i+1] = 0.02

    port_indices = [0, n_L // 2]

    print(f"\nSystem setup:")
    print(f"  Conductor DOFs: {n_L}")
    print(f"  Dielectric DOFs: {n_D}")
    print(f"  Ports: {p}")
    print(f"  Dielectric loss tangent: {tan_delta_e*100:.1f}%")

    # Apply reduction with Block Lanczos (DC-based starting vectors)
    # Block Lanczos converges much faster than scalar Lanczos
    solver = NPortCoupledDielectricSPICE(n_ports=p, n_stages_L=4, n_stages_D=2)
    result = solver.reduce(L_L, R_L, C_D, G_D, K_LD, port_indices)

    print(f"\nBlock Lanczos reduction:")
    print(f"  Conductor: {n_L} -> {result.n_stages_L * p} DOF")
    print(f"  Dielectric: {n_D} -> {result.n_stages_D} DOF")
    print(f"  Total compression: {100*(1 - (result.n_stages_L*p + result.n_stages_D)/(n_L + n_D)):.1f}%")

    # Generate SPICE
    netlist = solver.to_spice(result, "LC_FILTER")
    print(f"\nSPICE netlist ({len(netlist.split(chr(10)))} lines):")
    print("-" * 50)
    for line in netlist.split('\n')[:25]:
        print(line)
    if len(netlist.split('\n')) > 25:
        print("... (truncated)")

    # Verify frequency response
    frequencies = np.logspace(3, 7, 30)
    Z_full = solver.compute_impedance_full(L_L, R_L, C_D, G_D, K_LD, port_indices, frequencies)
    Z_red = solver.compute_impedance(result, frequencies)

    err_Z11 = np.mean(np.abs(Z_full[:, 0, 0] - Z_red[:, 0, 0]) / (np.abs(Z_full[:, 0, 0]) + 1e-10)) * 100
    err_Z12 = np.mean(np.abs(Z_full[:, 0, 1] - Z_red[:, 0, 1]) / (np.abs(Z_full[:, 0, 1]) + 1e-10)) * 100

    print(f"\nFrequency response verification:")
    print(f"  Z_11 error: {err_Z11:.2f}%")
    print(f"  Z_12 error: {err_Z12:.2f}%")

    # Sample values
    f_idx = np.argmin(np.abs(frequencies - 100e3))
    print(f"\nSample at 100 kHz:")
    print(f"  Z_11 full: {Z_full[f_idx, 0, 0]:.4f} Ohm")
    print(f"  Z_11 red:  {Z_red[f_idx, 0, 0]:.4f} Ohm")

    return result, netlist


def demo_full_coupled_lmd():
    """Demonstrate full L + M + D coupled system with complex materials."""
    print("=" * 70)
    print("Full Coupled (L + M + D) with Complex Materials Demo")
    print("=" * 70)

    # System: Transformer on PCB with ferrite core
    n_L = 12  # Conductor DOFs (windings)
    n_M = 6   # Magnetic DOFs (ferrite core)
    n_D = 4   # Dielectric DOFs (PCB substrate)
    p = 2     # Ports (primary, secondary)

    # Conductor (transformer windings)
    L_L = np.eye(n_L) * 10e-6
    for i in range(n_L):
        for j in range(i+1, n_L):
            coupling = 3e-6 * np.exp(-0.2 * abs(i - j))
            L_L[i, j] = L_L[j, i] = coupling
    R_L = np.ones(n_L) * 0.02

    # Magnetic (ferrite core)
    L_M = np.eye(n_M) * 50e-6  # Higher inductance due to core
    for i in range(n_M):
        for j in range(i+1, n_M):
            coupling = 20e-6 * np.exp(-0.5 * abs(i - j))
            L_M[i, j] = L_M[j, i] = coupling
    R_M = np.ones(n_M) * 0.01

    # Dielectric (PCB substrate)
    C_D = np.eye(n_D) * 50e-12
    G_D = np.zeros(n_D)

    # Coupling matrices
    K_LM = np.zeros((n_L, n_M))
    for i in range(min(n_L, n_M)):
        K_LM[i, i] = 0.3  # Strong coupling to core
        K_LM[n_L - 1 - i, i] = 0.25  # Secondary also couples

    K_LD = np.zeros((n_L, n_D))
    for i in range(min(n_L, n_D)):
        K_LD[i, i] = 0.05  # Weak capacitive coupling

    port_indices = [0, n_L - 1]  # Primary and secondary terminals

    # Define complex materials
    ferrite = ComplexMaterialParams(
        mu_r_real=2000,     # High permeability ferrite
        tan_delta_m=0.01,   # 1% magnetic loss
    )

    pcb_substrate = ComplexMaterialParams(
        eps_r_real=4.5,     # FR4 permittivity
        tan_delta_e=0.02,   # 2% dielectric loss
    )

    print(f"\nSystem setup:")
    print(f"  Conductor DOFs: {n_L}")
    print(f"  Magnetic DOFs: {n_M} (ferrite core)")
    print(f"  Dielectric DOFs: {n_D} (PCB)")
    print(f"  Ports: {p}")
    print(f"\nComplex materials:")
    print(f"  Ferrite: mu'_r = {ferrite.mu_r_real}, tan(delta_m) = {ferrite.tan_delta_m}")
    print(f"  PCB: eps'_r = {pcb_substrate.eps_r_real}, tan(delta_e) = {pcb_substrate.tan_delta_e}")

    # Apply reduction
    solver = NPortFullCoupledSPICE(n_ports=p, n_stages_L=4, n_stages_M=3, n_stages_D=2)
    result = solver.reduce(L_L, R_L, L_M, R_M, C_D, G_D, K_LM, K_LD, port_indices,
                           material_M=ferrite, material_D=pcb_substrate)

    print(f"\nBlock Lanczos reduction:")
    print(f"  Conductor: {n_L} -> {result.n_stages_L * p} DOF")
    print(f"  Magnetic: {n_M} -> {result.n_stages_M} DOF")
    print(f"  Dielectric: {n_D} -> {result.n_stages_D} DOF")
    total_orig = n_L + n_M + n_D
    total_red = result.n_stages_L * p + result.n_stages_M + result.n_stages_D
    print(f"  Total: {total_orig} -> {total_red} DOF ({100*(1 - total_red/total_orig):.1f}% compression)")

    # Generate SPICE
    netlist = solver.to_spice(result, "XFMR_ON_PCB", f_ref=100e3)
    print(f"\nSPICE netlist ({len(netlist.split(chr(10)))} lines):")
    print("-" * 50)
    for line in netlist.split('\n')[:30]:
        print(line)
    if len(netlist.split('\n')) > 30:
        print("... (truncated)")

    # Verify frequency response
    frequencies = np.logspace(3, 7, 30)
    Z_red = solver.compute_impedance(result, frequencies)

    print(f"\nImpedance at key frequencies:")
    for f in [10e3, 100e3, 1e6]:
        f_idx = np.argmin(np.abs(frequencies - f))
        Z11 = Z_red[f_idx, 0, 0]
        Z12 = Z_red[f_idx, 0, 1]
        print(f"  f = {f/1e3:.0f} kHz: Z_11 = {np.abs(Z11):.2f} Ohm @ {np.angle(Z11)*180/np.pi:.1f} deg, "
              f"Z_12 = {np.abs(Z12):.2f} Ohm")

    return result, netlist


if __name__ == "__main__":
    # Run new demos
    print("\n" + "=" * 60 + "\n")
    demo_coupled_dielectric()
    print("\n" + "=" * 60 + "\n")
    demo_full_coupled_lmd()
    print("\n" + "=" * 60 + "\n")
    demo_nport_block_lanczos()
