"""
PRIMA-based Model Order Reduction for Loop-Star-Magnetic System

This module implements PRIMA (Passive Reduced-order Interconnect Macromodeling
Algorithm) with Lanczos process for SPICE-compatible circuit extraction.

Key features:
1. PRIMA Lanczos with re-orthogonalization (higher accuracy than plain Lanczos)
2. Tridiagonal structure yields RL ladder network (no Arnoldi, direct circuit)
3. ACA+ compression for material blocks (magnetic, dielectric, conductor)
4. Schur complement for port impedance extraction
5. PyKAN integration for complex material properties (mu", eps")

The goal is to reduce:
    [Z_LL    Z_LS    Z_LM  ] [I_L]   [V_L  ]
    [Z_SL    Z_SS    Z_SM  ] [Q_S] = [0    ]
    [Z_ML    Z_MS    Z_MM  ] [M  ]   [H_ext]

To a sparse equivalent circuit with KAN-based material models.

Note: This implementation uses PRIMA (not CLN/Cauer), avoiding patent issues
while providing passive, stable reduced-order models.

References:
- A. Odabasioglu, M. Celik, L.T. Pileggi, "PRIMA: Passive Reduced-order
  Interconnect Macromodeling Algorithm," IEEE TCAD, 1998.
- B. Gustavsen, A. Semlyen, "Rational approximation of frequency domain
  responses by Vector Fitting," IEEE TPWRD, 1999.

Author: Radia Development Team
Date: 2026-01-16
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


class ACAPlus:
    """
    ACA+ (Adaptive Cross Approximation with pivoting).

    Computes low-rank approximation: A ≈ U @ V.T

    Algorithm:
    1. Select pivot row/column based on maximum residual
    2. Update U and V
    3. Repeat until ||residual|| < tol * ||A||
    """

    def __init__(self, tol: float = 1e-4, max_rank: int = None):
        self.tol = tol
        self.max_rank = max_rank

    def compress(self, A: np.ndarray) -> LowRankBlock:
        """
        ACA+ compression of dense matrix.

        Parameters:
            A: Dense matrix [m x n]

        Returns:
            LowRankBlock approximation
        """
        m, n = A.shape
        max_rank = min(m, n) if self.max_rank is None else min(self.max_rank, m, n)

        # Working copy for residual tracking
        R = A.copy()

        U_list = []
        V_list = []

        norm_A = np.linalg.norm(A, 'fro')
        if norm_A < 1e-15:
            # Zero matrix
            return LowRankBlock(
                U=np.zeros((m, 1)),
                S=np.array([0.0]),
                V=np.zeros((n, 1)),
                rank=0
            )

        for k in range(max_rank):
            # Find pivot (maximum absolute value in residual)
            idx = np.unravel_index(np.argmax(np.abs(R)), R.shape)
            i_pivot, j_pivot = idx

            pivot_val = R[i_pivot, j_pivot]
            if np.abs(pivot_val) < self.tol * norm_A / max(m, n):
                break

            # Extract row and column
            u = R[:, j_pivot].copy()
            v = R[i_pivot, :].copy() / pivot_val

            U_list.append(u)
            V_list.append(v)

            # Update residual
            R -= np.outer(u, v)

            # Check convergence
            norm_R = np.linalg.norm(R, 'fro')
            if norm_R < self.tol * norm_A:
                break

        if len(U_list) == 0:
            return LowRankBlock(
                U=np.zeros((m, 1)),
                S=np.array([0.0]),
                V=np.zeros((n, 1)),
                rank=0
            )

        # Convert to U @ S @ V.T form via SVD of small matrix
        U = np.column_stack(U_list)
        V = np.column_stack(V_list)

        # Compact SVD
        C = U.T @ A @ V  # Small matrix [r x r]
        Uc, Sc, Vct = linalg.svd(C, full_matrices=False)

        # Truncate small singular values
        r = np.sum(Sc > self.tol * Sc[0])
        r = max(1, r)

        return LowRankBlock(
            U=U @ Uc[:, :r],
            S=Sc[:r],
            V=V @ Vct[:r, :].T,
            rank=r
        )


class HierarchicalReducer:
    """
    Hierarchical model order reduction combining Lanczos and ACA+.

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

    Step 3: ACA+ on coupling blocks
        Z_LM' ≈ U_LM @ S_LM @ V_LM.T
        Z_LS' ≈ U_LS @ S_LS @ V_LS.T

    Step 4: Schur complement
        Z_port = Z_LL' - Z_LM' @ Z_MM^{-1} @ Z_ML' - Z_LS' @ Z_SS^{-1} @ Z_SL'
    """

    def __init__(self, lanczos_tol: float = 1e-6, aca_tol: float = 1e-4,
                 max_lanczos_rank: int = None, max_aca_rank: int = None):
        self.lanczos = LanczosReducer(tol=lanczos_tol, max_rank=max_lanczos_rank)
        self.aca = ACAPlus(tol=aca_tol, max_rank=max_aca_rank)

        # Stored results
        self._lanczos_result = None
        self._reduced_blocks = {}
        self._aca_blocks = {}

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

    def apply_aca_compression(self) -> dict:
        """
        Step 3: ACA+ compression of transformed coupling blocks.

        Returns:
            Dictionary of LowRankBlock approximations
        """
        self._aca_blocks = {}

        for name, block in self._reduced_blocks.items():
            if block is not None and block.size > 0:
                self._aca_blocks[name] = self.aca.compress(block)

        return self._aca_blocks

    def compute_schur_complement(self, s: complex,
                                  Z_SS: np.ndarray = None,
                                  Z_MM: np.ndarray = None,
                                  use_aca: bool = False) -> np.ndarray:
        """
        Step 4: Schur complement to eliminate Star and Magnetic DOFs.

        PEEC-MMM System (physically correct formulation):
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
            use_aca: Use ACA+ compressed blocks (experimental)

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

        # Use dense reduced blocks (more accurate)
        if not use_aca:
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

        else:
            # Use ACA+ compressed (experimental, less accurate)
            if 'Z_LM' in self._aca_blocks and Z_MM is not None:
                lr_LM = self._aca_blocks['Z_LM']
                try:
                    Z_MM_inv = np.linalg.inv(Z_MM)
                    temp = (lr_LM.S[:, None] * lr_LM.V.T) @ Z_MM_inv @ (lr_LM.V * lr_LM.S)
                    Z_schur += s * lr_LM.U @ temp @ lr_LM.U.T
                except np.linalg.LinAlgError:
                    pass

            if 'Z_LS' in self._aca_blocks and Z_SS is not None:
                lr_LS = self._aca_blocks['Z_LS']
                lr_SL = self._aca_blocks['Z_SL']
                try:
                    Z_SS_inv = np.linalg.inv(Z_SS)
                    temp = (lr_LS.S[:, None] * lr_LS.V.T) @ Z_SS_inv @ (lr_SL.U * lr_SL.S)
                    Z_schur -= s**3 * lr_LS.U @ temp @ lr_SL.V.T
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


def symmetrize_mmm_block(Z_MM: np.ndarray) -> np.ndarray:
    """
    Symmetrize MMM (Magnetic Moment Method) impedance matrix.

    For PEEC-MMM coupling, the MMM block must be symmetric to ensure
    proper energy conservation and passivity of the coupled system.

    The symmetrization follows the approach:
        Z_MM_sym = 0.5 * (Z_MM + Z_MM.T)

    Parameters:
        Z_MM: MMM impedance matrix (may be slightly asymmetric)

    Returns:
        Symmetrized Z_MM matrix
    """
    return 0.5 * (Z_MM + Z_MM.T)


def build_peec_mmm_system(L: np.ndarray, R: np.ndarray, L_LM: np.ndarray,
                           Z_MM: np.ndarray, symmetrize: bool = True) -> dict:
    """
    Build PEEC-MMM coupled system matrices.

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
        Z_MM: MMM reluctance/impedance matrix [n_M x n_M]
        symmetrize: Whether to symmetrize Z_MM (default: True)

    Returns:
        dict with keys:
            'L': PEEC inductance
            'R': PEEC resistance
            'K': Coupling matrix (= L_LM)
            'Z_MM': Symmetrized MMM impedance
            'L_coupling': K @ Z_MM^{-1} @ K.T (coupling contribution)
            'L_eff': L + L_coupling (effective inductance, INCREASES with mu > 1)
    """
    if symmetrize:
        Z_MM = symmetrize_mmm_block(Z_MM)

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
    """Demonstrate hierarchical reduction on a physically realistic PEEC-MMM system."""
    print("=" * 60)
    print("PRIMA-style Hierarchical Model Order Reduction Demo")
    print("=" * 60)
    print("(PEEC-MMM coupled system with physically correct formulation)")

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

    # Generate MMM impedance matrix with proper scaling
    Z_MM = generate_magnetic_impedance(n_M, L[0, 0], off_diag_ratio=0.1)

    # Build PEEC-MMM system with physically correct formulation
    system = build_peec_mmm_system(L, R, K, Z_MM, symmetrize=True)

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
    aca_tol: float              # ACA tolerance for this group
    n_lanczos: int = 0          # Lanczos order (0 = no reduction)
    material_type: str = 'linear'  # 'linear', 'nonlinear', 'kan'
    kan_model: Optional['KANMaterialInterface'] = None  # KAN model if applicable


@dataclass
class SPICEExtractionConfig:
    """Configuration for SPICE circuit extraction."""
    # PRIMA Lanczos settings
    n_lanczos_loop: int = 20        # Lanczos order for Loop (conductor) DOFs
    n_lanczos_star: int = 10        # Lanczos order for Star (capacitive) DOFs

    # Element group ACA tolerances
    aca_tol_magnetic: float = 1e-4    # ACA tolerance for magnetic elements
    aca_tol_dielectric: float = 1e-4  # ACA tolerance for dielectric elements
    aca_tol_conductor: float = 1e-4   # ACA tolerance for conductor (shield) elements

    # Port settings
    port_indices: List[int] = field(default_factory=list)
    port_names: List[str] = field(default_factory=list)

    # Frequency range for fitting
    f_min: float = 1.0              # Minimum frequency [Hz]
    f_max: float = 1e6              # Maximum frequency [Hz]
    n_freq: int = 50                # Number of frequency points


class PRIMASchurExtractor:
    """
    PRIMA-based Schur complement extractor with per-group ACA control.

    This class implements a complete workflow for SPICE netlist extraction:

    1. PRIMA Lanczos Reduction:
       - Apply block Lanczos to Loop and Star DOFs
       - Specify Lanczos order (ladder stages) for each block
       - Result: Tridiagonal (PRIMA ladder) structure

    2. Schur Complement:
       - Eliminate internal DOFs, keep port DOFs
       - Port impedance: Z_port = Z_PP - Z_PI @ Z_II^{-1} @ Z_IP

    3. ACA+ Compression (per element group):
       - Magnetic: Z_MM block with aca_tol_magnetic
       - Dielectric: Z_DD block with aca_tol_dielectric
       - Conductor/Shield: Z_CC block with aca_tol_conductor

    4. SPICE Netlist Generation:
       - Each Lanczos stage -> RL ladder element
       - ACA modes -> coupled subcircuits

    Usage:
    ------
    config = SPICEExtractionConfig(
        n_lanczos_loop=20,
        n_lanczos_star=10,
        aca_tol_magnetic=1e-3,
        aca_tol_dielectric=1e-4,
        aca_tol_conductor=1e-5,
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
            Extraction configuration with Lanczos orders and ACA tolerances
        """
        self.config = config

        # ACA extractors for each element group
        self.aca_magnetic = ACAPlus(tol=config.aca_tol_magnetic)
        self.aca_dielectric = ACAPlus(tol=config.aca_tol_dielectric)
        self.aca_conductor = ACAPlus(tol=config.aca_tol_conductor)

        # Lanczos reducer
        self.lanczos = LanczosReducer(tol=1e-10, max_rank=max(
            config.n_lanczos_loop, config.n_lanczos_star, 100))

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
            'aca_magnetic': ACA result for magnetic block
            'aca_dielectric': ACA result for dielectric block
            'aca_conductor': ACA result for conductor block
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
        print("PRIMA-Schur SPICE Extraction")
        print("=" * 70)
        print(f"\nSystem dimensions:")
        print(f"  Loop (conductor):    {n_L} DOFs -> Lanczos order {self.config.n_lanczos_loop}")
        print(f"  Star (capacitive):   {n_S} DOFs -> Lanczos order {self.config.n_lanczos_star}")
        print(f"  Magnetic:            {n_M} DOFs -> ACA tol {self.config.aca_tol_magnetic:.1e}")
        print(f"  Dielectric:          {n_D} DOFs -> ACA tol {self.config.aca_tol_dielectric:.1e}")
        print(f"  Conductor (shield):  {n_C} DOFs -> ACA tol {self.config.aca_tol_conductor:.1e}")
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
        # Step 2: ACA compression for each material group
        # ================================================================
        print("\n[Step 2] ACA compression per element group...")

        aca_results = {}

        # Magnetic ACA
        if n_M > 0:
            print(f"  Magnetic: ", end="")
            aca_M = self._aca_compress_block(
                Z_MM, self.aca_magnetic, "magnetic", kan_models.get('magnetic'))
            aca_results['magnetic'] = aca_M
            print(f"rank {aca_M['rank']}/{n_M} (compression {100*(1-aca_M['rank']/n_M):.1f}%)")

        # Dielectric ACA
        if n_D > 0:
            print(f"  Dielectric: ", end="")
            aca_D = self._aca_compress_block(
                Z_DD, self.aca_dielectric, "dielectric", kan_models.get('dielectric'))
            aca_results['dielectric'] = aca_D
            print(f"rank {aca_D['rank']}/{n_D} (compression {100*(1-aca_D['rank']/n_D):.1f}%)")

        # Conductor/Shield ACA
        if n_C > 0:
            print(f"  Conductor: ", end="")
            aca_C = self._aca_compress_block(
                Z_CC, self.aca_conductor, "conductor", kan_models.get('conductor'))
            aca_results['conductor'] = aca_C
            print(f"rank {aca_C['rank']}/{n_C} (compression {100*(1-aca_C['rank']/n_C):.1f}%)")

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
                aca_results, port_proj, kan_models
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
            aca_results, K_LM_red, K_LD_red, K_LC_red, K_LS_red,
            port_indices, frequencies, Z_port_func
        )

        print(f"  Generated {circuit['n_elements']} circuit elements")
        print(f"  Netlist: {len(netlist.split(chr(10)))} lines")

        return {
            'netlist': netlist,
            'circuit': circuit,
            'lanczos_loop': lanczos_L,
            'lanczos_star': lanczos_S,
            'aca_results': aca_results,
            'Z_port_func': Z_port_func,
            'config': self.config,
            'dimensions': {
                'n_L': n_L, 'n_S': n_S, 'n_M': n_M, 'n_D': n_D, 'n_C': n_C,
                'n_lanczos_L': n_lanczos_L, 'n_lanczos_S': n_lanczos_S,
                'n_ports': n_ports
            }
        }

    def _aca_compress_block(self,
                            Z: np.ndarray,
                            aca: 'ACAPlus',
                            group_name: str,
                            kan_model: Optional['KANMaterialInterface'] = None) -> dict:
        """Apply ACA compression to a material block."""
        if Z is None or Z.size == 0:
            return {'rank': 0, 'U': None, 'V': None, 'kan': kan_model}

        n = Z.shape[0]
        U, V, rank = aca.decompose(Z)

        return {
            'rank': rank,
            'U': U,
            'V': V,
            'Z_original': Z,
            'kan': kan_model,
            'group': group_name
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
                                aca_results: dict,
                                port_proj: np.ndarray,
                                kan_models: dict) -> np.ndarray:
        """Compute port impedance via Schur complement."""
        omega = np.abs(s.imag) if np.abs(s.imag) > 1e-10 else 1.0

        # Build reduced system matrix
        n_total = n_L + n_S + n_M + n_D + n_C
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

        # Magnetic block
        if n_M > 0:
            idx_M = idx
            aca_M = aca_results.get('magnetic', {})
            kan_mu = kan_models.get('magnetic') if kan_models else None
            if kan_mu:
                mu_eff = kan_mu.evaluate_frequency([omega/(2*np.pi)])[0]
            else:
                mu_eff = 1.0
            if aca_M.get('U') is not None:
                # Use low-rank: Z_MM ≈ U @ V.T
                Z_MM_approx = aca_M['U'] @ aca_M['V'].T * mu_eff
            else:
                Z_MM_approx = aca_M.get('Z_original', np.eye(n_M)) * mu_eff
            Z[idx:idx+n_M, idx:idx+n_M] = Z_MM_approx
            idx += n_M
        else:
            idx_M = idx

        # Dielectric block
        if n_D > 0:
            idx_D = idx
            aca_D = aca_results.get('dielectric', {})
            kan_eps = kan_models.get('dielectric') if kan_models else None
            if kan_eps:
                eps_eff = kan_eps.evaluate_frequency([omega/(2*np.pi)])[0]
            else:
                eps_eff = 1.0
            if aca_D.get('U') is not None:
                Z_DD_approx = aca_D['U'] @ aca_D['V'].T * eps_eff
            else:
                Z_DD_approx = aca_D.get('Z_original', np.eye(n_D)) * eps_eff
            Z[idx:idx+n_D, idx:idx+n_D] = Z_DD_approx
            idx += n_D
        else:
            idx_D = idx

        # Conductor/Shield block
        if n_C > 0:
            idx_C = idx
            aca_C = aca_results.get('conductor', {})
            kan_sigma = kan_models.get('conductor') if kan_models else None
            if kan_sigma:
                sigma_eff = kan_sigma.evaluate_frequency([omega/(2*np.pi)])[0]
            else:
                sigma_eff = 1.0
            if aca_C.get('U') is not None:
                Z_CC_approx = aca_C['U'] @ aca_C['V'].T * sigma_eff
            else:
                Z_CC_approx = aca_C.get('Z_original', np.eye(n_C)) * sigma_eff
            Z[idx:idx+n_C, idx:idx+n_C] = Z_CC_approx
            idx += n_C
        else:
            idx_C = idx

        # Coupling blocks
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
                                aca_results: dict,
                                K_LM: np.ndarray,
                                K_LD: np.ndarray,
                                K_LC: np.ndarray,
                                K_LS: np.ndarray,
                                port_indices: List[int],
                                frequencies: np.ndarray,
                                Z_port_func: Callable) -> Tuple[str, dict]:
        """Generate SPICE netlist from reduced system."""
        lines = []
        lines.append("* PRIMA-Schur Extracted SPICE Netlist")
        lines.append("* Generated by Radia PEEC + Lanczos MOR")
        lines.append(f"* Lanczos order (Loop): {self.config.n_lanczos_loop}")
        lines.append(f"* Lanczos order (Star): {self.config.n_lanczos_star}")
        lines.append(f"* ACA tol (Magnetic):   {self.config.aca_tol_magnetic:.1e}")
        lines.append(f"* ACA tol (Dielectric): {self.config.aca_tol_dielectric:.1e}")
        lines.append(f"* ACA tol (Conductor):  {self.config.aca_tol_conductor:.1e}")
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
        # Material subcircuits (from ACA)
        # ================================================================
        for group_name, aca_data in aca_results.items():
            if aca_data.get('rank', 0) > 0:
                lines.append(f"* === {group_name.capitalize()} (ACA rank {aca_data['rank']}) ===")

                U = aca_data.get('U')
                V = aca_data.get('V')
                rank = aca_data['rank']

                # Each ACA mode becomes a subcircuit
                for k in range(rank):
                    n_in = next_node()
                    n_out = next_node()

                    # Mode impedance (simplified: R + sL approximation)
                    # In practice, fit from frequency response
                    sigma_k = np.sum(U[:, k] * V[:, k])
                    R_mode = float(np.abs(sigma_k)) * 1e-3
                    L_mode = float(np.abs(sigma_k)) * 1e-9

                    lines.append(f"* {group_name} mode {k+1}")
                    lines.append(f"R_{group_name[0].upper()}{k+1} {n_in} {n_in}a {R_mode:.6e}")
                    lines.append(f"L_{group_name[0].upper()}{k+1} {n_in}a {n_out} {L_mode:.6e}")

                    circuit_elements.append({
                        'type': 'mode',
                        'group': group_name,
                        'mode': k,
                        'R': R_mode,
                        'L': L_mode
                    })

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

    # Configuration
    config = SPICEExtractionConfig(
        n_lanczos_loop=10,          # Reduce 30 -> 10 stages
        n_lanczos_star=5,           # Reduce 10 -> 5 stages
        aca_tol_magnetic=1e-3,      # Coarser for magnetic
        aca_tol_dielectric=1e-4,    # Medium for dielectric
        aca_tol_conductor=1e-5,     # Finer for conductor (shield)
        port_indices=[0, 15],       # Two ports
        port_names=['IN', 'OUT'],
        f_min=100,                  # 100 Hz
        f_max=10e6,                 # 10 MHz
        n_freq=50
    )

    print(f"\nExtraction configuration:")
    print(f"  Lanczos (Loop):  {n_L} -> {config.n_lanczos_loop} stages")
    print(f"  Lanczos (Star):  {n_S} -> {config.n_lanczos_star} stages")
    print(f"  ACA tol (Mag):   {config.aca_tol_magnetic:.1e}")
    print(f"  ACA tol (Diel):  {config.aca_tol_dielectric:.1e}")
    print(f"  ACA tol (Cond):  {config.aca_tol_conductor:.1e}")

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

    for name, aca in result['aca_results'].items():
        n_orig = {'magnetic': n_M, 'dielectric': n_D, 'conductor': n_C}[name]
        print(f"  {name.capitalize()}: rank {aca['rank']}/{n_orig}")

    print(f"\nSPICE netlist ({len(result['netlist'].split(chr(10)))} lines):")
    print("-" * 50)
    for line in result['netlist'].split('\n')[:40]:
        print(line)
    if len(result['netlist'].split('\n')) > 40:
        print("... (truncated)")

    return result


if __name__ == "__main__":
    demo_hierarchical_reduction()
    print("\n" + "=" * 60 + "\n")
    demo_kan_material()
    print("\n" + "=" * 60 + "\n")
    demo_aca_circuit_extraction()
    print("\n" + "=" * 60 + "\n")
    demo_prima_schur_extraction()
