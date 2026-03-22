"""Type stubs for sparsesolv_ngsolve — SparseSolv iterative solvers for NGSolve.

Provides IC/SGS preconditioners, ICCG/SGSMRTR iterative solvers,
Compact AMS preconditioners for HCurl eddy-current problems,
and COCR/GMRES Krylov solvers.
"""

from ngsolve import BaseMatrix, BaseVector, BitArray

__all__ = [
    "SparseSolvResult",
    "ICPreconditionerD",
    "ICPreconditionerC",
    "SGSPreconditionerD",
    "SGSPreconditionerC",
    "SparseSolvSolverD",
    "SparseSolvSolverC",
    "ICPreconditioner",
    "SGSPreconditioner",
    "SparseSolvSolver",
    "COCRSolverD",
    "COCRSolverC",
    "COCRSolver",
    "GMRESSolverD",
    "GMRESSolverC",
    "GMRESSolver",
    "CompactAMSPreconditionerImpl",
    "ComplexCompactAMSPreconditionerImpl",
    "CompactAMSPreconditioner",
    "ComplexCompactAMSPreconditioner",
    "CompactAMGPreconditioner",
    "has_compact_ams",
]

# =============================================================================
# Result type
# =============================================================================

class SparseSolvResult:
    """Result of a SparseSolv iterative solve."""

    converged: bool
    """Whether the solver converged within tolerance."""
    iterations: int
    """Number of iterations performed."""
    final_residual: float
    """Final relative residual (or best residual if save_best_result enabled)."""
    residual_history: list[float]
    """Residual at each iteration (empty unless save_residual_history enabled)."""

# =============================================================================
# IC Preconditioner types
# =============================================================================

class ICPreconditionerD(BaseMatrix):
    """Incomplete Cholesky preconditioner for real (double) matrices."""

    shift: float
    """IC shift parameter."""
    use_abmc: bool
    """Enable ABMC ordering for parallel triangular solves."""
    abmc_block_size: int
    """Rows per block in ABMC aggregation."""
    abmc_num_colors: int
    """Target number of colors for ABMC coloring."""
    diagonal_scaling: bool
    """Diagonal scaling for improved conditioning."""

    def Update(self) -> None:
        """Recompute the IC factorization with current shift."""
        ...

class ICPreconditionerC(BaseMatrix):
    """Incomplete Cholesky preconditioner for complex matrices."""

    shift: float
    """IC shift parameter."""
    use_abmc: bool
    """Enable ABMC ordering for parallel triangular solves."""
    abmc_block_size: int
    """Rows per block in ABMC aggregation."""
    abmc_num_colors: int
    """Target number of colors for ABMC coloring."""
    diagonal_scaling: bool
    """Diagonal scaling for improved conditioning."""

    def Update(self) -> None:
        """Recompute the IC factorization with current shift."""
        ...

# =============================================================================
# SGS Preconditioner types
# =============================================================================

class SGSPreconditionerD(BaseMatrix):
    """Symmetric Gauss-Seidel preconditioner for real (double) matrices."""

    def Update(self) -> None:
        """Recompute the SGS factorization."""
        ...

class SGSPreconditionerC(BaseMatrix):
    """Symmetric Gauss-Seidel preconditioner for complex matrices."""

    def Update(self) -> None:
        """Recompute the SGS factorization."""
        ...

# =============================================================================
# Solver types
# =============================================================================

class SparseSolvSolverD(BaseMatrix):
    """Iterative solver for real (double) matrices."""

    method: str
    """Solver method: ``"ICCG"``, ``"SGSMRTR"``, ``"CG"``, or ``"COCR"``."""
    tol: float
    """Convergence tolerance (relative residual)."""
    maxiter: int
    """Maximum iterations."""
    shift: float
    """IC shift parameter."""
    auto_shift: bool
    """Auto-increase shift on IC breakdown (default: False)."""
    diagonal_scaling: bool
    """Diagonal scaling for improved conditioning (default: False)."""
    save_best_result: bool
    """Track and return the best iterate found (default: True)."""
    save_residual_history: bool
    """Record residual at each iteration (default: False)."""
    printrates: bool
    """Print convergence info to stdout (default: False)."""
    conjugate: bool
    """Use conjugated inner product for Hermitian systems (default: False)."""
    divergence_check: bool
    """Enable stagnation-based divergence detection (default: True)."""
    divergence_threshold: float
    """Threshold for divergence detection (default: 1.0)."""
    divergence_count: int
    """Iterations before declaring divergence (default: 50)."""
    use_abmc: bool
    """Enable ABMC ordering for parallel triangular solves (default: False)."""
    abmc_block_size: int
    """Rows per block in ABMC aggregation (default: 4)."""
    abmc_num_colors: int
    """Target number of colors for ABMC coloring (default: 4)."""
    abmc_reorder_spmv: bool
    """Reorder SpMV in ABMC space (default: False)."""
    abmc_use_rcm: bool
    """Use RCM ordering (default: False)."""

    @property
    def last_result(self) -> SparseSolvResult:
        """Result from the last Solve() call."""
        ...

    def Solve(self, rhs: BaseVector, sol: BaseVector) -> SparseSolvResult:
        """Solve Ax = b.

        Args:
            rhs: Right-hand side vector.
            sol: Solution vector (input: initial guess, output: solution).

        Returns:
            SparseSolvResult with convergence info.
        """
        ...

class SparseSolvSolverC(BaseMatrix):
    """Iterative solver for complex matrices."""

    method: str
    tol: float
    maxiter: int
    shift: float
    auto_shift: bool
    """Auto-increase shift on IC breakdown (default: False)."""
    diagonal_scaling: bool
    """Diagonal scaling for improved conditioning (default: False)."""
    save_best_result: bool
    save_residual_history: bool
    printrates: bool
    conjugate: bool
    divergence_check: bool
    divergence_threshold: float
    divergence_count: int
    use_abmc: bool
    abmc_block_size: int
    abmc_num_colors: int
    abmc_reorder_spmv: bool
    abmc_use_rcm: bool

    @property
    def last_result(self) -> SparseSolvResult: ...

    def Solve(self, rhs: BaseVector, sol: BaseVector) -> SparseSolvResult:
        """Solve Ax = b.

        Args:
            rhs: Right-hand side vector.
            sol: Solution vector (input: initial guess, output: solution).

        Returns:
            SparseSolvResult with convergence info.
        """
        ...

# =============================================================================
# Factory functions (auto-dispatch real/complex based on mat.IsComplex())
# =============================================================================

def ICPreconditioner(
    mat: BaseMatrix,
    freedofs: BitArray | None = None,
    shift: float = 1.05,
) -> ICPreconditionerD | ICPreconditionerC:
    """Incomplete Cholesky (IC) Preconditioner.

    Auto-dispatches to real/complex based on ``mat.IsComplex()``.
    Calls ``Update()`` automatically on construction.

    Args:
        mat: SPD sparse matrix (real or complex).
        freedofs: Free DOFs. Constrained DOFs treated as identity.
        shift: Shift parameter for stability (default: 1.05).
    """
    ...

def SGSPreconditioner(
    mat: BaseMatrix,
    freedofs: BitArray | None = None,
) -> SGSPreconditionerD | SGSPreconditionerC:
    """Symmetric Gauss-Seidel (SGS) Preconditioner.

    Auto-dispatches to real/complex based on ``mat.IsComplex()``.
    Calls ``Update()`` automatically on construction.

    Args:
        mat: SPD sparse matrix (real or complex).
        freedofs: Free DOFs. Constrained DOFs treated as identity.
    """
    ...

def SparseSolvSolver(
    mat: BaseMatrix,
    method: str = "ICCG",
    freedofs: BitArray | None = None,
    tol: float = 1e-10,
    maxiter: int = 1000,
    shift: float = 1.05,
    save_best_result: bool = True,
    save_residual_history: bool = False,
    printrates: bool = False,
    conjugate: bool = False,
    use_abmc: bool = False,
    abmc_block_size: int = 4,
    abmc_num_colors: int = 4,
    abmc_reorder_spmv: bool = False,
    abmc_use_rcm: bool = False,
) -> SparseSolvSolverD | SparseSolvSolverC:
    """Iterative solver (ICCG / SGSMRTR / CG / COCR).

    Can be used as a BaseMatrix inverse operator (``solver * rhs``) or
    via ``Solve()`` for detailed convergence results.

    Auto-dispatches to real/complex based on ``mat.IsComplex()``.

    Args:
        mat: SPD sparse matrix (real or complex).
        method: ``"ICCG"``, ``"SGSMRTR"``, ``"CG"``, or ``"COCR"``.
        freedofs: Free DOFs.
        tol: Convergence tolerance (default: 1e-10).
        maxiter: Maximum iterations (default: 1000).
        shift: IC shift parameter (default: 1.05).
        save_best_result: Track and return best iterate (default: True).
        save_residual_history: Record residual per iteration (default: False).
        printrates: Print convergence info to stdout (default: False).
        conjugate: Conjugated inner product for Hermitian systems (default: False).
        use_abmc: Enable ABMC parallel ordering (default: False).
        abmc_block_size: Rows per ABMC block (default: 4).
        abmc_num_colors: Target ABMC colors (default: 4).
        abmc_reorder_spmv: Reorder SpMV in ABMC space (default: False).
        abmc_use_rcm: Use RCM ordering (default: False).
    """
    ...

# =============================================================================
# COCR Solver (C++ native, NGSolve BaseMatrix interface)
# =============================================================================

class COCRSolverD(BaseMatrix):
    """COCR solver for real (double) matrices."""

    @property
    def iterations(self) -> int:
        """Number of iterations performed in last solve."""
        ...

class COCRSolverC(BaseMatrix):
    """COCR solver for complex matrices."""

    @property
    def iterations(self) -> int:
        """Number of iterations performed in last solve."""
        ...

def COCRSolver(
    mat: BaseMatrix,
    pre: BaseMatrix,
    freedofs: BitArray | None = None,
    maxiter: int = 500,
    tol: float = 1e-8,
    printrates: bool = False,
) -> COCRSolverD | COCRSolverC:
    """COCR (Conjugate Orthogonal Conjugate Residual) solver.

    For complex-symmetric systems (A^T = A, NOT Hermitian).
    Uses unconjugated inner products. Minimizes ||A r~||_2.

    Auto-dispatches to real/complex based on ``mat.IsComplex()``.

    Usage (same as NGSolve CGSolver):
        inv = COCRSolver(mat, pre, maxiter=500, tol=1e-8)
        gfu.vec.data = inv * rhs.vec

    For COCG, use ``CGSolver(mat, pre, conjugate=False)`` instead.

    Args:
        mat: System matrix (real or complex).
        pre: Preconditioner (BaseMatrix).
        maxiter: Maximum iterations (default: 500).
        tol: Relative convergence tolerance (default: 1e-8).
        printrates: Print convergence info (default: False).
    """
    ...

# =============================================================================
# GMRES Solver (C++ native, NGSolve BaseMatrix interface)
# =============================================================================

class GMRESSolverD(BaseMatrix):
    """GMRES solver for real (double) non-symmetric matrices."""

    @property
    def iterations(self) -> int:
        """Number of iterations performed in last solve."""
        ...

class GMRESSolverC(BaseMatrix):
    """GMRES solver for complex non-symmetric matrices."""

    @property
    def iterations(self) -> int:
        """Number of iterations performed in last solve."""
        ...

def GMRESSolver(
    mat: BaseMatrix,
    pre: BaseMatrix,
    freedofs: BitArray | None = None,
    maxiter: int = 500,
    tol: float = 1e-8,
    restart: int = 0,
    printrates: bool = False,
) -> GMRESSolverD | GMRESSolverC:
    """Left-preconditioned GMRES solver for non-symmetric linear systems.

    Optimal for AMS-preconditioned eddy current problems where COCR
    cannot be used (non-symmetric preconditioner).

    Auto-dispatches to real/complex based on ``mat.IsComplex()``.

    Args:
        mat: System matrix (real or complex).
        pre: Preconditioner (BaseMatrix).
        freedofs: Free DOFs mask.
        maxiter: Maximum iterations (default: 500).
        tol: Relative convergence tolerance (default: 1e-8).
        restart: GMRES restart (0 = full, default: 0).
        printrates: Print convergence info (default: False).
    """
    ...

# =============================================================================
# CompactAMG Preconditioner
# =============================================================================

def CompactAMGPreconditioner(
    mat: BaseMatrix,
    freedofs: BitArray | None = None,
    theta: float = 0.25,
    max_levels: int = 25,
    min_coarse: int = 50,
    num_smooth: int = 1,
    print_level: int = 0,
) -> BaseMatrix:
    """Compact AMG (Algebraic Multigrid) preconditioner for H1 Poisson-type systems.

    Header-only implementation, no external dependency.
    Uses classical AMG with PMIS coarsening and extended+i interpolation.

    Args:
        mat: Real SPD sparse matrix.
        freedofs: Free DOFs mask.
        theta: Strength threshold (default: 0.25).
        max_levels: Maximum AMG levels (default: 25).
        min_coarse: Minimum coarsest-level DOFs (default: 50).
        num_smooth: Smoother sweeps per level (default: 1).
        print_level: Verbosity (default: 0).
    """
    ...

# =============================================================================
# Compact AMS Preconditioner types (with Update() for nonlinear solvers)
# =============================================================================

class CompactAMSPreconditionerImpl(BaseMatrix):
    """Compact AMS preconditioner for real HCurl systems.

    Supports ``Update()`` for Newton iteration: geometry (G, Pi matrices)
    is preserved, only matrix-dependent parts are rebuilt.
    """

    def Update(self, new_mat: BaseMatrix | None = None) -> None:
        """Rebuild with current or new matrix values (geometry preserved).

        Args:
            new_mat: If provided, replaces the system matrix before rebuilding.
                     If None, rebuilds using the current matrix.
        """
        ...

def CompactAMSPreconditioner(
    mat: BaseMatrix,
    grad_mat: BaseMatrix,
    freedofs: BitArray | None = None,
    coord_x: list[float] = ...,
    coord_y: list[float] = ...,
    coord_z: list[float] = ...,
    cycle_type: int = 1,
    print_level: int = 0,
    subspace_solver: int = 0,
    num_smooth: int = 1,
) -> CompactAMSPreconditionerImpl:
    """Compact AMS (Auxiliary-space Maxwell Solver) Preconditioner.

    For real HCurl curl-curl + mass systems. No external dependency.
    Supports ``Update()`` for Newton iteration.

    Args:
        mat: Real HCurl system matrix.
        grad_mat: Discrete gradient G (HCurl -> H1).
        freedofs: Free DOFs mask.
        coord_x: Vertex x-coordinates (length = H1 DOFs).
        coord_y: Vertex y-coordinates.
        coord_z: Vertex z-coordinates.
        cycle_type: AMS cycle type (1=01210, 7=0201020, default=1).
        print_level: Verbosity (default: 0).
        subspace_solver: 0=CompactAMG (default), 1=SparseCholesky.
        num_smooth: Smoother sweeps (default: 1).
    """
    ...

class ComplexCompactAMSPreconditionerImpl(BaseMatrix):
    """Complex Compact AMS preconditioner with fused Re/Im operations.

    Supports ``Update()`` for Newton iteration: geometry is preserved,
    only matrix-dependent parts are rebuilt.
    """

    def Update(self, new_a_real: BaseMatrix | None = None) -> None:
        """Rebuild with current or new real auxiliary matrix (geometry preserved).

        Args:
            new_a_real: If provided, replaces the real auxiliary matrix before rebuilding.
                        If None, rebuilds using the current matrix.
        """
        ...

def ComplexCompactAMSPreconditioner(
    a_real_mat: BaseMatrix,
    grad_mat: BaseMatrix,
    freedofs: BitArray | None = None,
    coord_x: list[float] = ...,
    coord_y: list[float] = ...,
    coord_z: list[float] = ...,
    ndof_complex: int = 0,
    cycle_type: int = 1,
    print_level: int = 0,
    correction_weight: float = 1.0,
    subspace_solver: int = 0,
    num_smooth: int = 1,
) -> ComplexCompactAMSPreconditionerImpl:
    """Complex Compact AMS preconditioner with fused Re/Im operations.

    For complex eddy-current systems ``A = K + jw*sigma*M``.
    Uses CompactAMG (header-only, no external dependency).
    Supports ``Update()`` for Newton iteration.

    Use with ``COCRSolver`` (complex symmetric) or ``GMRESSolver``.

    Args:
        a_real_mat: Real SPD auxiliary matrix (K + eps*M + |omega|*sigma*M).
        grad_mat: Discrete gradient G (HCurl -> H1).
        freedofs: Free DOFs mask.
        coord_x: Vertex x-coordinates.
        coord_y: Vertex y-coordinates.
        coord_z: Vertex z-coordinates.
        ndof_complex: Complex DOF count (0 = auto-derive from matrix).
        cycle_type: AMS cycle type (default: 1).
        print_level: Verbosity (default: 0).
        correction_weight: Correction weight (default: 1.0).
        subspace_solver: Subspace solver type (default: 0 = CompactAMG).
        num_smooth: Number of smoothing steps (default: 1).
    """
    ...

def has_compact_ams() -> bool:
    """Returns True if Compact AMG/AMS support is available."""
    ...
