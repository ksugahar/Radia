"""Type stubs for sparsesolv_ngsolve — SparseSolv iterative solvers for NGSolve.

Provides IC preconditioner, ICCG iterative solver,
Compact AMS preconditioners for HCurl eddy-current problems,
and COCR/GMRES Krylov solvers.
"""

from ngsolve import BaseMatrix, BaseVector, BitArray, FESpace

__all__ = [
    "SparseSolvResult",
    "ICPreconditionerD",
    "ICPreconditionerC",
    "SparseSolvSolverD",
    "SparseSolvSolverC",
    "ICPreconditioner",
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
    """Whether the recursive relative residual fell below tol."""
    iterations: int
    """Number of iterations performed."""
    best_iteration: int
    """Iteration of the returned iterate (0 = initial guess)."""
    final_residual: float
    """Recursive relative residual of the returned iterate (scaled system when
    diagonal_scaling is on)."""
    true_residual: float
    """||b - A x|| / ||b|| of the returned x on the original (free-DOF) system."""
    actual_shift: float
    """IC shift used (0 when no IC factor was applied)."""
    residual_history: list[float]
    """[initial, iteration 1, ...] (empty unless save_residual_history enabled)."""

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
# Solver types
# =============================================================================

class SparseSolvSolverD(BaseMatrix):
    """Iterative solver for real (double) matrices."""

    method: str
    """Solver method: ``"ICCG"``, ``"CG"``, or ``"COCR"``."""
    tol: float
    """Stop when the recursive relative residual < tol (default: 1e-8)."""
    maxiter: int
    """Iteration limit; 0 means 2*n (default: 0)."""
    shift: float
    """IC shift >= 1; start value of the auto-shift search (default: 1.0)."""
    auto_shift: bool
    """Raise the shift by 0.01 while a pivot has Re(d) < 1e-6|a_ii| (limit 5, then error) (default: True)."""
    diagonal_scaling: bool
    """Solve the scaled system S A S (default: True)."""
    save_best_result: bool
    """Return the best iterate, initial guess included (default: True)."""
    save_residual_history: bool
    """Record residual at each iteration (default: False)."""
    printrates: bool
    """Print convergence info to stdout (default: False)."""
    conjugate: bool
    """Conjugated inner product for Hermitian systems; CG only (default: False)."""
    divergence_check: bool
    """Enable the stagnation stop (default: True)."""
    divergence_threshold: float
    """A residual below best*threshold resets the counter (default: 10.0)."""
    divergence_count: int
    """Stop when the counter exceeds this value (default: 10)."""
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
    """See SparseSolvSolverD.auto_shift (default: True)."""
    diagonal_scaling: bool
    """See SparseSolvSolverD.diagonal_scaling (default: True)."""
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

def SparseSolvSolver(
    mat: BaseMatrix,
    method: str = "ICCG",
    freedofs: BitArray | None = None,
    tol: float = 1e-8,
    maxiter: int = 0,
    shift: float = 1.0,
    save_best_result: bool = True,
    save_residual_history: bool = False,
    printrates: bool = False,
    conjugate: bool = False,
    use_abmc: bool = False,
    abmc_block_size: int = 4,
    abmc_num_colors: int = 4,
    abmc_reorder_spmv: bool = False,
    abmc_use_rcm: bool = False,
    auto_shift: bool = True,
    diagonal_scaling: bool = True,
    divergence_check: bool = True,
    divergence_threshold: float = 10.0,
    divergence_count: int = 10,
) -> SparseSolvSolverD | SparseSolvSolverC:
    """Iterative solver (ICCG / CG / COCR).

    Can be used as a BaseMatrix inverse operator (``solver * rhs``) or
    via ``Solve()`` for detailed convergence results.

    Auto-dispatches to real/complex based on ``mat.IsComplex()``.

    Args:
        mat: SPD sparse matrix (real or complex), full (non-symmetric) storage.
        method: ``"ICCG"``, ``"CG"``, or ``"COCR"``.
        freedofs: Free DOFs.
        tol: Stop when the recursive relative residual < tol (default: 1e-8).
        maxiter: Iteration limit, 0 means 2*n (default: 0).
        shift: IC shift >= 1, auto-shift start value (default: 1.0).
        save_best_result: Return best iterate incl. initial guess (default: True).
        save_residual_history: Record residual per iteration (default: False).
        printrates: Print convergence info to stdout (default: False).
        conjugate: Conjugated inner product for CG only (default: False).
            ICCG and COCR reject True.
        use_abmc: Enable ABMC parallel ordering (default: False).
        abmc_block_size: Rows per ABMC block (default: 4).
        abmc_num_colors: Target ABMC colors (default: 4).
        abmc_reorder_spmv: Reorder SpMV in ABMC space (default: False).
        abmc_use_rcm: Use RCM ordering (default: False).
        auto_shift: Auto IC shift search (default: True).
        diagonal_scaling: Solve the scaled system S A S (default: True).
        divergence_check: Stagnation stop (default: True).
        divergence_threshold: Counter reset factor (default: 10.0).
        divergence_count: Stop when the counter exceeds it (default: 10).
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

    @property
    def beta_zero(self) -> bool: ...

    @property
    def reuse_hierarchy(self) -> bool: ...

    @property
    def mixed_precision(self) -> bool: ...

    @property
    def hierarchy_refreshes(self) -> int: ...

    @property
    def in_place_updates(self) -> int:
        """Updates whose Galerkin products ran numerically on the previous patterns."""
        ...

    @property
    def setup_workers(self) -> int:
        """Workers observed in strength construction; zero if no coarsening."""
        ...

    def Update(self, new_mat: BaseMatrix | None = None) -> None:
        """Rebuild with current or new matrix values (geometry preserved).

        Args:
            new_mat: If provided, replaces the system matrix before rebuilding.
                     If None, rebuilds using the current matrix.
        """
        ...

def TaskManagerActive() -> bool:
    """True inside an ngsolve.TaskManager context (where AMS setup refuses to run)."""
    ...


def AMSCoarseStats() -> dict:
    """Counters of the most recent BDDC wirebasket AMS (``coarsetype="sparsesolv_ams"``).

    Importing this module registers ``"sparsesolv_ams"`` with NGSolve's
    preconditioner classes: ``Preconditioner(a, "bddc", coarsetype="sparsesolv_ams",
    coarseflags={...})`` on an HCurl form replaces the direct wirebasket inverse
    by Compact AMS on the lowest-order edge block (built inside Assemble, which
    may run in TaskManager). coarseflags: ``cycles`` (k stationary AMS steps on
    the wirebasket system, default 1), ``lean_coarse`` (default 1: the
    auxiliary AMGs stop coarsening where it stalls and solve a coarsest level
    of at most 1024 rows with a dense inverse, accepted only if it reproduces a
    test vector; 0 keeps sparse Cholesky), ``cycle_type``, ``num_smooth``,
    ``print_level``, ``eps`` (relative diagonal shift of the AMS surrogate),
    ``beta_zero`` and ``mixed_precision`` (real systems only). A complex
    wirebasket matrix S uses the real surrogate Re S + Im S. Keys: ``builds``,
    ``n_edges``, ``n_extra`` (non-edge wirebasket dofs, smoothed only),
    ``cycles``, ``n_free``, ``n_vertices``, ``complex``, ``extract_s``,
    ``setup_s``, ``applies``, ``apply_s`` and, for a complex system while the
    preconditioner is alive, ``cycle_stage_s``, ``cycles_run`` and
    ``amg_levels`` (per auxiliary AMG: rows, nonzeros, seconds per level and
    the coarsest solve, 1 sparse / 2 dense).
    """
    ...


def LowestOrderCurlSystem(fes: FESpace, coefficient: object = None) -> dict:
    """Element data of a lowest-order HCurl space on straight tetrahedra, in one pass.

    Keys: ``dofs`` (ne, 6) int32 ascending, ``curl`` (ne, 3, 6) element-constant
    basis curls, ``volume`` (ne,), ``matrix`` (element-graph SparseMatrix; with
    ``coefficient`` (ne,) it holds sum_e c_e vol_e curl_e^T curl_e, exactly
    symmetric, else zeros) and ``positions`` (ne, 36) int32 value indices.
    """
    ...


def ClosedCoilCurrentPhi(mesh: object, materials: list[int], current_A: float, origin: object,
                         normal: object, radius: float, inverse: str = "sparsecholesky") -> dict:
    """A-phi DC current of a closed conductor with one thick cut (conductor-only, native).

    Used by radia.meshed_current.solve_closed_coil_current_phi. ``inverse`` is
    "iccg" (IC(0)-CG on the scaled system to 1e-12, fails loudly) or an NGSolve inverse type. Returns
    ``elements``, ``density`` (n, 3), ``vertices``, ``phi``, the cut statistics,
    ``relative_weak_divergence``, ``cut_face_flux_A`` and ``timing``.
    """
    ...


class NativePCG:
    """Preconditioned CG stopped on the true relative residual over free dofs.

    ``Solve(b, x, tolerance, maxiter)`` starts from x = 0 and returns
    ``(iterations, true_relative_residual, converged)``; one product and one
    preconditioner application per iteration, the true residual computed only to
    confirm convergence; raises when p.Ap <= 0 or r.z <= 0.
    """

    def __init__(self, mat: BaseMatrix, pre: BaseMatrix, freedofs: BitArray | None) -> None: ...
    def Solve(self, b: BaseVector, x: BaseVector, tolerance: float, maxiter: int) -> tuple: ...


class LowestOrderCurlResidual:
    """Element flux, material state and residual of a lowest-order HCurl problem.

    ``Evaluate(x, source_mean, grid, nu, dhdb, load, residual)`` returns
    ``(b, magnitude, nu, q)`` on the nonlinear elements and writes
    ``residual`` = sum_e vol_e C_e^T (nu_e c_e + [iron] (nu_e - nu0) Bs_e) - load.
    """

    def __init__(self, dofs: object, curl: object, volume: object, ndof: int, iron: object,
                 nu0: float) -> None: ...
    def Evaluate(self, x: object, source_mean: object, grid: object, nu: object, dhdb: object,
                 load: object, residual: object) -> tuple: ...


class LowestOrderCurlJacobian:
    """Newton Jacobian refresh on the elements given (LowestOrderCurlSystem data).

    ``Refresh(nu, q, b)`` rewrites every row the elements touch as the saved
    constant part plus sum_e vol_e C_e^T (nu_e I + q_e b_e b_e^T) C_e, gathered
    per row in ascending element order (exactly symmetric).
    """

    def __init__(self, matrix: BaseMatrix, dofs: object, curl: object, volume: object,
                 positions: object) -> None: ...
    def Refresh(self, nu: object, q: object, b: object) -> None: ...
    @property
    def elements(self) -> int: ...
    @property
    def rows(self) -> int: ...


def LowestOrderGradient(fes: FESpace) -> BaseMatrix:
    """Discrete gradient H1(order 1) -> lowest-order HCurl from the edge table.

    Equals ``fes.CreateGradient()[0]`` for ``HCurl(order=1, nograds=True)``
    (or order 0); raises for any other dof layout.
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
    beta_zero: bool = False,
    reuse_hierarchy: bool = False,
    mixed_precision: bool = False,
) -> CompactAMSPreconditionerImpl:
    """Compact AMS (Auxiliary-space Maxwell Solver) Preconditioner.

    For real HCurl curl-curl + mass systems, or compatible pure curl-curl
    systems with ``beta_zero=True``. No external dependency.
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
        beta_zero: Skip gradient correction and its hierarchy (default: False).
        reuse_hierarchy: Update() keeps the first AMG coarsening and interpolation
            and refreshes only the Galerkin coarse matrices (default: False).
        mixed_precision: residual products inside the cycle read float32 value
            copies (vectors and sums double); requires beta_zero (default: False).
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
