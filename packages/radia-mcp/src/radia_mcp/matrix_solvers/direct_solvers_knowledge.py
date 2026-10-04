"""Direct-solver guidance for Radia's supported FE and legacy dense paths."""

SOLVER_CAPACITY = r"""
# Solver capacity and measured diagnosis

Direct solvers can be fast when the factor fits. Capacity and speed are separate:
element count or total runtime alone does not establish a factorization bottleneck.
There is no universal DOF/element cutoff or guaranteed memory estimate.

Sparse factors acquire fill-in depending on graph, ordering, FE space and order.
Matrix NNZ is not factor-entry storage. Real versus complex entries, indices,
workspace and temporary allocations all matter. On Windows, record system commit
headroom (commit limit minus committed bytes), not only free physical RAM;
process peak memory and system commit are different quantities. Neither a small
matrix nor available RAM guarantees that factorization will fit.

Report the actual chain: outer solver -> preconditioner -> coarse solver.
An iterative outer solve with BDDC can still use a direct SparseCholesky coarse
factorization. Explicit BDDC + AMS coarse selection changes that stage; do not
infer it from the word BDDC or claim that all auxiliary levels are factor-free.
For validated HCurl applications select AMS or BDDC+AMS explicitly within their
space/order/periodicity restrictions. Maxwell AMS is not an H1 prescription.
CG needs a positive-definite compatible system, not arbitrary complex matrices.
Native ICCG with conjugate=True supports Hermitian positive-definite systems,
checks Hermitian structure and rejects encountered non-positive curvature;
complex-symmetric eddy current is a different contract (COCR in validated HCurl
paths). IC shifts act on the preconditioner, not permission to alter the operator.
Iterative methods are not guaranteed faster, convergent or memory-safe.

## Diagnostic report: measured, estimated, unavailable

Label each quantity measured (including method/units), estimated (including
assumptions), or unavailable. Record space, order, free DOFs versus total DOFs,
real/complex arithmetic, matrix NNZ, actual solver chain, backend/version and
threads. Report factor/setup/solve times separately only if actually instrumented;
include factor entries, peak memory, Windows commit headroom, iterations and
true relative residual when available. An aggregate solve time is not a measured
factorization time. Missing telemetry is unavailable, never zero. An extrapolated
factor count is estimated, never measured; this guidance adds no memory predictor.

For IH calc_fem_kelvin JSON, linear_solver_requested and linear_solver distinguish
request from selection; bddc_ams_coarse_cycles describes that configured coarse
route. ndof is total DOFs, ne is element count, linear_krylov_iterations records outer
iterations, and linear_true_relative_residual / linear_true_residual_limit report the
linear gate. t_solve_s and t_total_s are aggregates, not isolated factor timings.
Unreported free DOFs, factor entries, peak memory or commit headroom remain
unavailable unless independently measured. Preserve the existing result schema.

The shared true linear residual gate remains 1e-6; keep stricter application,
nonlinear and physical acceptance criteria unchanged. No alternative direct backend and no silent
fallback after allocation failure or nonconvergence. Residual checks validate
returned solutions, not allocation/index safety before factorization.

## SparseCholesky evidence boundary

See docs/solver/SPARSECHOLESKY_LIMITS.md in the source repository. The 6.2.2607
Windows long counter is a source-level risk, not an established crash cause.
The reproduction was memory-confounded; no universal threshold was established.
A production pre-factorization symbolic-storage guard is not yet implemented.
This guidance does not implement or validate such a guard.
"""

OVERVIEW = r"""
# Radia direct-solver policy

For supported real symmetric and complex symmetric FE systems, explicitly select
`inverse="sparsecholesky"`. Do not select PARDISO or silently switch backends.
For HCurl applications, prefer the validated AMS or BDDC+AMS iterative path;
SparseCholesky is the direct reference and an explicit user choice, not an
assumption that direct factorization is always faster.

Factor reuse helps multiple right-hand sides only while the matrix is unchanged.
A frequency sweep generally changes K + i omega M, so it does not in general
reuse one factorization. Measure setup, solve, memory and true residual for the
actual problem; there is no universal DOF cutoff for direct versus iterative.

Keep prescribed boundary values and check the original free-row residual.
For constrained mixed systems, use the application's verified exact reduction
and multiplier recovery. Do not regularize the physical operator or change the
acceptance tolerance merely to make a factorization succeed. Unsupported matrix
structure, singularity, memory exhaustion or a failed residual check must be
reported as failures.

MKL remains a dependency of Radia's dense BLAS/LAPACK and FFT kernels. Removing
PARDISO selection does not make the entire Radia package MKL-independent.
"""

SPARSECHOLESKY = r"""
# SparseCholesky: Radia's FE direct backend

Use the built-in NGSolve backend explicitly. The caller owns TaskManager and
assembly. In this example a, f, fes and gfu already exist; gfu holds prescribed
boundary data before solving the free-DoF correction.

```python
from math import isfinite
from ngsolve import Norm, Projector

free = Projector(fes.FreeDofs(), True)
inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
rhs = f.vec.CreateVector()
rhs.data = f.vec - a.mat * gfu.vec
gfu.vec.data += inv * rhs
r = f.vec.CreateVector()
r.data = free * (f.vec - a.mat * gfu.vec)
scale = Norm(free * rhs)
rnorm = Norm(r)
relative_residual = rnorm / scale if scale else (0.0 if rnorm == 0 else float("inf"))
from radia._residual_gate import RELATIVE_LIMIT
if not isfinite(relative_residual) or relative_residual > RELATIVE_LIMIT:
    raise RuntimeError(f"direct solve failed: residual={relative_residual}")
```

This example is for a supported symmetric system. It does not authorize applying
SparseCholesky to arbitrary nonsymmetric or unhandled saddle-point matrices.
A failed solve must not fall back to PARDISO, alter a gauge/penalty, or accept a
nonfinite result. Applications can impose a stricter residual contract.

MATLAB `radia.ngsolve.Matrix.inverse()` explicitly selects this same native
backend. No Python process is needed by that MEX path.
"""

PARDISO = r"""
# Retired Radia FE solver selection: PARDISO

This legacy documentation topic remains discoverable to explain migration.
Radia does not select PARDISO for FE direct solves. Use SparseCholesky explicitly
for supported symmetric systems, or the validated AMS / BDDC+AMS application
path. Do not install another backend or silently fall back after a failure.

NGSolve's `ngsolve-openblas` dependency and Radia's `mkl>=2026,<2027` dependency
serve different native components. MKL remains for Radia dense BLAS/LAPACK and FFT;
its presence is not permission to choose PARDISO.
"""

MUMPS = r"""
# MUMPS: separate external integration

This topic is background, not an automatic Radia fallback. Radia's supported
symmetric FE direct path explicitly uses SparseCholesky. A workflow requiring a
different matrix structure or external solver needs its own implementation and
validation; it must not silently change the current problem or backend.

MUMPS is a distributed-memory multifrontal solver with real and complex,
symmetric and unsymmetric variants. It is relevant to separately configured MPI
workflows and indefinite systems. NGSolve exposes it when built with
`-DUSE_MUMPS=ON`; do not assume it is present in an installed wheel.

In such a separately validated integration, the NGSolve selector is
`a.mat.Inverse(fes.FreeDofs(), inverse="mumps")`. Initialize the MPI environment
before using that integration, retain the matrix's actual symmetry contract,
and verify the original residual. This selector is not a Radia FE fallback.
"""

LU_RADIA = r"""
# Legacy Radia dense LU (`method=0`)

The retained C++ relaxation ABI accepts dense LU:

```python
import radia as rad
rad.Solve(legacy_object, 1e-4, 1000, 0)
```

This compatibility path uses LAPACK through MKL. Factorization is not a guarantee
of a valid answer: singularity, conditioning and the true residual still matter.
This is not the HDiv-VIM solver. New soft-iron models must be mesh-backed and
use `radia.vim.Solve` or `radia.vim.HDivSolver.Solve`, whose `linear_solver`,
`preconditioner`, `gram_eps`, `leaf`, and `eta` options are named independently
of the legacy method integer. Do not recommend retired `rad.Solve` methods 1 or 2.
"""


OVERVIEW += SOLVER_CAPACITY
SPARSECHOLESKY += SOLVER_CAPACITY

def get_direct_solvers_knowledge(topic: str = "overview") -> str:
    """Return current direct guidance, or explain a retired topic's migration."""
    topic = topic.lower().strip()
    if topic in ("overview", "general", "summary"):
        return OVERVIEW
    if topic in ("sparsecholesky", "sparse_cholesky", "direct"):
        return SPARSECHOLESKY
    if topic in ("capacity", "memory", "diagnostics", "bddc"):
        return SOLVER_CAPACITY
    if topic == "pardiso":
        return PARDISO
    if topic == "mumps":
        return MUMPS
    if topic in ("lu", "lu_radia", "radia_lu"):
        return LU_RADIA
    if topic == "all":
        return "\n\n".join([OVERVIEW, SPARSECHOLESKY, PARDISO, MUMPS, LU_RADIA])
    return (f"Unknown topic '{topic}'. Available: overview, sparsecholesky, "
            "capacity (memory, diagnostics, bddc), pardiso (migration), mumps, lu_radia, all.")
