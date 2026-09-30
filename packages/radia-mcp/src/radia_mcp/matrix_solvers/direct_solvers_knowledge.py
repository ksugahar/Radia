"""Direct-solver guidance for Radia's supported FE and legacy dense paths."""

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


def get_direct_solvers_knowledge(topic: str = "overview") -> str:
    """Return current direct guidance, or explain a retired topic's migration."""
    topic = topic.lower().strip()
    if topic in ("overview", "general", "summary"):
        return OVERVIEW
    if topic in ("sparsecholesky", "sparse_cholesky", "direct"):
        return SPARSECHOLESKY
    if topic == "pardiso":
        return PARDISO
    if topic == "mumps":
        return MUMPS
    if topic in ("lu", "lu_radia", "radia_lu"):
        return LU_RADIA
    if topic == "all":
        return "\n\n".join([OVERVIEW, SPARSECHOLESKY, PARDISO, MUMPS, LU_RADIA])
    return (f"Unknown topic '{topic}'. Available: overview, sparsecholesky, "
            "pardiso (migration), mumps, lu_radia, all.")
