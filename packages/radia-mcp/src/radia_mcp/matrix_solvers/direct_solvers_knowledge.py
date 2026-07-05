"""Direct sparse solvers: LU, PARDISO, MUMPS, SuperLU."""

OVERVIEW = r"""
# Direct sparse solvers

Direct solvers factor the matrix A = LU (or LDL^T for symmetric) and
back-substitute.  Cost: O(N^{1.5}) memory for 2D, O(N^2) memory for 3D
elliptic problems — limits practical N to ~500k DOF in 3D.

| Solver | Distributed? | Lic | Best for |
|--------|--------------|-----|----------|
| PARDISO (Intel MKL) | shared mem (OpenMP) | Intel MKL | Default in NGSolve; fast on Intel |
| MUMPS | MPI distributed | CeCILL-C | Large clusters, indefinite |
| SuperLU_DIST | MPI | BSD | Non-symmetric, sparse direct |
| UMFPACK | serial | LGPL | Small/medium, used by scipy |
| KLU | serial | LGPL | Circuit simulation (very sparse) |
| LU built-in (radia method=0) | serial | (lab) | N < 500 |

## When direct beats iterative

1. **Many RHS sweeps with same A** (frequency sweep, port S-parameters):
   factor once, back-substitute N_rhs times. Cost amortizes.
2. **Ill-conditioned matrices** where Krylov stalls (κ > 10^10).
3. **Small N** (< 10k): iterative setup overhead dominates.
4. **Saddle-point / indefinite**: many Krylov + AMG fail without
   sophisticated block preconditioning.

## When direct loses

1. **Single solve, large N (>500k 3D)**: memory blows up.
2. **Memory-limited GPU**: no GPU direct solver matches MKL/HYPRE iter.
3. **Time-stepping with changing A**: refactor cost per step kills.
"""


PARDISO = r"""
# PARDISO (Schenk-Gärtner 2003)

Reference: Schenk-Gärtner, "Solving unsymmetric sparse systems of
linear equations with PARDISO", FGCS 2003.  DOI: 10.1016/S0167-739X(03)00188-2.

## What it is

Multifrontal direct solver with **METIS** reordering and **OpenMP**
parallelism.  Bundled in Intel MKL as `mkl_pardiso` — ships with
NGSolve binary distributions automatically.

## How to invoke from NGSolve

```python
from ngsolve import *
mesh = Mesh(...); fes = H1(mesh, order=2)
a = BilinearForm(fes); a += ...; a.Assemble()
f = LinearForm(fes); f += ...; f.Assemble()

inv = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso")   # or "sparsecholesky"
u = GridFunction(fes); u.vec.data = inv * f.vec
```

Other choices: `inverse="sparsecholesky"` (built-in symmetric direct),
`inverse="umfpack"` (built-in non-symmetric), `inverse="mumps"`
(if compiled with MUMPS).

## Performance characteristics

- 8-core Xeon, 3D Poisson 1M DOF tet mesh: ~30 s setup, ~5 s back-sub.
- Memory: ~12x sparse matrix size for elliptic 3D (cubic fill-in).
- Threading: scales well to ~8 threads, plateaus beyond.

## Pitfalls

- Complex symmetric: PARDISO supports `mtype=6` but NGSolve `Inverse`
  defaults to LU on complex.  Use `inverse="pardiso"` explicitly.
- Singular matrices (Periodic with ungauged Omega-reduction): PARDISO
  may return garbage without error.  Always check with one CG iteration
  on the same matrix.
- MKL license: included with NGSolve PyPI wheel via `mkl>=2024.2.0`
  dependency.  No separate license needed.
"""


MUMPS = r"""
# MUMPS (Amestoy-Duff-L'Excellent 2000)

Reference: Amestoy et al., "Multifrontal Parallel Distributed Symmetric
and Unsymmetric Solvers", CMAME 184:501, 2000.

## What it is

**Distributed-memory** multifrontal direct solver via MPI.  Handles:
- Real / complex
- Symmetric / unsymmetric
- Definite / indefinite (Bunch-Kaufman pivoting for symmetric indefinite)

## When to use over PARDISO

| Need | Pick |
|------|------|
| Shared-memory single node | PARDISO (simpler, no MPI) |
| Distributed cluster (>1 node) | MUMPS |
| Symmetric **indefinite** (saddle-point) | MUMPS (PARDISO mtype handling is brittle) |
| Out-of-core (huge matrices) | MUMPS |

## NGSolve invocation

```python
inv = a.mat.Inverse(fes.FreeDofs(), inverse="mumps")
```

Requires NGSolve compiled with `-DUSE_MUMPS=ON`.  The PyPI wheel does
NOT include MUMPS — build from source if needed.

## Pitfalls

- MUMPS uses Fortran ordering — make sure CSR conversion is correct.
- MPI initialization must come BEFORE NGSolve import.
- Output verbosity defaults to chatty; set `-1, -1, -1, -1` for silence.
"""


LU_RADIA = r"""
# Built-in LU (Radia method=0)

For Radia HDiv-VIM interaction matrices (N < 500 elements), the built-in
LAPACK LU is appropriate:

```python
import radia as rad
rad.SolverConfig(...)        # not needed for LU
rad.Solve(container, 1e-4, 1000, 0)    # method=0 -> LU
```

- Uses LAPACK `dgetrf` / `dgetri` from MKL.
- Guaranteed convergence (it's not iterative).
- Cost: O(N^3) factor, O(N^2) back-sub per RHS.
- Threading: MKL handles internally via `MKL_NUM_THREADS`.

For N > 500, switch to BiCGSTAB (method=1) or HACApK (method=2).  See
`HMATRIX_EVALUATION.md` in the Radia docs.
"""


def get_direct_solvers_knowledge(topic: str = "overview") -> str:
    """Dispatch direct-solver topics.

    Topics:
        overview     - When direct beats iterative (DEFAULT)
        pardiso      - Intel MKL PARDISO (NGSolve default)
        mumps        - Distributed MUMPS
        lu_radia     - Built-in LU for Radia HDiv-VIM
        all          - Everything
    """
    topic = topic.lower().strip()
    if topic in ("overview", "general", "summary"):
        return OVERVIEW
    if topic == "pardiso":
        return PARDISO
    if topic == "mumps":
        return MUMPS
    if topic in ("lu", "lu_radia", "radia_lu"):
        return LU_RADIA
    if topic == "all":
        return "\n\n".join([OVERVIEW, PARDISO, MUMPS, LU_RADIA])
    return (f"Unknown topic '{topic}'. Available: overview, pardiso, mumps, "
            "lu_radia, all.")
