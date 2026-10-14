"""Matrix-solver overview: lab stack, decision tree, history."""

from radia_mcp.matrix_solvers.direct_solvers_knowledge import SOLVER_CAPACITY

LAB_STACK = r"""
# Radia matrix-solver stack

★ Production solver library = **`radia.sparsesolv_ngsolve`**
  (built from `src/ext/sparsesolv/`, shipped inside the radia wheel)

| Component | Origin paper | Year | Where in code |
|-----------|--------------|------|----------------|
| `CGSolver` (real spd) | Hestenes-Stiefel | 1952 | `cg_solver.hpp` |
| `COCRSolver` (complex sym) ★ | Sogabe-Zhang | 2007 | `cocr_solver.hpp` |
| `ICPreconditioner` | Meijerink-vanderVorst | 1977 | `ic_preconditioner.hpp` |
| `CompactAMG` ★ | Ruge-Stuben 1987 / Henson-Yang 2002 | — | `compact_amg.hpp` |
| `CompactAMSPreconditioner` ★ | Hiptmair-Xu | 2007 | `compact_ams.hpp` |
| `ComplexCompactAMS` ★ | Hiptmair-Xu + Re/Im split | (lab) | `complex_compact_ams.hpp` |
| ABMC ordering | parallel triangular solve | (lab custom) | `abmc_ordering.hpp` |
| RCM ordering | Cuthill-McKee | 1969 | `rcm_ordering.hpp` |

The starred methods are available application building blocks. Choose by matrix
structure and verified residual, not by an assumed universal size threshold.
For direct FE solves, Radia explicitly uses SparseCholesky; PARDISO is not selected.

## Application priorities

| Goal | Choice | Why |
|------|--------|-----|
| Real HCurl systems | CG + CompactAMS | Use the validated application setup and check the original residual |
| Complex symmetric eddy current | COCR + ComplexCompactAMS | Check the complex system, not only the real auxiliary problem |
| Higher-order HCurl | BDDC + AMS coarse solver | Validated paths include order 2 and 3; do not infer higher-order or mixed-space support |
| Direct symmetric FE reference | SparseCholesky | Explicit backend, boundary lifting and true-residual validation |

## What `radia.sparsesolv_ngsolve` does NOT solve

- **Non-symmetric real** systems need a separately validated nonsymmetric
  method. Do not silently apply SparseCholesky or switch to another backend.
- **Helmholtz** / full Maxwell at >GHz: out of scope (lab is Laplace-kernel
  MQS/Darwin only — see Radia CLAUDE.md).
- **Indefinite** symmetric (Stokes-like): use MINRES + Schur-complement
  preconditioner (not in lab; defer to NGSolve).
"""


DECISION_TREE = r"""
# Decision tree: which solver+preconditioner combination

```
1. Matrix structure
   ├── Real, symmetric positive definite (HCurl mass, magnetostatic Phi)
   │   ├── HCurl             → validated CG + CompactAMS path
   │   ├── Other spaces      → validated space-appropriate preconditioner
   │   └── Direct FE choice  → SparseCholesky + true-residual check
   │
   ├── Real, symmetric indefinite (Stokes-like, saddle-point)
   │   └── MINRES + block preconditioner (defer to NGSolve)
   │
   ├── Real, non-symmetric (advection-diffusion)
   │   ├── Mild non-symmetry → BiCGSTAB + ILU
   │   └── Severe + restart  → GMRES(m) + AMG / SOR
   │
   ├── Complex symmetric (eddy current)
   │   ├── HCurl, order 1     → validated COCR + ComplexCompactAMS path
   │   ├── HCurl, order 2/3   → validated BDDC + AMS coarse path
   │   └── Direct FE choice   → SparseCholesky + true-residual check
   │
   └── Complex non-symmetric  → GMRES(m) + ILU(0) (rarely needed for MQS)

2. Special problem features
   ├── Multiply connected (loops in conductor)
   │   └── Tree-cotree gauging or A-V formulation (Biro-Preis 2000)
   ├── Air region σ=0 with conductor σ>0
   │   └── Use the application's validated nullspace/auxiliary setup
   ├── Frequency-stable from DC to MHz
   │   └── Two-Step Maxwell + tree-cotree (Ostrowski-Hiptmair 2021)
   └── Open boundary
       └── Kelvin transformation (NOT solver-related; geometry-level)
```

## Failure and comparison contract

Do not silently fall back to PARDISO or another backend. Preserve the original
operator, right-hand side, boundary conditions and acceptance tolerance. Check
the true free-row residual and report nonconvergence or factorization failure.
Benchmark AMS/BDDC+AMS and SparseCholesky on the same problem before making a
speed claim. Availability of HYPRE or PETSc must be checked explicitly; it is not
an automatic replacement for the supported Radia path.
"""


HISTORY = r"""
# Genealogy: how the lab's stack came to be

```
1952  Hestenes-Stiefel CG (NBS J Res 49:409)
   │  — original CG paper, conjugate directions
   │
1977  Meijerink-vanderVorst IC (Math Comp 31:148)
   │  — Incomplete Cholesky preconditioning for symmetric M-matrices
   │
1978  Kershaw ICCG (J Comp Phys 26:43)
   │  — ICCG popularization for plasma physics
1978  Gustafsson MICCG (BIT 18:142)
   │  — modified IC: better cond. number for elliptic problems
   │
1986  Saad-Schultz GMRES (SIAM SISC 7:856)
1987  Ruge-Stuben classical AMG (SIAM Frontiers vol 3 ch 4)
   │  — algebraic multigrid based on strength-of-connection
1989  Biro-Preis 3D eddy current A-V (IEEE IntMag)
1990  vanderVorst COCG (IEEE TMAG) — complex symmetric
1992  vanderVorst BiCGSTAB (SIAM SISC 13:631) — non-symmetric Krylov
1998  Hiptmair Multigrid for H(curl) (SIAM SINUM 36:204)
   │  — first multigrid for edge elements
2002  Henson-Yang BoomerAMG (Appl Num Math 41:155)
   │  — parallel AMG; basis for HYPRE
2003  Schenk-Gartner PARDISO (FGCS) — direct sparse defacto
2007  Hiptmair-Xu Nodal Auxiliary Space (SIAM SINUM 45:2483) ★
   │  — H(curl)/H(div) preconditioning via auxiliary nodal space
   │  — basis of NGSolve `Pi`-correction approach
2007  Sogabe-Zhang COCR (J Comp Appl Math 199:297) ★
   │  — modern complex symmetric Krylov, alternative to COCG
2009  Sonneveld-vGijzen IDR(s) (SIAM SISC 31:1035)
   │  — short-recurrence non-symmetric Krylov
2021  Ostrowski-Hiptmair Two-Step Maxwell (SIAM SISC 43)
   │  — frequency-stable formulation with tree-cotree gauge
```

★ = the two pillars of `radia.sparsesolv_ngsolve` (Hiptmair-Xu 2007 +
Sogabe-Zhang 2007).

## Cross-references

- `radia_ngsolve` MCP: `sparsesolv` tool — concrete code usage of the
  lab stack (CompactAMS, COCR, examples).
- `matrix_solvers.preconditioners` (this subpackage) — theory of the
  preconditioners shipped in sparsesolv.
- `motor_hollaus_eddy` MCP — application of CompactAMS to laminated cores.
"""


DECISION_TREE += SOLVER_CAPACITY

def get_overview_knowledge(topic: str = "lab_stack") -> str:
    """Dispatch overview topics.

    Topics:
        lab_stack       - Production solver components + which to pick
        decision_tree   - Solver+preconditioner decision tree (DEFAULT)
        history         - Genealogy from 1952 CG to Hiptmair-Xu 2007
        all             - Everything
    """
    topic = topic.lower().strip()
    if topic in ("lab_stack", "stack", "production"):
        return LAB_STACK
    if topic in ("decision_tree", "decision", "choose", "tree"):
        return DECISION_TREE
    if topic in ("capacity", "memory", "diagnostics", "bddc"):
        return SOLVER_CAPACITY
    if topic in ("history", "genealogy", "lineage"):
        return HISTORY
    if topic == "all":
        return "\n\n".join([LAB_STACK, DECISION_TREE, HISTORY])
    return (f"Unknown topic '{topic}'. Available: lab_stack, decision_tree, "
            "history, all.")
