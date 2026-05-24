"""Krylov subspace methods: CG, BiCGSTAB, GMRES, COCG, COCR, IDR(s)."""

CATALOG = r"""
# Krylov subspace methods catalog

All Krylov methods build approximations in K_m = span{r0, Ar0, A^2 r0, ...}.
What differs: HOW the m-th basis is orthogonalized, and WHICH matrix
property is exploited.

| Method | Year | Matrix req | Short rec? | Lab use |
|--------|------|------------|------------|---------|
| **CG** (Hestenes-Stiefel) | 1952 | spd | yes (3-term) | ★ lab CGSolver |
| MINRES | 1975 | symmetric (indef ok) | yes | (NGSolve) |
| GMRES(m) | 1986 | general | no (m-term restart) | (NGSolve fallback) |
| BiCG | 1976 | general | yes | rarely (BiCGSTAB instead) |
| **BiCGSTAB** (vdVorst) | 1992 | general | yes | (NGSolve fallback) |
| QMR | 1991 | general | yes | rarely |
| **COCG** (vdVorst) | 1990 | complex symmetric | yes | legacy |
| **COCR** (Sogabe-Zhang) | 2007 | complex symmetric | yes | ★ lab COCRSolver |
| IDR(s) (Sonneveld-vGijzen) | 2008 | general | quasi short | future |

## "Short recurrence" matters

Short recurrence = constant memory per iteration, no restart cost.
Long recurrence (GMRES) = m-vector basis kept in memory; must restart.

For HCurl eddy current problems (lab core), the matrix is **complex
symmetric** so COCR (short rec, no restart) wins over restarted
GMRES.  The shipped `radia.sparsesolv_ngsolve.COCRSolver` is exactly
this.
"""


CG_HESTENES_STIEFEL = r"""
# Conjugate Gradient — Hestenes-Stiefel 1952

Reference: M. Hestenes, E. Stiefel, "Methods of Conjugate Gradients
for Solving Linear Systems", J. Res. NBS 49(6):409, 1952.
[LOCAL] 03_krylov/01_Hestenes_Stiefel_1952_CG_original.pdf

## Algorithm (simplified)

```
r0 = b - A x0;  p0 = r0
for k = 0, 1, 2, ...
    alpha_k = (r_k, r_k) / (p_k, A p_k)
    x_{k+1} = x_k + alpha_k * p_k
    r_{k+1} = r_k - alpha_k * A * p_k
    if ||r_{k+1}|| < tol: stop
    beta_k = (r_{k+1}, r_{k+1}) / (r_k, r_k)
    p_{k+1} = r_{k+1} + beta_k * p_k
```

3-term recurrence — only 4 vectors in memory regardless of iteration count.

## Convergence

For spd A with condition number κ = λ_max/λ_min:
```
||x_k - x*||_A / ||x_0 - x*||_A ≤ 2 * ((√κ-1)/(√κ+1))^k
```
→ iterations ~ √κ to reach tolerance.

## Preconditioned CG (PCG)

```
PCG solves M^{-1} A x = M^{-1} b
where M ≈ A (the preconditioner).
Effective κ becomes κ(M^{-1} A) instead of κ(A).
```

Typical M choices:
- M = diag(A) → Jacobi PCG (cheap, weak)
- M = IC factor → ICCG (Meijerink-vdVorst 1977)  ★ lab
- M = AMG V-cycle → AMG-PCG (Ruge-Stuben)
- M = AMS cycle → AMS-PCG for H(curl)  ★ lab
"""


GMRES = r"""
# GMRES — Saad-Schultz 1986

Reference: Y. Saad, M. Schultz, "GMRES: A Generalized Minimal Residual
Algorithm for Solving Nonsymmetric Linear Systems", SIAM SISC
7(3):856-869, 1986.  DOI: 10.1137/0907058.

## What it is

Builds an orthonormal basis V_m of K_m via Arnoldi process, then solves
the least-squares problem min ||β e_1 - H_m y|| for y.

**Strength**: works for ANY non-singular A; minimizes residual in
2-norm at every iteration.

**Weakness**: requires storing m basis vectors → O(m·N) memory.
Practical restart: GMRES(m), restart every m iterations (usually
m=20..50).  Restarting can stagnate on hard problems.

## Use in NGSolve

```python
from ngsolve.solvers import GMRes
GMRes(a.mat, f.vec, freedofs=fes.FreeDofs(), pre=c.mat,
      maxsteps=200, tol=1e-8, restart=30)
```

## When GMRES wins

- Non-symmetric problems with **moderate** non-symmetry (advection-diffusion).
- Need monotonic residual reduction (BiCGSTAB can oscillate).

## When GMRES loses

- Very long iteration count (memory blowup).
- Symmetric problems (CG / MINRES are 2-3x faster and constant memory).
- Severe non-symmetry / large complex part: try IDR(s) or BiCGSTAB(2).
"""


BICGSTAB = r"""
# BiCGSTAB — van der Vorst 1992

Reference: H. A. van der Vorst, "Bi-CGSTAB: A Fast and Smoothly
Converging Variant of Bi-CG", SIAM SISC 13(2):631-644, 1992.
DOI: 10.1137/0913035.

## What it is

Combines BiCG (Lanczos for non-symmetric) with a 1-step GMRES residual
smoothing.  Short recurrence (~7 vectors).  No restart needed.

## When to pick

- General **non-symmetric** large systems (N > 100k).
- Not too far from symmetric (otherwise BiCGSTAB(2) / IDR(s) are better).

## NGSolve

```python
from ngsolve.solvers import BVP, CGSolver
inv = ksp.GMRes  # or BiCGStab if exposed; for ngsolve 6.2.2603+:
from radia.sparsesolv_ngsolve import SparseSolvSolver  # wraps CG/COCR
```

## Convergence quirks

- **Bumpy residual**: BiCGSTAB residual can spike up before going down
  (unlike GMRES which is monotonic).  This is normal.
- **Breakdown**: if (r_tilde, r_k) ≈ 0, BiCGSTAB stalls.  Restart with
  different r_tilde or switch to IDR(s).
- **Round-off**: residual can decouple from true residual on very long
  runs; check ||b - A x|| every 50 iter.
"""


COCG_COCR = r"""
# Complex symmetric Krylov: COCG and COCR

For eddy current MQS systems, the discretized Maxwell operator in HCurl
is **complex symmetric** (A = A^T, complex entries):

```
A = ν · curl curl + jω σ · I    on H(curl)
```

CG cannot be used (not Hermitian), but a Lanczos-like 3-term recurrence
DOES exist for complex symmetric matrices.  Two families:

## COCG — Conjugate Orthogonal CG (van der Vorst-Melissen 1990)

[LOCAL] 03_krylov/COCG/02_COCG_vanderVorst.pdf
Reference: vanderVorst-Melissen, IEEE TMAG 1990.

Direct analogue of CG with the bilinear (·, ·) replaced by the
**unconjugated** inner product (no complex conjugate).  Same 3-term
recurrence, same low memory.

**Issue**: can break down when (p_k, A p_k) → 0 even though A is
non-singular.  Residual can stagnate or jump.

## COCR — Conjugate Orthogonal Conjugate Residual (Sogabe-Zhang 2007) ★

Reference: T. Sogabe, S.-L. Zhang, "A COCR method for solving complex
symmetric linear systems", J. Comp. Appl. Math. 199(2):297-303, 2007.
DOI: 10.1016/j.cam.2005.07.032.

★ This is the solver shipped as `radia.sparsesolv_ngsolve.COCRSolver`.

Algorithm (concise):
```
r0 = b - A x0;  p0 = r0;  s0 = A r0
for k = 0, 1, 2, ...
    alpha = (s_k, r_k)_uncon / (s_k, s_k)_uncon
    x_{k+1} = x_k + alpha * p_k
    r_{k+1} = r_k - alpha * s_k
    if ||r_{k+1}|| < tol: stop
    beta = (s_{k+1}, r_{k+1})_uncon / (s_k, r_k)_uncon       (1-step lookback)
    p_{k+1} = r_{k+1} + beta * p_k
    s_{k+1} = A r_{k+1} + beta * s_k
```

Where `(u, v)_uncon = sum(u_i * v_i)` (no conjugate).

## Why COCR > COCG for production

| Property | COCG | COCR |
|----------|------|------|
| Breakdown risk | high near (p, Ap)=0 | low (residual-based) |
| Residual smoothness | bumpy | smoother |
| Cost / iter | 1 matvec + 6 axpy | 1 matvec + 7 axpy |
| Lab CompactAMS compatibility | yes | yes (shipped) |
| Default in radia.sparsesolv_ngsolve | no | **yes** |

## Use

```python
import radia.sparsesolv_ngsolve as ssn
from radia.sparsesolv_ngsolve import COCRSolver, ComplexCompactAMSPreconditioner

prec = ComplexCompactAMSPreconditioner(a.mat, ...)
solver = COCRSolver(a.mat, prec, tol=1e-10, maxiter=500)
u_vec = solver.Solve(f.vec)
```

See `radia_ngsolve` MCP tool `sparsesolv('compact_ams')` for the full
code recipe.
"""


IDR_S = r"""
# IDR(s) — Sonneveld-van Gijzen 2009

Reference: P. Sonneveld, M. van Gijzen, "IDR(s): A Family of Simple
and Fast Algorithms for Solving Large Nonsymmetric Systems of Linear
Equations", SIAM SISC 31(2):1035-1062, 2009.
DOI: 10.1137/070685804.

## Idea

Generate residuals in a sequence of shrinking subspaces of decreasing
dimension N, N-s, N-2s, ...  At each shrink, project onto a
**fixed** shadow space of dimension s.

## Properties

- **Short recurrence** (memory ~ s+5 vectors, typically s=4 → ~9 vectors).
- Often **faster than BiCGSTAB** on hard non-symmetric problems.
- Has the **finite termination** property: at most N matvecs (vs GMRES
  which also terminates but at exponential memory cost).

## When IDR(s) wins

- Strongly non-symmetric, complex non-Hermitian.
- Convection-dominated convection-diffusion.
- Helmholtz with absorbing BC (frequency-domain wave).

## Status in lab

Not yet shipped in `radia.sparsesolv_ngsolve`.  For non-symmetric
problems use BiCGSTAB or GMRES from NGSolve.  IDR(s) is on the future
TODO list for cases where COCR breaks down (very low-frequency complex
systems where eddy current term dominates over ν·curl·curl).
"""


def get_krylov_knowledge(topic: str = "catalog") -> str:
    """Dispatch Krylov-method topics.

    Topics:
        catalog       - All Krylov methods at a glance (DEFAULT)
        cg            - Hestenes-Stiefel 1952 CG
        gmres         - Saad-Schultz 1986 GMRES
        bicgstab      - vdVorst 1992 BiCGSTAB
        cocg_cocr     - Complex symmetric: COCG (vdVorst 1990) and COCR (Sogabe-Zhang 2007) ★ lab core
        idr           - Sonneveld-vGijzen 2009 IDR(s)
        all           - Everything
    """
    topic = topic.lower().strip()
    if topic in ("catalog", "overview"):
        return CATALOG
    if topic in ("cg", "pcg", "hestenes"):
        return CG_HESTENES_STIEFEL
    if topic == "gmres":
        return GMRES
    if topic == "bicgstab":
        return BICGSTAB
    if topic in ("cocg_cocr", "cocg", "cocr", "complex_symmetric"):
        return COCG_COCR
    if topic in ("idr", "idr_s", "idrs"):
        return IDR_S
    if topic == "all":
        return "\n\n".join([CATALOG, CG_HESTENES_STIEFEL, GMRES, BICGSTAB,
                            COCG_COCR, IDR_S])
    return (f"Unknown topic '{topic}'. Available: catalog, cg, gmres, "
            "bicgstab, cocg_cocr, idr, all.")
