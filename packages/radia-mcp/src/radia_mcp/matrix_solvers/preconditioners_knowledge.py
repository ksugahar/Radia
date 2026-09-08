"""Preconditioner catalog: classical, AMG, AMS (Hiptmair-Xu)."""

CATALOG = r"""
# Preconditioner catalog

A preconditioner M approximates A^{-1}.  Solving M^{-1} A x = M^{-1} b
converges in O(√κ(M^{-1} A)) iterations.  Goal: M should be cheap to
apply AND make M^{-1} A well conditioned.

| Family | Examples | Use case | Cost / apply |
|--------|----------|----------|--------------|
| **Diagonal** | Jacobi | Always works, weak | O(N) |
| **Triangular** | SOR / SSOR | Better than Jacobi for spd | O(nnz(A)) |
| **Incomplete factor** | IC, MIC, ILU(0), ILU(k) | Classical default | O(nnz(L)) |
| **Sparse approximate inverse** | SPAI, FSAI, AINV | GPU-friendly | O(nnz(M^{-1})) |
| **Algebraic multigrid** | Ruge-Stuben, Smoothed Aggregation, BoomerAMG | Large elliptic | O(N) |
| **Auxiliary space** | AMS (Hiptmair-Xu), ADS | H(curl), H(div) ★ | O(N) |
| **Domain decomposition** | BDDC, FETI-DP, Schwarz | Distributed | O(N/p) |
| **Block** | Block-Jacobi, Block-Gauss-Seidel | Saddle-point | O(N) |
"""


CLASSICAL = r"""
# Classical preconditioners — Jacobi / SOR / IC / MIC / ILU

## Jacobi (diagonal)
M = diag(A).  Very cheap (O(N) apply), but only works for
strongly-diagonally-dominant A.  Default fallback when nothing else
applies.

## (S)SOR — Symmetric Successive Over-Relaxation
M_SSOR = (D + ω L) D^{-1} (D + ω L^T)  with ω ∈ (0, 2).
Slightly better than Jacobi for spd; convergent for ω close to optimal.

## IC / MICCG genealogy (the lab's classical baseline)

[LOCAL] 04_preconditioner_classical/ICCG/01_Meijerink_vanderVorst_1977_original_IC.pdf
Meijerink-vanderVorst 1977 — the original IC paper.  Symmetric M-matrices.

For A = L_A D_A L_A^T spd, the IC(0) factorization computes L̃ with the
**sparsity pattern of A** by dropping any fill-in:

```
L̃[i,j] = 0  for (i,j) ∉ pattern(A)
otherwise standard Cholesky
```

Then M = L̃ L̃^T ≈ A.  Apply: forward + back substitute on triangular
L̃ (sparse).

[LOCAL] 04_preconditioner_classical/ICCG/02_Kershaw_1978_ICCG.pdf
Kershaw 1978 — popularized ICCG in plasma physics.

[LOCAL] 04_preconditioner_classical/ICCG/03_MICCG_modified.pdf
Gustafsson 1978 — Modified IC: subtract row-sum of dropped entries
from diagonal.  Improves condition number of M^{-1} A from O(h^{-2})
to O(h^{-1}) for 2nd-order elliptic problems.

## Shift trick: when IC factorization fails (negative diagonal)

If A is indefinite (or has zero diagonal entries — see Periodic-Omega
in CLAUDE.md), Cholesky breaks down.  Standard fix: factor (A + α·M)
instead, with M = mass and α small.

```python
import radia.sparsesolv_ngsolve as ssn
prec = ssn.IC_AutoShift(a.mat, max_tries=8)   # auto-bisect shift
```

## ILU(k) — Incomplete LU with level-of-fill

For non-symmetric A.  ILU(0) keeps sparsity of A; ILU(k) allows
fill-in up to level k.  Trade-off: larger k = better preconditioner
but more memory and apply cost.

NGSolve default for non-symmetric problems via:
```python
c = Preconditioner(a, "local")           # block Jacobi
c = Preconditioner(a, "bddc")            # BDDC (best)
# ILU is not directly exposed; via SparseSolvSolver wrapper
```

## When classical works / fails

| Problem | Classical works? | Better alternative |
|---------|------------------|---------------------|
| 2D Poisson, structured grid | ✓ ICCG | (no need to upgrade) |
| 3D Poisson, unstructured | △ ICCG (slow convergence O(h^{-1})) | AMG |
| 3D HCurl curl-curl | ✗ ICCG (κ blows up at high p) | **AMS (Hiptmair-Xu)** ★ |
| 3D anisotropic | △ (depends on alignment) | AMG with smoothed aggregation |
"""


AMG = r"""
# Algebraic Multigrid (AMG)

Geometric multigrid (GMG) requires a hierarchy of meshes; AMG builds
the hierarchy **from the matrix entries alone** — no geometry needed.
Critical when meshes are unstructured.

## Ruge-Stuben classical AMG (1987)

Reference: Ruge-Stuben, "Algebraic Multigrid", in McCormick (ed)
"Multigrid Methods", SIAM Frontiers vol 3, ch 4, 1987.

Steps:
1. **Coarsen**: partition DOFs into C (coarse) and F (fine) sets via
   strength-of-connection threshold θ (usually 0.25):
   `i strongly connected to j ⇔ |a_ij| ≥ θ · max_k |a_ik|`
2. **Interpolate**: build P (F → C+F) using strong connections.
3. **Coarse operator**: A_c = P^T A P (Galerkin).
4. **Smooth + V-cycle**: pre-smooth (Jacobi/GS), restrict residual,
   solve coarsest exactly, prolong, post-smooth.

Works well for **M-matrices** (Poisson, diffusion).  Struggles with
strong anisotropy and 3D HCurl curl-curl (need AMS instead).

## BoomerAMG (Henson-Yang 2002)

Reference: V. Henson, U. Yang, "BoomerAMG: A parallel algebraic
multigrid solver and preconditioner", Appl. Num. Math. 41(1):155-177,
2002.  DOI: 10.1016/S0168-9274(01)00115-5.

Industrial-strength parallel implementation of Ruge-Stuben (in C with
MPI).  Ships in **HYPRE** (LLNL).  De-facto AMG reference.

Key options:
- `coarsen_type`: 6 (Falgout, default), 8 (PMIS, parallel), 10 (HMIS)
- `interp_type`: 0 (classical), 6 (extended+i, default), 13 (FF1)
- `relax_type`: 6 (hybrid sym-GS, default), 18 (l1-Jacobi, parallel-safe)
- `theta`: 0.25 (strong connection threshold)

## Smoothed Aggregation AMG (Vanek-Mandel-Brezina 1996)

Alternative coarsening: group DOFs into aggregates, smooth tentative
prolongator with one Jacobi step.  Better for elasticity (vector PDEs)
and anisotropic problems.

Available in NGSolve via `Preconditioner(a, "h1amg")` (only for H1
scalar problems).

## CompactAMG ★ (Sugahara lab, src/ext/sparsesolv)

Reference implementation: `src/ext/sparsesolv/compact_amg.hpp`.

Lightweight Ruge-Stuben variant designed for NGSolve TaskManager:
- PMIS coarsening (parallel, deterministic given hash seed)
- Classical interpolation
- l1-Jacobi smoother (fully parallel, no triangular solve)
- V-cycle, ω = 1.0

Performance (mesh1_3.5T, 197k DOFs HCurl):
- CompactAMG (Compact HX subspace solver): 25 iterations to 1e-10
- HYPRE BoomerAMG (subspace): 25 iterations (matches)
- Memory: ~0.3 GB (vs HYPRE ~0.8 GB including MPI buffers)

## Anisotropy and weak scaling

Pure AMG (Ruge-Stuben) coarsens by strength-of-connection.  In
anisotropic problems (e.g. thin laminates aligned with one axis),
the strong direction is over-coarsened and the weak direction
under-coarsened → poor convergence.

Mitigations:
- Smoothed Aggregation with vector node aggregation.
- Compatible relaxation coarsening.
- Use AMG only on the **scalar subspace** (nodes), AMS for vector.
"""


AMS_HIPTMAIR_XU = r"""
# Auxiliary Space Maxwell Solver (AMS / HX) — Hiptmair-Xu 2007 ★

★ THIS IS THE LAB CORE PRECONDITIONER FOR HCurl PROBLEMS.

[LOCAL] 05_AMG_AMS/01_Hiptmair_Xu_2007_Nodal_Auxiliary_Space_Hcurl_Hdiv.pdf

Reference: R. Hiptmair, J. Xu, "Nodal Auxiliary Space Preconditioning
in H(curl) and H(div) Spaces", SIAM J. Numer. Anal. 45(6):2483-2509,
2007.  DOI: 10.1137/060660588.

## The problem AMS solves

The HCurl curl-curl system
```
A x = b   with A = ν · curl_h^T curl_h + α · M_h     (HCurl edge elements)
```
has a **huge** null space:  curl of any H1 scalar gradient is zero.
Classical preconditioners (Jacobi, IC, AMG-on-edges) fail because
they don't know about this kernel.

Hiptmair 1998 (predecessor) — first multigrid for H(curl), but
requires geometric multigrid.

## Hiptmair-Xu 2007 — the auxiliary space trick

The HCurl space has two **auxiliary** scalar subspaces:
1. **Gradient subspace**: G : H1 → HCurl,  G·φ = ∇φ
2. **Nodal coordinate**: Π : (H1)^3 → HCurl, Π·v = Π_h(v)  (nodal interpolation of vector)

Preconditioner application:
```
M_HX^{-1} v = S_h v                                 (fine smoother)
            + G * S_∇^{-1} * G^T v                  (gradient subspace correction)
            + Π * (S_x^{-1}, S_y^{-1}, S_z^{-1}) * Π^T v   (nodal subspace correction)
```

Each S is itself a standard preconditioner (Jacobi, AMG, etc.) on the
**scalar** auxiliary problem.  This means **all the AMG machinery
ported to scalar problems automatically lifts to HCurl** — no
HCurl-specific multigrid needed.

## Why this is so powerful

| Property | Classical IC on HCurl | Hiptmair-Xu AMS |
|----------|----------------------|------------------|
| Handles kernel | ✗ (κ blows up) | ✓ (gradient subspace = kernel) |
| Scales with order p | ✗ degrades fast | ✓ stable up to p=10 |
| Requires geometric MG | n/a | ✗ (uses AMG on scalar) |
| Implementation effort | low | medium (need G, Π matrices) |

## Lab implementation: CompactAMS ★

[CODE] `src/ext/sparsesolv/compact_ams.hpp`  (HYPRE-free, TaskManager-native)
[CODE] `radia.sparsesolv_ngsolve.CompactAMSPreconditioner`

Validated 2026-04 on complex eddy current @ 30 kHz, 155k-1.44M DOFs:

| Parameter | Value | Note |
|-----------|-------|------|
| Cycle type | 1 (01210) | pre-smooth, G-correct, Pi-correct, G-correct, post-smooth |
| Outer solver | COCR (complex) / CG (real) | from radia.sparsesolv_ngsolve |
| Fine smoother | l1-Jacobi | TaskManager-parallel |
| Subspace solver | CompactAMG | own AMG, no HYPRE |
| Pi mode | Separate Pix / Piy / Piz | multiplicative |
| AMG theta | 0.25 | strength-of-connection |
| Correction weight | 1.0 | no damping |

Performance vs HYPRE AMS + BoomerAMG (mesh1_3.5T, 197k DOFs):
- CompactAMS + CompactAMG: 25 BiCGStab iterations to 1e-10
- HYPRE AMS + BoomerAMG: 25 iterations (matches)

## When to use — MEASURED, and it depends on the ORDER p (2026-09-08)

**AMS is for p=1 (lowest-order Nedelec).  For p=2 use NGSolve `bddc`.**
Measured on mdx1 (idle, 38 cores) on the hiruma complex eddy-current problem
(conductor + core + air, σ=0 in air, 30 kHz, `nograds=True` tets, COCR, 1e-8):

| order | use | why (measured) |
|-------|-----|----------------|
| **p=1** | **AMS** ★ | 197k dof: AMS 523 MB / 2.9 s setup + 2.0 s solve, versus a direct factorisation 2029 MB / 6.4 s.  Memory scales **O(N^0.81)** versus **O(N^1.72)** |
| **p=2** | **NGSolve `bddc`** | 1.46M dof in 67 iterations / 5.5 s / 4.35 GB.  Iterations are flat in mesh refinement (95 → 96 → 112 over 680k → 1.46M) |
| any | NOT `multigrid` | p=2, 865k dof: 154 iterations and **205 s of solve versus bddc's 7.2 s — 28x slower** |

**AMS cannot be applied directly at p≥2**: its Π interpolation is built from
vertex coordinates and the lowest-order discrete gradient.

**Do not "extend AMS to p=2" — it is not needed.**  p=2 does not suffer a
convergence collapse; `bddc` handles it with flat iteration counts.

**Using AMS as the BDDC coarse solver was built and measured — it is a
TRADE-OFF, not a win.**  The p=2 BDDC wirebasket coarse space is 56% lowest-order
edge dofs (i.e. exactly the p=1 space, where AMS applies).  Substituting AMS
there gives **memory −63% but solve time +58%** at 1.46M dof.  Reason: a direct
method puts its cost in the SETUP (done once) and has a cheap APPLY; a Krylov
preconditioner applies ~100 times per solve, so apply cost dominates and the
direct factorisation wins.  AMS wins only where it is the SOLVER (1 apply) —
which is the p=1 row above.  Full data: `memory/hcurl_preconditioner_measurements_2026_09_08.md`.

**ICCG / IC is NOT usable as the coarse solver**: 1000 iterations without
converging (it cannot see the gradient kernel).  Compressed factorisations
(MUMPS BLR, STRUMPACK HSS, H-LU) are also a poor fit here — low-rank
compressibility requires a conductive medium between clusters, and a σ=0 AIR
region destroys it.

| Problem | Use CompactAMS |
|---------|----------------|
| Real spd HCurl curl-curl + mass, **p=1** | ✓ + CG |
| Complex sym HCurl (eddy current MQS), **p=1** | ✓ + COCR ★ |
| Same but **p=2 or higher** | ✗ — use NGSolve `bddc` |
| HCurl with air region (σ=0) | ✓ + shifted preconditioner (see em_specific) |
| HDiv (flux variable) | Use ADS — same paper, dual construction |
| HCurl helmholtz (high freq) | NOT a lab use case (Laplace kernel only) |

## ⚠ Build it OUTSIDE `with TaskManager():`

Constructing `CompactAMSPreconditioner` / `ComplexCompactAMSPreconditioner`
inside an active TaskManager region kills the process with `0xC0000409`
(`__fastfail`) — **no Python exception, no traceback**.  Verified 2026-09-08 on
the identical matrix: outside = OK, inside = dead; reproduced on two
independently built binaries.  The unit tests pass only because they do not
wrap.  This contradicts the CLAUDE.md "TaskManager-Only" policy, which tells
callers to wrap the whole NGSolve block — so wrap the mesh/space/forms/assembly
and the SOLVE, but build the preconditioner between them, outside the region.
See `memory/ams_constructor_dies_inside_taskmanager.md`.

## Code recipe

The real signature takes the discrete gradient and the vertex COORDINATES —
it does not take the FESpace.  It also needs a REAL surrogate matrix even for
a complex system.  Verified against the shipped binary 2026-09-08.

```python
import numpy as np
from ngsolve import HCurl, BilinearForm, LinearForm, GridFunction, TaskManager
from radia.sparsesolv_ngsolve import ComplexCompactAMSPreconditioner, COCRSolver

# order=1: AMS's Pi is built from the LOWEST-ORDER gradient. See "When to use".
kw = dict(order=1, nograds=True, dirichlet="dirichlet")

with TaskManager():
    fes = HCurl(mesh, complex=True, **kw)          # the complex system
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += nu*curl(u)*curl(v)*dx + eps*nu*u*v*dx + 1j*omega*sigma*u*v*dx("cond")
    f = LinearForm(fes); f += nu*CF((0, 0, 1))*v*dx("cond")
    a.Assemble(); f.Assemble()

    fes_r = HCurl(mesh, complex=False, **kw)       # REAL surrogate for the AMS
    ur, vr = fes_r.TnT()
    ar = BilinearForm(fes_r)
    ar += nu*curl(ur)*curl(vr)*dx + eps*nu*ur*vr*dx
    ar += abs(omega)*sigma*ur*vr*dx("cond")        # |omega|, not 1j*omega
    ar.Assemble()
    G, _ = fes_r.CreateGradient()
    pts = mesh.ngmesh.Points()
    cx, cy, cz = ([pts[i + 1][k] for i in range(mesh.nv)] for k in range(3))

# BUILD OUTSIDE the TaskManager region -- inside it dies with 0xC0000409
prec = ComplexCompactAMSPreconditioner(
    a_real_mat=ar.mat, grad_mat=G, freedofs=fes_r.FreeDofs(),
    coord_x=cx, coord_y=cy, coord_z=cz,
    ndof_complex=fes.ndof, cycle_type=1, print_level=0)

with TaskManager():
    gfu = GridFunction(fes)
    inv = COCRSolver(a.mat, prec, freedofs=fes.FreeDofs(), tol=1e-8, maxiter=2000)
    gfu.vec.data = inv * f.vec
```

For the REAL system use `CompactAMSPreconditioner(a.mat, G, freedofs=...,
coord_x=..., coord_y=..., coord_z=...)` and NGSolve's `CGSolver`.

See `radia_ngsolve` MCP tool `sparsesolv('compact_ams')` and
`sparsesolv('example_compact_ams')` for full recipes.
"""


SHIFTED_PRECONDITIONER = r"""
# Shifted Preconditioner for HCurl with air+conductor regions

CLAUDE.md POLICY (`Shifted Preconditioner for Air+Conductor Problems`):
For HCurl eddy current with air regions (σ=0), use **Shifted
Preconditioner** instead of system regularization.  Add ε·mass to the
PRECONDITIONER only, not the system matrix.

## The problem

In the system
```
A = ν · curl_h^T curl_h + jω·σ_h M_h
```
with σ = 0 in air, the mass term M_air = 0.  The system has a huge
null space on the air region (any solenoidal field).  CompactAMS still
needs the nodal subspace correction to be NON-SINGULAR.

## Fix: shift the preconditioner ONLY

```python
# Preconditioner: shifted (non-singular on whole domain)
a_shifted = BilinearForm(fes, symmetric=True)
a_shifted += SymbolicBFI(ν * curl(u) * curl(v))
a_shifted += SymbolicBFI(1j * ω * σ_cf * u * v, definedon=mesh.Materials("cond"))
a_shifted += SymbolicBFI(eps * ν * u * v)              # eps = 1e-6 * ν
a_shifted.Assemble()

# System: original (singular in air, but physically correct)
a = BilinearForm(fes, symmetric=True)
a += SymbolicBFI(ν * curl(u) * curl(v))
a += SymbolicBFI(1j * ω * σ_cf * u * v, definedon=mesh.Materials("cond"))
# NO eps here
a.Assemble()

# Build preconditioner from shifted, but solve original
prec = ComplexCompactAMSPreconditioner(a_shifted.mat, fes, ...)
solver = COCRSolver(a.mat, prec, tol=1e-10)
```

Verified (sparsesolv repo `examples/hiruma/shifted_ams_experiment.py`):
- ε from 1e-4 to 1e-8 give **identical** ||B||² → solution independent of ε
- Without shift: diverges (nan)

## Common mistake

Adding ε·M to **both** the system AND the preconditioner is the WRONG
fix.  It modifies the physics — you get a fictitious displacement
current in air.  Always shift the PRECONDITIONER only.
"""


def get_preconditioners_knowledge(topic: str = "catalog") -> str:
    """Dispatch preconditioner topics.

    Topics:
        catalog        - Preconditioner family overview (DEFAULT)
        classical      - Jacobi/SOR/IC/MIC/ILU (Meijerink-vdVorst 1977 family)
        amg            - Algebraic multigrid (Ruge-Stuben, BoomerAMG, CompactAMG)
        ams_hiptmair_xu - ★ AMS (Hiptmair-Xu 2007) — lab CompactAMS
        shifted        - Shifted preconditioner for air+conductor HCurl
        all            - Everything
    """
    topic = topic.lower().strip()
    if topic in ("catalog", "overview", "family"):
        return CATALOG
    if topic in ("classical", "iccg", "ic", "miccg", "ilu"):
        return CLASSICAL
    if topic in ("amg", "boomeramg", "ruge_stuben", "compact_amg"):
        return AMG
    if topic in ("ams", "ams_hiptmair_xu", "hiptmair_xu", "compact_ams", "hx"):
        return AMS_HIPTMAIR_XU
    if topic in ("shifted", "shifted_prec", "air_conductor"):
        return SHIFTED_PRECONDITIONER
    if topic == "all":
        return "\n\n".join([CATALOG, CLASSICAL, AMG, AMS_HIPTMAIR_XU,
                            SHIFTED_PRECONDITIONER])
    return (f"Unknown topic '{topic}'. Available: catalog, classical, amg, "
            "ams_hiptmair_xu, shifted, all.")
