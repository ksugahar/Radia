"""EM-specific solver / formulation choices.

Tree-cotree gauging, Biro-Preis A-V, ungauged formulations, eddy-current
low-frequency stabilization.
"""

GAUGING = r"""
# Gauging electromagnetic potentials

For HCurl-A formulations of static / quasi-static Maxwell, A is
defined only up to a gradient: A' = A + ∇φ gives the same B.  This
gauge freedom shows up as a **null space** in the discrete operator.

| Strategy | What it does | When to use |
|----------|--------------|-------------|
| **Tree-cotree** (Albanese-Rubinacci 1988) | Set A to zero on a spanning tree of edges | Topologically connected conductor with prescribed boundary current |
| **Coulomb gauge** ∇·A = 0 | Add penalty (∇·u, ∇·v) | Simple but adds positive operator |
| **A-V gauge** (Biro-Preis 1989) | Add scalar V; solve coupled (A, V) | Multiply-connected, mixed conductor / nonconductor |
| **Ungauged + Krylov** | Live with null space, use AMS that knows the kernel | Hiptmair-Xu AMS handles this automatically |
| **Two-Step** (Ostrowski-Hiptmair 2021) | Generalized tree-cotree, frequency-stable | Full Maxwell DC to MHz |

## Why the modern lab uses ungauged + AMS

[CODE] `radia.sparsesolv_ngsolve.CompactAMSPreconditioner`

The Hiptmair-Xu AMS preconditioner has a **gradient subspace
correction**: it applies a scalar AMG solve on H1 to handle the
gradient null space.  Result: no gauge is needed at the system level
— the preconditioner absorbs the kernel.

Tradeoff:
- Gauged formulation: smaller (no extra scalar), system is SPD/definite.
- Ungauged + AMS: easier to set up, no tree-finding code, but iterations
  cost ~10% more.

**Default**: ungauged + CompactAMS unless geometry-specific reason.
"""


TREE_COTREE = r"""
# Tree-cotree gauge (Albanese-Rubinacci 1988)

For HCurl edge basis on a mesh with N_v vertices and N_e edges:
- Spanning tree T has N_v - 1 edges
- Cotree T̄ has N_e - (N_v - 1) edges

**Gauge**: set A_h on tree edges to zero.  Remaining DOFs are on cotree
edges.  System is now SPD and non-singular.

## Implementation in NGSolve

```python
from ngsolve import HCurl, Periodic, GridFunction
fes_full = HCurl(mesh, order=p)

# Build tree (use NetworkX or Netgen's spanning tree)
# Convert to Dirichlet DOFs:
tree_dofs = compute_spanning_tree_dofs(mesh, fes_full)
fes = HCurl(mesh, order=p, dirichletx_bbnd=tree_dofs)
```

NGSolve does NOT ship a built-in tree finder — must compute externally.
Hiptmair has an "ngscxx" example.

## Two-Step / Frequency-Stable formulation (Ostrowski-Hiptmair 2021)

[LOCAL] 06_EM_specific/01_TwoStep_Maxwell_TreeCotree_LowFreq_Stability.pdf

Reference: J. Ostrowski, R. Hiptmair, "Frequency-Stable Full Maxwell
in Electro-quasistatic Gauge", SIAM J. Sci. Comput. 43(6), 2021.
DOI: 10.1137/20M1356300.

Problem: standard A-V formulation has ω→0 stability issues — the mass
term jωσM vanishes, and the curl-curl part becomes singular.

Fix: split the Maxwell A-V system into TWO steps:
1. Solve electro-quasistatic gauge: yields V exactly at ω=0
2. Solve full Maxwell on the gauged A: stable across all ω

Uses generalized tree-cotree to identify which DOFs live in each step.
"""


BIRO_PREIS_AV = r"""
# A-V Reduced Potential Formulation (Bíró-Preis)

Reference: O. Bíró, K. Preis, "An edge finite element eddy current
formulation using a reduced magnetic and a current vector potential",
IEEE Trans. Magn. 36(5):3128-3130, 2000.  DOI: 10.1109/20.908708.

Earlier: Bíró-Preis 1989 IEEE IntMag (multiply-connected eddy current).

## The formulation

For eddy current with multiply-connected conductors:
- In conductor: solve for **A** (HCurl) and **V** (H1) coupled
- In non-conductor: solve for **A_r** (reduced HCurl, A_r = A - A_source)
- Continuity of B·n at interface enforced

System (frequency domain, time-harmonic ω):
```
[ S_νν + jω σ M_AA    jω σ G_AV ] [A_r]   [F]
[                                 ] [   ] = [ ]
[ jω σ G_AV^T          σ K_VV    ] [ V ]   [0]
```

where S_νν = (ν · curl·, curl·), M_AA = (σ ·, ·), G_AV = (σ ∇·, ·),
K_VV = (σ ∇·, ∇·).

## Why use it (over pure A or pure T-Ω)

- **Pure A**: needs gauge (tree-cotree or AMS).
- **T-Ω**: requires cohomology basis (cuts) for multiply connected — hard.
- **A-V**: V handles the loop-current ambiguity naturally; no cut-finding.

## Implementation status

Not directly shipped in `radia.sparsesolv_ngsolve`.  NGSolve users
write the bilinear form themselves; the system is solved with
CompactAMS on the A-block + direct on V (Schur complement) or with
HYPRE AMS-aware preconditioner.

For pure A formulation (most lab use), ungauged + CompactAMS is
the default — see `preconditioners_knowledge.AMS_HIPTMAIR_XU`.
"""


EDDY_CURRENT_STABILIZATION = r"""
# Low-frequency stabilization for eddy current MQS

The MQS HCurl system
```
A = ν · curl·curl + jω · σ · M
```
has two stability limits:
1. **ω → 0**: mass term vanishes, A becomes singular (curl·curl has
   gradient null space).
2. **σ → 0 in air regions**: M·u = 0 there, system is singular.

## Fix #1 (ω → 0): use a gauged or AMS-aware formulation

- Tree-cotree gauge: forces A = 0 on tree, eliminates kernel.
- A-V Bíró-Preis: V absorbs the kernel.
- Ungauged + CompactAMS: AMS gradient subspace handles it.

## Fix #2 (σ → 0 in air): shifted preconditioner

See `preconditioners_knowledge.SHIFTED_PRECONDITIONER`:
add ε·M to the PRECONDITIONER only (not the system).

## Combined practice (lab default)

```python
# Original physics (correct, no perturbation)
a += ν * curl(u) * curl(v) * dx
a += 1j * ω * σ_cf * u * v * dx("cond")

# Preconditioner with shift for ε > 0 over WHOLE domain
a_shifted += ν * curl(u) * curl(v) * dx
a_shifted += 1j * ω * σ_cf * u * v * dx("cond")
a_shifted += eps * ν * u * v * dx       # eps = 1e-6 * ν

# Compact AMS on shifted; COCR on original
prec = ComplexCompactAMSPreconditioner(a_shifted.mat, fes, ...)
solver = COCRSolver(a.mat, prec, tol=1e-10)
```

This is the **production lab configuration** for IH workpiece SIBC,
PEEC+FEM coupling, and the radia_motor lamination effective-material
solver.  See `radia_ngsolve.sparsesolv('compact_ams')`.
"""


def get_em_specific_knowledge(topic: str = "gauging") -> str:
    """Dispatch EM-specific solver topics.

    Topics:
        gauging              - Gauge choices overview (DEFAULT)
        tree_cotree          - Albanese-Rubinacci + Ostrowski-Hiptmair Two-Step
        biro_preis_av        - A-V reduced potential formulation
        eddy_current_stab    - Low-frequency + air-region stabilization
        all                  - Everything
    """
    topic = topic.lower().strip()
    if topic in ("gauging", "gauge", "overview"):
        return GAUGING
    if topic in ("tree_cotree", "tree", "cotree", "two_step", "twostep"):
        return TREE_COTREE
    if topic in ("biro_preis_av", "biro_preis", "av_gauge", "av", "biropreis"):
        return BIRO_PREIS_AV
    if topic in ("eddy_current_stab", "eddy_current", "stabilization",
                 "low_frequency", "lowfreq"):
        return EDDY_CURRENT_STABILIZATION
    if topic == "all":
        return "\n\n".join([GAUGING, TREE_COTREE, BIRO_PREIS_AV,
                            EDDY_CURRENT_STABILIZATION])
    return (f"Unknown topic '{topic}'. Available: gauging, tree_cotree, "
            "biro_preis_av, eddy_current_stab, all.")
