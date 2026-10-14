"""
Finite Element Exterior Calculus (FEEC): the modern unified framework
that generalizes Whitney elements to all polynomial orders and shows
why they make FEM stable.

Distilled from:
  - Arnold-Falk-Winther 2006, "Finite element exterior calculus,
    homological techniques, and applications", Acta Numerica 15
    (the canonical reference)
  - Bossavit's later papers (1998, 2001) reformulated in FEEC language
  - Whitney 1957 (the historical origin of the W^k basis)
"""

FEEC_OVERVIEW = """
# Finite element exterior calculus (FEEC) — Arnold-Falk-Winther 2006

## Central thesis
PDE stability is governed by **structural** properties of the
underlying differential complex — geometry, topology, and homology.
Numerical methods that preserve these structures (the "cochain
projection condition") are stable; those that don't are not.

Many seemingly-different FEM choices (Lagrange P_k, Raviart-Thomas
RT_k, Nedelec N_k I and II, Brezzi-Douglas-Marini BDM_k) fall into a
unified family of **finite-element differential forms**.  Two families
exist for each polynomial degree r ≥ 1:

    P_r Λ^k(T_h)        and        P_r^- Λ^k(T_h)

(read: "P r Lambda k", and "P r minus Lambda k").  Each is a space of
piecewise polynomial differential k-forms on the mesh T_h.

## Classical correspondences

| k=0 (scalars)      | k=1 (vectors)         | k=2 (vectors)            | k=3 (scalars) |
|--------------------|------------------------|--------------------------|---------------|
| P_r Λ^0 = Lagrange P_r | P_r Λ^1 = Nedelec II of degree r | P_r Λ^2 = BDM_r | P_r Λ^3 = discontinuous P_r |
| P_r^- Λ^0 = Lagrange P_r | P_r^- Λ^1 = Nedelec I of degree r | P_r^- Λ^2 = RT_r | P_r^- Λ^3 = discontinuous P_{r-1} |

So Whitney elements (P_1^- Λ^k) and Raviart-Thomas (P_r^- Λ^2 in 2D)
and Nedelec edge elements (P_r^- Λ^1) all live in the same family.

## The two families differ by HOMOGENEITY
- P_r Λ^k:  piecewise polynomials of full degree r.
- P_r^- Λ^k:  "trimmed" — a subspace such that  dim(P_r^- Λ^k) = dim P_r Λ^k − dim P_r Λ^{k+1}.  Combinatorial reduction by Koszul complex.

This trimming is crucial for the discrete de Rham complex:

    P_r Λ^0 →^d  P_{r-1} Λ^1 →^d  P_{r-2} Λ^2 →^d  ...        (lose 1 degree per d)
    P_r^- Λ^0 →^d  P_r^- Λ^1 →^d  P_r^- Λ^2 →^d  ...           (stays at "−r")

The first sequence is conventional polynomial reduction.  The second
gives an **isomorphic copy** of the de Rham complex at fixed
nominal degree.  This isomorphism is what makes FEEC powerful.
"""

CHAIN_PROJECTION = """
# Cochain projections: the discrete-continuous bridge

## What is a cochain projection?
A linear map  Π_h : Λ^k(Ω) → V_h^k (the finite-element subspace)
such that

    d Π_h = Π_h d           (commuting square)
    Π_h v = v   for v ∈ V_h  (projection / fixed-point on V_h)

Pictorially:
    Λ^k   ──d──→   Λ^{k+1}
     │ Π_h           │ Π_h
     ▼               ▼
    V_h^k ──d──→  V_h^{k+1}

The horizontal d is the **same** operator on both rows (with the
discrete d restricted to the FE subspace giving the relevant matrix).

## Why it matters
- Stability: a cochain projection gives a uniform Fortin operator,
  giving inf-sup conditions (LBB) for free.
- Convergence: ||v - V_h v|| → 0 as h → 0 at the right rate, AND
  ||d v - d V_h v|| → 0 at the same rate — without loss.
- Spectral accuracy: discrete eigenmodes converge to continuous
  eigenmodes without spurious modes (provided the chain projection
  exists).

## Constructions of cochain projections
- Canonical projection: integrate v over the appropriate moment
  (node value for P_r Λ^0, edge moment for P_r Λ^1, etc.).
  Works for smooth v but NOT bounded uniformly on H^k spaces.
- Schöberl 2005 / Christiansen 2005: **smoothed** projections.
  First mollify v (regularize) and then apply the canonical
  projection.  Uniformly bounded.
- Arnold-Falk-Winther 2006: combine both approaches.

## Why classical FEM choices SOMETIMES fail
- Nodal (P_1) elements for the curl-curl problem with mixed
  boundary conditions on a non-convex domain converge to the
  **wrong solution**.  Reason: the FE subspace V_h is in H^1
  but the true solution is in H(curl) but NOT H^1.  No cochain
  projection exists.
- Whitney elements (P_1^- Λ^1) for the same problem converge to the
  correct solution.  The cochain projection exists because Whitney
  elements properly approximate H(curl), not the strict subspace H^1.

This is the **fundamental answer** to "why edge elements" — and it
took 25 years (1980 → 2002) to be fully formulated.
"""

KOSZUL_COMPLEX = """
# Koszul complex: the algebraic engine

## Setup
Let Λ^k = Λ^k(R^n) be the space of polynomial k-forms on R^n.  The
**Koszul differential** is a map κ : Λ^k → Λ^{k-1} defined on a
monomial form:

    κ(x_{i_1} ... x_{i_m} dx^{j_1} ∧ ... ∧ dx^{j_k})
        = ∑_p (-1)^{p-1} x_{j_p} · x_{i_1} ... x_{i_m}
                 · dx^{j_1} ∧ ... ∧ dx̂^{j_p} ∧ ... ∧ dx^{j_k}

In coordinate-free language:  κ ω = X⌟ω  where X is the radial vector
field  x · ∂/∂x  (i.e., the Euler vector field).

## Properties
- κ² = 0 (Koszul-nilpotent)
- κ d + d κ = (r + k) id   (on polynomial forms of degree r)

The latter identity (the Cartan magic formula) is the algebraic core.

## The P_r^- Λ^k space
    P_r^- Λ^k(R^n) := { ω ∈ Λ^k : ω is polynomial of degree ≤ r, AND
                         κ ω is polynomial of degree ≤ r }

Equivalent: P_r^- Λ^k = κ(Λ^{k+1}_{r-1}) ⊕ d(Λ^{k-1}_{r-1})  (a
direct decomposition, by the magic formula).

This characterizes the "minus" family.  The full family is just
P_r Λ^k = polynomial k-forms of degree ≤ r.

## Concrete consequences
- P_1^- Λ^k = Whitney elements (lowest order).
- dim P_r^- Λ^k = (n + r choose n − k)(r + k choose k) / (k + 1)
- The Koszul complex is the homological device that proves
  P_r^- Λ^k is **affinely invariant** (a finite-element transferable
  space).

## Why "affinely invariant" matters for FEM
On a mesh, each tetrahedron has its own coordinate system.  The
reference tetrahedron map is an affine transformation.  Spaces that
are invariant under affine pullback can be assembled element-by-
element on the mesh in a consistent way.  Arnold-Falk-Winther
showed that  P_r Λ^k  and  P_r^- Λ^k  are **essentially the only**
spaces of polynomial differential forms with this property.

So **the design space of stable finite elements is finite** and
classified.  No surprises lurking in undiscovered families.
"""

HODGE_LAPLACIAN = """
# The Hodge Laplacian and the abstract FEM template

## Abstract setup
Let  V^0 → V^1 → ... → V^n  be the Sobolev de Rham complex
H^1 / H(curl) / H(div) / L^2.  The Hodge Laplacian on level k is

    L = d δ + δ d                  on Λ^k

where δ = d^* is the formal adjoint of d (= ±⋆d⋆ depending on sign
conventions).

## Source problem (mixed form)
Given f ∈ V^k, find σ ∈ V^{k-1}, u ∈ V^k, p ∈ harmonic_k such that

    (σ, τ) − (d τ, u) = 0           for all τ ∈ V^{k-1}
    (d σ, v) + (d u, d v) + (p, v) = (f, v)   for all v ∈ V^k
    (u, q) = 0                       for all q ∈ harmonic_k

This is the **mixed variational formulation** of the Hodge Laplacian.
It captures many standard problems in a unified framework:

| k | What this is |
|---|--------------|
| 0 | The standard Poisson problem in mixed (σ = grad u) form |
| 1 | The vector Poisson problem  −Δ u = f  in H(curl)/H¹ mixed form |
| 2 | The vector Poisson problem in H(div)/H(curl) mixed form |
| 3 | Trivial (just L²) |

## Galerkin discretization
Pick FE subspaces V_h^{k-1} ⊂ V^{k-1} and V_h^k ⊂ V^k that form a
discrete subcomplex (d V_h^{k-1} ⊂ V_h^k) and admit cochain
projections.  Discretize the variational form by Galerkin
projection.

**Result (Arnold-Falk-Winther 2006, Theorem 7.4)**: the discrete
problem is uniformly stable and convergent at the rate determined by
the approximation error of the spaces.  No inf-sup checking needed
case-by-case — the cochain projection condition alone suffices.

## Practical instances
- Maxwell eigenvalue problem: avoid spurious modes by choosing
  V_h^1 ⊂ V_h^0 ⊕ V_h^2 in a balanced way (Whitney elements +
  discontinuous L² for the constraint).
- Mixed Poisson: RT_k for σ, discontinuous P_{k-1} for u (RT-DG).
- Vector Laplacian: Nedelec for σ, Lagrange for u.

## Connection to preconditioning
Arnold-Falk-Winther 2006 §10 shows that effective preconditioners
for the Hodge Laplacian can be built using the de Rham subcomplex
structure: auxiliary-space preconditioners (Hiptmair-Xu 2007 AMS,
known to Radia users as "Compact HX") leverage exactly this
structure.  See `compact_ams.hpp` in src/ext/sparsesolv/.
"""


HIPTMAIR_ROBUST_MAXWELL = """
# Robust Maxwell formulation for all frequencies (Hiptmair-Krämer-Ostrowski 2008)

The classical Coulomb-gauged A-V formulation fails at low frequency:
the variational system loses stability as ω → 0, even though the
physics (electrostatics + magnetostatics) is perfectly well-posed.

Hiptmair-Krämer-Ostrowski 2008 (IEEE Trans. Mag. 44(6):682) propose
a **generating-system approach** that fixes this.

## Low-frequency instability

Standard A-V formulation in time-harmonic form:
    curl(ν curl A) + jω σ (A + grad V) = j^impressed   in Ω
    div(σ(A + grad V)) = 0                              in Ω_c
    div(ε(A + grad V)) = ρ                              in Ω_a

When ω → 0, Gauss's law div(ε grad V) = ρ becomes an INDEPENDENT
equation but is treated as a consequence in the variational form.
Conditioning blows up; roundoff dominates the solution in air.

Empirical (Hiptmair Fig. 2): below 1 kHz, the imaginary part of V
in non-conducting region is randomly disturbed, even with a
direct sparse solver.

## Generating-system fix

Introduce a redundant scalar q on the non-conducting domain with
the redundant equation
    test the 1st PDE with grad q', ∀ q' ∈ H¹ (q'=0 on Ω_c)

Redundant for any ω ≠ 0; at ω = 0 it BECOMES Gauss's law in air.

Result: a **singular but consistent** linear system that iterative
solvers (BiCGstab) handle naturally.

## Stationary limit (ω = 0)

The augmented system decouples PERFECTLY:
- Electrostatics in air:        -div(ε grad V) = ρ
- DC conduction in conductor:   -div(σ grad V) = 0
- Magnetostatics:               curl(ν curl A) = j^impressed (Coulomb gauge)

All well-posed by classical theory.

## Iterative solution

Hiptmair Fig. 5: BiCGstab + preliminary preconditioner converges in
3-5 iterations to machine precision for ALL frequencies including
ω = 0.  Robust over 0 Hz to 1 MHz.

## Why this matters for Radia + NGSolve

- IH problems span DC to ~500 kHz (6 orders of magnitude)
- Accelerator magnets are quasi-static but may need DC limit
- PEEC at low frequency has similar issues; same singular-iterative
  fix applies

## Connection to FEEC

This is a "redundant discretization" — adding a variable that's
mathematically dependent on existing ones but numerically separates
the ill-conditioning.  Recurring theme:
- Whitney elements have built-in redundancies (kernel of curl is
  non-trivial unless gauge-fixed) — iterative solvers exploit
- Arnold-Falk-Winther FEEC shows redundant / null-space-aware
  discretizations are systematically more stable
- Hiptmair 2008 = prototype for industrial-scale Maxwell solvers

Reference: R. Hiptmair, F. Krämer, J. Ostrowski, "A Robust Maxwell
Formulation for All Frequencies", IEEE Trans. Mag. 44(6):682-685,
June 2008.
"""


def get_feec_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "overview"          — P_r Λ^k vs P_r^- Λ^k families, classical equivalents
      "cochain_projection" — commuting diagram, why Whitney is correct
      "koszul"            — Koszul complex, algebraic engine, affine invariance
      "hodge_laplacian"   — abstract Hodge problem, mixed FEM, applications
      "hiptmair_robust"   — Hiptmair-Krämer-Ostrowski 2008 robust low-frequency
                            Maxwell formulation via generating systems
    """
    topic = topic.lower().strip()
    if topic == "overview":
        return FEEC_OVERVIEW
    if topic in ("cochain_projection", "projection", "commuting"):
        return CHAIN_PROJECTION
    if topic in ("koszul", "koszul_complex", "algebra"):
        return KOSZUL_COMPLEX
    if topic in ("hodge_laplacian", "hodge_lap", "laplacian", "mixed"):
        return HODGE_LAPLACIAN
    if topic in ("hiptmair_robust", "robust_maxwell", "low_frequency",
                 "all_frequencies", "hiptmair"):
        return HIPTMAIR_ROBUST_MAXWELL
    if topic == "all":
        return "\n\n".join([
            FEEC_OVERVIEW, CHAIN_PROJECTION, KOSZUL_COMPLEX,
            HODGE_LAPLACIAN, HIPTMAIR_ROBUST_MAXWELL,
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, overview, cochain_projection, koszul, hodge_laplacian, "
        "hiptmair_robust."
    )
