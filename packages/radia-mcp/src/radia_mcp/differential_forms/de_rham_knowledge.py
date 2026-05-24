"""
de Rham complex, Sobolev spaces of differential forms, the "Maxwell house"
diagram, and the connection between weak operators and integration by parts.

Distilled from:
  - Bossavit 1998 Ch.5 §5.1 (weak grad/rot/div, L^2_grad / H(curl) / H(div),
    Maxwell's house figure 5.1)
  - Arnold-Falk-Winther 2006 §2.2-2.3 (de Rham complex, Hodge theory)
  - Kameari 2011 figs. 3-4 (the same de Rham complex with the 3-D layout)
"""

DE_RHAM_COMPLEX = """
# The de Rham complex

## Setup
On an open domain Ω ⊂ R^n (or a smooth n-manifold), let Λ^k(Ω)
denote the space of smooth k-forms.  The exterior derivative
d : Λ^k → Λ^{k+1}  satisfies  d² = 0.  The de Rham complex is:

    0 → Λ^0(Ω) →^d Λ^1(Ω) →^d Λ^2(Ω) →^d ... →^d Λ^n(Ω) → 0.

In 3-D this reads:
                grad           rot            div
    0 → C^∞(Ω) ─→ C^∞(Ω)^3 ──→ C^∞(Ω)^3 ──→ C^∞(Ω) → 0.

## de Rham cohomology
The cohomology groups
    H^k_dR(Ω) := ker(d : Λ^k → Λ^{k+1}) / im(d : Λ^{k-1} → Λ^k)
measure the failure of exactness, and their dimensions are the
**Betti numbers** b_k(Ω) — topological invariants of Ω.

For a contractible Ω, all H^k_dR are trivial (Poincaré lemma).
For a torus T^3:  b_0=1, b_1=3, b_2=3, b_3=1.

This is the foundational result of de Rham (1931): differential
forms detect topology.

## Why is this useful for EM?
The kernel of rot (curl-free fields) = image of grad (gradients)
ONLY when the domain is contractible.  In a domain with loops
(e.g. the interior of a coil), curl-free fields that are NOT
gradients exist — these are precisely the "cuts" needed for scalar
potential formulations of magnetostatics with current loops.

The number of independent cuts = b_1(D) = number of independent
loops you need to break to make D contractible.

Concrete example:  inside a torus the field of "1/(2π r) e_φ"
(the magnetic field of an infinitely thin wire on the torus axis)
is curl-free but NOT a gradient — you cannot define a single-valued
magnetic scalar potential.  You must cut the torus with a surface.
"""

WEAK_OPERATORS = """
# Weak operators and Sobolev spaces of forms

## Why weakening is needed
Strong  grad/rot/div  require differentiability of the input.  But
many physical fields have jumps (across material interfaces) and
distributional sources (point charges, line currents).  The strong
operators cannot handle these.

## Definition (3-D illustration)
A vector field b ∈ L^2(Ω) has **weak divergence**  div b ∈ L^2(Ω)
if there exists a scalar field f ∈ L^2(Ω) such that

    ∫_Ω b · grad φ + ∫_Ω f φ = 0       for all  φ ∈ C^∞_0(Ω).

(The IBP boundary term vanishes because φ has compact support in Ω.)

Then by definition  div b := f.  Same trick for weak grad and weak rot.

## Functional spaces (Bossavit notation ←→ Sobolev notation)

| Sobolev | Bossavit | Definition |
|---------|----------|------------|
| H^1(Ω)         | L²_grad(Ω) | { φ ∈ L²(Ω) : grad φ ∈ L²(Ω)³ } |
| H(curl;Ω)      | IL²_rot(Ω) | { u ∈ L²(Ω)³ : rot u ∈ L²(Ω)³ } |
| H(div;Ω)       | IL²_div(Ω) | { u ∈ L²(Ω)³ : div u ∈ L²(Ω) }  |
| L²(Ω)          | L²(Ω)      | (square-integrable functions)   |

Each has its **graph norm** (sum of input + weak-derivative L² norms);
they are Hilbert spaces.

## The sequence with weak operators
                  grad             rot              div
    0 → H¹(Ω) ─→ H(curl;Ω) ──→ H(div;Ω) ──→ L²(Ω) → 0

The Poincaré lemma extends:  for a bounded regular contractible Ω
this sequence is **exact** at all interior levels.  This is what
makes the variational FE method well-defined for Maxwell.

## Trace operators
Each space has a "boundary trace" operator:
- H¹(Ω) → H^{1/2}(∂Ω):  trace.  (Function values on the boundary.)
- H(curl;Ω) → tangential trace γ_τ:  τ-component on ∂Ω.
- H(div;Ω) → normal trace γ_n:  n-component on ∂Ω.

These are continuous, surjective.  They are the discrete-form
versions of "Dirichlet BC" (H¹), "essential E_τ BC" (H(curl)),
"essential B_n BC" (H(div)).

## What the FEM does
Whitney elements W^k are **conforming finite-dimensional subspaces**:
- W^0 ⊂ H^1
- W^1 ⊂ H(curl)
- W^2 ⊂ H(div)
- W^3 ⊂ L²

The discrete sequence
    W^0 →^G W^1 →^R W^2 →^D W^3
mirrors the continuous one, and **commutes** with the inclusion.
This is the modern Arnold-Falk-Winther "cochain projection"
condition — necessary for stable approximation.
"""

MAXWELL_HOUSE = """
# Maxwell's house: the functional framework for Maxwell's equations
(Bossavit 1998 Fig. 5.1)

The full Maxwell equations live in a 4-pillar structure with two
"primal" columns (left, for e/a/φ) and two "dual" columns (right,
for h/d/j/ρ).  Each column is a copy of the de Rham complex.

    PRIMAL (inner orientation, "straight" forms)
    ─────────────────────────────────────────────
    ψ ∈ H¹      ←  d  ←   φ ∈ H¹        scalar potentials
    │ grad                │
    ▼                     ▼
    e ∈ H(curl) ←  d  ←   a ∈ H(curl)  vector potential (e = -∂_t a - grad φ)
    │ rot                 │ rot
    ▼                     ▼
    b ∈ H(div)            b ∈ H(div)   magnetic flux (b = rot a)
    │ div                 │ div
    ▼                     ▼
    0                     0

    DUAL (outer orientation, "twisted" forms)
    ─────────────────────────────────────────────
    h ∈ H(curl) ←  d  ←   T ∈ H(curl)  (Bossavit's electric vector potential)
                                       (j = rot T, electric vector potential)
    ↑ rot                 ↑ rot
    d, j ∈ H(div) ←  d  ←
    ↑ div                 ↑ div
    ρ, q ∈ L²

The horizontal arrows between PRIMAL and DUAL columns are the
**constitutive laws**:

    d  = ε e         (Hodge ⋆ scaled by ε)
    b  = μ h         (Hodge ⋆ scaled by μ)
    j  = σ e         (Ohm's law: Hodge ⋆ scaled by σ; only in conductors)

Time derivative ∂_t connects layers vertically (Faraday and Ampère):
    ∂_t b + rot e = 0          ← reading down the primal column
    −∂_t d + rot h = j         ← reading up the dual column.

## Why two columns?
Maxwell's equations have a beautiful symmetry between (e, b) and
(h, d).  The primal column hosts the variables which integrate over
inner-oriented domains (lines, surfaces); the dual hosts those that
integrate over outer-oriented (crossing) domains.  Ohm's law and the
constitutive ⋆ are the only places this symmetry is broken (because
they involve material parameters).

## Discrete realization
Replace each pillar by W^0 → W^1 → W^2 → W^3 on a primal simplicial
mesh and (its dual on the dual mesh, or via mass matrices).  Then
Maxwell's equations become finite linear systems — and DR=0, RG=0
hold EXACTLY, no spurious eigenmodes, no charge conservation
violations.

This is why the modern view (Bossavit, Hiptmair, Arnold-Falk-Winther)
calls Whitney elements "structure-preserving" or "compatible"
discretizations.

## Three vs four columns
Some authors (Tonti) draw a 4-cell × 4-row chart with explicit
twisted/primal labels.  Bossavit draws Maxwell's house with the
constitutive ⋆ as a horizontal arrow.  Arnold's FEEC draws the
abstract Hodge-Laplacian complex.  All three are the same picture
in different notation.
"""


def get_de_rham_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "complex"        — de Rham complex, exact sequence, Poincaré lemma
      "weak"           — weak operators, H¹/H(curl)/H(div), traces
      "maxwell_house"  — Bossavit's 4-pillar diagram, constitutive arrows
    """
    topic = topic.lower().strip()
    if topic in ("complex", "de_rham", "deRham", "exact"):
        return DE_RHAM_COMPLEX
    if topic in ("weak", "weak_operators", "sobolev", "h_curl", "h_div"):
        return WEAK_OPERATORS
    if topic in ("maxwell_house", "house", "diagram", "pillars"):
        return MAXWELL_HOUSE
    if topic == "all":
        return "\n\n".join([DE_RHAM_COMPLEX, WEAK_OPERATORS, MAXWELL_HOUSE])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, complex, weak, maxwell_house."
    )
