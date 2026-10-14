"""
Whitney elements (Whitney forms): the discrete realization of the de Rham
complex on a simplicial mesh.

Distilled from:
  - Bossavit 1998 Ch.5 §5.2 (Whitney complex, w_n / w_e / w_f / w_T,
    continuity, DR=0 and RG=0, Betti numbers)
  - Kameari 2011 §4-5 (Whitney's elementary forms φ_i dφ_j − φ_j dφ_i,
    higher-order edge elements, hexahedral Serendipity 2nd order,
    Kameari's symmetric 2nd-order element)
  - Arnold-Falk-Winther 2006 §5 (modern FEEC view: Whitney forms =
    P_1^- Λ^k(T_h), the lowest-order subspace)
  - Whitney 1957, "Geometric Integration Theory", IV.27
"""

WHITNEY_OVERVIEW = """
# Whitney elements: discrete differential forms on a simplicial mesh

## The four basis families (3-D tetrahedral mesh)
For a simplicial mesh with sets N (nodes), E (edges), F (faces),
T (tetrahedra), and barycentric coordinates  w_n  (the standard "hat"
function for node n):

| Degree | Name | Formula | Continuity | DoF |
|--------|------|---------|------------|-----|
| 0 | Node element  w_n | the hat function | C^0 (continuous) | value at n |
| 1 | Edge element  w_e | w_e = w_m ∇w_n − w_n ∇w_m   for e={m,n} | tangential only | circulation along e |
| 2 | Face element  w_f | w_f = 2(w_l ∇w_m × ∇w_n + w_m ∇w_n × ∇w_l + w_n ∇w_l × ∇w_m) for f={l,m,n} | normal only | flux through f |
| 3 | Cell element  w_T | 1/vol(T) on T, 0 elsewhere | discontinuous | integral over T |

## Key duality (defining property)
- ∫_n  w_{n'}    = δ_{n,n'}    (value at the dual node)
- ∫_e  τ·w_{e'}  = δ_{e,e'}    (circulation along the dual edge)
- ∫_f  n·w_{f'}  = δ_{f,f'}    (flux through the dual face)
- ∫_T  w_{T'}    = δ_{T,T'}    (integral over the dual cell)

So **the degrees of freedom are integrals of the field over the mesh
entity of matching dimension** — exactly the kind of degrees of
freedom that physical measurement gives you (voltmeter reads a line
integral; flux probe reads a surface integral; etc.).

## Inclusion in Sobolev spaces
- W^0 ⊂ H^1(D)      ( = L^2_grad(D) in Bossavit's notation)
- W^1 ⊂ H(curl; D)  ( = IL^2_rot(D))
- W^2 ⊂ H(div; D)   ( = IL^2_div(D))
- W^3 ⊂ L^2(D)

So Whitney elements are the **conforming** finite-element subspaces
of these natural Sobolev spaces.

## The Whitney complex
With the incidence matrices G, R, D (node→edge, edge→face, face→cell):

       grad             rot              div
   W^0 ───→ W^1  ───→  W^2  ───→  W^3
   ↕         ↕          ↕          ↕
   W₀ ───→  W₁  ───→  W₂   ───→   W₃
       G                R                D

Top row: spaces of Whitney forms (functional spaces).
Bottom row: DoF vector spaces (combinatorial / discrete).
Vertical arrows: bijections (each form ↔ its DoF tuple).

Exact relations (this is the heart of the construction):
- grad(W^0) ⊂ W^1   ↔   image of G is in DoF-edges
- rot(W^1)  ⊂ W^2   ↔   image of R is in DoF-faces
- div(W^2)  ⊂ W^3   ↔   image of D is in DoF-cells

And dually:
- R G = 0   (discrete  rot grad = 0)
- D R = 0   (discrete  div rot  = 0)

These are EXACT (not approximate).  No numerical error introduced.
"""

EDGE_ELEMENTS = """
# Edge elements (Whitney 1-forms) — the workhorse for H(curl) FEM

## Why edge DoFs, not nodal?
Classical nodal FEM expresses a vector field a = (a_x, a_y, a_z) by
three independent scalar interpolants — one Cartesian component each.
But the curl-curl operator mixes components, and across material
interfaces the **normal** component of a may jump while the **tangential**
part stays continuous.  Cartesian-component interpolation forces ALL
three components to be continuous, which is wrong physics.

Edge elements fix this:
- Only tangential continuity is enforced (matches the physics of E
  and H across material interfaces).
- Normal jump is naturally allowed (matches discontinuity of normal
  H or normal E_n at material interface).

## Whitney's original formula (1957)
For an edge from node i to node j with barycentric coordinates
φ_i, φ_j on a simplex:

    w_{ij} = φ_i d φ_j − φ_j d φ_i.

The 3-D vector-calculus form is what Bossavit and Nedelec use:

    w_e = w_m ∇w_n − w_n ∇w_m   (edge e from m to n)

Per-tetrahedron analytic shape:  w_e(x) = α × x + β
(α a fixed vector parallel to the opposite edge of the tetrahedron,
β a fixed vector).  This  α × x + β  form is Nedelec's "mixed FEM
in R^3" (1980).

## Degrees of freedom: circulation
For h ∈ W^1, the DoF on edge e is the **edge circulation**

    h_e := ∫_e τ · h        (line integral, tangent vector τ).

The interpolation w_e is designed so that
    ∫_{e'} τ · w_e = δ_{e,e'},
so h = ∑ h_e w_e is uniquely determined by its tangential-integral
DoFs.

## How rot maps to the incidence matrix R
If h ∈ W^1 and j = rot h ∈ W^2, the face flux of j on face f is

    j_f = ∑_e R_{f,e} h_e

where R is the face-edge incidence matrix (sign ±1 if e is on
boundary of f, 0 otherwise).  This is **exact** by Stokes' theorem:
∫_f n·rot h = ∮_∂f τ·h, and ∂f decomposes into a signed sum of edges.

## Tree-cotree gauge
For solving curl-curl problems  rot(ν rot a) = j  on a contractible
domain, the kernel of rot in W^1 is exactly  grad W^0  (the gauge
freedom a → a + grad φ).  To fix the gauge, pick a **spanning tree**
of edges T ⊂ E: a maximal subset with no cycles.  Set the DoFs
a_e = 0 for all e ∈ T.  The remaining "cotree" edges carry the true
degrees of freedom — they are in bijection with face fluxes of
b = rot a.

In non-contractible domains (loops, holes), the spanning tree must
be augmented with one cotree edge per loop (a "belt fastener"), giving
a "belted tree".  The number of belts = Betti number b_1.

## Curved high-order edge elements (Nedelec hierarchy)
The Nedelec 1980 paper introduces tetrahedral edge elements of
arbitrary polynomial degree k.  Number of DoFs:

    N_k = k(k+2)(k+3)/2

- k=1: 6 (Whitney): one DoF per edge, six edges per tet.
- k=2: 20.  Split: 2 per edge (12) + 2 per face (8).  Several
  different bases are equivalent — see Kameari 2011 §5 for
  Kameari's symmetric variant.
- k=3: 45.  And so on.

The Arnold-Falk-Winther 2006 framework calls these P_k^- Λ^1 and
P_k Λ^1 (two families) — see `feec`.
"""

FACE_ELEMENTS = """
# Face elements (Whitney 2-forms) — for H(div) FEM

## Definition
For a face f = {l, m, n} on a tetrahedron:

    w_f = 2 (w_l ∇w_m × ∇w_n + w_m ∇w_n × ∇w_l + w_n ∇w_l × ∇w_m)

This is the 3-D analogue of the Raviart-Thomas 2-D element (1977),
discovered independently for 3-D by Nedelec (1980) and by Schaubert-
Wilton-Glisson (1984) for scattering problems (RWG functions).

Per-tetrahedron shape:  w_f(x) = α x + β  with scalar α ∈ R and
vector β ∈ R^3.

## Continuity: normal across faces
- Normal component  n·w_f  is continuous across the face f.
- Tangential component can jump.

This matches the physics of magnetic flux density b (normal continuous,
tangential jumps where H jumps) or current density j.

## Degrees of freedom: face flux
For b ∈ W^2 and face f:

    b_f := ∫_f n · b        (signed flux through f).

div(W^2) is captured by the cell-face incidence matrix D:
    (div b)_T = ∑_f D_{T,f} b_f      (Gauss theorem, exactly).

## Use cases in EM
- Vector potential a in the H(div) formulation (b = a, treat as
  primal).
- Magnetic flux density b (when b is the primary unknown).
- Current density j in mixed-method PEEC.
- Hybrid H-A formulations where eddy currents and air share interfaces.

## Connection to the dual mesh
Whitney 1-forms on the **dual** mesh are isomorphic to twisted
2-forms on the primal mesh, by the duality of Kameari/Tonti.
This is why "edge elements for H" and "face elements for D" are
mirror images on the two sides of the Maxwell house diagram.

## Mass matrix (Hodge ⋆ in the discrete setting)
M^{(2)}_{f,f'} = ∫_D α(x) w_f(x) · w_{f'}(x) dx

where α = 1/μ for the magnetic energy term, etc.  This is the
**only** place where the material coefficient (the Hodge star) enters
the discrete system — D, R, G are entirely metric-free.
"""

PROPERTIES_PROOFS = """
# Key combinatorial-topological properties of the Whitney complex

## Exactness modulo topology
For a **contractible** meshed domain (simply connected, with connected
boundary):

    grad(W^0) = ker(rot ; W^1)     ←  curl-free 1-forms are gradients
    rot(W^1)  = ker(div ; W^2)     ←  divergence-free 2-forms are curls

These are EXACT identities on the discrete Whitney complex,
**proven combinatorially** from the corresponding continuous fact
(Poincaré lemma).  No mesh-refinement or limit is needed.

## Failure of exactness = topology
If the meshed domain has b_1 loops (e.g., a torus), then
    dim[ker(rot ; W^1) / grad(W^0)] = b_1.
Each "missing gradient" corresponds to a non-bounding 1-cycle
(a loop that goes around a hole).

If the meshed domain has b_2 holes (e.g., a hollow sphere shell),
then
    dim[ker(div ; W^2) / rot(W^1)] = b_2.
Each "missing curl" corresponds to a non-bounding 2-cycle (a closed
surface around a hole).

## Euler-Poincaré identity
    χ(D) = N − E + F − T = b_0 − b_1 + b_2 − b_3

with b_0 = number of connected components, b_3 = 0 for bounded D.

For:
- contractible D: χ = 1
- simple torus:   χ = 0  (b_0=1, b_1=1, b_2=0)
- hollow sphere: χ = 2  (b_0=1, b_1=0, b_2=1)

## Mass matrices
For each degree k, define M^{(k)}(α) with entries

    M^{(k)}(α)_{ss'} = ∫_D α(x) w_s(x) · w_{s'}(x) dx       (k=1,2)
                    = ∫_D α(x) w_s(x)   w_{s'}(x) dx        (k=0,3)

α is the material property (μ, ε, σ, or their inverses).  Anisotropic
α (a 3×3 tensor) gives anisotropic media — supported by Whitney
elements without any modification of the basis.

## The discrete Hodge star
Hodge star ⋆ : Λ^k → Λ^{n−k} discretizes as the matrix combination

    ⋆_h^{(k)} := M^{(n−k)}(1/α) · [conversion between primal and dual]

There is **no canonical discrete Hodge** — different choices give
different but consistent schemes (FDTD's Yee grid, FIT, mass-matrix
Galerkin, etc.).  The choice affects accuracy and conditioning, not
exactness of the conservation laws.

## Why DR=0 and RG=0 matter for FEM
- div(rot a) = 0 is automatic in the **discrete** problem.  No need
  to impose div = 0 as a constraint.
- The kernel of curl-curl is exactly the gradient image: classical
  gauge fixing applies.
- Spurious modes (Chapter 9 of Bossavit) are eliminated as long as
  the FEM space is a true sub-complex (grad W^0 ⊂ W^1, rot W^1 ⊂ W^2).
  Mixing element types (e.g., nodal H and edge E) breaks this and
  produces spurious modes.
"""


def get_whitney_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "overview"        — four basis families, complex diagram, DoF nature
      "edge"            — edge (W^1) elements: formula, gauge, Nedelec hierarchy
      "face"            — face (W^2) elements: formula, RWG, mass matrix
      "properties"      — DR=0/RG=0, Betti numbers, Euler-Poincaré, mass matrices
    """
    topic = topic.lower().strip()
    if topic == "overview":
        return WHITNEY_OVERVIEW
    if topic in ("edge", "edge_elements", "1form", "1-form", "w_e"):
        return EDGE_ELEMENTS
    if topic in ("face", "face_elements", "2form", "2-form", "w_f"):
        return FACE_ELEMENTS
    if topic in ("properties", "exactness", "betti"):
        return PROPERTIES_PROOFS
    if topic == "all":
        return "\n\n".join([
            WHITNEY_OVERVIEW,
            EDGE_ELEMENTS,
            FACE_ELEMENTS,
            PROPERTIES_PROOFS,
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, overview, edge, face, properties."
    )
