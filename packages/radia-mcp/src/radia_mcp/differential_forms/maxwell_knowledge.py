"""
Maxwell's equations in differential-form language: which physical
quantity is which type of form, the Hodge star operator, twisted vs
straight forms, and the Tonti / Kameari diagram.

Distilled from:
  - 微分形式.pdf §6-7 (Maxwell's equations as 1-/2-/3-forms, Hodge in
    Minkowski space, constitutive relations via Hodge)
  - Bossavit 1998 Ch.1 §1.1 (field equations, units)
  - Kameari 2011, fig. 6 (the canonical diagram: 0-/1-/2-/3-form
    on primal and dual mesh, straight vs twisted)
  - Arnold-Falk-Winther 2006, Table 2.1 (vector calculus ↔ exterior
    algebra correspondences in R^3)
"""

MAXWELL_AS_FORMS = """
# Maxwell's equations as differential forms

## Field type by integration domain (3+1-D split)

| Quantity | Symbol | Type   | Integrates over | Units |
|----------|--------|--------|-----------------|-------|
| Electric potential        | φ        | 0-form  | points     | V        |
| Electric field            | e        | 1-form  | lines      | V/m      |
| Magnetic vector potential | a        | 1-form  | lines      | V·s/m    |
| Magnetic field intensity  | h        | 1-form (twisted) | lines | A/m   |
| Magnetic flux density     | b        | 2-form  | surfaces   | T = Wb/m^2 |
| Electric flux density     | d        | 2-form (twisted) | surfaces | C/m^2 |
| Current density           | j        | 2-form (twisted) | surfaces | A/m^2 |
| Charge density            | ρ        | 3-form (twisted) | volumes  | C/m^3 |

(Energy density, action, Lagrangian densities are also forms — 3-form
or 4-form — but we focus on 3-D quasi-statics here.)

## Why the integration domain
- "Voltage between two points" = line integral of e → 1-form
- "Magnetic flux through a surface" = surface integral of b → 2-form
- "Charge in a volume" = volume integral of ρ → 3-form

Vector-calculus identification:
- 1-form  ω = E_1 dx + E_2 dy + E_3 dz  ↔  vector E = (E_1, E_2, E_3)
- 2-form  ω = B_1 dy∧dz + B_2 dz∧dx + B_3 dx∧dy  ↔  vector B
- 3-form  ω = ρ dx∧dy∧dz  ↔  scalar ρ

## Faraday's law:  ∂_t b + d e = 0
    1-form e → 2-form de  (this is rot E)
    Then de = −∂_t b, the differential form  Faraday's law.

## Ampère's law:  d h = j + ∂_t d
    1-form h → 2-form dh  (this is rot H)
    Equals current density 2-form  j  plus displacement 2-form ∂_t d.

## Continuity: d j + ∂_t ρ = 0
    2-form j → 3-form dj  (this is div J)
    Charge conservation = differential identity dd = 0 applied to (dh).

## 4-D covariant form (special relativity)
Combine e and b into a 4-D 2-form  F = b − (1/c) e ∧ dt :

    F = − (E_x/c) dx⁰∧dx¹ − (E_y/c) dx⁰∧dx² − (E_z/c) dx⁰∧dx³
         + B_x dx²∧dx³ + B_y dx³∧dx¹ + B_z dx¹∧dx²

Combine d, h into G:
    G = h ∧ dt − (1/c) d  (with appropriate signs).

Then the two pairs of Maxwell equations become **two** 4-D form equations:

    dF = 0         (homogeneous: Faraday + ∇·B = 0)
    dG = J         (inhomogeneous: Ampère + Gauss)

Plus the constitutive law  G = ⋆ F  (Hodge star in Minkowski space).
This is what 微分形式.pdf §6 derives explicitly.
"""

HODGE_STAR = """
# Hodge star ⋆ — the metric in disguise

## Definition (in an oriented inner-product space)
Let V be an n-dimensional oriented inner-product space.  For a
k-form ω there is a unique (n−k)-form ⋆ω satisfying

    ω ∧ μ = ⟨⋆ω, μ⟩ vol      for all (n−k)-forms μ,

where vol is the volume form and ⟨·,·⟩ the induced inner product on
forms.

## Properties
- Linear isometry  Λ^k → Λ^{n−k}.
- Double Hodge:  ⋆⋆ω = (−1)^{k(n−k)} ω
   In E^3:  ⋆⋆ = id  (no sign);  in Minkowski R^{1,3}: ⋆⋆ = ±id
   depending on signature.
- Reverses inner/outer orientation (straight ↔ twisted forms).

## Explicit action on basis (R^3, Euclidean metric)
- ⋆1     = dx ∧ dy ∧ dz       (0-form → 3-form)
- ⋆dx    = dy ∧ dz            (1-form → 2-form)
- ⋆(dy ∧ dz) = dx             (2-form → 1-form)
- ⋆(dx ∧ dy ∧ dz) = 1         (3-form → 0-form)

## Constitutive relations via ⋆
In a linear isotropic medium:

    d = ε e                  ↔   d = ε ⋆ e   (1-form → 2-form via ⋆)
    b = μ h                  ↔   b = μ ⋆ h   (1-form → 2-form via ⋆)

In a **conducting** medium (Ohm's law):
    j = σ e                  ↔   j = σ ⋆ e   (1-form → 2-form via ⋆)

The constitutive relations are **the place** where geometry (the
metric) enters Maxwell's theory.  The de Rham operator d itself is
**metric-free** (purely topological).

## In Minkowski 4-D
⋆ takes 2-forms to 2-forms.  Then  G = ⋆ F  becomes:

    constitutive relation = single 4-D form equation,
    G   = (1/μ₀)(⋆F)  in vacuum.

## Why FEM cares
The Hodge star is the source of the **mass matrix**
   M^{(k)}_{ij} = ∫_D ⟨α w_i, w_j⟩  (basis-mass for k-forms)
in the Whitney-form FEM (see `whitney`).  The metric/material
properties of the problem enter **only** through these mass matrices.
The "topological" matrices G, R, D (grad/rot/div incidence) are
metric-free.
"""

TWISTED_VS_STRAIGHT = """
# Inner vs outer orientation (straight vs twisted forms)

## Inner orientation
Intrinsic to the manifold.  Examples:
- For an edge: a "forward" direction (arrow along the edge).
- For a face: a sense of rotation (clockwise vs counter-clockwise
  when viewed from one side).
- For a volume: an ordering of vertices that defines "positive" cell.

## Outer (external) orientation
Requires an embedding in a larger space:
- For an edge: a "crossing direction" — which way charge flows through
  the edge.
- For a face: a normal vector — which side is "outside".
- For a volume: a sign  +/−  (whether the cell is on the "positive" or
  "negative" side of some reference).

In an oriented ambient space (E^3 with its standard right-hand rule),
**outer orientation of a k-cell = inner orientation of its complement**
(an (n−k)-cell).  So in E^3 the two notions coincide given the
ambient orientation.

## Straight forms vs twisted forms (densities)
- **Straight forms** depend on inner orientation.  Reversing the
  orientation of the integration domain flips the sign of the integral.
  Examples in EM:  e (1-form), a (1-form), b (2-form).
- **Twisted forms** (also called **densities**, or **pseudo-forms**)
  depend on outer orientation.  Their integration over a domain
  needs a crossing-direction; under reversal of inner orientation the
  sign does NOT flip.  Examples in EM:  h (1-form twisted),
  j (2-form twisted), d (2-form twisted), ρ (3-form twisted).

## Why pseudovectors are awkward — and forms are clean
Vector calculus distinguishes "polar vectors" (E, A) from "axial
vectors" (H, B, J) — the latter flip sign under coordinate inversion.
This is the same distinction as straight vs twisted forms, but
hidden in coordinate machinery.  In differential-form language each
quantity has its OWN type and the difference is **structural**, not
a convention.

## Kameari diagram (canonical layout, see fig. 6 of SA-11-013)

Straight (primal mesh)             Twisted (dual mesh)
─────────────────────                ────────────────────
0-form  φ  (nodes)                   3̃-form  ρ        (dual cells)
1-form  e  (edges)        ⋆          2̃-form  d, j     (dual faces)
2-form  b  (faces)        ⋆          1̃-form  h, T     (dual edges)
3-form  Ω* (cells)                   0̃-form  --       (dual nodes)

   ↓ grad                                  ↓ div
   ↓ rot                                   ↓ rot
   ↓ div                                   ↓ grad

Left column: d acts as grad → rot → div (top-to-bottom).
Right column: d acts going UP (twisted forms run on the dual mesh).
The horizontal arrows linking them are the constitutive ⋆ (material
properties ε, μ, σ).

This diagram is the same as the FDTD "primal-dual lattice" picture
and as Tonti's "classification diagram".

## Practical consequence for FEM
Choosing the **right kind** of finite element for each variable:
- e, a  → edge elements (Whitney 1-forms)  on primal mesh
- b     → face elements (Whitney 2-forms)  on primal mesh
- h, j  → edge/face elements on the **dual** mesh — OR equivalently,
  use the constitutive Hodge ⋆ to map them back to primal-mesh
  elements via the mass matrix.

The "Maxwell house" of Bossavit (see `de_rham`) makes this
two-column structure of primal and dual quantities explicit.
"""


def get_maxwell_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "maxwell"             — equations as forms, 4-D form, F and G
      "hodge"               — ⋆ operator, constitutive laws, metric role
      "twisted"             — twisted vs straight forms, Kameari diagram
    """
    topic = topic.lower().strip()
    if topic in ("maxwell", "equations", "field_eq"):
        return MAXWELL_AS_FORMS
    if topic in ("hodge", "hodge_star", "constitutive"):
        return HODGE_STAR
    if topic in ("twisted", "twisted_forms", "orientation", "kameari"):
        return TWISTED_VS_STRAIGHT
    if topic == "all":
        return "\n\n".join([MAXWELL_AS_FORMS, HODGE_STAR, TWISTED_VS_STRAIGHT])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, maxwell, hodge, twisted."
    )
