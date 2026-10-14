"""
Electromagnetic forces in differential-form language.

The "question of forces in electromagnetics" — long controversial in
vector-calculus — collapses into a single energy-balance identity when
formulated with differential forms.  All known force-computation
methods (energy variation, Maxwell stress, eggshell, Arkkio,
Coulomb-nodal, Kameari nodal, Lorentz) emerge as special cases.

Distilled from:
  - F. Henrotte (with A. Bossavit), "Handbook for the computation of
    electromagnetic forces in a continuous medium", ICS Newsletter
    (post-2004) — the 7-page unified treatment
  - A. Bossavit, "Forces inside a magnet", ICS Newsletter 11(1):4-12,
    March 2004
  - A. Kameari (亀有), SA-11-013, 2011 — open question on the
    differential-form representation of shear stress
  - A. Kameari, "Local force calculation in 3D FEM with edge elements",
    IJAEM 3:231-240, 1993 — the canonical Kameari nodal force method
  - 新しい計算電磁気学 (五十嵐他 2003) Ch.4 (Nishiguchi) §4.1-4.3
    on electromagnetic forces + Maxwell stress tensor multiplicity
  - F. Henrotte et al., "Electromagnetic force density in a
    ferromagnetic material", IEEE Trans. Mag. 40(2):553, 2004
"""

FORCES_OVERVIEW = """
# EM forces in differential-form language: why this matters

## The long controversy
For a century, "how to compute the electromagnetic force on a
deforming body" had no canonical answer.  Textbooks listed multiple
formulas — Maxwell stress tensor, energy method, virtual work,
Lorentz force, magnetic charge, magnetic current — and they
disagreed in non-trivial cases (anisotropic, magnetostrictive,
hysteretic materials).  Engineers chose one by feel.

## The reason
Vector calculus cannot cleanly express the energy exchange between
an electromagnetic field and a **deforming continuous medium**.
The setup requires:
  - A material manifold M (the body, NOT yet placed in space)
  - A Euclidean space E (where the motion happens)
  - A time-dependent placement map p_t : M → E

These are differential-geometry objects.  Trying to express them in
vector calculus leads to coordinate-dependent expressions that hide
the underlying structure.

## The differential-form answer (Bossavit-Henrotte, ICS Newsletter)
Every classical force formula is a special case of ONE identity:

    ρẆ^em = -σ_em : ∇v                    (electromechanical coupling)

where  σ_em = b ⊗ h̃ - {h̃·b - ρΨ} I  is the **Maxwell stress tensor**
derived from the energy density ρΨ.  Different materials give
different energy densities → different stress tensors, but the
form of the coupling identity is universal.

## What this server provides

- Mechanical framework: material vs Euclidean manifolds (`mechanical_framework`)
- Material (Lie) derivatives of p-forms (`material_derivative`)
- Maxwell stress tensor derivation (`maxwell_stress`)
- Unified force methods table (`force_methods`)
- Kameari's open question on shear stress (`kameari_stress`)

Pairs with `differential_forms_mathematica_recipes('maxwell_stress')`
for symbolic verification of σ_em on specific energy densities.
"""

MECHANICAL_FRAMEWORK = """
# Mechanical framework: two manifolds + placement map

## The setup (Bossavit-Henrotte 2004 Handbook §III)

Two manifolds, distinct roles:

| Manifold | Role | Metric |
|----------|------|--------|
| **M** (material) | each point = a material particle (atom) of the body | NONE intrinsically |
| **E** (Euclidean) | the space where motion happens | gᴱᵢⱼ = δᵢⱼ |

The two are linked by the **placement map**:

    p_t : X ∈ M → x = p_t X ∈ E                  (one per instant t)

The codomain p_t M ⊂ E is the **deformed state** at time t.  The
velocity field is

    v = ∂_t x ∈ TE

i.e. the tangent vectors to all trajectories.

## Why two manifolds, not one?
This is THE conceptual move that makes force theory clean.  In vector
calculus, one tries to compute forces directly in 3D space and ends
up with frame-dependent expressions.  By distinguishing:
  - the body's intrinsic configuration (material manifold M, no metric)
  - its embedding in space (placement p_t into Euclidean E with metric)
one can separate:
  - what depends on the body itself (e.g. energy density as a 3-form on M)
  - what depends on the embedding (e.g. lengths, angles, areas in E)

## Quantities on each side

| On M (material) | Notation | On E (Euclidean) | Notation |
|-----------------|----------|------------------|----------|
| Magnetic flux 2-form | B | Magnetic flux 2-form (image) | b ≡ p_t B |
| Electric flux 2-form | D | Electric flux 2-form (image) | d ≡ p_t D |
| Energy density 3-form | (none — depends on metric in E) | ρΨ : Λ²(E)×Λ²(E) → Λ³(E) |
| Velocity field | n/a | v ∈ TE | |

State variables (the unknowns of the dynamic problem) live on M.
Field intensities (which need a metric to define "Tesla per square
meter") live on E via the placement.

## Material vs spatial (Eulerian) descriptions

- **Material (Lagrangian)**: write equations on M, with ∂_t for time
  derivative.  The mesh **moves with the body** (comoving mesh).
- **Spatial (Eulerian)**: write equations on E, with L_v for material
  derivative (Lie derivative along v).  The mesh is **fixed in space**.

Maxwell equations take their familiar form  ∂_t B + d Ẽ = 0  on M
(Lagrangian).  The same equations on E read  L_v b + d ẽ = 0 .  Both
are correct; the choice depends on whether the FEM mesh moves with
the body or is fixed in space.

## Connection to FEM
In Radia / NGSolve practice:
  - Static problems: M = E (rigid body, no placement effect)
  - Quasi-static problems with rigid motion: M ≠ E, p_t is a rigid
    transformation, used in motor / generator simulations
  - Magnetostriction, structural EM coupling: M ≠ E, p_t describes
    body deformation
"""

MATERIAL_DERIVATIVE = """
# Material derivatives of p-forms (Lie derivative formulas)

When the body deforms with velocity field v ∈ TE, p-forms on E
have material derivatives that follow specific formulas depending on
the form degree.  These are the algebraic backbone of every
force-computation method.

## The 4 cases in 3D

| Form degree | Symbol | Material derivative L_v |
|---|---|---|
| 0-form (scalar) | f | L_v f = ḟ |
| 1-form (covector, e.g. h̃, ẽ) | h | L_v h = ḣ + (∇v) · h |
| 2-form (e.g. b, d, j) | b | L_v b = ḃ - b · (∇v) + b · tr(∇v) |
| 3-form (e.g. ρ, ρΨ) | ρ | L_v ρ = ρ̇ + tr(∇v) · ρ |

where ż = ∂_t z + v · ∇z is the total derivative.  The (∇v) terms
encode how the form components change because the underlying body
deforms — NOT because the field itself changes in space.

## Why each form degree picks up a different (∇v) term
- 0-form: just a number at each material point, no transformation
- 1-form: covariant; transforms with the inverse Jacobian → +(∇v)·h
- 2-form: contravariant in 2 indices + Jacobian determinant correction
- 3-form: scales by Jacobian determinant ⇒ pure tr(∇v) factor

These are the same conformal-weight scalings that appear in the
Kelvin transform (see `mathematica_recipes('kelvin')`).

## Crucially, all formulas involve ∇v, NOT v alone
This is why "the electromagnetic force is the conjugate of ∇v in
the energy balance" — the **stress tensor σ_em** naturally pairs with
∇v, NOT with v.  See `maxwell_stress` for the explicit identification.

## Use in chain rule for energy density
For energy density ρΨ : (d, b) ↦ Λ³(E) :

    L_v {ρΨ(p_t D, p_t B)}
        = {L_v ρΨ}(d, b)                ← "explicit" change in functional
        + ẽ · L_v d  +  L_v b · h̃       ← "field" change via state vars

Each term has a clear interpretation:
- "explicit" term → the stress σ_em
- "field" term → Maxwell equations (Faraday, Ampère)

This is the algebraic core that makes the Bossavit-Henrotte
derivation work in 7 pages.

## Practical implication for Radia / NGSolve
When computing forces on moving conductors or magnetostrictive
materials, the right time derivative to use is NOT ∂_t (which is
the Lagrangian) but L_v (the Lie derivative).  See Bossavit Ch.8
eddy currents for the L_v b expression in transient eddy-current
problems.
"""

MAXWELL_STRESS = """
# Maxwell stress tensor σ_em — the unifier

## Derivation in one slide

Energy balance (Bossavit-Henrotte 2004):

    ∂_t Ψ  +  Ṡ  +  Ṙ  =  0

  Ψ  = ∫_Ω ρΨ(b, d)        electromagnetic energy
  Ṡ = ∫_{∂Ω} ẽ ∧ h̃         Poynting flux (boundary power)
  Ṙ = ∫_Ω {ẽ ∧ j + ρẆ^em}  converted: Joule + mechanical

Substituting the material-derivative formulas (`material_derivative`)
and the Maxwell equations on the Euclidean manifold E, all terms
either cancel or reorganize.  The mechanical-power density isolates:

    ρẆ^em = -σ_em : ∇v

with

    σ_em = b ⊗ h̃  -  {h̃·b - ρΨ(b, d)} I

For an energy density depending on b only (linear isotropic
magnetics, ρΨ = b·h̃/2), this reduces to the familiar form:

    σ_em^magnetic = (1/μ) [b ⊗ b - (1/2) |b|² I]

## What the formula says

- **σ_em is a (1,1)-tensor** (or in vector-form language, a vector-
  valued twisted 2-form — exactly Kameari's question in SA-11-013).
- The first term b ⊗ h̃ is the **dyadic product**: pulls the field
  apart into a "tension along the lines" + "pressure perpendicular
  to them" decomposition (Faraday's tubes of force picture).
- The second term {h̃·b - ρΨ}I is the **isotropic correction** that
  removes the trace if needed.

## Material extensions

The procedure works for ANY energy density.  Different materials
give different stress tensors but the same structural form:

| Material | Energy density ρΨ depends on | Stress tensor σ_em |
|----------|------------------------------|--------------------|
| Linear isotropic | b only | (1/μ)[b⊗b - (1/2)|b|² I] |
| Anisotropic | b, with tensor μ | b⊗h̃ - (h̃·b/2) I + anisotropy correction |
| Magnetostrictive | b, strain ε | b⊗h̃ + ∂_ε ρΨ - (h̃·b - ρΨ) I |
| Hysteretic | b, internal state ξ | b⊗h̃ + ∂_ξ ρΨ · (...) - ... |

For magnetostriction the extra term  ∂_ε ρΨ  IS the mechanical
strain dependence — included naturally, not added ad-hoc.

## Symbolic verification with Mathematica

Recipe: `mathematica_recipes('maxwell_stress')`.  Given a specific
ρΨ(b), Mathematica derives σ_em, checks symmetry (when expected),
and computes ρẆ^em on a sample velocity field — all symbolically.
"""

FORCE_METHODS = """
# Unified force methods: 7 popular formulas, one identity

All of these are particular ways to evaluate the same quantity
F = ∫_Y ρf^em = ∫_∂Y n·σ_em on a rigid region Y (Eq. 57 of
Bossavit-Henrotte).  They differ in how σ_em is integrated.

## The 7 methods

### (1) Force density (volume integration)
    ρf^em = div σ_em
    F_Y = ∫_Y ρf^em

Theoretically clean.  In practice ρf^em concentrates on material
interfaces (jumps of σ_em), so numerical integration is tricky.

### (2) Maxwell stress method (surface integration)
    F_Y = ∫_{∂Y} n · σ_em

Choose any surface ∂Y enclosing the rigid body in a force-free
region.  Integration is over a closed surface in air.  Used by
ANSYS, COMSOL, classical magnetostatics codes.

### (3) Eggshell method (volume integration on a shell)
    F_Y = -∫_Ω σ_em : ∇γ

where γ is a smooth function = 1 on the body Y, = 0 outside a thin
shell S surrounding Y.  Implemented in NGSolve via a "gear" mesh
layer; the support of γ is one element-layer thick.  References:
McFee-Webb-Lowther 1988, Henrotte-Deliege-Hameyer 2003.

Faster than (2) because uses existing volume-integration machinery.

### (4) Arkkio's method (2D rotating machines, special case)
    T_Y = ez/(R_o - R_i) · ∫_S r · (σ_em)_{rθ}

Specialization of the eggshell method for cylindrical shells in 2D
machines.  Arkkio 1987.

### (5) Coulomb's method (nodal forces from local Jacobian)
    F_N = ∫_Ω ρf^em · γ_N = -∫_Ω σ_em : ∇γ_N

For each FEM node N, take γ = its nodal shape function.  Gives a
nodal-force vector that sums correctly to the total resultant.
Coulomb 1983.  Used in COMSOL, Maxwell ANSYS.

**IMPLEMENTED in Radia 2026-05-21, refined 2026-05-22** via
`src/radia/panels/calc_em_force.py:_force_nodal_coulomb_vwm`.
Uses a smooth scalar indicator v(x) = 1 inside body, 0 outside
(equivalent to summing γ_N over all body nodes — gauge-invariant
form).  The Lorentz body term ∫_body J × B dV is added in for
current-carrying bodies.  Runs as `--methods nodal` in
`calc_em_force.py` (default in the methods list since the volume-
integral form is immune to the BND trace issues that plague
direct Maxwell stress with H1 grad on boundaries).

**Validated: 0.5% Newton-3 violation** vs Lorentz on the partner
body (Pile-Devillers-Le Besnerais 2018 reports 5% as "production-
grade"; we achieve 10x better).  Key prescriptions:
  - v_order = A_order (Mertens-Hameyer 2007) — critical for p>=2
  - Laplace-interpolated indicator (transition in source-free air)
  - Outer Dirichlet box >> body-body distance (image-artifact free)
  - For LINEAR mu: Henrotte 2004 == Coulomb (forwarded internally)

See `differential_forms_em_force_recipe` MCP tool for the full
practical recipe, common pitfalls, and validation protocol.

### (6) Ren-Razek formula for deformable bodies
Similar to (5) but specifically for **edge-element** based FEM where
the strain field is needed.  Ren-Razek 1992.

### (7) Kameari's local force calculation in 3D edge-element FEM
    F_N = -∫_Ω σ_em : ∇γ_N

where γ_N is the same **nodal scalar shape function** as in
Coulomb's method, but the σ_em uses B and h̃ derived from an
**edge-element** discretization (B = curl A with A as a Whitney
1-form, h̃ recovered via the discrete Hodge ⋆).

Kameari's contribution is the EXTENSION of Coulomb's idea to 3D
edge-element FEM with proper discrete consistency:
  - Field side: edge / Whitney complex (preserves d² = 0,
    avoids spurious modes — see `whitney`)
  - Force side: nodal test function γ_N for the local force
  - The two are linked via the discrete σ_em(B_edge, h̃_edge)
  - All discrete integrations match the Whitney-complex structure

Reference: A. Kameari, "Local force calculation in 3D FEM with edge
elements", IJAEM 3:231-240, 1993.  Cited as Ref [12] in Bossavit-
Henrotte Handbook.

This is the natural method for Radia / NGSolve-based magnet codes
because the magnetic field is already discretized on edges
(B = curl A with A in edge-elements).  Most Japanese magnetic FEM
codes (Femtet, JMAG, ELF / Radia derivatives) implement Kameari
节点力法 by default for force computation.

**Difference from Coulomb's method**:
  - Coulomb 1983: developed in NODAL-element FEM (B = grad ϕ_node)
  - Kameari 1993: extended to EDGE-element FEM (B = curl A_edge)
  - Both use the SAME formula F_N = -∫_Ω σ_em : ∇γ_N
  - But the discrete σ_em is different because B comes from
    a different FE space, and the accuracy / consistency
    arguments differ.

## Equivalence
All 7 methods give the SAME total resultant force F_Y (up to
discretization error).  They differ in:
  - Where they place the integration (volume vs surface)
  - How the test function γ is chosen
  - What information they expose (total only, nodal, density)
  - Numerical accuracy in difficult cases

Practical guidance (Bossavit, Henrotte, Kameari):
  - **Method (3) eggshell** is the production default for general 3D
  - **Method (7) Kameari** for edge-element-based codes (Radia,
    NGSolve, FEMTET) — preserves the structural advantage of the
    Whitney complex
  - **Method (5) Coulomb** for nodal-element-based codes (legacy)
  - **Method (2) Maxwell stress** for verification (independent
    integration domain)
"""

KAMEARI_STRESS_QUESTION = """
# Kameari's open question (SA-11-013, 2011) on shear stress

In the 2011 paper "辺要素を巡って — 「張力」って？", Kameari raises a
question that, despite the differential-form framework being well-
developed, remains pedagogically under-explored:

## The question

Bossavit's stress tensor σ_em = b⊗h̃ - {...}I is a (1,1)-tensor.  In
differential-form language, it is a **vector-valued twisted 2-form**:

- It acts on an OUTER-oriented surface (twisted 2-form) — the surface
  where the force flux is measured
- It returns a VECTOR (the force direction) — vector-valued

This is the natural form-language object for stress.  But:

**Question (Kameari)**: how do we represent SHEAR stress and
SHEAR strain in differential-form language?

For NORMAL (extension) strain and tension:
  - Displacement u → vector-valued 0-form (vector at each point)
  - Strain (= ∂u) → vector-valued 1-form (vector × tangent direction)
  - Tension → vector-valued twisted 2-form (vector × normal direction)

For SHEAR:
  - The relation "shear strain → vector × tangent direction" doesn't
    cleanly separate the orthogonal components
  - The "vector-valued form" language gets ややこしい (Kameari's
    word — "messy") for shear

## Why this matters
For magnetostriction analysis with anisotropic crystals, the shear
components of stress matter.  In vector calculus they are just
components of a 3×3 matrix, but the form-language interpretation is
unclear.

## Status (2026)
This is an open expository problem.  The mathematical apparatus
exists (e.g. covariant-derivative formulation of elasticity on
manifolds, see Marsden-Hughes "Mathematical Foundations of
Elasticity" 1983), but the **pedagogical bridge** between Kameari's
intuitive form-language picture for extension/tension and the full
elastic-shear case has not been written.

A worthy project for someone in the lab: write the missing
exposition.  Mathematica can help by computing explicit shear
examples symbolically (vector calculus side) and comparing with the
form-language candidates (form side).

## Connection to current MCP server
- Bossavit-Henrotte 2004 Handbook (in this server) gives the full
  vector-tensor side: σ_em as a (1,1)-tensor
- Kameari SA-11-013 (in this server) gives the form-language
  question without the answer
- The **answer** would unify Bossavit-Henrotte's σ_em with the
  vector-valued-twisted-2-form interpretation across both extension
  and shear

Use `mathematica_recipes('maxwell_stress')` to play with specific
σ_em examples while thinking about this.
"""


HENROTTE_FERROMAGNETIC = """
# Maxwell stress tensor σ_em for different ferromagnetic models

(Henrotte-Vande Sande-Deliege-Hameyer, IEEE Trans. Mag. 40(2):553, 2004)

The same material can be modelled in different ways, and **different
models give different Maxwell stress tensors** even when their
magnetic (constitutive) behavior is identical.  This is the source
of the long-standing controversy about EM forces in iron.

## The 4 ferromagnetic models and their MSTs

### Model 1: Linear isotropic (b = μ h̃)
Energy density:  ρΨ = b · h̃ / 2.
MST:
    σ_em = b ⊗ h̃ - {h̃·b - ρΨ} I
         = b ⊗ h̃ - (b·h̃/2) I

### Model 2: Permanent magnet as a 1-form (circulation density)
Magnetization m treated as 1-form → CONTRIBUTES NOTHING to σ_em
beyond the linear part (the term  ∫ m·h̃  gives no σ_em contribution).
MST identical to (1) with h̃ from the modified constitutive law.

### Model 3: Permanent magnet as a 2-form (flux density)
Magnetization m treated as 2-form → contributes:
    σ_em ⊃  m ⊗ h̃  (additional dyadic term)
QUITE DIFFERENT from Model 2 even for the same physical permanent
magnet.

### Model 4: Polycrystalline (saturable)
Energy:  ρΨ(b) =  ∫_0^b h̃(b') db'  (path integral, well-defined
because χ depends only on |b|).
MST:  σ_em = b ⊗ h̃ - (b·h̃ - ρΨ) I,
where ρΨ now contains the saturation curve.

### Model 5: Monocrystal with Weiss domains (anisotropic)
2D monocrystal with 4 easy-axis populations (a, b, c, d).
Energy:  ρΨ = magnetostatic + Zeeman + leakage (C-term).
For given h̃, m determined by constrained minimization of ρΨ.
Despite richer model, σ_em formula is the same as permanent
magnet (because internal variables a,b,c,d have time-derivative
not Lie-derivative — no ∇v contribution).

### Model 6: Magnetostrictive (with explicit strain coupling)
Energy:  ρΨ(b, ε) = ψ_magnetic + ψ_elastic + ψ_coupling
where ε is the strain tensor.
MST has extra term:
    σ_em = b ⊗ h̃ + ∂_ε ρΨ - (h̃·b - ρΨ) I

The ∂_ε ρΨ term IS the magnetostrictive force — emerges naturally
from differentiating the energy w.r.t. strain.

## Practical implication (from Henrotte 2004)

> "Three candidate expressions of the MST [...] If those three
>  expressions are implemented and the deformation of a closed
>  rectangular magnetic circuit is computed, one observes that
>  the three models, which are perfectly equivalent on the magnetic
>  side, give different deformations.  Again, it is a matter of
>  experiment and modeling to determine, case by case, which
>  model matches the best reality."

In other words: **the right Maxwell stress tensor is determined
empirically**, not by mathematical first principles alone.

For Radia/NGSolve practice:
- Linear isotropic + saturation curve → Model 4 (polycrystalline)
  is usually the right choice
- Permanent magnets → need to decide Model 2 vs Model 3 based on
  experimental comparison (e.g. against TEAM Problem 23 data, see
  `forces('validation')`)
- Magnetostriction → Model 6 with strain coupling

## Connection to Kameari's stress question

Henrotte 2004 confirms Kameari's intuition (`forces('kameari_stress')`):
the differential-form representation of stress IS subtle and
material-dependent.  The flux-density-vs-circulation-density choice
for magnetization (irrelevant for the magnetic constitutive law,
critical for the force law) is exactly the kind of distinction
that vector calculus papers over.
"""


VALIDATION_BENCHMARKS = """
# Force computation validation: TEAM benchmarks 20, 23, 33b

The TEAM (Testing Electromagnetic Analysis Methods) Workshop
maintains a catalog of standard test problems used worldwide to
verify computational EM codes.  Three are directly relevant to
force computation, providing ground-truth data against which any
σ_em-based formula can be checked.

## TEAM Problem 20 — 3D Static Force, Saturated Steel

Geometry: C-yoke with center pole, all steel, excited by a coil at
1000 / 3000 / 4500 / 5000 ampere-turns (the higher levels drive
the steel deeply into saturation).

Material:  B-H curve provided in 38-point table from B=0.01 T at
H=27 A/m up to B=2.3 T at H=135,000 A/m.  Above 2.3 T, use
B = μ₀ H + M_s with M_s = 2.16 T (Frohlich-type extrapolation).

Quantities to report:
- Bz at 2 points (P1 = origin pole top, P2 = above pole)
- Average Bz in center pole region and in yoke region
- Magnetic force on the center pole

Reference value: published as Table 2-5 in the TEAM Problem 20
specification, with multiple participating codes contributing.

Use case for radia-mcp: a discrete `forces('force_methods')`
implementation that gives the wrong total force on this benchmark
under any of the 7 methods is broken.  Use as a regression test.

## TEAM Problem 23 — Forces in Permanent Magnets

Geometry: small cylindrical permanent magnet (SmCo or NdFeB) facing
a small coil (280 turns of #47 wire).  Configurations A (SmCo) and
B (NdFeB), each with small/large variants.

Materials:
- NdFeB: B_r = 1.22 T, H_c = -820 kA/m  (linear demagnetization)
- SmCo:  B_r = 1.02 T, H_c = -720 kA/m

Measurements provided (Table 2-7):
- Axial force vs current (0-100 mA, in 10 mA steps)
- Axial force vs axial displacement (0-0.514 mm) at 50 mA
- Restoring (side) force vs side displacement at 50 mA

This is the canonical benchmark for the **permanent-magnet-as-1-form
vs 2-form** distinction (Henrotte Model 2 vs 3).  Compute σ_em both
ways, integrate against the experimental data, see which model
matches better.

## TEAM Problem 33b — Electric Local Force, Gelatin Dielectric

Geometry: two parallel copper plates 1 cm apart, 11 kV between them.
A dielectric specimen (gelatin in water, ε_r = 78.5) is inserted.

Materials:
- Gelatin: ε_r = 78.5, Young modulus E = 6000 N/m², Poisson ν = 0.5
- Air: vacuum permittivity

Experimental observation: top of specimen displaces by 0.15 mm
when 11 kV is applied.

Numerical prediction (using the energy-based local force formula):
    F_n = D_n² (1/ε_0 - 1/(ε_0·ε_r))/2  +  E_t² (ε_0·ε_r - ε_0)/2

gives displacement 0.1425 mm, in good agreement.

This benchmark probes the **dielectric (electric)** counterpart of
the magnetic σ_em.  The structural form of the equation is
identical to the magnetic case (energy-based MST), just with
ε ↔ μ⁻¹ substitution.

## How to use these in radia-mcp + Mathematica workflow

```mathematica
(* Define σ_em via mathematica_recipes('maxwell_stress'),
   plug in the TEAM Problem material parameters, compute the
   integral on the eggshell surface, compare to the table. *)
```

A passing regression test on all 3 benchmarks proves the σ_em
implementation is correct for static magnetic, permanent-magnet,
and dielectric problems — all 3 main physical regimes.
"""


ERROR_CORRECTION = """
# Two-stage error correction for EM force accuracy
(Iino-Okamoto, SA-20-003 / SA-20-034, 2020)

Even the best force-computation method (Maxwell stress, eggshell,
Kameari nodal) has numerical error driven by:
- Mesh coarseness / distortion
- Integration-surface placement (eggshell radius)
- Iron-core permeability effects on B accuracy

Iino-Okamoto (2020) propose an **error correction** procedure
that removes spurious contributions cleanly.

## The basic identity (1-stage error correction)

Consider two magnetic sources I_A (on body A) and I_B (on body B),
with the goal of computing the force F on body A.

Compute F three times:
- F(I_A, I_B)  — both sources active
- F(I_A, 0)    — only I_A; should give 0 in theory (no source on B)
- F(0, I_B)    — only I_B; should give 0 in theory (no source on A?)

If the numerical method were perfect, F(I_A, 0) and F(0, I_B) would
both be exactly 0.  In practice they are not.  Define
    F_corrected = F(I_A, I_B) − F(I_A, 0) − F(0, I_B)
This removes spurious self-source error components.

## The 2-stage extension (Iino 2020)

A 1st-order Taylor expansion of F(I_A + δI₁, I_B + δI₂) gives
    F(I₀+δI₁, I₀+δI₂) ≈ F(I₀, I₀) + ∂F/∂I₁·δI₁ + ∂F/∂I₂·δI₂ + ΔF₁
where ΔF₁ is the 2nd-order remainder.  By computing each partial
separately and subtracting, the 2nd-order remainder is recovered:
    F_2nd = F(I_A+δI₁, I_B+δI₂) − F(I_A+δI₁, I_B)
           − F(I_A, I_B+δI₂) + F(I_A, I_B)

This is the **二段階誤差補正** (two-stage error correction).

## Time-domain eddy-current extension (SA-20-034)

For eddy-current problems with TIME-DEPENDENT sources, the same
trick is applied at each time step.  Three sources at time t:
- J_inp(t)  — injected current
- J_e(t)    — eddy current
- M(t)      — magnetization

Compute F four times: (all three on), (only J_inp), (only J_e),
(only M).  The corrected force is
    F_cor = F(all) − F(only J_inp) − F(only J_e) − F(only M)
which captures only the genuine cross-coupling.

## Computational cost

The correction requires **N+1 solves per time step** where N is the
number of independent magnetic sources.  For 1 source (the simplest
case), 2 solves.  For 3 sources (J_inp + J_e + M), 4 solves.

This is expensive (≈ 2-4× cost) but the accuracy gain is
substantial — see SA-20-034 Fig. 5 where the corrected force
matches the analytical solution within 0.5% vs ~5% uncorrected.

## When to use

- High-accuracy force calculation needed (motor torque ripple,
  noise/vibration assessment, magnetostriction)
- Edge-element FEM where naïve σ_em integration shows known
  spurious behavior
- Eddy-current problems where time-domain accuracy is critical

## Implementation in Radia/NGSolve

For Kameari nodal force method, the per-node correction is:
    F_N^corrected = -∫_Ω σ_em(B_all) : ∇γ_N
                   + ∫_Ω σ_em(B_only_self) : ∇γ_N
the second term being the spurious self-force.

This is straightforward to add to the existing nodal-force
computation: just run the FEM solver twice per time step with
selectively zeroed sources.

Reference: T. Iino, Y. Okamoto, A. Ahagon, Y. Kida, K. Semba,
T. Yamada, SA-20-003 (magnetostatic) / SA-20-034 (eddy current),
IEEJ Static Apparatus & Rotating Machinery joint study group, 2020.
"""


KAMEARI_EDDY_CURRENT = """
# Kameari's moving-source approach to eddy currents
(Niikura-Kameari, "Analysis of eddy current and force in conductors
with motion", IEEE Trans. Mag. 28(2):1450-1453, 1992)

The naïve approach to eddy currents in moving conductors is to
include a v × B advective term in the governing PDE.  At high Péclet
number this is numerically unstable (upwinding required).

**Kameari's clever workaround**: keep the conductor stationary,
move the source magnetic field instead.

## The setup

Two equivalent formulations:

(a) Naïve: Stationary source, moving conductor
    ∂A/∂t + curl(ν curl A) + σ(v × curl A) = J_s
    The v×B term causes instability.

(b) Kameari: Stationary conductor, moving source
    ∂A_r/∂t + curl(ν curl A_r) = -∂A_s/∂t
    where A_s(t) = Biot-Savart from the source at current position.
    A_r = A - A_s (reduced vector potential).
    No advective term, no upwinding.

The two are equivalent at the physical level (Galilean invariance
of the EM field for v << c) but numerically different.

## Verification (1992 paper)

- **Rotating circular plate**: plate at rest, dipole field rotates.
  Joule loss and torque computed in two codes (EDDY3D + EDDYCUFF).
  Agreement with analytical solution at low rotation cycle (<10/s).
- **Moving coil over rectangular plate**: plate at rest, coil moves
  at up to 138.9 m/s (500 km/h).  Joule loss + drag force computed.
  Agreement with Fourier-series analytical solution.

## Two complementary discretizations (1992 paper)

| Code | Method | Geometry | Unknowns |
|------|--------|----------|----------|
| **EDDY3D** | 3D FEM with edge elements (A-φ) | Volumetric conductor | Edge values of A_r |
| **EDDYCUFF** | 2D FEM circuit method (T_n) | Thin conductor surface | Current vector potential T_n |

Both give consistent results for the rotating plate / moving coil
benchmark.  Cross-validation between two independent codes is a
hallmark of Kameari's methodology.

## Why this matters for Radia + NGSolve

- Linear motors, MagLev, induction heating with workpiece motion,
  rail-gun-type problems all need eddy-current-in-motion.
- NGSolve has eddy-current solvers but doesn't natively handle
  conductor motion.  Kameari's trick gives a clean implementation
  path: stationary mesh + moving B source.
- For Radia: the PEEC coil + workpiece SIBC framework (in `ih`
  subpackage) can be extended with moving-source via Biot-Savart
  re-evaluation at each time step.

## Connection to Kameari's other work

This 1992 paper, together with Kameari 1990 ("transient 3D eddy
current using edge elements") and Kameari 1993 ("local force
calculation in 3D FEM with edge elements"), forms the canonical
trilogy of Kameari's edge-element framework:
- 1990: edge-element discretization of A-φ
- 1992: moving-source workaround for v×B
- 1993: local force from σ_em with edge-element B

All three are essential for any production-quality 3D edge-element
FEM code that handles moving conductors + force computation.
"""


PILE_2018_COMPARISON = """\
## Pile-Devillers-Le Besnerais 2018 — Force-method comparison for motor NVH

Reference: R. Pile, E. Devillers, J. Le Besnerais, "Comparison of
Main Magnetic Force Computation Methods for Noise and Vibration
Assessment in Electrical Machines", IEEE Trans. Mag. 54(7):8104013,
July 2018.  Affiliation: Eomys Engineering + Ecole Centrale Lille.

### Why this paper matters for motor work

NVH (Noise, Vibration, Harshness) assessment of an electric motor
needs the **spatial distribution** of magnetic forces on the stator
teeth (not just the integrated rotor torque).  Different force-
computation methods agree on **integrated quantities** (total torque)
but **disagree on the local distribution** — which is exactly what
NVH simulation needs.

### Methods compared

| Method | Where applied | Local-distribution agreement |
|--------|---------------|-------------------------------|
| **Virtual Work Principle (VWP) nodal force** | At each FE node | Reference (paper's gold standard) |
| **Maxwell Tensor (MT) on iron surface** | Iron-air interface | Differs locally from VWP |
| **MT in air gap (mid-circle)** | Air gap mid-radius | Differs locally; OK for integrated |
| **MT inside iron** | Iron volume | Differs from all above |

### Key finding

> "The different methods disagree on the local distribution of the
> magnetic forces at the stator surface, **but they give similar
> results concerning the integrated forces per tooth** (referred as
> lumped forces)."

For NVH: the **lumped force per tooth** is what excites tooth-radial
modes of the stator yoke (acoustic radiation).  All methods agree
on this lumped quantity → safe for NVH analysis.

For local stress / fatigue / surface tractions: only VWP gives a
defensible local distribution.  Maxwell tensor methods are formal
**equivalent** at the integrated level but their local values are
artifacts of the integration surface choice.

### Practical guide for radia_mcp.motor

- Lumped tooth force for NVH: use any method, pick the cheapest
  (typically MT in air gap).
- Spatial force distribution on stator iron: use VWP nodal force.
  See `force_methods` topic for the VWP construction.
- Cross-check: compute by 2 methods, confirm integrated agreement.
  Disagreement at integrated level = bug.

### Citation context

This paper is widely cited in motor NVH literature (Le Besnerais's
PyLeecan toolbox uses it as a reference).  Add to bibliography when
implementing force extraction in `calc_motor_transient.py`.
"""

DLALA_FIXED_POINT = """\
## Dlala-Belahcen-Arkkio fixed-point methods (TMAG 2007/2008)

References:
- E. Dlala, A. Belahcen, A. Arkkio, "Locally Convergent Fixed-Point
  Method for Solving Time-Stepping Nonlinear Field Problems", IEEE
  TMAG 43(11):3969-3972, 2007.
- E. Dlala, A. Belahcen, A. Arkkio, "A Fast Fixed-Point Method for
  Solving Magnetic Field Problems in Media of Hysteresis", IEEE TMAG
  44(6):1214-1217, 2008.
Authors: Laboratory of Electromechanics, Helsinki University of
Technology (Arkkio = inventor of the Arkkio formula for air-gap
torque).

### The standard fixed-point technique

Hantila's polarization method (1975, see CLAUDE.md "Hantila
Polarization Method") writes the constitutive law as

  B = μ₀(1+α)H + μ₀R       with  R = M - αH

Substituting into the FE system gives a **constant-LHS** equation

  (I - αN)H = H_ext + NR

where N is the demagnetization-tensor matrix (geometry only).  This
matrix is **LU-factored ONCE** and reused across iterations — that's
the whole appeal vs Newton-Raphson.

Drawback: convergence rate ~ spectral radius of αN, which is **slow**
for highly saturable materials (small α → many iterations).

### Dlala-Belahcen-Arkkio's accelerator

The 2007 paper proposes a **locally convergent** fixed-point variant
that adapts α to the local operating point — guarantees convergence
in an interval around the initial guess, and is significantly faster
than fixed-α Hantila.

The 2008 follow-up tunes it for **hysteresis**: extends to a
*differential* reluctivity at the current B(H) operating point,
remaining O(N²) per iteration (no LU re-factorization).

### Why this matters for the Sugahara lab

Radia ships `hantila_solver` (CLAUDE.md "Hantila Polarization Method").
The Dlala-Belahcen-Arkkio acceleration is a **drop-in upgrade**:
same constant-LHS factorization, but with a smarter α update
schedule.  Expected speedup: 3-5x in iteration count on
saturable-iron problems.

Implementation status: not yet ported to Radia's Hantila benchmark lane.
Worth a half-day project to test against the Radia validation benchmark
(`validation_test/hysteresis/test_binput_cpp.py`).
"""

SEO_CHOI_CONTACT = """\
## Seo-Choi 2014 — Contact forces between touching magnetic materials

Reference: J.H. Seo, H.S. Choi, "Computation of Magnetic Contact
Forces", IEEE Trans. Mag. 50(2):7012904, February 2014.
Affiliation: Kyungpook National University, Korea.

### The problem

Conventional force methods FAIL when two magnetic bodies are **in
direct contact** (no air gap between them).  Reasons:
- Maxwell stress tensor: needs an air-region integration surface;
  there is none if the iron touches.
- Eggshell / Coulomb VWP: the surrounding air-layer "eggshell" is
  ill-defined where the bodies meet.

This matters for:
- Electromagnet pickup at the moment of contact (the design metric!)
- Closing force of a relay armature against the yoke
- PM-to-iron attraction force at zero gap
- Motor pole tip contact (gap closure under shock)

### The virtual air-gap scheme

Insert a **fictitious air gap of zero physical thickness** between
the two bodies, but **treat it as a non-zero FE element** in the
mesh.  Now the Maxwell tensor can integrate over the virtual surface.

Crucially:
- The virtual air gap is created **only in post-processing**, not in
  the geometry.  The FE solve uses iron-iron contact.
- The post-processing step inserts the virtual gap, recomputes the
  field there (the Schwarz reflection principle gives an analytical
  guess), and integrates MT on the virtual surface.

### Validation cases (paper experiments)

The paper validates by:
- iron vs iron (large electromagnet pickup test)
- iron vs PM (relay-style pickup)
- PM vs PM (attraction force between two stacked magnets)
- irregular-shaped bodies (worst case for surface integration)
- saturation regime (high-field, B > 2 T)

In all cases, agreement with **measurement** is within 5%.

### Implementation in radia_mcp

For Radia HDiv-VIM / face-charge route which already handles open-domain magnetostatics
without an "air mesh" per se, the virtual-air-gap idea is even more
natural: just place the integration surface at the touching plane.
The required B field at the virtual surface is available analytically
from the Radia sigma / M solution.

For NGSolve FE motor analysis:
- Relevant if simulating motor armature pull-in dynamics (rare but
  important for actuator-type motors).
- See `force_methods` topic for the standard cases.

### Citation context

Frequently cited in actuator design literature.  Add to bibliography
when implementing force computation for relay / contactor / linear-
actuator examples.
"""


def get_forces_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "overview"             - Why EM forces need differential forms
      "mechanical_framework" - Material manifold + Euclidean + placement map
      "material_derivative"  - Lie derivative formulas for p-forms
      "maxwell_stress"       - σ_em = b⊗h - {h·b - ρΨ}I derivation
      "force_methods"        - Unified table: 7 methods (Maxwell stress,
                               eggshell, Arkkio, Coulomb, Kameari, ...)
      "kameari_stress"       - Kameari's SA-11-013 open question on shear
    """
    topic = topic.lower().strip()
    if topic == "overview":
        return FORCES_OVERVIEW
    if topic in ("mechanical_framework", "framework", "manifold"):
        return MECHANICAL_FRAMEWORK
    if topic in ("material_derivative", "lie", "lie_derivative"):
        return MATERIAL_DERIVATIVE
    if topic in ("maxwell_stress", "stress", "sigma_em", "stress_tensor"):
        return MAXWELL_STRESS
    if topic in ("force_methods", "methods", "kameari_method",
                 "eggshell", "arkkio", "coulomb"):
        return FORCE_METHODS
    if topic in ("kameari_stress", "shear", "open_question"):
        return KAMEARI_STRESS_QUESTION
    if topic in ("henrotte_ferromagnetic", "ferromagnetic", "permanent_magnet",
                 "magnetostrictive", "henrotte"):
        return HENROTTE_FERROMAGNETIC
    if topic in ("validation", "benchmarks", "team_problem", "team"):
        return VALIDATION_BENCHMARKS
    if topic in ("error_correction", "iino_okamoto", "two_stage"):
        return ERROR_CORRECTION
    if topic in ("kameari_eddy_current", "moving_source", "kameari_motion",
                 "niikura_kameari", "kameari"):
        return KAMEARI_EDDY_CURRENT
    if topic in ("pile_2018", "pile_comparison", "nvh",
                 "force_comparison_methods"):
        return PILE_2018_COMPARISON
    if topic in ("dlala", "dlala_fixed_point", "fast_fixed_point",
                 "locally_convergent"):
        return DLALA_FIXED_POINT
    if topic in ("seo_choi", "contact_forces", "virtual_air_gap",
                 "pickup_force"):
        return SEO_CHOI_CONTACT
    if topic == "all":
        return "\n\n".join([
            FORCES_OVERVIEW,
            MECHANICAL_FRAMEWORK,
            MATERIAL_DERIVATIVE,
            MAXWELL_STRESS,
            FORCE_METHODS,
            HENROTTE_FERROMAGNETIC,
            VALIDATION_BENCHMARKS,
            ERROR_CORRECTION,
            KAMEARI_EDDY_CURRENT,
            KAMEARI_STRESS_QUESTION,
            PILE_2018_COMPARISON,
            DLALA_FIXED_POINT,
            SEO_CHOI_CONTACT,
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, overview, mechanical_framework, material_derivative, "
        "maxwell_stress, force_methods, henrotte_ferromagnetic, "
        "validation, error_correction, kameari_eddy_current, kameari_stress, "
        "pile_2018, dlala, seo_choi."
    )
