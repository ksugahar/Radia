"""
Mathematica (Wolfram Language) recipes for symbolic verification of
differential-form identities, Whitney element formulas, and Maxwell
equations.

These recipes pair this server (`radia_mcp.differential_forms`, which
holds the theory) with the sibling MCP server `radia_mcp.mathematica`
(which exposes wolframscript as `mathematica_evaluate`, etc.).  The
intended workflow is:

  1. Read theory in `differential_forms_basics` / `whitney` / `de_rham`.
  2. Pick a verifiable identity (e.g. d² = 0, Stokes' theorem on a
     tetrahedron face, curl of an edge element).
  3. Use the Wolfram recipes below as input to `mathematica_evaluate`
     via the `mathematica` MCP server.
  4. Get a symbolic confirmation (or counter-example).

Distilled from:
  - 微分形式.pdf + Bossavit Ch.5 (theory)
  - Mathematica V14.2 (Curl, Div, Grad, FullSimplify, Series, TeXForm)
  - Internal experience in symbolic-FEM verification (Kameari-style)
"""

VERIFY_DSQUARED = """
# Recipe: Verify d² = 0 (exterior derivative is nilpotent)

The most fundamental identity in differential forms.  In 3D vector
calculus this is two separate facts:
  - rot(grad f) = 0  for any scalar field f
  - div(rot a) = 0  for any vector field a

Both follow from d² = 0.

## Wolfram Language input
Copy-paste into `mathematica_evaluate` (radia_mcp.mathematica):

```mathematica
(* rot(grad f) = 0 — use direct Grad[f[x,y,z], ...], DO NOT do
   f = f[x,y,z] which is self-referential and recurses *)
FullSimplify[Curl[Grad[f[x, y, z], {x, y, z}], {x, y, z}]]
(* Expected: {0, 0, 0} *)

(* div(rot a) = 0 *)
aVec = {a1[x, y, z], a2[x, y, z], a3[x, y, z]};
FullSimplify[Div[Curl[aVec, {x, y, z}], {x, y, z}]]
(* Expected: 0 *)
```

**Gotcha**: Mathematica's pattern matcher treats  `f = f[x, y, z]`
as a recursive redefinition (infinite loop, RecursionLimit error).
Always use the form `Grad[f[x, y, z], {x, y, z}]` directly, with
`f` as the symbolic head — never assign `f = f[x, y, z]` first.

## What this proves
At the symbolic level, Mathematica verifies that smooth fields satisfy
d² = 0 identically.  This is the algebraic identity behind:
  - RG = 0 (discrete rot-grad)
  - DR = 0 (discrete div-rot)
in the Whitney complex.

## Caveats
For NON-smooth fields (with jumps across material interfaces), d²
still holds in the **distributional / weak** sense, but symbolic
verification needs a different setup.  Use Mathematica's
`DiracDelta` and `HeavisideTheta` for that.

See `differential_forms_basics('exterior_derivative')` for theory.
"""

VERIFY_STOKES = """
# Recipe: Verify Stokes' theorem on a tetrahedron face

Stokes' theorem  ∫_S (rot a) · n dS = ∮_∂S a · t dl  is the foundation
of edge-element discretization.  Verify it symbolically on a single
triangular face of a tetrahedron.

## Setup
Triangular face with vertices p_l, p_m, p_n.  Parametrize it with
barycentric coordinates (λ_l, λ_m, λ_n), λ_l + λ_m + λ_n = 1.

## Wolfram Language input

```mathematica
(* Tetrahedron with simple vertices for verification *)
pl = {0, 0, 0}; pm = {1, 0, 0}; pn = {0, 1, 0};

(* Generic linear vector field — any vector field with linear components *)
a[x_, y_, z_] := {a1 + b1 x + c1 y + d1 z,
                  a2 + b2 x + c2 y + d2 z,
                  a3 + b3 x + c3 y + d3 z};

(* Curl of a (constant in 3D since coefficients are linear) *)
rotA = Curl[a[x, y, z], {x, y, z}];

(* Normal to the triangle (z direction here) and area = 1/2 *)
nVec = {0, 0, 1};
fluxSurface = (rotA . nVec) * (1/2);   (* uniform field × area *)

(* Boundary integral: 3 edges  pl→pm, pm→pn, pn→pl *)
edge[p1_, p2_] := Module[{tvec = p2 - p1},
  Integrate[a[p1[[1]] + s tvec[[1]],
              p1[[2]] + s tvec[[2]],
              p1[[3]] + s tvec[[3]]] . tvec,
            {s, 0, 1}]
];
boundary = edge[pl, pm] + edge[pm, pn] + edge[pn, pl];

(* Check equality *)
FullSimplify[fluxSurface - boundary]
(* Expected: 0 *)
```

## Why this matters
This is how the discrete relation  (rot h)_f = ∑_e R_{f,e} h_e  is
PROVED to be exact (not approximate) when h is in W^1 (edge element
space).  The Stokes identity on each face is the building block.

See `differential_forms_basics('stokes')` and
`differential_forms_whitney('edge')`.
"""

WHITNEY_EDGE_ELEMENT_FORMULA = """
# Recipe: Build Whitney edge element  w_e = w_m ∇w_n − w_n ∇w_m

Verify the formula symbolically, check tangential continuity, compute
its curl (which gives 2 ∇w_m × ∇w_n — a constant per tetrahedron).

## Setup
Tetrahedron with vertices p_1, p_2, p_3, p_4 (use generic vertices
for full generality).

## Wolfram Language input

```mathematica
(* Barycentric coordinates as PLAIN EXPRESSIONS, NOT function defs.
   This is crucial: defining lam1[x_,y_,z_]:=... and then computing
   Grad[lam1[s,0,0],...] gives 0 because Grad sees a constant after
   substitution.  Keep them as symbolic expressions in x, y, z and
   substitute coordinates AFTER taking Grad. *)
lam1 = 1 - x - y - z;
lam2 = x;
lam3 = y;
lam4 = z;

(* Whitney edge element from node 1 to node 2 *)
we12 = lam1 * Grad[lam2, {x, y, z}] - lam2 * Grad[lam1, {x, y, z}];
Print["we12 = ", we12];
(* Expected: {1 - y - z, x, x}   on the unit tetrahedron *)

(* Verify circulation along edge 1-2 (= 1).  Edge goes (0,0,0)->(1,0,0),
   parametrize by s ∈ [0,1] and integrate after substitution. *)
edge12Vec = {1, 0, 0};
circulation12 = Integrate[
  (we12 /. {x -> s, y -> 0, z -> 0}) . edge12Vec, {s, 0, 1}
];
Print["circulation along edge 1-2: ", circulation12];
(* Expected: 1 *)

(* Verify circulation along edge 1-3 is zero *)
circulation13 = Integrate[
  (we12 /. {x -> 0, y -> s, z -> 0}) . {0, 1, 0}, {s, 0, 1}
];
Print["circulation along edge 1-3: ", circulation13];
(* Expected: 0 *)

(* Compute curl — should be 2 ∇lam1 × ∇lam2 (constant on tet) *)
rotWe12 = FullSimplify[Curl[we12, {x, y, z}]];
Print["curl: ", rotWe12];
(* Expected: {0, -2, 2}  (= 2 × (∇lam1 × ∇lam2)) *)
```

## What this gives
The recipe demonstrates Bossavit's defining duality:
  ∫_{e'} τ · w_e  =  δ_{e,e'}

and shows that the curl is **constant per tetrahedron**, which is
why face elements w_f = 2(...) capture the curl space exactly.

See `differential_forms_whitney('edge')`.
"""

VERIFY_MAXWELL_IDENTITY = """
# Recipe: Verify Maxwell identity  d(F) = 0  (Faraday + ∇·B = 0)

In differential-form language, two of Maxwell's equations are a
single identity:  dF = 0, where F is the 2-form

   F = (B_x dy∧dz + B_y dz∧dx + B_z dx∧dy)
       − (1/c) (E_x dx∧dt + E_y dy∧dt + E_z dz∧dt)

Then dF = 0 reads:
   ∂_t B + rot E = 0           (Faraday)
   div B = 0                   (no magnetic monopoles)

## Wolfram Language input

```mathematica
(* Vector potential a and B = rot a, E = -grad phi - ∂_t a.
   Use names aVec, phiF, EE — NOT a, phi, E — to avoid:
     - self-referential assignments (RecursionLimit)
     - shadowing Wolfram's built-in symbol  E  (Euler's number) *)
aVec = {a1[x, y, z, t], a2[x, y, z, t], a3[x, y, z, t]};
phiF = phi[x, y, z, t];

B = Curl[aVec, {x, y, z}];
EE = -Grad[phiF, {x, y, z}] - D[aVec, t];

(* Faraday: ∂_t B + rot E ?= 0 *)
faradayLHS = D[B, t] + Curl[EE, {x, y, z}];
Print["Faraday: ", FullSimplify[faradayLHS]];
(* Expected: {0, 0, 0} *)

(* Gauss for magnetism: div B ?= 0 *)
gaussB = Div[B, {x, y, z}];
Print["div B: ", FullSimplify[gaussB]];
(* Expected: 0 *)
```

## What this confirms
The two homogeneous Maxwell equations are AUTOMATIC consequences of
B = rot a, E = -grad ϕ - ∂_t a.  They do not need to be imposed.

This is why FEM using edge elements for A is **exact** with respect
to ∇·B = 0: structure-preserving discretization.

See `differential_forms_maxwell('maxwell')` and
`differential_forms_de_rham('maxwell_house')`.
"""

HODGE_STAR_VERIFICATION = """
# Recipe: Hodge star ⋆ — symbolic verification

Verify the Hodge star formulas:
  ⋆⋆ω = (-1)^{k(n-k)} ω  (in Euclidean R^3:  ⋆⋆ = id)
  ⋆1 = dx∧dy∧dz
  ⋆dx = dy∧dz
  ⋆(dx∧dy) = dz

and the constitutive relations d = εe, b = μh as Hodge-star scaled
maps.

## Wolfram Language: tensor approach (rank-k antisymmetric arrays)

```mathematica
(* For R^3 Euclidean: ε^{ijk} is the Levi-Civita symbol *)
levi = LeviCivitaTensor[3];

(* 1-form ω = (ω1, ω2, ω3): act of Hodge ⋆ gives a 2-form *)
omega1 = {w1, w2, w3};
hodge1to2[v_] := Table[Sum[levi[[i, j, k]] v[[k]], {k, 3}],
                       {i, 1, 3}, {j, 1, 3}] / 2;
(* This gives the antisymmetric matrix corresponding to a 2-form *)
H = hodge1to2[omega1];
Print["⋆(1-form): ", H];

(* Double Hodge — should return original *)
omega1back = {H[[2, 3]], H[[3, 1]], H[[1, 2]]};
Print["⋆⋆(1-form) - original: ", FullSimplify[omega1back - omega1]];
(* Expected: {0, 0, 0} *)
```

## Wolfram Language: differential-form approach (component notation)

```mathematica
(* B = 2-form components on dy∧dz, dz∧dx, dx∧dy *)
B = {Bx, By, Bz};

(* Constitutive law:  H = (1/μ) ⋆ B  *)
mu0 = 4 Pi * 10^-7;
mu = mu0 * mur;       (* relative permeability *)
H = B / mu;           (* trivial in Cartesian Euclidean *)

(* In Minkowski (1+3) signature, sign matters: *)
metric = {{-1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
detg = Det[metric];   (* should be -1 *)
```

## What this teaches
The Hodge star is **purely algebraic** in component form.  All the
"material" content of Maxwell's equations is encoded in the scalar
prefactor (ε, μ, σ).  Geometry-vs-material separation, which is the
philosophical thrust of FEEC.

See `differential_forms_maxwell('hodge')`.
"""

TEX_OUTPUT_FOR_PAPER = """
# Recipe: TeXForm output for paper writing

Get LaTeX strings for differential-form expressions, ready to paste
into a paper.

## Wolfram Language input

```mathematica
(* Whitney edge element formula *)
expr1 = lambdaI * Grad[lambdaJ] - lambdaJ * Grad[lambdaI];
Print["TeX: ", ToString[TeXForm[expr1]]];

(* Maxwell stress tensor *)
T = mu0 * Outer[Times, H, H] - (1/2) * mu0 * (H . H) * IdentityMatrix[3];
Print["TeX (T_ij): ", ToString[TeXForm[T]]];

(* Skin depth *)
delta = Sqrt[2 / (omega * mu * sigma)];
Print["TeX (delta): ", ToString[TeXForm[delta]]];

(* Hodge Laplacian: d δ + δ d *)
(* (symbolic — represent d and δ as operators via Inactive[]) *)
HodgeLap = Inactive[d][Inactive[delta][omega]] +
           Inactive[delta][Inactive[d][omega]];
Print["TeX (Hodge Lap): ", ToString[TeXForm[HodgeLap]]];
```

## Tip via the high-level helper
The radia_mcp.mathematica server exposes `mathematica_to_tex` which
wraps this for you:

```python
result = mathematica_to_tex("lambdaI * Grad[lambdaJ] - lambdaJ * Grad[lambdaI]")
# → result['tex'] = "\\lambda_I \\,\\nabla \\lambda_J - \\lambda_J \\,\\nabla \\lambda_I"
```

(Use `mathematica_to_tex` from `radia_mcp.mathematica`.)
"""

LORENTZ_TRANSFORMATION = """
# Recipe: Symbolic Lorentz transformation and Maxwell covariance

Verify that Maxwell's equations are Lorentz-invariant in 4D
differential-form language.

## Wolfram Language input

```mathematica
(* Lorentz boost along x with speed v, gamma = 1/sqrt(1 - beta^2) *)
beta = v/c;
gam = 1/Sqrt[1 - beta^2];
Lambda = {
  {gam,      -gam beta, 0, 0},
  {-gam beta, gam,      0, 0},
  {0,         0,        1, 0},
  {0,         0,        0, 1}
};
Print["Det[Lambda] = ", FullSimplify[Det[Lambda]]];
(* Expected: 1 — proper Lorentz transformation *)

(* Electromagnetic 2-form F_{μν} as antisymmetric tensor *)
F = {{0,       -Ex/c,  -Ey/c,  -Ez/c},
     {Ex/c,     0,      Bz,    -By  },
     {Ey/c,   -Bz,      0,      Bx  },
     {Ez/c,    By,    -Bx,      0  }};

(* Transformed F' = Lambda . F . Lambda^T *)
Fprime = Lambda . F . Transpose[Lambda];
ExPrime = -c * Fprime[[1, 2]];
ByPrime = -Fprime[[2, 4]];
BzPrime = Fprime[[2, 3]];

Print["E_x' = ", FullSimplify[ExPrime]];
Print["B_y' = ", FullSimplify[ByPrime]];
Print["B_z' = ", FullSimplify[BzPrime]];
(* Expected:
   E_x' = E_x
   B_y' = gam (B_y + v E_z / c^2)
   B_z' = gam (B_z - v E_y / c^2)   *)
```

## What this confirms
The "tensor of electromagnetism" F transforms cleanly under Lorentz,
giving the familiar E ↔ B mixing.  This is why the 4D form
dF = 0, dG = J is a manifestly covariant formulation.

See `differential_forms_maxwell('maxwell')` 4-D section.
"""

MAXWELL_STRESS_TENSOR = """
# Recipe: Maxwell stress tensor σ_em — symbolic derivation

The Bossavit-Henrotte 2004 Handbook gives:
    σ_em = b ⊗ h̃  -  {h̃·b - ρΨ(b)} I

For a linear isotropic magnetic material (energy density ρΨ = b·h̃/2 =
|b|²/(2μ)), this reduces to the familiar Maxwell stress.  Verify it
symbolically with Mathematica.

## Wolfram Language input

```mathematica
(* B field as a 3-vector (proxy for the 2-form b) *)
B = {Bx, By, Bz};

(* Linear isotropic material:  H = B/μ *)
H = B / mu;

(* Energy density:  ρΨ = (1/2) μ |H|^2 = (1/2) |B|^2 / μ *)
rhoPsi = (B . B) / (2 mu);

(* Outer/dyadic product b ⊗ h̃ *)
bOuterH = Outer[Times, B, H];

(* Identity tensor *)
Id3 = IdentityMatrix[3];

(* Maxwell stress tensor *)
sigma = bOuterH - (H . B - rhoPsi) * Id3;
Print["σ_em = ", MatrixForm[Simplify[sigma]]];

(* Check symmetry — for linear isotropic, σ_em should be symmetric *)
Print["σ - σ^T = ",
      MatrixForm[Simplify[sigma - Transpose[sigma]]]];
(* Expected: zero matrix *)

(* Trace — the "isotropic pressure" component *)
Print["tr(σ_em) = ", Simplify[Tr[sigma]]];
(* Expected: -|B|^2 / (2 μ)  (negative of magnetic energy density —
   the isotropic component is a *pressure*, not a tension) *)

(* Special case B = (0, 0, B0): σ should be diag(-B0²/(2μ),
   -B0²/(2μ), +B0²/(2μ)). The (+) on the diagonal along B is the
   FARADAY TENSION (field-line tension). The (-) perpendicular is
   the FARADAY PRESSURE (field-line pressure). *)
sigmaZ = sigma /. {Bx -> 0, By -> 0, Bz -> B0};
Print["σ_em with B along z: ", MatrixForm[Simplify[sigmaZ]]];
(* The diagonal (-, -, +) is the canonical Maxwell-stress signature *)
```

## Force on a rigid body from σ_em

```mathematica
(* Test: force on a permanent magnet block in uniform external H_ext.
   Energy density picks up an "applied" term  -B · H_ext.
   The force comes from the gradient of energy, i.e. div σ_em. *)

(* Set up a uniform external B = {0, 0, B0} field *)
B = {0, 0, B0};
H = B / mu0;
rhoPsi = (B . B) / (2 mu0);
sigma = Outer[Times, B, H] - (H . B - rhoPsi) IdentityMatrix[3];
Print["σ in vacuum (uniform field): ", MatrixForm[Simplify[sigma]]];

(* Force density ρf_em = div σ_em.  For a UNIFORM field, this is 0
   in the volume (Newton's third law internal forces cancel), and
   the force on a body arises only from σ_em jumps at interfaces. *)
```

## Symbolic check of force on a wire in a uniform field

```mathematica
(* Lorentz force: F = ∫ J × B dV  — for a wire with current I, length L *)
(* Compare with the Maxwell-stress integral over a surface around the wire *)

(* Wire: straight along ẑ, current I along ẑ *)
(* Field around wire (cylindrical, B_φ = μ_0 I / (2π r))  *)
(* External uniform field B_ext = B0 ẑ *)

(* Total field is superposition. The Maxwell stress integral on a
   cylinder of radius R around the wire gives F = I L × B_ext.       *)
(* — direct computation via Mathematica needs cylindrical coordinates  *)

(* This is the textbook "F = IL × B" derivation, but here it emerges
   from the Maxwell stress tensor, not from a postulate.            *)
```

## Anisotropic / magnetostrictive extension

Mathematica scales naturally to anisotropic materials.  Replace
`H = B/mu` with `H = MuInv.B` where `MuInv` is a 3×3 inverse-tensor.
The stress tensor automatically picks up the right anisotropy
correction.  For magnetostriction, set
  rhoPsi = ψ_mag(B) + ψ_elastic(strain) + ψ_coupling(B, strain)
and the stress tensor σ_em automatically includes the
magnetostrictive coupling term ∂_strain ψ_coupling.

## Why use this

Hand-derivation of σ_em for anisotropic / magnetostrictive materials
is error-prone.  Mathematica produces the correct stress tensor in
closed form, ready for codegen into NGSolve/Radia C++ integrators.
This is Bossavit-Henrotte's framework made executable.

See `differential_forms_forces('maxwell_stress')` for theory and
`differential_forms_forces('force_methods')` for the 7-method table
including Kameari's local force method.
"""

KELVIN_TRANSFORM = """
# Recipe: Kelvin transformation — `kelvin_factor` from conformal weight

The Kelvin transform maps an unbounded exterior domain to a bounded
inverted-coordinate domain via inversion in a sphere of radius R:

    T:  x  →  y = R² x / |x|²

Conformal factor:  λ(y) = R² / |y|²   (T is conformal).

The "kelvin_factor" that everyone struggles to derive by hand is just
the conformal weight of a k-form under λ in dimension n:

    factor_k  =  λ^((n − 2k) / 2)        (n = 3)

| k | physical entity            | factor (R = 1) |
|---|----------------------------|----------------|
| 0 | scalar potential ϕ         | 1/|y|          |
| 1 | vector potential A         | |y|            |
| 2 | magnetic flux B            | |y|³           |

This factor is forced by the requirement that the bilinear-form
integrand (e.g. ∫|∇ϕ|² dV) be invariant under the Kelvin map.  Vector
calculus would require chain-rule on the 3×3 Jacobian and the
two-step Laplacian (27 terms) — but the conformal-weight formula
gives all three factors instantly.

## Wolfram Language input — verify the k=0 case (scalar potential)

```mathematica
(* Take a harmonic test function: dipole  phi = x1 / |x|^3 *)
phiFn[xx_, yy_, zz_] := xx / (xx^2 + yy^2 + zz^2)^(3/2);

(* Verify phi is harmonic *)
laplPhi = FullSimplify[
  D[phiFn[xx, yy, zz], {xx, 2}]
  + D[phiFn[xx, yy, zz], {yy, 2}]
  + D[phiFn[xx, yy, zz], {zz, 2}]
];
Print["Laplacian(phi) = ", laplPhi];   (* Expected: 0 *)

(* Kelvin transform: psi(y) = (1/|y|) phi(y/|y|^2) *)
rsq = yx^2 + yy^2 + yz^2;
psi = (1/Sqrt[rsq]) * phiFn[yx/rsq, yy/rsq, yz/rsq];
psiSimp = FullSimplify[psi];
Print["psi(y) = ", psiSimp];           (* Expected: yx *)

(* Verify psi is harmonic — this is the Kelvin theorem *)
laplPsi = FullSimplify[
  D[psiSimp, {yx, 2}] + D[psiSimp, {yy, 2}] + D[psiSimp, {yz, 2}]
];
Print["Laplacian(psi) = ", laplPsi];   (* Expected: 0 *)
```

## Verify the factor is unique

Try different exponents k in the prefactor (1/|y|)^k and check which
preserve harmonicity:

```mathematica
phiFn[xx_, yy_, zz_] := xx / (xx^2 + yy^2 + zz^2)^(3/2);
rsq = yx^2 + yy^2 + yz^2;

Do[
  psi = (1/rsq^(k/2)) * phiFn[yx/rsq, yy/rsq, yz/rsq];
  lapl = FullSimplify[
    D[psi, {yx, 2}] + D[psi, {yy, 2}] + D[psi, {yz, 2}]
  ];
  Print["k = ", k, ":  Laplacian = ", lapl],
  {k, 1, 3}
];
```

Result:
  k = 1:  Laplacian = 0       ← (n-2)/2 = 1/2 for n=3 ⇒ factor = 1/|y|^1
  k = 2:  Laplacian ≠ 0
  k = 3:  Laplacian ≠ 0

So k=1 (i.e. factor = R/|y| with R=1) is FORCED by the conformal
weight formula — no other exponent works.

## Wolfram Language — verify the k=1 case (vector potential A)

For curl-curl problems, the kelvin_factor for A is |y|/R, NOT R/|y|.
This is harder to derive by hand but follows the same formula.

```mathematica
(* Simple harmonic vector field — take A = (0, 0, x) *)
aPhys[xx_, yy_, zz_] := {0, 0, xx};
bPhys = Curl[aPhys[xx, yy, zz], {xx, yy, zz}];
Print["B in physical space = ", bPhys];

(* Kelvin transform of a 1-form picks up factor lambda^(-1/2) = |y|/R *)
(* (formal pullback computation; for paper-quality derivation see
   Bossavit Ch.7 + Arnold-Falk-Winther 2006 §9) *)
```

## Why this matters for Radia + NGSolve

- `radia.axifem` Kelvin: axisymmetric Henrotte + Kelvin shell
  combines TWO conformal effects (axi r² + Kelvin λ).  Mathematica
  can verify the combined factor symbolically.

- `radia_mcp.radia_ngsolve` Kelvin identification: dirichlet vs
  periodic vs Kelvin BCs each pick a different conformal weight.
  Mathematica gives the right one without hand-derivation.

- C++ codegen: the kelvin_factor that NGSolve evaluates at runtime
  in `kelvin_source.py` and the C++ bilinear form integrators can be
  compared bit-for-bit against the symbolic ground truth derived here.

## Pedagogical note

In vector-calculus textbooks, the Kelvin transform is presented as
"a clever trick".  In differential-form language, it is **the
inversion T applied as a pullback**, and the factor is just the
conformal weight λ^((n-2k)/2) of a k-form.  Single formula → all
cases.

See `differential_forms_de_rham('weak')` for the Sobolev-space
framework underlying the Kelvin variational formulation, and
`differential_forms_homology('chain_complex')` for the chain-level
view of the inversion's effect on the de Rham complex.

## Committed, self-tested realization (lab .wls)

The weak-form / Hodge-as-material-modulation half of this recipe is a
committed, self-testing `.wls` suite (run via `mathematica_evaluate` or
`wolframscript -file`), in
`packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry/`:

- `weakform_hodge.wls` — derives the **material modulation**
  `nu' = nu |det P| (P^T P)^{-1}` (P = pullback on the orthonormal 1-form
  basis), shows it is the SAME formula as the scalar Dirichlet weight
  `W = |det J|(J^T J)^{-1}`, and lands the **Nagamine-Yamaguchi-Sugahara
  pullback Kelvin** law: spherical (conformal) `nu' = (r'/R)^2 nu`
  isotropic; cylindrical (non-conformal) `nu' = diag(1,1,(rho'/R)^4) nu`
  anisotropic; plus the nonlinear `Star_nu` tangent `nu I + nu'(B(x)B)/|B|`
  is SPD (elliptic, no fold).  20 assertions, ALL PASS.
- `hodograph.wls` — the 3-axis backbone: cochain-map split (dB=0 metric-free),
  2-D Kelvin inversion conformal => weight-free, Clebsch + helicity
  obstruction (Moffatt), `A_z`-as-Clebsch-potential.  21 assertions, ALL PASS.
- `canonical.wls` — the Hamiltonian / Legendre reading: flux lines are
  canonical equations with `A_z` as Hamiltonian, `H = dw/dB` is a Legendre
  transform, and the hodograph is a canonical transformation.  12 assertions,
  ALL PASS.
- `surface_derham.wls` — the HOIBC / surface-de-Rham reading: global HOIBC
  modes are surface harmonic 1-forms counted by `b1`, while the exact exterior
  DtN has eigenvalues `Lambda_n=(n+1)/R` and HOIBC is a local Pade/Laplace-
  Beltrami approximation; the vector case splits into two distinct de Rham
  Steklov ladders.  12 assertions, ALL PASS.
- `dtn_geometry.wls` — the DtN / Steklov operator as the exterior Hodge star
  condensed onto the boundary: Schur complement, radial Riccati fixed point,
  shifted square-root of the surface Laplacian, self-adjoint positivity, and
  the nonlocal reason HOIBC has an approximation knee.  9 assertions, ALL PASS.

Note these use the **k-form conformal weight as the material tensor**
`nu' = nu |det P|(P^T P)^{-1}` (the energy/cofactor route), which is the
same fact as the `factor_k = lambda^((n-2k)/2)` table above written as a
metric weight; the spherical k=2 (flux B) case there is the `(r'/R)^2`
isotropic `nu'` here.  The topology half is locked by the companion
`basis_functions/cohomology.wls`, which verifies that harmonic 1-form
dimension is the Betti number and is invariant under SPD Hodge metrics.
See `docs/differential_geometry/PREMETRIC_PICTURE.md` and
`docs/clebsch_hodograph/HODOGRAPH_BACKBONE.md`.
"""

DIFFERENTIAL_GEOMETRY_WLS_RUNBOOK = """
# Recipe: Differential-geometry `.wls` runbook for hodograph work

This is the MCP-facing index for the committed Wolfram Language scripts in
`packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry/`.
Use it when the question is "what did the hodograph commits establish, and
which symbolic file should I run before changing the FEM formulation?"

## What is covered

| file | question it answers | assertions |
|---|---|---:|
| `weakform_hodge.wls` | Is the weak form a Hodge pairing, and does a coordinate map move only the material/Hodge weight? | 20 |
| `curve_surface_basics.wls` | What curve/surface differential-geometry identities support hodograph charts and stream-function surfaces? | 23 |
| `hodograph.wls` | How do Kelvin, Clebsch, Chaplygin, and the `A_z` potential fit into one hodograph diagram? | 21 |
| `canonical.wls` | What is the Hamiltonian / Legendre / canonical-transform reading of the weak form and hodograph? | 12 |
| `surface_derham.wls` | How do HOIBC global modes split into surface topology plus the analytic exterior DtN? | 12 |
| `dtn_geometry.wls` | What is the DtN/Steklov operator: exterior Hodge condensation, Schur/Riccati flow, shifted `sqrt(-Delta_Gamma)`, and boundary `H^{1/2}` metric? | 9 |

Every file is self-contained and prints `ALL PASS` when the symbolic
assertions hold.  The companion topology guard is
`../basis_functions/cohomology.wls` (16 assertions, ALL PASS): it proves that
the harmonic-loop count lives in `d`, not in the Hodge/material tensor.

## Run locally

```powershell
cd packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry
wolframscript -file weakform_hodge.wls
wolframscript -file curve_surface_basics.wls
wolframscript -file hodograph.wls
wolframscript -file canonical.wls
wolframscript -file surface_derham.wls
wolframscript -file dtn_geometry.wls
wolframscript -file ../basis_functions/cohomology.wls
```

## Run through MCP

Use the sibling `mcp-server-mathematica` bridge:

```python
mathematica_evaluate(code=<contents of the .wls file>, timeout=120)
```

Use this `mcp-server-differential-forms` tool for the map:

- `differential_forms_mathematica_recipes("weakform_hodge")`
- `differential_forms_mathematica_recipes("curve_surface")`
- `differential_forms_mathematica_recipes("hodograph")`
- `differential_forms_mathematica_recipes("canonical")`
- `differential_forms_mathematica_recipes("surface_derham")`
- `differential_forms_mathematica_recipes("dtn_geometry")`

Cross-reference:

- `docs/clebsch_hodograph/DIFFERENTIAL_GEOMETRY_WLS.md` is the
  human-readable docs twin.
- `packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry/curve_surface_reading_guide.md`
  maps the local Mathematica curve/surface differential-geometry textbook to
  the hodograph and stream-function anchors without copying book text.
"""

CURVE_SURFACE_DG_READING_GUIDE = """
# Recipe: curve_surface_basics.wls -- curve/surface DG for hodograph and SF

Executable anchor:
`packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry/curve_surface_basics.wls`

Local background source: `Mathematica 曲線と曲面の微分幾何.pdf`.
The PDF is NOT required for `pip install radia-mcp` users.  The shipped
knowledge is the `.wls` file and this summary; the book is only a local study
source for improving the assertions.

Packaged guide:
`packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry/curve_surface_reading_guide.md`

## What the `.wls` proves

`curve_surface_basics.wls` is self-contained and prints `ALL PASS`.  It checks
23 assertions:

- curve speed, curvature, and torsion for a circle and helix;
- surface first fundamental form, area density, unit normal, second fundamental
  form, Gaussian curvature, and mean curvature for sphere/cylinder;
- the chart identity `g = J^T J`;
- surface gradient on a sphere;
- stream-function current `K = n x grad_Gamma psi` is tangential and
  `div_Gamma K = 0` on sphere and plane;
- `Delta_Gamma psi = div_Gamma grad_Gamma psi` and
  `Delta_Gamma cos(theta) = -2 cos(theta)/R^2` on a sphere.

## Run

```powershell
cd packages/radia-mcp/src/radia_mcp/mathematica/differential_geometry
wolframscript -file curve_surface_basics.wls
```

Use it before changing:

- stream-function contour extraction on curved surfaces;
- `K = n x grad_Gamma psi` sign/orientation conventions;
- curved-surface current-density norms;
- surface Laplace-Beltrami / DtN / HOIBC logic;
- hodograph chart geometry.

## How it connects to hodograph

Read curve/surface coordinate transformations through the split:

- `d` / pullback / closedness are topology and commute with smooth coordinate
  changes;
- metric / area / inner product / Hodge star are geometry and become explicit
  weights in the weak form.

Then run `hodograph.wls`, `canonical.wls`, and `weakform_hodge.wls`.

## How it connects to stream-function

The surface identity is `K = n x grad_Gamma psi`.  The first fundamental form
controls `grad_Gamma`; the oriented normal controls the sign; the area density
controls norms and integration.  `curve_surface_basics.wls` verifies the core
contracts symbolically before the numerical SF code uses them.

Related anchors:
`surface_derham.wls`, `dtn_geometry.wls`,
`radia_mcp.streamfunction.streamfunction("single_stroke")`, and
`radia_mcp.streamfunction.streamfunction("regularized")`.

## Local search terms

When improving this `.wls` from the local OCR text, search for:

```text
曲線, 曲率, 接ベクトル, 法線, 曲面, 接平面, 第一基本形式,
第二基本形式, 測地線, ガウス曲率, 平均曲率, パラメータ表示,
面積要素, 座標変換, Mathematica
```

Every useful result should be turned into a Radia-owned assertion before it
becomes MCP guidance.
"""

WEAKFORM_HODGE_RECIPE = """
# Recipe: Weak form = Hodge pairing; material = pulled-back Hodge weight

`weakform_hodge.wls` makes the central FEM statement operational:

```text
weak form energy = integral <d test, star_nu d trial>
```

Under a coordinate map, the exterior derivative / cochain side is natural:
`k^* d = d k^*`.  The part that changes is the Hodge star, hence the
material tensor.

The script verifies:

- The scalar Dirichlet energy weight is `W = |det J| (J^T J)^-1`.
- The magnetic material modulation is the same cofactor/Hodge rule:
  `nu' = nu |det P| (P^T P)^-1`.
- The pullback Kelvin laws follow from that one rule:
  spherical inversion gives isotropic `nu' = (r'/R)^2 nu`;
  cylindrical inversion gives anisotropic `diag(1,1,(rho'/R)^4) nu`.
- The nonlinear tangent of `star_nu` is SPD for ordinary saturation, so
  the magnetostatic operator stays elliptic and the material law itself
  does not create a hodograph fold.
- Metric is not curvature: Kelvin and polar metrics are flat; a round sphere
  control case is curved, so the symbolic curvature check is real.
- Transformation optics is the same statement in premetric EM language:
  push the geometry onto the Hodge/material tensor; topology is not in the
  Hodge.

Use this before changing Kelvin weights, bidirectional-coordinate transforms,
or nonlinear-material linearizations.
"""

HODOGRAPH_WLS_RECIPE = """
# Recipe: Hodograph symbolic backbone

`hodograph.wls` is the executable version of the three-axis hodograph map in
`docs/clebsch_hodograph/HODOGRAPH_BACKBONE.md`.

It verifies:

- The topological half (`dB = 0`) is a cochain-map statement and is preserved
  by any smooth coordinate transform.
- A 2-D Kelvin inversion is conformal, so the Hodge weight on 1-forms is
  unchanged; this is why the 2-D hodograph/Kelvin cell is weight-free.
- The complex-potential `(Phi, A_z)` chart is the Cauchy-Riemann net: the
  potential-plane transform is classical, bounded-domain prior art.
- The 3-D Clebsch representation `B = grad(alpha) x grad(beta)` is blocked
  globally by helicity when `integral A.B dV` is nonzero.
- In 2-D, `A_z` is the Clebsch potential and also the flux-line Hamiltonian.
- The Chaplygin field-plane variables `(theta, q=|H|)` linearize the nonlinear
  saturable 2-D equation into a self-adjoint linear elliptic equation:

```text
d/dq((q/mu) A_q) + ((mu q)'/(mu^2 q)) A_thetatheta = 0
```

The ellipticity condition is `(mu q)' > 0`, equivalently
`nu + nu' |B| > 0`.  Ordinary magnetic saturation therefore has no gas-dynamics
style limiting line.
"""

CANONICAL_HODOGRAPH_RECIPE = """
# Recipe: Canonical / Hamiltonian reading of the hodograph

`canonical.wls` explains the hodograph from the variational side:
the weak form is `delta(action)=0`, and the Hamiltonian picture is its
Legendre shadow.

It verifies:

- 2-D flux lines satisfy Hamilton's canonical equations with `A_z` as the
  Hamiltonian:

```text
dx/ds =  dA_z/dy
dy/ds = -dA_z/dx
```

  Therefore `A_z` is conserved along the line and the flow is area-preserving.
- The constitutive law `H = dw/dB` is a Legendre transform.  The coenergy
  `w* = H.B - w` satisfies `B = dw*/dH`.
- The tangent `dH/dB` is the Hessian of `w`; convex magnetic energy gives an
  SPD tangent, the same "no fold" condition seen in `weakform_hodge.wls`.
- The hodograph is the Legendre transform of the potential: the `(x,H)` graph
  is Lagrangian, and the dual potential satisfies `x = d chi / dH`.
- The Poisson bracket `{x,y}=1` gives the compact formula `df/ds={f,A_z}`.

Use this topic when connecting flux-line tracking, Chaplygin variables,
coenergy, and the A-formulation into one canonical structure.
"""

SURFACE_DERHAM_HOIBC_RECIPE = """
# Recipe: Surface de Rham + HOIBC / analytic DtN

`surface_derham.wls` splits higher-order impedance boundary conditions into
two independent pieces:

1. Surface topology.  The boundary surface has its own de Rham complex with
   `d o d = 0`.  Global HOIBC modes are surface harmonic 1-forms, counted by
   the surface Betti number `b1`.  A sphere has `b1=0`; an annulus/cylinder
   boundary has `b1=1`.  The count is topological and cannot be changed by a
   material/Hodge weight.
2. Surface geometry.  For the exterior sphere, the exact Dirichlet-to-Neumann
   operator is diagonal in spherical harmonics:

```text
Lambda_n = (n + 1) / R
```

   Kelvin / infinite elements reproduce this nonlocal operator exactly.
   HOIBC is a local Laplace-Beltrami polynomial/Pade approximation: order 1
   matches low modes and then has a knee; higher order matches more modes but
   cannot become exact because `Lambda_n` is not a polynomial in the surface
   Laplacian eigenvalue `n(n+1)`.

Use this before treating HOIBC as "just another material law": the global
loop modes live in the surface de Rham complex, while the impedance values live
in the surface Hodge / Steklov operator.
"""

DTN_GEOMETRY_RECIPE = """
# Recipe: DtN / Steklov operator as condensed exterior Hodge

`dtn_geometry.wls` is the capstone of the differential-geometry suite.  It
answers: what geometric object do Kelvin transforms, infinite elements, and
HOIBC approximate or realize?

It verifies:

- The DtN operator is the boundary condensation (Schur complement) of the
  exterior Hodge-Laplace weak form:

```text
Lambda = A_BB - A_BI A_II^{-1} A_IB
```

- For a sphere, the condensed exterior solution `(R/r)^(n+1)` has conormal
  trace `(n+1)/R`, and the same quantity is the stable fixed point of the radial
  Riccati equation (the infinitesimal thin-shell Schur step).
- The operator is a shifted square root of the surface Laplacian:

```text
Lambda = (1/2 + sqrt(R^2 (-Delta_Gamma) + 1/4)) / R
```

  so it is pseudodifferential and nonlocal.  A finite local HOIBC polynomial can
  match only finitely many modes; the knee is unavoidable.
- Green's identity / Hodge-pairing symmetry makes the DtN self-adjoint and
  positive.  It is the boundary `H^{1/2}` metric, and its spectrum is the
  Steklov spectrum.

Use this before replacing a Kelvin / infinite-element exterior with a local
boundary model, or before interpreting HOIBC coefficients as material constants.
"""

WHITNEY_HEXAHEDRAL = """
# Recipe: Whitney-like basis for hexahedral / polyhedral elements (DGA)

For non-simplicial primal grids (hex, polyhedra), the simple Whitney
formulas don't apply.  Codecasa-Specogna-Trevisan 2010 ("A new set
of basis functions for the discrete geometric approach", JCP 229)
introduces vector basis functions piecewise-uniform on subregions
of the polyhedron, satisfying:

  Property 1:  ∫_{r_j} v^r_i · dr = δ_{ij}      (duality)
  Property 2:  exact reproduction of uniform fields by DoF interpolation
  Property 3:  ∫_v v^r_i dv = ~r_i              (consistency)

The basis functions take the form (for edge type, with dual face ~e_i):

  v^e_i(p) = (~e_j / t_j) δ_{ij}  +  ((I − T_j/t_j) · ~e_i / |v|)
                                                   for  p ∈ s^e_j

where T_j = ~e_j ⊗ e_j, t_j = ~e_j · e_j, and s^e_j is the subregion.

## Wolfram Language input (symbolic version)

```mathematica
(* For a single hexahedron, compute the 12-edge basis matrix v^e *)
(* This is illustrative — full hex needs barycentric subdivision *)

(* Edge vectors e_i (i = 1..12) and dual face vectors ~e_i *)
(* Define for a unit cube *)
edges = {
  {1, 0, 0}, {0, 1, 0}, {0, 0, 1},   (* axis-aligned *)
  ...
};

(* Tensor product helper *)
TensorProd[a_, b_] := Outer[Times, a, b];

(* Constitutive matrix entries via energetic approach *)
constitutiveEntry[i_, j_, mu_, basis_, region_] :=
  Integrate[basis[i] . mu . basis[j], region];

(* The matrix is symmetric, positive definite, and CONSISTENT *)
(* — exactly reproduces uniform fields *)
```

## Why this matters for Radia
Radia C++ uses face-charge extraction on hexahedra with
6 DoFs (one per face).  The mass-matrix construction follows the
same energetic principle: ensure positive-definite stable + consistent
constitutive matrices.  See `rad_interaction.cpp` in the Radia tree
and the Codecasa-Specogna-Trevisan reference in
`differential_forms_bibliography`.

For polyhedra: the basis is **piecewise uniform** (constant per
sub-region), unlike Whitney's polynomial basis on simplices.
"""


def get_mathematica_recipes_documentation(topic: str = "all") -> str:
    """Return Wolfram Language recipes for symbolic verification.

    Topics:
      "all"
      "dsquared"           - rot(grad) = 0 and div(rot) = 0 via Mathematica
      "stokes"             - Stokes' theorem on a tetrahedron face
      "whitney_edge"       - Verify  w_e = w_m ∇w_n - w_n ∇w_m  symbolically
      "maxwell"            - Faraday + ∇·B = 0 from B = rot a, E = -∂_t a - grad ϕ
      "hodge"              - Hodge star ⋆ in R^3 and Minkowski
      "weakform_hodge"     - Committed .wls: weak form as Hodge pairing,
                              pulled-back material/Kelvin weights
      "hodograph"          - Committed .wls: Kelvin/Clebsch/Chaplygin
                              hodograph backbone
      "canonical"          - Committed .wls: Hamiltonian/Legendre reading
                              of the hodograph
      "surface_derham"     - Committed .wls: HOIBC surface de Rham + DtN
      "dtn_geometry"       - Committed .wls: DtN/Steklov as condensed exterior
                              Hodge star / shifted sqrt(-Delta_Gamma)
      "differential_geometry" - Index for the full committed .wls suite
      "curve_surface"      - Committed .wls: curve/surface DG basics for
                              hodograph and stream-function surfaces
      "tex"                - LaTeX (TeXForm) output for paper writing
      "lorentz"            - Lorentz boost of the EM 2-form F
      "hex_dga"            - Codecasa DGA basis functions for polyhedra
    """
    topic = topic.lower().strip()
    if topic in ("dsquared", "d_squared", "d2"):
        return VERIFY_DSQUARED
    if topic == "stokes":
        return VERIFY_STOKES
    if topic in ("whitney_edge", "whitney", "edge"):
        return WHITNEY_EDGE_ELEMENT_FORMULA
    if topic == "maxwell":
        return VERIFY_MAXWELL_IDENTITY
    if topic == "hodge":
        return HODGE_STAR_VERIFICATION
    if topic in ("weakform_hodge", "hodge_material", "pullback_hodge",
                 "transformation_optics", "metric_curvature"):
        return WEAKFORM_HODGE_RECIPE
    if topic in ("hodograph", "clebsch_hodograph", "chaplygin",
                 "chaplygin_hodograph"):
        return HODOGRAPH_WLS_RECIPE
    if topic in ("canonical", "hamiltonian", "legendre", "poisson"):
        return CANONICAL_HODOGRAPH_RECIPE
    if topic in ("surface_derham", "surface_de_rham", "hoibc", "surface_dtn",
                 "vector_hoibc"):
        return SURFACE_DERHAM_HOIBC_RECIPE
    if topic in ("dtn_geometry", "dtn", "steklov", "boundary_metric",
                 "schur_dtn"):
        return DTN_GEOMETRY_RECIPE
    if topic in ("differential_geometry", "diffgeo", "wls", "runbook",
                 "hodograph_runbook"):
        return DIFFERENTIAL_GEOMETRY_WLS_RUNBOOK
    if topic in ("curve_surface", "curve_surface_geometry", "curves_surfaces",
                 "surface_geometry", "streamfunction_geometry",
                 "mathematica_curve_surface_textbook"):
        return CURVE_SURFACE_DG_READING_GUIDE
    if topic in ("tex", "latex"):
        return TEX_OUTPUT_FOR_PAPER
    if topic == "lorentz":
        return LORENTZ_TRANSFORMATION
    if topic in ("hex_dga", "hex", "dga", "polyhedral", "codecasa"):
        return WHITNEY_HEXAHEDRAL
    if topic in ("kelvin", "kelvin_transform", "kelvin_factor", "inversion"):
        return KELVIN_TRANSFORM
    if topic in ("maxwell_stress", "stress", "stress_tensor", "sigma_em"):
        return MAXWELL_STRESS_TENSOR
    if topic == "all":
        return "\n\n".join([
            VERIFY_DSQUARED,
            VERIFY_STOKES,
            WHITNEY_EDGE_ELEMENT_FORMULA,
            VERIFY_MAXWELL_IDENTITY,
            HODGE_STAR_VERIFICATION,
            DIFFERENTIAL_GEOMETRY_WLS_RUNBOOK,
            WEAKFORM_HODGE_RECIPE,
            HODOGRAPH_WLS_RECIPE,
            CANONICAL_HODOGRAPH_RECIPE,
            SURFACE_DERHAM_HOIBC_RECIPE,
            DTN_GEOMETRY_RECIPE,
            CURVE_SURFACE_DG_READING_GUIDE,
            TEX_OUTPUT_FOR_PAPER,
            LORENTZ_TRANSFORMATION,
            WHITNEY_HEXAHEDRAL,
            KELVIN_TRANSFORM,
            MAXWELL_STRESS_TENSOR,
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, dsquared, stokes, whitney_edge, maxwell, hodge, tex, "
        "weakform_hodge, hodograph, canonical, surface_derham, dtn_geometry, "
        "differential_geometry, curve_surface, lorentz, hex_dga, kelvin, "
        "maxwell_stress."
    )
