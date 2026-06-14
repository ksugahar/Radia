# Polynomial-charge field kernel (the order≥2 field)

Status: **Step 1 done** (external field, `reconstruct_field_polynomial`, golden-locked,
**tetrahedral AND hexahedral**); **the analytic charge-field kernel is done at degree 0, 1, 2 (closed
form) AND arbitrary degree (general assembler)** — `flat_/linear_/quadratic_triangle_charge_field`,
`tet_volume_field_linear/quadratic`, and the general `polynomial_triangle_charge_field` /
`tet_volume_field_polynomial` (validated vs Gauss for cubic/quartic, machine precision; 27 golden tests).
Step 3 (wiring the fast internal field into the nonlinear `set_field`) and the C++ port remain.

The kernel is **element-type agnostic**: the quadrature points / weights / normals come from NGSolve's
own `mesh.GetTrafo` + `IntegrationRule(el.type)` + `specialcf.normal`, so the *same* code handles tet
**and hex** (and prism) meshes, flat or curved (`mip.measure` carries the curved Jacobian,
`specialcf.normal` the curved outward normal). A hex box and a tet box of the same body give the same
external field to ~machine precision (`test_hex_matches_tet_*`).

## Why

A genuine order ≥ 2 nonlinear HDiv-VIM solve needs the magnetic field of a **polynomial**
magnetization `M(x)` (HDiv order *p*) — both for the engineering deliverable (the stray field
around a soft-iron part) and for the constitutive law `M = χ(H) H` inside the body.

The committed `reconstruct_field` uses the per-element **centroid** `M` (piecewise-constant), i.e.
the **surface charge `σ = M·n` only** — it silently **drops the volume charge `ρ = -div M`**. That
is exact for uniform `M` (`div M = 0`) but, wherever `div M ≠ 0`, omitting `ρ` is a **90–230 % error**
(measured, `tests/feec/test_hdiv_vim_poly_field.py::test_volume_charge_is_essential_for_div_M`).

This is also why the old order ≥ 2 nonlinear solve diverged: `M_mass⁻¹ N m` (the weak demag field)
has a solenoidal nullspace at order ≥ 2, and even the centroid field reconstruction misses `ρ`.

## The kernel

The field of a magnetization is the field of its magnetic charges (`H = -∇φ_M`):

```
ρ = -div M   (volume charge, L2 order p-1)
σ =  M·n     (surface charge, SurfaceL2 order p)

           1   ⎡        r-r'                        r-r'             ⎤
H(r) =  ─────  ⎢ ∫_V ρ ─────── dV'  +  ∫_S σ ─────────── dS' ⎥
         4π    ⎣       |r-r'|³                   |r-r'|³            ⎦
```

These are exactly the charges the HDiv-VIM already forms in `build_charge_gram` (the `B` map:
`ρ = -div M ∈ L2(p-1)`, `σ = M·n ∈ SurfaceL2(p)`). The charge **Gram** `G` is the charge–charge
**energy** `∫∫ q q'/|r-r'|`; this kernel is its **field-at-a-point** companion
`∫ q (r-r')/|r-r'|³`. The C++ already has the *constant*-charge potential building blocks
(`_hdiv_phi_tet`, `_hdiv_tri_potential` = Wilton); the field/polynomial generalisation is the work.

## Staging

| Step | Scope | Singular? | Status |
|------|-------|-----------|--------|
| **1** | **External** points (stray field of a polynomial-M body), **tet + hex** | no (`r` clear of the body) | **done** — `reconstruct_field_polynomial`, element-agnostic Python reference, golden-locked |
| **2** | Internal / near points (the field at the body's own quadrature points), **tet** | yes (`1/r²` at `r'→r`) | **assembled** — `reconstruct_field_internal` (self-volume spherical + far-volume + analytic surface), golden-locked; polynomial surface σ / curved faces / C++ remain |
| 3 | Wire Step 2 into the per-element nonlinear Newton (`set_field` ⇐ polynomial field, not `M_mass⁻¹ N m`) | — | designed: genuine order ≥ 2 nonlinear M, golden vs a finer-mesh / Radia MMM reference |

## Step 2 — internal/near singular field

For a query point `r` inside element `e`, split the charge sum by proximity:

- **far** (elements/faces not containing or adjacent to `r`): non-singular → the Step-1 quadrature.
- **self / near** (the element holding `r`, and `r`'s own faces): singular → handled analytically.

**Kernel A — self-element VOLUME charge, spherical ray-trace** (`tet_self_volume_field`).
Substituting `r' = r + s·ŝ` (`dV' = s² ds dΩ`, `r-r' = -s ŝ`, `|r-r'|³ = s³`) cancels the kernel:

```
                 1                     ⌠       ⌠ smax(ŝ)
H_self(r) = ───── INT_e ρ ... dV' = - ─── ⎮  ŝ ⎮  ρ(r+s ŝ) ds dΩ      (NON-singular)
                4π                    4π ⌡S² ⌡0
```

`smax(ŝ)` = ray distance from `r` to the element boundary; the inner `∫ρ ds` is closed-form for a
polynomial `ρ`. Same substitution as the self-energy spherical method. Golden
(`test_tet_self_volume_field_vs_phitet_gradient`): constant `ρ` on the unit tet vs `-(ρ/4π)∇(phi_tet)`
(central FD of the exact analytic Newtonian potential `_hdiv_phi_tet`) at three interior points → **rel
7e-4 … 3e-3**.

**Kernel B — surface charge, analytic uniform-triangle field** (`flat_triangle_charge_field`).
The exact `INT_T (r-r')/|r-r'|³ dS'` for a flat triangle (Wilton/Graglia: solid-angle normal term +
per-edge log tangential term), valid at any `r` (near/far/on-face PV). Golden
(`test_flat_triangle_charge_field_exact`): matches a fine Gauss reference to ~machine precision.

**Assembly — `reconstruct_field_internal` (DONE, tet).** Per obs point `r`: locate the self tet
(barycentric), self-volume via Kernel A, far-volume via Gauss-Duffy (`r` outside those elements →
non-singular), surface via Kernel B over all boundary triangles (constant σ = M·n per face). Golden
`test_internal_field_assembly_uniform_sphere`: uniform sphere → **center = −M/3** (pins the assembly
factors) and **near-surface (0.95R) the analytic surface is ~24× better than Step-1** (Step-1 plain
Gauss-Duffy is ~58 % off there). The residual ~2.4 % at 0.95R is **faceting** (the flat mesh's
near-surface field genuinely differs from the smooth −M/3; Kernel B gives the faceted body's field
*exactly*), removed only by curving.

**C++ kernels (DONE).** The analytic constant-charge field is now in C++ (`rad_hdiv::TriField`,
`rad_hdiv::TetField`; probes `_hdiv_tri_field` / `_hdiv_tet_field`), exact near AND far, NO quadrature:
- `TriField` = the Wilton triangle field (`-grad TriPotential`) — matches Python
  `flat_triangle_charge_field` to **machine precision** (`test_cpp_tri_field_matches_python`).
- `TetField` = the tet volume-charge field (`-grad PhiTet` via the divergence theorem,
  `0.5 Σ_faces[n·TriPotential + d·TriField]`) — matches `-grad(_hdiv_phi_tet)` (FD) to **~1e-9** and the
  Python spherical ray-trace to ~1e-3 (`test_cpp_tet_field_matches_grad_phitet`).
This is the speed enabler for a practical order≥2 nonlinear solve (the Python reference field is too
slow in the loop — the self-volume sphere-integral × mesh-location is prohibitive).

## The polynomial volume-charge field — degree 1 (DONE, closed form)

The order-2 volume charge `ρ = −div M` is **linear** per cell, and its field is now a **closed form**,
exact to machine precision at any point (interior, surface, exterior), with **no quadrature** — the
exact, ~orders-faster replacement for the ~1e-3 spherical `tet_self_volume_field`. The key identity is
`(r-r')/R³ = ∇'(1/R)`, so by the product rule + divergence theorem the field of a polynomial volume
charge reduces to **lower-degree potential integrals** (one differential order down):

```
∫_V ρ (r-r')/R³ dV'  =  Σ_faces n_f ∫_face ρ/R dS'  −  ∫_V (∇ρ)/R dV'
```

For linear `ρ = ρ0 + g·r'` (`∇ρ = g` const) this is `Σ_f n_f [ρ0 I0_f + g·M1_f] − g·PhiTet`, needing only:

- `triangle_potential_const` — `I0 = ∫_T 1/R dS'` (Wilton; pure-Python, **bit-identical** to the C++ `_hdiv_tri_potential`);
- `triangle_potential_moment` — `M1 = ∫_T r'/R dS'`, first moment via the **surface** divergence theorem `∫_T (r'−r_p)/R dS' = Σ_edges m_e ∫_edge R dl` (closed-form edge integrals);
- `tet_newtonian_potential` — `PhiTet = ∫_V 1/R dV' = −½ Σ_f h_f I0_f` (from `1/R = ½∇'²R`; pure-Python, matches the C++ `_hdiv_phi_tet` to machine precision).

`tet_volume_field_linear(verts, r, rho0, grho)` assembles these. Validated four ways
(`tests/feec/test_hdiv_vim_poly_field.py`): the pure-Python building blocks vs the C++ probes
(**machine**), the constant case vs the **independently-derived C++ `TetField`** (`−grad PhiTet`,
**1e-12**), the linear case vs far tet Gauss (**1e-10**), and the interior linear case vs the spherical
ray-trace (**~1e-3**, the spherical method's own accuracy — confirming the closed form *is* the exact
value the spherical method converges to). All building blocks are pure-Python (no debug-probe runtime
dependency); the C++ probes are used only as test oracles.

## The polynomial surface-charge field — degree 1 (DONE, closed form)

The order-2 surface charge `σ = M·n` is **linear** per face. Since `(r-r')/R³ = −∇_r(1/R)`, the field
of a linear σ is exactly `−∇_r φ_σ` with `φ_σ = ∫_T σ/R dS' = σ0·I0 + s·M1` — the **degree-1 triangle
potential** we already have. Differentiating in closed form:

```
∫_T (σ0 + s·r')(r-r')/R³ dS'  =  (σ0 + s·r_p) F_const  −  Σ_edges (s·m_e) G_e  −  I0 · s_∥
```

needing only `F_const` = `flat_triangle_charge_field` (the constant-σ field), `I0` =
`triangle_potential_const`, `s_∥` = in-plane part of `s`, and one new elementary building block
`G_e = ∫_edge (r-r')/R dl` (`_edge_field_integral`, closed-form `asinh`/`sqrt`). `linear_triangle_charge_field`
assembles these — validated vs off-plane Gauss to **machine precision** (`test_linear_triangle_charge_field_vs_gauss`),
and `s = 0` reproduces `σ0·flat_triangle_charge_field` bit-identically. This is the surface companion of
`tet_volume_field_linear` — together they are the **complete degree-1 (linear) charge field**, exact and
closed-form for both the volume `−div M` and the surface `M·n` terms.

## The polynomial volume-charge field — degree 2 (DONE, closed form)

The quadratic volume charge `ρ = ρ0 + g·r' + r'ᵀQr'` (Q symmetric) field, via the same
divergence-theorem recursion `∫_V ρ(r-r')/R³ = Σ_f n_f ∫_face ρ/R − ∫_V (∇ρ)/R`:

```
∫_V ρ (r-r')/R³ dV'  =  Σ_faces n_f [ρ0 I0_f + g·M1_f + Q:M2_f]  −  (g·PhiTet + 2 Q·V1)
```

adds two degree-2 moment building blocks, each from the **same identities one degree up**:

- `triangle_potential_moment2` — surface second moment `M2 = ∫_T r'⊗r'/R dS'`, from the Hessian identity
  `∇'_s∇'_s(R³) = 3(ξ⊗ξ/R + R·P)` ⟹ `∫_T ξ⊗ξ/R = Σ_e (∫_edge R ξ dl)⊗m_e − P·∫_T R dS'` (with
  `∫_T R dS' = ⅓[Σ_e m_e·∫_edge R ξ dl + h² I0]`), then shifted by `r_p`. Symmetric, exact.
- `tet_newtonian_moment` — volume first moment `V1 = ∫_V r'/R dV' = ⅓[r·PhiTet − Σ_f h_f·M1_f]`, from
  `1/R = ½∇'²R` weighted by `r'_k`.

`tet_volume_field_quadratic` assembles these — validated vs far tet Gauss to **machine precision**
(`test_tet_volume_field_quadratic_vs_gauss_far`), `M2`/`V1` each vs Gauss to machine precision, and
`Q = 0` reduces to `tet_volume_field_linear` bit-identically.

## The polynomial surface-charge field — degree 2 (DONE, closed form)

The quadratic surface charge `σ = σ0 + s·r' + r'ᵀSr'` (S symmetric) field, via the systematic
in-plane/normal split `(r-r')/R³ = ∇'_s(1/R) + h·n/R³`:

```
∫_T σ(r-r')/R³ dS'  =  [Σ_e m_e ∫_edge σ/R dl − (P·s·I0 + 2·P·S·M1)]  +  h·n·[σ0 J3_0 + s·J3_1 + S:J3_2]
                        └──────────── in-plane ──────────────┘            └─────── normal ────────┘
```

with the `1/R³` moments `J3_0 = ∫_T 1/R³ = (n·F_const)/h` (reusing the constant-σ field for the solid
angle), `J3_1 = ∫_T r'/R³`, `J3_2 = ∫_T r'⊗r'/R³` (from `ξ⊗ξ/R³ = P/R − Hess_s R`), plus the quadratic
edge integrals `∫_edge l^k/R dl` (`_edge_monomial_over_R`, closed-form). `quadratic_triangle_charge_field`
assembles these — validated vs off-plane Gauss to **machine precision**
(`test_quadratic_triangle_charge_field_vs_gauss`), and `S = 0` matches `linear_triangle_charge_field`
(two independent derivations — in-plane/normal vs `−∇φ` — agreeing). It subsumes the constant and linear
cases. Together with `tet_volume_field_quadratic` this is the **complete degree-2 (quadratic) charge
field**, exact and closed-form for both terms.

## Arbitrary degree (DONE, the general assembler)

The degree-0/1/2 closed forms are fast, hand-derived special cases. The **general assembler**
(`polynomial_triangle_charge_field`, `tet_volume_field_polynomial`) handles **any polynomial degree**
via the general moment recursion, and reduces to the closed forms at degree ≤ 2 (verified — two code
paths agree to machine precision).

**Surface moments** `A_k = ∫_T ξ^⊗k/R`, `B_k = ∫_T ξ^⊗k/R³` (in an in-plane basis) from the master
recursion `(2+p+k)∫_T ξ^α R^p − p h² ∫_T ξ^α R^{p-2} = ∮ ξ^α(ξ·m) R^p dl`:
- `A_k = (E⁻¹_k − h² B_k)/(k+1)` with the gradient relation `B_k = A_{k-2}(deriv) − edge`, so **`A_k`
  comes from `A_{k-2}` + edge integrals alone** — and it is **h-safe** (fold `h² B_k = h²[(a-1)A_{k-2}
  − edge]`, finite even when `r` lies in a face plane). `triangle_inplane_moments(P, r, degree)`.
- Edge integrals `∫_edge l^n/R dl` via the `∫u^n/R` reduction formula (any `n`).

**Surface field** = the in-plane/normal split, contracting the charge's in-plane monomial coefficients
(extracted by sampling at barycentric nodes — exact for a degree-`d` polynomial) with `A_k`, `B_k`, and
the edge integrals.

**Volume potential moments** `∫_V r'^α/R = 1/(|α|+2)[ −Σ_f h_f ∫_face r'^α/R + Σ_i r_i α_i ∫_V
r'^{α-e_i}/R ]` (from `1/R = ½∇'²R` + Euler), bottoming at `PhiTet`, reducing to **surface** potentials.
**Volume field** = the divergence-theorem recursion `Σ_f n_f ∫_face ρ/R − ∫_V (∇ρ)/R`.

Validated vs Gauss to **machine precision** for cubic (surface + volume) and quartic (surface), and
against the degree-2 closed forms. All pure-Python, no debug-probe runtime dependency.

**Remaining work:** (a)
**curved faces** (flat-triangle only) — note ngsolve.bem's `LaplaceSL` is curved+triangle but gives the
**potential**, not the field gradient (`grad G` is a documented ngsolve.bem gap, so it cannot be reused
here); (b) **hex** internal; (c) a **C++ port** of the moment building blocks + a
**charge-coefficient assembly** (sum the kernels weighted by the polynomial charge coeffs from
`build_charge_gram`'s B, H-matrix accelerated). Then **Step 3** wires the (fast) internal field into the
nonlinear `set_field`.

## Step 1 — validated (`reconstruct_field_polynomial`)

- **Uniform-M sphere** (`div M = 0`): center `H = -M/3` to 1.2e-3; external = analytic dipole to ~6 %
  (= the flat-mesh faceting at `h=0.4`, removed by `mesh.Curve` — *not* a kernel error).
- **Linear M** (`M = (0,0,M₀(1+z))`, `div M = M₀`): the full kernel is coarse→fine self-convergent,
  while dropping `ρ` (surface-only) is 90–230 % wrong — the volume-charge term is essential.

Reference (Python) implementation: charges sampled once over (elements × Gauss-Duffy) and summed
vectorised over the observation points. Cost is independent of the number of observation points; a
C++/H-matrix-accelerated version is part of Step 2.

## Notes / pitfalls (verify-first record)

- Use NGSolve's own geometry, not a hand-rolled affine map: `trafo = mesh.GetTrafo(ElementId(...))`,
  `for ip in IntegrationRule(mesh[ei].type, order)`, `mip = trafo(ip)`. Physical point = `mip.point`,
  physical quadrature weight = `ip.weight * mip.measure` (the MIP exposes `.point`, `.measure`,
  `.jacobi` — **not** `GetMeasure()` / `GetJacobiDet()` / `.weight`). This makes the kernel handle
  tet **and** hex (and curved) for free.
- Surface charge `M·n`: build it as a CF `InnerProduct(gfM.Trace(), specialcf.normal(dim))` and
  evaluate at the boundary MIP — `specialcf.normal` is the correct outward (and curved) normal; do not
  hand-roll the face normal from vertices (wrong sign / no curving).
- `div(gfM)(mip)` may return a 1-tuple — extract the scalar.
- The field integrand is `1/r²`-singular, so Step 1 is **external only**; an internal/near point needs
  the Step-2 singular-aware kernel. Do not use Step 1 inside the body.
