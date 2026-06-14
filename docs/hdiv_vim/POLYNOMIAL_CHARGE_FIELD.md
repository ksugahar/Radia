# Polynomial-charge field kernel (the order≥2 field)

Status: **Step 1 done** (external field, `reconstruct_field_polynomial`, golden-locked,
**tetrahedral AND hexahedral**); Steps 2–3 designed, not yet implemented.

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
| 2 | Internal / near points (the field at the body's own quadrature points) | yes (`1/r²` at `r'→r`) | designed: analytic Wilton/PhiTet base + smooth-remainder Gauss-Duffy (mirror the self-Gram subtraction), likely C++ |
| 3 | Wire Step 2 into the per-element nonlinear Newton (`set_field` ⇐ polynomial field, not `M_mass⁻¹ N m`) | — | designed: genuine order ≥ 2 nonlinear M, golden vs a finer-mesh / Radia MMM reference |

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
