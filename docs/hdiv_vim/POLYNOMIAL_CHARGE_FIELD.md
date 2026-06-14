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

**Remaining Step-2 work:** (a) **polynomial surface charge** σ (the C++/Python triangle field is exact
for the constant per-face part; the polynomial remainder needs the Graglia linear/higher triangle
field); (b) **curved faces** (flat-triangle only; a curved boundary face needs a curved-element
near-field) — note ngsolve.bem's `LaplaceSL` is curved+triangle but gives the **potential**, not the
field gradient (the `grad G` kernel is a documented ngsolve.bem gap, so it cannot be reused here); (c)
**hex** internal (self-volume ray-trace needs hex face planes; quad boundary faces split into 2
triangles); (d) a **C++ charge-coefficient assembly** (the field-version of the charge Gram, using
`TriField`/`TetField` + the polynomial charge coeffs from `build_charge_gram`'s B). Then **Step 3**
wires the (fast) internal field into the nonlinear `set_field`.

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
