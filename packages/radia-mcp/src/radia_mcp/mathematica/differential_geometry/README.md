# Differential geometry, operational — weak form, Hodge star, hodograph, canonical structure, HOIBC/DtN

Symbolic (Mathematica / `.wls`) reference for the **lab policy**:

> Use differential geometry to make two things easy to SEE —
> **(1)** the FEM **weak (energy) form** is a **Hodge pairing**, and
> **(2)** the **Hodge star carries the material constant**, so a coordinate
> map (mesh deformation / Kelvin open boundary / a hodograph chart) **modulates
> the material** as an explicit weight.

This is the differential-geometry **half** of `examples/clebsch_hodograph/`
(the spine is [`docs/clebsch_hodograph/HODOGRAPH_BACKBONE.md`](../../../../../../docs/clebsch_hodograph/HODOGRAPH_BACKBONE.md)),
re-expressed symbolically. Companion to the FEEC shape-function bases in
[`../basis_functions/`](../basis_functions/) (the de Rham complex those weak
forms are assembled on).

**Prose synthesis (the one picture):**
[`docs/differential_geometry/PREMETRIC_PICTURE.md`](../../../../../../docs/differential_geometry/PREMETRIC_PICTURE.md)
— the single `d` (topology / basis) vs `⋆` (geometry / material / DtN) split,
premetric EM, with the volume × exterior × surface master table and the
`.wls`-evidence map for weak form / Kelvin / IE / PML / HOIBC / hodograph.

## The one idea

Magnetostatics is a **closed 2-form with a Hodge-star constitutive law**:

```
  B in Omega^2 ,  dB = 0          (de Rham: B = dA, "no monopoles")   [TOPOLOGY]
  H in Omega^1 ,  H = Star_nu B   (the material Hodge star carries nu) [METRIC]
```

Under **any** diffeomorphism `phi` the two halves behave completely differently
— and that split is the whole point:

| half | transforms | consequence |
|---|---|---|
| `dB = 0` (topology) | **naturally** — the pullback is a cochain map (`kStar d = d kStar`), no metric data | preserved by *any* `phi` for free |
| `H = Star_nu B` (metric) | the Hodge star pulls back with a **weight** | this is the **only** thing that moves |

The weight is one formula, used twice:

```
  scalar Dirichlet (0-form) energy weight :  W   =      |det J| (J^T J)^{-1}      (J = coord Jacobian)
  material constant modulation            :  nu' = nu |det P| (P^T P)^{-1}        (P = pullback on the orthonormal 1-form basis)
```

`nu' = nu |det P| (P^T P)^{-1}` is derived here (energy invariance + the 2-form
cofactor pullback) and reproduces the **Nagamine–Yamaguchi–Sugahara pullback
Kelvin law** (CEFC 2026):

- **spherical** inversion (conformal) → `nu' = (r'/R)^2 nu` **isotropic**;
- **cylindrical** inversion (non-conformal) → `nu' = diag(1,1,(rho'/R)^4) nu` **anisotropic**
  (anisotropy is the *signature* of non-conformality).

## Files

| file | content | self-test |
|------|---------|-----------|
| `weakform_hodge.wls` | weak form = Hodge pairing; `nu' = nu \|det P\|(P^T P)^{-1}` = the scalar weight `W` (the unification); 2-D conformal weight-free, 3-D scale survives; pullback Kelvin `nu'` (spherical isotropic / cylindrical anisotropic); nonlinear `Star_nu` tangent `nu I + nu'(B⊗B)/\|B\|` is **SPD** (elliptic, no fold); **metric ≠ curvature** — the Kelvin & polar metrics are FLAT (`Riemann = 0`; a map can't curve flat space), sphere control genuinely curved (`R = 2`); **transformation optics = premetric "push everything onto the Hodge"** — our `ν' = χ(g) = √(det g) g⁻¹` is Pendry's `ε' = ΛΛᵀ/det Λ` (the one metric→material map, on the pullback vs pushforward metric), and any coordinate map (conformal or not) keeps `Riemann = 0`: the geometry moves onto the Hodge, the topology can't | 20 assertions, ALL PASS |
| `hodograph.wls` | the 3-axis backbone: **[1]** the topology/metric split (pullback = cochain map); **[2]** axis 1 — 2-D Kelvin inversion is conformal ⇒ weight-free, complex-potential `(Phi, A_z)` Cauchy–Riemann net; **[3]** axis 2 — Clebsch `B = grad(alpha)×grad(beta)` div-free with `A = alpha grad(beta)`, helicity `A.B = 0` vs ABC Beltrami `A.B = \|B\|^2 != 0` (Moffatt obstruction); **[4]** axis 3 — in 2-D `A_z` IS the Clebsch potential / flux-line Hamiltonian; **[5]** axis 3 deep — the **Chaplygin full nonlinear linearisation**: `div(nu(\|grad A\|) grad A)=0` → the LINEAR self-adjoint elliptic `d/dq((q/mu) A_q) + ((mu q)'/(mu^2 q)) A_thth = 0` (`mu(q)` a coefficient, not a nonlinearity), elliptic ⇔ `(mu q)'>0` ⇔ `nu+nu'\|B\|>0` (no limiting line) | 21 assertions, ALL PASS |
| `canonical.wls` | the **Hamiltonian / canonical reading** of the weak form (weak form = `δ(action)=0`; the Hamiltonian is its Legendre shadow): **[1]** flux lines ARE Hamilton's canonical equations with `A_z` the Hamiltonian (`dx/ds=A_z,y`, `dy/ds=-A_z,x`; `A_z` conserved; area-preserving) — the repo's symplectic flux-line tracking; **[2]** the constitutive law `H=∂w/∂B` is a **Legendre transform** (`B↔H` conjugate, coenergy `w*=H·B-w`, `B=∂w*/∂H`; tangent = Hessian, SPD ⇔ `w` convex ⇔ ⋆_ν "no fold"); **[3]** the hodograph IS the Legendre transform of the potential = a **canonical transformation** (`H=∇Φ` ⇔ the `(x,H)` graph is Lagrangian); **[4]** Poisson bracket `{x,y}=1`, flow `ḟ={f,A_z}` | 12 assertions, ALL PASS |
| `surface_derham.wls` | **HOIBC** read via **topology + the analytic DtN**: **[A]** the boundary's surface de Rham (`d∘d=0`); HOIBC global modes = surface harmonic 1-forms = surface Betti `b1` (sphere `b1=0` ⇒ no global loop; cylinder `b1=1` ⇒ one loop; the count is topological, the impedance values geometric); **[B]** the exact exterior DtN ladder `Λ_n=(n+1)/R` (diagonal in `Y_n`; Kelvin/IE reproduce it exactly), and HOIBC = a **local Laplace–Beltrami-polynomial (Padé) approximation** — order-1 matches `Λ_n` at `n=0,1`, deviates at `n=2` (the knee); order-2 at `n≤2`, deviates at `n=3` (higher order = more matched modes, never exact, because `Λ(n)=(n+1)/R` is not a polynomial in the Laplace–Beltrami eigenvalue `n(n+1)`); **[C]** the **vector** case (cf. act7_30) — the tangential field splits (surface Hodge–Helmholtz) into `grad_Γ α` + `curl_Γ β`, with **two distinct de Rham Steklov ladders**: H(curl) tangential `n/R` (derived from the toroidal field's curl) vs 0-form/H(div) normal `(n+1)/R` (differ by `1/R`), so a vector HOIBC is a 2-channel (tensor) impedance, each a local Padé of its ladder | 12 assertions, ALL PASS |
| `dtn_geometry.wls` | the **DtN / Steklov operator** in DG language — the exterior Hodge `⋆` **condensed** onto the boundary: **[1]** DtN = **boundary condensation (Schur complement)** of the bulk Hodge–Laplace form — the conormal trace of the condensed exterior solution `(R/r)^{n+1}`, and the decaying fixed point of the radial **Riccati** `m'=m²-(2/r)m-n(n+1)/r²` (the thin-shell Schur step); **[2]** DtN = a **shifted square root of the surface Laplacian** `Λ=(½+√(R²(-Δ_Γ)+¼))/R = (n+1)/R`, pseudodifferential / **nonlocal** ⇒ HOIBC's local-polynomial knee is unavoidable (the `√` is not a polynomial); **[3]** **self-adjoint & positive** (Green's identity = Hodge-pairing symmetry + Legendre orthogonality; the DtN is the boundary `H^{1/2}` metric, its spectrum the Steklov spectrum) | 9 assertions, ALL PASS |

The **dual half** of "push everything onto the Hodge" — *topology is NOT in the
Hodge* — is verified in [`../basis_functions/cohomology.wls`](../basis_functions/cohomology.wls):
the harmonic-1-form count (= Betti `b1`) is **invariant under any SPD metric**
(Hodge theorem), so a coordinate map / material / Hodge star moves the geometry
but never the holes (`dB=0` survives; only `H=⋆B` moves).

## Run

```powershell
wolframscript -file weakform_hodge.wls
wolframscript -file hodograph.wls
wolframscript -file canonical.wls
wolframscript -file surface_derham.wls
wolframscript -file dtn_geometry.wls
```

(or via the `mathematica_evaluate` MCP tool / `radia_mcp.mathematica`). Each file
is self-contained and self-testing (`ALL PASS` on success), in the same clean
`.wls` discipline as `../basis_functions/`.

## Source of truth / prior art

- **Weak form / Hodge / pullback weight**: `docs/clebsch_hodograph/HODOGRAPH_BACKBONE.md` §1–2;
  `examples/clebsch_hodograph/bidirectional_coordinate_transform_2d.py`
  (`symbolic_pullback_check`, golden) — these `.wls` are the Mathematica twin/extension.
- **Pullback Kelvin material law**: H. Nagamine, T. Yamaguchi, K. Sugahara,
  "A Pullback-Based Formulation of Kelvin Transformation in Electromagnetic Field
  Analysis," CEFC 2026 (lineage Wong–Ciric 1985, Freeman–Lowther 1988, Sugahara 2022).
- **Clebsch / helicity**: P. Robert, *IEEE Trans. Magn.* 27(5) 1991; H. K. Moffatt,
  *J. Fluid Mech.* 35, 1969.
- **Chaplygin / von Mises hodograph**: classical gas dynamics (the field-plane
  linearization and the limiting line the magnetic `Star_nu` does NOT have).

These are classical; the files are a working symbolic map of the machinery, not an
originality claim (see HODOGRAPH_BACKBONE.md §4).
