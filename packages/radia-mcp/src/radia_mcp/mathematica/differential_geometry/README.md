# Differential geometry, operational — weak form, Hodge star, hodograph, canonical structure

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
| `weakform_hodge.wls` | weak form = Hodge pairing; `nu' = nu \|det P\|(P^T P)^{-1}` = the scalar weight `W` (the unification); 2-D conformal weight-free, 3-D scale survives; pullback Kelvin `nu'` (spherical isotropic / cylindrical anisotropic); nonlinear `Star_nu` tangent `nu I + nu'(B⊗B)/\|B\|` is **SPD** (elliptic, no fold) | 13 assertions, ALL PASS |
| `hodograph.wls` | the 3-axis backbone: **[1]** the topology/metric split (pullback = cochain map); **[2]** axis 1 — 2-D Kelvin inversion is conformal ⇒ weight-free, complex-potential `(Phi, A_z)` Cauchy–Riemann net; **[3]** axis 2 — Clebsch `B = grad(alpha)×grad(beta)` div-free with `A = alpha grad(beta)`, helicity `A.B = 0` vs ABC Beltrami `A.B = \|B\|^2 != 0` (Moffatt obstruction); **[4]** axis 3 — in 2-D `A_z` IS the Clebsch potential / flux-line Hamiltonian; **[5]** axis 3 deep — the **Chaplygin full nonlinear linearisation**: `div(nu(\|grad A\|) grad A)=0` → the LINEAR self-adjoint elliptic `d/dq((q/mu) A_q) + ((mu q)'/(mu^2 q)) A_thth = 0` (`mu(q)` a coefficient, not a nonlinearity), elliptic ⇔ `(mu q)'>0` ⇔ `nu+nu'\|B\|>0` (no limiting line) | 21 assertions, ALL PASS |
| `canonical.wls` | the **Hamiltonian / canonical reading** of the weak form (weak form = `δ(action)=0`; the Hamiltonian is its Legendre shadow): **[1]** flux lines ARE Hamilton's canonical equations with `A_z` the Hamiltonian (`dx/ds=A_z,y`, `dy/ds=-A_z,x`; `A_z` conserved; area-preserving) — the repo's symplectic flux-line tracking; **[2]** the constitutive law `H=∂w/∂B` is a **Legendre transform** (`B↔H` conjugate, coenergy `w*=H·B-w`, `B=∂w*/∂H`; tangent = Hessian, SPD ⇔ `w` convex ⇔ ⋆_ν "no fold"); **[3]** the hodograph IS the Legendre transform of the potential = a **canonical transformation** (`H=∇Φ` ⇔ the `(x,H)` graph is Lagrangian); **[4]** Poisson bracket `{x,y}=1`, flow `ḟ={f,A_z}` | 12 assertions, ALL PASS |

## Run

```powershell
wolframscript -file weakform_hodge.wls
wolframscript -file hodograph.wls
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
