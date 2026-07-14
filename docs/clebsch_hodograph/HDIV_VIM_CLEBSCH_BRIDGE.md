# The HDiv-VIM ↔ Clebsch bridge (the de Rham capstone) and the saturation frontier

This note maps how the two research lines in `docs/clebsch_hodograph/demos/` —
the **HDiv-VIM demag solver** (`radia.vim`, the FEEC `H(div)` element Radia
is migrating to) and the **Clebsch hodograph** design line — are one structure
seen from two sides, and pins down exactly what is *verified* versus what is the
*open* hard part.

The verified linear (kinematic) level is
[`docs/clebsch_hodograph/demos/hdiv_vim_clebsch_loopstar.py`](demos/hdiv_vim_clebsch_loopstar.py)
(golden `test_hdiv_vim_clebsch_loopstar`).

## One structure, four faces

| face | object | role |
|---|---|---|
| **operator** | HDiv-VIM `N = BᵀGB` | discrete demag; loop modes field-null |
| **potential** | Clebsch `(α,β)`, `M = ∇α×∇β` | coordinates of the field-null subspace |
| **conditioning** | convex B-input energy `∫W(\|B\|)` | well-posed where the M-/H-Picard is not |
| **2-D reduction** | Chaplygin hodograph | *linearises* the saturation |

### 1. Operator ↔ potential (verified)

A magnetization `M` is a **2-form** (flux). Its magnetic charge is `ρ = -div M
= -dM`. The HDiv-VIM demag operator `N = BᵀGB` (`B` = the charge map
`M ↦ (ρ, σ=M·n)`) annihilates every `M` with `dM = 0` (closed 2-forms) — *"loop
modes are field-null by construction"*. A **Clebsch** field
`M = ∇α×∇β = d(α dβ)` is an **exact**, hence closed, 2-form ⇒ charge-free ⇒
field-null. So **the loop modes ARE the Clebsch fields**, and the Hodge split

```
M  =  ∇φ        (+)   ∇α × ∇β
      "star"           "loop" = Clebsch
      carries charge   charge-free → field-null
      pole-forming     flux-guiding (yoke return)
```

is the loop–star decomposition the VIM is built on. Verified on a unit sphere
(HDiv order 1): `D(Clebsch (y,−x,0)) = 1.6e-3 ≈ 0` (‖div M‖ = 1e-14, machine
zero) vs `D(gradient) = 0.997`, `D(uniform) = 1/3`; the Clebsch field makes
~no external field; and adding `t·Clebsch` to a charged magnetization changes
the external field by `< 1.3 %` (gauge invariance).

### 2. Conditioning: the convex B-input form is the de Rham dual (verified in pieces)

The nonlinear demag solve is **ill-conditioned at the saturation knee in the
H-/M-input form** — a naive field-exact Picard on the magnetization *diverges*
(spectral radius ≥ 1; reproduced here), and the reduced-Ω `μ(|H|)` Picard finds
spurious fixed points (`saturation_loop_2d.py`). The de Rham/Hodge structure says
*which* formulation is convex: the one where the **flux 2-form `B`** is primary
and the energy is `∫W(|B|)` with `W` convex (`ν(|B|)` monotone). That is the
**B-input A-formulation — the de Rham dual of the HDiv-VIM** — and it converges
*to machine precision* (`saturation_loop_2d.py`, the `nu(|B|)` A-formulation).
The well-conditioned weak form `A⁺ = (1/χ)M_mass + N_weak` used in the HDiv-VIM
linear solve is the same idea: solve in the flux/charge pairing, not the M
fixed point.

### 3. 2-D reduction: the hodograph linearises saturation (verified) — and A_z IS the Clebsch potential

In 2-D, `B = ∇A_z × ẑ`, i.e. **`A_z` is exactly the Clebsch potential `α`**
(with `β = z`). So the **Chaplygin hodograph** that linearises the saturable
`A_z` problem (`chaplygin_hodograph_2d.py`: a 1-shot quadrature reproducing the
full nonlinear loop) *is* the linearisation of the HDiv-VIM **loop** saturation.
The 2-D frontier is therefore essentially already in hand — it had simply not
been *named* as the HDiv-VIM connection. In 2-D the single Clebsch scalar `A_z`
+ the interchange `(x,y) ↔ (θ, q=|B|)` turns `div(ν(|∇A_z|)∇A_z)=0` into a
**linear** variable-coefficient Chaplygin equation.

## The dynamical face — flux-line closure (Noguchi / Sugahara 2020)

A flux line is an integral curve of `B` (`dx/ds = B`); in 2-D this **is** Hamilton's
equations with **`A_z` as the Hamiltonian**, so `A_z` (= the Clebsch potential) is
conserved and flux lines close. `flux_line_closure_symplectic.py` isolates the two
requirements: the field must be a **closed 2-form** (`div B = 0` — the de Rham /
edge-`H(curl)` requirement, Noguchi) *and* the integrator must be **symplectic**
(`A_z`-conserving — accelerator beam tracking, Sugahara 2020). `A_z` = Clebsch
potential = flux-line-flow Hamiltonian is the single object behind both the field
(FEEC) face and the dynamical (tracking) face. This connects the solver's field
reconstruction quality directly to whether its flux lines close — a leaky
(`M_mass⁻¹ N m`) reconstruction is *not* a closed 2-form and its flux lines spiral.

That diagnostic is made concrete on a **real solved FE field** in
`flux_line_realfield_ngsolve.py` (golden `test_flux_line_realfield_ngsolve`): a 2-D
magnetostatic solve `−div(grad A_z) = J`, one flux line traced for three
reconstructions of the *same* solve with the *same* RK4 integrator. The de Rham
field `rot(grad A_z)` (the edge-FE `B = curl A`) is exactly tangent to the flux
surfaces (`B·∇A_z = 0`, misalignment `0.0`) and its line **closes**; a
nodal-averaged reconstruction (`3.7e-2`) and an explicit charge admixture
`B + ε∇A_z` (`6e-2`) both leak off the flux surface and **spiral** (`A_z` drift
`21×` / `~2000×` worse). So *"trace a flux line; if it spirals, the reconstruction
leaked solenoidal content"* is a directly usable field-quality test for the
HDiv-VIM migration.

**Order vs representation** (`derham_closure_order_sweep.py`, golden
`test_derham_closure_order_sweep`) answers two old questions: can Noguchi's de Rham
flux line be *extended* to a symplectic tracker, and was Kameari's "2nd-order elements
don't close" a de Rham issue? Measured: `B = curl A` is divergence-free for **any**
conforming `A` — edge `H(curl)` *or* nodal `[H1]³` — at **every order** (`~1e-15`),
so "edge vs nodal `A`" is *not* the closure discriminator. The closure-breaker is
*leaving* the de Rham representation — nodally **smoothing** `B` — whose spurious
divergence falls with order but is **never zero** (`1.2e-2 → 1.3e-3 → 6e-5` at
`p = 1,2,3`): the flux lines do not close *even at 2nd order*, while the native
de Rham `B` closes at every order. Closure is governed by the **representation**, not
the order; de Rham is the closed-2-form precondition that makes the symplectic /
volume-preserving tracker meaningful (the Noguchi extension; Bossavit / Nédélec edge
elements).

## The open frontier (3-D helicity) — now characterised

In 3-D the Clebsch representation needs **two** potentials `(α,β)` for three
coordinates, and a *global* Clebsch pair need not exist — the obstruction is the
**helicity** `h = ∫A·B` (the topological linking of field lines, Moffatt 1969).
This is now **verified** (`clebsch_3d_closing_condition.py`): a Clebsch field is
helicity-free pointwise (`h = −7e-17`), while the ABC Beltrami field has
`h = 3(2π)³ ≠ 0` and chaotic, never-closing flux lines (Poincaré occupancy 0.55).
So the 3-D closing condition is **`h = 0` ⟺ a global Clebsch pair exists ⟺ flux
lines lie on flux surfaces**, and the 3-D saturable HDiv-VIM solve admits the
Chaplygin/Clebsch linearisation **only where the helicity vanishes** (integrable /
foliated fields). The concrete open questions:

1. **Conditioning (tractable next rung).** Does solving the saturable 3-D
   HDiv-VIM in the convex B-input / `A⁺` form (the de Rham-dual Newton) stay
   well-conditioned at the knee where the M-Picard diverges? Validate against
   the scalar sphere fixed point `H_int = 3H₀/(μ_r(H_int)+2)` and then a
   `div M ≠ 0` body. (The existing `solve_nonlinear_newton` already gives the
   right answer via a damped Newton; the question is whether the convex form
   removes the damping/warmstart heuristics.)

2. **Linearisation (the real prize, open).** Is there a 3-D analogue of the
   Chaplygin linearisation — a coordinate system (local Clebsch / flux
   coordinates) in which the 3-D saturation becomes a *linear* variable-
   coefficient problem on a fixed domain — valid wherever the helicity
   obstruction vanishes (e.g. integrable / foliated fields, flux guides with a
   global flux function)? This is genuinely open.

## Status

- **Verified + committed:**
  - the linear bridge — Clebsch = loop modes
    (`hdiv_vim_clebsch_loopstar.py`, 3-D);
  - the 2-D unification — `A_z` *is* the Clebsch potential
    (`hdiv_vim_clebsch_2d_az.py`: a loop field `∇A_z×ẑ` is machine-zero
    divergence *and* boundary-tangential → fully charge-free; `A_z` recovered
    from `B` to `4e-7`);
  - the convex B-input conditioning cure and the 2-D hodograph linearisation
    (`saturation_loop_2d.py` / `chaplygin_hodograph_2d.py`);
  - the dynamical face — flux-line closure needs a closed 2-form *and* a
    symplectic integrator (`flux_line_closure_symplectic.py`; Noguchi / Sugahara
    2020), and the same diagnostic **on a real solved NGSolve field**: the de Rham
    `rot(grad A_z)` reconstruction's flux line closes while nodal-averaged / charge-
    admixed reconstructions spiral (`flux_line_realfield_ngsolve.py`); and the
    **order-vs-representation** sweep — closure is set by the de Rham representation,
    not the polynomial order, so a de Rham 2nd-order field closes while a nodally-
    smoothed one does not, even at 2nd order (`derham_closure_order_sweep.py`);
  - the 3-D closing condition — **helicity** is the obstruction to a global
    Clebsch pair (`clebsch_3d_closing_condition.py`: Clebsch `h = −7e-17`, ABC
    `h = 3(2π)³ ≠ 0` + chaotic Poincaré).
  (Each locked by a golden in `validation_test/feec/test_clebsch_hodograph_research.py`.)
- **Rung-1 finding (conditioning).** The convex `A⁺` form *is* well-conditioned
  at the knee (the linear field-exact solve gives the sphere `D = 1/3`; the
  committed `solve_nonlinear` outer-χ loop on `A⁺ = (1/χ)M_mass + N` converges).
  The remaining obstacle to a *practical* field-exact nonlinear solve is **not
  conditioning but cost**: the exact-field operator is assembled by an
  element-by-element charge field, which is `O(N²)`.  Production solve results
  now own a persistent C++ RT1 field evaluator: NumPy target buffers, one-pass
  IMA, TaskManager observation parallelism, exact analytic tet leaves, and a
  direct-probed quadrupole source tree for large non-IMA maps.  IMA automatic
  evaluation remains direct to preserve the reduced/full roundoff contract.
- **Open (the hard prize):** the 3-D saturation linearisation — a 3-D
  Chaplygin-analogue valid where the helicity obstruction vanishes.

Per the repository-first policy this note records the *map* — what is connected
and what is not — so the frontier is dug in the right place rather than
re-derived.

## References

The Clebsch-potential / flux-line / helicity structure used throughout this note
is classical; the examples cite:

- **P. Robert**, "Clebsch Potentials and the Visualization of Three-Dimensional
  Solenoidal Vector Fields," *IEEE Trans. Magn.* **27**(5), Sept. 1991. The
  foundational reference for this bridge, in the *same differential-forms
  language*: field intensity `H ∈ Ω¹`, flux density `B ∈ Ω²`, and the **global
  Clebsch existence as the de Rham question** "*given `β ∈ Ω²` with `dβ = 0`, do
  (possibly multivalued) `ξ, η` exist with `β = dξ ∧ dη`?*" — obstructed by the
  **helicity / Hopf invariant** ("tangled flux lines"). Also gives the
  micromagnetics gyrovector `g = −∇(cosθ)×∇φ` (Clebsch potentials of a
  magnetization map `R → S²`) and the **stream-function inverse** (current
  synthesis), i.e. the field ⇄ current de Rham duality.
- **M. Hirahatake, S. Noguchi, H. Igarashi, H. Yamashita**, "Analytical
  Calculation of Magnetic Flux Line from Hexahedral **Edge** Finite Element
  Analysis and its Allocation using the **Bubble System** for Visualization,"
  *IEEJ Trans.* (pp. 1205–1212). Computes flux lines from edge (`H(curl)`) FE
  results by solving `dξ/Bξ = dη/Bη = dζ/Bζ` on the normalised hexahedron — the
  practical edge-FE flux-line method whose closure is exactly the de Rham
  "closed 2-form" requirement (`B = curl A` is divergence-free by construction).
- **H. K. Moffatt**, "The degree of knottedness of tangled vortex lines,"
  *J. Fluid Mech.* **35**, 1969 — helicity = topological linking; the obstruction
  to a global Clebsch pair.
- Clebsch / Euler-potential foundations: *Reviews of Plasma Physics* (the
  flux-coordinate / Euler-potential literature); Bretherton (1970) on the
  Hamiltonian/Clebsch formulation; Arnold; Dombre et al., *JFM* **167**, 1986
  (ABC chaos).
