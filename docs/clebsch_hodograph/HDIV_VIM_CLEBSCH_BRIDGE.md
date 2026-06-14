# The HDiv-VIM ↔ Clebsch bridge (the de Rham capstone) and the saturation frontier

This note maps how the two research lines in `examples/clebsch_hodograph/` —
the **HDiv-VIM demag solver** (`radia.hdiv_vim`, the FEEC `H(div)` element Radia
is migrating to) and the **Clebsch hodograph** design line — are one structure
seen from two sides, and pins down exactly what is *verified* versus what is the
*open* hard part.

The verified linear (kinematic) level is
[`examples/clebsch_hodograph/hdiv_vim_clebsch_loopstar.py`](../../examples/clebsch_hodograph/hdiv_vim_clebsch_loopstar.py)
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

## The open frontier (3-D helicity)

In 3-D the Clebsch representation needs **two** potentials `(α,β)` for three
coordinates, and a *global* Clebsch pair need not exist — the obstruction is
**helicity** (`∫ A·B`). So the 3-D saturable HDiv-VIM solve does **not**
auto-linearise the way the 2-D hodograph does. The concrete open questions:

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

- **Verified + committed:** the linear bridge (Clebsch = loop modes), the 2-D
  hodograph linearisation, the convex B-input conditioning cure (each locked by
  a golden in `tests/feec/test_clebsch_hodograph_research.py`).
- **Open:** the 3-D conditioning rung (tractable) and the 3-D saturation
  linearisation (the hard prize, helicity-obstructed in general).

Per the repository-first policy this note records the *map* — what is connected
and what is not — so the frontier is dug in the right place rather than
re-derived.
