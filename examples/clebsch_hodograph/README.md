# Clebsch hodograph — research examples

Research demonstrations of the **Clebsch / hodograph (potential-coordinate)**
formulation of magnetostatics, extending the verified forward **"Clebsch
hodograph"** mode shipped in the `radia-electromagnet` panel
(`src/radia/panels/calc_clebsch_hodograph.py`).

> The design methodology these examples support (accelerator magnet design
> via the A–φ duality + reduced potential) is in
> [`docs/clebsch_hodograph/DESIGN_METHODOLOGY.md`](../../docs/clebsch_hodograph/DESIGN_METHODOLOGY.md).

These scripts are **research examples, not panel modes.** Per the repo-first
"Don't Publish the Unfinished" policy, only the verified *forward* mode is in
the panel; the *inverse* pole-design (3-D geometry-as-unknown hodograph PDE)
stays on the research branch until it is verified on a real device.

## Background

Clebsch / Euler potentials: `B = ∇φ × ∇ψ`, and a vector potential is
`A = φ ∇ψ`. The two scalar potentials are conjugate "potential coordinates";
the **hodograph** swaps them with the spatial coordinates (Hamilton 1846 —
the dependent↔independent variable swap, *not* a Legendre transform).

The panel's forward mode solves the **reduced scalar potential Ω** as primary
and recovers the flux `ψ` and scalar `Φ` as conjugates, reporting the
**bidirectional consistency** `‖B(from ψ) − B(from Φ)‖ / ‖B‖` — an
a-posteriori field-quality metric. These examples push two axes the panel
mode does not cover.

## Scripts

### `a_method_clebsch_2d.py` — A-method (vector potential primary)

The **dual** of the panel mode: the **vector potential `A_z` is the primary
unknown** (`div((1/μ)∇A_z) = 0`), and the conjugate scalar `V` is recovered.
`(A_z, μ₀V)` are Cauchy–Riemann conjugates — `W = A_z + iμ₀V` is analytic —
so the orthogonal nets `{A_z = const}` (flux lines) and `{V = const}`
(equipotentials) **are** the hodograph grid (the classical complex-potential /
conformal pole-design picture, here as an FEM realisation).

Verified (μ_r = 1000, order 3, maxh 0.08): interior `B = 2μ_r/(μ_r+1)·B₀`
to ~3e-5; `A_z`/`V` accuracy ~5e-5; hodograph consistency ~4e-4.
Figure: `a_method_clebsch_2d.png` (the orthogonal hodograph net).

### `hodograph_kelvin_axisym.py` — exact open boundary via Kelvin

Replaces the panel mode's *far-truncated* sphere boundary (`r/a = 6`,
~3e-3 field error) with the lab's verified **axisymmetric Kelvin
transformation** (open boundary, `μ_ext = μ₀(R/ρ′)²`, periodic
`kelvin_int ↔ kelvin_ext`), then computes the hodograph (flux `ψ`,
anchored `ψ = 0` on the axis; algebraic scalar `Φ = H₀z + Ω_r`) on the
exact field.

Verified (μ_r = 100, order 3, maxh 0.03): interior `Hz = 3/(μ_r+2)·H₀`,
**field_error ~1e-7** (vs ~3e-3 truncated — the Kelvin open boundary is the
win); hodograph consistency ~2e-5. Figure: `hodograph_kelvin_axisym.png`.

### `hodograph_kelvin_2d.py` — Kelvin in the hodograph (2-D Cartesian, no air box)

**rung 1 of "Kelvin in the hodograph."** The dipole FEM rung truncates the
exterior with a Dirichlet air box; this replaces it with the **2-D Cartesian
Kelvin transformation**, which is especially clean: in 2-D the Dirichlet energy
`∫|∇u|²` is **conformally invariant** and the Kelvin inversion `z → R²/z̄` is
conformal, so the Kelvin exterior is **weight-free** (`μ′ = μ₀`, *no* `(R/ρ′)²`
factor — unlike the axisymmetric case above). The open boundary at infinity maps
to a single interior point (`GND`). In a current-free 2-D region the flux
function `A_z` and scalar `V` are **Cauchy–Riemann conjugates** (`W = A_z + iμ₀V`
analytic), so "Kelvin in the hodograph" is literally **analytic continuation of
`W` across the conformal inversion** — the classical complex-potential pole-design
picture with an exact open boundary.

A magnetisable cylinder (`μ_r`) in a uniform field, reduced-Ω, the uniform
background carried into the Kelvin exterior by the **2-D** reduced-potential rule
(`-F_s`, `radia.kelvin_material.make_reduced_potential_background_cf(dim=2)`).
Verified (μ_r = 100, a = 0.3, Kelvin R = 1.0, order 3, maxh 0.04): interior
`B = 2μ_r/(μ_r+1)·B₀`, **field_error ~2e-8** (machine precision) vs **~2.8e-2**
for a truncated air box at `r/a = 6` — **~10⁶× more accurate**; `Bx ~ 1e-17`
(no transverse field by symmetry); hodograph consistency `B(A_z)` vs `B(V)`
~3e-5. Figure: `hodograph_kelvin_2d.png` (the orthogonal flux/equipotential net
on the exact field).

### `bidirectional_coordinate_transform_2d.py` — the coordinate transform, unified (Tampere-inspired)

**The exterior-calculus statement behind every map in this folder.** Inspired by
Dervisha, Marjamäki, **Rasilo** & Tarhasaari, *"Bidirectional Coordinate
Transformation and Its Application to 2-D Magnetic Field Problems"* (CEFC 2026,
WA-O1, Tampere): exterior calculus gives the governing equation in **any**
coordinate system, so a magnetostatic problem may be solved in whichever domain is
simplest — the physical one (carries the geometry) or a transformed one (a complex
map to a simple computational domain). For a map `F` with Jacobian `J`, the 2-D
magnetostatic Dirichlet energy pulls back with the **metric weight**
`W = |det J|·(JᵀJ)⁻¹`:

```
∫_phys ∇u·∇v dx dy  =  ∫_comp (∇u)ᵀ W (∇v) dξ dη .
```

So assembling on the **deformed physical mesh** and assembling on the **simple
computational mesh with `W`** give the **identical** stiffness matrix — verified
here to `~1e-16` for a *conformal* map (`F = z + 0.3z²`) **and** a *non-conformal*
stretch+shear, with the two solutions and the exact harmonic field `Re(z³)` all
agreeing. This is the **bidirectional** statement: solve where it is convenient,
the metric carries the geometry; and it works for **arbitrary** maps, not only
conformal ones.

**It unifies the lab's coordinate transforms.** A 2-D **conformal** map has
`J = sR`, so `W = s²·(s²I)⁻¹ = I` — **weight-free**. That is *exactly* why
`hodograph_kelvin_2d.py` carries `μ′ = μ₀` with no `(R/ρ′)²` factor (the Kelvin
inversion is conformal) and why the hodograph/Cauchy–Riemann nets need no metric:
the 2-D Kelvin and hodograph maps are the *same* object — a weight-free pullback.
The example also **derives this symbolically** (`symbolic_pullback_check`, sympy):
2-D conformal `W = I`, a generic 2-D map `W ≠ I`, and **3-D / axisymmetric keep a
weight** (`s` for a 3-D scaling, `2πr` for the meridian plane — the reason
`clebsch_kelvin_3d.py` and `hodograph_kelvin_axisym.py` are *not* weight-free).
Figure: `bidirectional_coordinate_transform_2d.png` (the physical coordinate grids:
(a) conformal — angles preserved, `W=I` | (b) non-conformal — skewed, `W≠I`).
Golden `test_bidirectional_coordinate_transform_2d` (fast).

### `weakform_pullback_kata.py` — the weak form in differential forms (the kata)

**The "kata" companion to the bidirectional file.** Works the *same* exterior-calculus
pullback for the *specific* magnetostatic conjugate map — the hodograph. Writes the
2-D Laplace weak form in forms (`d(μ⋆du)=0` → `a(u,v)=∫ μ dv∧⋆du = ∫ μ ∇u·∇v`, one line
because `α∧⋆β = ⟨α,β⟩vol`), then pulls it into the `(A,φ)` potential chart and shows the
material survives as the **anisotropic** weight `W = diag(μ, 1/μ)` — exactly the Tampere
bidirectional governing equation (8) (the unknown having swapped to the geometry
`x(A,φ), y(A,φ)`). **Derived symbolically** (sympy): the line element
`g = diag(1/(μq)², 1/q²)` gives `W = diag(μ, 1/μ)` with the field magnitude `q`
**cancelling**, and the conjugate-map pullback for `F(z) = z, z², z³, 1/z` gives
`W = diag(μ, 1/μ)` *exactly* in every case (the conformal factor cancels). That
`q`-independence is why the hodograph linearises (for `μ(q)`, `W = diag(μ(q), 1/μ(q))`
is a known function of the coordinate — Chaplygin). Golden `test_weakform_pullback_kata`
(fast, sympy-only).

### `clebsch_kelvin_3d.py` — Kelvin in the hodograph (3-D, the real formulation)

**rung 2 of "Kelvin in the hodograph."** The genuine 3-D case: the **Clebsch
potentials** (`B = ∇ψ × ∇χ`) on the **3-D Cartesian Kelvin two-sphere domain**
(Sugahara 2022; `radia.kelvin_geometry.add_kelvin_exterior_domain`), the air box
abolished. The 3-D Kelvin pullback (de Rham / Nagamine CEFC 2026): the inversion
`T: x′ = R²x/|x|²` maps the physical exterior to the offset Kelvin sphere; the
Clebsch potentials `ψ, χ` are **0-forms** so they pull back trivially
(`ψ′ = ψ∘T`); the field `B` is a **2-form** so it carries `−(R/ρ′)⁴`; and — unlike
2-D — the 3-D Jacobian is **not** conformally invariant, so the reduced-Ω material
gets `μ′ = (R/ρ′)²·μ₀` (the standard 3-D Kelvin weight, *not* weight-free).

Verified on the canonical 3-D test (a magnetisable sphere in a uniform field,
exact interior `H = 3/(μ_r+2)·H₀`): reduced-Ω on the two-sphere domain (periodic
`kelvin_int ↔ kelvin_ext`, GND at the Kelvin centre = ∞, background via
`make_reduced_potential_background_cf(dim=3)`); the Clebsch net `χ = atan2(y,x)`
(the exact azimuthal coordinate) + `ψ` recovered as the Stokes flux, checked
`B = ∇ψ × ∇χ`. Result (μ_r=100, a=0.2, R_K=0.5, order 3, maxh 0.05): interior
**field_error ~1.5e-5** vs **~8e-3** for a truncated air box at `r/a=5`
(~10³× more accurate); `Hx ~7e-9`; Clebsch consistency `B(ψ,χ)` vs `B` ~6e-4
(off-axis, away from the `χ` branch cut). **The STRONG test** (the interior is
boundary-*insensitive*, so it does not prove the open boundary is right): the
**exterior** field matches the **exact uniform + induced-dipole** to **~1.4e-3**
(mesh-limited), confirming the Kelvin transform is interpreted correctly — the
material weight `μ′=(R/ρ′)²μ₀` is the exact conformal change-of-variables, and
the reduced-potential background (an *engineering* convention, not the
naive covariant pullback which is singular at ∞) reproduces the right exterior.
Figure: `clebsch_kelvin_3d.png` (the
Clebsch flux `ψ` field lines on a meridional slice). *Later rungs (both below
and DONE): 1.5 = 2-D Chaplygin (the hodograph **linearises** the saturation
nonlinearity — Molenbroek–Chaplygin); 3 = the 3-D merged geometry+material
single Picard.*

### `saturation_loop_2d.py` — the nonlinear saturation loop (Chaplygin reference)

The **nonlinear half** of the program (the user's "非線形のループ"), and the
reference the **Chaplygin** rung will reproduce *without* iteration. A
magnetisable cylinder whose permeability **saturates**, solved by a nonlinear
FEM loop. Formulation = the standard **stable** one: the **A-formulation with a
B-input reluctivity** `ν(|B|)` (`B = (∂A_z/∂y, −∂A_z/∂x)`, `|B| = |∇A_z|`),
Fröhlich `μ_r(B) = 1 + (μ_r0−1)/(1+(B/B_k)²)`. Plain Picard + a `B₀`
continuation. Result (μ_r0=10, B_k=1 T): the interior ratio `B_in/B₀` falls from
the **unsaturated demag value** `2μ_r0/(μ_r0+1) ≈ 1.82` toward **1** (saturated
cylinder, "transparent"), monotone, respecting the demag limit `[1, 1.82]`.

**Why this formulation** (the lesson worth keeping): the *reduced-Ω* `μ(|H|)`
loop is **ill-conditioned** for a saturable cylinder (`μ_r` steepest at small
H, exactly at the cylinder boundary) — its Picard converges to a **spurious,
non-physical** fixed point (interior B *exceeding* the demag limit `2B₀`,
impossible). The B-input `ν(|B|)` A-formulation has a **convex** energy, so
plain Picard converges cleanly **to machine precision**: the converged iterate
is its own frozen re-solve (`|blend − frozen| ~ 1e-12`) — the diagnostic that
the loop found the *true* solution, not a false fixed point. Open boundary: a
large domain (`R/a = 10`, < 1% truncation) stands in for the exact one
(saturation is a *local* effect, orthogonal to the rung 1-2 open-boundary
treatment; the two compose). Figure: `saturation_loop_2d.png` (the saturation
S-curve). *Next rung — `chaplygin_hodograph_2d.py` below — reproduces the
saturation curve in one shot, no loop.*

### `chaplygin_hodograph_2d.py` — the hodograph linearises saturation (1-shot vs loop)

**The Chaplygin rung** (the user's *"hodograph 変換をすると どーせ反復が
必要なので 非線形ループと混ぜて短くできないか?"*). 2-D current-free
nonlinear magnetostatics **is** 2-D steady irrotational compressible gas flow
(`H↔v`, `B↔ρv`, `μ(|H|)↔ρ(q)`, scalar potential `Φ↔φ`, flux function `A↔ψ`).
Taking `q=|H|, θ=arg(H)` as the **independent** coordinates (the hodograph
plane; Molenbroek 1890, Chaplygin 1902) turns the nonlinear PDE into a
**linear, variable-coefficient** one: `μ(q)` is now a **coefficient** (a known
function of the independent coordinate `q`), **not** a nonlinearity in an
unknown field. The whole nonlinearity is absorbed — *one* linear solve, no
outer Picard over the constitutive law.

**Honest caveat (a correction to the naive "Chaplygin = transonic" framing):**
the operator `div(μ(|H|)∇Φ)` is **elliptic everywhere** for ordinary
saturation — `μ` falls but `μq = |B|` still rises, so `d(μq)/dq > 0` always.
There is **no limiting line, no shock, no hyperbolic hard part** (that
type-change needs falling `|B|` vs `|H|`, unphysical for a passive magnet). So
the linear Chaplygin equation is an elliptic boundary-value problem on the
**hodograph image** of the region — a clean 1-shot **iff that image is a
simple, fixed domain** (special geometries); general geometry gives a
free-boundary (unknown image), the genuinely hard case, deferred.

**The special geometry** (simple image): a slender saturable **flux guide**
pinched to a throat. The field is nearly axial, so the hodograph image is a
thin **band** about `θ=0` (a "simple wave" → segment), on which the Chaplygin
system collapses to a single quadrature:
`ΔΦ(Ψ) = ∫₀ᴸ ν(|B|(x)) |B|(x) dx` with `|B|(x) = Ψ/w(x)` (flux conservation,
geometry-only) — evaluated **once, no iteration**. Compared to the full 2-D
nonlinear FEM loop `div(ν(|∇A|)∇A)=0` on the actual curved-wall guide
(μ_r0=200, B_k=1 T): the 1-shot matches the loop to **1.0–1.5 %** across
`|B|_throat` = 0.25 → 2 T, and the saturation **bend** of drive/flux agrees
(**×2.78** loop vs **×2.77** 1-shot). The image is a band (max **21°** off-axis
at the throat shoulders, *not* a tautological 0 nor a 2-D blob); the 1-shot
holds because the drive integral is throat-**centre** dominated, where the
field is axial. And it is the **slender limit**: rel.err **8.7e-3 → 2.1e-3** as
the throat is made more gradual (slenderness 3 → 9). Figure:
`chaplygin_hodograph_2d.png` (1-shot vs loop + the rel.err). *Convergence
needs **under-relaxation** — undamped Picard oscillates under the throat's
strong reluctivity contrast even though the problem is convex.* **Next rung: a
guide that *turns* (`θ` varies) → a real 2-D hodograph image needing a linear
PDE solve; rung 3 = the 3-D merged geometry+material single Picard (next section).*

### `chaplygin_design_sweep_2d.py` — nonlinear-as-linear, MAXIMALLY: a saturable design space at linear cost

The payoff of "the hodograph makes saturation a **coefficient**": a whole
saturable-magnet **design space** is computed by **one linear quadrature per
point**. A grid of (throat width × throat field) — 36 nonlinear operating
points — is the **design map** `ΔΦ(w_throat, |B|_throat)`, each point the
mesh-free 1-shot drive, in **~2 ms**; the drive **bends up ×3.4–4.0** as the
throat saturates (the design content). Inverting it gives the **flux-regulator /
field-limiter** transfer `Ψ(drive)`, which **clamps** (concave) at saturation,
and a **narrower throat clamps at a lower flux** (`knee Ψ ≈ B_k·w_throat`) — the
throat width is the clamp-level knob. Validated against the full **nonlinear FEM
Picard** at sample points: the 1-shot agrees to **~1–2 %** (slenderness-limited),
and each FEM point costs **~12 Picard iterations**. So the whole map = **36
quadratures** vs the equivalent **~440 curved-mesh linear solves** the nonlinear
FEM would need — the **entire nonlinear design space explored at linear cost**,
where the Picard sweep would be prohibitive. **Inverse design (spec → geometry)**:
`size_for_clamp(target)` sizes the throat so the regulator **clamps at a target
flux** by bisecting the throat depth (the clamp knee is monotone in depth) — a
handful of mesh-free 1-shot quadratures, **sub-second**, hitting the target to
machine precision (e.g. **12 mWb/m → w = 14.1 mm**); the **sized design is then
confirmed by the full nonlinear FEM** (1-shot 477 A vs FEM 483 A, **1.2 %**) — the
**design loop closes**. This is the engineer's actual step (spec → geometry), not
just a forward sweep. Figure: `chaplygin_design_sweep_2d.png` (design map | the
clamp curves with the sized design ★ | the cost bar). Goldens
`test_chaplygin_design_sweep_*` + `test_chaplygin_inverse_design_sizing` (sizing
fast, FEM-validate+cost+loop-close slow).
*Honest scope: the 1-shot is the slender-guide limit (its error is the
slenderness error); a guide that **turns** needs the linear hodograph PDE
(`chaplygin_turning_guide_2d.py`) — still one **linear** solve per case, so the
same "design space at linear cost" exploit applies one rung up.*

### `clebsch_kelvin_nonlinear_3d.py` — rung 3: the 3-D merge (saturation + exact open boundary, one Picard)

**The capstone.** Rungs 1–2 put the Kelvin open boundary in the
hodograph/Clebsch frame (linear); rung 1.5b showed the 2-D hodograph
*linearises* saturation (1-shot, no loop). **3-D does NOT auto-linearise** —
the Clebsch pair (ψ,χ) is two potentials for three coordinates (a gauge/helicity
freedom), so there is no clean (x,y,z)↔hodograph interchange. So the honest 3-D
form of *"混ぜて短く"* is a **single Picard that updates the geometry (the exact
Kelvin open boundary, rung 2) and the material (μ(|H|) saturation, rung 1.5a)
TOGETHER** — instead of nesting an air-box-truncation loop inside a material
loop. The geometry (open boundary) lives in the *fixed* Kelvin weight
`μ'=(R/ρ')²μ₀`; only the magnet's μ changes per step. There is no air-box loop
to converge — it is exact, once.

A saturable **sphere** in a uniform field has a *uniform* interior, so its
nonlinear solution is exactly a linear sphere with the self-consistent
`μ_r_eff` — the exact reference is the scalar demag fixed point
`H_int = 3H₀/(μ_r(H_int)+2)`. The single reduced-Ω Picard on the Kelvin
two-sphere (μ_r0=20, H_k=0.25, order 3, maxh 0.06, ne 14301) converges in
**42 iters** to:
- interior `field_error` **2.5e-4** vs the exact fixed point (saturated to
  `μ_r_eff`≈**11.9** < μ_r0=20);
- **EXTERIOR** `field_error` **3.1e-3** vs the equivalent induced dipole (the
  strong Kelvin test, now on the *nonlinear* solution) — vs **6.7e-3** for a
  truncated r/a=5 air box (the Kelvin win, ~26×);
- self-consistency **3.2e-9** (the iterate IS its own frozen re-solve = a true
  fixed point; the residual 2.5e-4 is *mesh* error, not iteration error).

**Honest scope:** the STABLE saturation regime. The H-input reduced-Ω fixed
point `H=3H₀/(μ_r(H)+2)` is contractive only while `|dRHS/dH|<1` (here ρ≈0.67);
for a *very steep* B-H curve that exceeds 1 and the H-input loop goes unstable —
and **no under-relaxation cures a positive-slope unstable fixed point** (the
general cure is the convex B-input ν(|B|) form, rung 1.5a's lesson, in 3-D a
vector-potential formulation). The figure reports ρ so the limit is explicit.
Figure: `clebsch_kelvin_nonlinear_3d.png` (Picard convergence to the exact
fixed point + the self-consistent BH operating point). *This closes the
program: 1–2 (exact open boundary), 1.5a (nonlinear loop), 1.5b (2-D 1-shot),
3 (3-D merge).*

### `chaplygin_turning_guide_2d.py` — the Chaplygin frontier: a TURNING field, ONE linear PDE (not a quadrature)

One step past rung 1.5b. There the slender guide's hodograph image was a thin
**segment** (θ≈const) → the Chaplygin system collapsed to a 1-shot **quadrature**.
Here the field **turns**, so θ varies over a genuine 2-D range and the hodograph
image is a 2-D **region** → the Chaplygin equation is a genuine **linear,
variable-coefficient elliptic PDE that must be SOLVED** (one linear solve), not
integrated. For the flux function `A(q,θ)` it is the self-adjoint form

```
d/dq( (q/μ) A_q ) + ((μq)′/(μ²q)) A_θθ = 0,
```

with `μ(q)` a **known coefficient** (the linearised nonlinearity). The μ_r=const
limit is `A_qq+(1/q)A_q+(1/q²)A_θθ=0` = **Laplace in (ln q, θ)**, so `A=ln(q)·θ`
is an exact solution — the solver verification: the FEM reproduces it to
**6.3e-8**. The saturating case (μ_r0=20, q spanning q_k) is the **same single
linear solve** (residual ~2e-15) with the saturating coefficient; the solution
bends ~34% away from the linear harmonic (genuine 2-D content). The result
**back-maps** to physical space (integrate `dx,dy`, recover `Φ` from the
first-order pair `Φ_q=-q(μq)′/(μq)²·A_θ`, `A_q=(μ/q)Φ_θ`) **single-valued to
1.1e-2** — a realisable, genuine 2-D nonlinear **turning** field from one linear
solve. Figure: `chaplygin_turning_guide_2d.png` (the 2-D hodograph image + the
back-mapped physical turning region).

**Honest scope / remaining frontier:** a turning flux *guide* (iron walls = flux
lines) on a **fixed** hodograph rectangle is exactly the **constant-width** bend,
which is θ-independent (1-D, self-linearising — `|H|~1/r` is forced by geometry,
μ only reshapes B). A genuinely turning+tapering guide has θ-**dependent** wall
positions → a **θ-dependent hodograph image = a free boundary**. So this file
does the *forward* construction (prescribe 2-D hodograph data → one linear solve
→ back-map to *a* physical patch); the **free-boundary closure** (prescribe the
physical guide, solve for its hodograph image) is the genuine open frontier.

### `chaplygin_turning_design_sweep_2d.py` — design-at-linear-cost, one rung up: TURNING fields, one linear solve each

`chaplygin_design_sweep_2d.py` exploited the **slender** guide (segment image →
1-shot quadrature). This is the rung up: a field that **turns** has a 2-D
hodograph image needing a genuine **linear elliptic PDE** — but still **linear**
(`μ(q)` a coefficient), so each nonlinear turning field is **one direct linear
solve, no Picard**. So a whole design space of TURNING saturable fields is at
linear cost too. Sweeping **material** `μ_r0` (5→20→80) and **operating depth**
`q₁/q_k` (1.2→2.0→3.0), each entry is one solve with **machine-zero residual**
(`~8e-15` — it *is* one direct solve); the **saturation bend** (deviation from the
linear `ln(q)·θ` harmonic) grows monotonically (`0.12→0.47→2.6` with `μ_r0`); and
each **back-maps single-valued** to a realisable physical turning field (closure
`≤1.4e-2`). The `μ_r=1` Laplace limit reproduces the exact harmonic to **3e-8**
(the solver is correct). **The linear-cost win:** 7 turning designs in **~450 ms**
(7 direct linear solves) vs the equivalent physical-space nonlinear solve at
**~12 Picard iterations each** (measured on the same Froehlich material via the
slender-guide reference) ≈ **88 linear solves** — and the turning 2-D case is at
least as expensive. Figure: `chaplygin_turning_design_sweep_2d.png` (bend vs
design | a back-mapped physical turning field | the cost bar). Goldens
`test_chaplygin_turning_design_sweep_*` (sweep fast, Picard-cost slow). *Honest
scope: the forward construction (free-boundary closure deferred, as above); each
of those Newton steps is also **linear**, so the exploit carries.*

### `chaplygin_free_boundary_2d.py` — the turning-guide free boundary, made concrete

The achievable, honest piece of that open frontier. The turning-guide forward
construction solved the Chaplygin PDE on a *fixed* hodograph rectangle; the
**inverse** direction — prescribe the physical guide, find its hodograph image —
is a **free-boundary** problem. This file *computes* the hodograph image of a
prescribed turning flux guide (a first-quadrant annular-sector bend) and shows
concretely:
- **constant-width circular bend → the image is a RECTANGLE** in (q,θ): the
  field turns (θ_B spans the bend, here ~56°) but the q-extent is the *same* at
  every position angle (free-boundary measure **0.05** ≈ 0). This is the **1-D
  self-linearising** case (`|H|∝1/r` forced by geometry).
- **tapering bend (outer wall spirals inward) → a θ-DEPENDENT image** (free-
  boundary measure **0.18**, ~3.6×): the gap width varies along the bend, so the
  q-extent varies with θ — a genuine **free boundary** (the image is unknown a
  priori).

The image's boundary shape is a *geometric* property of the field, so a linear
flux solve already exhibits it (saturation shifts the q-values, not the
rectangle-vs-free-boundary distinction). Figure: `chaplygin_free_boundary_2d.png`
(the two hodograph images, coloured by position angle). *Recovering the image
from the prescribed physical boundary — the inverse hodograph solve / free-
boundary iteration — remains the open frontier (a Newton/mixed free-boundary
scheme); this file establishes its structure.*

### `chaplygin_inverse_vonmises_2d.py` — dissolving the free boundary (von Mises inverse)

The inverse direction, attacked. The free boundary is an artefact of the
**coordinate choice**: take the potential `Φ` and the flux function `A` as the
**independent** variables (von Mises). A flux guide is bounded by two flux lines
(`A=0,Ψ` = walls) and two equipotentials (`Φ=const` = ports), so in `(Φ,A)` the
domain is **always the fixed rectangle** `[0,Φ1]×[0,A1]` — the free boundary is
**gone**. One solves instead for the physical map `(x,y)(Φ,A)` from the same
first-order Chaplygin pair, here `x_A=-y_Φ/μ(q)`, `y_A=x_Φ/μ(q)`,
`q=1/|∇_Φ(x,y)|`, in **least squares**. For `μ=1` this is Cauchy–Riemann
(`x+iy` analytic in `Φ+iA`), so the exact annular-bend map is the conformal
`f=e^{i(Φ+iA)}` — and the solver **recovers it to 3.6e-9** (residual `J→0`): the
free boundary is dissolved into a fixed-domain solve. Figure:
`chaplygin_inverse_vonmises_2d.png` (the fixed `(Φ,A)` rectangle + the recovered
physical map).

**The nonlinear wall — now CLOSED in `chaplygin_inverse_nonlinear_2d.py`.** With
`μ=μ(q)` the A-spacing `1/(μq)` changes, so the consistent `Φ,A` distribution
*along the walls* is μ-dependent, and the correct BC is a **slip condition** (the
boundary point lies *on* the wall curve, tangential position free — a nonlinear
constraint for curved walls). The remaining obstacle was an **over-determination**:
the von Mises rectangle has two dimensions, `Φ1` (the MMF) and `A1` (the flux `Ψ`),
and a guide of fixed geometry+material cannot have **both** prescribed (its
permeance fixes one once the other is chosen). This file solves the linear case;
the companion closes the nonlinear case.

### `chaplygin_inverse_nonlinear_2d.py` — the nonlinear free-boundary inverse, CLOSED

The fix for the over-determination: **prescribe `Φ1=ψmax` (the drive/MMF) and
free `A1=λ`** — the μ-dependent saturable flux — as a single **global
`NumberSpace` unknown**. Polar physical unknowns `(r,ψ)` make the arc walls
coordinate lines, so the slip BCs are clean (`r=r_in` inner Dirichlet, `ψ=0,ψmax`
ports, the tapered outer wall via an on-curve penalty `β(r−r_out(ψ))²`). Damped
Newton with continuation in `Ms`. Verified:

| case | `J` (PDE residual) | wall fit | `free_measure` | map |
|---|---|---|---|---|
| const-width | **2.6e-18** (machine zero) | — | 0.04 (rectangle) | valid (det>0) |
| tapered 30% | **1.4e-07** | 2.3e-09 | 1.18 (θ-dependent) | valid (det>0) |
| extreme 50% | 1e-3 | — | — | **folds** (throat limit) |

The saturable flux grows far above its linear value (`λ: 0.69→9.4`, ~14× — a
high-μ guide carries much more flux for the same MMF). Three independent checks:
the PDE residual `J→0`, the on-curve **wall fit ~1e-9** (the map's outer edge
*does* lie on the prescribed tapered wall), and `free_measure≈1.18` vs `0.04`
(the hodograph image is genuinely **θ-dependent = the free boundary recovered**,
not collapsed to a rectangle). The previous `J~0.24` plateau is gone (6 orders
better). Only the extreme 50% taper folds — the throat narrows until the
saturable guide can no longer carry monotone flux, a **geometric** limit
reported honestly, not a formulation failure. Figure:
`chaplygin_inverse_nonlinear_2d.png` (const-width annular bend + tapered spiral,
both with `det>0`).

### `chaplygin_taper_design_sweep_2d.py` — closing Task 1's free boundary: a TAPER design sweep, Newton steps all linear

The connection of the **turning design sweep** (`chaplygin_turning_design_sweep_2d.py`,
forward, one linear solve each) to the **free-boundary closure** above. A guide that
turns *and* **tapers** has θ-dependent walls → a θ-dependent hodograph image = a **free
boundary**; the **von Mises** change of variables dissolves it onto a **fixed
rectangle**, so each tapered design is the nonlinear inverse — a **Newton continuation
whose every step is one linear solve, on a fixed mesh (no remesh)**. Sweeping the
**taper** (the free-boundary lever, `solve_inverse(taper)`): each closes (`J ≤ 3e-6`)
with a **valid, non-folded map** (`jac_min ≥ 0.05`), and the **free-boundary measure
grows** `0.04 → 0.86 → 1.19` as the taper goes `0 → 0.15 → 0.30` — constant width is a
rectangle image (self-linearising), tapering is θ-dependent = the free boundary
recovered (the saturable flux `λ` falls `5.9 → 3.6` as the throat tightens). So even
the **turning+tapering (free-boundary) saturable design space** is reachable, the
nonlinearity handled by Newton with **linear** steps — 3 designs in ~4 s, fixed mesh,
where a physical-space free-boundary solver would re-mesh and re-solve a nonlinear
problem at every shape iteration. Figure: `chaplygin_taper_design_sweep_2d.png` (the
free-boundary measure + λ vs taper | the constant-width annular bend | the tapered
free-boundary guide). Golden `test_chaplygin_taper_design_sweep`. *Honest scope: the
slip-wall von Mises formulation; extreme taper eventually **folds** at the geometric
throat limit (`chaplygin_inverse_nonlinear_2d.py`) — this sweep stays in the valid
regime.*

### `cohomology_hodograph_currentlink.py` — when the hodograph needs cohomology

Answers *"is cohomology needed for the hodograph?"*. The hodograph's **scalar
coordinate** `V` (`H = ∇V`) is single-valued **iff the field has zero period
(zero linked current)** around every hole — *not* merely iff the region is
simply connected:

| region | period `∮H·dl` | scalar coordinate `V` |
|---|---|---|
| simply connected (`b₁=0`) | — | single-valued (no cohomology) |
| multiply connected, **zero** period | `0` | single-valued (no cut) |
| multiply connected, **non-zero** period | `= I ≠ 0` | **multi-valued → cohomology** |

So the cylinder/sphere demos (no current) need no cohomology; a real magnet
with a **coil window** (current threads the air) does. Demonstrated on a
washer (`b₁=1`) with the lab's **gmsh-free `radia.cohomology`**: the
unit-circulation harmonic generator `h` is curl-free but **not a gradient**
(`∮h·dl = 1`), so the best single-valued grad-fit residual `‖F−∇V‖/‖F‖` is
**~0.95 for the current-linking field** (no single-valued `V` → cohomology
required) vs **~1e-11 for a zero-period field** (single-valued `V` exists).
The solid cylinder (`b₁=0`) has no generator at all.

In de Rham terms: cohomology = the **periods** that single-valued potentials
miss. Which conjugate coordinate carries the obstruction depends on the
physics — the **scalar** `V/Φ` for *current-linking* (`∮H·dl≠0`, the T-Ω cut),
the **flux** `ψ/A` for *net flux* through a hole. In genuine 3-D there is an
additional **helicity** obstruction (a global Clebsch pair may not exist) —
that is the research frontier, not covered here.

> **Architectural note.** `radia.cohomology` is the *interim* home. NGSolve /
> Netgen ship the de Rham complex (the natural fit) but **no turnkey homology
> generator**; gmsh has one (Pellikka 2013). Native NGSolve support
> (a harmonic-1-form extractor in `ngsolve.comp`) would be the proper place —
> `radia.cohomology` is structured to be upstreamable.

### `hdiv_vim_clebsch_loopstar.py` — the de Rham capstone: HDiv-VIM loop modes ARE Clebsch fields

The bridge between **this directory's Clebsch line** and the **HDiv-VIM demag
solver** (`radia.vim`, the FEEC `H(div)` element that Radia is migrating
to). The Hodge / Helmholtz split of a magnetization is

```
M  =  ∇φ            (+)     ∇α × ∇β
      "star"                "loop" = Clebsch
      carries the charge    charge-free → FIELD-NULL in N = BᵀGB
      pole-forming          flux-guiding (yoke return)
```

The HDiv-VIM is built so its demag operator `N = BᵀGB` (B = the magnetic-charge
map `M ↦ (ρ=−div M, σ=M·n)`) annihilates the divergence-free RT modes — *"loop
modes are field-null by construction"* (de Rham: a curl carries no charge). That
property **is** the statement that a Clebsch magnetization `∇α×∇β = d(α dβ)` (an
exact, hence closed, 2-form) makes **no demagnetizing field**. So the Clebsch
potentials `(α,β)` are the coordinates of the HDiv-VIM kernel: the loop–star
split *is* the Clebsch–gradient (Hodge) split *is* the pole / flux-guide split.

Verified on a unit sphere (HDiv order 1, exact analytic charge Gram):

| magnetization | form | demag factor `D` | external field |
|---|---|---|---|
| **Clebsch** `M=(y,−x,0)` | `∇(½r⊥²)×∇z` | **1.6e-3 ≈ 0** (‖div M‖=1e-14) | **2.4e-4** (≈0) |
| gradient `M=(x,y,z)` | `∇(½‖r‖²)` | 0.997 | — |
| uniform `M=(0,0,1)` | `∇z` | 0.333 ✓ (=1/3) | 9.2e-2 (dipole) |

`D_Clebsch / D_gradient = 1.6e-3`; and adding `t·Clebsch` to the uniform star
changes the external field by **<1.3 % even at t=5** (the residual is the
faceting `M·n`, → 0 under `mesh.Curve`). Three checks: (1) field-null, (2) no
stray field, (3) gauge — flux-circulation is invisible from outside.

This is the **linear (kinematic)** level. The nonlinear payoff — the saturable
HDiv-VIM solve reformulated in Clebsch/hodograph coordinates, where the
Chaplygin hodograph (`saturation_loop_2d.py` / `chaplygin_hodograph_2d.py`)
linearises saturation — is the next frontier.

### `hdiv_vim_clebsch_2d_az.py` — the 2-D unification: A_z IS the Clebsch potential

The 2-D companion of the capstone, and the explicit link to the Chaplygin
hodograph. In 2-D every in-plane divergence-free field is
`B = ∇A_z × ẑ = (∂A_z/∂y, −∂A_z/∂x)`, so the **flux function `A_z` is exactly the
Clebsch potential `α`** (with `β = z`). The 3-D Clebsch *pair* `(α,β)` collapses
to the **single scalar `A_z`** — which is precisely why the Chaplygin hodograph
linearises the 2-D saturation: `div(ν(|∇A_z|)∇A_z)=0` has one scalar unknown,
the Clebsch potential.

Verified (unit square, `H1` order 3): a loop field `B = ∇A_z × ẑ` built from a
known `A_z = sin πx sin πy` is **machine-zero divergence** (`‖div B‖ = 0`) *and*
**tangential on the boundary** (`‖B·n‖ = 4e-16`) — fully charge-free →
field-null (the 2-D loop–star, matching the 3-D capstone); a gradient field
`(x,y)` carries the charge (`‖div‖ = 2`). And `A_z` is **recovered from `B`** via
the stream-function weak form `∫∇A_z·∇w = ∫B·rot(w)` to `4e-7`, with
`∇A_z(rec)×ẑ` reproducing `B` to `2e-5` — `A_z` *is* the Clebsch potential.

### `flux_line_closure_symplectic.py` — the dynamical face: flux-line closure

A magnetic flux line is an integral curve of `B`: `dx/ds = B`. In 2-D this is
**Hamilton's equations with `A_z` as the Hamiltonian** (`ẋ = ∂A_z/∂y`,
`ẏ = −∂A_z/∂x`), so `A_z` (= the Clebsch potential) is conserved along a flux
line → flux lines are `A_z` level sets and **close**. Closure has two
requirements, isolated here:

- **the field must be a closed 2-form** (`div B = 0`): on a closed field the
  line closes (`A_z` drift `1e-2`, returns to start to `2e-4`); adding a
  charge/star part `ε∇φ` (div ≠ 0) gives no global `A_z` and the line **spirals**
  (drift `1.7e3`, `1.6e5×` worse) — *the de Rham / edge-FE requirement: flux
  lines computed from an edge (`H(curl)`) potential `B = curl A` are exactly
  divergence-free.*
- **the integrator must be symplectic** (`A_z`-conserving): same 1st-order
  forward Euler **spirals out** (drift `26`, `2440×` worse than symplectic) — *the
  accelerator-tracking requirement.*

`A_z` = the Clebsch potential = the flux-line-flow Hamiltonian is the single
object behind both the **field** (FEEC/Clebsch) face and the **dynamical**
(symplectic) face. *Refs: Noguchi, "Flux-line computation from hexahedral
edge-finite-element results / bubble-system placement", IEEJ (JP); Sugahara
2020, "Implicit symplectic flux-line tracking" (noting the circular-accelerator
beam-orbit-tracking analogy).*

### `flux_line_realfield_ngsolve.py` — the diagnostic on a REAL solved FE field

The *field* face of the same idea carried onto an actual NGSolve solve (the
script above isolated the *integrator* face on an analytic field). A 2-D
magnetostatic solve `−div(grad A_z) = J` for a current dipole gives `A_z ∈ H1`;
one flux line is traced for **three reconstructions of the same solve with the
same RK4 integrator**, so the only variable is the reconstruction:

- **de Rham `rot(grad A_z)` = `(∂A_z/∂y, −∂A_z/∂x)`** — the edge-FE `B = curl A`
  (Noguchi): a **closed 2-form** exactly tangent to the flux surfaces
  (`B·∇A_z = 0` pointwise, misalignment `0.0`) → the flux line **closes** (`A_z`
  drift `7e-4`, returns to `3e-4`, integrator floor).
- **nodal-averaged (`VectorH1` order 1)** — the realistic, parameter-free
  edge-vs-nodal leak: forcing both components continuous breaks the de Rham
  structure (misalignment `3.7e-2`) → **spirals** (`A_z` drift `21×` worse).
- **explicit charge leak `B + ε∇A_z` (`ε = 0.06`)** — the controlled
  de Rham-complement admixture (misalignment `6e-2`) → **spirals** dramatically
  (`A_z` drift `1.34`, `~2000×` worse, never returns).

This turns the de Rham picture into the **field-reconstruction-quality
diagnostic for the HDiv-VIM migration**: trace a flux line; if it spirals, the
reconstruction has leaked solenoidal content — the visible symptom of the
`M_mass⁻¹ N m` leak. The 3-D generalisation is helicity-obstructed (next file).
*Same refs as the symplectic example (Noguchi edge-FE flux lines).*

### `derham_closure_order_sweep.py` — order vs representation: which one closes the flux line?

Two old questions answered with an order sweep. **(Q1, extending Noguchi)** can the
edge-FE flux line (de Rham, exactly `div B = 0`) be *extended* — "if the field is
de Rham, a symplectic tracker can be built on it"? **(Q2, Kameari's remark)** "2nd-order
elements — the flux lines don't close" — is that because they are *not de Rham*?
Would a de Rham 2nd-order element close?

The measured answer (this script reproduces the table):

- **`B = curl A` is divergence-free for ANY conforming `A`** — edge `H(curl)` *or*
  nodal Lagrange `[H1]³` — at **every order** (weak interior divergence `~1e-15` at
  `p = 1, 2, 3`). So "edge vs nodal `A`" is **not** the closure discriminator (it is
  the spurious-eigenmode / interface-continuity one — the classical edge-element
  motivation). A first guess that "nodal curl leaks through normal jumps" is **wrong**:
  `(curl A)·n` depends only on the *tangential* trace of `A`, continuous for nodal `A`
  too.
- **The closure-breaker is leaving the de Rham representation — nodally *smoothing* `B`**
  (the legacy "evaluate at nodes and interpolate" post-processing). The smoothed field
  acquires a spurious divergence that **decreases with order but is never zero**
  (`1.2e-2 → 1.3e-3 → 6e-5` at `p = 1, 2, 3`) — so the flux lines do **not** close
  *even at 2nd order* (Q2: yes, matches Kameari). The de Rham `B`, kept native, is
  exactly div-free at every order — **a de Rham 2nd-order field closes** (Q2: yes).
- On a 2-D solve the de Rham `rot(grad A_z)` is exactly tangent to the flux surfaces
  (misalignment `0.0`) and closes at every order; the smoothed reconstruction's
  misalignment falls with order (`1.5e-1 → 2.5e-2 → 4.5e-3`) but stays far above
  de Rham — it does not close even at 2nd order.

So **closure is governed by the representation (de Rham vs nodal-smoothed), not the
order**. de Rham is the *closed-2-form precondition* that makes a symplectic (2-D,
Hamiltonian `A_z`) / volume-preserving (3-D, `div B = 0`) tracker meaningful (Q1: yes
— Noguchi supplies the field, structure-preserving integration supplies the tracker;
their union is the extension). *Refs: Robert 1991; Noguchi edge-FE flux lines;
Bossavit / Nédélec (edge elements); Moffatt (helicity).*

### `clebsch_3d_closing_condition.py` — the 3-D frontier: helicity is the obstruction

The 2-D flux line always has a conserved `A_z` (1.5 DOF → always integrable, always
closes). In **3-D** a *global* Clebsch pair `B = ∇α×∇β` need not exist, and the
fundamental obstruction is the **helicity** `h = ∫A·B` (the topological
linking/knotting of field lines — Moffatt 1969). A Clebsch field is helicity-free
**pointwise** (`A = α∇β` ⟹ `A·B = α∇β·(∇α×∇β) = 0`), so **`h ≠ 0` ⟹ no global
Clebsch pair** ⟹ flux lines are linked/chaotic and do **not** close. This is the
3-D "closing condition": *flux lines close (lie on flux surfaces) iff a global
Clebsch pair exists iff the helicity vanishes.*

Verified on the 3-torus: a Clebsch field `∇α×∇β` has `h = −7e-17` (machine zero,
rel `4e-19`) → global Clebsch exists; the **ABC** Beltrami field
(`∇×B = B`, `A = B`) has `h = ∫|B|² = 744.15 = 3(2π)³` exactly → no global Clebsch,
and a single ABC flux line's **Poincaré section fills a 2-D region** (occupancy
`0.55` — the famous ABC chaos, never closing onto a curve). This pins the open
frontier: the saturable 3-D HDiv-VIM solve admits the Chaplygin/Clebsch
linearisation only where the helicity obstruction vanishes (integrable / foliated
fields). *Ref (concept): Moffatt, "The degree of knottedness of tangled vortex
lines", JFM 1969; Arnold; Dombre et al., JFM 1986 (ABC chaos).*

## Design track — accelerator pole (the "iron face = equipotential" lever)

These quantify the design methodology (§3–5 of `DESIGN_METHODOLOGY.md`): the
iron pole face is a scalar-potential **equipotential**, and *deviation from it
= field harmonics*. They are analytic (no FEM), so each is machine-exact.

- **`accel_pole_design.py`** — the **multipole analyzer** `multipoles(B, r_ref)`
  (FFT → normal `b_n` / skew `a_n`) + the hodograph quad pole geometry
  (`Φ = G·xy` ⇒ pole `xy = r₀²/2`, hyperbola). Verified: a pure quad gives only
  `b₂` (spurious ~4e-17); an injected octupole is recovered to ~1e-15.
- **`accel_pole_harmonics.py`** — the **2-D cross-section lever**: the ideal
  hyperbola pole → pure quad; a shim off the equipotential → sextupole `a₆`
  grows ~linearly with the shim.
- **`accel_pole_ends_3d.py`** — the **3-D ENDS** (the July main goal, §3.2). The
  **integrated** transverse field `∫B⊥ dz` is *always* a 2-D multipole (it is
  2-D div- and curl-free), so the integrated analyzer is what beam optics sees.
  A **Maxwellian** (symmetry-preserving, equipotential-following) end has two
  honest properties: **(a)** the gradient's radial Maxwell corrections
  `∝ G″(z), G⁗(z)` are total z-derivatives → they **integrate to zero** → the
  integrated quad strength is *exact and radially undistorted*
  (`b̄₂ = (∫G)·r_ref` to ~5e-7); **(b)** the preserved m=2 symmetry generates no
  azimuthal `b₆` at all (`|b̄₆/b̄₂| ~ 3e-17`). A **non-equipotential** end
  **breaks the symmetry**, injecting a genuine one-signed `b₆` whose integral
  is nonzero → a spurious integrated `b̄₆` growing **linearly** with the
  deviation (slope spread 0.0%). So "follow the 3-D equipotential at the end" =
  *keep the symmetry and let the radial fringe corrections cancel, so the
  integrated field stays the pure designed multipole.* Figure:
  `accel_pole_ends_3d.png` (left: the radial correction ∝ G″(z) oscillates and
  integrates to zero; right: `b₆(z)` is flat for the equipotential end, a
  one-signed bump for the defect). **Rung 2 (FEM) below** replaces the analytic
  field with a real reduced-Ω solve.
- **`accel_pole_ends_fem.py`** — the **FEM rung** (§3.2, rung 2). Replaces the
  analytic Maxwellian field with a real **reduced-Ω + CoilBuilder** solve of a
  finite-length dipole — an x-symmetric H-frame iron yoke (netgen.occ, *no
  Cubit*) + a CoilBuilder racetrack pair as a Biot-Savart source (*no coil
  mesh*); `∫μ∇Ω·∇v = ∫μ Hₛ·∇v`, `H = Hₛ − ∇Ω` — and feeds the SAME integrated
  analyzer. Verified (finer mesh, ne ≈ 32 k): a **clean flat-top dipole**
  `B_z(body) ≈ 0.14 T` with `B_x/B_z ≈ 0.3 %` at centre (the x-symmetric H-frame;
  the residual shrinks with mesh refinement — 5.6 % → 0.3 %), effective length
  `L_eff ≈ 153 mm` (iron 120 mm + ~17 mm fringe each end), a small
  (mesh-sensitive ~3 %) pole-end enhancement, integrated dipole `b̄₁ ≈ 0.026
  T·m`, integrated spurious (n=3,5) ≈ 8 % (the ends + finite pole width — the
  part the next rung's re-shaping targets). Figure left: the on-axis `B_z(y)`
  (flat body + fringe).
  Two engineering notes baked in: **RadiaField (Biot-Savart) is not thread-safe
  under TaskManager** — the source LinearForm + field readout assemble serially,
  only the stiffness+solve are wrapped; and the CoilBuilder racetrack arcs curve
  toward −x (start at +R to centre it). It also **reads the equipotential as the
  end-iron contour**: in the current-free gap `H = −∇Ψ`, so along z the total
  scalar potential is `Ψ(y,z) = −∫₀ᶻ H_z dz'`; the iron face is the equipotential
  `Ψ_pole = Ψ(0, g/2)`, and the curve `z_p(y)` where `Ψ(y,z)=Ψ_pole` is the ideal
  end edge. Verified: in the body `z_p = g/2` *exactly* (self-consistency — the
  equipotential sits on the iron face), and it **lifts ~10 mm past the iron
  end** (the field bows out → the chamfer to follow). The figure's right panel is
  this contour. It then **closes the §3.2 design loop** (`end_shaping_sweep`):
  re-shape the pole END (chamfer it, following the equipotential lift) →
  re-solve → the longitudinal **pole-end enhancement is driven through zero**
  (optimal chamfer ~4–6 mm), shown in `*_loop.png`. **Honest two-lever finding:**
  the chamfer controls the *longitudinal* end bump, but the integrated
  *transverse* spurious `b₃,₅` (~9 %) barely moves — it is **body/pole-width
  dominated, not an end effect** (that is a Rogowski pole-shape problem, a
  separate lever). So end shaping is the right knob for the end bump, the wrong
  knob for `b₃,₅`.
  Finally it **follows `z_p(y)` *exactly*** with a **CURVED** end chamfer
  (`curved_chamfer_study`, figure `*_curved.png`) instead of the linear taper:
  the measured equipotential bow-out is a **parameter-free convex shape**
  `Ĝ(s)` (rises faster than linear near the body, then saturates); cutting the
  pole end along it drives the end bump from `+` through zero to negative, so
  the curved profile zeros it. Honest caveat: the **naive single-pass depth**
  (the raw equipotential lift, ~10 mm) **over-corrects** — the shape is right,
  the depth needs one knob (~25 % of the naive lift) or one design iteration
  (the straight-pole equipotential is read in the already-bumped field). The
  precise zero is mesh-noise-limited; the integrated transverse spurious stays
  body-dominated throughout (the end shape is not its lever).
  It also answers the **open-boundary** question (`open_boundary_convergence`,
  figure `*_openbnd.png`): growing the Dirichlet air box leaves the integrated
  dipole `b̄₁` **stable to 0.17 %** across 200 → 450 mm (< 1 % even out to
  600 mm), and *non-monotonic* (mesh noise, not a converging trend) — the spread
  *shrinks* with finer mesh (5 % at `maxh=0.035` → 0.17 % at `maxh=0.03`), the
  signature of a noise floor, not a real effect. The iron H-frame **flux return
  contains the field**,
  so the air-box truncation is **below the mesh-noise floor** — the open
  boundary is *not* the limiting error for this iron-dominated dipole, and an
  exact Kelvin would not change the answer. (`B_z`(body), a near-point value, is
  more mesh-sensitive, so its larger spread is mesh noise too.) The **exact**
  open boundary — which *does* matter for a flux-return-**free** magnet — is the
  lab's Kelvin transform, demonstrated exactly in `hodograph_kelvin_axisym.py`
  (above: `field_error ~1e-7` vs `~3e-3` truncated). A full 3-D *Cartesian*
  Kelvin for this *coil-driven* dipole would need the localized Biot–Savart
  source mapped into the Kelvin exterior (the closed-form
  reduced-potential-background helper covers uniform / dipole-at-infinity
  backgrounds, not a localized coil); it is deferred because the convergence
  here shows it would not change the result.
- **`accel_pole_dipole_body_2d.py`** — the **BODY lever** (the *other* half of
  the two-lever split): the transverse `b₃,₅` are a pole-**shape** knob, not an
  end effect. A 2-D cross-section solve (the body is translationally invariant
  → no ends) of the pole-as-equipotential Laplace problem. A finite **flat**
  pole droops at its edges (`b₃ < 0`); two body levers fix it: **(1) width** —
  a wider flat pole flattens the field over a fixed aperture (`|b₃|` falls
  ~100× from `w=30→60 mm`); **(2) curvature** — the dipole "shim"/Rogowski
  analog `z_face = g/2 − δ(x/w)²` (narrow the gap toward the edges to boost the
  drooped field) drives `b₃` *through zero* at `δ ≈ 0.4 mm`, leaving `|b₅|` as
  the residual (spurious `8e-4 → 1e-4`). This is the finite-aperture analog of
  the quad hyperbola: the ideal finite-width dipole pole is slightly concave,
  not flat — and it is the lever the end chamfer (above) cannot move.
- **`accel_quad_ends_fem.py`** — the **QUADRUPOLE FEM rung** (the analyzer
  handles *any* multipole). The dipole work built the method; this is its
  quadrupole instance — a real finite-length **4-pole hyperbola** quad solved by
  the forward engine and fed to the **same** integrated analyzer. Excitation =
  the 3-D scalar-potential high-μ model (the 3-D form of the 2-D body lever
  above): the four hodograph poles `xy = ±r₀²/2` (the `n=2` equipotential, from
  `accel_pole_design.py`) held at alternating `±Φ₀`, Laplace-solved in the
  current-free aperture; the finite iron bars give real ends. Verified (FEM,
  `ne ≈ 70 k`, refined axis core for clean harmonics): a **clean integrated
  quadrupole** — main `b₂`; the symmetry-**forbidden** normals `n=1,3,5`
  suppressed to `~6e-4` (≈ 5× below) the first **allowed** spurious, the 12-pole
  `b₆ ≈ 3e-3`. A length comparison (`b₆/b₂` changes only ~16 % over 2× length)
  shows `b₆` is **body-dominated** (the finite pole angular width), not the ends
  — the **same** two-lever picture the dipole found for `b₃,₅` (the ends drive
  the longitudinal profile; the body pole shape drives the transverse harmonics;
  the analytic theorem in `accel_pole_ends_3d.py` shows an ideal Maxwellian end
  contributes *zero* integrated `b₆`). Figure: `accel_quad_ends_fem.png` (the
  4-pole hyperbola geometry + the integrated harmonic spectrum).
- **`one_turn_coil_streamfunction.py`** — the A-side **(B) track**: a 1-turn
  coil is the coarsest stream-function discretization (one contour = one wire);
  the task is the single best wire path, and the script shows the honest 1-turn
  limit vs the full multi-turn stream-function current.
- **`clebsch_dipole_design_workflow.py`** — the **end-to-end workflow** that ties
  the track together: a 3-D dipole designed by a **Clebsch / scalar-potential
  LEVEL SET** carried from the 2-D cross-section into the 3-D magnet. One
  geometric object at three places — the iron pole **surface** is the
  equipotential `Ω = const` (the scalar potential conjugate to the flux function
  `A_z` of the flux-line work).
  **Stage A** (2-D, fast): the cross-section pole face `z_p(x)` is the level set;
  a finite flat pole droops (`b₃ < 0`), and a curvature shim
  `z_p = g/2 − δ(x/w)²` drives `b₃` through zero. The **width knob**: a wider
  pole needs less shim (`δ_opt` 1.12 → 0.04 mm as `w` 30 → 60 mm). This is the
  cheap, *exact* instrument for the transverse harmonics.
  **Stage B**: reflect the level set into the 3-D pole **surface** — body = the
  2-D contour **extruded** along the beam (body field = the 2-D field, `B =
  ∇A_z × ŷ`); end = the **equipotential / Maxwellian** end so the fringe
  integrates away.
  **Stage C** (3-D reduced-Ω FEM, slow): a real finite-length H-frame solve
  confirms a clean flat-top dipole (`B_z ≈ 0.16 T`, `B_x/B_z ≈ 7 %` coarse), a
  clean integrated dipole (`b̄₁ ≈ 0.027 T·m`, `L_eff ≈ 159 mm`), and reads back
  the equipotential **end contour** `z_p(y)` (= `g/2` in the body, lifting ~9 mm
  past the iron end = the Maxwellian end to cut). **Honest scope**: the
  transverse harmonics are designed *and* verified in 2-D (cheap, exact); the 3-D
  integrated transverse spurious at a golden-feasible mesh is **mesh-noise-limited
  (~8e-2)**, well above the cross-section's intrinsic `b₃` — so 2-D is the
  instrument for it, exactly how accelerator dipoles are designed (2-D
  cross-section optimisation + 3-D end correction). The design chooses the width
  so the residual shim is negligible (`δ_opt < 0.05 mm`), so the flat-body 3-D IS
  the reflected design; reflecting a non-negligible body shim and resolving its
  sub-mesh-noise integrated improvement is a fine-mesh study, deferred. Figure:
  `clebsch_dipole_design_workflow.png` (the cross-section level set | the width
  knob | the level set carried into 3-D). Goldens `test_clebsch_dipole_workflow_*`
  (design fast, FEM slow).
- **`clebsch_pole_shape_optimization_2d.py`** — **3-D Clebsch pole-face shape
  OPTIMIZATION**: null `b3` **and** `b5` together. The pole face is a magnetic-scalar
  **equipotential = a Clebsch level set**, and the 3-D body pole is its extrusion, so
  optimizing the 2-D contour optimizes the 3-D pole. A finite flat pole droops at its
  edges (`B_z = b1 + b3 x² + b5 x⁴ + …`); a **single** quadratic shim has one knob, so
  it nulls `b3` but **leaves `b5`** (`≈1.1e-4`). A **two-parameter** Clebsch contour
  `z = g/2 − d2(x/w)² − d4(x/w)⁴`, optimized by a **2-D Newton on `(b3, b5) = 0`** (the
  harmonics respond near-linearly to the shim coefficients, so the Jacobian solve lands
  close and one refinement tightens it), nulls **BOTH** leading spurious harmonics
  simultaneously: `flat 8.3e-4 → 1-shim 1.1e-4 (b5 stuck) → 2-shim 1.8e-6` (a **≈470×
  cleaner** field), at `d2 = −0.33 mm, d4 = 1.65 mm`. A genuine multi-parameter shape
  optimization, each eval a 2-D Laplace solve. Figure:
  `clebsch_pole_shape_optimization_2d.png` (the three contours | the `b3`/`b5` bars).
  Golden `test_clebsch_pole_shape_optimization_2d`.
- **`clebsch_dipole_saturation_2d.py`** — **saturation** added to the dipole design,
  at **linear cost**: the level-set workflow above assumed an *unsaturated* iron
  equipotential; at high excitation the iron **saturates** and the gap field stops
  rising linearly. The iron flux **return path is a saturable Chaplygin guide**, so the
  operating curve `B_gap(NI)` is a **magnetic-circuit 1-shot**
  `NI = (B_gap/μ₀)·gap + Σ_seg ν(B_seg)·B_seg·L_seg` — the iron sum *is* the Chaplygin
  reluctance integral, and a deliberately **necked segment (a throat,
  `A_throat<A_gap`)** is the Chaplygin throat embedded in the magnet: it carries
  `B_throat=B_gap·A_gap/A_throat` and **saturates first**, setting the knee. Given `NI`
  the circuit root-finds `B_gap` in microseconds (no mesh, no Picard). A **120-point**
  operating-curve design map (3 throat widths × 40 drives) is **~10 ms**; a **thinner
  iron throat saturates the magnet EARLIER** (knee drive `6.3→19.8 kA·t` as the throat
  goes `6→16 mm`) — the saturation design knob the (linear) pole-surface design cannot
  see. **Validated against a real 2-D nonlinear FEM** (an A-formulation Froehlich-iron
  window-frame electromagnet, under-relaxed Picard): the lumped circuit matches the FEM
  to **~6–20 %** (the lumped-circuit error — gap fringing + corner crowding), each FEM
  point costing **~22 Picard iterations**. So the whole map = 120 root-finds (10 ms) vs
  the equivalent **~2670 nonlinear FEM linear solves**. Figure:
  `clebsch_dipole_saturation_2d.png` (the `B_gap(NI)` curves + FEM points | the knee
  design knob | the cost bar). Goldens `test_clebsch_dipole_saturation_*` (circuit fast,
  FEM slow). *Honest scope: the lumped series circuit is the standard ~10 %-accurate
  engineering model; the iron-path reluctance is the Chaplygin 1-shot whose slender-guide
  form is validated in `chaplygin_hodograph_2d.py`, and the 2-D FEM here is the
  end-to-end check.*
- **`clebsch_dipole_saturation_3d.py`** — saturation in a **3-D** dipole, done **right**:
  the **B-input A-formulation**, the documented **cure** for the reduced-Ω
  ill-conditioning. The obvious 3-D choice — the reduced scalar potential
  `H = H_s − ∇Ω`, `μ(|H|)` Picard — is ill-conditioned at high permeability: the
  saturation knee sits in the **low-drive (unsaturated, high-μ)** regime, and there the
  reduced-Ω Picard **stalls** (spurious `|H|` in the iron from an under-resolved `Ω`).
  The fix is the convex **`ν(|B|)` B-input** form. The **reduced vector potential**
  avoids meshing the coil / any `div J`: `B = B_s + curl A_r` with `B_s = μ₀H_s` the
  coil's Biot–Savart field (Radia), `A_r ∈ H(curl)` the iron reaction, and the weak form
  `∫ν(|B|) curl A_r·curl v + ∫_iron (ν−ν₀) B_s·curl v = 0` (the coil drives the iron
  through its known `B_s`; a tiny `ε∫A_r·v` fixes the gauge). Verified on the accel
  H-frame + CoilBuilder coil (reused from `accel_pole_ends_fem.py`): **the cure** — at a
  low drive (`NI=3 kA·t`, iron `⟨μ_r⟩≈1200`) the B-input A-formulation converges to
  `resid 8e-6` in **15 iters**, while the reduced-Ω `μ(|H|)` Picard **stalls at
  `resid ~2e-2`**; **the saturation** — swept over drive, the A-formulation converges at
  *every* point (~14 iters, `resid<1e-5`) including the low-drive regime the reduced-Ω
  cannot reach, and the iron `⟨μ_r⟩` falls `1228 → 161` (it saturates). Figure:
  `clebsch_dipole_saturation_3d.png` (the convergence histories — the cure | `B_gap(NI)`
  + iron `⟨μ_r⟩(NI)`). Golden `test_clebsch_dipole_saturation_3d_aform` (slow). *Honest
  scope: this **large-gap** dipole is gap-reluctance-dominated (`R_gap/R_iron ~ 160`
  unsaturated), so `B_gap` softens only mildly with iron saturation (a large-gap dipole
  is intrinsically robust to it; a strong `B_gap` knee needs a small gap / necked iron,
  as in the 2-D throat geometry above) — the point demonstrated is the **solver**: the
  B-input A-formulation resolves the high-μ regime the reduced-Ω cannot.*

- **`clebsch_dipole_saturation_3d_throat.py`** — the **strong 3-D saturation knee** that
  B(a)'s large-gap dipole could not show, by **throat flux-concentration**. A
  rectangular iron **window frame** with a thin bottom-limb **throat** in series with
  the gap: series flux conservation gives `B_throat = B_gap·(A_gap/A_throat)`, so the
  throat (here `A_gap/A_throat = 2`) reaches `J_sat` **first**; its permeability then
  collapses (`μ_r→1`), its reluctance explodes, and the **circulating flux clamps**. The
  operating curve is a **magnetic-circuit 1-shot** (a scalar reluctance root-find per
  drive, no mesh) — the design map at **linear cost**, exactly as in the 2-D case. *The
  B(b) bug and its fix:* a throat only concentrates flux if the drive forces
  **circulating** flux around the closed window. A racetrack coil over the gap drives
  only the **local** gap flux, so the throat sees nothing (`B_throat ≈ B_gap`, no knee —
  the failed H-frame/C-core attempts). The cure is a coil that **links** the window — a
  `rad.ObjFlmCur` loop threading one limb — driving flux all the way around, through the
  throat. The convergent **B-input A-formulation** FEM (reused from B(a), with **adaptive
  under-relaxation** to cross the steep knee) validates it: a **linear** solve measures
  `B_throat/B_gap ≈ 3.8` (> the geometric area ratio 2 — the throat genuinely **funnels**
  flux, the fix proven), and the nonlinear sweep shows the **throat field** (the
  bottleneck) **clamp at ≈2.84 T** (`J_sat=2 T`) as the drive grows 8× — `dB_throat/dNI`
  drops **≈157×** (a strong, converged knee, `resid<1e-4`). *Honest 3-D caveat (a real
  result):* the **raw gap field does NOT clamp** — it keeps rising, because (1) the
  linking coil sits beside the gap, so `B_gap` carries the coil's local un-saturable
  field, and (2) once the throat saturates the circulating flux **leaks around** it in
  3-D (alternative return paths). So a single saturated throat clamps the field
  **locally** (at the bottleneck), not everywhere; the lumped circuit's full
  flux-conservation is the design-level approximation, the FEM is the truth. Figure:
  `clebsch_dipole_saturation_3d_throat.png` (left: the circuit design map — both fields
  bend at the knee | right: the FEM — `B_throat` clamps at ≈2.84 T while `B_gap` keeps
  rising, the honest 3-D leakage). Golden
  `test_clebsch_dipole_saturation_3d_throat_circuit` (fast, circuit) +
  `test_clebsch_dipole_saturation_3d_throat_fem` (slow). *This is the iron-bottleneck
  counterpart to B(a): B(a) demonstrated the **solver** on a gap-dominated dipole; B(b)
  demonstrates the **strong channeled-flux knee** on a throat-necked circuit, and the
  honest 3-D limit of throat clamping.*
  **3-D inverse design (spec → geometry):** `size_throat_for_knee(target)`
  sizes the **throat thickness** so the magnet's **saturation onset** (the knee drive)
  hits a target — the knee drive is **monotone in throat thickness**, so the 1-shot
  circuit bisects it **mesh-free, to machine precision** (e.g. 4000 A → `t = 7.1 mm`),
  the 3-D analog of the 2-D flux-limiter sizing. `run_inverse_3d(target, with_fem=True)`
  then **closes the 3-D loop**: the adaptive B-input A-formulation FEM at the *sized*
  throat shows `B_throat` crossing `J_sat` (a real knee), with the FEM onset **earlier
  than the circuit target by the flux-funneling factor** (`≈1.7` — the lumped-circuit
  correction, since the linking coil funnels more flux than the geometric area ratio;
  reported, not hidden). Goldens `..._inverse_sizing` (fast) + `..._inverse_fem` (slow).

### `scaling_ffag_pole_2d.py` — achromatic (scaling-FFAG) gantry pole, certified by A/phi bounds

The **engineering target the whole hodograph line serves**: a fixed-field
(scaling-FFAG) **proton gantry** delivers a *range* of beam momenta without
re-exciting the magnets (fast energy switching, a compact gantry). The enabling
property is **achromaticity** — momentum-independent tune — which for a scaling
field `B_y(r) = B0 (r/r0)^k` (orbits geometrically similar, `p ~ r^{k+1}`) is
the *single* condition that the **field index** `k(r) = d log B_y / d log r` be
**constant** across the aperture. The hodograph is the natural language: in
`u = log r`, momentum change is a *translation* and the scaling field is a
*straight line* (`log B` vs `log r`, slope `k`), so **achromatic ⇔ that line
stays straight**.

The design is **shape optimization** (fixed topology — the pole face is the
`Phi = const` equipotential / Clebsch level set), so the hodograph chart is
single-valued (no fold). For a **super-ferric** magnet (SC coils set the NI,
the iron pole shapes the field, saturation is the wall) the value is to hold
`k` constant *into saturation* — pushing the field for a smaller gantry.

**Step 1 here (linear):** the field-index metric (exact on `B ~ r^k`,
`~7e-14`); the proton 70-250 MeV band → radial aperture ratio `~1.12` (`k=5`);
and the achromaticity is **certified by the complementary A/phi bracket** — the
same gap solved as the scalar potential (`phi`, Dirichlet on the equipotential
poles, `B=-grad phi`) AND as the flux function (`A`, the dual BCs, `B=curl A`),
the energy↔co-energy (Legendre) conjugate pair. The two field indices converge
from discretisation-*complementary* sides; their **gap < 1e-6** certifies that
the measured `k(r)` is physics, not mesh. Result: the naive `g ~ r^{-k}` pole
gives a bulk index `k ≈ 4.88` (just under the design `k=5`) — a small, certified
2-D-fringing **deficit** that the (upcoming) Step 3 reshape closes; *honest*:
the bracket certifies the deficit is real. Figure
`scaling_ffag_pole_2d.png` (left: `|B_y|(r)` log-log vs ideal `r^k` | right:
`k_phi(r)`, `k_A(r)` and the A/phi bracket vs `k_design`). Golden
`test_scaling_ffag_pole_2d_step1`.

**Step 2 (saturation, `--step2`):** a super-ferric iron pole (SC coils set the
NI, the iron shapes the field) with a Froehlich `mu(B)` (knee `Bk=1.2 T`).
Driving the scalar potential through the air gap (`B = -mu0 mu_r(|B|) grad phi`,
Picard on `mu_r(B)`) and sweeping the drive: referenced to the lowest
(unsaturated) level, the field-index loss `Dk(r) = k - k_ref` isolates the
saturation. As the high-r edge (highest `B`) crosses the iron knee, `Dk`
develops a **negative dip THERE** — the achromaticity degrades at the
high-energy edge of the momentum acceptance (the **super-ferric operating
wall**). Verified monotone: at gap fields up to `~3.3 T` the high-r index drops
`~0.26` below baseline and the aperture **tilts** (`dk_tilt ~ -0.13`: the high-r
edge droops more than the low-r edge). Figure
`scaling_ffag_pole_2d_saturation.png` (left: `Dk(r)` deepening at the high-r
edge | right: `|B_gap|(r)` crossing `Bk`). Golden
`test_scaling_ffag_pole_2d_step2_saturation`.

**Step 3 (reshape, `--step3`):** restore a flat (achromatic) `k(r)` INTO
saturation by reshaping the pole face. In the log chart the gap correction is a
polynomial `g(r) = g0 exp(-k u - gamma/2 u^2 - gamma2/6 u^3)` (`u = log r/r0`),
so the local geometric index `k_geom(u) = k + gamma u + gamma2/2 u^2` has two
knobs — `gamma` cancels the **tilt** of `k(r)`, `gamma2` the **curvature** (the
single-valued von Mises / log-chart reshape; cf. the 2-param Newton in
`clebsch_pole_shape_optimization_2d` that nulls `b3` AND `b5`). A 2-D Newton on
`(tilt, curv)=0` runs on the **nonlinear (saturated)** solve, so it corrects the
ACTUAL operating field index. Result at the super-ferric design excitation
(`B@r0=1.8 T`): the naive pole's saturated `k(r)` varies `ptp ~ 0.089`
(tilt `+0.084`, curv `+0.039`); **one Newton step** (`gamma,-0.85`,
`gamma2,+42.9`) nulls both to `~0.002` and flattens the field index **7.2x**
(`ptp 0.089 -> 0.012`). Figure `scaling_ffag_pole_2d_reshape.png` (left: `k(r)`
naive vs reshaped vs `k_design` | right: the pole-face shape change). Golden
`test_scaling_ffag_pole_2d_step3_reshape`. *Together Steps 1-3 are the
saturation-robust achromatic scaling pole: certify the index (A/phi bracket),
measure the saturation droop, and reshape it flat — all in the hodograph
(log/von Mises) chart, single-valued shape optimization.*

## Run

```bash
python a_method_clebsch_2d.py                 # A-method net figure
python hodograph_kelvin_axisym.py             # Kelvin exact-open-boundary flux figure
python hodograph_kelvin_2d.py                 # Kelvin in the hodograph (2-D Cartesian, no air box)
python bidirectional_coordinate_transform_2d.py  # the coordinate transform unified: pullback W (Tampere) (--fig)
python weakform_pullback_kata.py              # weak form in forms -> hodograph weight W=diag(mu,1/mu) (kata)
python clebsch_kelvin_3d.py                   # Kelvin in the hodograph (3-D Clebsch, no air box)
python saturation_loop_2d.py                  # the nonlinear saturation loop (Chaplygin reference)
python chaplygin_design_sweep_2d.py           # nonlinear-as-linear: a saturable design space at linear cost (--fem)
python chaplygin_turning_design_sweep_2d.py   # ... one rung up: TURNING fields, one linear solve each (--fem)
python chaplygin_taper_design_sweep_2d.py     # ... free-boundary closure: TAPER sweep, Newton steps all linear
python cohomology_hodograph_currentlink.py    # when the hodograph needs cohomology
python accel_pole_design.py                   # multipole analyzer + quad pole geometry
python accel_pole_harmonics.py                # 2-D equipotential lever (shim → harmonics)
python accel_pole_ends_3d.py                  # 3-D ends: integrated analyzer + end rule
python accel_pole_ends_fem.py                 # FEM rung: reduced-Omega + CoilBuilder dipole
python accel_pole_dipole_body_2d.py           # the BODY lever: pole width + curvature -> b3,5
python accel_quad_ends_fem.py                 # the QUADRUPOLE FEM rung (any multipole)
python one_turn_coil_streamfunction.py        # (B) the 1-turn stream-function limit
python clebsch_dipole_design_workflow.py      # end-to-end: 2-D level set -> 3-D dipole (--fem for Stage C)
python clebsch_pole_shape_optimization_2d.py  # 3-D Clebsch pole shape opt: null b3 AND b5 (2-param Newton) (--fig)
python scaling_ffag_pole_2d.py                # achromatic scaling-FFAG gantry pole: k(r) index + A/phi bracket (--step2 saturation, --step3 reshape; --fig)
python clebsch_dipole_saturation_2d.py        # saturation in the dipole: iron flux-path B_gap(NI) at linear cost (--fem)
python clebsch_dipole_saturation_3d.py        # saturation in 3-D, done right: the B-input A-formulation (the cure)
python clebsch_dipole_saturation_3d_throat.py # B(b): the STRONG 3-D B_gap knee via throat flux-concentration (--fem)
python chaplygin_free_boundary_2d.py          # Frontier 2: the turning-guide free boundary (image)
python chaplygin_inverse_vonmises_2d.py       # Frontier 2 inverse: von Mises dissolves it (linear)
python chaplygin_inverse_nonlinear_2d.py      # Frontier 2 CLOSED: nonlinear inverse, flux (lambda) freed
python hdiv_vim_clebsch_loopstar.py           # de Rham capstone: HDiv-VIM loop modes ARE Clebsch fields
python hdiv_vim_clebsch_2d_az.py              # 2-D unification: A_z IS the Clebsch potential
python flux_line_closure_symplectic.py        # dynamical face: flux-line closure (de Rham + symplectic)
python flux_line_realfield_ngsolve.py         # dynamical face on a REAL FE field: de Rham closes, leaky spirals
python derham_closure_order_sweep.py          # order vs representation: de Rham closes at every order, smoothed does not
python clebsch_3d_closing_condition.py        # 3-D frontier: helicity obstructs the global Clebsch pair
```

Locked by `tests/feec/test_clebsch_hodograph_research.py` (26 tests; the nine
heavy FEM rungs — forward+contour, the design loop, the curved chamfer, the
open-boundary convergence, the quadrupole, the Chaplygin 1-shot-vs-loop, the
3-D nonlinear Kelvin merge, the turning-guide hodograph PDE, and the
free-boundary image — are `@pytest.mark.slow`; the von Mises inverses, linear
and nonlinear, are fast least-squares solves).

## Prior art (honest)

The 2-D bidirectional "potential coordinate" map is published — Dervisha,
Marjamäki, Rasilo, Tarhasaari, *"Bidirectional Coordinate Transformation and
Its Application to 2-D Magnetic Field Problems"*, CEFC 2026 (Bossavit
school). The 2-D complex-potential / conformal pole design lineage is older
still (Rogowski 1923, Halbach). These scripts are an independent NGSolve
realisation coupled to the lab's reduced-potential + Kelvin machinery — **not
a novelty claim.** See `memory/derham_field_current_duality.md` for the full
prior-art verdict.
