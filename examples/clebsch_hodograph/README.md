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

## Run

```bash
python a_method_clebsch_2d.py                 # A-method net figure
python hodograph_kelvin_axisym.py             # Kelvin exact-open-boundary flux figure
python hodograph_kelvin_2d.py                 # Kelvin in the hodograph (2-D Cartesian, no air box)
python clebsch_kelvin_3d.py                   # Kelvin in the hodograph (3-D Clebsch, no air box)
python saturation_loop_2d.py                  # the nonlinear saturation loop (Chaplygin reference)
python cohomology_hodograph_currentlink.py    # when the hodograph needs cohomology
python accel_pole_design.py                   # multipole analyzer + quad pole geometry
python accel_pole_harmonics.py                # 2-D equipotential lever (shim → harmonics)
python accel_pole_ends_3d.py                  # 3-D ends: integrated analyzer + end rule
python accel_pole_ends_fem.py                 # FEM rung: reduced-Omega + CoilBuilder dipole
python accel_pole_dipole_body_2d.py           # the BODY lever: pole width + curvature -> b3,5
python accel_quad_ends_fem.py                 # the QUADRUPOLE FEM rung (any multipole)
python one_turn_coil_streamfunction.py        # (B) the 1-turn stream-function limit
```

Locked by `tests/feec/test_clebsch_hodograph_research.py` (19 tests; the eight
FEM rungs — forward+contour, the design loop, the curved chamfer, the open-
boundary convergence, the quadrupole, the Chaplygin 1-shot-vs-loop, the
3-D nonlinear Kelvin merge, and the turning-guide hodograph PDE — are
`@pytest.mark.slow`).

## Prior art (honest)

The 2-D bidirectional "potential coordinate" map is published — Dervisha,
Marjamäki, Rasilo, Tarhasaari, *"Bidirectional Coordinate Transformation and
Its Application to 2-D Magnetic Field Problems"*, CEFC 2026 (Bossavit
school). The 2-D complex-potential / conformal pole design lineage is older
still (Rogowski 1923, Halbach). These scripts are an independent NGSolve
realisation coupled to the lab's reduced-potential + Kelvin machinery — **not
a novelty claim.** See `memory/derham_field_current_duality.md` for the full
prior-art verdict.
