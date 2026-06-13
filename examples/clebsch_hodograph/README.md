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
Clebsch flux `ψ` field lines on a meridional slice). *Next rungs: 1.5 = 2-D
Chaplygin (the hodograph **linearises** the saturation nonlinearity —
Molenbroek–Chaplygin); 3 = the 3-D merged geometry+material single Picard.*

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
python cohomology_hodograph_currentlink.py    # when the hodograph needs cohomology
python accel_pole_design.py                   # multipole analyzer + quad pole geometry
python accel_pole_harmonics.py                # 2-D equipotential lever (shim → harmonics)
python accel_pole_ends_3d.py                  # 3-D ends: integrated analyzer + end rule
python accel_pole_ends_fem.py                 # FEM rung: reduced-Omega + CoilBuilder dipole
python accel_pole_dipole_body_2d.py           # the BODY lever: pole width + curvature -> b3,5
python accel_quad_ends_fem.py                 # the QUADRUPOLE FEM rung (any multipole)
python one_turn_coil_streamfunction.py        # (B) the 1-turn stream-function limit
```

Locked by `tests/feec/test_clebsch_hodograph_research.py` (15 tests; the five
FEM rungs — forward+contour, the design loop, the curved chamfer, the open-
boundary convergence, and the quadrupole — are `@pytest.mark.slow`).

## Prior art (honest)

The 2-D bidirectional "potential coordinate" map is published — Dervisha,
Marjamäki, Rasilo, Tarhasaari, *"Bidirectional Coordinate Transformation and
Its Application to 2-D Magnetic Field Problems"*, CEFC 2026 (Bossavit
school). The 2-D complex-potential / conformal pole design lineage is older
still (Rogowski 1923, Halbach). These scripts are an independent NGSolve
realisation coupled to the lab's reduced-potential + Kelvin machinery — **not
a novelty claim.** See `memory/derham_field_current_duality.md` for the full
prior-art verdict.
