# Accelerator magnet design methodology — A–φ duality + reduced potential

A principled, open-source design framework for accelerator electromagnets
that gives **design insight into where to place the iron — especially at the
magnet ENDS** — and unifies the *iron-pole* and *coil* realizations of the
same target field.

This document fixes the overall picture (the "establish the method" step).
It cleanly separates what is **verified today** from the **research program**
(per the repo-first policy: the verified core is shipped; the frontier is
named, not claimed done).

---

## 1. Foundation — the A–φ duality on the de Rham complex

Magnetostatics has two complementary potential formulations, dual on the
de Rham complex `H1 → HCurl → HDiv → L2`:

| | φ-method (scalar) | A-method (vector) |
|---|---|---|
| potential | `Ω ∈ H1` (0-form), `H = H_s − ∇Ω` | `A ∈ HCurl` (1-form), `B = ∇×A` |
| conforming | H-conforming (curl-free H exact) | B-conforming (div-free B exact) |
| physical realization | **iron pole** = equipotential `Ω = const` | **coil/current** = flux / stream function |
| iron/coil BC | pole face: `Ω` Dirichlet | flux/current: `A` Neumann (dual BC) |

The Hodge star `★` (carrying `µ`) connects the two sides; in 2-D the **hodograph
net** `Φ ⟂ A_z` is exactly the conjugate pair (complex potential
`W = A_z + i µ₀ V`).

**Verified core:** `a_method_clebsch_2d.py` (the A-method dual, vector
potential primary), the panel's Ω-reduced solver, `hodograph_kelvin_axisym.py`
(exact open boundary).

---

## 2. Forward engine — reduced potential + CoilBuilder

The forward solve uses the **reduced potential**, with the coil supplied by
**CoilBuilder** (Biot–Savart source, **no coil mesh**):

- **Coil** → CoilBuilder wire path → `H_s` (Biot–Savart). The source enters
  the reduced potential; the coil is never meshed.
- **Iron** → reduced potential `Ω` (or `A`) FEM + **Kelvin open boundary**
  (no air-box truncation).
- **Exists today:** `src/radia/panels/calc_accel_magnet.py`
  (`--formulation omega|a`), driven from the `radia-electromagnet` panel.

The forward engine is therefore an existing, validated asset; the new work is
the *design loop* (§3) and the *end-field insight* (§3.2), not the solver.

---

## 3. Design principle — "the iron face is an equipotential"

On a high-permeability iron surface `H_t ≈ 0`, so the **iron face is an
equipotential of the total scalar potential `Ω`**. Design is the
self-consistent fixed point:

```
place iron  →  reduced-potential solve (CoilBuilder source)
            →  read the equipotential surfaces of Ω
            →  move the iron face onto the target equipotential
            →  iterate to the fixed point
fixed point ⇒ iron face = equipotential ⇒ clean field
```

The **deviation of the actual iron from the equipotential is the source of
field harmonics** — this is the quantitative design lever.

### 3.1 The 2-D cross-section (solved)

The equipotential is the classic pole profile, known in closed form:

| magnet | `Ω` (scalar) | pole face `Ω = const` |
|---|---|---|
| dipole (n=1) | `(B₀/µ₀) y` | `y = const` (flat) |
| quad (n=2) | `G x y` | `x y = r₀²/2` (hyperbola) |
| sextupole (n=3) | `Im(z³)·c` | `Im(z³) = const` |

Generated + measured by `accel_pole_design.py` (hodograph pole geometry +
the multipole analyzer).

### 3.2 The magnet ENDS — the research frontier (the main goal)

In 3-D at the magnet end the **equipotential surface curves**: the 2-D pole
profile no longer holds, and the *integrated* field `∫B dl` picks up
end harmonics. The methodology's central claim:

> **The 3-D equipotential surface of the reduced-potential solution tells you
> the optimal end-iron contour (chamfer / end shim).** Shape the end iron to
> follow the equipotential ⇒ the integrated-field harmonics are minimized.

This turns the conventionally *empirical / iterative* end-shimming into a
**principled equipotential-following design**. It is the part this method
aims to **establish** (track A, §5).

**Established analytically (rung 1, `accel_pole_ends_3d.py`, golden-tested).**
The claim has an exact backbone, provable with no FEM:

1. The **integrated** transverse field `B̄⊥(x,y) = ∫B⊥ dz` is *always* a 2-D
   multipole field — in the current-free aperture `∫(∂ₓBx+∂_yBy)dz = −∫∂_zBz dz
   = 0` and `∫(∇×B)_z dz = 0`, so `B̄⊥` is 2-D div- **and** curl-free. The
   *integrated* multipole analyzer (what beam optics actually sees) is therefore
   well defined regardless of the per-slice fringe.
2. A **Maxwellian** (symmetry-preserving, equipotential-following) end:
   **(a)** the gradient's radial Maxwell corrections `∝ G″(z), G⁗(z)` (the
   pseudo-multipoles, in the *same* m=2 channel) are **total z-derivatives** ⇒
   they integrate to zero ⇒ the integrated quad is *exact and radially
   undistorted* (`b̄₂ = (∫G)·r_ref` to ~5e-7); **(b)** the preserved m=2
   symmetry generates **no** azimuthal `b₆` (`|b̄₆/b̄₂| ~ 3e-17`, machine).
3. A **non-equipotential** end **breaks the m=2 symmetry**, injecting a genuine,
   one-signed `b₆` whose z-integral is nonzero ⇒ a spurious integrated `b̄₆`
   growing **linearly** with the deviation.

So "follow the 3-D equipotential at the end" = *keep the symmetry and let the
radial fringe corrections cancel, so the integrated field stays the pure
designed multipole.*

### 3.3 The radial field index — achromatic (scaling) vs isochronous, into saturation

§3.1–3.2 shape the iron so the *transverse / longitudinal* field is the
designed multipole. A **circular** machine (FFAG / cyclotron) adds an
orthogonal design axis: the **radial field index**
`k(r) = d log B_y / d log r` that controls how the field grows with radius.
This is a hodograph-native, single-valued *shape* design (no topology change),
golden in `scaling_ffag_pole_2d.py`.

**Two achromaticities — and they are different.** "Achromatic" (momentum-
independent) means one of two distinct, relativistically **mutually exclusive**
things:

| | scaling FFAG | isochronous (cyclotron / non-scaling FFAG) |
|---|---|---|
| invariant | betatron **tune** | revolution **time** |
| field law | `B_y(r) = B0 (r/r0)^k` | `<B>(r) = B0 γ(r)` |
| field index | `k = const` (rigid) | `k_iso(r) = (β γ)² = β²/(1−β²)` (**rising**) |
| chart | `u = log r`: log B straight, slope `k` | `u = log r`: a **convex rising** curve |

In `u = log r` a momentum scaling `r → λr` is a **translation**; the scaling
field is translation-covariant (hence `k = const`), while the isochronous field
deliberately **breaks** that symmetry (rising `k_iso`). The pole gap is
`g(r) ∝ 1/B(r)` in either case (thin-gap, `B ∝ 1/g`).

**The super-ferric wall = the nonlinear END PACK.** With a Froehlich `µ(B)`
iron pole, the high-`r` edge carries the highest `B`, saturates first, and the
achieved field index **droops** there — degrading achromaticity at the
high-energy edge of the momentum acceptance. For the **isochronous** magnet
this is most acute: the high-`r` end must deliver the *steepest* rise exactly
where the iron gives out.

**The fix is the same hodograph machinery in both cases.** A **2-parameter
pole reshape** in the log / von Mises chart (`g = g0 exp(−k u − γ/2 u² −
γ2/6 u³)`, local index `k_geom(u) = k + γ u + γ2/2 u²`) — single-valued (the
full 2-variable hodograph folds once `µ = µ(q)`, so the von Mises single-
variable chart is used) — drives the **saturated** index back onto the target:

- **scaling** (`run_step3`): target `k = const`; a 2-D Newton on
  (tilt, curvature) = 0 flattens the saturated index **~7.2×** in one step
  (`test_scaling_ffag_pole_2d_step3_reshape`).
- **isochronous** (`run_isochronous`): target the **rising** `k_iso(r)`;
  the *same* Newton drives the saturated `<B>(r)` back onto `B0 γ(r)`,
  restoring isochronism **3.1×** (field-shape residual `|<B>/(B0 γ) − 1|`
  2.3 % → 0.73 % at `B_gap ≈ 1.33 T > Bk = 1.2 T`)
  (`test_scaling_ffag_pole_2d_isochronous`).

**Certified into saturation (the A–φ bracket of §1, nonlinear).** The same
operating point is solved both ways — φ (Dirichlet on the poles) and A
(Dirichlet on the flux walls, driven to the φ-solve's median flux so both sit
at the **same** saturation state). Monotone `BH ⇒ convex energy ⇒` the energy
bracket survives into saturation (Synge hypercircle / Rikabi–Bryant–Freeman);
`k_φ(r)` and `k_A(r)` converge from discretisation-complementary sides, so a
tight gap (≈5e-4) certifies the saturated index is **physics, not mesh** — for
the scaling *and* the isochronous reshaped pole.

**Hodograph AS the solver (no remesh).** For the linear pole, `run_pullback`
solves on a **fixed** computational mesh with the pole shape entering as a
pullback deformation (`mesh.SetDeformation`, weight `W = |det J|(JᵀJ)⁻¹`), so
a reshape is a new *weight* on the same mesh — Netgen runs **once** for the
whole shape sweep (the genuine no-remesh win; reproduces the physical-remesh
`k(r)` to ~5e-4).

**Honest scope.** The reshape residual (0.73 % isochronous) is the higher-order
mismatch a *2-parameter* quadratic reshape leaves against a ~5× rising `k_iso`
— more shape DOF closes it; the Newton itself converges in one step. This is
the **radial `<B>(r)` isochronism only**; the AVF flutter (vertical focusing),
the betatron tunes, and the orbit↔field self-consistency are separate problems,
not modeled here. The no-remesh pullback is shown for the *linear* pole; wiring
it through the nonlinear saturated Newton is the next rung.

### 3.4 Foliate-and-perturb: when is the body 2-D and the end a perturbation?

§3.1 designs the 2-D cross-section; §3.2 designs the 3-D end. The **quantitative
bridge** is: slice the magnet into 2-D `(x,z)` **leaves** along the beam `y`,
**stack** the cross-section solution (0th order), and **connect** leaves by a
beam-direction **perturbation** (the ends). When does this land? Measured on a
real reduced-Omega + CoilBuilder dipole, parametrised by the iron length so the
aspect ratio `L/gap` can be swept (`leaf_coupling_perturbation_3d.py`):

- the **0th-order leaf-stacking error** `delta(y) = ||B_perp(.,y) -
  B_perp(.,body)|| / ||B_perp(.,body)||` (for a straight constant-gap magnet the
  body slice IS the 2-D infinite-long leaf) is ~0 in the body and grows at the
  ends;
- the **inter-leaf perturbation parameter** `eps(y) = (g/2)|dBz/dy|/|Bz_body|`
  (transverse / beam-variation scale; **not** an operator-norm ratio, which is
  trapped at 1 by `grad_perp^2 = -d^2/dy^2` in current-free air) is ~0 in the
  body, O(1) at the ends;
- the integrated **fringe excess** `(L_eff - L_iron)/L_iron` is the 1st-order
  correction, and it **scales as ~ gap/L** (log-log slope **-0.95**):

| L/gap | 2 | 3 | 5 | 8 |
|---|---|---|---|---|
| fringe excess | +180 % | +111 % | +70 % | +48 % |

**Consequence for the design.** A **compact** magnet (the §3.2 end-study dipole,
`L/gap = 3`) is **non-perturbative** (+111 % fringe, the 0th-order stack misses
~40 %, the 3-D-ness is *not* end-localised) -- you cannot foliate it, the ends
are the whole magnet. Foliate-and-perturb lands only for **long** magnets, the
fringe dropping to ~10 % near `L/gap ~ 40` (typical beamline dipole). **Where it
lands, the body is a 2-D cross-section design and only the ends need the 3-D
treatment -- and that 3-D end treatment is exactly §3.2 (follow the
beam-referenced equipotential surface).** An equipotential-following end removes
the fringe's *harmonic* contamination (§3.2 theorem) but **not** the fringe
itself (`L_eff > L_iron` is a free-space effect, the table above) -- so the end
fixes the *integrated strength*, the body 2-D design fixes the *field quality*.

**Established (rung 2, FEM, `accel_pole_ends_fem.py`, golden-tested).** The
analytic field is replaced by a real **reduced-Ω + CoilBuilder** forward solve
of a finite-length dipole — x-symmetric H-frame iron (netgen.occ, no Cubit) + a
CoilBuilder racetrack Biot-Savart source (no coil mesh); `∫μ∇Ω·∇v = ∫μ Hₛ·∇v`,
`H = Hₛ − ∇Ω` — fed to the **same** integrated analyzer. It reproduces a clean
flat-top dipole (`B_z ≈ 0.14 T`, `B_x/B_z ≈ 0.3 %` at centre on a refined mesh),
the effective magnetic length `L_eff > L_iron` (the two fringes), a small
(mesh-sensitive) pole-end enhancement, and the integrated dipole + spurious
(n=3,5 ≈ 8 %, the ends + finite pole width). This is the
forward-engine bridge "analytic ⇒ FEM" and exposes the ends the design step
acts on. It also **reads the equipotential as the end-iron contour**: in the
current-free gap `H = −∇Ψ`, so `Ψ(y,z) = −∫₀ᶻ H_z dz'`; the iron face is
`Ψ_pole = Ψ(0, g/2)`, and `z_p(y)` with `Ψ(y,z)=Ψ_pole` is the ideal end edge —
verified to recover `z_p = g/2` *exactly* in the body (self-consistency) and to
**lift ~8.5 mm past the iron end** (the field bows out ⇒ the chamfer to follow).
And it **closes the loop**: re-shaping the pole END (a chamfer following the
equipotential lift) and re-solving drives the *longitudinal* pole-end
enhancement **through zero** (optimal chamfer ~4–6 mm). **Honest two-lever
result:** the chamfer controls the longitudinal end bump, but the integrated
*transverse* harmonics `b₃,₅` (~9 %) are **body/pole-width dominated** and barely
move — end shaping is the right lever for the end bump, a Rogowski body-pole
shape is the lever for `b₃,₅`. The §3.2 loop is closed for the end-field; the
transverse-harmonic lever is a separate (body) problem.

---

## 4. Dual realization — iron (φ) ⟷ coil (A)

The same target field has two realizations, bridged by the A–φ duality:

- **φ-side → iron pole** = the equipotential surface (hodograph).
- **A-side → coil/current** = the flux function / stream function
  (e.g. the quad `A_z = (G/2) r² cos 2θ` ⇒ a `cos 2θ` current sheet; the
  superconducting `cos nθ` coil).

So a field designed once yields **iron-pole, coil, or hybrid** realizations —
choose iron for warm DC high-field, the `cos nθ` coil for
superconducting / tunable, or split the field between them (the A–φ
decomposition) for a trimmed hybrid.

This **unifies the lab's two design lines**:

- φ-side iron pole — the hodograph (`examples/clebsch_hodograph/`).
- A-side coil current — the stream function
  (`src/radia/stream_function.py`, `examples/vim/foliated_solenoid_wires.py`,
  verified against Radia to 3.4e-10).

Multiply-connected / current-linking cases (a coil window) bring in the dual
**cohomology** roles (scalar needs a cut, vector needs a harmonic form) —
handled by `radia.cohomology` (gmsh-free), see
`cohomology_hodograph_currentlink.py`.

---

## 5. Two parallel tracks

The method advances on two **independent** tracks (the user's "(A) and (B)
are parallel"):

- **(A) Iron placement at the ends** — reduced potential + CoilBuilder +
  the equipotential-following design loop (§3.2); quantify with an
  **integrated** multipole analyzer. Deliverable: the 3-D equipotential map
  as the end-iron design rule.
- **(B) 1-turn coil via the stream function** — the A-side. A single-turn
  coil is the *coarsest* stream-function discretization (one contour = one
  wire), so a pure field is impossible; the task is the **single best wire
  path**, solved by field-RMS minimization (the `single-stroke` /
  `foliated_solenoid_wires` machinery). Deliverable: a designed 1-turn coil
  and the honest limit of what one turn can achieve.

Both tracks read out of the **same** A–φ Clebsch structure (iron = φ-side,
coil = A-side), so the framework is one method, not two.

---

## 6. Status (honest)

**Verified today** (golden-tested in `tests/feec/`):
- hodograph 2-D / axisym (cylinder, sphere) — bidirectional consistency;
- A-method dual (vector-potential primary) — Cauchy–Riemann conjugate net;
- Kelvin exact open boundary (field_error ~1e-7);
- cohomology current-linking (radia.cohomology, gmsh-free);
- the multipole harmonic analyzer (machine-precision);
- the **3-D end-field *integrated* multipole theorem** (§3.2 rung 1,
  `accel_pole_ends_3d.py`): the integrated field is always a 2-D multipole; a
  symmetry-preserving (equipotential) end ⇒ exact integrated quad (radial
  corrections ∝G″ integrate away) with no `b̄₆` (~3e-17); a non-equipotential
  end ⇒ spurious `b̄₆` linear in the deviation;
- the **3-D end-field FEM rung** (§3.2 rung 2, `accel_pole_ends_fem.py`): a real
  reduced-Ω + CoilBuilder finite-length dipole (netgen.occ, no Cubit) reproduces
  a clean flat-top dipole (`B_x/B_z ≈ 0.3 %`) + `L_eff > L_iron` through the same
  integrated analyzer, reads the solved equipotential as the end-iron contour
  (body `z_p = g/2` exactly; lifts ~10 mm past the iron end), and **closes the
  design loop** — a chamfer following that lift drives the longitudinal pole-end
  enhancement through zero (the transverse `b₃,₅` stays body-dominated);
- the stream-function coil (A-side), vs Radia to 3.4e-10;
- the reduced-potential + CoilBuilder forward engine (the panel);
- the **radial field-index design** (§3.3, `scaling_ffag_pole_2d.py`): the
  achromatic *scaling* pole (`k ≈ 4.88`, A/φ bracket ~9e-7) and its saturation
  droop + 2-param reshape (flat `k` restored ~7.2×); the **isochronous** variant
  (rising `k_iso(r) = (β γ)²`) whose nonlinear END PACK is driven back onto
  `B0 γ(r)` (saturation-broken 2.3 % → 0.73 %, **3.1×**, A/φ-certified ~5e-4);
  and the hodograph-as-solver pullback (fixed mesh, Netgen runs once);
- the **foliate-and-perturb scaling** (§3.4, `leaf_coupling_perturbation_3d.py`):
  the inter-leaf coupling (fringe excess) of a straight dipole decays as ~ gap/L
  (log-log slope -0.95 over L/gap = 2..8), so a compact magnet (L/gap=3) is
  non-perturbative (+111 % fringe) and the body-2-D + end-perturbation scheme
  lands only for long magnets (~10 % fringe at L/gap ~ 40).

**Research program (named, not claimed done):**
- the end-design loop is **closed for the longitudinal end-field** (§3.2 rung 2);
  what remains is the **transverse** integrated-harmonic lever — a Rogowski /
  body-pole-shape problem (a separate, body-not-end design knob), plus a
  *curved* (not linear) end chamfer that follows the equipotential contour exactly;
- the 1-turn coil stream-function design (§5 B);
- nonlinear µ(B) (saturation) inside the potential framework — the elegant
  "design sophistication", done within the reduced potential / hodograph
  (NOT via the volume-integral HDiv-MMM, which is the wrong tool here:
  it solves the iron's volume magnetization, a different part of the
  de Rham complex, and is 3-D-tet where this design is 2-D-native +
  thin-end).

---

## 7. Assets

| Role | Code |
|---|---|
| φ-side iron pole (hodograph) | `examples/clebsch_hodograph/` |
| A-method dual | `examples/clebsch_hodograph/a_method_clebsch_2d.py` |
| exact open boundary | `examples/clebsch_hodograph/hodograph_kelvin_axisym.py` |
| cohomology (current-linking) | `src/radia/cohomology.py` |
| pole geometry + multipole analyzer | `examples/clebsch_hodograph/accel_pole_design.py` |
| 3-D ends: integrated analyzer + end rule | `examples/clebsch_hodograph/accel_pole_ends_3d.py` |
| 3-D ends: FEM rung (reduced-Ω + CoilBuilder) | `examples/clebsch_hodograph/accel_pole_ends_fem.py` |
| forward (reduced potential + CoilBuilder) | `src/radia/panels/calc_accel_magnet.py` |
| radial field index (scaling + isochronous, saturation) | `examples/clebsch_hodograph/scaling_ffag_pole_2d.py` |
| foliate-and-perturb scaling (leaf coupling ~ gap/L) | `examples/clebsch_hodograph/leaf_coupling_perturbation_3d.py` |
| A-side coil (stream function) | `src/radia/stream_function.py`, `examples/vim/foliated_solenoid_wires.py` |
