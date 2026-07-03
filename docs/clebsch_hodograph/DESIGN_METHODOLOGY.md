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

### 3.5 The two-plane → 3-D method (the FFAG sector cell)

§3.1–3.4 are fragments of one method, stated plainly: **design the magnet in
two orthogonal 2-D planes and REFLECT the two designs into one 3-D pole.** For a
circular (FFAG / cyclotron sector) magnet the beam orbit is an azimuthal arc, so
the two planes are:

| plane | what it sets | engine |
|---|---|---|
| **transverse `(r, z)`** (⊥ orbit) | the field the beam **sees** — the scaling field index `B_z(r) ∝ r^k`, gap `g(r)=g₀(r/r₀)^{−k}` (§3.3) | `scaling_ffag_pole_2d.py` (Plane A) |
| **azimuthal `(s, z)`** (along the orbit, `s=r₀θ`) | how the field **turns on/off** — the sector ENDS, the effective magnetic length `L_eff=∫B_z ds / B_z(body)` | `ffag_sector_two_plane.py::solve_azimuthal_end` (Plane B) |

The **3-D pole is the `(r,z)` gap profile SWEPT around the sector arc**
`θ∈[−Δθ/2,Δθ/2]` and truncated at the azimuthal ends shaped by Plane B. The two
2-D designs are **exact in the body** and couple only at the ends; the §3.4
leaf-coupling law (`coupling ~ gap/L`) is therefore the **validity lever of the
reflection** — now on the AZIMUTHAL plane.

**Verified (rung 1, `ffag_sector_two_plane.py`, ngsolve only, golden-tested).**

- Plane A: the scaling index `k ≈ 4.88` (vs design 5; the naive pole's droop,
  reshaped by §3.3), A/φ bracket `~2e-6` (physics, not mesh).
- Plane B: a finite-length iron pole of gap `g(r₀)` over the sector arc
  `L_sector=r₀Δθ`, scalar-potential reduced solve. **Each sector end adds
  ~`0.75 g` of effective length**, i.e. `L_eff = L_sector + ~1.5 g`, so the
  **fringe excess `(L_eff−L_sector)/L_sector` falls as `~gap/L` (log-log slope
  −0.98)** — the same leaf-coupling law as §3.4's straight magnet, now azimuthal:

  | L_sector/g | 2 | 3.5 | 6 | 10 |
  |---|---|---|---|---|
  | fringe excess | +73% | +43% | +25% | +15% |

- **Validity:** the aspect `L_sector/g(r₀)` is the design lever. A **compact**
  cell (`L/g=3` → +50% fringe) is non-perturbative: the sector ENDS are a genuine
  3-D problem, not a body-stack correction. The two-plane reflection is exact only
  as `L/g → ∞`.

**Established (rung 2, the 3-D reflection — `ffag_sector_two_plane.py --rung2`,
ngsolve, golden-tested).** Sweep `g(r)` around the sector arc (**revolve** the
`(r,z)` gap cross-section about the bend axis) into a 3-D iron pole; drive it as
an **iron-pole equipotential** (upper-half model, median `z=0` the up-down
antisymmetry plane, the high-μ pole back at `Ψ=mmf` — the 3-D form of Plane A/B's
scalar potential). The orbit sees **`B_z(r) ∝ r^k` recovered: field index mean
`≈ 4.88`** (range `[4.6, 5.2]`, vs design 5 — the same naive-pole droop as Plane A
+ mesh scatter), confirming the swept `g(r) ∝ r^{−k}` pole reproduces the designed
radial field in full 3-D; the azimuthal sector ends add `L_eff/L_sector − 1 ≈
+39%` (the Plane-B fringe, cross-checked in 3-D). The field index is set by the
pole **geometry** (a high-μ equipotential forces `B ∝ 1/g(r)`), so it is
**drive-independent** — a CoilBuilder + reduced-Ω coil drive would set the field
**amplitude** (a further step), not the index.

**Honest scope.** Plane B is the linear-iron azimuthal-END geometry (the
magnetic-length excess; saturation is Plane A's §3.3 lever, composed
orthogonally). The radial profile is `<B>(r)` only — AVF flutter (vertical
focusing), the betatron tunes, and the orbit↔field self-consistency are separate
(as in §3.3). This is the **method scaffold**: two 2-D hodograph-native designs +
a measured `L/g` reflection criterion, now verified end-to-end (rung 1 the two 2-D
planes, rung 2 the 3-D reflection); the curved-orbit twist (combined-function, the
beam-referenced equipotential surface rotating along a bent orbit) is the next
axis.

### 3.6 The beam-referenced equipotential surface as the design primitive — and the twist

§3.1–3.5 read the multipoles OUT of a solved field. The reframing makes the
**beam-referenced equipotential surface the design SPEC**: in the Frenet frame of
the orbit `s`, the iron pole face is

$$\Omega(r,\theta;s) = \sum_n r^n\, b_n(s)\,\sin\!\big(n\theta + \phi_n(s)\big)$$

(high-μ ⇒ `H_tangential = 0` ⇒ `Ω = const`), so the **multipole `(b_n,φ_n)(s)`
IS the surface's angular Fourier mode.** Design = prescribe `(b_n,φ_n)(s)`, sweep
the equipotential surface along the orbit, place iron there — *not* solve-then-expand.

**The twist (the curved-orbit / combined-function axis).** The genuinely 3-D
content is that the transverse multipole **rotates** along `s` — the Frenet frame
turns with the bend, or a rotating-gradient magnet turns the pole on purpose. The
key fact is the **n-fold law**: rotating the equipotential SURFACE by `φ` rotates
the order-`n` multipole PHASE by `nφ`. For the quadrupole (`n=2`):

$$\text{rotate the pole by }\phi \iff (b_2,a_2)\to|b_2|(\cos 2\phi,\ \sin 2\phi),$$

so a quad twisted by 45° becomes a pure **skew** quad.

**Verified (`twisting_quadrupole_pole.py`, ngsolve only, golden-tested).** The
quad pole face is the hyperbola `xy = ±r0²/2` (the `Ω=const` equipotential); a 2-D
Laplace solve in the aperture with the 4 hyperbola poles at alternating `±Ω0`
recovers a **clean quad** (skew `a_2/|c_2| ~ 5e-6`, forbidden `n=1,3,5` at the
`~5e-5` floor, the leading allowed spurious the finite-pole 12-pole `b_6 ~ 5.6e-3`).
**Rotating the poles by `φ` rotates the recovered pole orientation by exactly `φ`**
— `α = −½ atan2(a_2,b_2)` tracks the prescribed twist to **slope 1.000, max error
0.00°** — and `b_6` is rotation-invariant. The twist is the surface angular mode,
measured.

**Honest scope.** This is the per-station (Frenet cross-section) 2-D design — the
**slow-twist (adiabatic) limit** `dφ/ds → 0`, where the magnet is a stack of 2-D
leaves (the §3.4 foliate-and-perturb picture, now twisting). A fast twist / tight
bend couples adjacent leaves (a longitudinal-field correction) — the twist rate
`dφ/ds` is a leaf-coupling perturbation parameter (the next rung). The
combined-function (dipole + quad together = a shifted+rotated hyperbola) and the
genuine curved-orbit Frenet sweep are the extensions; the quad here establishes
the n-fold twist law that governs them.

### 3.7 The confluence — a combined-function magnet on its curved orbit (the Frenet sweep is the twist)

§3.5 bent the beam with a pure **dipole** sector; §3.6 twisted a pure **quad** on
a fixed station. **`combined_function_frenet_sweep.py` merges them**: a
**combined-function** magnet (dipole `b1` + quad gradient `b2` in ONE
cross-section) swept along the **curved orbit it bends**.

In the **Frenet frame** of the orbit the cross-section is FIXED (the design spec:
`b1` bends, `b2` focuses). But the Frenet frame **rotates** with the bend — by the
bend angle `θ(s) = s/ρ`, `ρ = (Bρ)/b1` — so in the **lab frame** the whole
combined-function pole **twists by `θ(s)`**. §3.6's n-fold law then gives:

$$\text{geometric roll }\theta \;\Rightarrow\; \text{dipole phase }\psi_1=\theta,\quad \text{quad phase }\psi_2=2\theta,$$

so **both** the dipole and quad orientations track the *same* Frenet angle `θ`
(the rigid roll), while their multipole **phases** differ by the factor `n`
(`ψ_2 = 2 ψ_1`).

**The combined-function cross-section** is a **tilted-gap** dipole (`z = ±(g/2 −
t x)`): the gap narrows toward `+x`, so `B_z(x) ∝ 1/g(x)` carries a dipole `b1`
plus a gradient `b2` (the quad) — §3.5's flat gap, *tilted*.

**Verified (`combined_function_frenet_sweep.py`, ngsolve only, golden-tested).**
The 2-D Laplace solve recovers the combined function (`b1` dipole + `b2/b1 ≈ 6%`
gradient + a small `b3 ~ 3e-3`, the `1/g` curvature a real magnet shims out);
rolling the magnet by `θ`, **both the dipole and quad orientations track `θ`
(slope 1.000, error 0.00°)** and the **quad multipole-phase change is exactly
2× the dipole's** (`ψ_2/ψ_1 = 2.000`, the n-fold law) — the design-primitive
surface (fixed in the Frenet frame) reflected into a lab pole that twists by `θ(s)`.

**Honest scope.** A pure sector (rigid Frenet roll; no spiral edge, no s-varying
gradient), per-station 2-D = the **slow-bend** limit. A spiral sector (the pole
twist `φ ≠` the orbit bend `θ`) and an s-ramped `(b1(s),b2(s))` are extensions;
**when the per-station 2-D breaks** — the fast-twist `dφ/ds` leaf coupling
(the §3.4 perturbation parameter, now on the twist) — is the next rung.

### 3.8 When does the per-station 2-D twist break? — the fast-twist leaf coupling

§3.6–3.7 design a twisting magnet as a STACK of independently-rotated 2-D
cross-sections — the **slow-twist (adiabatic)** limit `dφ/ds → 0`. **When does
that break?** A magnet whose order-`n` multipole twists at rate `k = dφ/ds` has
helical symmetry, and the exact current-free harmonic is the **helical multipole**
(standard for helical undulators / Siberian snakes / twisted quads):

$$\Phi_n = C\, I_n(n k r)\,\sin\!\big(n(\theta - k s)\big),$$

`I_n` the modified Bessel function. As `k → 0`, `I_n(nkr) → (nkr/2)^n/n!` and the
field reduces to the pure 2-D multipole rotated by `φ(s) = ks` — **exactly the
per-station 2-D stack** (§3.6–3.7). The leaf coupling is the deviation from that
stack, controlled by the dimensionless **twist-per-aperture** `ka = 2π a/P`
(`a` aperture, `P` pitch).

**Measured (`twist_rate_leaf_coupling.py`, analytic, golden-tested).** On the
aperture circle:
- the **transverse** focusing error `ε = ‖B_⊥(3D) − B_⊥(2D stack)‖/‖B_⊥(2D)‖`
  scales as **`(ka)²`** (slope 2.04 — 2nd order, the quality the optics sees);
- the **longitudinal** field `B_s/B_⊥` scales as **`ka`** (slope 0.97 — 1st order,
  the genuinely-3-D component absent in any 2-D stack);
- the **threshold**: `ε = 1%` at `ka* ≈ 0.14`, i.e. **pitch/aperture `P/a ≈ 46`**.

**The bridge to §3.4 (rung-1).** The per-station 2-D twist design holds when the
**pitch exceeds the aperture by ~ several tens** — the *same* "longitudinal scale
≫ transverse scale by ~40×" rule as the straight magnet's foliate-and-perturb
(`L/gap ~ 40`), with the twist replacing `gap/L` by `a/P`. So the whole twist axis
closes: §3.6 (the twist + n-fold law), §3.7 (the combined-function confluence on a
curved orbit), §3.8 (its validity threshold) — a per-station 2-D design with a
*measured* fast-twist coupling, exactly mirroring the straight-magnet rung-1.

### 3.9 The END PACK in two planes — x-y cross-section + s-y end → 3-D

§3.5 applied the two-plane method to a *curved* FFAG sector **cell** (the whole
bend cell). The **end pack** — the magnet's longitudinal **termination**, the
genuinely-3-D hard part of a *straight* magnet — is the **same two-plane thought
on the literal `(x-y)` cross-section + `(s-y)` longitudinal planes**:

| plane | what it sets | engine |
|---|---|---|
| **`x-y` cross-section** (⊥ beam) | the transverse multipole: a finite flat pole droops (`b₃<0`), the **shim** `z=g/2−δ(x/w)²` zeroes it | `accel_pole_dipole_body_2d.solve` (Plane 1) |
| **`s-y` longitudinal** (beam `s`, gap `y`) | the **end chamfer**: a *standalone* 2-D Laplace fringe → the Rogowski bow-out `ĝ(s)` + the effective length `L_eff` | `endpack_two_plane.solve_sy_endpack` (Plane 2) |

The distinction from §3.2/§3.5: **both planes are cheap 2-D DESIGN solves done
FIRST**; the 3-D solve only **reflects + verifies** (it does *not* extract the end
profile out of itself). The reflection drives the 3-D pole as an **equipotential**
(`Ψ=mmf` on the pole, `Ψ=0` on the median plane — the high-μ limit, the same drive
as §3.5 rung-2, *pure Laplace, no coil*), then sweeps the chamfer **depth** to
drive the pole-tip corner field through its body value.

**Verified (`endpack_two_plane.py`, ngsolve only, golden-tested).**

- **Plane 1 (`x-y`):** the flat finite pole (half-width 60 mm) droops
  `b₃/b₁ ≈ −3.6×10⁻⁵`; a `δ ≈ 0.41 mm` concave shim zeroes it (residual
  `~1×10⁻⁴`).
- **Plane 2 (`s-y`):** a standalone 2-D Laplace fringe gives the Rogowski end
  bow-out `ĝ(s)` and `L_eff ≈ 151 mm` — a **+26 %** excess over the 120 mm iron
  (each end ~`0.75 g`, the same fringe law as §3.2/§3.5).
- **Reflection (3-D):** the hard-cut pole tip **over-fields by ~+11 %** (the corner
  flux concentration); reflecting the `ĝ(s)` shape and sweeping the depth drives
  that corner over-field **through zero at ~2.1 mm**, while the integrated
  transverse `b̄₃,₅` stays **~0.3 %** — body/Plane-1 dominated (the END shape is the
  wrong lever for it, the honest two-lever split of §3.2).
- **Cross-check:** the *cheap 2-D* `s-y` chamfer SHAPE predicts the *expensive 3-D*
  end equipotential bow-out to **~7 % rms** — the two-plane reflection is sound.

So the end-pack realizes the design thought literally: **design the END in `x-y`
(the transverse harmonic) and `s-y` (the longitudinal termination), each a cheap
2-D plane, then loft into one 3-D pole** — the transverse and longitudinal levers
cleanly separated (Plane-1 owns `b₃,₅`; Plane-2 owns the end taper / `L_eff`).

### 3.10 The SPECTROMETER end pack, NONLINEAR — the pole-tip corner is a saturable throat

A large **bending spectrometer** dipole runs near the iron knee, so the end pack
must be designed **with saturation**. The linear §3.9 result hands this to the
saturation framework directly: the hard-cut pole END concentrates flux at the tip
corner by `κ = tip_enhancement ≈ 1.11`, so in saturating iron the **corner reaches
the knee FIRST** — it saturates at a gap field

$$B_{gap}^{\text{corner-knee}} = B_K/\kappa \approx 1.5/1.13 \approx 1.33\ \text{T},$$

**~12 % BELOW the bulk iron knee** `B_K = 1.5 T`. Above that the corner `μ_r`
collapses, the flux can no longer follow the pole edge, and the **EFB (effective
field boundary `≈ L_eff`) drifts with excitation** — fatal for a spectrometer,
whose pole-edge **edge focusing** `tan β/ρ` depends on the EFB (the optics would
change with the field setting).

**The corner is exactly a Chaplygin saturable THROAT** (§3.3,
`clebsch_dipole_saturation_2d.py`): `κ` is its inverse cross-section, and it
saturates first. **The Rogowski end chamfer (§3.9's `s-y` design) is the
corner-throat width knob** — it lowers `κ`, raising the corner knee `B_K/κ` toward
the bulk knee. The **SAME chamfer that zeroes the linear corner over-field removes
the premature corner saturation**: linear (cosmetic field quality) and nonlinear
(avoid early saturation + EFB drift) levers POINT THE SAME WAY, and saturation gives
the chamfer its **hard engineering justification**.

**Verified (`endpack_spectrometer_saturation.py`, golden-tested, ngsolve only).**

- **The map at LINEAR cost.** Reuse the §3.9 equipotential corner concentration
  `κ(chamfer)` (the depth sweep) and overlay the Froehlich iron BH (`B_K=1.5 T`,
  `μ_r0=2000`, from `clebsch_dipole_saturation_3d.py`): the corner knee
  `B_K/κ(chamfer)` — flat `1.33 T` → `2.4 mm` chamfer `1.55 T` → `5 mm` `1.75 T`.
  To clear a `B_op = 1.45 T` operating field without corner saturation needs
  `κ ≤ 1.034`, a `≈ 1.4 mm` chamfer (the linear cosmetic optimum `κ=1` is `≈ 1.9 mm`
  — the same lever). The **whole nonlinear end-pack map = 4 linear equipotential
  solves + a BH overlay** — the Chaplygin "nonlinear analysis done linearly" applied
  to the END corner.
- **Design-grade, components validated.** The lumped `κ`-throat overlay is the same
  lumped-magnetic-circuit class as `clebsch_dipole_saturation_2d` (~10 % vs FEM), and
  its two ingredients are independently verified elsewhere: the corner concentration
  `κ` is the LINEAR equipotential `tip_enhancement` (§3.9, golden — geometry-only,
  drive-agnostic), and the Froehlich BH + the well-conditioned **B-input A-formulation**
  that backs the iron saturation are the §3.3 `clebsch_dipole_saturation_3d.py`
  (the documented cure: the reduced-Ω `μ(|H|)` Picard STALLS at high `μ`, the
  A-formulation does not). So the composition `B_corner = κ·B_gap` until `B_K` is
  well-founded.

**Honest scope.** A fully coil-driven 3-D corner-saturation FEM is the documented
expensive extension (the equipotential/MMF drive forces flux across the gap — a
*uniform* applied field does NOT reproduce the corner concentration `κ`, and the
coil's Biot-Savart `B_s` projection is the serial bottleneck), not run here.
The corner-`κ` softens before the hard knee (a real FEM is the truth). Curved/rotated
-EFB edge focusing (the horizontal `x-s` edge contour) and the fully-saturating sector
body remain the spectrometer extensions (§3.9 honest scope + the sector §3.5).

### 3.11 The two planes CO-BAKED into one pole — `z(x,s)=g/2−δ(x/w)²+lift(s)`

§3.9 named its own next refinement: its 3-D reflection carried the `s-y` chamfer at a
fixed body width, so the `x-y` shim `δ` was *verified as the transverse lever* but not
yet baked into the same 3-D pole. This **co-bakes both** into one gap face

$$z_{\text{face}}(x,s) = g/2 \;-\; \underbrace{\delta\,(x/w)^2}_{x\text{-}y\ \text{shim}} \;+\; \underbrace{\text{lift}(s)}_{s\text{-}y\ \text{chamfer}},$$

and shows BOTH levers act **at once**: the co-baked pole achieves a clean integrated
transverse `b̄₃,₅` AND a rounded pole-tip corner.

**Verified (`endpack_cobake.py`, ngsolve only, golden-tested).** The 4 cases (the
same equipotential-pole drive + integrated analyzer as §3.9) on the co-baked face:

| case | corner tip | transverse `b̄₃,₅` |
|---|---|---|
| baseline (flat cut) | `1.16` (over-fields) | `0.7 %` |
| shim only (`δ`) | `1.23` | **`0.07 %`** (x-y lever cleans it) |
| chamfer only (`ĝ`) | **`0.96`** (s-y lever rounds it) | `0.3 %` |
| **BOTH (co-baked)** | **`1.02`** | **`0.07 %`** |

i.e. **one pole face delivers a clean transverse harmonic AND a rounded corner** — the
two cleanly-separated two-plane levers composed in 3-D (`δ ≈ 0.41 mm` from §3.9 Plane 1;
`ĝ` the Rogowski shape from Plane 2).

**Honest scope of the staircase build (`endpack_cobake.py`).** The exact `δ(x/w)²`
shim needs an *x-varying* face, built there as an x-prism STAIRCASE (per-slab shim
offset), so the no-shim cases mesh coarser than the shim cases — the per-case *absolute*
numbers are research-grade, not precision. The locked claim is the **co-existence** of
both levers in the (well-resolved) BOTH pole; the per-lever *causation* is golden-locked
separately (the `x-y` shim zeroes `b₃`: `accel_pole_dipole_body_2d` / §3.9 Plane 1; the
`s-y` chamfer drives the corner over-field through `1`: §3.9's depth sweep).

**Precision construction (`endpack_cobake_loft.py`, ngsolve only, golden-tested).** The
clean construction the staircase pointed to is now built: the gap face is a SMOOTH OCC
`ThruSections` LOFT through per-x-station cross-section wires (each carrying its shim
offset `δ(xᵢ/w)²` + the chamfer `lift(s)`), so the surface is smooth in `x` (no facets).
The headline is **mesh consistency**: the smooth loft meshes the baseline (`δ=0`) and the
shim (`δ>0`) cases at the *same* density —

| build | `ne(shim)/ne(baseline)` |
|---|---|
| smooth LOFT (`endpack_cobake_loft.py`) | **`≈ 0.97`** (same density) |
| x-prism STAIRCASE (`endpack_cobake.py`) | `≈ 36` (merges `δ=0` slabs → coarse; steps `δ>0` → fine) |

so the loft **RESOLVES the staircase artifact** and the co-baked pole's `b̄₃,₅` + corner
become a *precision* claim. On that consistent mesh both levers still act: the chamfer
rounds the corner (`tip 0.99`), the shim removes the transverse content the chamfer
introduces (`both 0.47 % < chamfer-only 0.84 %`) and returns it to the baseline
mesh-noise floor (`~0.5 %`). *(Absolute `b̄₃,₅` differs from the staircase table above
because the meshes differ; the loft's value is the precision one.)*

**Rotated-EFB edge focusing — characterized NEGATIVE, not shipped (the field-EFB slope
is the wrong observable).** A first attempt read the `∫B` effective-field-boundary angle
out of the *equipotential* drive and it attenuated to `~0.47·β_cut`, recorded as an
unattributed open step. This session retried it the better way — a rigidly-rotated
*whole magnet* (iron **and** the CoilBuilder coil rotated together) driven by the genuine
reduced-Ω + Biot-Savart **source**, with the EFB read both as `∫B_z dy / B_z(body)` and as
the half-field crossing `y_half(x)` — and across all variants (parallelogram vs rigid
rotation, integral vs half-field EFB, narrow vs wide pole) the result is robust and
**negative**: at `β=0` the EFB slope is cleanly `≈ 0` (the method is unbiased), but for
`β>0` the field-EFB slope does **not** recover the geometric edge angle (it comes out
wrong-sign and many times `tan β`, and at larger `β` the per-line `B_z(body)` normalization
passes through zero). The attribution is now clear: **the per-beam-line field integral
`∫B_z dy` through a tilted finite magnet is not a local edge tracker** — the compact
fringe is fully 3-D and the field-EFB slope simply is not the edge angle (so the prior
`~0.47` was the same surrogate failing, not the drive). The genuine edge-focusing strength
is a **trajectory** quantity — the vertical kick `∫(v×B)` along particle orbits through the
fringe — so it needs **particle tracking**, not a field-EFB slope. Per the repository's
honest-results policy this is kept as a characterized open problem (recorded here), not
shipped as a working example.

### 3.12 The saturating sector body — the two planes respond OPPOSITELY to saturation

§3.5's scaling-FFAG **sector** is designed in two planes — the radial `(r,z)` field index
`k(r) ~ r^k` (the achromaticity) and the azimuthal `(s,z)` sector end (the effective length
`L_eff`). `scaling_ffag_pole_2d.py` showed the RADIAL plane droops `k(r)` at the high-r edge
under iron saturation (Step 2, the achromaticity wall); the azimuthal end was solved only
*linearly*. This rung drives the **sector body** into saturation by making the azimuthal
end solve NONLINEAR (Froehlich `μ_r(|B|)`, the same knee as Step 2) and solving it at the
high-r aperture edge (smallest gap, highest `B`) and the low-r body.

**Verified (`scaling_ffag_sector_saturation.py`, ngsolve only, golden-tested).** The honest
result for a scaling (large-gap) pole is a **contrast** between the two planes:

| plane | quantity | under saturation |
|---|---|---|
| azimuthal `(s,z)` | effective length `L_eff` | **ROBUST** — drift `< 0.1 %` even where the high-r iron `⟨μ_r⟩` collapses `×0.3` |
| radial `(r,z)` | field index `k(r)` | **FRAGILE** — droops `Δk ≈ −0.26` at the high-r edge |

The sector END is **gap-reluctance-dominated**: the high-r iron saturates hardest (its
`⟨μ_r⟩` collapses `1574 → 481`, `×0.3`), yet `L_eff` barely moves (the fringe stays `~1.5`
gaps) — the same honest scope as §3.6 / `clebsch_dipole_saturation_3d` (a large-gap magnet's
gap field softens only mildly with iron saturation). So **saturation degrades the radial
field SHAPE (achromaticity) but not the azimuthal end LENGTH** — a real design insight: the
high-r achromaticity needs the radial reshape (§3.5 Step 3), while the sector ENDS are
saturation-robust and need no nonlinear end correction.

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

- φ-side iron pole — the hodograph (`docs/clebsch_hodograph/demos/`).
- A-side coil current — the stream function
  (`src/radia/stream_function.py`, `validation_test/feec/vim_legacy/foliated_solenoid_wires.py`,
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
  lands only for long magnets (~10 % fringe at L/gap ~ 40);
- the **two-plane → 3-D method** (§3.5, `ffag_sector_two_plane.py`): the FFAG
  scaling sector designed in two orthogonal 2-D planes — transverse `(r,z)`
  scaling index (Plane A) + azimuthal `(s,z)` sector ends (Plane B) — with the
  reflection's validity the `L/g` leaf-coupling: each sector end adds ~`0.75 g`,
  `L_eff = L_sector + 1.5 g`, fringe ~ gap/L (slope −0.98); and **rung 2** (the
  3-D reflection, `--rung2`) revolves `g(r)` into a 3-D iron pole and **recovers
  the field index `B(r) ∝ r^k` in full 3-D** (mean `≈ 4.88`, design 5) with the
  azimuthal end fringe (`+39%` at `L/g=3`);
- the **beam-referenced equipotential surface as the design primitive + the
  twist** (§3.6, `twisting_quadrupole_pole.py`): the quad pole = the hyperbola
  equipotential; rotating the surface by `φ` rotates the recovered orientation by
  exactly `φ` (**slope 1.000, err 0.00°**) — the n-fold law (surface twist `φ` ⟷
  multipole phase `2φ`), with a clean quad (skew `~5e-6`, `b_6 ~ 5.6e-3`
  rotation-invariant);
- the **combined-function magnet on its curved orbit** (§3.7,
  `combined_function_frenet_sweep.py`): a tilted-gap dipole+quad cross-section,
  fixed in the Frenet frame, **twists by `θ(s)` in the lab** — both the dipole and
  quad orientations track `θ` (**slope 1.000, err 0.00°**) with the multipole
  phase-change ratio `ψ_2/ψ_1 = 2.000` (the n-fold law) — the confluence of §3.5
  (the dipole sector) and §3.6 (the twist);
- the **fast-twist leaf coupling** (§3.8, `twist_rate_leaf_coupling.py`): the
  per-station 2-D twist design breaks when the twist is fast — the exact helical
  multipole deviates from the 2-D stack as `ε ~ (ka)²` (transverse, slope 2.04)
  with a longitudinal `B_s ~ ka` (slope 0.97); the threshold `ε = 1%` is at
  **pitch/aperture `~46`**, the twist analogue of rung-1's `L/gap ~ 40`
  (`gap/L → a/P`).
- the **end pack in two planes** (§3.9, `endpack_two_plane.py`): the magnet END
  designed in the `x-y` cross-section (the shim `δ≈0.41 mm` zeroes `b₃≈−3.6e-5`)
  **and** the `s-y` longitudinal plane (a standalone 2-D Laplace fringe →
  Rogowski `ĝ(s)`, `L_eff` **+26 %**), reflected into one 3-D pole
  (equipotential-pole drive, pure Laplace): the chamfer-depth sweep drives the
  pole-tip corner over-field (`+11 %` hard-cut) **through zero at ~2.1 mm**, the
  integrated `b̄₃,₅` stays `~0.3 %` (body lever), and the **cheap 2-D `s-y` chamfer
  shape predicts the 3-D end equipotential to ~7 % rms**.
- the **spectrometer end pack, NONLINEAR** (§3.10,
  `endpack_spectrometer_saturation.py`): the pole-tip corner is a Chaplygin saturable
  THROAT — `κ ≈ 1.13` means it reaches the iron knee FIRST (corner saturates at
  `B_K/κ ≈ 1.33 T`, **~12 % below the bulk knee** `B_K=1.5 T`), drifting the EFB /
  edge focusing; the Rogowski chamfer is the throat-width knob (`κ↓` ⇒ corner knee
  `B_K/κ↑`), the SAME lever as the linear over-field zero, now with a hard saturation
  justification. The whole nonlinear map = the linear equipotential `κ(chamfer)` sweep
  + a BH overlay (the §3.3 Chaplygin "nonlinear-done-linearly" on the END corner) —
  design-grade (the lumped-circuit class of `clebsch_dipole_saturation_2d`), with `κ`
  from the §3.9 linear equipotential (golden) and the BH + A-formulation from
  `clebsch_dipole_saturation_3d` (committed) as its independently-validated components.
- the **two planes CO-BAKED into one pole** (§3.11, `endpack_cobake.py` staircase +
  `endpack_cobake_loft.py` precision LOFT): both the `x-y` shim `δ` and the `s-y`
  Rogowski chamfer `ĝ` baked into one gap face `z(x,s)=g/2−δ(x/w)²+lift(s)` — the
  co-baked pole delivers a clean integrated transverse `b̄₃,₅` AND a rounded pole-tip
  corner at once, the §3.9 two-plane levers composed in 3-D. The PRECISION construction
  (`endpack_cobake_loft.py`, OCC `ThruSections`) builds the gap face as a SMOOTH loft so
  the baseline + shim cases mesh at the SAME density (`ne(shim)/ne(baseline) ≈ 0.97`, vs
  the staircase's `≈ 36`) — it RESOLVES the documented staircase artifact, making the
  co-baked `b̄₃,₅` + corner a precision claim. The rotated-EFB **edge-focusing** extension
  remains **not shipped** — retried this session with a whole-magnet rotation + reduced-Ω
  source, the field-EFB slope does NOT recover the geometric edge angle (it is the wrong
  observable; the genuine focusing is a particle-tracking quantity), kept as a characterized
  negative (§3.11).
- the **saturating sector body** (§3.12, `scaling_ffag_sector_saturation.py`): the
  scaling-FFAG sector's azimuthal `(s,z)` end made NONLINEAR (Froehlich `μ_r(|B|)`) — the
  two sector planes respond OPPOSITELY to iron saturation: the azimuthal effective length
  `L_eff` is ROBUST (gap-reluctance-dominated; drift `< 0.1 %` even where the high-r iron
  `⟨μ_r⟩` collapses `×0.3`), while the radial field index `k(r)` is FRAGILE (`Δk ≈ −0.26`,
  the §3.5 achromaticity wall). Saturation degrades the field SHAPE, not the end LENGTH.

**Research program (named, not claimed done):**
- the end-design loop is **closed in two planes** (§3.9): the longitudinal
  end-field (the `s-y` Rogowski chamfer, depth tuned to zero the corner over-field)
  AND the transverse harmonic (the `x-y` shim that zeroes `b₃`) — cleanly
  separated levers; what remains is the FULL co-baked 3-D loft carrying BOTH the
  `δ`-shim curvature and the `ĝ(s)` chamfer in one pole surface (here the 3-D
  reflection carries the `s-y` chamfer at fixed body width; the `x-y` shim is
  verified as the transverse lever, baked-together loft is the next refinement);
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
| φ-side iron pole (hodograph) | `docs/clebsch_hodograph/demos/` |
| A-method dual | `docs/clebsch_hodograph/demos/a_method_clebsch_2d.py` |
| exact open boundary | `docs/clebsch_hodograph/demos/hodograph_kelvin_axisym.py` |
| cohomology (current-linking) | `src/radia/cohomology.py` |
| pole geometry + multipole analyzer | `docs/clebsch_hodograph/demos/accel_pole_design.py` |
| 3-D ends: integrated analyzer + end rule | `docs/clebsch_hodograph/demos/accel_pole_ends_3d.py` |
| 3-D ends: FEM rung (reduced-Ω + CoilBuilder) | `docs/clebsch_hodograph/demos/accel_pole_ends_fem.py` |
| forward (reduced potential + CoilBuilder) | `src/radia/panels/calc_accel_magnet.py` |
| radial field index (scaling + isochronous, saturation) | `docs/clebsch_hodograph/demos/scaling_ffag_pole_2d.py` |
| foliate-and-perturb scaling (leaf coupling ~ gap/L) | `docs/clebsch_hodograph/demos/leaf_coupling_perturbation_3d.py` |
| two-plane → 3-D method (FFAG sector: transverse + azimuthal) | `docs/clebsch_hodograph/demos/ffag_sector_two_plane.py` |
| beam-referenced equipotential surface + the twist (n-fold law) | `docs/clebsch_hodograph/demos/twisting_quadrupole_pole.py` |
| combined-function on a curved orbit (the Frenet sweep = twist) | `docs/clebsch_hodograph/demos/combined_function_frenet_sweep.py` |
| fast-twist leaf coupling (the per-station 2-D validity threshold) | `docs/clebsch_hodograph/demos/twist_rate_leaf_coupling.py` |
| end pack in two planes (x-y cross-section + s-y end → 3-D) | `docs/clebsch_hodograph/demos/endpack_two_plane.py` |
| spectrometer end pack NONLINEAR (corner saturable throat) | `docs/clebsch_hodograph/demos/endpack_spectrometer_saturation.py` |
| two planes co-baked into one pole (δ shim + ĝ chamfer) | `docs/clebsch_hodograph/demos/endpack_cobake.py` |
| co-bake as a PRECISION tensor loft (OCC ThruSections) | `docs/clebsch_hodograph/demos/endpack_cobake_loft.py` |
| saturating sector body (azimuthal L_eff robust, radial k fragile) | `docs/clebsch_hodograph/demos/scaling_ffag_sector_saturation.py` |
| A-side coil (stream function) | `src/radia/stream_function.py`, `validation_test/feec/vim_legacy/foliated_solenoid_wires.py` |
