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
  analyzer. Verified: a flat-top dipole `B_z(body) ≈ 0.15 T`, effective length
  `L_eff ≈ 165 mm` (iron 120 mm + ~23 mm fringe each end), a **+11 % pole-END
  enhancement** (the flux concentration at the iron ends — the design-relevant
  end effect), integrated dipole `b̄₁ ≈ 0.026 T·m`, integrated spurious
  (n=3,5) ≈ 9 %. Figure: the on-axis `B_z(y)` (flat body + end bumps + fringe).
  Two engineering notes baked in: **RadiaField (Biot-Savart) is not thread-safe
  under TaskManager** — the source LinearForm + field readout assemble serially,
  only the stiffness+solve are wrapped; and the CoilBuilder racetrack arcs curve
  toward −x (start at +R to centre it). **Next rung:** extract the solved 3-D
  equipotential surface `Ψ = Ψ_pole` at the end as the end-iron contour, then
  re-shape → re-solve to drive the integrated spurious down.
- **`one_turn_coil_streamfunction.py`** — the A-side **(B) track**: a 1-turn
  coil is the coarsest stream-function discretization (one contour = one wire);
  the task is the single best wire path, and the script shows the honest 1-turn
  limit vs the full multi-turn stream-function current.

## Run

```bash
python a_method_clebsch_2d.py                 # A-method net figure
python hodograph_kelvin_axisym.py             # Kelvin exact-open-boundary flux figure
python cohomology_hodograph_currentlink.py    # when the hodograph needs cohomology
python accel_pole_design.py                   # multipole analyzer + quad pole geometry
python accel_pole_harmonics.py                # 2-D equipotential lever (shim → harmonics)
python accel_pole_ends_3d.py                  # 3-D ends: integrated analyzer + end rule
python accel_pole_ends_fem.py                 # FEM rung: reduced-Omega + CoilBuilder dipole
python one_turn_coil_streamfunction.py        # (B) the 1-turn stream-function limit
```

Locked by `tests/feec/test_clebsch_hodograph_research.py` (8 tests; the FEM rung
is `@pytest.mark.slow`).

## Prior art (honest)

The 2-D bidirectional "potential coordinate" map is published — Dervisha,
Marjamäki, Rasilo, Tarhasaari, *"Bidirectional Coordinate Transformation and
Its Application to 2-D Magnetic Field Problems"*, CEFC 2026 (Bossavit
school). The 2-D complex-potential / conformal pole design lineage is older
still (Rogowski 1923, Halbach). These scripts are an independent NGSolve
realisation coupled to the lab's reduced-potential + Kelvin machinery — **not
a novelty claim.** See `memory/derham_field_current_duality.md` for the full
prior-art verdict.
