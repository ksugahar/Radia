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
designed multipole.* **Next rung:**
replace the analytic Maxwellian field with the **reduced-potential + CoilBuilder
FEM** end and read the solved 3-D equipotential surface as the end-iron contour
(then close the loop: re-shape, re-solve, converge the integrated harmonics).

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
  (`src/radia/stream_function.py`, `examples/feec_vim/foliated_solenoid_wires.py`,
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
- the stream-function coil (A-side), vs Radia to 3.4e-10;
- the reduced-potential + CoilBuilder forward engine (the panel).

**Research program (named, not claimed done):**
- the self-consistent *iron = equipotential* end-design loop (§3.2) — the FEM
  rung: read the solved 3-D equipotential as the end-iron contour and converge;
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
| forward (reduced potential + CoilBuilder) | `src/radia/panels/calc_accel_magnet.py` |
| A-side coil (stream function) | `src/radia/stream_function.py`, `examples/feec_vim/foliated_solenoid_wires.py` |
