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

## Run

```bash
python a_method_clebsch_2d.py                 # A-method net figure
python hodograph_kelvin_axisym.py             # Kelvin exact-open-boundary flux figure
python cohomology_hodograph_currentlink.py    # when the hodograph needs cohomology
```

Locked by `tests/feec/test_clebsch_hodograph_research.py` (3 tests).

## Prior art (honest)

The 2-D bidirectional "potential coordinate" map is published — Dervisha,
Marjamäki, Rasilo, Tarhasaari, *"Bidirectional Coordinate Transformation and
Its Application to 2-D Magnetic Field Problems"*, CEFC 2026 (Bossavit
school). The 2-D complex-potential / conformal pole design lineage is older
still (Rogowski 1923, Halbach). These scripts are an independent NGSolve
realisation coupled to the lab's reduced-potential + Kelvin machinery — **not
a novelty claim.** See `memory/derham_field_current_duality.md` for the full
prior-art verdict.
