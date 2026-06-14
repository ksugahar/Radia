# 3D Stream Function + Cohomology — Clebsch / vector-T coil design on the de Rham complex

Research demonstration of a **volume-source 3D stream-function coil-design method**, and an
honest map of **where it stops working**. Everything here — the stream function, the Clebsch
potentials, the cohomology cuts, the single-stroke deformation — lives on **one NGSolve de Rham
(FEEC) complex**:

```
H1  --grad-->  H(curl)  --curl-->  H(div)  --div-->  L2
0-form          1-form             2-form          3-form
```

| object | space | role |
|---|---|---|
| stream / Clebsch scalars `ψ, χ`; scalar potential `φ`; 2D flux function `A` | **H1** | the potentials you design with |
| cohomology generators `h_k`; `grad` of a multivalued `χ`; the `H` field | **H(curl)** | H¹ cohomology = curl-free-not-gradient (net current / cuts) |
| `B`, and `J = grad(ψ)×grad(χ) = grad(λ)×grad(μ) = curl T` | **H(div)** | H² cohomology = div-free-not-curl (flux / cavities) |
| single-stroke wire deformation field | **VectorH1** | shape derivative for manufacturing |

Because all of these are native NGSolve/Netgen FEEC objects, the whole stack is **Gmsh-free**
(cohomology generators are the Hodge-Laplacian zero-eigenspace on H(curl); no `computeHomology`).

> **Origin (CEFC 2026, WA-O1, Thessaloniki).** Talk #2 (Tampere — Dervisha, Marjamäki, Rasilo,
> **Tarhasaari**) = the exterior-calculus *bidirectional potential-coordinate map* (`bidirectional_map_2d.py`).
> Talks #1 & #4 (Zhou/Yan/Xu/**Ren**, IEE-CAS + GeePs) = the *cohomology dual / cut* formulation
> (the multivalued-`χ` / net-current layer). Cohomology-cut citation stack: **Kotiuga** (theory)
> → **Pellikka 2013** (the algorithm `cohomology_cut.py` wraps) → Ren 2002 (T-Ω application).

## The method in one paragraph

A divergence-free current `J` (DC continuity) is the design unknown. There are **three genuine
"volume stream function" forms**, all native to the complex:

| form | representation | inverse design | wires? | helicity |
|---|---|---|---|---|
| **vector-T** | `J = curl T`, `T ∈ H(curl)` | **CONVEX** (linear in `T`; ACA+TSVD, gauge auto-truncated) | not directly (a current *distribution*) | no limit (any `J`) |
| **Clebsch, `μ` free** | `J = grad(λ)×grad(μ)` | non-convex (bilinear) | yes — `{λ=c}∩{μ=c}` level-set intersections | forces `H = 0` |
| **Clebsch, `μ` fixed (foliated)** | `J = grad(λ)×grad(μ)`, `μ` = a chosen foliation | **CONVEX** (linear in `λ`; existing `radia.stream_function` ACA+TSVD **unchanged**) | yes — `λ`-contours per `μ`-leaf | forces `H = 0` |

The surface stream-function method (`K = n̂×grad(ψ)`, `radia.stream_function`) is the
**single-surface special case** of the foliated form (`μ` = one winding surface); the volume
version stacks it into a true 3D bulk winding. **The whole story is verified by `sympy`** (helicity
identity, the unified Jacobian, gauge null space) before any FEM is run.

## Stage A — the convex baselines (it works)

| # | script | does | validated against | numbers |
|---|---|---|---|---|
| 1 | `bidirectional_map_2d.py` | Tampere potential-coordinate map `(A,φ)`; forward geometry→potentials, inverse potentials→geometry (`x,y` harmonic in `(A,φ)` at `μ=1`) | **ANALYTIC** `w=ζ²` conjugate pair | forward `A,φ` `1e-15`; orthogonality `grad A·grad φ = 2.6e-14`; **Jacobian `det = |B||H|` to `1.6e-16`**; inverse `1.6e-8`; round-trip `6e-8` |
| 1b | `bidirectional_map_axisym.py` | **axisymmetric** (helicity-0) map `(ψ=rA_φ, Φ)`; forward GS + axisym-Laplace, inverse pointwise Newton. **Finding**: the `1/r` Stokes metric breaks conformality → map is only *r-conformal* (`\|∇ψ\|=r\|∇Φ\|`), so the inverse is **not** a plain harmonic image solve (`z(ψ,Φ)` is cubic, non-harmonic) | **ANALYTIC** `l=2` zonal gradient field `ψ=r²z, Φ=z²−r²/2` | forward `ψ,Φ` `1e-15`; axisym CR conj. `1.2e-14`; **Jacobian `det = r\|∇Φ\|² = r\|B\|²` to `9.8e-15`**; round-trip `7e-16` |
| 2 | `vector_t_inverse.py` | `J = curl T` (HCurl) **convex** inverse via ACA+TSVD | self-consistency + **Radia** recompute | fit `6e-16`; `div J = 5.7e-16` (curl T exact); **gauge `T=grad χ → field 2.2e-20`** (TSVD truncates; `k_aca=9 ≤ M ≪ ndof=3360`) |
| 3 | `foliated_clebsch_solenoid.py` | foliated `J=grad(λ)×grad(μ)`, `μ=r` fixed → linear in `λ` → existing ACA+TSVD **unchanged** | **ANALYTIC** uniform `Bz` + weak div | uniform `Bz` to `3.3e-5`; independent Biot-Savart recompute matches; weak `div J = 8.5e-6` |

**The key reframe (verified):** the volume inverse design is **not** fundamentally non-convex — that
is *specific to the two-scalar Clebsch* parameterisation. The full vector `T` gives a **convex**
least-norm (`#2`), with the gauge null space (`curl grad χ = 0`) truncated by TSVD automatically —
**no tree-cotree gauge needed**.

## Stage B — wire extraction (distribution → windable wires)

| # | script | does | validated against | numbers |
|---|---|---|---|---|
| B1 | `helicity_diagnostic.py` | `H_rel = |∫T·curl T| / (‖T‖‖curl T‖) ∈ [0,1]` — gates whether a clean Clebsch/level-set extraction is even possible | **ANALYTIC** axial vs Beltrami | axial (Clebsch-type) `9e-5`; ABC/Beltrami (`curl T = kT`) `1.000` |
| B2 | `foliated_solenoid_wires.py` | foliated `λ` → per-cylinder equal-Δλ contours → equal-current wires → Biot-Savart | target `Bz` + **two-codebase** (Radia `ObjFlmCur`+`Fld`) | 59 equal-current wires (`I=Δλ·Δμ`); `Bz` to `4.7%` (discretisation); **agrees with Radia to `3.4e-10`** |
| B3 | `line_topology_gate.py` | the **second** extraction gate — do the lines *close*? rotation-number closure-defect of a torus winding (robust closure metric, **not** a fragile chaotic-Poincaré 2-D-fill) | **ANALYTIC** rational vs irrational `ι` | `ι=1/5, 2/5` **close at turn 5** (defect `0.068` minor-radii → windable); `ι=1/φ` **never closes** (`0.572`); **8.4× margin** |

`H_rel ~ 0` for `#2`/`#3` (so Stage A is extractable); the two-codebase invariant (`3.4e-10`) means
the golden cannot pass on a wrong Biot-Savart formula.

## Frontier — where it stops working (honest walls, golden-locked)

Each frontier was **adversarially attacked** (a genuine attempt to crack it) before conceding, then
demonstrated with a runnable probe. **All three verdicts = `fundamental_wall`** — honest
partly-negative results (repo-first: a clearly-demonstrated wall is a real result).

| # | script | the wall | honest remedy |
|---|---|---|---|
| F1 | `foliation_choice_wall.py` | choosing the foliation `μ` for a general target = **bilevel convex-inner / non-convex-outer**. Field is fit to `~1e-15` at every tilt → the non-convexity lives in **coil complexity** `‖J‖` (the REGCOIL measure), *invisible to the field residual*. single-axis target = 1 min (convex); two-axis = **2 minima + 18% barrier** (mesh-stable). = the stellarator winding-surface problem (Merkel/Landreman/**Kaptanoglu arXiv:2408.08267**). | **REMEDY IMPLEMENTED**: `cmaes_foliation_gauge.py` — CMA-ES (Optuna `CmaEsSampler`) over a 2-param foliation normal; finds the global basin (4-basin landscape) where a local optimizer sticks at **1.7×** higher `‖J‖`. **or** accept vector-T |
| F2 | `clebsch_recovery_wall.py` | recovering `(λ,μ)` **from** a general `J` is **ill-posed**: WALL1 non-unique (transport `J·grad μ = 0`: seeds `μ₀=x` vs `x²` both valid, `μ_b=μ_a²` → infinite-dim **area-preserving-diffeo gauge orbit**); WALL2 conditioning `~1/|J|²`, log-log slope **−4.00** at field nulls (**Yoshida 2009**, Qin 2018) | local, seed-dependent recovery only |
| F3 | `helical_current_no_clebsch.py` | **nonzero helicity = topological obstruction**. ABC/Beltrami target `H_rel~1`; any convex Clebsch fit is FORCED to `H_rel~0` → the **helicity gap `~1.000`** (foliation- AND mesh-independent; the L2 residual is the *wrong* metric). Clebsch streamline closes (`2e-3`), ABC does **not** (`143` ring-lengths, ergodic; **Moffatt 1969**, Enciso–Peralta-Salas 2020) | vector-T bulk distribution **or** a multi-patch Clebsch atlas + cohomology cuts |

## Domain of validity (the bonsai shape)

The method works **cleanly** exactly in the **zero-helicity, fixed-foliation, closed-streamline**
regime (Stage A/B). F1/F2/F3 mark precisely where it stops being clean — not bugs, but the true
shape of the method. There is **no fake "general 3D solver"** here, by design.

```
zero helicity + closed lines + chosen foliation   →   clean equal-current wires   (Stage A/B ✅)
general foliation / arbitrary target               →   non-convex outer loop        (F1 wall → CMA-ES ✅)
recover Clebsch from an arbitrary J                 →   ill-posed (non-unique/null)  (F2 wall)
nonzero helicity (linked/knotted lines)            →   no clean wires at all        (F3 wall)
```

## Net-current (cohomology) + the F1 remedy — now implemented

The two remaining "finish" pieces of the coil designer are shipped + golden-locked:

| piece | what | files |
|---|---|---|
| **net current (cohomology)** | a single-valued `ψ`/`λ` carries zero net current around a winding-surface cycle; the secular term `Σ c_k h_k` (H¹, dim = b1) supplies it, INTEGRATED into the design solve | `../stream_function/cohomology_net_current.py` (+ golden) |
| **harmonic generators from scratch** | the canonical (min-norm) generators `h_k` computed via the **matched Whitney complex** (`H1₁ →grad→ HCurl₀`, `curl∘grad=0` to `1.4e-14`) — no analytic ansatz; 2 forms drop out of the Hodge-Laplacian kernel (gap `1.2e13`), class-matched to the net current (`2.3e-3`) | `../stream_function/cohomology_generators_whitney.py` (+ golden) |
| **F1 gauge choice (CMA-ES)** | global optimization over the foliation gauge escapes the non-convex wall; a local optimizer sticks at `1.7×` higher `‖J‖` | `cmaes_foliation_gauge.py` (+ golden) |
| **single-stroke (rigorous)** | the manufacturing single-stroke = the **level-set of `Ψ = ψ + secular`** (one connected helix, vs N disconnected rings); the connection IS the cohomology secular term (net axial `I`), **distributed** axisymmetrically — vs the heuristic **localised** rung (near-join stray `2.94×` less). Puts the `single-stroke-chain` skill's heuristic on a rigorous cohomology footing. | `single_stroke_clebsch.py` (+ golden) |
| **multi-layer single-stroke** | multi-layer SF = the **foliated Clebsch `J=∇λ×∇μ`** (μ=layer); the single-stroke = a `(λ,μ)` path (within-layer λ-spiral + layer-step μ) = **boustrophedon** (alternating axial sense). One wire, valid M-layer solenoid (`Bz≈M×`), and the alternating μ-steps **cancel** the per-layer axial-lead secular currents → `~121×` less far-field stray than naively stacking same-handed layers. | `multilayer_single_stroke_clebsch.py` (+ golden) |

So the cohomology-and-gauge "finish" of the stream-function coil designer is complete; F2/F3 remain genuine frontier walls (by design).  The single-stroke is the **clean solenoid case** of the rigorous construction; arbitrary-pattern routing is the F1-like open part (rigorous parameterisation, non-convex routing).

## Run + golden tests

```bash
# any single example prints its diagnostics and a PASS line
python examples/vim/bidirectional_map_2d.py
python examples/vim/vector_t_inverse.py
python examples/vim/foliated_clebsch_solenoid.py
python examples/vim/helicity_diagnostic.py
python examples/vim/foliated_solenoid_wires.py
python examples/vim/foliation_choice_wall.py
python examples/vim/clebsch_recovery_wall.py
python examples/vim/helical_current_no_clebsch.py

# the 38 golden tests (Stage A 14 + Stage B 9 + frontier 15)
python -m pytest tests/feec/test_bidirectional_map_2d.py tests/feec/test_vector_t_inverse.py \
  tests/feec/test_foliated_clebsch_solenoid.py tests/feec/test_helicity_diagnostic.py \
  tests/feec/test_foliated_solenoid_wires.py tests/feec/test_foliation_choice_wall.py \
  tests/feec/test_clebsch_recovery_wall.py tests/feec/test_helical_current_no_clebsch.py
```

All examples wrap NGSolve work in `with TaskManager():` (Caller-Wraps policy); each `main()`
returns `(dict, ok)` so the goldens import it directly. Reuses `src/radia/stream_function.py`
(`aca_tsvd`, `pseudo_inverse_solve`) and Radia (`ObjFlmCur`+`Fld`) ground truth — **no new solver
code**, only new basis-current callbacks on the existing engine.

## Shipped as a panel (examples -> panels)

Stage A/B (the foliated-Clebsch volume designer) is promoted to a Layer-3 panel:
the **`radia_streamfunction` "Volume 3D" mode** -- a hollow-cylinder conductor
`.vol` + a target axial `Bz` -> equal-current windable wires + a GMSH wire overlay.

| tier | file |
|---|---|
| shipped pipeline | `src/radia/streamfunction_volume.py` (`design_volume_coil`) |
| headless calc (Stage 2) | `src/radia/panels/calc_streamfunction_volume.py` |
| PySide6 mode (Stage 3) | `src/radia/radia_streamfunction.py` (`_Volume3DPanel`) |
| golden lock | `tests/panels/test_streamfunction_volume_golden.py` (+ frozen `.vol` fixture) |

The golden reproduces this example's numbers (n_wires=59, field residual 4.7%,
two-codebase Biot-Savart `3.4e-10`).  The panel mode covers the **clean regime
only** (Stage A/B); the F1/F2/F3 frontier walls stay research demos, not panel
modes, by design.

## Detailed home

- MCP knowledge: `streamfunction("clebsch_3d")` (and `cohomology` / `fusion` for the surface
  net-current / `b1=2` case already shipped).
- `memory/clebsch_cohomology_streamfunction_unification.md` — the full record, verified-math
  provenance, and the honest-frontier delineation.
- Related: `src/radia/clebsch_potential.py` (the axisymmetric `ClebschSolver`),
  `src/radia/cohomology_cut.py` (gmsh-free cohomology cut via `radia.cohomology`, the Ren/Pellikka layer).
