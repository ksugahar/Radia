# Radia Open-Boundary Methods — A Unified Map

Radia closes an unbounded exterior several ways — **Kelvin** transformation,
**BEM** coupling, **PML**, **IABC** (improvised absorbing shells), the
**asymptotic Robin** BC, and the **CLN** (Cauer-ladder) open boundary. They look
unrelated, but they are all discretizations or reductions of **one** object: the
exterior Dirichlet-to-Neumann operator. This document maps them onto the few axes
that actually select one in practice, with measured anchors reproduced in-repo.
See **The selector** below for the radia policy (which method when; **IABC is
retired** as a radia method — kept only as the comparison record). The **measured,
quantitative** companion to this map is
[`DTN_SPECTRUM_COMPARISON.md`](DTN_SPECTRUM_COMPARISON.md) — every closure (Kelvin / BEM /
PML / CFS-PML / ballooning / Robin) on **one yardstick, the per-multipole DtN-spectral
defect**, across the static / eddy / high-freq regimes (reproduced by `act7_22`).

It is the cross-cutting companion to the two single-method references:
[`docs/kelvin/DTN_SPECTRUM_COARSE_MESH.md`](../kelvin/DTN_SPECTRUM_COARSE_MESH.md)
(Kelvin, the **spatial** side) and
[`docs/cln/CAUER_LADDER_NETWORK.md`](../cln/CAUER_LADDER_NETWORK.md) (CLN, the
**temporal** side).

---

## 0. The selector (radia / Sugahara-lab open-boundary policy, refined 2026-06-21)

**LAB POLICY: the lab closes the exterior with the EXACT-OPERATOR → CLN route, NOT PML.**
PML / CFS-PML are kept only as the benchmark **foil** (measured on the DtN yardstick —
[`DTN_SPECTRUM_COMPARISON.md`](DTN_SPECTRUM_COMPARISON.md), `act7_22` + the real NGSolve FEM+PML
`act7_24`), never as a solving method. Pick the method by **what the boundary operator IS** — and
note the two axes are different: **Kelvin = the spatial boundary; CLN = the temporal `s` axis**:

| boundary operator | regime | method |
|---|---|---|
| **static Laplace, `ω`-independent** — air truncation (the usual outer boundary) | non-radiating | **Kelvin ALONE** (`−(n+1)/R`; no `s` ⇒ **no CLN job**) |
| **`√s`-dependent** — conductor / diffusive SIBC (a semi-infinite conductor interior) | evanescent | **multipole-Zs → CLN** (analytic `√s`; **Kelvin-built → CLN only for a non-separable conductor**) |
| **complex, `s`-dependent** — Helmholtz wave DtN | radiation | **Grote-Keller → CLN** (exact rational NRBC = pole/Cauer ladder; prior art, not PML) |

- **CLN earns its keep ONLY where the boundary operator is `s`-dependent** — a conductor boundary
  or radiation. The air open boundary is static ⇒ **Kelvin alone**; "Kelvin → CLN" *combined* is
  essentially never needed in the non-radiating regime (only an exotic non-separable conductor).
  Kelvin (spatial) and CLN (temporal/`s`) are **different boundaries / different axes** — do not
  conflate them.
- **Honest scope:** principled for the lab's MQS-magnetics class (compact / quasi-spherical).
  Exact for separable geometry; an arbitrary non-separable boundary ⇒ build (Kelvin-FEM / BEM)
  then CLN band-fit. PML's generality (antennas / scattering / arbitrary radiating geometry) is for
  problem classes outside the lab — "no PML" is domain-appropriate, **not** "PML universally bad."
- **Grote-Keller → CLN is prior art** (Grote-Keller 1995/96/98; Hagstrom-Warburton; Guddati; Birk;
  Cauer-Foster; Kameari-Sugahara 2018) — adopted as a repository capability, not a novelty claim.

The homogeneous-exterior-vs-material-body refinement below is *within* the exact-operator route
(read "transient → CLN" as the **conducting / diffusive** case — a static air boundary stays
Kelvin-alone even in a transient problem):

Radia is **Laplace-kernel / MQS–Darwin** (no full-wave). Within that scope the
open boundary is picked on two axes — **is the exterior homogeneous, or does it
contain a material body?** (does the DtN come in closed form, or must Kelvin-FEM
build it?) × **do you solve once, or evolve in time / couple to a circuit?** (do
you need CLN temporal reduction on top?):

| | **solve once** (static / single shot) | **transient / repeated / circuit-coupled** |
|---|---|---|
| **homogeneous exterior** (free space / layered) | **Kelvin** — build the DtN, solve once (closed-form `−(n+1)/R`) | **exact DtN → CLN** — `√s`-Cauer, exact at `n+1` stages, every multipole `n=1..6` (`act7_20`: beats the impedance-shell route on all 4 axes) |
| **exterior contains a body** (iron shield, 2nd magnet) | **Kelvin** — mesh the exterior body, solve once | **Kelvin builds a material-aware DtN → CLN reduces it** (Track-B SF-with-iron) |

- **Exterior body ⇒ Kelvin** is *the* reason Kelvin exists: the closed-form DtN
  symbol assumes a homogeneous / separable exterior, so a material body in the
  exterior breaks it, and only Kelvin (a genuine compactified exterior FEM) carries
  it. If the exterior is homogeneous you usually prefer the closed-form DtN.
- **There is no "high-freq → PML" branch inside radia.** The MQS eddy exterior is
  *evanescent for every ω* (`√s` diffusion — no propagating regime), so raising ω
  stays in the Kelvin / DtN-CLN home (`act6_09`: CLN beats even CFS-PML there).
  **PML is only for genuine wave radiation (Helmholtz) = outside radia's
  Laplace-kernel scope → that is NGSolve's job, not radia's.** "I need PML" ⇒ "I
  have left radia's scope."
- **IABC (absorbing impedance shells) is RETIRED as a radia method** (2026-06-20):
  the impedance boundary is realised better as *exact DtN → CLN* (top-right cell),
  and IABC's only niche — high-freq radiation (HOIBC) — is outside radia's scope.
  There is **no `iabc()` MCP tool**; the exact-impedance / `Zs` → DtN → CLN
  knowledge moved to `dtn_coarse_mesh(topic="dtn_to_cln")`. The IABC investigation
  is kept as the **negative-result / comparison** corpus only: the `act7_*` demos +
  `dtn_coarse_mesh(topic="method_map")`. (The historical map below still analyses
  IABC as one fixed-error surrogate — that analysis is *why* it was retired.)
- **BOTH paths are ADOPTED as `radia.open_boundary`** (2026-06-20):
  - **`dtn_cln`** — the exact closed-form separable `Zs`-DtN-CLN: the exact
    eddy/diffusion DtN per multipole, its Cauer ladder (**exact at `n+1` stages**,
    well-conditioned, passive), the Grote-Keller companion auxiliary-ODE form.
  - **`kelvin_dtn`** — the **Kelvin-BUILT material-aware / non-separable** DtN: a
    Kelvin-FEM (`kelvin_fem_radial_dtn`, pure numpy) BUILDS the DtN with no DC floor;
    `kelvin_dtn_matrix` + `steklov_spectrum` (NGSolve) build an arbitrary-shape /
    iron-exterior DtN whose Steklov ladder is point-group-split (cube O_h, square C4v),
    reduced by `band_cln_fit`. **Honest provenance (3 layers):** Kelvin-with-air-exterior
    = Freeman-Lowther 1988/89 (FEMM); transforming `σ` under a conformal map =
    transformation optics / Ward-Pendry 1996 (DC-`σ` cloaks) — both **classical**.
    Sugahara's OWN *validated fusions* are the contributions: (i) Kelvin inversion as an
    exact open boundary with the **`σ`-conformal** transform so a **conductor crosses the
    truncation**, ECT-validated (IEEE Magnetics 2022 — the basis for this module's
    `(a/r)⁴ σ` / `(a/r)² μ` weights); (ii) the Kelvin material-aware DtN as an
    inverse-design kernel (SF-with-iron, conf ~0.83). Full record: MCP
    `kelvin_transformation(topic="material_exterior")`. Goldens `tests/open_boundary/`, usage `docs/open_boundary/demo_dtn_cln_usage.py`.

### Is it "better than PML"? — three senses of *superior* (be precise)

`Zs`/Kelvin-DtN-CLN vs CFS-PML has **no single winner**. Separate the axes:

| axis | winner | why |
|---|---|---|
| **numerical, on the island** (separable / quasi-spherical, MQS evanescent) | **Kelvin / `Zs`-DtN-CLN** (`act6_09`: 8 vs 128 online DOF, cond `~1` vs `~5e4`) | it has the **exact closed-form operator** (Kelvin **builds** it for material / non-separable interiors) — exact, DC-well-conditioned, + a compact passive circuit ROM |
| **arbitrary geometry** | **CFS-PML** (by default) | the exact DtN does **not exist** there → `Zs`-DtN-CLN can't even play; Kelvin pays a **spherical-truncation** waste (Liouville) a box CFS-PML avoids |
| **general adoption** (all shapes / all physics / every solver) | **PML** | PML is *approximate* (has a floor) but **local, trivial, conforms to box/sphere/cylinder shells, ships everywhere** — it won on **generality, not accuracy** |

So: **inside the island Kelvin/`Zs`-DtN-CLN genuinely beats CFS-PML on capability**
(exact + DC-conditioning + passive ROM); **outside it (elongated / arbitrary / waves)
CFS-PML wins**. "PML is superior" means *generally applicable*, **not** *more accurate*.
Do **not** state "`Zs`-DtN-CLN > PML" as a general claim — only the island claim is true.
(Separable ≠ only the sphere: it is the family of separable systems — sphere, cylinder /
circle, half-space, … — but **Kelvin** itself is sphere-locked (3D) / circle-locked (2D)
by Liouville, tighter than separable.)

**Prior art (cite — NOT novel).** The exact rational radiation DtN + local
auxiliary-ODE realisation is **Grote-Keller** (SIAM J. Appl. Math. 1995) /
**Hagstrom-Warburton** (complete radiation BCs / continued-fraction ABCs; **Guddati**
2006, **Birk** 2012); the `√s` (Warburg) diffusion impedance as a **Cauer ladder** is
classical network synthesis (e.g. *PCCP* 18 (2016) 9498); the **Cauer Ladder Network**
MOR is Kameari-Ebrahimi-Sugahara-Shindo-Matsuo (*IEEE T-Magn* 54(3):7201804, 2018).
`radia.open_boundary` is a **verified, reusable operator + the reverse-Bessel
wave↔diffusion unification in CLN form** — *not* a novelty claim. It earns its place by
being correct and useful in the repository (the bonsai), **no paper required**.

---

## 1. One operator behind all of them

Truncate the exterior at a surface `Γ` and impose the exact transparent condition

```
    ∂u/∂n |_Γ  =  Λ_ext u |_Γ ,
```

where `Λ_ext` is the exterior Steklov–Poincaré (Dirichlet-to-Neumann) operator.
Every method realizes some discrete `S_h` on `Γ` targeting this single `Λ_ext`,
in two fidelity classes:

| class | members | behaviour |
|---|---|---|
| **convergent discretizations** | Kelvin, BEM, exact-DtN / CLN, rich infinite elements | `S_h → Λ_ext` for every mode as the discretization refines |
| **fixed-error surrogates** | PML, asymptotic Robin (`S_h=−1/R`), single IABC shell | fixed, mesh-independent modal error (Robin: exact `n=0` only) |

On a sphere `Λ_ext` diagonalizes in spherical harmonics; each mode `n` has a
**scalar frequency symbol** `Λ_n(s)`. That symbol is **reverse-Bessel rational**
— its poles are the Bessel/Thomson filter poles (verified `act6_10_iabc_time_domain`,
`besselap` match `9.9e-16`). The whole map below follows from this one fact.

---

## 2. Three axes that select the method

### 2.1 The kR axis (frequency regime)

The exterior wavenumber sets the character of `Λ_n`:

- **Quasi-static / low-freq (`kR → 0`)** — `Λ_n` is the **real** ladder
  `−(n+1)/R` (3-D), `−n/R` (2-D). Closed by **Kelvin** (spatial) or its
  **CLN** time-domain realization. *This is the SA / Hachinohe regime.*
- **Radiating / wave (finite `kR`)** — `Λ_n(kR)` turns **complex** (spherical
  Hankel log-derivative; `Im` = radiation). Closed by **IABC / PML /
  extended-Kelvin** (`act7_01_highfreq_spectrum_comparison`, `act6_10_iabc_time_domain`).

The same reverse-Bessel structure carries both: **wave** is rational in `s`,
**diffusion** is rational in `√s`, with the **same** poles `roots(θ_n)`
(`act6_02_cln_dtn_cauer`, machine-ε). So one datasheet, two variables.

### 2.2 The geometry axis (truncation shape)

- **Kelvin is sphere-locked.** Its inversion `r ↦ R²/r` is a spherical (3-D) /
  circular (2-D) conformal map — Liouville's theorem leaves no non-sphere Kelvin.
  The **body inside** the sphere is arbitrary (cube, elongated, iron); only the
  **truncation** must be a sphere.
- **CLN / IABC / BEM are surface-free** at the operator level. CLN is
  **finite-exact only on a separable surface** (sphere → spherical Bessel →
  *rational* → `n+1`-stage exact); on a cylinder / cube / arbitrary surface the
  symbol is transcendental / non-separable → CLN is a *convergent* approximation.
- **Aspect ratio** is where the sphere lock bites: an elongated body needs a big
  enclosing sphere (wasted air). A **cylinder** truncation hugs an axial device
  (CLN can use it; Kelvin cannot). In 2-D a conformal pre-map *relocates* the
  cost to spectral order rather than removing it (`memory/conformal_kelvin_2d_no_free_lunch.md`).

### 2.3 The spatial vs temporal axis (division of labor)

The deepest split — **Kelvin and CLN are not substitutes, they are on different
axes**:

- **Kelvin = a SPATIAL operator FACTORY.** It *builds* the exterior DtN as a
  sparse SPD volume FEM — no Green's function, arbitrary body, exterior material
  (iron) — for the case nothing else can (`act1_05_assemble_dtn_matrix/w/q/t/bb`).
- **CLN = a TEMPORAL operator REDUCER.** It *compresses* a DtN's `s`-dependence
  into a small Cauer ladder (auxiliary ODEs) for transient / repeated / control
  solves. It needs a DtN to reduce.
- **They compose:** *separable* exterior → CLN works on the analytic symbol
  directly (no mesh). *Arbitrary geometry + iron, transient* → **Kelvin builds
  the DtN, CLN reduces it**.

---

## 3. The methods side by side

| method | class | freq regime | matrix | geometry | role |
|---|---|---|---|---|---|
| **Kelvin** | convergent | quasi-static | **sparse** SPD (`∝N`) | sphere truncation, **any body + material** | spatial **factory** |
| **BEM** | convergent | DC → wave | **dense** (`∝N²`) | any `Γ`, free-space/Helmholtz kernel | dense factory |
| **PML** | surrogate | radiating (weak DC) | sparse | box / conformal | absorber; **low-freq cond. degrades** |
| **IABC** | surrogate→exact | radiating | sparse shell / Robin+ODE | sphere | electrical-image absorber; exact-DtN route |
| **CLN** | convergent (reduced) | DC → radiating | tiny ladder | surface-free; **finite-exact on sphere** | temporal **reducer** |
| **asymptotic Robin** | surrogate | static | sparse (`−1/R`) | sphere | cheap, exact `n=0` only |

---

## 4. The no-free-lunch axis (finite multipole truncation)

All of them truncate the multipole ladder at some finite `n_max` — this cost is
**universal** and unavoidable:

- **Kelvin**: `p ≥ n` captures mode `n` (the order is the multipole reach).
- **BEM**: surface Nyquist `√N_Γ`.
- **CLN**: `n_max` modes × the per-mode stage count. CLN's "exact" means **exact
  per-mode impedance synthesis over a *finite* multipole set** — not global
  exactness.

The **source's multipole content sets `n_max`**: compact source → low `n_max` →
every method is cheap; sharp / elongated source → high `n_max` → every method
pays (the aspect-ratio penalty). No method escapes this axis; they differ only in
*how they represent the modes they keep* (spatial polynomial / dense surface /
circuit impedance).

---

## 5. Measured anchors (audit-verified in-repo, 2026-06-19)

### Kelvin — spatial, quasi-static (six `mesh_control` pillars)
Exterior volume is free (`‖u_h−P_n‖ ≈ 1e-15`); floor = Curve order (`k=1→3`:
`1.3e-2 → 1.3e-5`); `p ≥ n` & `p`-vs-`h` (~1000× per DOF); optimal `R/a ≈ 2.78`;
re-entrant corner `α_h=0.357 / α_p=0.661`; `ΔDoF ≈ 58 → 1/45` of the interior FEM
error. Holds in 2-D on the `−n/R` ladder. (See `kelvin_transformation(topic="mesh_control")`.)

### IABC / exact-DtN — wave, radiating
- `act6_10_iabc_time_domain`: per-mode DtN = degree-`n` reverse-Bessel rational; its poles ARE the
  **Bessel/Thomson filter poles** (`besselap` `9.9e-16`); the transient auxiliary
  network converges to the exact DtN (`9.8e-11`), all rates `Re<0` (passive).
- `act6_11_exact_dtn_fetd`: exact-DtN FETD reflection falls as **O(h²)** — `5.70e-4 → 3.57e-5`
  (`l=1`, `N=100→400`) — **×2471** better than a 1st-order Sommerfeld ABC, and
  passive (`E_final/E_peak = 1.7e-16`). Separable-geometry (Grote–Keller class).

### CLN — diffusion, temporal
- `act6_02_cln_dtn_cauer`: the diffusion DtN = a Cauer continued fraction **in `√s`**, EXACT at
  `n+1` stages (`NRMSE 1.1e-16`, well-conditioned) **for EVERY multipole, not just the
  dipole** — verified `n=1..6` (each `~1e-16`, spread `1→582`); whereas a Foster fit
  **in `s`** floors at `1.7e-3` and ill-conditions (spread `1e5`) at every `n` — the
  structural win of the natural variable. (The full multipole field on the sphere is a
  **bank** of per-mode exact ladders; non-separable bodies are a band approximation.)
- `act6_04_cln_mor_radial_eddy`: a genuine lab CLN (Lanczos/PVL) reduces a **700-DOF** radial eddy FEM
  to an **`N=16`-stage** integer-order Cauer ladder (reduction `1.2e-5` / total
  `3.6e-4` at the FEM floor) — **~43× state reduction**; `T_N` SPD → real-negative
  poles → unconditionally stable, directly time-domain.
- `act6_05_cln_fetd_reflection`: in a **transient** eddy-current FETD, the CLN open boundary
  (`N=16`, a 16-DOF exterior vs 699 full) gives reflection **`9.9e-7`** (`n=1`) —
  **62659×** better than Dirichlet truncation (`6.19e-2`).
- `act6_09_cln_vs_pml`: CLN vs **CFS-PML** in the DC-to-evanescent eddy-current band — at
  matched accuracy (`NRMSE ≈ 1.46e-4`) CLN uses **8 online DOF vs CFS-PML's 128**
  (~16× fewer) **and** is **`53189×` better-conditioned at DC** (`cond` `1.00` vs
  `5.32e4`); CLN is **DC-exact** (`4.2e-5` vs vanilla-PML `1.83e-2`, whose stretch
  `~1/√s` blows up). CFS-PML's conditioning fix *evaporates* when pushed to CLN's
  accuracy (a thick/strong layer). **Non-claim:** propagating waves are PML's home
  — this is the *diffusion* operator; arbitrary geometry still needs Kelvin + CLN.
- `act7_20_impedance_vs_kelvin_dtn_cln`: **impedance-boundary-DtN-CLN vs Kelvin-DtN-CLN, head-to-head** (4 axes). The
  eddy DtN is `√s`-native (Warburg); Kelvin's `√s`-Cauer is **exact** (`~1e−16`, `n+1`
  stages, well-conditioned, passive) while an impedance **`s`-network** (passive RLC) hits the
  **Warburg wall** (floors `~1.7e−3` + ill-conditions `~1e5`) and the IABC `N`-shell
  transformer ill-conditions `~1`–`1.5` decade/shell. **Kelvin wins on all four axes in the
  MQS/eddy scope**; the impedance route's home is high-freq radiation (outside scope). For the
  sphere the exact DtN is one operator (exact-DtN-Robin `==` Kelvin → same CLN) — the gap is
  the impedance *approximation* vs Kelvin's exact operator.

### Kelvin **builds** + CLN **reduces** — separable → non-separable (the DtN→CLN arc)
The arc that joins the two columns: realize the operator (don't discretize the air),
then compress it.
- `act6_03_dtn_to_cln_wideband`: **separable, band-UNLIMITED.** The 3-D sphere eddy DtN is *exactly*
  rational in `q=√s` (reverse-Bessel), so its CLN (Cauer in `√s`) is exact for all
  `s`: NRMSE **`1.5e-16` over 11 decades** at **2 online DOF**, while FEM-MOR-CLN and
  CFS-PML degrade at the band edge (and adding DOF does not save them).
- `act6_01_kelvin_fem_eddy_dtn`: the Kelvin-FEM **BUILDS** the eddy DtN where there is no closed form
  (radial sphere = analytic-checkable proof). The 3-D ball carries the scalar
  `(R/ρ')²` weight (3-D Dirichlet energy is **not** conformally invariant); with it
  DC error `3e-5` and **no** `(R0/Rfar)^(2n+1)` closure floor (truncated FEM `5.8e-2`).
- `act6_06_square_eddy_dtn_to_cln` / `act6_07_cube_eddy_dtn_to_cln`: the build on a genuinely **NON-separable** body — 2-D
  **square** (C4v) / 3-D **cube** (O_h) — verified **analytic-free by the symmetry
  splitting** of the static Steklov ladder `(S, Mg)`: C4v splits the `m=2` quadrupole
  (square); O_h splits the `l=2` quintet **2+3 = E_g+T_2g** (cube) while the dipole
  stays a degenerate triplet; mesh-convergent; the dipole DtN interpolates
  DC→evanescent and a few-stage CLN-in-`√s` reduces it (6 stages → `3.5e-4`–`6e-4`).
- `act6_08_disk2d_kelvin_eddy_dtn`: the **2-D Kelvin disk is CONFORMAL → NO weight** (2-D in-plane
  Dirichlet energy IS invariant = the Nagamine 2-D-cylindrical in-plane identity
  slots). Hits `−m` to `1e-5` with no DC floor; applying the 3-D `(R/ρ')²` weight in
  2-D **MISSES** the ladder (`m=1` by ~30%) — proving 2-D needs no weight. (2-D `K_m`
  is not exactly rational in `√s`, so its 2-D CLN is a band approximation.)

### Review-hardening: honest refinements (SA / Hachinohe paper Q&A, 2026-06-22)

Twenty reviewer questions + seven new demos sharpened claims this map (and the SA paper)
had stated loosely. The headline is the conditioning correction:

- **Separate the SOLVE cond from the per-mode RATIO.** The "`cond ~1`" / "`1.30`" quoted
  above and in `act6_09` is the **per-mode DtN eigenvalue ratio** (a spectral property),
  **not** the assembled linear-system condition number one actually solves. Measured
  (`act2_14`, 3-D scalar ball): the Kelvin ball's **assembled-stiffness cond is `1e3–1e4`
  and GROWS with refinement**, its `λmax` amplified `~(R/h)²` by the singular centre weight
  `(R/ρ')²`. The PML solve cond (`act7_40`, radial, matched mesh) **grows `~1/k`** toward
  static. So *neither is "1.30"*: **Kelvin pays a fixed (frequency-independent) centre
  penalty; PML loses on frequency-robustness** — Kelvin is frequency-robust but not cheap.
  (3-D only — the 2-D Kelvin disk is conformal / weight-free, no centre penalty.)
- **The `2·min(p,k)` error law is FORM-dependent** (`act3_06`, H(curl)): it extends to edge
  elements **for curved geometry** (`k=2` → the `2k=4` superconvergence, = scalar), but the
  **flat-facet `k=1` is DEGRADED** (`q~1.4 <` the scalar `h²`) and the **vector dipole needs
  higher `p`** than scalar `p=1`. "Kelvin inherits the de Rham family for free" = free of
  *bespoke per-coordinate construction*, **not** of the *curved high-order machinery* the
  vector form genuinely needs.
- **Edge-element `A` solves SPARSE on the Kelvin ball** (`act3_07`): `nnz ∝ N`, no dense DtN
  formed (`(ρ'/R)² → 0` at the centre is integrable + a gauge). The "centre-less DtN
  (FEM-BEM)" is a **niche alternative, not a forced fallback** — "do not form the DtN and it
  stays sparse". (Advantage over dense DtN: asymptotic in storage + qualitative — no
  Green-kernel / BEM assembly, sparse solve, carries iron.)
- **AC / MQS open boundary is FREQUENCY-FLAT** (`act6_13`): the air (`σ=0`) exterior is
  Laplace, so its DtN is the static `−(n+1)/R` at **every** ω; the eddy / skin physics is
  entirely in the conductor (matches analytic DC→skin→evanescent, no DC floor). **Every
  static open-boundary result applies to the MQS (static-apparatus / rotating-machine)
  problem verbatim** — only the interior block carries `jω`. (Measures the "static air ⇒
  Kelvin alone" policy of §0.)
- **The 5–6-digit floor is GEOMETRY, not centre-quadrature** (`act2_14`): it drops with the
  Curve order and is `0.0%`-sensitive to extra centre quadrature (the smooth image `P_n`
  keeps the eigenvalue's centre integral finite; conditioning is the separate matrix effect).
- **Validated beyond symmetry** (`act2_15`): on a genuinely asymmetric 3-source field Kelvin
  recovers the analytic exterior to `~1e-9` and agrees with an independent large-box solve
  (`~3e-8`) — the open-boundary conclusions are not a symmetry artifact.
- **`p`-selection chicken-and-egg, resolved** (`act2_13`): a cheap coarse-`p` solve recovers
  `d_max/R` from the low multipoles; an adaptive loop (estimate the tail rate from the top
  resolved modes, raise `p_c` until it stabilises) sizes the production `p*` — single body
  `p_c=2` nails it, a multi-body mix needs the loop (a 2-mode peek under-sizes).

---

## 6. Selection map (the decision table)

| situation | use | why |
|---|---|---|
| separable + **transient** | **CLN alone** | analytic per-mode symbol → finite-exact ladder; beats PML (`act6_05_cln_fetd_reflection/xx6`) |
| separable + static | analytic DtN **or** Kelvin | either works; Kelvin's exactness here is mostly elegance |
| **arbitrary body + iron + static** | **Kelvin** | the only sparse, Green-function-free DtN factory — the SA / Hachinohe paper |
| arbitrary + iron + **transient** | **Kelvin builds + CLN reduces** | no analytic symbol → Kelvin makes it, CLN compresses it (`act6_06_square_eddy_dtn_to_cln`/`act6_07_cube_eddy_dtn_to_cln`: non-separable square/cube build, symmetry-verified) |
| **radiating** (finite `kR`) | **IABC / PML / extended-Kelvin** | the DtN is complex; static Kelvin is pinned to the real axis |

---

## 7. Runnable layer (source of truth)

- **Module (production API):** [`radia.open_boundary`](../../src/radia/open_boundary/) — the adopted operator, two paths:
  - [`dtn_cln`](../../src/radia/open_boundary/dtn_cln.py) (exact separable): `eddy_dtn` /
    `wave_dtn`, `cauer_ladder` / `eval_ladder` (exact at `n+1` stages), `companion_poles`,
    `sqrt_s_passive_ladder`.
  - [`kelvin_dtn`](../../src/radia/open_boundary/kelvin_dtn.py) (Kelvin-built, material-aware /
    non-separable): `kelvin_fem_radial_dtn` (pure numpy), `kelvin_dtn_matrix` + `steklov_spectrum`
    (NGSolve), `band_cln_fit`.
  - Goldens [`tests/open_boundary/`](../../tests/open_boundary/) (`test_dtn_cln.py`, `test_kelvin_dtn.py`);
    usage [`docs/open_boundary/demo_dtn_cln_usage.py`](demo_dtn_cln_usage.py).
- **Docs:** [`DTN_SPECTRUM_COMPARISON.md`](DTN_SPECTRUM_COMPARISON.md) (**the measured
  method × regime × multipole comparison** — reproduced by `act7_22_dtn_spectrum_consolidated`),
  [`docs/kelvin/DTN_SPECTRUM_COARSE_MESH.md`](../kelvin/DTN_SPECTRUM_COARSE_MESH.md)
  (Kelvin spectral datasheet), [`docs/cln/CAUER_LADDER_NETWORK.md`](../cln/CAUER_LADDER_NETWORK.md)
  (CLN), [`docs/kelvin/KELVIN_TRANSFORMATION.md`](../kelvin/KELVIN_TRANSFORMATION.md).
- **MCP knowledge:** `dtn_coarse_mesh(topic=...)` — including `topic="dtn_to_cln"`
  (the exact-impedance / `Zs` → DtN → CLN realisation, relocated 2026-06-20 from the
  retired `iabc` tool) and `topic="method_map"` (the IABC comparison record) —
  `kelvin_transformation(topic="mesh_control")`, `mor_cln(...)`.
- **Maintained implementation/validation:** `src/radia/open_boundary/{dtn_cln.py,kelvin_dtn.py}`,
  `src/radia/infinite_element.py`, and `validation_test/open_boundary/`.
  The former wave/diffusion research corpus is routed through
  [`docs/kelvin/ARCHIVE_RETIREMENT.md`](../kelvin/ARCHIVE_RETIREMENT.md);
  production docs should cite the maintained API/tests rather than the retired
  full-source ledger.

---

## 8. The map in one line

> **One operator (`Λ_ext`, reverse-Bessel rational), selected by three axes:**
> *frequency* (quasi-static → Kelvin / radiating → PML/BEM, out of radia's MQS scope), *geometry* (sphere-locked
> Kelvin / surface-free CLN), and *space-vs-time* (Kelvin **builds** the operator,
> CLN **reduces** it). The "exactness" of any of them is always relative to a
> **finite multipole truncation** — there is no free lunch on that axis; the
> methods differ in *capability* (Kelvin builds arbitrary-body+material) and
> *deployment* (CLN compresses to a time-domain ladder), not in escaping the
> modal cost.
