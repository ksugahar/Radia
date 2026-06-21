# Open-boundary methods, measured on the DtN spectrum

The quantitative companion to [`OPEN_BOUNDARY_MAP.md`](OPEN_BOUNDARY_MAP.md): every
open-boundary closure compared on **one yardstick — the per-multipole Dirichlet-to-Neumann
(DtN) spectral defect**

```
    d_n  =  | λ_h(n) − λ_exact(n) |  /  | λ_exact(n) |
```

across the three regimes, with the measured numbers reproduced in-repo by
[`examples/kelvin_transformation/DtN_spectrum/act7_22_dtn_spectrum_consolidated.py`](../../examples/kelvin_transformation/DtN_spectrum/act7_22_dtn_spectrum_consolidated.py)
(data: `act7_22_dtn_spectrum_consolidated.json`; figure: `…​.pdf`). It consolidates the
scattered per-regime measurements — static [`act7_21`](../../examples/kelvin_transformation/DtN_spectrum/act7_21_lowfreq_openbc_4way.py),
eddy [`act6_09`](../../examples/kelvin_transformation/DtN_spectrum/act6_09_cln_vs_pml.py),
high-frequency `act7_01`/`act7_07` — into one table.

## Why the DtN-spectral defect (the lens)

A **field-error** comparison conflates the *interior FEM error* with the *open-boundary
error*. The **per-multipole DtN defect isolates the boundary operator's accuracy, mode by
mode** — a diagnostic the usual comparisons cannot give. The exact eigenvalue is a property
of the *continuous* exterior operator `Λ_ext`, regime by regime:

| regime | exact `λ_exact(n)` | operator |
|---|---|---|
| static (Laplace) | `−(n+1)/a` (3D), `−n/a` (2D) | real ladder |
| eddy (diffusion) | rational in `q=√s` (reverse-Bessel) | `radia.open_boundary.eddy_dtn` |
| high-freq (Helmholtz) | `z h_n^(1)′(z)/h_n^(1)(z)`, `z=ka` | `radia.open_boundary.wave_dtn` (complex/radiating) |

## The measured comparison (a = 1)

| method | accuracy (per-mode `d_n`) | convergent? | parameter-free? | DC conditioning | cost |
|---|---|---|---|---|---|
| **Kelvin** | **exact** (static `~2e−6`, eddy converges; **extended/radiating Kelvin** carries high-freq — HOIBC `~3e−2`, exact-Z `~6e−6`) | **YES** | **YES** | **flat** (`~(R/h)²` solve) | sparse |
| **BEM** (`ngsolve.bem`) | **exact** (all regimes; the real `ngsolve.bem` Helmholtz BEM reproduces `wave_dtn` to `~1e−5`, `act7_23`) | **YES** | **YES** | — | **DENSE** (`N²`) |
| **PML** (`NGSolve` FEM+PML) | accurate per-mode (`~1e−4` in its home; the real 3-D NGSolve FEM+PML reproduces `wave_dtn` to `~1e−3`, `act7_24`) | no (a tuned **layer**) | no (`σ, L`) | **blows up** (`2.4e4`@DC) | sparse |
| **CFS-PML** | accurate per-mode | no (a tuned layer) | no (`σ, α, L`) | **fixed** (`2.3e3`@DC) | sparse |
| **ballooning** (truncation wall) | **fails the LOW modes** (`~(a/R)^{2n+1}`) | no (finite **reach** `R`) | `R` | — | sparse |
| **infinite element** (Bettess) | **exact `n ≤ P−1`**, degrades `n ≥ P` (decay-basis; `act7_25`) | **YES** (p in decay order) | `P` | — | sparse |
| **Robin** (`λ=−1/a`) | **exact `n=0` only**, fails HIGH modes | no (a fixed floor) | YES | — | sparse |

> **Conditioning — the honest distinction (review-hardening 2026-06-22).** The "DC conditioning"
> column is the per-mode / DC-frequency BEHAVIOUR (does it blow up at DC?). The assembled
> linear-system SOLVE cond is a separate *magnitude*: Kelvin is frequency-flat but at `~(R/h)²`
> (the singular centre weight `(R/ρ')²` amplifies `λ_max`; assembled cond `1e3–1e4`, `act2_14`),
> PML grows `~1/k` toward static (`act7_40`). So PML loses on frequency-robustness, but Kelvin is
> frequency-robust **not cheap** — and neither figure is the per-mode `1.30`. (3-D only; the 2-D
> Kelvin disk is conformal / weight-free.)

### Per-regime numbers (from `act7_22`)

- **static** (`λ_n=−(n+1)`): Kelvin `≤ 2.1e−6` every mode (converged); **ballooning** `n=0..4 =
  0.33 → 0.024 → 1.6e−3 → 1.1e−4 → 6.9e−6` (fails the *low* slow-decaying modes, shrinks with
  `n` and with `R`); **Robin** `0 → 0.5 → 0.67 → 0.75 → 0.80` (exact `n=0`, fails the *high*
  modes — the opposite failure); the **infinite element** (decay order `P`) is **exact for
  `n ≤ P−1`** (its `(a/r)^k` basis contains the exact `r^{−(n+1)}`) and degrades for `n ≥ P` — the
  OPPOSITE failure to the wall, so it is the *convergent* member of the ballooning / infinite-element
  family (`act7_25`; NGSolve has no infinite element — implemented directly).
- **eddy** (`s=i·1`): all three closures resolve the mode (`d_n < 5e−2`); the **distinguishing**
  behaviour is **convergence** — Kelvin-built `n=2` defect drops `1.2e−3 → 2.2e−5` under
  `(h, R_mid)` refinement (parameter-free), while PML/CFS-PML are a tuned absorbing layer; and
  **conditioning** — vanilla PML's interior matrix conditioning **blows up toward DC**
  (`cond ≈ 2.4e4` at `s=i·1e−4`) while **CFS-PML removes it** (`cond ≈ 2.3e3`) — *which is why
  CFS-PML exists*.
- **high-freq** (`z=ka=2`, radiating) — **a studied regime** (the radiating extended-Kelvin /
  HOIBC / PML track, `act7_01`–`act7_07`; Sugahara, *IEICE Trans. C* 2024). The DtN goes
  **COMPLEX** (`Im` = radiation). The **static** Kelvin is only the `kR→0` limit (real axis); the
  regime is carried by the **extended (radiating) Kelvin** — transformation-optics medium +
  matched HOIBC — which reproduces the complex DtN to **`~6e−6` with the exact impedance** and
  **`~3e−2` with the 2nd-order HOIBC** (the radiating-band knee at `n≈ka`, measured in `act7_22`),
  competitive with **PML** (`d_n ~ 1e−4` in its home) and exact **BEM-FEM**. *(The Laplace-kernel /
  MQS limit is on Radia's CORE field solver — MMM/MSC — not on this open-boundary study.)*

## The headline (the two classes, measured)

- **Convergent + parameter-free**: **Kelvin** (static / eddy; the **extended / radiating Kelvin**
  carries the complex DtN at high-freq via the matched HOIBC — `act7_05`/`act7_07`/`act7_22`) and
  **BEM** (all regimes, but DENSE) — the discrete operator `S_h → Λ_ext` for every mode under
  refinement.
- **Fixed-error surrogates**: **PML** (accurate per-mode but **DC-ill-conditioned** + tuned),
  **CFS-PML** (fixes the conditioning at modest accuracy + tuned), **Robin** (exact `n=0` only).
- **Finite-reach**: **ballooning** (cheap, but a finite wall — fails the low modes; shrinks with
  `R`).
- **High-freq / radiating is studied, not excluded**: the DtN goes **complex**; carried by the
  **extended (radiating) Kelvin** (matched HOIBC, Sugahara *IEICE Trans. C* 2024), **PML**, and
  **BEM-FEM**. The Laplace-kernel / MQS limit is on Radia's *core field solver* (MMM/MSC), not on
  this open-boundary comparison.

> The comparison is **not a single number**: on the *accuracy* axis Kelvin, BEM and PML are all
> good per-mode; the methods separate on **convergence** (Kelvin/BEM vs tuned layer),
> **conditioning** (PML blows up at DC; CFS-PML fixes it; Kelvin frequency-flat but at `~(R/h)²`), **cost** (BEM dense vs
> the rest sparse), and **reach** (ballooning finite). Reporting all axes is the honest result —
> not "method X is best".

## Relation to the conventional reflection coefficient (`d_n ≡ reflection`)

The community grades an open boundary by its **reflection coefficient** `R_n` (Bérenger;
Engquist–Majda; Bayliss–Turkel), not by a DtN defect. The two are the **same quantity**: for a
mode at the truncation the spurious "wrong" solution is the *growing* mode (static / evanescent)
or the *incoming* wave (high-freq), with DtN `λ_other`, and a boundary DtN `λ_h` admits

```
    R_n  =  | λ_h − λ_exact |  /  | λ_h − λ_other |
```

— the **same numerator** as `d_n = |λ_h − λ_exact| / |λ_exact|`. So `R_n` is the
physically-measured *face* of the DtN defect; they carry identical information (only the
normalisation differs). Measured in `act7_22` (reflection view): static **Kelvin** is
reflectionless (`R ~ 1e−13…1e−6`) while **ballooning** reflects the low modes (`R[n=0] = 0.25`,
shrinking with `n`); high-freq propagating (`n ≤ ka`) **extended-Kelvin-HOIBC** and **PML** are
both low-reflection (`R ≲ 3e−2`); and the real **`ngsolve.bem`** (`act7_23`, `R ~ 1e−5`) and the
real **3-D NGSolve FEM+PML** (`act7_24`, `R ~ 1e−3`) are ~reflectionless because they reproduce the
exact outgoing DtN. **The DtN-spectrum view adds nothing physically new
over reflection** — its only convenience is being defined *uniformly across regimes* (static /
evanescent have no propagating wave to "reflect"). This is a measurement / comparison, framed in
the standard reflection language, **not a new metric**.

## De Rham / vector extension, and the two infinite-element families

The scalar infinite element (`act7_25`, the `H1` / 0-form end) extends to a **de Rham
(exact-sequence) infinite element** — the `H1 → H(curl) → H(div) → L2` complex on the exterior, with
the radial decay families shifted by **+1 per form degree** (`S0 = {n+1..n+P}`, `S1 = S0+1`,
`S2 = S0+2`, `S3 = S0+3`) chosen so that grad / curl / div **commute** (the commuting diagram).
[`act7_26_derham_infinite_element`](../../examples/kelvin_transformation/DtN_spectrum/act7_26_derham_infinite_element.py)
verifies the diagram **exactly with sympy** (`grad(V0) ⊂ V1`, `curl(V1) ⊂ V2`, `div(V2) ⊂ V3` with
explicit structure constants; `curl∘grad = 0`, `div∘curl = 0`), the toroidal / div closure resting on
the Legendre angular eigenvalue `−n(n+1)`. This is the construction of **Demkowicz & Pal**, *An
infinite element for Maxwell's equations* (CMAME 164, 1998) — a known construction, implemented and
verified here (no novelty claimed). "High-order" is the same decay order `P` (a radial p-method,
exact for `n ≤ P−1`).

**Two distinct "infinite element" families** — they reach de Rham by different routes:

| family | radial treatment | de Rham source | example |
|---|---|---|---|
| **decay-basis / shape-function** | special radial functions `(a/r)^k` | hand-built **commuting radial families** (exact-sequence) | Bettess scalar (`act7_25`); Demkowicz–Pal vector (`act7_26`) |
| **coordinate-mapping / coordinate-scaling** | a real coordinate stretch on a standard FE layer | **standard Nédélec / RT elements** on the mapped region — inherited **for free** | the **Kelvin transformation** (conformal map); the coordinate-scaling infinite-element domain shipped in commercial FE |

The **practical / shipped** de Rham open boundary is the coordinate-mapping kind (the same family as
the Kelvin transformation, and the real-stretch cousin of PML); the decay-basis exact-sequence IE is
the academic alternative. For a VECTOR problem handled by a **volume integral** method instead
(MMM / MSC, or the H(div) charge-Gram VIM), the open boundary lives in the free-space kernel — there
is no exterior to mesh, so the IE-vs-Kelvin choice does not arise. That integral route is an
**alternative to**, not a combination with, the mapped-exterior FEM.

**Trefftz framing.** The whole *convergent* class = "represent the exact exterior" = the **Trefftz
umbrella**: Kelvin (a conformal map of the exact exterior), BEM (the exact kernel), and the
decay-basis IE (whose basis `(a/r)^k` contains the one Trefftz function `r^{−(n+1)}`, hence exact for
`n ≤ P−1`). A *pure* Trefftz element (exact multipoles as the basis) reproduces the exact DtN
trivially. The trade is **no free lunch**: sparse ⇒ approximate (Kelvin's FE floor; the IE's order-`P`
cutoff), exact ⇒ dense (BEM / pure-Trefftz multipole).

**Build decision (Gate 1, [`act7_27`](../../examples/kelvin_transformation/DtN_spectrum/act7_27_ie_vs_kelvin_vs_pml_gate1.py)).**
Before committing to a C++ infinite element, the honest rival is **not Kelvin but box-PML** —
NGSolve's `pml.Cartesian` / `BrickRadial` *also* escapes the Liouville sphere-lock. Measured:
on exterior-DOF vs aspect ratio `AR = L/d`, **Kelvin scales `AR²`** (sphere-lock — it must enclose
the body in a sphere) while **box-PML and IE both scale `AR¹`** (the IE a bit leaner). So Kelvin is
out for elongated / planar bodies, but the IE only **ties** box-PML on geometry. The IE's one proven
unique edge is **spectral exactness** (`n ≤ P−1`); however the naive reciprocal-power basis `(a/r)^k`
is **Hilbert-ill-conditioned** (`cond` ≈ `10 → 5e3 → 4e6 → 4e9` for `P = 2,4,6,8`), so a production IE
needs an *orthogonalized* basis, and a cheap 1-D proxy does **not** establish an IE conditioning win
over box-PML. **Verdict: no clean GO** — prefer box-PML (complement NGSolve, do not reimplement),
*unless* the IE basis is orthogonalized **and** IE-vs-box-PML is settled on a real 3-D model first.

**Fair re-test on the DtN yardstick ([`act7_28`](../../examples/kelvin_transformation/DtN_spectrum/act7_28_ie_vs_kelvin_fair_dtn.py)).**
Gate 1 above compared the IE against box-PML but never put IE and **Kelvin** on the same per-mode
`d_n`. Doing so reveals they are the **same method on the sphere**: the Kelvin inversion `ξ = a²/r`
maps `r^{−(n+1)}` to the polynomial `ξ^{n+1}`, and the IE decay-matrix `A_kl = a(kl+n(n+1))/((k+l)−1)`
*is* the Kelvin-mapped monomial energy (measured identical). Consequences, measured: IE and Kelvin
give **identical DtN at matched DOF** (a TIE on accuracy-per-DOF; both exact for `n ≤ P−1`); the IE's
*only* deficit is that the monomial basis `(a/r)^k` is the Hilbert-ill-conditioned coordinate system
for that space (`cond ≈ 10 → 4e9` for `P=2..8`), whereas an orthogonal/nodal basis for the **same**
space is well-conditioned (`cond ≈ 2 → 339`) — so the deficit is **fixable by orthogonalization, not
intrinsic**. Net: an *orthogonalized* IE has Kelvin-grade accuracy + conditioning **plus** the
geometry edge (no Liouville sphere-lock) — so the honest build decision reads in the IE's favour
*provided the orthogonalized basis is used*. (The Gate-1 "prefer box-PML" was an unfair-comparison
artifact; box-PML stays relevant only as the tuned, wave-native alternative for the eddy/radiating
regime.) See [`INFINITE_ELEMENT_SOTA.md`](INFINITE_ELEMENT_SOTA.md) for the full infinite-element
state of the art (acoustic conjugated / Jacobi conditioning / Hardy-space; the de Rham Maxwell IE)
and the Gate-2/3 build spec.

## Provenance (cite — this is a measurement / comparison, not a new method)

Kelvin / inversion open boundary: Freeman-Lowther 1988/89; Brunotte-Meunier-Imhoff 1992.
Ballooning / infinite elements: Silvester-Hsieh 1971; Bettess. PML: Bérenger 1994;
Collino-Monk. CFS-PML: Kuzuoğlu-Mittra 1996. BEM-FEM open boundary: standard. The contribution
here is the **unified method × regime × multipole DtN-spectral comparison on one yardstick**.
